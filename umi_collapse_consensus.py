#!/usr/bin/python3
"""
UMI collapse and consensus (z-score filter, k = 1).

- Group reads by UMI (5′ or 3′, user-specified).
- If a UMI has exactly 1 sequence originally, output that sequence unchanged.
- Keep reads with |len - mean| < 1 * SD (strict).
- For exactly 2 reads: if their length span > --length-threshold, drop both.
- Keep UMI only if it has at least --min-reads before and after filtering,
  except for original single-sequence UMIs, which are written as-is.
- Skip MAFFT if kept reads share the same length; otherwise run MAFFT.
- Majority-rule consensus per column (ties -> 'N'); remove gaps.
- Drop consensus if it contains more than --max-Ns Ns.
- Log all drop reasons to UMI.dropouts.txt.
"""

from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
import os
import argparse
from io import StringIO
import contextlib
from concurrent.futures import ProcessPoolExecutor
from statistics import mean, pstdev


def parse_args():
    p = argparse.ArgumentParser(
        description="UMI collapse + consensus (z-score, k=1) with configurable UMI side/length"
    )
    p.add_argument("-f", "--fastafile", required=True,
                   help="Input FASTA (UMI extracted from 5' or 3' end)")
    p.add_argument("-o", "--filename", required=True,
                   help="Output consensus FASTA")
    p.add_argument("-s", "--summaryFile", required=True,
                   help="UMI → read count summary")
    p.add_argument("-N", "--NumberWorkers", dest="Nworkers", type=int, default=os.cpu_count(),
                   help="Parallel workers (default: all CPUs)")
    p.add_argument("-m", "--min-reads", dest="min_reads", type=int, default=2,
                   help="Minimum reads per UMI (checked before and after filtering; original singleton UMIs are still kept)")
    p.add_argument("-t", "--length-threshold", dest="len_thresh", type=int, default=5,
                   help="Two-read guard: if span > this (nt), drop both")
    p.add_argument("--max-Ns", dest="max_Ns", type=int, default=10,
                   help="Maximum allowed 'N' characters in consensus; if exceeded, drop (default: 10)")
    p.add_argument("--umi-side", dest="umi_side", choices=["5p", "3p"], default="3p",
                   help="Side to extract UMI from (5p or 3p) [default: 3p]")
    p.add_argument("--umi-len", dest="umi_len", type=int, default=14,
                   help="Length of UMI region [default: 14]")
    p.add_argument("-v", "--verbose", action="count", default=0,
                   help="-v, -vv for more logs")
    args = p.parse_args()
    if args.min_reads < 1:
        p.error("--min-reads must be ≥ 1")
    if args.len_thresh < 0:
        p.error("--length-threshold must be ≥ 0")
    if args.max_Ns < 0:
        p.error("--max-Ns must be ≥ 0")
    if args.umi_len < 1:
        p.error("--umi-len must be ≥ 1")
    return vars(args)


def write_umi_summary(filename, groups):
    with open(filename, "w") as fh:
        for umi, seqs in groups.items():
            fh.write(f"{umi}: {len(seqs)} sequences\n")


def log_dropout(umi, reason):
    with open("UMI.dropouts.txt", "a") as fh:
        fh.write(f"{umi}: {reason}\n")


def filter_by_length_zscore_k1(sequences, len_thresh):
    if not sequences:
        return [], []

    lengths = [len(s) for s in sequences]
    span = max(lengths) - min(lengths)

    if len(sequences) == 2 and span > len_thresh:
        removed = list(zip(sequences, lengths))
        return [], removed

    mu = mean(lengths)
    sigma = pstdev(lengths)
    if sigma == 0:
        return sequences[:], []

    kept, removed = [], []
    for s, L in zip(sequences, lengths):
        if abs(L - mu) < sigma:
            kept.append(s)
        else:
            removed.append((s, L))
    return kept, removed


def calculate_consensus(align):
    aln_len = align.get_alignment_length()
    out = []
    for i in range(aln_len):
        symbols = [rec[i] for rec in align]
        counts = {}
        for x in symbols:
            counts[x] = counts.get(x, 0) + 1
        base = max(counts, key=counts.get)
        if counts[base] > len(symbols) / 2:
            out.append(base)
        else:
            out.append("N")
    return Seq("".join(out).upper().replace("-", ""))


def align_and_consensus(task):
    umi, sequences, min_reads, len_thresh = task

    if len(sequences) == 1:
        return (umi, sequences[0])

    if len(sequences) < min_reads:
        log_dropout(umi, f"fewer than min_reads before filtering (n={len(sequences)} < {min_reads})")
        return None

    kept, removed = filter_by_length_zscore_k1(sequences, len_thresh)
    if removed:
        removed_lens = ",".join(str(L) for (_, L) in removed)
        log_dropout(umi, f"removed {len(removed)} read(s) by length (k=1); removed lengths: {removed_lens}")

    if len(kept) < min_reads:
        log_dropout(umi, f"fewer than min_reads after filtering (n={len(kept)} < {min_reads})")
        return None

    kept_lengths = [len(s) for s in kept]
    if max(kept_lengths) == min(kept_lengths):
        seq_string = "".join(f">seq{i}\n{s}\n" for i, s in enumerate(kept))
        align = AlignIO.read(StringIO(seq_string), "fasta")
    else:
        tmp = f"aligned_seqs_{os.getpid()}.fasta"
        with open(tmp, "w") as fh:
            for i, s in enumerate(kept):
                fh.write(f">Seq_{i}\n{s}\n")
        mafft = MafftCommandline(input=tmp)
        stdout, stderr = mafft()
        with contextlib.suppress(FileNotFoundError):
            os.remove(tmp)
        align = AlignIO.read(StringIO(stdout), "fasta")

    consensus = calculate_consensus(align)
    return (umi, str(consensus))


def executor_map(executor, fn, iterable, chunksize=1):
    return executor.map(fn, iterable, chunksize=chunksize)


def main(args):
    fasta = args["fastafile"]
    out_fa = args["filename"]
    umi_summary = args["summaryFile"]
    max_workers = args["Nworkers"]
    min_reads = args["min_reads"]
    len_thresh = args["len_thresh"]
    max_Ns = args["max_Ns"]
    verbose = args["verbose"]
    umi_side = args["umi_side"]
    umi_len = args["umi_len"]

    if os.path.exists("UMI.dropouts.txt"):
        os.remove("UMI.dropouts.txt")

    groups = {}
    with open(fasta, "r") as fh:
        for i, line in enumerate(fh):
            if verbose and i % 2000 == 0:
                print(f"Processing line {i}")
            if line.startswith(">"):
                continue
            seq = line.strip()
            if not seq:
                continue
            if len(seq) < umi_len:
                log_dropout(seq, f"sequence shorter than UMI length ({umi_len})")
                continue
            if umi_side == "5p":
                umi = seq[:umi_len]
            else:
                umi = seq[-umi_len:]
            groups.setdefault(umi, []).append(seq)

    if verbose:
        print(f"Writing UMI summary to {umi_summary}")
    write_umi_summary(umi_summary, groups)

    tasks = [(umi, seqs, min_reads, len_thresh) for umi, seqs in groups.items()]

    if verbose:
        print(
            f"Consensus with {max_workers} worker(s) | min_reads={min_reads} | z-score k=1 | "
            f"2-read span guard={len_thresh} nt | umi_side={umi_side} | umi_len={umi_len} | max_Ns={max_Ns}"
        )

    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        results = list(executor_map(ex, align_and_consensus, tasks))

    kept = 0
    written = 0
    if verbose:
        print(f"Writing consensus FASTA to {out_fa}")
    with open(out_fa, "w") as out:
        for r in results:
            if not r:
                continue
            kept += 1
            umi, seq = r
            n_count = seq.upper().count("N")
            if n_count > max_Ns:
                log_dropout(umi, f"consensus has >{max_Ns} 'N's (N={n_count}); removed")
                continue
            out.write(f">UMI_{umi}\n{seq}\n")
            written += 1

    if verbose:
        print(f"UMIs passing filters (pre-N filter): {kept} / {len(groups)}")
        print(f"UMIs written (post N filter): {written} / {len(groups)}")


if __name__ == "__main__":
    args = parse_args()
    main(args)
