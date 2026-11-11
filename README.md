# UMI Consensus Sequence Generator

This script groups sequencing reads by unique molecular identifier (UMI), filters reads within each UMI group by length, and generates one consensus sequence per retained UMI group.

## Overview

Reads are grouped according to a UMI extracted from either the 5′ or 3′ end of each sequence. For each UMI group, the script:

1. Checks that the group has at least the required minimum number of reads.
2. Filters reads by sequence length:
   - Reads are retained only when their length is within 1 standard deviation of the group mean.
   - For groups containing exactly two reads, both reads are discarded if their length difference exceeds the user-defined threshold.
3. Aligns retained reads:
   - If all retained reads have the same length, MAFFT is skipped.
   - If retained reads differ in length, MAFFT is used for multiple-sequence alignment.
4. Builds a consensus sequence using a majority-rule approach:
   - At each alignment position, the most frequent symbol is selected if it represents more than half of the reads.
   - If no symbol has a strict majority, `N` is assigned.
   - Gap characters are removed from the final consensus.
5. Removes consensus sequences containing more than the allowed number of `N` characters.
6. Writes dropout reasons to `UMI.dropouts.txt`.

## Requirements

- Python 3
- Biopython
- MAFFT

## Arguments
| Argument                   | Description                                                                                       | Default            |
| -------------------------- | ------------------------------------------------------------------------------------------------- | ------------------ |
| `-f`, `--fastafile`        | Input FASTA file                                                                                  | Required           |
| `-o`, `--filename`         | Output FASTA file containing consensus sequences                                                  | Required           |
| `-s`, `--summaryFile`      | Output file containing read counts for each UMI                                                   | Required           |
| `-N`, `--NumberWorkers`    | Number of parallel worker processes                                                               | All available CPUs |
| `-m`, `--min-reads`        | Minimum number of reads required per UMI before and after filtering                               | `2`                |
| `-t`, `--length-threshold` | For exactly two-read UMI groups, discard both reads if their length difference exceeds this value | `5`                |
| `--max-Ns`                 | Maximum number of `N` bases allowed in a consensus sequence                                       | `10`               |
| `--umi-side`               | Side from which to extract the UMI: `5p` or `3p`                                                  | `3p`               |
| `--umi-len`                | Length of the UMI region                                                                          | `14`               |
| `-v`, `--verbose`          | Enable verbose logging; use `-v` or `-vv` for more output                                         | Off                |


##Usage/Example

python umi_collapse_consensus.py \
  -f reads.fasta \
  -o umi_consensus.fasta \
  -s umi_read_counts.txt \
  --umi-side 3p \
  --umi-len 14 \
  --max-Ns 10 \
  -m 2 \
  -t 5 \
  -N 8 \
  -v
 
##Dropout Log

Dropped UMIs and removal reasons are written to: UMI.dropouts.txt
Possible dropout reasons include:
sequence shorter than the configured UMI length
fewer than the minimum required reads before filtering
reads removed by length filtering
fewer than the minimum required reads after filtering
consensus sequence containing more than the allowed number of N bases
