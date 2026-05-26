#!/usr/bin/env python3
"""Score Dorado FASTQ output from make_r10_truthset_fast5.py."""

from __future__ import annotations

import argparse
import csv
import pathlib


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Score R10 truth-set basecalls.")
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--expected-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--fastq", type=pathlib.Path, required=True)
    parser.add_argument("--out-tsv", type=pathlib.Path, required=True)
    parser.add_argument("--k", type=int, default=15)
    return parser.parse_args()


def read_fasta(path: pathlib.Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def read_fastq(path: pathlib.Path) -> dict[str, str]:
    reads: dict[str, str] = {}
    with path.open() as handle:
        while True:
            header = handle.readline().rstrip()
            if not header:
                break
            seq = handle.readline().rstrip().upper()
            handle.readline()
            handle.readline()
            read_id = header[1:].split()[0]
            reads[read_id] = seq
    return reads


def kmerset(seq: str, k: int) -> set[str]:
    return {seq[i : i + k] for i in range(0, max(0, len(seq) - k + 1))}


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def containment(query: str, target: str, k: int) -> float:
    query_kmers = kmerset(query, k)
    if not query_kmers:
        return 0.0
    target_kmers = kmerset(target, k)
    return len(query_kmers & target_kmers) / len(query_kmers)


def main() -> None:
    args = parse_args()
    expected = read_fasta(args.expected_fasta)
    reads = read_fastq(args.fastq)
    rows = []
    with args.manifest.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            basecall = reads.get(row["read_id"], "")
            expected_seq = expected[row["expected_name"]]
            original_seq = expected["original"]
            scores = {
                "expected_forward_kmer_containment": containment(basecall, expected_seq, args.k),
                "expected_revcomp_kmer_containment": containment(
                    basecall, revcomp(expected_seq), args.k
                ),
                "original_forward_kmer_containment": containment(basecall, original_seq, args.k),
                "original_revcomp_kmer_containment": containment(
                    basecall, revcomp(original_seq), args.k
                ),
            }
            rows.append(
                {
                    **row,
                    "basecalled": "yes" if basecall else "no",
                    "basecall_len": len(basecall),
                    **{key: f"{value:.4f}" for key, value in scores.items()},
                }
            )

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    ranked = sorted(
        rows,
        key=lambda row: float(row["original_forward_kmer_containment"])
        + float(row["original_revcomp_kmer_containment"]),
        reverse=True,
    )
    print(f"Wrote scores to {args.out_tsv}")
    print("Top variants by original-sequence k-mer containment:")
    for row in ranked[:8]:
        print(
            row["variant"],
            "len=" + str(row["basecall_len"]),
            "orig_fwd=" + row["original_forward_kmer_containment"],
            "orig_rc=" + row["original_revcomp_kmer_containment"],
            "expected_fwd=" + row["expected_forward_kmer_containment"],
            "expected_rc=" + row["expected_revcomp_kmer_containment"],
        )


if __name__ == "__main__":
    main()
