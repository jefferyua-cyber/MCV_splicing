#!/usr/bin/env python3

import gzip
import os
import glob
import csv
from itertools import islice

FASTQ_FOLDER = "/Users/shudalab/Library/CloudStorage/OneDrive-UniversityofPittsburgh/RNAseq/260107 nuclear vs cytoplasmic/Shuda, M_01_07_2026_Total+mRNA_fastq/mRNA"

HALF_LEADER_SEQUENCE = "GGGAGGGAGACGGAAGACTCTTAACTTTTTTTCAACAAGG"
JUNCTION_SEQUENCE    = "GGGAGGAGGACATTAAAAGAGGTTTCTGGCATTGACTCAT"

MAX_MISMATCHES = 2

OUTPUT_CSV_PATH = "sequence_abundance_summary_junction_half_leader.csv"

R1_PATTERN = "*_R1_*.fastq.gz"

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def fastq_sequences(path):
    with gzip.open(path, "rt") as fh:
        while True:
            lines = list(islice(fh, 4))
            if not lines or len(lines) < 4:
                break
            yield lines[1].strip()


def hamming_leq(a: str, b: str, max_mismatches: int) -> bool:
    mismatches = 0
    for x, y in zip(a, b):
        if x != y and x != "N" and y != "N":
            mismatches += 1
            if mismatches > max_mismatches:
                return False
    return True


def contains_approx(read: str, target: str, max_mismatches: int) -> bool:
    k = len(target)
    if len(read) < k:
        return False
    for start in range(len(read) - k + 1):
        if hamming_leq(read[start:start + k], target, max_mismatches):
            return True
    return False


def read_matches(seq, target, target_rc, max_mismatches):
    seq = seq.upper()
    return (contains_approx(seq, target, max_mismatches)
            or contains_approx(seq, target_rc, max_mismatches))


def count_pair_file(r1_path, r2_path, targets, max_mismatches):
    total_pairs = 0
    match_counts = {label: 0 for label in targets}

    for seq1, seq2 in zip(fastq_sequences(r1_path), fastq_sequences(r2_path)):
        total_pairs += 1
        for label, (target, target_rc) in targets.items():
            m1 = read_matches(seq1, target, target_rc, max_mismatches)
            m2 = read_matches(seq2, target, target_rc, max_mismatches)
            if m1 or m2:
                match_counts[label] += 1

    return total_pairs, match_counts


def find_r1_files(folder, pattern):
    return sorted(glob.glob(os.path.join(folder, pattern)))


def main():
    targets = {
        "HALF_LEADER": (HALF_LEADER_SEQUENCE.upper(), reverse_complement(HALF_LEADER_SEQUENCE.upper())),
        "JUNCTION":    (JUNCTION_SEQUENCE.upper(),    reverse_complement(JUNCTION_SEQUENCE.upper())),
    }

    r1_files = find_r1_files(FASTQ_FOLDER, R1_PATTERN)

    if not r1_files:
        print(f"No R1 files matching {R1_PATTERN} found in:\n  {FASTQ_FOLDER}")
        return

    print(f"Found {len(r1_files)} R1 file(s) in:\n  {FASTQ_FOLDER}\n")

    rows = []
    for r1_path in r1_files:
        r2_path = r1_path.replace("_R1_", "_R2_")
        r1_name = os.path.basename(r1_path)
        r2_name = os.path.basename(r2_path)

        if not os.path.exists(r2_path):
            print(f"WARNING: no matching R2 found for {r1_name} (expected {r2_name}) -- skipping.\n")
            continue

        print(f"Processing pair: {r1_name} + {r2_name} ...")
        total_pairs, match_counts = count_pair_file(r1_path, r2_path, targets, MAX_MISMATCHES)

        half_leader_count = match_counts["HALF_LEADER"]
        junction_count = match_counts["JUNCTION"]
        half_leader_cpm = (half_leader_count / total_pairs * 1e6) if total_pairs else 0.0
        junction_cpm = (junction_count / total_pairs * 1e6) if total_pairs else 0.0
        ratio = (junction_count / half_leader_count) if half_leader_count else float("nan")

        rows.append({
            "R1_filename": r1_name,
            "R2_filename": r2_name,
            "fragments_scanned": total_pairs,
            "HALF_LEADER_matches": half_leader_count,
            "JUNCTION_matches": junction_count,
            "HALF_LEADER_CPM": round(half_leader_cpm, 3),
            "JUNCTION_CPM": round(junction_cpm, 3),
            "JUNCTION_to_HALF_LEADER_ratio": round(ratio, 4) if half_leader_count else "NA",
        })

        print(
            f"  fragments: {total_pairs} | HALF_LEADER: {half_leader_count} (CPM {half_leader_cpm:.3f}) "
            f"| JUNCTION: {junction_count} (CPM {junction_cpm:.3f}) "
            f"| JUNCTION:HALF_LEADER ratio: {rows[-1]['JUNCTION_to_HALF_LEADER_ratio']}\n"
        )

    if not rows:
        print("No valid R1/R2 pairs were processed.")
        return

    csv_path = os.path.join(FASTQ_FOLDER, OUTPUT_CSV_PATH)
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Done. Summary written to:\n  {csv_path}")


if __name__ == "__main__":
    main()
