#!/usr/bin/env python3
"""Validate paired FASTQ files for HMP read-processing handoffs."""

from __future__ import annotations

import argparse
from pathlib import Path


VALID_BASES = set("ACGTNacgtn.")


def read_record(handle):
    header = handle.readline()
    if header == "":
        return None
    seq = handle.readline()
    plus = handle.readline()
    qual = handle.readline()
    if seq == "" or plus == "" or qual == "":
        return (header, seq, plus, qual, "truncated_record")
    return (header, seq, plus, qual, "")


def clean_id(header: str) -> str:
    text = header.strip()
    if text.startswith("@"):
        text = text[1:]
    text = text.split()[0]
    if text.endswith("/1") or text.endswith("/2"):
        text = text[:-2]
    return text


def validate_record(record, mate: str) -> list[str]:
    if record is None:
        return [f"{mate}:missing"]
    header, seq, plus, qual, status = record
    reasons: list[str] = []
    if status:
        reasons.append(f"{mate}:{status}")
        return reasons
    if not header.startswith("@"):
        reasons.append(f"{mate}:header_not_at")
    if not plus.startswith("+"):
        reasons.append(f"{mate}:plus_not_plus")
    seq_text = seq.rstrip("\n\r")
    qual_text = qual.rstrip("\n\r")
    if len(seq_text) != len(qual_text):
        reasons.append(f"{mate}:seq_qual_len:{len(seq_text)}:{len(qual_text)}")
    if "\x00" in seq_text or "\x00" in qual_text:
        reasons.append(f"{mate}:nul_byte")
    if not set(seq_text) <= VALID_BASES:
        reasons.append(f"{mate}:invalid_base")
    return reasons


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--max-examples", type=int, default=20)
    parser.add_argument("--expected-pairs", type=int)
    parser.add_argument("--required-read-id", action="append", default=[])
    parser.add_argument("--track-read-id", action="append", default=[])
    args = parser.parse_args()

    args.report.parent.mkdir(parents=True, exist_ok=True)
    total_pairs = 0
    error_pairs = 0
    examples: list[str] = []
    required_ids = {clean_id(item): False for item in args.required_read_id}
    tracked_ids = {clean_id(item): False for item in args.track_read_id}

    with args.r1.open(errors="replace") as r1, args.r2.open(errors="replace") as r2:
        while True:
            rec1 = read_record(r1)
            rec2 = read_record(r2)
            if rec1 is None and rec2 is None:
                break
            total_pairs += 1
            reasons = validate_record(rec1, "r1")
            reasons.extend(validate_record(rec2, "r2"))
            if rec1 is not None and rec2 is not None and not reasons:
                id1 = clean_id(rec1[0])
                id2 = clean_id(rec2[0])
                if id1 != id2:
                    reasons.append(f"pair:id_mismatch:{id1}:{id2}")
                if id1 in required_ids:
                    required_ids[id1] = True
                if id1 in tracked_ids:
                    tracked_ids[id1] = True
            if reasons:
                error_pairs += 1
                if len(examples) < args.max_examples:
                    name = rec1[0].strip() if rec1 is not None else ""
                    examples.append(f"{total_pairs}\t{name}\t{';'.join(reasons)}")

    with args.report.open("w") as report:
        report.write("metric\tvalue\n")
        report.write(f"total_pairs\t{total_pairs}\n")
        report.write(f"error_pairs\t{error_pairs}\n")
        if args.expected_pairs is not None:
            report.write(f"expected_pairs\t{args.expected_pairs}\n")
            report.write(f"pair_count_matches_expected\t{int(total_pairs == args.expected_pairs)}\n")
        for read_id, found in required_ids.items():
            report.write(f"required_read_found:{read_id}\t{int(found)}\n")
        for read_id, found in tracked_ids.items():
            report.write(f"tracked_read_found:{read_id}\t{int(found)}\n")
        report.write("\nexample_pair_index\tread_name\treasons\n")
        for example in examples:
            report.write(example + "\n")

    failures = []
    if error_pairs:
        failures.append(f"{error_pairs} error pairs")
    if args.expected_pairs is not None and total_pairs != args.expected_pairs:
        failures.append(f"expected {args.expected_pairs} pairs, observed {total_pairs}")
    missing_ids = [read_id for read_id, found in required_ids.items() if not found]
    if missing_ids:
        failures.append("missing required reads: " + ",".join(missing_ids))
    if failures:
        raise SystemExit("FASTQ validation failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
