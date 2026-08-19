#!/usr/bin/env python3
"""Repair paired FASTQ files by dropping malformed read pairs."""

from __future__ import annotations

import argparse
from pathlib import Path


VALID_BASES = set("ACGTNacgtn.")
MAX_REASON_EXAMPLES = 20


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


def record_id(header: str) -> str:
    text = header.strip()
    if text.startswith("@"):
        text = text[1:]
    text = text.split()[0]
    if text.endswith("/1") or text.endswith("/2"):
        text = text[:-2]
    return text


def validate_record(record, mate: str) -> list[str]:
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


def write_record(handle, record) -> None:
    header, seq, plus, qual, _status = record
    handle.write(header)
    handle.write(seq)
    handle.write(plus)
    handle.write(qual)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--out-r1", required=True, type=Path)
    parser.add_argument("--out-r2", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    args.out_r1.parent.mkdir(parents=True, exist_ok=True)
    args.out_r2.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    total = kept = dropped = 0
    examples: list[str] = []
    with args.r1.open(errors="replace") as r1, args.r2.open(errors="replace") as r2, \
            args.out_r1.open("w") as out1, args.out_r2.open("w") as out2:
        while True:
            rec1 = read_record(r1)
            rec2 = read_record(r2)
            if rec1 is None and rec2 is None:
                break
            total += 1
            reasons: list[str] = []
            if rec1 is None:
                reasons.append("r1:missing_record")
            if rec2 is None:
                reasons.append("r2:missing_record")
            if rec1 is not None:
                reasons.extend(validate_record(rec1, "r1"))
            if rec2 is not None:
                reasons.extend(validate_record(rec2, "r2"))
            if rec1 is not None and rec2 is not None and not reasons:
                id1 = record_id(rec1[0])
                id2 = record_id(rec2[0])
                if id1 != id2:
                    reasons.append(f"pair:id_mismatch:{id1}:{id2}")
            if reasons:
                dropped += 1
                if len(examples) < MAX_REASON_EXAMPLES:
                    read_name = ""
                    if rec1 is not None:
                        read_name = rec1[0].strip()
                    elif rec2 is not None:
                        read_name = rec2[0].strip()
                    examples.append(f"{total}\t{read_name}\t{';'.join(reasons)}")
                continue
            kept += 1
            write_record(out1, rec1)
            write_record(out2, rec2)

    with args.report.open("w") as report:
        report.write("metric\tvalue\n")
        report.write(f"total_pairs\t{total}\n")
        report.write(f"kept_pairs\t{kept}\n")
        report.write(f"dropped_pairs\t{dropped}\n")
        report.write("\nexample_pair_index\tread_name\treasons\n")
        for item in examples:
            report.write(item + "\n")


if __name__ == "__main__":
    main()
