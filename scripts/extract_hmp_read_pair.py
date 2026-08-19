#!/usr/bin/env python3
"""Extract an HMP paired FASTQ from a sample tarball, optionally subsetting reads."""

from __future__ import annotations

import argparse
import tarfile
from pathlib import Path


def copy_fastq(member_file, output: Path, read_pair_limit: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    line_limit = read_pair_limit * 4 if read_pair_limit > 0 else None
    with output.open("wb") as out:
        for i, line in enumerate(member_file):
            if line_limit is not None and i >= line_limit:
                break
            out.write(line)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tarball", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--read-pair-limit", type=int, default=0)
    args = parser.parse_args()

    targets = {
        f"{args.sample}/{args.sample}.denovo_duplicates_marked.trimmed.1.fastq": args.r1,
        f"{args.sample}/{args.sample}.denovo_duplicates_marked.trimmed.2.fastq": args.r2,
    }
    found: set[str] = set()
    with tarfile.open(args.tarball, "r:bz2") as tar:
        for member in tar:
            output = targets.get(member.name)
            if output is None:
                continue
            handle = tar.extractfile(member)
            if handle is None:
                raise SystemExit(f"Could not extract {member.name} from {args.tarball}")
            with handle:
                copy_fastq(handle, output, args.read_pair_limit)
            found.add(member.name)
            if len(found) == len(targets):
                break

    missing = sorted(set(targets) - found)
    if missing:
        raise SystemExit(f"Missing FASTQ member(s) in {args.tarball}: {', '.join(missing)}")


if __name__ == "__main__":
    main()
