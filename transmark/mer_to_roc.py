#!/usr/bin/env python3
"""
Convert a .mer file to a 3-column ROC file.

Output format:
    pct_found  pct_false  tool

One data point is emitted per false positive hit encountered.
pct_found = cumulative true positives / total positives in benchmark
pct_false = cumulative false positives / total negatives in database

Usage:
    python mer_to_roc.py <merfile> <tool_name> <total_negatives> [--no-header]

Arguments:
    merfile          Path to the .mer file
    tool_name        Label to use in the tool column (e.g. "bathsearch")
    total_negatives  Total number of negative (decoy) sequences in the
                     search database -- used as the denominator for pct_false
    --no-header      Omit the "pct_found pct_false tool" header line

The total number of positives is read from the "*summary*" line in the
.mer file. If no summary line is found, the script falls back to counting
all "+" hits in the file.
"""

import sys
import argparse


def parse_total_pos(lines):
    """Extract total positives from the *summary* line."""
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("*summary*"):
            fields = stripped.split()
            # format: *summary* <mer> <npos> <fn> <fp> <thresh> <time>
            if len(fields) >= 3:
                return int(fields[2])
    return None


def mer_to_roc(merfile, tool_name, total_neg, write_header):
    with open(merfile) as fh:
        lines = fh.readlines()

    total_pos = parse_total_pos(lines)

    # Fallback: count all "+" data lines
    if total_pos is None:
        total_pos = sum(
            1 for line in lines
            if line.startswith("=") and " + " in line
        )
        sys.stderr.write(
            f"Warning: no *summary* line found; counted {total_pos} positives from '+' lines.\n"
        )

    if total_pos == 0:
        sys.stderr.write("Error: could not determine total positives.\n")
        sys.exit(1)

    num_pos = 0
    num_neg = 0

    if write_header:
        print("pct_found pct_false tool")

    for line in lines:
        if not line.startswith("="):
            continue
        fields = line.split()
        if len(fields) < 6:
            continue

        label = fields[5]

        if label == "+":
            num_pos += 1
        elif label == "-":
            num_neg += 1
            pct_found = num_pos / total_pos
            pct_false = num_neg / total_neg
            print(f"{pct_found:.6f} {pct_false:.6f} {tool_name}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a .mer file to a pct_found/pct_false ROC file."
    )
    parser.add_argument("merfile", help="Input .mer file")
    parser.add_argument("tool_name", help="Tool label for the third column")
    parser.add_argument(
        "total_negatives",
        type=float,
        help="Total number of negative sequences in the search database",
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Omit the header line from the output",
    )
    args = parser.parse_args()

    mer_to_roc(
        args.merfile,
        args.tool_name,
        args.total_negatives,
        write_header=not args.no_header,
    )


if __name__ == "__main__":
    main()
