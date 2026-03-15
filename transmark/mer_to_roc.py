#!/usr/bin/env python3
"""
Convert a .mer file (rmark-mer.pl output) to a 3-column ROC file.

Output format:
    pct_found  pct_false  tool

One data point is emitted per false positive hit encountered.
  pct_found = (TPs with E-value strictly less than this FP's E-value) / total_pos
  pct_false = (cumulative FPs before this point) / num_families

Hits are grouped by E-value.  Within a group, FPs are emitted first using
the TP count accumulated before that group; TPs in the group are added to
the running total only after all FPs in the group have been emitted.  This
ensures a TP with the same E-value as a FP is not counted as "before" it.

The first row has pct_false = 0 (TPs found before the first FP).
The last row has pct_false = 1.0 (one average FP per family).
num_families is derived from the unique family names in the "=" lines.

Usage:
    python mer_to_roc.py <merfile> <tool_name> [--no-header]

Arguments:
    merfile     Path to the .mer file
    tool_name   Label for the third column (e.g. "bathsearch")
    --no-header Omit the "pct_found pct_false tool" header line
"""

import sys
import argparse
from itertools import groupby


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


def mer_to_roc(merfile, tool_name, write_header):
    with open(merfile) as fh:
        lines = fh.readlines()

    total_pos = parse_total_pos(lines)

    # Collect data lines and count unique families
    data_lines = []
    families = set()
    for line in lines:
        if not line.startswith("="):
            continue
        fields = line.split()
        if len(fields) < 6:
            continue
        data_lines.append(fields)
        families.add(fields[2])   # family column

    num_families = len(families)

    # Fallback: count "+" lines for total_pos
    if total_pos is None:
        total_pos = sum(1 for f in data_lines if f[5] == "+")
        sys.stderr.write(
            f"Warning: no *summary* line found; "
            f"counted {total_pos} positives from '+' lines.\n"
        )

    if total_pos == 0:
        sys.stderr.write("Error: could not determine total positives.\n")
        sys.exit(1)

    if num_families == 0:
        sys.stderr.write("Error: no '=' data lines found in .mer file.\n")
        sys.exit(1)

    if write_header:
        print("pct_found pct_false tool")

    num_pos = 0
    num_neg = 0
    done = False

    # Group consecutive hits by E-value (field[4]).
    # Within each group, emit FPs first (using TPs accumulated before this
    # group), then add the group's TPs to the running total.
    for _evalue, group in groupby(data_lines, key=lambda f: f[4]):
        if done:
            break
        tp_in_group = 0
        fp_in_group = 0

        for fields in group:
            if fields[5] == "+":
                tp_in_group += 1
            elif fields[5] == "-":
                fp_in_group += 1

        # Emit one data point per FP, using num_pos from before this group
        for _ in range(fp_in_group):
            pct_found = num_pos / total_pos
            pct_false = num_neg / num_families
            print(f"{pct_found:.6f} {pct_false:.6f} {tool_name}")
            num_neg += 1
            if pct_false >= 1.0:
                done = True
                break

        # Now credit TPs from this group
        num_pos += tp_in_group


def main():
    parser = argparse.ArgumentParser(
        description="Convert a .mer file to a pct_found/pct_false ROC file."
    )
    parser.add_argument("merfile", help="Input .mer file")
    parser.add_argument("tool_name", help="Tool label for the third column")
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Omit the header line from the output",
    )
    args = parser.parse_args()

    mer_to_roc(args.merfile, args.tool_name, write_header=not args.no_header)


if __name__ == "__main__":
    main()
