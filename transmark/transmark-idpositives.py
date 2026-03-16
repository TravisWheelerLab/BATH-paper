#!/usr/bin/env python3
#
# Given a positive file (.pos) and an output file of a transmark benchmark,
# first remove all overlapping hits from the file, as defined below.
# Then given the set of non-overlapping hits, create a positive-annotated
# output file that lists the non-overlapping hits in the benchmark output file,
# along with two extra fields indicating if each hit overlaps a positive or not
# and whether the overlap is on the correct strand.
#
# The <transmark outfile> MUST be sorted properly by score (E-value or bit score),
# with better scores preceding worse scores.
#
# Usage:   python transmark-idpositives.py <posfile> <sorted transmark output>
# Example: ./transmark-idpositives.py transmark-00.pos sorted-cmsearch.out
#
# Faster rewrite of transmark-idpositives.pl:
#   - Reads the outfile once instead of three times
#   - Groups kept hits by (query, target_name, strand) so overlap removal is
#     O(hits_per_group) instead of O(total_kept_hits) — eliminates the O(N^2) loop
#   - Uses binary search (bisect) in CheckIfPositive to skip positives whose
#     start position is already past the hit's end
#
import sys
import bisect

OVERLAP_THR = 0.5


def get_overlap(from1, to1, from2, to2):
    """Return fractional overlap of two intervals (each already normalized so from <= to)."""
    minlen = min(to1 - from1 + 1, to2 - from2 + 1)
    if from1 > from2:
        from1, from2 = from2, from1
        to1, to2 = to2, to1
    # from1 <= from2 now
    if to1 < from2:
        return 0.0           # case 1: no overlap
    if to1 < to2:
        return (to1 - from2 + 1) / minlen   # case 2
    return (to2 - from2 + 1) / minlen       # case 3


def check_if_positive(target_from, target_to, target_ori, overlap_thr,
                      pos_fam, pos_to_dict, pos_ori, pos_idx, pos_order):
    """
    Check whether a hit [target_from, target_to] (already normalized) overlaps
    any positive in this target sequence.

    pos_order is a sorted list of pos_from values.  We use binary search to skip
    all positives whose start position is already beyond target_to (they cannot
    overlap), then scan the remainder and skip those whose end is before target_from.

    Returns (fam, idx, strand).
    """
    fam = "decoy"
    strand = "decoy"
    idx = 0

    # Binary search: find index of first pos_from > target_to.
    # All positives at indices < right_idx have pos_from <= target_to and are candidates.
    right_idx = bisect.bisect_right(pos_order, target_to)

    for i in range(right_idx):
        pos_from = pos_order[i]
        pt = pos_to_dict[pos_from]
        if pt < target_from:
            continue   # this positive ends before our hit starts — no overlap
        overlap = get_overlap(target_from, target_to, pos_from, pt)
        if overlap > overlap_thr:
            fam = pos_fam[pos_from]
            strand = "same" if pos_ori[pos_from] == target_ori else "opposite"
            idx = pos_idx[pos_from]
            # keep going in case another positive also overlaps (matches original behavior)

    return fam, idx, strand


def main():
    if len(sys.argv) != 3:
        print("Incorrect number of command line arguments.")
        print("Usage: python transmark-idpositives.py <posfile> <transmark outfile>")
        sys.exit(1)

    posfile, outfile = sys.argv[1], sys.argv[2]

    import os
    if not os.path.exists(posfile):
        sys.exit(f"{posfile} doesn't exist")
    if not os.path.exists(outfile):
        sys.exit(f"{outfile} doesn't exist")

    # ------------------------------------------------------------------
    # Step 1: Read the positives file
    # ------------------------------------------------------------------
    pos_fam_HH = {}   # target_name -> {pos_from -> family}
    pos_to_HH  = {}   # target_name -> {pos_from -> pos_to}
    pos_ori_HH = {}   # target_name -> {pos_from -> ori}
    pos_idx_HH = {}   # target_name -> {pos_from -> idx}

    with open(posfile) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#'):
                continue
            fields = line.split(None, 4)
            family_and_idx = fields[0]
            target_name    = fields[1]
            target_from    = int(fields[3])
            target_to      = int(fields[4])

            slash = family_and_idx.rfind('/')
            if slash < 0:
                sys.exit(f"ERROR invalid family and idx field {family_and_idx}")
            idx    = family_and_idx[slash + 1:]
            family = family_and_idx[:slash]
            if '/' in family:
                sys.exit(f"ERROR family {family} from {family_and_idx}, contains two '/', it should only have 1!")
            if target_name == "decoy":
                sys.exit('ERROR a family named "decoy" exists in the benchmark dataset, this is not allowed.')

            if target_from > target_to:
                target_from, target_to = target_to, target_from
                target_ori = "-"
            else:
                target_ori = "+"

            if target_name not in pos_fam_HH:
                pos_fam_HH[target_name] = {}
                pos_to_HH[target_name]  = {}
                pos_ori_HH[target_name] = {}
                pos_idx_HH[target_name] = {}

            pos_fam_HH[target_name][target_from] = family
            pos_to_HH[target_name][target_from]  = target_to
            pos_ori_HH[target_name][target_from] = target_ori
            pos_idx_HH[target_name][target_from] = idx

    # Sorted arrays of start points per target (for binary search)
    pos_order_HA = {t: sorted(pos_to_HH[t].keys()) for t in pos_to_HH}

    # ------------------------------------------------------------------
    # Step 2: Read the outfile once into memory and verify sorting
    # ------------------------------------------------------------------
    out_lines = []
    with open(outfile) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#'):
                continue
            out_lines.append(line)

    # Check score sort order
    seen_sc = False
    sc_should_increase = False
    sc_should_decrease = False
    prv_sc = None
    for line in out_lines:
        sc = float(line.split(None, 1)[0])
        if seen_sc:
            if sc > prv_sc:
                if sc_should_decrease:
                    sys.exit("ERROR, results don't appear to be sorted by score")
                sc_should_increase = True
            elif sc < prv_sc:
                if sc_should_increase:
                    sys.exit("ERROR, results don't appear to be sorted by score")
                sc_should_decrease = True
        prv_sc = sc
        seen_sc = True

    # ------------------------------------------------------------------
    # Step 3: Greedy overlap removal.
    #
    # Key optimization over the Perl version: instead of checking every
    # new hit against ALL previously kept hits (O(N^2)), we group kept
    # hits by (query, target_name, strand).  A new hit only needs to be
    # compared against the (usually small) set of hits already kept in
    # the same group.
    # ------------------------------------------------------------------
    kept_groups = {}   # (query, target_name, ori) -> list of (from, to)
    hit_useme   = []   # parallel to out_lines

    for line in out_lines:
        fields       = line.split(None, 5)
        target_from  = int(fields[2])
        target_to    = int(fields[3])
        target_name  = fields[4]
        query        = fields[5]

        if target_from > target_to:
            target_from, target_to = target_to, target_from
            target_ori = "-"
        else:
            target_ori = "+"

        key       = (query, target_name, target_ori)
        intervals = kept_groups.get(key)

        use_hit = True
        if intervals:
            for ifrom, ito in intervals:
                if get_overlap(target_from, target_to, ifrom, ito) > OVERLAP_THR:
                    use_hit = False
                    break

        hit_useme.append(use_hit)

        if use_hit:
            if intervals is None:
                kept_groups[key] = []
            kept_groups[key].append((target_from, target_to))

    # ------------------------------------------------------------------
    # Step 4: Annotate each kept hit with positive/decoy status
    # ------------------------------------------------------------------
    matching_idx = 0   # matches Perl's uninitialized-but-numeric-context behavior

    pout_lines = []
    for h, line in enumerate(out_lines):
        if not hit_useme[h]:
            continue
        fields       = line.split(None, 5)
        target_from  = int(fields[2])
        target_to    = int(fields[3])
        target_name  = fields[4]

        if target_from > target_to:
            target_from, target_to = target_to, target_from
            target_ori = "-"
        else:
            target_ori = "+"

        if target_name not in pos_fam_HH:
            matching_fam    = "decoy"
            matching_strand = "decoy"
            # matching_idx retains previous value (matches Perl behavior)
        else:
            matching_fam, matching_idx, matching_strand = check_if_positive(
                target_from, target_to, target_ori, OVERLAP_THR,
                pos_fam_HH[target_name],
                pos_to_HH[target_name],
                pos_ori_HH[target_name],
                pos_idx_HH[target_name],
                pos_order_HA[target_name],
            )

        pout_lines.append(f"{line} {matching_fam}/{matching_idx} {matching_strand}\n")

    sys.stdout.writelines(pout_lines)


if __name__ == "__main__":
    main()
