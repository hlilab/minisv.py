"""Main interface for breakpoints caller"""

import argparse

from gaftools.io import load_gaf_to_grouped_reads


# function to return final node in graph path as this is where the BND is occuring
def get_contig(s, ch, location, node_end):
    # ignore if not in graph path
    # NOTE: hg19 do not have chr as start string
    #      we extended to hg19 with hard-coded conditions
    # print(s)
    if (
        s.startswith("chr")
        or s.isdigit()
        or s.startswith("GL")
        or s.startswith("NC")
        or s.startswith("hs")
        or s in ["MT", "X", "Y"]
    ):
        return s, str(location)
    else:
        nodes = [i for i, ltr in enumerate(s) if ltr in ch]
        if len(nodes) == 1:
            return s, str(location)
        # start of breakpoint, want end coord of read mapping, so last node
        if node_end:
            total = 0
            for i in range(1, len(nodes)):
                cur = s[nodes[i - 1] : nodes[i]]
                locs = cur.split(":")[-1].split("-")
                diff = int(locs[1]) - int(locs[0])
                total += abs(diff)
            return s[nodes[-1] :], str(location - total)
        # end of break want start coord, so first node
        else:
            return s[nodes[0] : nodes[1]], str(location)


def find_break(m1, m2):
    """
    Function to find the break point between two chunks of the read mapping
    Args:
    - m1 (list): First read mapping information.
    - m2 (list): Second read mapping information.
    Returns:
    - list: Breakpoint information.

    NOTE: m1 and m2 should share the same read name

    """
    assert m1[0] == m2[0], "Read names do not match"
    # translate strands into breakpoint notation
    strand_translation = {"+": ">", "-": "<"}
    # determine break locations based on strand
    b1 = 0
    # m1[5]: genome graph path
    if m1[4] == "+":
        # https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
        # m1[8]: end pos
        b1 = int(m1[8])
        contig_m1 = get_contig(m1[5], [">", "<"], b1, True)
    elif m1[4] == "-":
        # m1[8]: start pos
        b1 = int(m1[7])
        contig_m1 = get_contig(m1[5], [">", "<"], b1, False)
    b2 = 0
    if m2[4] == "+":
        b2 = int(m2[7])
        contig_m2 = get_contig(m2[5], [">", "<"], b2, False)
    elif m2[4] == "-":
        b2 = int(m2[8])
        contig_m2 = get_contig(m2[5], [">", "<"], b2, True)

    # return a more legible break point text format
    return [
        contig_m1[0],
        contig_m1[1],
        strand_translation[m1[4]] + strand_translation[m2[4]],
        contig_m2[0],
        contig_m2[1],
        str(min([int(m1[11]), int(m2[11])])),
        m1[0],
    ]


# function to convert node locations into chromosome/contig locations
def convert(contig, brk, direct):
    if not contig.startswith(("<", ">")):
        return contig, brk, direct
    locs = contig.split(":")[-1].split("-")
    # print(contig)
    if contig.startswith(">"):
        # print( contig.split(":")[0][1:], str(int(locs[0]) + int(brk)), contig.split(":")[0][0])
        return (
            contig.split(":")[0][1:],
            str(int(locs[0]) + int(brk)),
            contig.split(":")[0][0],
        )
    else:
        return (
            contig.split(":")[0][1:],
            str(int(locs[1]) - int(brk)),
            contig.split(":")[0][0],
        )


def adj_graph_paths(original):
    # swap "<<" to ">>" and "<>" to "><" based on symmetries
    #    print(original)
    first = convert(original[0], original[1], original[2].strip()[0])
    second = convert(original[3], original[4], original[2].strip()[1])
    arrows = first[2] + second[2]
    #    print(first)
    #   print(second)
    #  print(arrows)
    if arrows == "<<":
        temp = first
        first = second
        second = temp
        arrows = ">>"
    if arrows == "<>" or arrows == "><":
        if first[0] > second[0]:
            temp = first
            first = second
            second = temp

    original[0] = first[0]
    original[1] = first[1]
    original[2] = arrows
    original[3] = second[0]
    original[4] = second[1]
    return original


def call_breakpoints(read_cluster, min_mapQ=5, min_map_len=2000):
    """Main function to call breakpoints from a cluster of reads
    Args:
    - read_cluster (list): List of reads in a cluster.
    - min_mapQ (int): Minimum mapping quality threshold (default: 5).
    - min_map_len (int): Minimum mapping length threshold (default: 2000).
    Returns:
    - output (list): List of breakpoints in the cluster.
    """
    output = []
    for i in range(1, len(read_cluster)):
        breaks = find_break(read_cluster[i - 1], read_cluster[i])
        if breaks:
            output.append(adj_graph_paths(breaks))
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify Break Points from GAF input")
    parser.add_argument(
        "-i", metavar="<input.gaf>", required=True, help="input GAF file"
    )
    parser.add_argument(
        "-m",
        metavar="--min_mapping_quality",
        required=False,
        type=int,
        default=5,
        help="Minimum mapping quality (integer)",
    )
    parser.add_argument(
        "-a",
        metavar="--min_alignment_length",
        required=False,
        type=int,
        default=2000,
        help="Minimum alignment length (bps)",
    )

    args = parser.parse_args()
    gafFile = args.i
    min_mapQ = args.m
    min_map_len = args.a

    groups = load_gaf_to_grouped_reads(args.i, min_mapQ, min_map_len)
    all_breaks = []
    for group in groups:
        if len(group) > 1:
            brks = call_breakpoints(group, min_mapQ, min_map_len)
            for brk in brks:
                all_breaks.append("\t".join(brk))
    print("\n".join(all_breaks))
