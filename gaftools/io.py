""" A generalized input interface for both INDELs and breakpoints
We could put the output interface into this submodule as well.
"""
import argparse
import gzip
import sys
from dataclasses import dataclass

from gaftools.identify_breaks_v5 import call_breakpoints
from gaftools.merge_break_pts_v3 import merge_breaks


@dataclass
class GroupedResults:
    indels: str
    breakpoints: str


def load_gaf_to_grouped_reads(gafFile, min_mapQ=5, min_map_len=2000):
    """
    Load GAF file, filter based on mapping quality and minimum map length,
    and sort the lines by cluster location in read.

    Args:
    - gafFile (str): Path to the GAF file.
    - min_mapQ (int): Minimum mapping quality threshold (default: 5).
    - min_map_len (int): Minimum mapping length threshold (default: 2000).

    Returns:
    - sorted_lines (list): List of GAF lines sorted by cluster location in read.
    """
    # Read lines, split into fields, and filter based on conditions
    lines = []
    with open(gafFile, "r") as gafFileHandler:
        for line in gafFileHandler:
            fields = line.strip().split("\t")
            read_name = fields[0]

            # Skip low mapping quality, short aligned length
            # For minimap2 paf, we skipped supplementary alignment
            if (
                fields[16] == "tp:A:S"
                or int(fields[11]) < min_mapQ
                or int(fields[8]) - int(fields[7]) < min_map_len
            ):
                continue

            if len(lines) >= 1:
                if read_name == lines[-1][0]:
                    lines.append(fields)
                else:
                    # locally sort a group of reads
                    sorted_lines = sorted(lines, key=lambda x: (int(x[2])))
                    yield sorted_lines
                    lines = [fields]
            else:
                lines.append(fields)
    sorted_lines = sorted(lines, key=lambda x: (int(x[2])))
    yield sorted_lines
    del lines
    del sorted_lines


def merge_indel_breakpoints(prefix, break_point_file, indel_file):
    """
    Merge indel and breakpoint files into a single file.

    Args:
    - break_point_file (str): Path to the breakpoint file.
    - indel_file (str): Path to the indel file.

    Returns:
    - merged_file (str): Path to the merged file.
    """
    col_names = [
        "#chrom_contig1",
        "start1",
        "end1",
        "chrom_contig2",
        "start2",
        "end2",
        "indel_length",
        "tsd_length",
        "polyA_length",
        "indel_seq",
        "read_counts",
        "dot",
        "average_mapquality",
        "forward_strand_counts",
        "reverse_strand_counts",
        "VNTR_hits",
        "Centromere_hits",
        "L1_hits",
        "read_names",
        "sv_type",
    ]
    breakpoints_indexes = [
        0,
        1,
        -1,
        3,
        4,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        5,
        -1,
        -1,
        -1,
        -1,
        -1,
        6,
    ]
    indel_indexes = [0, 1, 2, -1, -1, -1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

    harmonized_output = gzip.open(f"{prefix}_harmonized.bed.gz", "wt")
    harmonized_output.write("\t".join(col_names) + "\n")

    with gzip.open(break_point_file, "rt") as break_point_file_handler:
        for line in break_point_file_handler:
            line = line.strip().split("\t")
            line = [line[i] if i != -1 else "." for i in breakpoints_indexes]
            out_str = "\t".join(line)
            harmonized_output.write(f"{out_str}\tbnd\n")

    with gzip.open(indel_file, "rt") as indel_file_handler:
        for line in indel_file_handler:
            line = line.strip().split("\t")
            line = [line[i] if i != -1 else "." for i in indel_indexes]
            out_str = "\t".join(line)
            harmonized_output.write(f"{out_str}\tindel\n")

    harmonized_output.close()


def merge_sv_vcf(prefix, indel_vcf, breakpoint_vcf):
    """Merge indel and breakpoint VCF files into a single file.
    Args:
    """
    return


if __name__ == "__main__":
    # test script locally
    parser = argparse.ArgumentParser(description="Identify Break Points from GAF input")
    parser.add_argument("-i", required=True, help="input GAF file")
    args = parser.parse_args()
    groups = load_gaf_to_grouped_reads(args.i)
    all_breaks = []
    for group in groups:
        if len(group) > 1:
            brks = call_breakpoints(group)
            for brk in brks:

              all_breaks.append(brk)  

    bnd_bed, bnd_vcf = merge_breaks(lines, 100, 5)
    #sys.stdout.write('\n'.join(bnd_bed))
    sys.stdout.write('\n'.join(bnd_vcf))
   
        








