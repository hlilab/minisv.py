""" A generalized input interface for both INDELs and breakpoints
We could put the output interface into this submodule as well.
"""
# import argparse
import gzip
import re
from dataclasses import dataclass

from .regex import re_info

# from .identify_breaks_v5 import call_breakpoints
# from gaftools.merge_break_pts_v3 import merge_breaks


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

    NOTE: This function is a generator. The PAF/GAF file already outputs sorted reads by names.
    """
    lines = []
    with open(gafFile, "r") as gafFileHandler:
        for line in gafFileHandler:
            fields = line.strip().split("\t")
            read_name = fields[0]

            # Skip low mapping quality, short aligned length
            # For minimap2 paf, we skipped secondary alignment marked by tp:A:S
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
    We use generic breakpoints as the output for both indel and breakpoints

    Args:
    - break_point_file (str): Path to the breakpoint file.
    - indel_file (str): Path to the indel file.

    Returns:
    - merged_file path (str): Path to the merged bedpe format file.
    """
    col_names = [
        "#chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",  # 5
        "name",  # 11
        "score",
        "strand1",
        "strand2",
        "indel_length",
        "tsd_length",
        "polyA_length",
        "indel_seq",
        "read_counts",
        "average_mapquality",
        "forward_strand_counts",
        "reverse_strand_counts",
        "VNTR_hits",
        "Centromere_hits",
        "L1_hits",
        "read_names",
        "strandness",
        "Centromere_dist",
        "sv_type",
    ]
    breakpoints_indexes = [
        0,
        1,
        -1,
        3,
        4,
        -1,  # 5
        -1,  # 11
        -1,  # score
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
        2,
        -1,
    ]
    indel_indexes = [
        0,
        1,
        2,
        -1,
        -1,
        -1,  # 5
        8,  # 11
        -1,  # score
        -1,
        -1,
        3,
        4,
        5,
        6,
        7,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        -1,
        16,
    ]

    harmonized_output = gzip.open(f"{prefix}_harmonized.bed.gz", "wt")
    harmonized_output.write("\t".join(col_names) + "\n")

    with gzip.open(break_point_file, "rt") as break_point_file_handler:
        for line in break_point_file_handler:
            line = line.strip().split("\t")
            # centromere/vntr annotation
            # print(line)
            # print(line[-1])
            cent_hit, cent_hit2, vntr_hit, vntr_hit2, cent_dist = line[-1].split(",")
            line = [line[i] if i != -1 else "." for i in breakpoints_indexes]
            line[15] = f"{vntr_hit},{vntr_hit2}"
            line[16] = f"{cent_hit},{cent_hit2}"
            line[-1] = str(cent_dist)
            out_str = "\t".join(line)
            harmonized_output.write(f"{out_str}\tbnd\n")

    with gzip.open(indel_file, "rt") as indel_file_handler:
        for line in indel_file_handler:
            line = line.strip().split("\t")
            line = [line[i] if i != -1 else "." for i in indel_indexes]
            out_str = "\t".join(line)
            harmonized_output.write(f"{out_str}\tindel\n")

    harmonized_output.close()


def merge_breakpoint_indel_vcf(prefix, indel_vcf, breakpoint_vcf):
    """Merge indel and breakpoint vcf files into a single file.

    Args:
    - prefix (str): Prefix for the output file.
    - indel_vcf (str): Path to the indel vcf file.
    - breakpoint_vcf (str): Path to the breakpoint vcf file.

    Returns:
    - merged_vcf (str): Path to the merged vcf file.
    """
    merged_vcf = f"{prefix}_merged.vcf"
    merged_vcf_handler = open(merged_vcf, "w")

    with open(indel_vcf, "r") as indel_vcf_handler:
        for line in indel_vcf_handler:
            merged_vcf_handler.write(line)

    with open(breakpoint_vcf, "r") as breakpoint_vcf_handler:
        for line in breakpoint_vcf_handler:
            if line.startswith("#"):
                continue
            merged_vcf_handler.write(line)

    merged_vcf_handler.close()
    return merged_vcf


def write_vcf(opt, input):
    print("##fileformat=VCFv4.3")
    print(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
    )
    print(
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
    )
    print(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'
    )
    print('##ALT=<ID=DEL,Description="Deletion">')
    print('##ALT=<ID=INS,Description="Insertion">')
    print('##ALT=<ID=DUP,Description="Duplication">')
    print('##ALT=<ID=INV,Description="Inversion">')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample")
    key = {"SVTYPE": 1, "SVLEN": 1}
    for line in input:
        t = line.strip().split()
        is_bp = re.match(r"[><]", t[2])
        off_info = 8 if is_bp else 6
        type = None
        info = ""
        for m in re_info.findall(t[off_info]):
            if m[0] in key:
                if len(info) > 0:
                    info += ";"
                info += f"{m[0]}={m[1]}"
            if m[0] == "SVTYPE":
                type = m[1]
        if type is None or type == "BND":
            continue
        info += f";END={t[4]}" if is_bp else f";END={t[2]}"
        print(
            t[0],
            t[1],
            ".",
            "N",
            f"<{type}>",
            t[off_info - 2],
            ".",
            info,
            "GT",
            "1/1",
            sep="\t",
        )


# if __name__ == "__main__":
#    # test script locally
#    parser = argparse.ArgumentParser(description="Identify Break Points from GAF input")
#    parser.add_argument("-i", required=True, help="input GAF file")
#    args = parser.parse_args()
#    groups = load_gaf_to_grouped_reads(args.i)
#    all_breaks = []
#    for group in groups:
#        if len(group) > 1:
#            brks = call_breakpoints(group)
#            for brk in brks:
#                all_breaks.append(brk)
#
#    # bnd_bed, bnd_vcf = merge_breaks(lines, 100, 5)
#    # sys.stdout.write('\n'.join(bnd_bed))
#    # sys.stdout.write("\n".join(bnd_vcf))
