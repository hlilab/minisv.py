"""An module for parsing GAF/PAF format from minigraph alignment and minimap2 alignment
"""

import functools
import gzip
import math
import os
import re
from datetime import datetime
from multiprocessing import Pool
from typing import Any, Optional

import intervaltree  # type: ignore
import mappy as mp
from intervaltree import Interval  # type: ignore

from .identify_breaks_v5 import call_breakpoints
from .io import GroupedResults, load_gaf_to_grouped_reads
from .merge_break_pts_v3 import merge_breaks

cigar_pattern = re.compile(r"(\d+)([=XIDM])")
path_seg_pattern = re.compile(r"([><])([^><:\s]+):(\d+)-(\d+)")
# ds:Z: tag example
#      *ag: SNV
#      :234: reference bp
#      -[gcgccgcgccggcgcaggcgcagagag]: TSD deletion, can be left or right end
#      +tggagggactgcccagt: insertion
ds_pattern = re.compile(r"([\+\-\*:])([A-Za-z\[\]0-9]+)")
path_seg_pattern = re.compile(r"([><])([^><:\s]+):(\d+)-(\d+)")
re_tsd = re.compile(r"(\[([A-Za-z]+)\])?([A-Za-z]+)(\[([A-Za-z]+)\])?")


class StartPosException(Exception):
    pass


class EndPosException(Exception):
    pass


def parse_one_line(line, min_mapq=30, min_len=100, verbose=False, lineno=0, ds=True):
    read = line.strip().split()
    parsed_read = AlignedRead(read)
    large_indels_one_read = parsed_read.get_indels(
        min_mapq=min_mapq, min_len=min_len, dbg=verbose, lineno=lineno, ds=ds
    )
    return large_indels_one_read


def parse_grouped_reads(
    grouped_reads, min_mapq=30, min_indel_len=100, verbose=False, lineno=0, ds=True
):
    all_indels = []
    for read in grouped_reads:
        parsed_read = AlignedRead(read)
        large_indels_one_read = parsed_read.get_indels(
            min_mapq=min_mapq, min_len=min_indel_len, dbg=verbose, lineno=lineno, ds=ds
        )
        all_indels += large_indels_one_read

    if len(grouped_reads) > 1:
        all_breaks = call_breakpoints(grouped_reads)
    else:
        all_breaks = []
    return GroupedResults(all_indels, all_breaks)


class GafParser(object):
    """Gaf file parser"""

    def __init__(
        self,
        gaf_paths: list[str],
        output: str,
        assembly: str,
        vntr: Optional[str] = None,
        cent: Optional[str] = None,
        l1: Optional[str] = None,
    ) -> None:
        """
        parameters
        ------------
        gaf_paths: a list of gaf files, restricted to tumor[0] and normal[1] pair
        """
        self.gaf_paths = gaf_paths
        assert len(self.gaf_paths) >= 1
        self.output = output
        self.vntr = vntr
        self.cent = cent
        self.l1 = l1
        self.assembly = assembly

    def parse_sv_on_group_reads(
        self,
        min_mapq: int = 5,
        min_indel_len: int = 50,
        min_map_len: int = 2000,
        verbose: bool = False,
        n_cpus: int = 4,
        ds: bool = True,
    ) -> None:
        self.ds = ds

        if os.path.exists(f"{self.output}.bed.gz"):
            return

        if len(self.gaf_paths) == 1:
            read_tags = [""]
        elif len(self.gaf_paths) == 2:
            read_tags = ["tumor", "normal"]

        lineno = 0
        indel_output = gzip.open(f"{self.output}_indel.bed.gz", "wt")
        breakpoint_output = gzip.open(f"{self.output}_brk.bed.gz", "wt")
        all_breaks = []
        for gaf_path, read_tag in zip(self.gaf_paths, read_tags):
            print(gaf_path, read_tag)
            if n_cpus == 1:
                for grouped_reads in load_gaf_to_grouped_reads(
                    gaf_path, min_mapq, min_map_len
                ):
                    if len(grouped_reads) > 1:
                        brks = call_breakpoints(grouped_reads)
                        for brk in brks:
                            brk[6] = (
                                f"{read_tag}_{brk[6]}" if read_tag != "" else brk[6]
                            )
                            all_breaks.append(brk)
                            brk_row_str = "\t".join(map(str, brk))
                            breakpoint_output.write(f"{brk_row_str}\n")
                    for read in grouped_reads:
                        lineno += 1
                        parsed_read = AlignedRead(read)
                        large_indels_one_read = parsed_read.get_indels(
                            min_mapq=min_mapq,
                            min_len=min_indel_len,
                            dbg=verbose,
                            lineno=lineno,
                            ds=ds,
                        )
                        for indel in large_indels_one_read:
                            indel[3] = (
                                f"{read_tag}_{indel[3]}" if read_tag != "" else indel[3]
                            )
                            indel_row_str = "\t".join(map(str, indel))
                            indel_output.write(f"{indel_row_str}\n")
                indel_output.close()
                breakpoint_output.close()
            else:
                fin = load_gaf_to_grouped_reads(gaf_path, min_mapq, min_map_len)
                with Pool(n_cpus) as pool:
                    parser = functools.partial(
                        parse_grouped_reads,
                        ds=ds,
                        min_mapq=min_mapq,
                        min_indel_len=min_indel_len,
                        verbose=verbose,
                    )
                    for group_results in pool.imap(parser, fin, chunksize=300000):
                        indel_results = group_results.indels
                        for indel in indel_results:
                            indel[3] = (
                                f"{read_tag}_{indel[3]}" if read_tag != "" else indel[3]
                            )
                            indel_row_str = "\t".join(map(str, indel))
                            indel_output.write(f"{indel_row_str}\n")

                        brk_results = group_results.breakpoints
                        for brk in brk_results:
                            brk[6] = (
                                f"{read_tag}_{brk[6]}" if read_tag != "" else brk[6]
                            )
                            all_breaks.append(brk)
                            brk_row_str = "\t".join(map(str, brk))
                            breakpoint_output.write(f"{brk_row_str}\n")
                indel_output.close()
                breakpoint_output.close()
        return all_breaks

    def merge_breakpts(self, all_breaks):
        self.breakpt_file = f"{self.output}_mergedbreaks.bed.gz"
        break_merged_file = gzip.open(f"{self.output}_mergedbreaks.bed.gz", "wt")
        break_merged_file.write("\n".join(merge_breaks(all_breaks)))
        break_merged_file.close()

    def parse_indel(
        self,
        min_mapq: int = 5,
        min_len: int = 100,
        verbose: bool = False,
        n_cpus: int = 4,
        ds: bool = True,
    ) -> None:
        self.ds = ds

        lineno = 0

        # if os.path.exists(f"{self.output}.bed.gz"):
        #    return

        output = gzip.open(f"{self.output}_indel.bed.gz", "wt")

        if len(self.gaf_paths) == 1:
            read_tags = [""]
        elif len(self.gaf_paths) == 2:
            read_tags = ["tumor", "normal"]

        for gaf_path, read_tag in zip(self.gaf_paths, read_tags):
            print(gaf_path, read_tag)
            with open(gaf_path) as fin:
                if n_cpus == 1:
                    for line in fin:
                        lineno += 1
                        read = line.strip().split()
                        parsed_read = AlignedRead(read)
                        large_indels_one_read = parsed_read.get_indels(
                            min_mapq=min_mapq,
                            min_len=min_len,
                            dbg=verbose,
                            lineno=lineno,
                            ds=ds,
                        )
                        for indel in large_indels_one_read:
                            indel[3] = (
                                f"{read_tag}_{indel[3]}" if read_tag != "" else indel[3]
                            )
                            indel_row_str = "\t".join(map(str, indel))
                            output.write(f"{indel_row_str}\n")
                else:
                    with Pool(n_cpus) as pool:
                        parser = functools.partial(
                            parse_one_line, ds=ds, min_mapq=min_mapq, min_len=min_len
                        )
                        for result in pool.imap(parser, fin, chunksize=10000):
                            for indel in result:
                                output.write(
                                    f"{indel[0]}\t{indel[1]}\t{indel[2]}\t{read_tag}_{indel[3]}\t{indel[4]}\t{indel[5]}\t{indel[6]}\t{indel[7]}\t{indel[8]}\t{indel[9]}\n"
                                )
        output.close()

    def bed2vcf(self, command: str = ""):
        """Convert merged indel bed file to a vcf file and index by bgzip and tabix"""
        current_time = datetime.now()
        formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")

        vcf_output = open(f"{self.output}_mergedindel.vcf", "w")
        # refer to severus and sniffle2
        hdr = ["##fileformat=VCFv4.2"]
        hdr.append("##source=gaftools1.0")
        hdr.append(f'##command="{command}"')
        hdr.append(f'##fileDate="{formatted_time}"')
        # remove lengths for both hg38 and chm13

        # chm13graph,chm13linear,grch37graph,grch37linear,grch38graph,grch38linear
        if self.assembly in [
            "chm13graph",
            "chm13linear",
            "grch38graph",
            "grch38linear",
        ]:
            hdr.append(
                """##contig=<ID=chr1>
##contig=<ID=chr2>
##contig=<ID=chr3>
##contig=<ID=chr4>
##contig=<ID=chr5>
##contig=<ID=chr6>
##contig=<ID=chr7>
##contig=<ID=chr8>
##contig=<ID=chr9>
##contig=<ID=chr10>
##contig=<ID=chr11>
##contig=<ID=chr12>
##contig=<ID=chr13>
##contig=<ID=chr14>
##contig=<ID=chr15>
##contig=<ID=chr16>
##contig=<ID=chr17>
##contig=<ID=chr18>
##contig=<ID=chr19>
##contig=<ID=chr20>
##contig=<ID=chr21>
##contig=<ID=chr22>
##contig=<ID=chrX>
##contig=<ID=chrY>
##contig=<ID=chrM>"""
            )
        elif self.assembly in ["grch37graph", "grch37linear"]:
            hdr.append(
                """##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
##contig=<ID=X>
##contig=<ID=Y>
##contig=<ID=MT>"""
            )

        hdr.append(
            """##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=BND,Description="Breakend; Translocation">"""
        )
        hdr.append('##FILTER=<ID=PASS,Description="All filters passed">')
        hdr.append(
            """##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="SV with precise breakpoints coordinates and length">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="SV with imprecise breakpoints coordinates and length">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=TSDLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=POLYALEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness">
##INFO=<ID=DETAILED_TYPE,Number=1,Type=Integer,Description="Detailed type of the SV">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of supporting reads">
##INFO=<ID=VNTR,Number=1,Type=Integer,Description="VNTR support">
##INFO=<ID=CENTROMERE,Number=1,Type=Integer,Description="Overlap with Centromere">
##INFO=<ID=SUPPREAD,Number=1,Type=Integer,Description="Number of supporting reads">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">
##INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotyping quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">"""
        )

        for h in hdr:
            vcf_output.write(f"{h}\n")
        # if len(self.gaf_paths) > 1:
        #     samples_header = list(map(lambda x: os.path.basename(x).replace(".gaf", ""), self.gaf_paths))
        # else:
        samples_header = ["SAMPLE"]
        columns = [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
        ] + samples_header
        columns_str = "\t".join(columns)
        vcf_output.write(f"{columns_str}\n")

        # merged indel bed column names
        # [ctg, start, end, len_str,  # average indel length
        #  tsd_len_str,  # average tsd length
        #  polyA_len_str,  # average polyA length
        #  indel_seq,  # randomly pick one
        #  n, ".", f"mq:i:{mapq}", f"cf:i:{nf}", f"cr:i:{nr}", vntr, centromere, f"rd:Z:+/-readname"]
        with gzip.open(f"{self.output}_mergedindel.bed.gz", "rt") as merge_indel_file:
            indel_num = 0
            for line in merge_indel_file:
                line = line.strip().split("\t")
                if not (
                    re.match(r"^[><HNC]", line[0]) or "#" in line[0] or "_" in line[0]
                ):  # remove non-linear genome
                    type = "INS" if int(line[3]) > 0 else "DEL"
                    start = int(line[1]) + 1
                    end = int(line[2])

                    # NOTE: shall we uppercase this base pair
                    # QUAL is assigned to . now, not sure how to compute it without PHRED score in bam since we input GAF/PAF, maybe not necessary for filter
                    # FILTER is set to PASS since we output only INDEL with supported read >= 3
                    # SV ID
                    ALT = line[6].upper() if type == "INS" else "<DEL>"
                    # Sniffle2 INFO: PRECISE;SVTYPE=INS;SVLEN=76;END=10872;SUPPORT=9;COVERAGE=8,9,9,11,10;STRAND=+-;AF=1.000;STDEV_LEN=8.319;STDEV_POS=4.930;SUPPORT_LONG=0
                    # Severus INFO: IMPRECISE;SVTYPE=INS;SVLEN=246;CHR2=chr1;END=34802253;DETAILED_TYPE=None;INSLEN=0;MAPQ=60;SUPPREAD=29;HVAF=0.00|0.00|0.50;CLUSTERID=severus_3
                    # How to tell if detailed type is BFB_foldback
                    # No need for CHR2 now since no breakpoint added yet
                    # Not sure how other tool define PRECISE/IMPRECISE
                    mapq = line[9].replace("mq:i:", "")
                    reads = line[-1].replace("rd:Z:", "")
                    VNTR = line[-3] == "True"
                    CENT = line[-2] == "True"
                    # GT:GQ:VAF:DR:DV 0|1:279:0.45:36:29
                    # Not sure how to decide GT, GQ and VAF, DR yet
                    vcf_output.write(
                        f"{line[0]}\t{start}\tgaftools1.{type}.{indel_num}\tN\t{ALT}\t.\tPASS\tEND={end}SVTYPE={type};SVLEN={line[3]};TSDLEN={line[4]};POLYALEN={line[5]};DETAILED_TYPE=None;MAPQ={mapq};SUPPREAD={line[7]};VNTR={VNTR};CENTROMERE={CENT};RNAMES={reads}\tGT:GQ:VAF:DR:DV\t.:.:.:.:{line[7]}\n"
                    )
                    indel_num += 1
        vcf_output.close()

    def parse_bed(self, bed) -> dict:
        bed_dict: dict[Any, Any] = {}  # type: ignore
        with open(bed) as bed_file:
            for line in bed_file:
                line = line.strip().split()
                if line[0] not in bed_dict:
                    bed_dict[line[0]] = intervaltree.IntervalTree()
                bed_dict[line[0]].add(
                    Interval(
                        int(line[1]),
                        int(line[2]),
                    )
                )
        return bed_dict

    def merge_indel(
        self,
        min_mapq: int = 5,
        min_cnt: int = 1,
        win_size: int = 100,
        max_diff: float = 0.05,
    ):
        """Merge the similar indels

        If the two indel overlap and have the same type and similar length, then merge together
        """

        if self.l1 is not None:
            minimap2 = mp.Aligner(self.l1)
            if not minimap2:
                raise Exception("Error: failed to load/build index")

        if self.vntr is not None:
            vntr_sites_dict = self.parse_bed(self.vntr)
        if self.cent is not None:
            cent_sites_dict = self.parse_bed(self.cent)

        self.indel_file = f"{self.output}_mergedindel.bed.gz"
        merged_output = gzip.open(f"{self.output}_mergedindel.bed.gz", "wt")

        def print_bed(ctg, t):
            # t is bj below: [start, end, ., mapq, ., indel length, [indel length], [mapq], [tsd length], [polyA length], [indel seq], [strandreadname]]
            n = len(t[-1])  # number of reads
            if n < min_cnt:
                return

            length, mapq, tsd_length, polyA_length = 0, 0, 0, 0
            # forward and reverse strand reads
            nf, nr = 0, 0
            for i in range(n):
                length += t[6][i]
                mapq += t[7][i]

                if self.ds:  # gaf file
                    tsd_length += t[8][i]
                    polyA_length += t[9][i]

                # 11th element
                if t[-1][i][0] == "+":
                    nf += 1
                else:
                    nr += 1

            indel_seq = "."
            # TODO: calculate concensus sequence in the next version
            if self.ds:  # gaf file, pick a arbitrary sequence
                indel_seq = t[10][0]

            mapq = math.floor(mapq / n + 0.499)
            if mapq < min_mapq:
                return

            length = math.floor(length / n + 0.499)
            len_str = f"+{length}" if length > 0 else str(length)

            if self.ds:
                tsd_length = math.floor(tsd_length / n + 0.499)
                # Check if TSD can be deletion
                tsd_len_str = f"+{tsd_length}" if tsd_length > 0 else str(tsd_length)
                polyA_length = math.floor(polyA_length / n + 0.499)
                # Check if polyA can be deletion
                polyA_len_str = (
                    f"+{polyA_length}" if polyA_length > 0 else str(polyA_length)
                )
            else:
                tsd_len_str = "0"
                polyA_len_str = "0"

            cent_hit = False
            vntr_hit = False
            if self.vntr is not None:
                if ctg in vntr_sites_dict:
                    vntr_hit = len(vntr_sites_dict[ctg].overlap(t[0], t[1])) > 0
            if self.cent is not None:
                if ctg in cent_sites_dict:
                    cent_hit = len(cent_sites_dict[ctg].overlap(t[0], t[1])) > 0

            if self.l1 is not None:
                if length > 0:  # only apply to insertion for L1
                    hits = minimap2.map(indel_seq)
                    hits_list = []
                    for hit in hits:
                        hits_list.append((hit.ctg, hit.mlen / hit.ctg_len))
                    hits_list.sort(key=lambda x: x[-1])
                    if len(hits_list) == 0:
                        best_hit = [".", 0]
                    else:
                        best_hit = hits_list[-1]
                else:
                    best_hit = [".", 0]
            else:
                best_hit = [".", 0]

            output_str = "\t".join(
                map(
                    str,
                    [
                        ctg,
                        t[0],
                        t[1],
                        len_str,  # average indel length
                        tsd_len_str,  # average tsd length
                        polyA_len_str,  # average polyA length
                        indel_seq,  # randomly pick one
                        n,
                        ".",
                        f"mq:i:{mapq}",
                        f"cf:i:{nf}",
                        f"cr:i:{nr}",
                        vntr_hit,
                        cent_hit,
                        f"{best_hit[0]}:{best_hit[1]}",
                        f"rd:Z:{','.join(t[-1])}",  # +/-readname
                    ],
                )
            )
            merged_output.write(f"{output_str}\n")

        # [path, start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
        merged_indel_dict = {}
        with gzip.open(f"{self.output}_indel.bed.gz", "rt") as indel_read_output:
            for line in indel_read_output:
                t = line.strip().split("\t")
                t[1] = int(t[1])
                t[2] = int(t[2])
                t[4] = int(t[4])
                t[6] = int(t[6])
                t[7] = int(t[7])
                t[8] = int(t[8])
                ctg = t.pop(0)
                if ctg not in merged_indel_dict:
                    merged_indel_dict[ctg] = []
                merged_indel_dict[ctg].append(t)

        # iterate through contig
        for ctg in merged_indel_dict:
            # key: contig, value:[start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
            # sort indels by start, end, indel length, and strand, read name
            # this will fix the merge indel issue
            merged_indel_dict[ctg].sort(key=lambda x: (x[0], x[1], x[5], x[4], x[2]))

            # sanity check reads before merging
            # merged_indel_dict[ctg].sort(key=lambda x: x[0])
            # for key in ctg:
            #     for indel in merged_indel_dict[ctg]:
            #         test='\t'.join(map(str, indel))
            #         print(f"{test}")

            # all indels
            a = merged_indel_dict[ctg]
            # candidate merged indel
            b = []
            # iterate through indel
            for i in range(len(a)):
                # ith indel or L1 element
                ai = a[i]

                # compare ai with existing merged indels
                while len(b) > 0:
                    # indel ai start position larger than indel end position
                    if ai[0] - b[0][1] > win_size:
                        t = b.pop(0)
                        print_bed(ctg, t)
                    else:
                        break

                merge_j = -1
                # where merge starts
                for j in range(len(b) - 1, -1, -1):
                    # start from the closest indel to ai
                    # ai: [start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
                    # bj: [start, end, ., mapq, ., indel length, [indel length], [mapq], [tsd length], [polyA length], [indel seq], [strandreadname]]
                    bj = b[j]

                    if bj[5] * ai[5] <= 0:
                        # bj and ai are not the same SV type,
                        # deletion cannot merge with insertion
                        # it is not break here...
                        # continue to compare with next candidate merged indel
                        continue

                    # ai and bj indel length
                    la = ai[5] if ai[5] > 0 else -ai[5]
                    # BUG!
                    # lb = bj[5] if bj[5] > 0 else -bj[5]

                    # insertion
                    if all([bjj > 0 for bjj in bj[6]]):
                        # conservative minimum distance between ai and a list of bj indels
                        diff = min([abs(la - bjj) for bjj in bj[6]])
                        # ai and bj indel maximum length
                        max_indel = max(bj[6] + [la])
                    # deletion
                    elif all([bjj < 0 for bjj in bj[6]]):
                        # minimum distance between ai and a list of bj indels
                        diff = min([abs(la - (-bjj)) for bjj in bj[6]])
                        # ai and bj indel maximum length
                        max_indel = max(-min(bj[6]), la)
                    else:
                        raise Exception("Bug!")

                    # two indels are different skip merge
                    # sorted order may change this result
                    if diff > max_indel * max_diff:
                        # continue to compare with next candidate merged indel
                        continue

                    # indel length
                    bj[6].append(ai[5])
                    # mapq
                    bj[7].append(ai[3])

                    if self.ds:  # gaf file
                        # tsd length
                        bj[8].append(ai[6])
                        # polyA length
                        bj[9].append(ai[7])
                        # indel seq
                        bj[10].append(ai[6])
                        # read name
                        bj[11].append(f"{ai[4]}{ai[2]}")
                    else:  # paf file
                        bj[8].append(f"{ai[4]}{ai[2]}")

                    # reassign the end point
                    # if bj contains the ai indel
                    # use bj end else ai end
                    bj[1] = bj[1] if bj[1] > ai[1] else ai[1]
                    merge_j = j
                    # ai indel only merge one time successfully
                    # with the last bj indel in the candidate merged list
                    # then break
                    break

                # empty indel b list, initialize with first ai indel
                # or ai is different from bj based on max_diff and sv_type
                # put in the ai into b
                if merge_j < 0:
                    # ai: [start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
                    # bj: [start, end, ., mapq, ., indel length, [indel length], [mapq], [tsd length], [polyA length], [indel seq], [strandreadname]]
                    if self.ds:  # gaf file with indel seq
                        b.append(
                            [
                                ai[0],
                                ai[
                                    1
                                ],  # always the most recent/largest coordinate indel
                                ".",
                                ai[3],
                                ".",
                                ai[
                                    5
                                ],  # NOTED issue: only a constant first indel length is used for max diff compute
                                [ai[5]],
                                [ai[3]],
                                [ai[6]],
                                [ai[7]],
                                [ai[8]],
                                [f"{ai[4]}{ai[2]}"],
                            ]
                        )
                    else:  # paf file
                        # paf output from minimap2 do not have ai[6-8]: tsd, polyA, indel seq yet
                        # ai [start, end, readname, mapq, strand, indel length]
                        # bj [start, end, ., mapq, ., indel length, [indel length], [mapq], [strandreadname]]
                        b.append(
                            [
                                ai[0],
                                ai[1],
                                ".",
                                ai[3],
                                ".",
                                ai[5],
                                [ai[5]],
                                [ai[3]],
                                [f"{ai[4]}{ai[2]}"],
                            ]
                        )

            # output indel that are skipped in the merging for loop
            while len(b) > 0:
                t = b.pop(0)
                print_bed(ctg, t)
        merged_output.close()

    def parse_tsd(
        self, min_mapq: int = 5, min_len: int = 100, verbose: bool = False
    ) -> None:
        output = open(f"{self.output}_tsd.bed", "w")

        with open(self.gaf_path) as fin:
            for line in fin:
                read = line.strip().split()
                parsed_read = AlignedRead(read)
                large_tsds = parsed_read.get_tsd(
                    min_mapq=min_mapq, min_len=min_len, dbg=verbose
                )
                if large_tsds:
                    for tsd in large_tsds:
                        tsd_row_str = "\t".join(map(str, tsd))
                        output.write(f"{tsd_row_str}\n")
        output.close()


class AlignedRead(object):
    """Parse one read from the gaf file"""

    def __init__(self, read: list):
        self.read = read

    @property
    def querylen(self):
        return int(self.read[1])

    @property
    def strand(self):
        return self.read[4]

    @property
    def path(self):
        return self.read[5]

    @property
    def mapq(self):
        return int(self.read[11])

    @property
    def cigar(self):
        # to be compatible with latest minigraph
        # https://github.com/lh3/minigraph/commit/54c74fc76b8b946e74aa6f44b559503b5821f12d
        if "ds:Z" in self.read[-1]:
            return self.read[-2]
        return self.read[-1]

    @property
    def path_start(self):
        """https://github.com/lh3/gfatools/blob/master/doc/rGFA.md"""
        return int(self.read[7])

    @property
    def path_end(self):
        return int(self.read[8])

    @property
    def ds_Z(self):
        """Parse the ds:Z: tag"""
        # to be compatible with latest minigraph
        # https://github.com/lh3/minigraph/commit/54c74fc76b8b946e74aa6f44b559503b5821f12d
        if "ds:Z" in self.read[-1]:
            return self.read[-1][5:]
        raise Exception("Upgrade your minigraph to 0.20-r574-dirty or later")

    def convert_type(self):
        """Convert coordinate column into integers"""
        # convert element types
        for i in range(1, 4):
            self.read[i] = int(self.read[i])
        for i in range(6, 12):
            self.read[i] = int(self.read[i])

    def get_indels(
        self,
        lineno: int,
        min_mapq: int = 5,
        min_len: int = 100,
        max_cnt: int = 5,
        min_frac: float = 0.7,
        ds: bool = True,
        polyA_pen: int = 5,
        polyA_drop: int = 100,
        dbg: bool = False,
    ) -> list:
        """Get indels from one aligned long read in GAF or PAF file

        Update script from https://github.com/lh3/minigraph/blob/master/misc/mgutils-es6.js#L232-L401
        """

        # at least 12 columns for one read
        # read should be minigraph or minimap2 output
        if len(self.read) < 12:
            return []

        self.convert_type()
        # low mapping quality
        if self.mapq < min_mapq:
            return []

        # mapped fraction is low
        if self.read[3] - self.read[2] < self.read[1] * min_frac:
            return []

        cg_segs = cigar_pattern.findall(self.cigar)

        # record cigar indels
        a = []
        x = self.path_start
        for length, op in cg_segs:
            length = int(length)
            if dbg:
                print("X0", x, length, op)
            if length >= min_len:
                if op == "I":
                    # [start, end, indel length, tsd length, polyA length, tsd seq, indel seq]
                    a.append([x - 1, x + 1, length, 0, 0, "", ""])
                elif op == "D":
                    a.append([x, x + length, -length, 0, 0, "", ""])
            if op in ["M", "D", "=", "X"]:
                x += length

        if dbg:
            print("X0", a)
            print("X0", self.path_start)
            print("X0", self.cigar, self.read[0])

        # No indel cigar or too many
        if len(a) == 0 or len(a) > max_cnt:
            return []

        if ds:
            i = 0  # indel index for one read
            x = self.path_start
            ds_Z = self.ds_Z
            ds_segs_iter = ds_pattern.findall(ds_Z)
            for op, ds_str in ds_segs_iter:
                seq = re.sub(r"[\]\[]", "", ds_str) if op in ["+", "-"] else ""
                if op == ":":
                    length = int(ds_str)
                elif op == "*":
                    length = 1
                elif op in ["+", "-"]:
                    length = len(seq)
                else:
                    length = -1
                assert length > 0, "not found ds:Z supported tag"
                if length >= min_len:
                    if op == "+":  # insertion
                        # a: [start, end, indel length, tsd length, polyA length, tsd seq, indel seq]
                        if a[i][0] != x - 1 or a[i][1] != x + 1 or a[i][2] != length:
                            raise Exception(
                                f"CIGAR and ds insertion not consistent line number {lineno}"
                            )
                        a[i][5] = ds_str
                        i += 1
                    elif op == "-":  # deletion
                        if a[i][0] != x or a[i][1] != x + length or a[i][2] != -length:
                            raise Exception(
                                f"CIGAR and ds deletion not consistent {lineno}"
                            )
                        a[i][5] = ds_str
                        i += 1
                if op == "*" or op == ":" or op == "-":
                    x += length

            for i in range(len(a)):
                # m[0]: [tsd_indel], m[1]: tsd_indel, m[2]: indel, m[3]: [tsd_indel], m[4]: tsd_indel
                m = re_tsd.search(a[i][5])
                m = m.groups()
                if m is None:
                    raise Exception("Bug!")
                left_tsd = m[1] if m[1] is not None else ""
                right_tsd = m[4] if m[4] is not None else ""

                # tsd sequences right + left is the reference sequences
                # due to aligner design
                # a: [start, end, indel length, tsd length, polyA length, tsd seq, indel seq]
                tsd = right_tsd + left_tsd
                a[i][3] = len(tsd)

                a[i][5] = tsd if len(tsd) > 0 else "."
                internal_seq = m[2]
                internal_seqlen = len(internal_seq)
                a[i][6] = internal_seq if internal_seqlen > 0 else "."
                if internal_seqlen > 0:  # internal insertion larger than 1bp
                    polyA_len, polyT_len = 0, 0
                    polyA_max, polyT_max = 0, 0
                    score, max = 0, 0
                    max_j = internal_seqlen
                    for j in range(internal_seqlen - 1, -1, -1):  # 3' -> 5' sense?
                        if internal_seq[j] in ["A", "a"]:
                            score += 1
                        else:
                            score -= polyA_pen
                        if score > max:
                            max = score
                            max_j = j
                        elif max - score > polyA_drop:
                            break
                    polyA_len = internal_seqlen - max_j
                    polyA_max = max
                    score, max = 0, 0
                    max_j = -1
                    for j in range(internal_seqlen):  # 5' -> 3' antisense?
                        if internal_seq[j] in ["T", "t"]:
                            score += 1
                        else:
                            score -= polyA_pen
                        if score > max:
                            max = score
                            max_j = j
                        elif max - score > polyA_drop:
                            break
                    polyT_len = max_j + 1
                    polyT_max = max
                    a[i][4] = polyA_len if polyA_max >= polyT_max else -polyT_len

        # path segments
        seg = []
        if re.match(r"[><]", self.path):
            assert self.strand == "+", "reverse strand on path"
            y = 0
            for m in path_seg_pattern.findall(self.path):
                st, en = int(m[2]), int(m[3])
                # [path, absolute path start, absolute path end, path orientation, relative seg start, relative seg end]
                seg.append([m[1], st, en, 1 if m[0] == ">" else -1, y, y + (en - st)])
                # accumulate seg start
                y += en - st
        else:
            # https://github.com/lh3/miniasm/blob/master/PAF.md
            # https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
            # [path, 0, path length, 1, 0, path length]
            seg.append([self.path, 0, self.read[6], 1, 0, self.read[6]])

        if dbg:
            print("X1", seg)

        # starts and ends points for the indels on the segments
        # also record start and end segments index
        # these are used to overlap between path segment and indels
        sts, ens = [], []
        for i in range(len(a)):
            k = 0
            # segment k relative end <= indel relative start
            while k < len(seg) and seg[k][5] <= a[i][0]:
                k += 1
            if k == len(seg):
                raise Exception("failed to find start position")
            # relative distance between indel start and segment start
            start_l = a[i][0] - seg[k][4]
            # graph path > or linear genome
            if seg[k][3] > 0:
                # [seg index, absolute path start + distance to indel]
                sts.append([k, seg[k][1] + start_l])
            # graph path <, ask Heng
            else:
                # [seg index, absolute path end - distance to indel]
                sts.append([k, seg[k][2] - start_l])

        for i in range(len(a)):
            k = 0
            # segment k relative end <= indel relative end
            while k < len(seg) and seg[k][5] <= a[i][1]:
                k += 1
            if k == len(seg):
                raise Exception("failed to find end position")
            # relative distance between indel end and segment start
            end_l = a[i][1] - seg[k][4]
            # graph path > or linear genome
            if seg[k][3] > 0:
                # [seg index, absolute path start + distance to indel]
                ens.append([k, seg[k][1] + end_l])
            # graph path <, ask Heng
            else:
                ens.append([k, seg[k][2] - end_l])

        output_indels = []
        for i in range(len(a)):
            if dbg:
                print("X2", a[i][0], a[i][1], sts[i][0], ens[i][0])

            # indel on the same segment of the path
            if sts[i][0] == ens[i][0]:
                s = seg[sts[i][0]]
                # seg starts with > or linear reference
                if s[3] > 0:
                    strand = self.strand
                # seg starts with <
                else:
                    strand = "-" if self.strand == "+" else "+"

                if sts[i][1] < ens[i][1]:
                    start = sts[i][1]
                else:
                    start = ens[i][1]

                if sts[i][1] > ens[i][1]:
                    end = sts[i][1]
                else:
                    end = ens[i][1]

                # [chrom, start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
                output_indels.append(
                    [
                        s[0],
                        start,
                        end,
                        self.read[0],
                        self.mapq,
                        strand,
                        a[i][2],
                        a[i][3],
                        a[i][4],
                        a[i][6],
                    ]
                )
            # indel on different segments
            else:
                path = []
                length = 0
                # iterate through segments
                # j from indel start seg index to ends seg index
                for j in range(sts[i][0], ens[i][0] + 1):
                    s = seg[j]
                    length += s[2] - s[1]
                    orientation = ">" if s[3] > 0 else "<"
                    # >chr1:1-2
                    path.append(orientation + s[0] + f":{s[1]}-{s[2]}")

                # relative seg start
                off = seg[sts[i][0]][4]
                # minigraph always use "+" in genome path
                # convert back to relative path coordinate
                # [chrom, start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
                output_indels.append(
                    [
                        "".join(path),
                        a[i][0] - off,
                        a[i][1] - off,
                        self.read[0],
                        self.mapq,
                        "+",
                        a[i][2],
                        a[i][3],
                        a[i][4],
                        a[i][6],
                    ]
                )
        return output_indels
