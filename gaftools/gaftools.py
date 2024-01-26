"""An module for parsing gfa alignment from minigraph
"""

import gzip
import math
import re
from datetime import datetime
from multiprocessing import Pool

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


def parse_one_line(line, min_mapq=30, min_len=100, verbose=False, lineno=0):
    read = line.strip().split()
    parsed_read = AlignedRead(read)
    large_indels_one_read = parsed_read.get_indels(
        min_mapq=min_mapq, min_len=min_len, dbg=verbose, lineno=lineno
    )
    return large_indels_one_read


class GafParser(object):
    """Gaf file parser"""

    def __init__(self, gaf_paths: list[str], output: str):
        """
        parameters
        ------------
        gaf_paths: a list of gaf files, restricted to tumor[0] and normal[1] pair
        """
        self.gaf_paths = gaf_paths
        assert len(self.gaf_paths) >= 1
        self.output = output

    def parse_indel(
        self,
        min_mapq: int = 5,
        min_len: int = 100,
        verbose: bool = False,
        n_cpus: int = 4,
    ) -> None:
        lineno = 0
        output = gzip.open(f"{self.output}.bed.gz", "wt")

        if len(self.gaf_paths) == 1:
            read_tags = [""]
        elif len(self.gaf_paths) == 2:
            read_tags = ["tumor", "normal"]

        for gaf_path, read_tag in zip(self.gaf_paths, read_tags):
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
                        )
                        for indel in large_indels_one_read:
                            indel[3] = f"{read_tag}_{indel[3]}"
                            indel_row_str = "\t".join(map(str, indel))
                            output.write(f"{indel_row_str}\n")
                else:
                    with Pool(n_cpus) as pool:
                        for result in pool.imap(parse_one_line, fin, chunksize=10000):
                            for indel in result:
                                # indel_row_str = "\t".join(map(str, indel))
                                # output.write(f"{indel_row_str}\n")
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
##INFO=<ID=SUPPREAD,Number=1,Type=Integer,Description="Number of supporting reads">
##INFO=<ID=READS,Number=1,Type=Integer,Description="Number of supporting reads">
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
        #  n, ".", f"mq:i:{mapq}", f"cf:i:{nf}", f"cr:i:{nr}", f"rd:Z:+/-readname"]
        with gzip.open(f"{self.output}_mergedindel.bed.gz", "rt") as merge_indel_file:
            indel_num = 0
            for line in merge_indel_file:
                line = line.strip().split("\t")
                if not re.match(r"[><HGN]", line[0]):  # remove non-linear genome
                    type = "INS" if int(line[3]) > 0 else "DEL"
                    if type == "DEL":  # deletion use the left end
                        pos = line[1]
                    else:  # insertion use the middle point
                        pos = (int(line[1]) + int(line[2])) // 2

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
                    # GT:GQ:VAF:DR:DV 0|1:279:0.45:36:29
                    # Not sure how to decide GT, GQ and VAF, DR yet
                    vcf_output.write(
                        f"{line[0]}\t{pos}\tgaftools1.{type}.{indel_num}\tN\t{ALT}\t.\tPASS\tSVTYPE={type};SVLEN={line[3]};TSDLEN={line[4]};POLYALEN={line[5]};DETAILED_TYPE=None;MAPQ={mapq};SUPPREAD={line[7]};READS={reads}\tGT:GQ:VAF:DR:DV\t.:.:.:.:{line[7]}\n"
                    )
                    indel_num += 1
        vcf_output.close()

    def merge_indel(
        self,
        min_mapq: int = 5,
        min_cnt: int = 1,
        win_size: int = 100,
        max_diff: float = 0.05,
    ):
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
                tsd_length += t[8][i]
                polyA_length += t[9][i]

                # 11th element
                if t[-1][i][0] == "+":
                    nf += 1
                else:
                    nr += 1

            # pick the first indel sequence for the merged indel
            # TODO: calculate concensus sequence in the next version
            indel_seq = t[10][0]
            mapq = math.floor(mapq / n + 0.499)
            if mapq < min_mapq:
                return

            length = math.floor(length / n + 0.499)
            len_str = f"+{length}" if length > 0 else str(length)
            tsd_length = math.floor(tsd_length / n + 0.499)
            # Check if TSD can be deletion
            tsd_len_str = f"+{tsd_length}" if tsd_length > 0 else str(tsd_length)
            polyA_length = math.floor(polyA_length / n + 0.499)
            # Check if polyA can be deletion
            polyA_len_str = (
                f"+{polyA_length}" if polyA_length > 0 else str(polyA_length)
            )

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
                        f"rd:Z:{','.join(t[-1])}",  # +/-readname
                    ],
                )
            )
            merged_output.write(f"{output_str}\n")

        # [path, start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
        merged_indel_dict = {}
        with gzip.open(f"{self.output}.bed.gz", "rt") as indel_read_output:
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
            # sort indels by start
            merged_indel_dict[ctg].sort(key=lambda x: x[0])
            # key: contig, value:[start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]

            a = merged_indel_dict[ctg]
            b = []
            # iterate through indel
            for i in range(len(a)):
                # ith indel or L1
                ai = a[i]

                while len(b) > 0:
                    # indel ai start position
                    # larger than indel end position
                    if ai[0] - b[0][1] > win_size:
                        t = b.pop(0)
                        print_bed(ctg, t)
                    else:
                        break

                merge_j = -1
                # where merge starts
                for j in range(len(b) - 1, -1, -1):
                    # ai: [start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
                    # bj: [start, end, ., mapq, ., indel length, [indel length], [mapq], [tsd length], [polyA length], [indel seq], [strandreadname]]
                    bj = b[j]

                    if bj[5] * ai[5] <= 0:
                        # bj and ai are not the same type, deletion or insertion
                        continue

                    # ai and bj indel length
                    la = ai[5] if ai[5] > 0 else -ai[5]
                    lb = bj[5] if bj[5] > 0 else -bj[5]

                    # ai and bj indel length difference
                    diff = la - lb if la > lb else lb - la
                    # ai and bj indel maximum length
                    max_indel = la if la > lb else lb

                    # two indels are different skip merge
                    if diff > max_indel * max_diff:
                        continue

                    # indel length
                    bj[6].append(ai[5])
                    # mapq
                    bj[7].append(ai[3])
                    # tsd length
                    bj[8].append(ai[6])
                    # polyA length
                    bj[9].append(ai[7])
                    # indel seq
                    bj[10].append(ai[6])
                    # read name
                    bj[11].append(f"{ai[4]}{ai[2]}")
                    # reassign the end point
                    bj[1] = bj[1] if bj[1] > ai[1] else ai[1]
                    merge_j = j
                    break

                # empty indel b list
                # initialize with first ai indel
                if merge_j < 0:
                    # ai: [start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
                    # bj: [start, end, ., mapq, ., indel length, [indel length], [mapq], [tsd length], [polyA length], [indel seq], [strandreadname]]
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
                            [ai[6]],
                            [ai[7]],
                            [ai[8]],
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

    def get_tsd(
        self,
        min_mapq: int = 5,
        min_len: int = 100,
        max_cnt: int = 5,
        min_frac=0.7,
        dbg: bool = False,
    ) -> list:
        """Get TSD from one aligned long read in GAF or PAF file"""
        # similar to cs tag
        # https://github.com/lh3/minimap2#cs
        # /(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)+/
        ds_Z = self.ds_Z
        assert len(self.read) >= 12, "incomplete GAF file"
        self.convert_type()

        # low mapping quality
        if self.mapq < min_mapq:
            return []

        # mapped fraction is low, default 70%
        if self.read[3] - self.read[2] < self.read[1] * min_frac:
            return []

        ds_segs_iter = ds_pattern.findall(ds_Z)
        # record TSD
        a = []
        x = self.path_start
        for op, ds_str in ds_segs_iter:
            seq = re.sub(r"[\]\[]", "", ds_str) if op in ["+", "-"] else ""
            if op == ":":
                length = int(ds_str)
            elif op == "*":
                length = 1
            elif op in ["+", "-"]:
                length = len(seq)
            else:
                raise Exception("not found ds:Z supported tag")
            if op == "-" or op == "+":
                if length >= min_len:
                    if op == "-":  # deletion
                        a.append([x, x + length, -length])
                    elif op == "+":  # insertion
                        a.append([x - 1, x + 1, length])
            if op == "*" or op == ":" or op == "-":
                x += length

        if dbg:
            print(self.path_start, ds_Z, self.read[0])

        # No indel cigar or too many
        if len(a) == 0 or len(a) > max_cnt:
            return []
        if dbg:
            print("X0", a)

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
                try:
                    raise StartPosException("failed to find start position")
                except StartPosException:
                    print("start", k, len(seg), seg, a, self.read[0])
                    return []
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
                try:
                    raise EndPosException("failed to find end position")
                except EndPosException:
                    print("end", k, len(seg), seg, a, self.read[0])
                    return []
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

                output_indels.append(
                    [s[0], start, end, self.read[0], self.mapq, strand, a[i][2]]
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
                output_indels.append(
                    [
                        "".join(path),
                        a[i][0] - off,
                        a[i][1] - off,
                        self.read[0],
                        self.mapq,
                        "+",
                        a[i][2],
                    ]
                )
        return output_indels
