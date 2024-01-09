"""An module for parsing gfa alignment from minigraph
"""

import gzip
import math
import re
from multiprocessing import Pool

import cython

cigar_pattern = re.compile(r"(\d+)([=XIDM])")
path_seg_pattern = re.compile(r"([><])([^><:\s]+):(\d+)-(\d+)")
ds_pattern = re.compile(r"([\+\-\*:])([A-Za-z\[\]0-9]+)")
path_seg_pattern = re.compile(r"([><])([^><:\s]+):(\d+)-(\d+)")
re_tsd = re.compile(r"(\[([A-Za-z]+)\])?([A-Za-z]+)(\[([A-Za-z]+)\])?")


class StartPosException(Exception):
    pass


class EndPosException(Exception):
    pass


def parse_one_line(line, min_mapq=30, min_len=100, verbose=False):
    read = line.strip().split()
    parsed_read = AlignedRead(read)
    large_indels_one_read = parsed_read.get_indels(
        min_mapq=min_mapq, min_len=min_len,
        max_cnt=5, min_frac=0.7, polyA_pen=5, polyA_drop=100
    )
    return large_indels_one_read


class GafParser:
    """Gaf file parser"""
    def __init__(self, gaf_path, output):
        self.gaf_path = gaf_path
        self.output = output

    def parse_indel(
        self,
        min_mapq: int = 5,
        min_len: int = 100,
        verbose: bool = False,
        n_cpus: int = 4,
        ):
        cdef:
            str line
            str indel_row_str
            list read
            list large_indels_one_read
            list indel

        output = gzip.open(f"{self.output}.bed.gz", "wt")
        with open(self.gaf_path) as fin:
            if n_cpus == 1:
                for line in fin:
                    read = line.strip().split()
                    parsed_read = AlignedRead(read)
                    large_indels_one_read = parsed_read.get_indels(
                        min_mapq=min_mapq, min_len=min_len, max_cnt=5, min_frac=0.7, polyA_pen=5, polyA_drop=100
                    )
                    for indel in large_indels_one_read:
                        indel_row_str = "\t".join(map(str, indel))
                        output.write(f"{indel_row_str}\n")
            else:
                with Pool(n_cpus) as pool:
                    for result in pool.imap(parse_one_line, fin, chunksize=10000):
                        for indel in result:
                            output.write(
                                f"{indel[0]}\t{indel[1]}\t{indel[2]}\t{indel[3]}\t{indel[4]}\t{indel[5]}\t{indel[6]}\t{indel[7]}\t{indel[8]}\t{indel[9]}\n"
                            )
        output.close()

    def merge_indel(
        self,
        min_mapq: int = 5,
        min_cnt: int = 1,
        win_size: int = 100,
        max_diff: float = 0.05,
    ):
        merged_output = gzip.open(f"{self.output}_mergedindel.bed.gz", "wt")

        def print_bed(ctg, t):
            # t [start, end, ., mapq, ., indel length, [indel length], [mapq], [strandreadname]]
            n = len(t[6])  # number of reads
            if n < min_cnt:
                return

            length, mapq = 0, 0
            nf, nr = 0, 0
            for i in range(n):
                length += t[6][i]
                mapq += t[7][i]
                if t[8][i][0] == "+":
                    nf += 1
                else:
                    nr += 1
            mapq = math.floor(mapq / n + 0.499)
            if mapq < min_mapq:
                return

            length = math.floor(length / n + 0.499)
            len_str = f"+{length}" if length > 0 else str(length)
            output_str = "\t".join(
                map(
                    str,
                    [
                        ctg,
                        t[0],
                        t[1],
                        len_str,
                        n,
                        ".",
                        f"mq:i:{mapq}",
                        f"cf:i:{nf}",
                        f"cr:i:{nr}",
                        f"rd:Z:{','.join(t[8])}",
                    ],
                )
            )
            merged_output.write(f"{output_str}\n")

        # [pathstr, start, end, readname, mapq, strand, indel length]
        merged_indel_dict = {}
        with gzip.open(f"{self.output}.bed.gz", "rt") as indel_read_output:
            for line in indel_read_output:
                t = line.strip().split("\t")
                t[1] = int(t[1])
                t[2] = int(t[2])
                t[4] = int(t[4])
                t[6] = int(t[6])
                ctg = t.pop(0)
                if ctg not in merged_indel_dict:
                    merged_indel_dict[ctg] = []
                merged_indel_dict[ctg].append(t)

        for ctg in merged_indel_dict:
            # sort indels by start
            merged_indel_dict[ctg].sort(key=lambda x: x[0])
            # key: contig, value: [start, end, readname, mapq, strand, indel length]

            a = merged_indel_dict[ctg]
            b = []
            for i in range(len(a)):
                # ith indel
                ai = a[i]

                while len(b) > 0:
                    # indel i start position
                    # larger than indel end position
                    if ai[0] - b[0][1] > win_size:
                        t = b.pop(0)
                        print_bed(ctg, t)
                    else:
                        break

                merge_j = -1
                # where merge starts
                for j in range(len(b) - 1, -1, -1):
                    bj = b[j]

                    if bj[5] * ai[5] <= 0:
                        # bj and ai are not the same type, deletion or insertion
                        continue

                    # ai and bj indel length
                    la = ai[5] if ai[5] > 0 else -ai[5]
                    lb = bj[5] if bj[5] > 0 else -bj[5]

                    diff = la - lb if la > lb else lb - la
                    max_indel = la if la > lb else lb

                    # two indels are different
                    if diff > max_indel * max_diff:
                        continue

                    # indel length
                    bj[6].append(ai[5])
                    # mapq
                    bj[7].append(ai[3])
                    # read name
                    bj[8].append(f"{ai[4]}{ai[2]}")
                    # reassign the end point
                    bj[1] = bj[1] if bj[1] > ai[1] else ai[1]
                    merge_j = j
                    break

                # empty indel b list
                # initialize with first ai
                if merge_j < 0:
                    # ai [start, end, readname, mapq, strand, indel length]
                    # bj [start, end, ., mapq, ., indel length, [indel length, mapq], [strandreadname]]
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
            while len(b) > 0:
                t = b.pop(0)
                print_bed(ctg, t)
        merged_output.close()
        return

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


cdef class AlignedRead:
    """Parse one read from the gaf file"""

    cdef list read

    def __init__(self, list read):
        self.read = read

    cdef void convert_type(self):
        """Convert coordinate column into integers"""
        # convert element types
        cdef unsigned int i

        for i in range(1, 4):
            self.read[i] = int(self.read[i])
        for i in range(6, 12):
            self.read[i] = int(self.read[i])

    cdef list get_indels(
        self,
        int min_mapq,
        int min_len,
        int max_cnt,
        double min_frac,
        int polyA_pen,
        int polyA_drop):
        """Get indels from one aligned long read in GAF or PAF file

        Update script from https://github.com/lh3/minigraph/blob/master/misc/mgutils-es6.js#L232-L401
        """

        cdef:
            list cg_segs
            list ds_segs_iter
            list a
            list seg
            list s
            list output_indels
            list sts
            list ens
            list path
            int x, length, y, internal_seqlen, max_j, k
            int i=0
            int start_l, end_l, start, end
            str length_str, op, ds_str, internal_seq, tsd, left_tsd, right_tsd

        # at least 12 columns for one read
        if len(self.read) < 12:
            return []

        self.convert_type()

        # low mapping quality
        if self.read[11] < min_mapq:
            return []

        # mapped fraction is low
        if self.read[3] - self.read[2] < self.read[1] * min_frac:
            return []

        cg_segs = cigar_pattern.findall(self.read[-2])

        # record cigar indels
        a = []
        x = self.read[7]
        for length_str, op in cg_segs:
            length = int(length_str)
            if length >= min_len:
                if op == "I":
                    # [start, end, indel length, tsd length, polyA length, tsd seq, indel seq]
                    a.append([x - 1, x + 1, length, 0, 0, "", ""])
                elif op == "D":
                    a.append([x, x + length, -length, 0, 0, "", ""])
            if op in ["M", "D", "=", "X"]:
                x += length

        # No indel cigar or too many
        if len(a) == 0 or len(a) > max_cnt:
            return []

        x = self.read[7]
        ds_Z = self.read[-1][5:]
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
                    a[i][5] = ds_str
                    i += 1
                elif op == "-":  # deletion
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
            tsd = right_tsd + left_tsd  # tsd sequences right + left?
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
        if re.match(r"[><]", self.read[5]):
            #assert self.strand == "+", "reverse strand on path"
            y = 0
            for m in path_seg_pattern.findall(self.read[5]):
                st, en = int(m[2]), int(m[3])
                # [path, absolute path start, absolute path end, path orientation, relative seg start, relative seg end]
                seg.append([m[1], st, en, 1 if m[0] == ">" else -1, y, y + (en - st)])
                # accumulate seg start
                y += en - st
        else:
            # https://github.com/lh3/miniasm/blob/master/PAF.md
            # https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
            # [path, 0, path length, 1, 0, path length]
            seg.append([self.read[5], 0, self.read[6], 1, 0, self.read[6]])

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
            # indel on the same segment of the path
            if sts[i][0] == ens[i][0]:
                s = seg[sts[i][0]]
                # seg starts with > or linear reference
                if s[3] > 0:
                    strand = self.read[4]
                # seg starts with <
                else:
                    strand = "-" if self.read[4] == "+" else "+"

                if sts[i][1] < ens[i][1]:
                    start = sts[i][1]
                else:
                    start = ens[i][1]

                if sts[i][1] > ens[i][1]:
                    end = sts[i][1]
                else:
                    end = ens[i][1]

                # [start, end, indel length, tsd length, polyA length, tsd seq, polyA seq]
                output_indels.append(
                    [
                        s[0],
                        start,
                        end,
                        self.read[0],
                        self.read[11],
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
                output_indels.append(
                    [
                        "".join(path),
                        a[i][0] - off,
                        a[i][1] - off,
                        self.read[0],
                        self.read[11],
                        "+",
                        a[i][2],
                        a[i][3],
                        a[i][4],
                        a[i][6],
                    ]
                )
        return output_indels
