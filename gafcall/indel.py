import re
from dataclasses import dataclass

from .annotation import cal_cen_dist
from .regex import cigar_pattern, ds_pattern, path_seg_pattern, re_tsd


@dataclass
class indel_coord:
    st: int
    en: int
    len: int
    indel_seq: str = "."
    tsd_len: int = 0
    tsd_seq: str = "."
    polyA_len: int = 0
    int_seq: str = "."

    stl: int = 0
    enl: int = 0
    str: int = 0
    enr: int = 0


# NOTE: why ignore n
dna_dict = {"a": "t", "t": "a", "g": "c", "c": "g", "n": ""}


def mg_revcomp(s):
    return "".join([dna_dict[seq] for seq in list(s)[::-1]])


def cal_polyA_len(opt, int_seq):
    polyA_len, polyT_len = 0, 0
    polyA_max, polyT_max = 0, 0
    score, max = 0, 0
    max_j = len(int_seq)

    for j in range(len(int_seq) - 1, -1, -1):
        if int_seq[j] in ["A", "a"]:
            score += 1
        else:
            score -= opt.polyA_pen
        if score > max:
            max = score
            max_j = j
        elif max - score > opt.polyA_drop:
            break

    polyA_len = len(int_seq) - max_j
    polyA_max = max
    score, max = 0, 0
    max_j = -1
    for j in range(len(int_seq)):
        if int_seq[j] in ["T", "t"]:
            score += 1
        else:
            score -= opt.polyA_pen
        if score > max:
            max = score
            max_j = j
        elif max - score > opt.polyA_drop:
            break
    polyT_len = max_j + 1
    polyT_max = max
    return polyA_len if polyA_max >= polyT_max else -polyT_len


def get_indel(opt, z):
    """ """
    if len(z) == 0:
        return

    for j in range(len(z)):
        y = z[j]

        # ignore short alignment
        if y.qen - y.qst < y.qlen * opt.min_frac:
            continue

        cg_segs = cigar_pattern.findall(y.cg)
        # record cigar indels
        x = y.tst
        a = []
        for length, op in cg_segs:
            length = int(length)
            if length >= opt.min_len:
                if op == "I":
                    # [start, end, indel length, tsd length, polyA length, tsd seq, indel seq]
                    a.append(
                        indel_coord(
                            st=x - 1,
                            en=x + 1,
                            len=length,
                            indel_seq=".",
                            tsd_len=0,
                            tsd_seq=".",
                            polyA_len=0,
                            int_seq=".",
                        )
                    )
                elif op == "D":
                    a.append(
                        indel_coord(
                            st=x,
                            en=x + length,
                            len=-length,
                            indel_seq=".",
                            tsd_len=0,
                            tsd_seq=".",
                            polyA_len=0,
                            int_seq=".",
                        )
                    )
            if op in ["M", "D", "=", "X"]:
                x += length

        # No indel cigar or too many
        if len(a) == 0 or len(a) > opt.max_cnt:
            continue

        # NOTE: what is this for?
        # set stl/enl and str/enr
        for i in range(len(a)):
            a[i].stl = a[i].str = a[i].st
            a[i].enl = a[i].enr = a[i].en

        # parse ds:Z tag
        if y.ds:
            # indel index for one read
            i = 0
            x = y.tst
            ds_segs_iter = ds_pattern.findall(y.ds)
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
                # extract INDEL sequence with check consistency with CIGAR
                if length >= opt.min_len:
                    if op == "+":  # insertion
                        if a[i].st != x - 1 or a[i].en != x + 1 or a[i].len != length:
                            raise Exception(
                                "CIGAR and ds insertion not consistent line number"
                            )
                        a[i].indel_seq = ds_str
                        i += 1
                    elif op == "-":  # deletion
                        if a[i].st != x or a[i].en != x + length or a[i].len != -length:
                            raise Exception("CIGAR and ds deletion not consistent")
                        a[i].indel_seq = ds_str
                        i += 1
                if op == "*" or op == ":" or op == "-":
                    x += length

            # compute TSD and polyA lengths
            for i in range(len(a)):
                m = re_tsd.search(a[i].indel_seq)
                # m[0]: [tsd_indel], m[1]: tsd_indel, m[2]: indel, m[3]: [tsd_indel], m[4]: tsd_indel
                m = m.groups()
                if m is None:
                    raise Exception("Bug!")

                left_tsd = m[1] if m[1] is not None else ""
                right_tsd = m[4] if m[4] is not None else ""
                # tsd sequences right + left is the reference sequences
                # due to aligner design
                tsd = right_tsd + left_tsd

                a[i].tsd_len = len(tsd)
                a[i].tsd_seq = tsd
                # internal sequencing
                int_seq = m[2]
                a[i].int_seq = int_seq
                if len(int_seq) > 0:
                    a[i].polyA_len = cal_polyA_len(opt, int_seq)

                # NOTE: why is that
                #       left TSD?
                llen = len(m[1]) if m[1] is not None else 0
                #       right TSD?
                rlen = len(m[4]) if m[4] is not None else 0

                # NOTE: adjust TSD length?
                a[i].stl = a[i].st - rlen
                a[i].enl = a[i].en - rlen
                # NOTE: for right point, why not choose the minimum point by subtracting llen
                #       Could CIGAR-based indel have l2 (qgap) < 0
                a[i].str = a[i].st + llen
                a[i].enr = a[i].en + llen

        # reference segments in the path
        seg = []
        # with <>: this is a path
        if re.match(r"[><]", y.path):
            assert y.strand == "+", "reverse strand on path"
            x = 0
            for m in path_seg_pattern.findall(y.path):
                st, en = int(m[2]), int(m[3])
                # [path, absolute path start, absolute path end, path orientation, relative seg start, relative seg end]
                seg.append([m[1], st, en, 1 if m[0] == ">" else -1, x, x + (en - st)])
                # accumulate seg start
                x += en - st
        else:
            # https://github.com/lh3/miniasm/blob/master/PAF.md
            # https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
            # [path, 0, path length, 1, 0, path length]
            seg.append([y.path, 0, y.tlen, 1, 0, y.tlen])

        # starts and ends points for the indels on the segments
        # also record start and end segments index
        # these are used to overlap between path segment and indels
        sts, ens = [], []
        for i in range(len(a)):
            k = 0
            # segment k relative end <= indel relative start
            while k < len(seg) and seg[k][5] <= a[i].st:
                k += 1
            if k == len(seg):
                raise Exception("failed to find start position")
            # relative distance between indel start and segment start
            start_l = a[i].st - seg[k][4]
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
            while k < len(seg) and seg[k][5] <= a[i].en:
                k += 1
            if k == len(seg):
                raise Exception("failed to find end position")
            # relative distance between indel end and segment start
            end_l = a[i].en - seg[k][4]
            # graph path > or linear genome
            if seg[k][3] > 0:
                # [seg index, absolute path start + distance to indel]
                ens.append([k, seg[k][1] + end_l])
            # graph path <, ask Heng
            else:
                ens.append([k, seg[k][2] - end_l])

        for i in range(len(a)):
            if opt.dbg:
                print("X2", a[i].st, a[i].en)

            # NOTE: this is a new update that for reverse complement tsd, why?
            #       reverse complement sequence if strand is reverse
            if sts[i][0] == ens[i][0] and seg[sts[i][0]][3] < 0:
                a[i].polyA_len = -a[i].polyA_len
                a[i].tsd_seq = mg_revcomp(a[i].tsd_seq)
                a[i].int_seq = mg_revcomp(a[i].int_seq)

            info1 = ("SVTYPE=INS" if a[i].len > 0 else "SVTYPE=DEL") + (
                f";SVLEN={a[i].len};tsd_len={a[i].tsd_len};polyA_len={a[i].polyA_len}"
            )
            info2 = f"source={opt.name};tsd_seq={a[i].tsd_seq if len(a[i].tsd_seq)>0 else '.'};insert={a[i].int_seq if len(a[i].int_seq)>0 else '.'}"

            # indel on the same segment of the path
            if sts[i][0] == ens[i][0]:
                s = seg[sts[i][0]]

                # seg starts with > or linear reference
                if s[3] > 0:
                    strand2 = y.strand
                # seg starts with <
                else:
                    strand2 = "-" if y.strand == "+" else "+"

                if s[0] in opt.cen:
                    dist_st = cal_cen_dist(opt, s[0], sts[i][1])
                    dist_en = cal_cen_dist(opt, s[0], ens[i][1])
                    info1 += f";cen_dist={dist_st if dist_st < dist_en else dist_en}"

                if sts[i][1] < ens[i][1]:
                    start = sts[i][1]
                else:
                    start = ens[i][1]

                if sts[i][1] > ens[i][1]:
                    end = sts[i][1]
                else:
                    end = ens[i][1]
                # [chrom, start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
                print(
                    s[0],
                    start,
                    end,
                    y.qname,
                    y.mapq,
                    strand2,
                    f"{info1};{info2}",
                    sep="\t",
                )
            # indel on different segments
            # NOTE: shall we split the indel into multiple ones? No.
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
                print(
                    "".join(path),
                    a[i].st - off,
                    a[i].en - off,
                    y.qname,
                    y.mapq,
                    "+",
                    f"{info1};{info2}",
                    sep="\t",
                )
