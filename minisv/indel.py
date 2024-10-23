import re
from dataclasses import dataclass

from .annotation import cal_cen_dist
from .graph_genome_coor import path2ctg
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

    qoff: int = 0
    qoff_l: int = 0
    qoff_r: int = 0


@dataclass
class path_segment:
    ctg: str
    ctg_st: int
    ctg_en: int
    strand: int
    path_st: int
    path_en: int


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

        is_rev = y.strand == "-"

        cg_segs = cigar_pattern.findall(y.cg)
        # record cigar indels
        x = y.tst
        q = 0
        a = []
        for length, op in cg_segs:
            length = int(length)
            if length >= opt.min_len:
                if op == "I":
                    qoff = y.qen - (q + length) if is_rev else q + y.qst
                    a.append(
                        # NOTE: to be more accurate in indels
                        indel_coord(
                            st=x,
                            en=x,
                            len=length,
                            indel_seq=".",
                            tsd_len=0,
                            tsd_seq=".",
                            polyA_len=0,
                            int_seq=".",
                            qoff=qoff,
                            qoff_l=qoff,
                            qoff_r=qoff + length,
                        )  # NOTE: difference of reverse q_off between ins and del?
                    )
                elif op == "D":
                    qoff = y.qen - q if is_rev else q + y.qst
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
                            qoff=qoff,
                            qoff_l=qoff,
                            qoff_r=qoff,
                        )
                    )
            if op in ["M", "D", "=", "X", "N"]:
                x += length
            # NOTE: no "N" in query offset?
            if op in ["M", "=", "X", "I", "S", "H"]:
                q += length

        # No indel cigar or too many indel
        # per 10k alleles cannot be too many
        # NOTE: considering ultralong and contig alignment
        if len(a) == 0 or len(a) > y.qlen * 1e-4 * opt.max_cnt_10k:
            continue

        # set stl/enl and str/enr
        for i in range(len(a)):
            a[i].stl = a[i].str = a[i].st
            a[i].enl = a[i].enr = a[i].en

        # parse ds:Z tag
        # for retrotransposon
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
                        if a[i].st != x or a[i].en != x or a[i].len != length:
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

                # NOTE: adjust TSD length
                a[i].stl = a[i].st - rlen
                a[i].enl = a[i].en - rlen
                # NOTE: for right point, why not choose the minimum point by subtracting llen
                #       Could CIGAR-based indel have l2 (qgap) < 0? Yes
                a[i].str = a[i].st + llen
                a[i].enr = a[i].en + llen

                # NOTE: adjust query length offset
                if is_rev:
                    a[i].qoff_l = a[i].qoff_l - llen
                    a[i].qoff_r = a[i].qoff_r + rlen
                else:
                    a[i].qoff_l = a[i].qoff_l - rlen
                    a[i].qoff_r = a[i].qoff_r + llen
        # end of y.ds

        if opt.dbg:
            print('X0', line)

        # genome graph path reference segments
        seg = []
        # with <>: this is a path
        if re.match(r"[><]", y.path):
            x = 0
            assert y.strand == "+", "reverse strand on path"
            for m in path_seg_pattern.findall(y.path):
                st, en = int(m[2]), int(m[3])
                # [path, absolute path start, absolute path end, path orientation, relative seg start, relative seg end]
                strand = 1 if m[0] == ">" else -1
                seg.append(
                    path_segment(
                        ctg=m[1],
                        ctg_st=st,
                        ctg_en=en,
                        strand=strand,
                        path_st=x,
                        path_en=x + (en - st),
                    )
                )
                # accumulate seg start
                x += en - st
        else:
            # https://github.com/lh3/miniasm/blob/master/PAF.md
            # https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
            # [path, 0, path length, 1, 0, path length]
            # NOTE: linear coordinate always strand == "+"
            seg.append(
                path_segment(
                    ctg=y.path,
                    ctg_st=0,
                    ctg_en=y.tlen,
                    strand=1,
                    path_st=0,
                    path_en=y.tlen,
                )
            )

        # starts and ends for the indels on the segments
        # also record start and end segments index
        # these are used to overlap between path segments and indels
        # left
        off_stl = []
        off_enl = []
        # right
        off_str = []
        off_enr = []
        for i in range(len(a)):
            off_stl.append(a[i].stl)
            off_enl.append(a[i].enl)
            off_str.append(a[i].str)
            off_enr.append(a[i].enr)

        global_qname = y.qname
        stl = path2ctg(seg, off_stl, False)
        enl = path2ctg(seg, off_enl, True)
        str = path2ctg(seg, off_str, False)
        enr = path2ctg(seg, off_enr, True)

        for i in range(len(a)):
            if not (
                stl[i].seg == str[i].seg
                and stl[i].seg == enl[i].seg
                and str[i].seg == enr[i].seg
            ):  # skip all indels st and end not on the same segment
                continue

            # find the corresponding segment
            s = seg[stl[i].seg]
            st = stl[i].pos
            en = enl[i].pos
            strand = y.strand

            # graph path reverse direction
            if s.strand < 0:
                # NOTE: this is a new update that for reverse complement tsd
                #       reverse complement sequence if strand is reverse
                # NOTE: How about linear genome - strand??
                a[i].polyA_len = -a[i].polyA_len
                a[i].tsd_seq = mg_revcomp(a[i].tsd_seq)
                a[i].int_seq = mg_revcomp(a[i].int_seq)
                # NOTE: why reverse start and end as well??
                st = enr[i].pos
                en = str[i].pos
                strand = "-" if strand == "+" else "+"

            info1 = ("SVTYPE=INS" if a[i].len > 0 else "SVTYPE=DEL") + (
                f";SVLEN={a[i].len};qoff_l={a[i].qoff_l};qoff_r={a[i].qoff_r};tsd_len={a[i].tsd_len};polyA_len={a[i].polyA_len}"
            )
            info2 = f"source={opt.name};tsd_seq={a[i].tsd_seq if len(a[i].tsd_seq)>0 else '.'};insert={a[i].int_seq if len(a[i].int_seq)>0 else '.'}"

            if s.ctg in opt.cen:
                dist_st = cal_cen_dist(opt, s.ctg, st)
                dist_en = cal_cen_dist(opt, s.ctg, en)
                info1 += f";cen_dist={dist_st if dist_st < dist_en else dist_en}"

            # 7 columns for indels
            print(
                s.ctg,
                st,
                en,
                global_qname,
                y.mapq,
                strand,
                f"{info1};{info2}",
                sep="\t",
            )
            # indel on different segments
            # NOTE: shall we split the indel into multiple ones? No.
            # else:
            #     path = []
            #     length = 0
            #     # iterate through segments
            #     # j from indel start seg index to ends seg index
            #     for j in range(sts[i][0], ens[i][0] + 1):
            #         s = seg[j]
            #         length += s[2] - s[1]
            #         orientation = ">" if s[3] > 0 else "<"
            #         # >chr1:1-2
            #         path.append(orientation + s[0] + f":{s[1]}-{s[2]}")

            #     # relative seg start
            #     off = seg[sts[i][0]][4]
            #     # minigraph always use "+" in genome path
            #     # convert back to relative path coordinate
            #     # [chrom, start, end, read_name, mapq, strand, indel length, tsd length, polyA length, indel seq]
            #     print(
            #         "".join(path),
            #         a[i].st - off,
            #         a[i].en - off,
            #         y.qname,
            #         y.mapq,
            #         "+",
            #         f"{info1};{info2}",
            #         sep="\t",
            #     )
