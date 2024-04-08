import re
from dataclasses import dataclass

path_seg_pattern = re.compile(r"([><])([^><:\s]+):(\d+)-(\d+)")


def get_breakpoint(opt, z):
    """ """

    if len(z) < 2:
        return

    # sort by start position on the read
    z.sort(key=lambda x: x.qst)

    # filter short alignment towards the end of the read
    zen = len(z)
    for j in range(len(z) - 1, -1):
        y = z[j]
        if y.qen - y.qst < opt.min_aln_len_end or y.mapq < opt.min_mapq_end:
            zen = j
        else:
            break

    if zen < 2:
        return
    # filter out short alignment towards the start of the read
    zst = 0
    for j in range(zen):
        y = z[j]
        if y.qen - y.qst < opt.min_aln_len_end or y.mapq < opt.min_mapq_end:
            zst = j + 1
        else:
            break

    if zen - zst < 2:
        return

    # construct the final alignment list
    zz = []
    for j in range(zst, zen):
        if z[j].qen - z[j].qst >= opt.min_aln_len_mid:
            zz.append(z[j])

    if len(zz) < 2:
        return

    # compute end coordinates of the breakpoints
    for j in range(len(zz)):
        r = get_end_coor(zz[j])
        zz[j].coor = r

    for j in range(1, len(zz)):
        y0 = zz[j - 1]
        y1 = zz[j]
        qgap = y1.qst - y0.qen
        # l2, if l2 < 0 means tsd

        c0 = y0.coor[1]
        c1 = y1.coor[0]
        strand2 = "+"
        ori = c0.ori + c1.ori  # example: >>, ><

        # breakpoints are symmetric
        # chr1 400 >> chr1 500
        # chr1 500 << chr1 400 => chr1 400 >> chr1 500
        if not ((c0.ctg < c1.ctg) or (c0.ctg == c1.ctg and c0.pos < c1.pos)):
            c0 = y1.coor[0]
            c1 = y0.coor[1]
            strand2 = "-"
            ori = ("<" if c1.ori == ">" else ">") + ("<" if c0.ori == ">" else ">")
            sv_info = infer_svtype(opt, c0, c1, ori, qgap)

            cen_str = ""
            if c0.ctg in opt.cen or c1.ctg in opt.cen:
                dist0 = cal_cen_dist(opt, c0.ctg, c0.pos)
                dist1 = cal_cen_dist(opt, c1.ctg, c1.pos)
                cen_str = f";cen_dist={dist0 if dist0 < dist1 else dist1}"

                # NOTE: is this condition equal to same chromosome sv?
                if sv_info.st >= 0 and sv_info.en >= sv_info.st:
                    ov = cal_cen_overlap(opt, c0.ctg, sv_info)
                    cen_str += f";cen_overlap={ov}"

            # NOTE: do we have long inserted L1 from breakpoints?
            print(
                c0.ctg,
                c0.pos,
                ori,
                c1.ctg,
                c1.pos,
                y0.qname,
                y0.mapq if y0.mapq < y1.mapq else y1.mapq,
                strand2,
                f"{sv_info.str};qgap={qgap};mapq={y0.mapq},{y1.mapq};aln_len={y0.qen-y0.qst},{y1.qen-y1.qst}{cen_str};source={opt.name}",
            )
    return None


def cal_cen_dist(opt, ctg, pos):
    if ctg not in opt.cen:
        return 1e9

    min = 1e9
    for i in range(len(opt.cen[ctg])):
        b = opt.cen[ctg][i]
        if pos < b[0]:
            d = b[0] - pos
        else:
            d = 0 if pos < b[1] else pos - b[1]
        min = min if min < d else d
    return min


def cal_cen_overlap(opt, ctg, st0, en0):
    """compute sv overlap with centromere"""
    if ctg not in opt.cen:
        return 0

    cov_st = 0
    cov_en = 0
    cov = 0
    for i in range(len(opt.cen[ctg])):
        b = opt.cen[ctg][i]

        if b[1] <= st0 or b[0] >= en0:  # not overlapping with [st0, en0)
            continue

        st1 = b[0] if b[0] > st0 else st0
        en1 = b[1] if b[1] < en0 else en0

        if st1 > cov_en:
            cov += cov_en - cov_st
            cov_st = st1
            cov_en = en1
        else:
            # NOTE: is it caused by overlapped centromeres
            cov_en = cov_en if cov_en > en1 else en1

    cov += cov_en - cov_st
    return cov


@dataclass
class svtype:
    st: int = -1
    en: int = -1
    str: str = "SVTYPE=BND"


def infer_svtype(opt, c0, c1, ori, qgap):
    """ """
    if c0.ctg != c1.ctg:
        return svtype()

    # l1: bp length
    l1 = c1.pos - c0.pos + 1

    if l1 < 0:
        raise Exception("Error: l1 < 0")

    # qgap: l2, l = l1-l2
    # long deletion from supplementary alignment
    if ori == ">>" and qgap < l1 and l1 - qgap >= opt.min_len:
        # NOTE: shall we use inner or outer part of the TSD for deletions?
        st = c0.pos + qgap if qgap < 0 else c0.pos
        en = c1.pos + 1 - qgap if qgap < 0 else c1.pos + 1
        return svtype(
            st=st,
            en=en,
            str=f"SVTYPE=DEL;SVLEN={-(l1-qgap)};sv_region={st},{en};tsd_len={-qgap if qgap < 0 else 0}",
        )

    # long insertion without TSD
    # NOTE: why not "<<"
    if ori == ">>" and l1 < qgap and qgap - l1 >= opt.min_len:
        # NOTE: why + 1
        return svtype(
            st=c0.pos,
            en=c1.pos + 1,
            str=f"SVTYPE=INS;SVLEN={qgap-l1};sv_region={c0.pos},{c1.pos+1}",
        )

    # long insertion with TSD
    # NOTE: why not ">>", only "-" strand?
    # l2 and l1 relative length condition? l2 > l1?
    if (
        ori == "<<"
        and qgap > 0
        and (l1 < c0.ql or l1 < c1.ql)
        and qgap + l1 >= opt.min_len
    ):
        # st: c1.pos? report TSD left end in reference coodinate?
        # en: c0.pos?
        return svtype(
            st=c0.pos,
            en=c1.pos + 1,
            str=f"SVTYPE=INS;SVLEN={qgap+l1};sv_region={c0.pos},{c1.pos+1};tsd_len={l1}",
        )

    # tandem duplication; similar to insertion with TSD
    # NOTE: why not ">>", only "-" strand?
    # l2 and l1 relative length condition?
    # coverage info?
    if ori == "<<" and qgap + l1 >= opt.min_len:
        if qgap < 0:
            st = c0.pos
        else:
            st = c0.pos - qgap if c0.pos > qgap else 0
        en = c1.pos + 1 if qgap < 0 else c1.pos + 1 + qgap
        return svtype(
            st=st, en=en, str=f"SVTYPE=DUP;SVLEN={qgap+l1};sv_region={st},{en}"
        )

    # NOTE: do we consider two reads?
    if ori == "<>" or ori == "><" and l1 >= opt.min_len:
        st = c0.pos + qgap if qgap < 0 else c0.pos
        en = c1.pos + 1 - qgap if qgap < 0 else c1.pos + 1
        return svtype(
            st=st, en=en, str=f"SVTYPE=INV;SVLEN={l1-qgap};sv_region={st},{en}"
        )

    return svtype()


@dataclass
class break_end_coord:
    ctg: str
    ori: str
    pos: str
    ql: int = -1


def get_end_coor(y):
    """ """
    if re.match(r"^[><]", y.path):  # a path
        if y.strand != "+":
            raise Exception("reverse strand on path")

        x = 0
        for m in path_seg_pattern.findall(y.path):
            st = int(m[2])
            en = int(m[3])
            len = en - st

            if y.tst >= x and y.tst < x + len:
                r1 = break_end_coord(ctg=m[1], ori=m[0], pos=-1)
                r1.pos = (
                    st + (y.tst - x) if m[0] == ">" else st + (x + len - y.tst) - 1
                )  # Why - 1 here
            if y.ten > x and y.ten <= x + len:
                r2 = break_end_coord(ctg=m[1], ori=m[0], pos=-1)
                r2.pos = (
                    st + (y.ten - x) - 1 if m[0] == ">" else st + (x + len - y.ten)
                )  # Why - 1 here
            x += len
    else:  # a contig
        r1 = break_end_coord(ctg=y.path, ori=">" if y.strand == "+" else "<", pos=-1)
        r1.pos = y.tst if y.strand == "+" else y.ten - 1
        r2 = break_end_coord(ctg=y.path, ori=">" if y.strand == "+" else "<", pos=-1)
        r2.pos = y.ten - 1 if y.strand == "+" else y.tst
    r1.ql = y.qen - y.qst
    r2.ql = y.qen - y.qst
    return [r1, r2]
