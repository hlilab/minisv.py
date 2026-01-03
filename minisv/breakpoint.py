import re
from dataclasses import dataclass

from .annotation import cal_cen_dist, cal_cen_overlap
from .regex import path_seg_pattern


def get_breakpoint(opt, z, file_handler=None):
    """
    opt: option dataclasses
    z: a list of reads in PAF/GAF/SAM
    """
    if len(z) < 2:
        return

    # sort by start position on the read
    z.sort(key=lambda x: x.qst)

    # filter short alignment towards the end of the read
    zen = len(z)
    for j in range(len(z) - 1, -1, -1):
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
        # l2, if l2 < 0 means a existing tsd
        qgap = y1.qst - y0.qen

        c0 = y0.coor[1]
        c1 = y1.coor[0]
        strand2 = "+"
        ori = c0.ori + c1.ori  # example: >>, ><

        # breakpoints are symmetric
        # chr1 400 >> chr1 500
        # chr1 500 << chr1 400 => chr1 400 >> chr1 500
        if not ((c0.ctg < c1.ctg) or (c0.ctg == c1.ctg and c0.pos < c1.pos)):
            # NOTE: why not c1.ori + c0.ori?
            ori = ("<" if c1.ori == ">" else ">") + ("<" if c0.ori == ">" else ">")
            c0 = y1.coor[0]
            c1 = y0.coor[1]
            # NOTE: why we still have <<? assign to negative strand.
            strand2 = "-"

        sv_info = infer_svtype(opt, c0, c1, ori, qgap)
        cen_str = ""
        if (c0.ctg in opt.cen) or (c1.ctg in opt.cen):
            dist0 = cal_cen_dist(opt, c0.ctg, c0.pos)
            dist1 = cal_cen_dist(opt, c1.ctg, c1.pos)
            cen_str = f";cen_dist={dist0 if dist0 < dist1 else dist1}"

            # NOTE: is this condition equal to same chromosome sv?
            #       we may need to add if c0.ctg == c1.ctg
            if sv_info.st >= 0 and sv_info.en >= sv_info.st:
                assert c0.ctg == c1.ctg, "sv contigs not the same"
                ov = cal_cen_overlap(opt, c0.ctg, sv_info.st, sv_info.en)
                cen_str += f";cen_overlap={ov}"

        # NOTE: do we have long inserted L1 from breakpoints?
        # NOTE: y0.qen for the qoff, do we need is_rev to pick between y0.qen and y0.qst?
        # visualize this part
        qoff_l = y0.qen if y0.qen < y1.qst else y1.qst
        qoff_r = y0.qen if y0.qen > y1.qst else y1.qst

        out = file_handler if file_handler is not None else None

        print(
            c0.ctg,
            c0.pos,
            ori,
            c1.ctg,
            c1.pos,
            y0.qname,
            y0.mapq if y0.mapq < y1.mapq else y1.mapq,
            strand2,
            f"{sv_info.str};qoff_l={qoff_l};qoff_r={qoff_r};qgap={qgap};mapq={y0.mapq},{y1.mapq};aln_len={y0.qen-y0.qst},{y1.qen-y1.qst}{cen_str};source={opt.name}",
            sep="\t",
            file=out
        )
    return None


@dataclass
class svtype:
    st: int = -1
    en: int = -1
    str: str = "SVTYPE=BND"


def infer_svtype(opt, c0, c1, ori, qgap):
    """
    inference of sv type from a pair of breakpoints

    # l1: bp length on the genome
    # l2: bp length on the query sequence
    # qgap: l2, l = l1-l2
    """
    if c0.ctg != c1.ctg:
        return svtype()

    l1 = c1.pos - c0.pos + 1

    if l1 < 0:
        # genome coordinate distance cannot be less than zero
        raise Exception("Error: l1 < 0 impossible")

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
    if ori == ">>" and l1 < qgap and qgap - l1 >= opt.min_len:
        return svtype(
            st=c0.pos,
            en=c1.pos + 1,
            str=f"SVTYPE=INS;SVLEN={qgap-l1};sv_region={c0.pos},{c1.pos+1}",
        )

    # long insertion with TSD
    # NOTE: only "-" strand has <<
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
    # NOTE: if qgap < 0, l = l1 - qgap?
    if ori == "<<" and qgap + l1 >= opt.min_len:
        if qgap < 0:
            st = c0.pos
        else:
            # NOTE: why c0.pos could < qgap?
            st = c0.pos - qgap if c0.pos > qgap else 0
        en = c1.pos + 1 if qgap < 0 else c1.pos + 1 + qgap
        return svtype(
            st=st, en=en, str=f"SVTYPE=DUP;SVLEN={qgap+l1};sv_region={st},{en}"
        )

    # NOTE: do we consider two reads? inversion
    if (ori == "<>" or ori == "><") and l1 >= opt.min_len:
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
    """get two end coordinate for one supp alignment"""
    # a path
    if re.match(r"^[><]", y.path):
        assert y.strand == "+", "reverse strand on path"

        x = 0
        for m in path_seg_pattern.findall(y.path):
            st = int(m[2])
            en = int(m[3])
            len = en - st

            # NOTE: visulize this compute
            if y.tst >= x and y.tst < x + len:
                r1 = break_end_coord(ctg=m[1], ori=m[0], pos=-1)
                r1.pos = (
                    st + (y.tst - x) if m[0] == ">" else st + (x + len - y.tst) - 1
                )
            if y.ten > x and y.ten <= x + len:
                r2 = break_end_coord(ctg=m[1], ori=m[0], pos=-1)
                r2.pos = st + (y.ten - x) - 1 if m[0] == ">" else st + (x + len - y.ten)
            x += len
    else:  # a contig
        r1 = break_end_coord(ctg=y.path, ori=">" if y.strand == "+" else "<", pos=-1)
        r1.pos = y.tst if y.strand == "+" else y.ten - 1
        r2 = break_end_coord(ctg=y.path, ori=">" if y.strand == "+" else "<", pos=-1)
        r2.pos = y.ten - 1 if y.strand == "+" else y.tst

    assert y.qen > y.qst, "query starting > end"
    r1.ql = r2.ql = y.qen - y.qst
    return [r1, r2]
