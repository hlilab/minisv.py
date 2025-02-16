import math
import re
from dataclasses import dataclass
from typing import Optional, List

from .regex import re_info, sv_region_regex


def normal_round(n):
    """replace original round which render error when
    encountering x.5"""
    if n - math.floor(n) < 0.5:
        return math.floor(n)
    return math.ceil(n)


@dataclass
class svinfo:
    ctg: str
    pos: int
    st: int = -1
    en: int = -1
    ori: Optional[str] = None
    ctg2: Optional[str] = None
    pos2: Optional[int] = None
    _mapq: int = 0
    strand: Optional[str] = None
    is_bp: bool = False
    info: Optional[str] = None
    SVTYPE: Optional[str] = ""
    source: Optional[str] = None
    SVLEN: Optional[int] = None
    cen_dist: Optional[int] = None
    cen_overlap: Optional[int] = None
    tsd_len: Optional[int] = None
    polyA_len: Optional[int] = None
    sv_region: Optional[str] = None
    name: Optional[str] = None
    vaf: Optional[float] = None
    inv: Optional[bool] = False
    count: Optional[int] = 0
    svid: Optional[str] = ""
    clean_svid: Optional[str] = ""
    readids: Optional[List[str]] = None
    asmreadids: Optional[List[str]] = None


@dataclass
class svobj:
    ctg: str
    pos_max: int
    SVTYPE: str
    v: list
    is_bp: bool = False


# NOTE: do we merge this with the eval function
# only for minisv gsv file
def parse_sv(t):
    v = svinfo(
        ctg=t[0],
        pos=int(t[1]),
        st=-1,
        en=-1,
        ori=None,
        ctg2=None,
        pos2=None,
        _mapq=0,
        strand=None,
        is_bp=False,
        name=None,
        info=None,
    )

    if re.match(r"[><]", t[2]):
        v.is_bp = True

    off = 6 if v.is_bp else 4
    v._mapq = int(t[off])
    v.strand = t[off + 1]
    v.info = t[off + 2]
    v.name = t[off - 1]

    for m in re_info.findall(v.info):
        setattr(v, m[0], m[1])

    if v.SVTYPE is None or v.source is None:
        raise Exception("missing SVTYPE or source")

    if v.SVLEN is not None:
        v.SVLEN = int(v.SVLEN)
    if v.cen_dist is not None:
        v.cen_dist = int(v.cen_dist)
    if v.cen_overlap is not None:
        v.cen_overlap = int(v.cen_overlap)

    if not v.is_bp:  # indel
        v.st = v.pos
        v.en = int(t[2])
    else:  # breakpoints
        v.ori = t[2]
        v.ctg2 = t[3]
        v.pos2 = int(t[4])
        if v.sv_region is not None:
            m = sv_region_regex.findall(v.sv_region)
            if len(m) > 0:
                v.st = int(m[0][0])
                v.en = int(m[0][1])
    return v


def splitmix32(a):
    def generate_random():
        nonlocal a
        a |= 0
        a = (a + 0x9E3779B9) & 0xFFFFFFFF
        t = a ^ (a >> 16)
        t = (t * 0x21F0AAAD) & 0xFFFFFFFF
        t = t ^ (t >> 15)
        t = (t * 0x735A2D97) & 0xFFFFFFFF
        t = t ^ (t >> 15)
        t = t >> 0
        return t / 4294967296.0

    return generate_random


def merge_sv(opt, input):
    # read lines from stdin
    sv = []
    rng = splitmix32(11)
    for line in input:
        t = line.strip().split("\t")
        v = parse_sv(t)

        while len(sv) > 0:
            # NOTE: so we don't merge SV with too many alleles
            if (
                sv[0].ctg != v.ctg
                or v.pos - sv[0].pos_max > opt.win_size
                or len(sv) > opt.max_allele
            ):
                out = sv.pop(0)
                write_sv(opt, out.v)
            else:
                break

        cnt_same = []
        for i in range(len(sv)):
            c = 0
            if sv[i].SVTYPE == v.SVTYPE or (
                sv[i].SVTYPE == "INS" and v.SVTYPE == "DUP"
            ):
                if len(sv[i].v) <= opt.max_check:
                    for j in range(len(sv[i].v)):
                        if same_sv(opt, sv[i].v[j], v):
                            c += 1
                else:
                    # Reservoir sampling for a subset of reads to reduce time
                    p = [-1] * len(sv[i].v)
                    for j in range(len(sv[i].v)):
                        k = j if j < opt.max_check else math.floor(j * rng())

                        if k < opt.max_check:
                            p[k] = j

                    for k in range(opt.max_check):
                        if same_sv(opt, sv[i].v[p[k]], v):
                            c += 1
                    c = math.floor(c / opt.max_check * len(sv[i].v) + 0.499)
            cnt_same.append(c)

        max = 0
        max_i = -1
        for i in range(len(sv)):
            if cnt_same[i] > max:
                max = cnt_same[i]
                max_i = i

        # NOTE: do we merge cross-chromosome svs?
        # yes, we select the median svs using >>1 index
        if max > 0 and max_i >= 0:
            sv[max_i].v.append(v)
            sv[max_i].pos_max = (
                sv[max_i].pos_max if sv[max_i].pos_max > v.pos else v.pos
            )
        else:
            sv.append(
                svobj(
                    ctg=v.ctg,
                    pos_max=v.pos,
                    SVTYPE=v.SVTYPE if v.SVTYPE != "DUP" else "INS",
                    is_bp=v.is_bp,
                    v=[v],
                )
            )

    while len(sv) > 0:
        out = sv.pop(0)
        write_sv(opt, out.v)
    return None


# NOTE: do we merge this with the eval function
def same_sv(opt, v, w):
    # not the same (breakpoint or indel type)
    if v.is_bp != w.is_bp:
        return False

    # not on the same contig
    # NOTE: do we merge cross-chrom svs? yes
    if v.ctg != w.ctg:
        return False

    # not the same type
    if v.SVTYPE != w.SVTYPE:
        return False

    # test inversions
    if v.is_bp and w.is_bp and v.ori != w.ori:
        # NOTE: clear
        if not ((v.ori == "><" and w.ori == "<>") or (v.ori == "<>" and w.ori == "><")):
            return False

    # pos differ too much
    if abs(v.pos - w.pos) > opt.win_size:
        return False

    if not v.is_bp:  # indel
        # check end position
        if abs(v.en - w.en) > opt.win_size:
            return False
    else:
        # breakpoint
        if v.ctg2 != w.ctg2:
            return False
        if abs(v.pos2 - w.pos2) > opt.win_size:
            return False

    # SVLEN differences
    if v.SVLEN is not None and w.SVLEN is not None:
        if v.SVLEN * w.SVLEN <= 0:  # redundant but doesnot hurt to check
            return False

        vl = abs(v.SVLEN)
        wl = abs(w.SVLEN)
        # NOTE: is this testing sequence identities?
        if abs(vl - wl) > 0.5 * (vl + wl) * opt.max_diff:
            return False

    return True


def write_sv(opt, s):
    """Selecting the median indexed SV"""
    if len(s) == 0:
        return

    # NOTE: this is for selecting the median position sv from a set of overlapped svs
    v = s[len(s) >> 1]

    # filter by centromere distance
    if opt.min_cen_dist > 0:
        if v.cen_overlap is not None and v.cen_overlap > 0:
            return
        if v.cen_dist is not None and v.cen_dist <= opt.min_cen_dist:
            return

    rt_len_arr = []
    rt_len = 0
    for i in range(len(s)):
        if s[i].tsd_len is not None and s[i].polyA_len is not None:
            # NOTE: why choose minimum length between TSD and polyA
            #       instead of filtering either one separately
            rt_len_arr.append(
                int(s[i].tsd_len)
                if int(s[i].tsd_len) < abs(int(s[i].polyA_len))
                else abs(int(s[i].polyA_len))
            )

    if len(rt_len_arr) > 0:
        # NOTE: choose the medium pos RT length 
        rt_len = rt_len_arr[len(rt_len_arr) >> 1]

    # count
    mapq = 0
    cnt = {}
    cnt_strand = [0, 0]
    name = []
    cnt_fr = 0
    cnt_rf = 0
    for i in range(len(s)):
        mapq += s[i]._mapq
        if s[i].source not in cnt:
            # NOTE: count the number of sample for each sv
            cnt[s[i].source] = [0, 0]
        j = 0 if s[i].strand == "+" else 1
        cnt[s[i].source][j] += 1
        cnt_strand[j] += 1
        name.append(s[i].name)
         
        if s[i].ori == "><":
            cnt_fr += 1
        elif s[i].ori == "<>":
            cnt_rf += 1

    mapq = normal_round(mapq / len(s))
    # filter by the count
    if opt.min_rt_len > 0 and rt_len >= opt.min_rt_len:
        if len(s) < opt.min_cnt_rt:
            return
    else:
        if len(s) < opt.min_cnt:
            return
        # NOTE: strand bias?
        if cnt_strand[0] < opt.min_cnt_strand or cnt_strand[1] < opt.min_cnt_strand:
            return

    cnt_arr = []
    for src in cnt:
        cnt_arr.append(f"{src}:{cnt[src][0]},{cnt[src][1]}")

    info = f"avg_mapq={mapq:.0f};"
    info += f"count={'|'.join(cnt_arr)};"
    info += f"rt_len={rt_len};"
    info += re.sub(r"(;?)source=[^;\s=]+", "", v.info)

    if cnt_fr + cnt_rf > 0:
        info += f";count_fr={cnt_fr};count_rf={cnt_rf}"
        if cnt_fr * cnt_rf == 0 and v.ctg == v.ctg2:
            info += ";foldback"

    info += f";reads={','.join(name)}"

    if not v.is_bp:
        print(v.ctg, v.st, v.en, ".", len(s), v.strand, info, sep="\t")
    else:
        print(
            v.ctg, v.pos, v.ori, v.ctg2, v.pos2, ".", len(s), v.strand, info, sep="\t"
        )
