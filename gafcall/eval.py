import math
import re
from dataclasses import dataclass
from typing import Optional

from .merge import svinfo


def eval(base, test, opt):
    """Evaluate the performance of the test SVs against the base SVs
    Args:
        base: list of base SVs
        test: list of test SVs
        opt: EvalConfig object
    Returns:
        None
    """
    base = gc_parse_sv(opt, base)
    test = gc_parse_sv(opt, test)
    tot_fn, fn = gc_cmp_sv(opt, test, base, "FN")
    tot_fp, fp = gc_cmp_sv(opt, base, test, "FP")
    print("RN", tot_fn, fn, (fn / tot_fn))
    print("RP", tot_fp, fp, (fp / tot_fp))


def gc_parse_sv(opt, file_path):
    min_len = math.floor(opt.min_len * opt.read_len_ratio + 0.499)
    sv = []
    ignore_id = {}
    with open(file_path) as f:
        for line in f:
            if line[0] == "#":
                continue
            t = line.strip().split("\t")

            if re.match(r"^\d+", t[1]) is None:
                continue

            # start pos
            t[1] = int(t[1])
            type = 0
            info = None

            if re.match(r"^[><][><]$", t[2]):
                # breakpoint type
                type = 3
                info = t[8]
            elif re.search(r";", t[7]):
                # VCF
                type = 1
                info = t[7]
            elif re.match(r"^\d+$", t[2]) and re.search(r";", t[6]):
                # NOTE: BED for CIGAR INDEL?
                type = 2
                info = t[6]
            # else:
            #    raise Exception('type == 0')

            if type == 0:
                continue

            svtype = None
            svlen = 0
            m = re.findall(r"SVTYPE=([^\s;]+)", info)
            if len(m) > 0:
                svtype = m[0]
            m = re.findall(r"SVLEN=([^\s;]+)", info)
            if len(m) > 0:
                svlen = int(m[0])

            # BED line
            if type == 2:
                t[2] = int(t[2])
                if t[1] > t[2]:
                    raise Exception("incorrect BED format")
                sv.append(
                    svinfo(
                        ctg=t[0],
                        pos=t[1],
                        ctg2=t[0],
                        pos2=t[2],
                        ori=">>",
                        SVTYPE=svtype,
                        SVLEN=svlen,
                    )
                )
            elif type == 3:
                # breakpoint line
                t[4] = int(t[4])
                sv.append(
                    svinfo(
                        ctg=t[0],
                        pos=t[1],
                        ctg2=t[3],
                        pos2=t[4],
                        ori=t[2],
                        SVTYPE=svtype,
                        SVLEN=svlen,
                    )
                )
            elif type == 1:
                # VCF line
                # reference allele
                rlen = len(t[3])
                en = t[1] + rlen - 1

                # initialize a sv info object
                s = svinfo(ctg=t[0], pos=t[1] - 1, ctg2=t[0], pos2=en, ori=">>")
                # indel
                if re.match(r"^[A-Z,]+$", t[4]):
                    # assume full allele sequence; override SVTYPE/SVLEN even if present
                    # multiple allele
                    alt = t[4].split(",")
                    for i in range(len(alt)):
                        a = alt[i]
                        # alt allele - reference allele
                        length = len(a) - rlen

                        if abs(length) < min_len:
                            continue

                        if length < 0:
                            # deletion
                            sv.append(
                                svinfo(
                                    ctg=s.ctg,
                                    pos=s.pos,
                                    ctg2=s.ctg,
                                    pos2=en,
                                    ori=">>",
                                    SVTYPE="DEL",
                                    SVLEN=length,
                                )
                            )
                        else:
                            # insertion
                            sv.append(
                                svinfo(
                                    ctg=s.ctg,
                                    pos=s.pos,
                                    ctg2=s.ctg,
                                    pos2=en,
                                    ori=">>",
                                    SVTYPE="INS",
                                    SVLEN=length,
                                )
                            )
                else:
                    if t[2] != ".":
                        # ignore previously visited ID
                        if t[2] in ignore_id:
                            continue
                        ignore_id[t[2]] = 1

                    m = re.findall(r"\bMATEID=(\d+)", info)
                    if len(m) > 0:
                        ignore_id[m[0]] = 1

                    if svtype is None:
                        # we don't infer SVTYPE from breakpoint for existing data
                        raise Exception("cannot determine SVTYPE")

                    s.SVTYPE = svtype

                    # NOTE: js might have a bug here
                    if svtype != "BND" and abs(svlen) < min_len:
                        continue  # too short

                    # ENCODE Deletion SVLEN consistently
                    if svtype == "DEL" and svlen > 0:
                        svlen = -svlen

                    s.SVLEN = svlen

                    m = re.findall(r"\bEND=([0-9]+)", info)
                    if len(m) > 0:
                        s.pos2 = int(m[0])  # NOTE: - 1
                    elif rlen == 1:
                        # ignore one-sided breakpoint
                        if opt.dbg:
                            print(t[4])
                        # NOTE: how this ignore one-sided breakpoint?
                        #       why < 6
                        if svtype == "BND" and len(t[4]) < 6:
                            continue
                        # NOTE: correct breakpoint end position?
                        if svtype == "DEL" or svtype == "DUP" or svtype == "INV":
                            s.pos2 = s.pos + abs(svlen)

                    # match VCF alt allele
                    # NOTE: ^\s in js?
                    m = re.findall(r"^[A-Z]+\[([^s:]+):(\d+)\[$", t[4])
                    # NOTE: ^\s in js?
                    m1 = re.findall(r"^\]([^s:]+):(\d+)\][A-Z]+$", t[4])
                    # [p[t: reverse comp piece extending right of p is joined before t
                    m2 = re.findall(r"^\[([^s:]+):(\d+)\[[A-Z]+$", t[4])
                    # t]p]: reverse comp piece extending left of p is joined after t
                    m3 = re.findall(r"^[A-Z]+\]([^s:]+):(\d+)\]$", t[4])
                    if len(m) > 0:
                        s.ctg2 = m[0]
                        s.pos2 = int(m[1])
                        s.ori = ">>"
                    elif len(m1) > 0:
                        s.ctg2 = m[0]
                        s.pos2 = int(m[1])
                        s.ori = "<<"
                    elif len(m2) > 0:
                        s.ctg2 = m[0]
                        s.pos2 = int(m[1])
                        # NOTE: why the direction? should be symmetric
                        #       maybe <>
                        s.ori = "><"
                    elif len(m3) > 0:
                        s.ctg2 = m[0]
                        s.pos2 = int(m[1])
                        # NOTE: why the direction?
                        #       maybe ><
                        s.ori = "<>"

                    if svtype != "BND" and s.ctg != s.ctg2:
                        raise Exception("different contigs for non-BND type")

                    if s.ctg == s.ctg2 and s.pos > s.pos2:
                        tmp = s.pos
                        s.pos = s.pos2
                        s.pos2 = tmp
                    sv.append(s)
    if opt.dbg:
        print(f"parsed {len(sv)} SVs")
        for i in range(len(sv)):
            s = sv[i]
            print(s.ctg, s.pos, s.ori, s.ctg2, s.pos2, s.SVTYPE, s.SVLEN, sep="\t")

    return sv


def gc_cmp_sv(opt, base, test, label):
    """Compare two lists of SVs
    Args:
        base: list of base SVs
        test: list of test SVs
        opt: EvalConfig object
    Returns:
        tuple of (total_fn, fn)
    """
    h = {}
    for i in range(len(base)):
        s = base[i]
        if s.ctg not in h:
            h[s.ctg] = []
        if s.ctg2 not in h:
            h[s.ctg2] = []
        h[s.ctg].append({"st": s.pos, "en": s.pos + 1, "data": s})
        h[s.ctg2].append({"st": s.pos2, "en": s.pos2 + 1, "data": s})

    for ctg in h:
        h[ctg] = iit_sort_copy(h[ctg])
        iit_index(h[ctg])

    tot = error = 0
    for j in range(len(test)):
        t = test[j]

        # Not long enough
        if t.SVTYPE == "BND" and abs(t.SVLNE) < opt.min_len:
            continue

        tot += 1
        n = eval1(opt, h, t.ctg, t.pos, t) + eval1(opt, h, t.ctg2, t.pos2, t)
        if n == 0:
            error += 1
            if opt.print_err:
                print(label, t.ctg, t.pos, t.ori, t.ctg2, t.pos2, t.SVTYPE, t.SVLEN)

    return [tot, error]


@dataclass
class iit_obj:
    st: int
    en: int
    max: int = 0
    data: Optional[svinfo] = None


# interval query
def iit_sort_copy(a):
    """Sort a list of SVs by start position and return a copy
    NOTE: what is iit sort? is it insertion sort?
    """
    a.sort(key=lambda x: x["st"])
    b = []
    for i in range(len(a)):
        b.append(iit_obj(st=a[i]["st"], en=a[i]["en"], max=0, data=a[i]["data"]))
    return b


def iit_index(a):
    """
    NOTE: what is iit index for?
    """
    if len(a) == 0:
        return -1

    for i in range(0, len(a), 2):
        last = a[i].en
        a[i].max = a[i].en
        last_i = i

    # NOTE: what is this loop for?
    k = 1
    while 1 << k <= len(a):
        i0 = (1 << k) - 1
        step = 1 << (k + 1)
        x = 1 << (k - 1)
        for i in range(i0, len(a), step):
            a[i].max = a[i].en
            if a[i].max < a[i - x].max:
                a[i].max = a[i - x].max
            e = a[i + x].max if i + x < len(a) else last
            if a[i].max < e:
                a[i].max = e
        last_i = last_i - x if (last_i >> k) & 1 else last_i + x
        if last_i < len(a):
            last = last if last > a[last_i].max else a[last_i].max
        k += 1
    return k - 1


def eval1(opt, h, ctg, pos, t):
    if ctg not in h:
        return False

    st = pos - opt.win_size if pos > opt.win_size else 0
    en = pos + opt.win_size
    a = iit_overlap(h[ctg], st, en)
    n = 0
    for i in range(len(a)):
        if same_sv1(opt, a[i].data, t):
            n += 1
    return n


def iit_overlap(a, st, en):
    """ """
    h = 0
    stack = []
    b = []

    h = 0
    while 1 << h <= len(a):
        h += 1

    h -= 1
    stack.append([(1 << h) - 1, h, 0])

    while len(stack) > 0:
        t = stack.pop()
        # NOTE: what does x, h, w represent?
        x = t[0]
        h = t[1]
        w = t[2]
        # NOTE: why this should be 3
        if h <= 3:
            # NOTE: why right shift then left shift bit
            i0 = x >> h << h
            i1 = i0 + (1 << (h + 1)) - 1
            if i1 >= len(a):
                i1 = len(a)
            i = i0
            while i < i1 and a[i].st < en:
                if st < a[i].en:
                    b.append(a[i])
                i += 1
        # NOTE: why this should be 0
        elif w == 0:
            stack.append([x, h, 1])
            y = x - (1 << (h - 1))
            if y >= len(a) or a[y].max > st:
                stack.append([y, h - 1, 0])
        elif x < len(a) and a[x].st < en:
            if st < a[x].en:
                b.append(a[x])
            stack.append([x + (1 << (h - 1)), h - 1, 0])
    return b


def same_sv1(opt, b, t):
    # check type
    if b.SVTYPE != t.SVTYPE:
        if (not (b.SVTYPE == "DUP" and t.SVTYPE == "INS")) and (
            not (b.SVTYPE == "INS" and t.SVTYPE == "DUP")
            and b.SVTYPE != "BND"
            and t.SVTYPE != "BND"
        ):  # special case for INS vs DUP
            return False

    # check length
    # NOTE: length should be similar to each other
    len_check = (
        abs(b.SVLEN) >= abs(t.SVLEN) * opt.min_len_ratio
        and abs(t.SVLEN) >= abs(b.SVLEN) * opt.min_len_ratio
    )

    if b.SVTYPE != "BND" and t.SVTYPE != "BND" and not len_check:
        return False

    match1 = match2 = 0
    # check the coordinates of end points within window
    # NOTE: what does match1 and match2 denote? two breakpoints?
    if (
        t.ctg == b.ctg
        and t.pos >= b.pos - opt.win_size
        and t.pos <= b.pos + opt.win_size
    ):
        match1 |= 1

    if (
        t.ctg == b.ctg2
        and t.pos >= b.pos2 - opt.win_size
        and t.pos <= b.pos2 + opt.win_size
    ):
        match1 |= 2

    if (
        t.ctg2 == b.ctg
        and t.pos2 >= b.pos - opt.win_size
        and t.pos2 <= b.pos + opt.win_size
    ):
        match2 |= 1
    if (
        t.ctg2 == b.ctg2
        and t.pos2 >= b.pos2 - opt.win_size
        and t.pos2 <= b.pos2 + opt.win_size
    ):
        match2 |= 2

    # NOTE: duplicated condition? seems to be a bug
    if b.SVTYPE == "DUP" and t.SVTYPE == "INS":
        return (match1 & 1) != 0
    elif b.SVTYPE == "INS" and t.SVTYPE == "DUP":
        return (match1 & 1) != 0
    elif b.SVTYPE == "BND" or t.SVTYPE == "BND":
        return ((match1 & 1) != 0 and (match2 & 2) != 0) or (
            (match1 & 2) != 0 and (match2 & 1) != 0
        )
    else:
        return (match1 & 1) != 0 and (match2 & 2) != 0
