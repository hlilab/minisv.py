"""TODO: Turn the graph genome conversion into here
"""

from dataclasses import dataclass


@dataclass
class path2seg:
    seg: int
    pos: int


def path2ctg(seg, path_off, is_end):
    """INDEL path segment coordinate calculation"""
    b = []
    # path_off length is equal to indel number
    k = 0
    for i in range(len(path_off)):
        # segment k relative end <= indel relative start
        # NOTE: it is the same for both ends. Ask this week
        if is_end:
            while k < len(seg) and seg[k].path_en < path_off[i]:
                k += 1
        else:
            # NOTE: shall we use seg[k].path_st here?
            while k < len(seg) and seg[k].path_en <= path_off[i]:
                k += 1

        if k == len(seg):
            raise Exception("failed to find start position")
        # relative distance between indel start and segment start
        start_l = path_off[i] - seg[k].path_st
        # graph path > or linear genome
        if seg[k].strand > 0:
            # [seg index, absolute path start + distance to indel]
            b.append(path2seg(seg=k, pos=seg[k].ctg_st + start_l))
        else:
            # [seg index, absolute path end - distance to indel]
            # NOTE: why not st + (x+len - tst)?
            b.append(path2seg(seg=k, pos=seg[k].ctg_en - start_l))
    return b
