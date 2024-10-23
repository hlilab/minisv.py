from dataclasses import dataclass


@dataclass
class path2seg:
    seg: int
    pos: int


def path2ctg(seg, path_off, is_end):
    """INDEL path segment coordinate calculation

    seg: list, graph path segments
    path_off: list, indels query offsets
    is_end: bool, path_off contains end points or not

    """
    b = []
    seg_num = len(seg)
    # path_off length is equal to indel number
    k = 0
    for i in range(len(path_off)):
        # segment k relative end <= indel relative start/end
        if is_end:
            while k < seg_num and seg[k].path_en < path_off[i]:
                k += 1
        else:
            while k < seg_num and seg[k].path_en <= path_off[i]:
                k += 1

        if k == seg_num:
            raise Exception("failed to convert path offset to contig offset")

        # relative distance between indel start and segment start
        start_l = path_off[i] - seg[k].path_st
        # graph path > or linear genome
        if seg[k].strand > 0:
            # [seg index, absolute path start + distance to indel]
            b.append(path2seg(seg=k, pos=seg[k].ctg_st + start_l))
        else:
            # [seg index, absolute path end - distance to indel]
            b.append(path2seg(seg=k, pos=seg[k].ctg_en - start_l))
    return b
