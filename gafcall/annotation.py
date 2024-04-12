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
