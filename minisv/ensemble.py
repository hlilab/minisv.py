import os
import numpy as np
from .regex import re_info


def insilico_truth(msv_union):

    group_ids = []
    file_ids = []
    res_list = []
    with open(msv_union) as inf:
        for line in inf:
            elements = line.strip().split('\t')
            # start
            elements[1] = int(elements[1])
            # end
            elements[4] = int(elements[4])
            group_id = ''
            file_id = ''
            for (key, val) in re_info.findall(elements[5]):
                 if key == 'group_id':
                     group_id = val
                 if key == 'file_id':
                     file_id  = val

            assert group_id != ''
            assert file_id != ''

            if len(group_ids) >= 1:
                if group_id == group_ids[-1]:
                    res_list.append(elements)
                    file_ids.append(file_id)
                else:
                    # output
                    sorted_lines = sorted(res_list, key=lambda x: (int(x[1]), int(x[4])))
                    med_line = '\t'.join(map(str, sorted_lines[len(sorted_lines) >> 1]))
                    print(med_line + '\t' + ','.join(file_ids))
                    group_ids.append(group_id)
                    file_ids = [file_id]
                    res_list = [elements]
            else:
                file_ids.append(file_id)
                group_ids.append(group_id)
                res_list.append(elements)


def double_strand_break(collapsed_msv_union):
    """ Identify four cases of double strand breaks
    """
    breakpts = []
    with open(collapsed_msv_union) as inf:
        for line in inf:
            elements = line.strip().split('\t')
            breakpts.append([elements[0], int(elements[1]), elements[2][0], elements[-1], "left(+)"])
            breakpts.append([elements[3], int(elements[4]), elements[2][1], elements[-1], "right(-)"])

    breakpts.sort(key = lambda x: (x[0], int(x[1])))
    dtype_list = []

    for index, pt0 in enumerate(breakpts):
        for pt1 in breakpts[(index+1):]:
            dsbtype = ""
            if pt0[0] != pt1[0]:
                break
            if pt0[0] == pt1[0] and abs(pt1[1] - pt0[1]) >= 25e3:
                break

            svtype0 = [ val for (key, val) in re_info.findall(pt0[-2]) if key == 'SVTYPE' ][0]
            svtype1 = [ val for (key, val) in re_info.findall(pt1[-2]) if key == 'SVTYPE' ][0]

            asmrid0 = [ val for (key, val) in re_info.findall(pt0[-2]) if key == 'asm_read_ids' ][0].split(',')
            asmrid1 = [ val for (key, val) in re_info.findall(pt1[-2]) if key == 'asm_read_ids' ][0].split(',')

            # ['chr10', 51935447, '<', '1,2,3', 'SVTYPE=BND;SVLEN=0;count=13;group_id=41;file_id=2', 'left(+)'] ['chr10', 51935447, '<', '0', 'SVTYPE=BND;SVLEN=0;count=13;group_id=935;file_id=0', 'right(-)'] BND BND
            # skip the same breakpoint coordinate, do not consider the orientation
            if pt0[0] == pt1[0] and pt0[1] == pt1[1]:
                break

            # case 4 from the same read
            if svtype0 == svtype1 and svtype0 == "INV" and (len(list(set(asmrid0) & set(asmrid1))) >= 1):
                dsbtype = "Parallel_Breakpoints_case4"
            # case 4 from diff read
            if pt1[-1] == pt0[-1]:
                dsbtype = "Parallel_Breakpoints_case4"

            # case 1 from the same read
            if svtype0 == svtype1 and svtype0 == "DEL" and (len(list(set(asmrid0) & set(asmrid1))) >= 1):
                dsbtype = "case1"
            # case 1 from the diff read
            if pt0[-1] == 'right(-)' and pt1[-1] == 'left(+)':
                dsbtype = "case1"
            # case 1 from the diff read
            if pt0[-1] == 'left(+)' and pt1[-1] == 'right(-)' and (svtype0 == "BND" or svtype1 == "BND"):
                dsbtype = "case1"

            # from same read or different reads
            if svtype0 == svtype1 and svtype0 == "DUP":
                dsbtype = "case2"

            if svtype0 == svtype1 and svtype0 == "INS" and (len(list(set(asmrid0) & set(asmrid1))) >= 1):
                dsbtype = "Insertions_case3"

            print('\t'.join(map(str, pt0)), '\t'.join(map(str, pt1)), svtype0, svtype1, dsbtype)
            dtype_list.append(dsbtype)
    val, cnt = np.unique(dtype_list, return_counts=True)
    print('#stat')
    print('\t'.join(list(val)))
    print('\t'.join(list(map(str, cnt))))
