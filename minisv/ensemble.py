import os
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
