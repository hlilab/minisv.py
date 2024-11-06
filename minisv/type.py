from dataclasses import dataclass
import re


@dataclass
class svisec:
    type: str = ""
    len: int = 0
    st: int = 0
    en: int = 0
    flag: int = 0


def get_type(t, col_info):
    info = t[col_info]
    m = re.findall(r"\b(SVTYPE|SVLEN|qoff_l|qoff_r)=([^\s;]+)", info)

    type = None
    qoff_l = -1
    qoff_r = -1
    flag = 0
    len = 0

    for element_type, element_val in m:
        if element_type == "SVTYPE":
            type = element_val
        elif element_type == "SVLEN":
            len = int(element_val)
        elif element_type == "qoff_l":
            qoff_l = int(element_val)
        elif element_type == "qoff_r":
            qoff_r = int(element_val)

    if type is None or qoff_l < 0 or qoff_r < 0:
        raise Exception("missing information")
    # print(type, len, qoff_l, qoff_r)

    # insertion
    if type == "INS" or type == "DUP":
        flag = 1
    # deletion
    elif type == "DEL":
        flag = 2
    # inversion
    elif type == "INV":
        flag = 4
    # translocation
    elif type == "BND" and col_info == 8 and t[0] != t[3]:
        flag = 8
    return svisec(type = type, len = len, st = qoff_l, en = qoff_r, flag = flag)


def simple_type(type, ctg1, ctg2):
    flag = 0

    # insertion
    if type == "INS" or type == "DUP":
        flag = 1
    # deletion
    elif type == "DEL":
        flag = 2
    # inversion
    elif type == "INV":
        flag = 4
    # translocation
    elif type == "BND" and ctg1 != ctg2:
        flag = 8
    return flag

