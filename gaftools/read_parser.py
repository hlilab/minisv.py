from dataclasses import dataclass
from typing import Optional

from .breakpoint import break_end_coord, get_breakpoint
from .indel import get_indel


@dataclass
class alignment:
    qname: str
    mapq: int
    qst: int
    qen: int
    qlen: int
    tlen: int
    tst: int
    ten: int
    path: str
    strand: str
    cg: Optional[str] = None
    ds: Optional[str] = None
    coor: Optional[break_end_coord] = None


def load_reads(input_file, opt):
    """
    Load GAF/PAF/SAM file, filter based on mapping quality and minimum map length,
    and sort the lines by cluster location in read.

    Args:
    - input_file (str): an input file
    - opt (dataclass): a object containing all configurations

    Returns:
    - sorted_lines (list): List of PAF/GAF/SAM lines sorted by cluster location in read.

    NOTE: This function is a generator. The PAF/GAF file already outputs sorted reads by names.
    """
    z = []
    with open(input_file, "r") as input_file_handler:
        for line in input_file_handler:
            t = line.strip().split("\t")
            if len(t) < 11:
                continue

            if len(z) > 0 and t[0] != z[0].qname:
                get_indel(opt, z)
                get_breakpoint(opt, z)
                z = []

            # parse GAF/PAF
            if len(t) >= 12 and (t[4] in ["+", "-"]):
                y = alignment(
                    qname=t[0],
                    mapq=int(t[11]),
                    qlen=int(t[1]),
                    qst=int(t[2]),
                    qen=int(t[3]),
                    strand=t[4],
                    path=t[5],
                    tlen=int(t[6]),
                    tst=int(t[7]),
                    ten=int(t[8]),
                )
                if y.mapq < opt.min_mapq:
                    continue
                for i in range(12, len(t)):
                    if t[i][:5] == "cg:Z:":
                        y.cg = t[i][5:]
                    elif t[i][:5] == "ds:Z:":
                        y.ds = t[i][5:]
                    elif t[i][:5] == "tp:A:":
                        tp = t[i][5:]

                if tp != "P":
                    continue  # only primary alignment, filter secondary alignment
                if y.cg is None:
                    continue
            else:  # parse SAM format
                continue
            z.append(y)
        get_indel(opt, z)
        get_breakpoint(opt, z)
