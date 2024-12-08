import math
import re
from dataclasses import dataclass
from typing import Optional

from .merge import svinfo
from .regex import re_info
from operator import itemgetter, attrgetter

from .eval import gc_eval_merge_sv, gc_parse_sv, gc_cmp_sv


def annot(inputfiles, opt):
    """Evaluate the performance of the test SVs against the base SVs
    Args:
        base: list of base SVs
        test: list of test SVs
        opt: EvalConfig object
    Returns:
        None
    """

    min_read_len = math.floor(opt.min_len * opt.read_len_ratio + 0.499)

    vcf = []
    for i in range(len(inputfiles)):
        vcf.append(
            gc_parse_sv(inputfiles[i], min_read_len, opt.min_count, opt.ignore_flt, opt.check_gt)
        )

    # first vcf is to be annotated
    for i in range(len(inputfiles)):
        other = []
        for j in range(len(inputfiles)):
            if i != j:
                other.append(vcf[j])

        merge = gc_eval_merge_sv(opt.win_size, opt.min_len_ratio, other)
        tot_fp, fp = gc_cmp_sv(opt, merge, vcf[i], 'merge')
        break
