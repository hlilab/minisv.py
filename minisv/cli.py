#!/usr/bin/env python

import re
import sys
from dataclasses import dataclass
from typing import Optional

import rich_click as click

from .eval import gc_read_bed
from .phase import extract_phase_HP, annotate_HP

##from .cygafcall import GafParser as cy_GafParser
from .minisv import GafParser
from .filtercaller import call_filterseverus, call_filtersnf, call_filtermsv
from .io import gc_cmd_view, merge_indel_breakpoints, parseNum, write_vcf

__version__ = "0.1.2"


@dataclass
class opt:
    cen: dict
    min_mapq: int = 5
    min_mapq_end: int = 30
    min_frac: float = 0.7
    min_len: int = 100
    min_aln_len_end: int = 2000
    min_aln_len_mid: int = 50
    max_cnt_10k: int = 5
    dbg: bool = False
    polyA_pen: int = 5
    polyA_drop: int = 100
    name: str = "foo"


@dataclass
class mergeopt:
    min_cnt: int = 3
    min_cnt_strand: int = 2
    min_cnt_rt: int = 1
    min_rt_len: int = 10
    win_size: int = 100
    max_diff: float = 0.05
    min_cen_dist: int = 5e5
    max_allele: int = 100
    max_check: int = 500


@dataclass
class EvalOpt:
    min_len: int = 100
    read_len_ratio: float = 0.8
    win_size: int = 500
    min_len_ratio: float = 0.6
    dbg: bool = False
    print_err: bool = False
    bed: Optional[str] = None
    min_vaf: float = 0
    ignore_flt: bool = False
    check_gt: bool = False


@click.group(help="minisv tool commands")
@click.version_option(__version__)
@click.pass_context
def cli(ctx):
    pass


#@cli.command()
#@click.option(
#    "-i",
#    "--input",
#    required=True,
#    type=click.Path(exists=True),
#    help="tumor or normal sample input gaf/paf file, this is required input, if paired normal sample provided, this should be input tumor sample",
#)
#@click.option(
#    "--vntr",
#    required=False,
#    type=click.Path(exists=True),
#    help="optional vntr sites for flagging potential false discovery",
#)
#@click.option(
#    "--cent",
#    required=False,
#    type=click.Path(exists=True),
#    help="optional centromere sites for flagging potential false discovery",
#)
#@click.option(
#    "--l1",
#    required=False,
#    type=click.Path(exists=True),
#    help="optional L1 sequence fasta for filtering L1 elements analysis",
#)
#@click.option(
#    "-r", "--support_read", required=True, type=int, help="supported read threshold"
#)
#@click.option("-m", "--mapq", required=True, type=int, help="read mapping quality")
#@click.option(
#    "-c", "--cpu", required=True, type=int, default=1, help="read mapping quality"
#)
#@click.option("-l", "--mlen", required=True, type=int, help="min indel length")
#@click.option("-a", "--maplen", required=True, type=int, help="min mapping length")
#@click.option(
#    "-n",
#    "--normal",
#    required=False,
#    type=str,
#    help="paired normal sample input gaf/paf file if tumor-only and normal-only mode, skip this option",
#)
#@click.option(
#    "-p", "--prefix", required=True, type=str, help="output prefix for table and figure"
#)
#@click.option(
#    "-s",
#    "--assembly",
#    required=True,
#    type=str,
#    help="assembly version, e.g., chm13graph,chm13linear,grch37graph,grch37linear,grch38graph,grch38linear",
#)
#@click.option("-d", "--ds", is_flag=True, help="support ds tag or not")
#@click.option("-v", "--verbose", is_flag=True, help="verbose option for debug")
#def getsv(
#    input: str,
#    vntr: str,
#    cent: str,
#    l1: str,
#    support_read: int,
#    mapq: int,
#    cpu: int,
#    mlen: int,
#    maplen: int,
#    normal: str,
#    prefix: str,
#    ds: bool,
#    assembly: str,
#    verbose: bool,
#):
#    """Get indel reads and merge the microhomology reads into large indels"""
#    print("get indel...")
#    command = " ".join(sys.argv)
#    print(command)
#
#    input_samples = [input, normal] if normal is not None else [input]
#    gafs = GafParser(input_samples, prefix, assembly, vntr, cent, l1)
#    all_breakpts = gafs.parse_sv_on_group_reads(
#        mapq, mlen, maplen, verbose, n_cpus=cpu, ds=ds
#    )
#
#    gafs.merge_indel(min_cnt=support_read, min_mapq=mapq)
#    merged_breakpt_vcf_list = gafs.merge_breakpts(all_breakpts)
#
#    # output bedpe format
#    merge_indel_breakpoints(prefix, gafs.breakpt_file, gafs.indel_file)
#
#    # output vcf format
#    gafs.bed2vcf(command=command, merged_breakpt_vcf_list=merged_breakpt_vcf_list)


@cli.command()
@click.option("-n", required=False, default="foo", type=str, help="sample name")
@click.option("-q", required=False, default=5, type=int, help="minimum mapping quality")
@click.option(
    "-b", required=False, type=click.Path(exists=True), help="centromere bed file"
)
@click.option(
    "-l", "--svlen", required=False, default="100", type=str, help="minimum sv length"
)
@click.option(
    "-f", required=False, default=0.7, type=float, help="min fraction of reads"
)
@click.option("-c", required=False, default=3, help="maximum count of sv per 10k")
@click.option("-a", required=False, type=int, default=5, help="polyA penalty")
@click.option("-d", is_flag=True, help="verbose option for debug")
@click.option(
    "-x",
    required=False,
    default=30,
    type=int,
    help="minimum mapping quality in the end",
)
@click.option(
    "-e",
    required=False,
    default=2000,
    type=int,
    help="minimum aligned length in the end",
)
@click.option(
    "-m",
    required=False,
    default=50,
    type=int,
    help="minimum aligned length in the middle",
)
@click.argument("filename", nargs=1)
def getsv(
    q: int,
    x: int,
    svlen: str,
    d: bool,
    f: float,
    c: int,
    a: int,
    e: int,
    m: int,
    n: str,
    b: str,
    filename: tuple,
):
    """Extract raw INDEL and breakpoints"""
    from .read_parser import load_reads

    options = opt(
        min_mapq=q,
        min_mapq_end=x,
        min_len=parseNum(svlen),
        dbg=d,
        min_frac=f,
        max_cnt_10k=c,
        polyA_pen=a,
        min_aln_len_end=e,
        min_aln_len_mid=m,
        name=n,
        cen={},
    )

    if b is not None:
        parse_centromere(b, options.cen)
    load_reads(filename, options)


@cli.command()
@click.option("-n", required=False, default="foo", type=str, help="sample name")
@click.option("-q", required=False, default=5, type=int, help="minimum mapping quality")
@click.option(
    "-b", required=False, type=click.Path(exists=True), help="centromere bed file"
)
@click.option(
    "-l", "--svlen", required=False, default="100", type=str, help="minimum sv length"
)
@click.option(
    "-f", required=False, default=0.7, type=float, help="min fraction of reads"
)
@click.option("-c", required=False, default=3, help="maximum count of sv per 10k")
@click.option("-a", required=False, type=int, default=5, help="polyA penalty")
@click.option("-d", is_flag=True, help="verbose option for debug")
@click.option(
    "-x",
    required=False,
    default=30,
    type=int,
    help="minimum mapping quality in the end",
)
@click.option(
    "-e",
    required=False,
    default=2000,
    type=int,
    help="minimum aligned length in the end",
)
@click.option(
    "-m",
    required=False,
    default=50,
    type=int,
    help="minimum aligned length in the middle",
)
@click.argument("filename", nargs=1)
def extract(
    q: int,
    x: int,
    svlen: str,
    d: bool,
    f: float,
    c: int,
    a: int,
    e: int,
    m: int,
    n: str,
    b: str,
    filename: tuple,
):
    """Extract raw INDEL and breakpoints"""
    from .read_parser import load_reads

    options = opt(
        min_mapq=q,
        min_mapq_end=x,
        min_len=parseNum(svlen),
        dbg=d,
        min_frac=f,
        max_cnt_10k=c,
        polyA_pen=a,
        min_aln_len_end=e,
        min_aln_len_mid=m,
        name=n,
        cen={},
    )

    if b is not None:
        parse_centromere(b, options.cen)
    load_reads(filename, options)


def parse_centromere(b, cen):
    with open(b) as centromere_file:
        for line in centromere_file:
            t = line.strip().split("\t")
            if t[0] not in cen:
                cen[t[0]] = []
            cen[t[0]].append([int(t[1]), int(t[2])])
        for ctg in cen:
            cen[ctg].sort(key=lambda x: x[0])


@cli.command()
@click.option("-w", required=False, default=100, type=int, help="window size")
@click.option(
    "-d",
    required=False,
    default=0.05,
    type=float,
    help="maximum allele length different ratio",
)
@click.option("-c", required=False, default=4, type=int, help="minimum sv counts")
@click.option(
    "-s", required=False, default=2, type=int, help="minimum count per strand"
)
@click.option(
    "-r",
    required=False,
    default=0,
    type=int,
    help="min min(TSD_len,polyA_len) to tag a candidate RT",
)
@click.option("-rr", required=False, default=1, type=int, help="minimum count for RT")
@click.option("-a", required=False, default=100, type=int, help="maximum allele number")
@click.option(
    "-cc",
    required=False,
    default=500,
    type=int,
    help="compare up to INT reads per allele",
)
@click.option(
    "-e",
    required=False,
    default="500k",
    type=str,
    help="minimum distance to centromere",
)  # NOTE: in mgutils, this is forced to be int
@click.argument("input", type=click.File("r"), required=True, nargs=1)
def merge(
    w: int, d: float, c: int, s: int, r: int, rr: int, a: int, cc: int, e: str, input
):
    """Usage: sort -k 1,1 -k2,2n sv.bed | gafcall merge [options] -"""
    from .merge import merge_sv

    options = mergeopt(
        win_size=w,
        max_diff=d,
        min_cnt=c,
        min_cnt_rt=rr,
        min_cen_dist=parseNum(e),
        min_rt_len=r,
        min_cnt_strand=s,
        max_allele=a,
        max_check=cc,
    )
    merge_sv(options, input)


@cli.command()
@click.option("-w", required=False, default="500", type=str, help="window size")
@click.option(
    "-svlen", required=False, default="100", type=str, help="sv minimum length"
)
@click.option(
    "-v", required=False, default=0, type=float, help="ignore VAF below FLOAT"
)
@click.option("-c", required=False, default=3, type=int, help="minimum sv counts")
@click.option(
    "-lenratio",
    required=False,
    default=0.6,
    type=float,
    help="minimum length ratio similarity between read SV and truthset SV",
)
@click.option(
    "-r",
    required=False,
    default=0.8,  # NOTE: in js this might be a bug, these ratio name looks confusing
    type=float,
    help="read SVs longer than svlen*FLOAT",
)
@click.option("-b", type=click.Path(exists=True), help="bed to restrict the comparison")
@click.option("-d", is_flag=True, help="verbose option for debug")
@click.option("-e", is_flag=True, help="print errors")
@click.option("-g", is_flag=True, help="check GT")
@click.option("-f", is_flag=True, help="ignore VCF filter")
@click.argument("filename", nargs=-1)
def eval(
    w: str,
    svlen: str,
    c: int,
    r: float,
    lenratio: float,
    d: bool,
    e: bool,
    v: float,
    b,
    g,
    f,
    filename,
):
    """Evaluation of SV calls"""
    from .eval import eval

    options = EvalOpt(
        min_len=parseNum(svlen),
        win_size=parseNum(w),
        read_len_ratio=r,
        min_len_ratio=lenratio,  # NOTE: the option are not input
        dbg=d,
        bed=gc_read_bed(b) if b is not None else None,
        print_err=e,
        min_vaf=v,
        check_gt=g,
        ignore_flt=f,
    )
    eval(filename, options)


@dataclass
class viewopt:
    min_read_len: int
    ignore_flt: bool
    check_gt: bool
    count_long: bool
    bed: Optional[str]


@cli.command()
@click.option(
    "-minlen", required=False, default="100", type=str, help="minimum sv length"
)
@click.option("-ignoreflt", is_flag=True, help="ignore FILTER field in VCF")
@click.option("-gt", is_flag=True, help="check GT in VCF")
@click.option("-c", is_flag=True, help="count 20kb,100kb,1Mb and translocations")
@click.option(
    "-b", required=False, type=click.Path(exists=True), help="restricted bed file"
)
@click.argument("input", nargs=-1)
def view(minlen: int, ignoreflt: bool, gt: bool, c: bool, b: str, input):
    opt = viewopt(
        min_read_len=minlen, ignore_flt=ignoreflt, check_gt=gt, count_long=c, bed=b
    )
    gc_cmd_view(opt, input)


@cli.command()
@click.option(
    "-minlen", required=False, default=100, type=int, help="minimum sv length"
)
@click.option("-ignoreflt", is_flag=True, help="ignore FILTER field in VCF")
@click.option("-gt", is_flag=True, help="check GT in VCF")
@click.option("-c", is_flag=True, help="count 20kb,100kb,1Mb and translocations")
@click.option(
    "-b", required=False, type=click.Path(exists=True), help="restricted bed file"
)
@click.argument("filter_input_gsv", type=click.File("r"), default=sys.stdin, nargs=1)
@click.argument("output_gsv", type=click.File("r"), default=sys.stdin, nargs=1)
def join(
    minlen: int,
    ignoreflt: bool,
    gt: bool,
    c: bool,
    b: str,
    filter_input_gsv,
    output_gsv,
):
    """Usage: gafcall join filter.gsv out.gsv"""

    h = {}

    def get_type(t, col_info):
        info = t[col_info]
        m = re.findall(r"\bSVTYPE=([^\s;])+", info)

        if len(m) > 0:
            if m[0] == "INS" or m[0] == "DUP":
                return 1
            elif m[0] == "DEL":
                return 2
            elif m[0] == "INV":
                return 4
            elif m[0] == "BND" and col_info == 8 and t[0] != t[3]:
                return 8
        return 0

    # file as filter
    for line in filter_input_gsv:
        t = line.strip().split()
        if re.match(r"[><]", t[2]):
            col_info = 8
        else:
            col_info = 6
        name = t[col_info - 3]
        if name not in h:
            h[name] = 0
        h[name] |= get_type(t, col_info)

    # file to be filtered
    for line in output_gsv:
        t = line.strip().split()
        if re.match(r"[><]", t[2]):
            col_info = 8
        else:
            col_info = 6
        name = t[col_info - 3]
        if name not in h:
            continue
        type = get_type(t, col_info)
        # NOTE: h[name] & 8 is redundant?
        #       shall we check chromosome?
        # TODO: add read query positions
        if type == 0 or (h[name] & type) or (h[name] & 8):
            print(line.strip())


@cli.command()
@click.argument("input", type=click.File("r"), default=sys.stdin)
def vcf(input):
    for line in input:
        if line.startswith("#CHROM"):
            print(
                """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotyping quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">"""
            )
            print(line.strip())
        elif line.startswith("##"):
            print(line.strip())
        else:
            line = line.split()
            line[6] = "PASS"
            line = "\t".join(line)
            # print(line.strip(), "GT:GQ:VAF:DR:DV", "0/1:.:.:.:.", sep="\t")
            print(line.strip(), sep="\t")


@cli.command()
@click.argument("input", type=click.File("r"), default=sys.stdin)
def formatvcf(input):
    options = None
    write_vcf(options, input)


@cli.command()
@click.argument("input", type=click.File("r"), default=sys.stdin)
def extracthp(input):
    options = None
    extract_phase_HP(input)


@cli.command()
@click.argument("hptagtsv", type=click.File("r"), default=sys.stdin, nargs=1)
@click.argument("msv", type=str, nargs=1)
def annotatehp(hptagtsv, msv):
    annotate_HP(hptagtsv, msv)


@cli.command()
@click.argument("severusvcf", type=str, nargs=1)
@click.argument("readidtsv", type=str, nargs=1)
@click.argument("msvasm", type=str, nargs=1)
@click.argument("outstat", type=str, nargs=1)
@click.argument("consensus_ids", type=str, nargs=1)
def filterseverus(severusvcf, readidtsv, msvasm, outstat, consensus_ids):
    """
    filter Severus results based on read ids overlap with graph alignment/self alignment
    """
    call_filterseverus(severusvcf, readidtsv, msvasm, outstat, consensus_ids)


@cli.command()
@click.argument("snfvcf", type=str, nargs=1)
@click.argument("msvasm", type=str, nargs=1)
@click.argument("outstat", type=str, nargs=1)
def filtersnf(snfvcf, msvasm, outstat):
    """
    filter Sniffles2 results based on read ids overlap with graph alignment/self alignment
    """
    call_filtersnf(snfvcf, msvasm, outstat)


@cli.command()
@click.argument("msvtg", type=str, nargs=1)
@click.argument("msvasm", type=str, nargs=1)
@click.argument("outstat", type=str, nargs=1)
def filtermsv(msvtg, msvasm, outstat):
    """
    filter Sniffles2 results based on read ids overlap with graph alignment/self alignment
    """
    call_filtermsv(msvtg, msvasm, outstat)


@cli.command()
@click.option(
    "-i",
    "--input",
    required=True,
    type=click.Path(exists=True),
    help="tumor or normal sample input gaf/paf file, this is required input, if paired normal sample provided, this should be input tumor sample",
)
@click.option(
    "--vntr",
    required=False,
    type=click.Path(exists=True),
    help="optional vntr sites for flagging potential false discovery",
)
@click.option(
    "--cent",
    required=False,
    type=click.Path(exists=True),
    help="optional centromere sites for flagging potential false discovery",
)
@click.option(
    "--l1",
    required=False,
    type=click.Path(exists=True),
    help="optional L1 sequence fasta for filtering L1 elements analysis",
)
@click.option(
    "-r", "--support_read", required=True, type=int, help="supported read threshold"
)
@click.option("-m", "--mapq", required=True, type=int, help="read mapping quality")
@click.option(
    "-c", "--cpu", required=True, type=int, default=1, help="read mapping quality"
)
@click.option("-l", "--mlen", required=True, type=int, help="min indel length")
@click.option("-a", "--maplen", required=True, type=int, help="min mapping length")
@click.option(
    "-n",
    "--normal",
    required=False,
    type=str,
    help="paired normal sample input gaf/paf file if tumor-only and normal-only mode, skip this option",
)
@click.option(
    "-p", "--prefix", required=True, type=str, help="output prefix for table and figure"
)
@click.option(
    "-s",
    "--assembly",
    required=True,
    type=str,
    help="assembly version, e.g., chm13graph,chm13linear,grch37graph,grch37linear,grch38graph,grch38linear",
)
@click.option("-d", "--ds", is_flag=True, help="support ds tag or not")
@click.option("-v", "--verbose", is_flag=True, help="verbose option for debug")
def getindel(
    input: str,
    vntr: str,
    cent: str,
    l1: str,
    support_read: int,
    mapq: int,
    cpu: int,
    mlen: int,
    maplen: int,
    normal: str,
    prefix: str,
    ds: bool,
    assembly: str,
    verbose: bool,
):
    """Get indel reads and merge the microhomology reads into large indels"""
    print("get indel...")
    command = " ".join(sys.argv)
    print(command)

    input_samples = [input, normal] if normal is not None else [input]
    gafs = GafParser(input_samples, prefix, assembly, vntr, cent, l1)
    gafs.parse_indel(mapq, mlen, verbose, n_cpus=cpu, ds=ds)
    gafs.merge_indel(min_cnt=support_read, min_mapq=mapq)


# @cli.command()
# @click.option(
#    "-i",
#    "--input",
#    required=True,
#    type=click.Path(exists=True),
#    help="input gaf/paf file",
# )
# @click.option(
#    "-r", "--support_read", required=True, type=int, help="supported read threshold"
# )
# @click.option("-m", "--mapq", required=True, type=int, help="read mapping quality")
# @click.option(
#    "-c", "--cpu", required=True, type=int, default=1, help="read mapping quality"
# )
# @click.option("-l", "--mlen", required=True, type=int, help="min indel length")
# @click.option(
#    "-n",
#    "--normal",
#    required=False,
#    type=str,
#    help="paired normal sample input gaf/paf file if tumor-only and normal-only mode, skip this option",
# )
# @click.option(
#    "-p", "--prefix", required=True, type=str, help="output prefix for table and figure"
# )
# @click.option("-v", "--verbose", is_flag=True, help="verbose option for debug")
# def getindel_cython(
#    input: str,
#    support_read: int,
#    mapq: int,
#    cpu: int,
#    mlen: int,
#    normal: str,
#    prefix: str,
#    verbose: bool,
# ):
#    """Get indel reads and merge the microhomology reads into large indels"""
#    print("get cython indel...")
#    command = " ".join(sys.argv)
#    input_samples = [input, normal] if normal is not None else [input]
#    gafs = cy_GafParser(input_samples, prefix)
#    gafs.parse_indel(mapq, mlen, verbose, n_cpus=cpu)
#    gafs.merge_indel(min_cnt=support_read, min_mapq=mapq)
#    gafs.bed2vcf(command=command)
