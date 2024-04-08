#!/usr/bin/env python

import sys
from dataclasses import dataclass

import rich_click as click

from .cygaftools import GafParser as cy_GafParser
from .gaftools import GafParser
from .io import merge_indel_breakpoints


@dataclass
class opt:
    cen: dict
    min_mapq: int = 5
    min_mapq_end: int = 30
    min_frac: float = 0.7
    min_len: int = 100
    min_aln_len_end: int = 2000
    min_aln_len_mid: int = 50
    max_cnt: int = 5
    dbg: bool = False
    polyA_pen: int = 5
    polyA_drop: int = 100
    name: str = "foo"


@click.group(help="Pangenome SV tool commands")
def cli():
    pass


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
def getsv(
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
    all_breakpts = gafs.parse_sv_on_group_reads(
        mapq, mlen, maplen, verbose, n_cpus=cpu, ds=ds
    )

    gafs.merge_indel(min_cnt=support_read, min_mapq=mapq)
    merged_breakpt_vcf_list = gafs.merge_breakpts(all_breakpts)

    # output bedpe format
    merge_indel_breakpoints(prefix, gafs.breakpt_file, gafs.indel_file)

    # output vcf format
    gafs.bed2vcf(command=command, merged_breakpt_vcf_list=merged_breakpt_vcf_list)


@cli.command()
@click.option("-q", required=False, default=5, type=int, help="minimum mapping quality")
@click.option(
    "-x",
    required=False,
    default=30,
    type=int,
    help="minimum mapping quality in the end",
)
@click.option(
    "-l", "--svlen", required=False, default=100, type=int, help="minimum sv length"
)
@click.option("-d", is_flag=True, help="verbose option for debug")
@click.option(
    "-f", required=False, default=0.7, type=float, help="min fraction of reads"
)
@click.option("-c", required=False, default=5, help="maximum count of sv in a read")
@click.option("-a", required=False, type=int, default=5, help="polyA penalty")
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
@click.option("-n", required=False, default="foo", type=str, help="sample name")
@click.option(
    "-b", required=False, type=click.Path(exists=True), help="centromere bed file"
)
@click.argument("filename", nargs=-1)
def sv(
    q: int,
    x: int,
    svlen: int,
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
    from .read_parser import load_reads

    options = opt(
        min_mapq=q,
        min_mapq_end=x,
        min_len=svlen,
        dbg=d,
        min_frac=f,
        max_cnt=c,
        polyA_pen=a,
        min_aln_len_end=e,
        min_aln_len_mid=m,
        name=n,
        cen={},
    )

    if b is not None:
        with open(b) as centromere_file:
            for line in centromere_file:
                t = line.strip().split("\t")
                if t[0] not in options.cen:
                    options.cen[t[0]] = []
                options.cen[t[0]].append([int(t[1]), int(t[2])])
            for ctg in options.cen:
                options.cen[ctg].sort(key=lambda x: x[0])
    click.echo(options)
    click.echo(filename)
    load_reads(filename[0], options)


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


@cli.command()
@click.option(
    "-i",
    "--input",
    required=True,
    type=click.Path(exists=True),
    help="input gaf/paf file",
)
@click.option(
    "-r", "--support_read", required=True, type=int, help="supported read threshold"
)
@click.option("-m", "--mapq", required=True, type=int, help="read mapping quality")
@click.option(
    "-c", "--cpu", required=True, type=int, default=1, help="read mapping quality"
)
@click.option("-l", "--mlen", required=True, type=int, help="min indel length")
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
@click.option("-v", "--verbose", is_flag=True, help="verbose option for debug")
def getindel_cython(
    input: str,
    support_read: int,
    mapq: int,
    cpu: int,
    mlen: int,
    normal: str,
    prefix: str,
    verbose: bool,
):
    """Get indel reads and merge the microhomology reads into large indels"""
    print("get cython indel...")
    command = " ".join(sys.argv)
    input_samples = [input, normal] if normal is not None else [input]
    gafs = cy_GafParser(input_samples, prefix)
    gafs.parse_indel(mapq, mlen, verbose, n_cpus=cpu)
    gafs.merge_indel(min_cnt=support_read, min_mapq=mapq)
    gafs.bed2vcf(command=command)
