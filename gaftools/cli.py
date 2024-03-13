#!/usr/bin/env python

import sys
import rich_click as click

from .cygaftools import GafParser as cy_GafParser
from .gaftools import GafParser


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
    help="optional vntr sites for flagging potential false discovery"
)
@click.option(
    "--cent",
    required=False,
    type=click.Path(exists=True),
    help="optional centromere sites for flagging potential false discovery"
)
@click.option(
    "--l1",
    required=False,
    type=click.Path(exists=True),
    help="optional L1 sequence fasta for filtering L1 elements analysis"
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
    "-s", "--assembly", required=True, type=str, help="assembly version, e.g., chm13graph,chm13linear,grch37graph,grch37linear,grch38graph,grch38linear"
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
    all_breakpts = gafs.parse_sv_on_group_reads(mapq, mlen, maplen, verbose, n_cpus=cpu, ds=ds)

    gafs.merge_indel(min_cnt=support_read, min_mapq=mapq)
    gafs.merge_breakpts(all_breakpts)

    gafs.bed2vcf(command=command)
    #gafs.bed2breakpoint_vcf(min_cnt=support_read, min_mapq=mapq)


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
