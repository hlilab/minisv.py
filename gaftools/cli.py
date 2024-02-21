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
def getindel(
    input: str,
    vntr: str,
    cent: str,
    support_read: int,
    mapq: int,
    cpu: int,
    mlen: int,
    normal: str,
    prefix: str,
    verbose: bool,
):
    """Get indel reads and merge the microhomology reads into large indels"""
    print("get indel...")
    command = " ".join(sys.argv)
    print(command)

    ds = False
    if input.replace(".paf", "") == input and normal.replace(".paf", "") == normal:
        # minimap2 output do not have ds:Z tag yet
        ds = True

    input_samples = [input, normal] if normal is not None else [input]
    gafs = GafParser(input_samples, prefix, vntr, cent)
    gafs.parse_indel(mapq, mlen, verbose, n_cpus=cpu, ds=ds)
    gafs.merge_indel(min_cnt=support_read, min_mapq=mapq)
    gafs.bed2vcf(command=command)


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
