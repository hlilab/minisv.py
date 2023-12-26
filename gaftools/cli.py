#!/usr/bin/env python

import rich_click as click

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
    "-p", "--prefix", required=True, type=str, help="output prefix for table and figure"
)
@click.option("-v", "--verbose", is_flag=True, help="verbose option for debug")
def getindel(
    input: str,
    support_read: int,
    mapq: int,
    cpu: int,
    mlen: int,
    prefix: str,
    verbose: bool,
):
    """Get indel reads and merge the microhomology reads into large indels"""
    print("get indel...")
    gaf = GafParser(input, prefix)
    gaf.parse_indel(mapq, mlen, verbose, n_cpus=cpu)
    # gaf.merge_indel(min_cnt=support_read, min_mapq=mapq)


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
@click.option("-l", "--mlen", required=True, type=int, help="min indel length")
@click.option(
    "-p", "--prefix", required=True, type=str, help="output prefix for table and figure"
)
@click.option("-v", "--verbose", is_flag=True, help="verbose option for debug")
def gettsd(
    input: str, support_read: int, mapq: int, mlen: int, prefix: str, verbose: bool
):
    """Get TSD reads"""
    print("get TSD...")
    gaf = GafParser(input, prefix)
    gaf.parse_tsd(mapq, mlen, verbose)
