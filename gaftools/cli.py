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
    help="input json file is generated through prepare",
)
@click.option(
    "-r", "--support_read", required=True, type=int, help="supported read threshold"
)
@click.option("-m", "--mapq", required=True, type=int, help="read mapping quality")
@click.option("-l", "--mlen", required=True, type=int, help="min indel length")
@click.option(
    "-p", "--prefix", required=True, type=str, help="output prefix for table and figure"
)
def parse(input: str, support_read: int, mapq: int, mlen: int, prefix: str):
    print(input)
    print(prefix)
    gaf = GafParser(input, prefix)
    gaf.parse_indel(mapq, mlen)
    gaf.merge_indel(min_cnt=support_read, min_mapq=mapq)
