#!/usr/bin/env python

import rich_click as click
from collections import defaultdict
from .gaftools import GafParser
import gzip


@click.group(help="Fusion Benchmarking analysis commands")
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
@click.option(
    "-p", "--prefix", required=True, type=str, help="output prefix for table and figure"
)
def parse(input: str, support_read:int, prefix: str):
    largedel_coords_dict = defaultdict(int)

    print(input)
    print(prefix)

    gaf = GafParser(input)
    output = gzip.open(f"{prefix}.bed.gz", 'wt')

    for read_blocks in gaf.parse_gaf():
        for del_read_seg in read_blocks:
            largedel_coords_dict[del_read_seg] += 1

    del_num = 0
    for key, value in largedel_coords_dict.items():
        if value >= support_read:
            output.write(f"{key[0]}\t{key[1]}\t{key[2]}\tlargedel{del_num}\t{value}\t{key[3]}\n")
            del_num += 1
    output.close()
