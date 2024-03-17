""" A generalized input interface for both INDELs and breakpoints
We could put the output interface into this submodule as well.
"""
import argparse
import sys 
from identify_breaks_v5 import call_breakpoints
from merge_break_pts_v3 import merge_breaks

def load_gaf_to_grouped_reads(gafFile, min_mapQ=5, min_map_len=2000):
    """
    Load GAF file, filter based on mapping quality and minimum map length,
    and sort the lines by cluster location in read.

    Args:
    - gafFile (str): Path to the GAF file.
    - min_mapQ (int): Minimum mapping quality threshold (default: 5).
    - min_map_len (int): Minimum mapping length threshold (default: 2000).

    Returns:
    - sorted_lines (list): List of GAF lines sorted by cluster location in read.
    """
    #print(gafFile)

    # Read lines, split into fields, and filter based on conditions
    lines = []
    with open(gafFile, 'r') as gafFileHandler:
        for line in gafFileHandler:        
            fields = line.strip().split('\t')
            read_name = fields[0]

            # Skip low mapping quality, short aligned length
            # For minimap2 paf, we skipped supplementary alignment
            if fields[16] == "tp:A:S" or int(fields[11]) < min_mapQ or int(fields[8]) - int(fields[7]) < min_map_len:
                continue

            if len(lines) >= 1:
                if read_name == lines[-1][0]:
                    lines.append(fields)
                else:
                    # locally sort a group of reads
                    sorted_lines = sorted(lines, key = lambda x: (int(x[2])))  
                    yield sorted_lines
                    lines = [fields]
            else:
                lines.append(fields)
    sorted_lines = sorted(lines, key = lambda x: (int(x[2])))  
    yield sorted_lines
    del lines
    del sorted_lines


if __name__ == '__main__':
    # test script locally
    parser = argparse.ArgumentParser(description='Identify Break Points from GAF input') 
    parser.add_argument('-i', required=True, help='input GAF file')
    args = parser.parse_args()
    groups = load_gaf_to_grouped_reads(args.i)
    all_breaks = [] 
    for group in groups:
        if len(group) > 1:
            brks = call_breakpoints(group) 
            for brk in brks:
              all_breaks.append(brk)  
    sys.stdout.write('\n'.join(merge_breaks(all_breaks))) 
        




