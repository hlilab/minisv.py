import argparse
import sys 

parser = argparse.ArgumentParser(description='Identify Break Points from GAF input') 
parser.add_argument('-i', metavar='<input.gaf>', required=True, help='input GAF file')
parser.add_argument('-m', metavar= '--min_mapping_quality', required= False, type=int, default = 10, help = 'Minimum mapping quality (integer)') 
parser.add_argument('-a', metavar= '--min_alignment_length', required= False, type=int, default = 2000, help = 'Minimum alignment length (bps)')

args = parser.parse_args()
gafFile = args.i
min_mapQ = args.m 
min_map_len = args.a


#read in lines and split into feilds
lines = [] 
with open(gafFile, 'r') as f: 
    for line in f: 
        #lines.append(line.strip()) 
        s = line.strip().split('\t')
        if int(s[11])  < min_mapQ:
            continue
        if int(s[8]) - int(s[7]) < min_map_len:
            continue 
        lines.append(s)

#sort lines by read ID and mapping start in read 

sorted_lines = sorted(lines, key = lambda x: (x[0], int(x[2])))

#for line in sorted_lines: 
 #   sys.stdout.write('\t'.join(line) + '\n') 
#sys.stdout.write('\n'.join(sorted_lines))
def get_contig(s, ch, location, is_start):
        if s.startswith('c'):
            return s, str(location) 
        else:            
            nodes = [i for i, ltr in enumerate(s) if ltr in ch]     
            if len(nodes) ==1 : 
                return s, str(location)

            if is_start: 
                total = 0
                for i in range(1, len(nodes)):
                    cur = s[nodes[i-1] :nodes[i]]
                    locs = cur.split(':')[-1].split('-')
                    diff = int(locs[1]) -int(locs[0])
                    total += abs(diff) 
                return s[nodes[-1]:], str(location - total) 
            else: 
                return s[nodes[0]:nodes[1]], str(location) 

def find_break(m1, m2):
    strand_translation = {'+':'>', '-':'<'}
    contig_m1 = get_contig(m1[5], ['>', '<'], int(m1[8]), True)
    contig_m2 = get_contig(m2[5], ['>', '<'], int(m2[7]), False) 
    return [contig_m1[0], contig_m1[1], strand_translation[m1[4]] + strand_translation[m2[4]], contig_m2[0], contig_m2[1], str(min([int(m1[11]), int(m2[11])])),m1[0]]

def convert(contig, brk, direct):
    if not contig.startswith(('<', '>')):
        return contig, brk , direct
    locs = contig.split(':')[-1].split('-')
    if contig.startswith('>'):
        return contig.split(':')[0][1:] ,  str(int(locs[0]) + int(brk)), contig.split(':')[0][0]
    else:
        return contig.split(':')[0][1:] ,  str(int(locs[1]) -  int(brk)), contig.split(':')[0][0]

def adj_graph_paths(original): 
    first = convert(original[0], original[1],original[2].strip()[0])
    second = convert(original[3], original[4], original[2].strip()[1])
    arrows = first[2] + second[2]
    if arrows == '<<':
        temp = first
        first = second
        second = temp
        arrows = '>>'
    if arrows == '<>':
        temp = first
        first = second
        second = temp
        arrows = '><'
    original[0] = first[0]
    original[1] = first[1]
    original[2] = arrows
    original[3] = second[0]
    original[4] = second[1]
    return '\t'.join(original) 

read_maps = []
prev_line = [0,0]
output = [] 
for line in sorted_lines:
    if line[0] == prev_line[0]:
        breaks = find_break(prev_line, line)
        if breaks:
            output.append(adj_graph_paths(breaks))
    prev_line = line
sys.stdout.write('\n'.join(output))

