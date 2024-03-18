import argparse
import sys 

#function to read in file and sort by read id
def load_gaf_for_breakpoints(gafFile, min_mapQ=10, min_map_len=2000):
    #read in lines and split into feilds
    lines = [] 
    #dont read in lines that are below mapping Q or min map length
    with open(gafFile, 'r') as f:
        for line in f: 
           
            s = line.strip().split('\t')
            if int(s[11])  < min_mapQ:
                continue
            if int(s[8]) - int(s[7]) < min_map_len:
                continue 
            lines.append(s)
    #sort lines by read ID and mapping start in read 
    
    sorted_lines = sorted(lines, key = lambda x: (int(x[2])))  
    return sorted_lines

#function to return final node in graph path as this is where the BND is occuring
def get_contig(s, ch, location, node_end):
    #ignore if not in graph path
    #NOTE: hg19 do not have chr as start string
    #      we extended to hg19 with hard-coded conditions
    if s.startswith('chr') or s.isdigit() or s.startswith("GL") or s.startswith('hs') or s in ['MT', 'X', 'Y']:
        return s, str(location) 
    else:
        nodes = [i for i, ltr in enumerate(s) if ltr in ch]     
        if len(nodes) ==1 : 
            return s, str(location)
        #start of breakpoint, want end coord of read mapping, so last node
        if node_end: 
            total = 0
            for i in range(1, len(nodes)):
                cur = s[nodes[i-1] :nodes[i]]
                locs = cur.split(':')[-1].split('-')
                diff = int(locs[1]) -int(locs[0])
                total += abs(diff) 
            return s[nodes[-1]:], str(location - total) 
        #end of break want start coord, so first node 
        else: 
            return s[nodes[0]:nodes[1]], str(location) 

#function to find the break point between two chuncks of the read mapping
def find_break(m1, m2):
    #translate strands into breakpoint notation
    strand_translation = {'+':'>', '-':'<'}
    # determine break locations based on strand 
    b1 = 0
    if m1[4] == '+':
        b1 =  int(m1[8]) 
        contig_m1 = get_contig(m1[5], ['>', '<'], b1, True)
    elif m1[4] == '-':
        b1 = int(m1[7])
        contig_m1 = get_contig(m1[5], ['>', '<'], b1, False)
    b2 = 0
    if m2[4] == '+':
        b2 = int(m2[7]) 
        contig_m2 = get_contig(m2[5], ['>', '<'], b2, False)
    elif m2[4]  == '-':
        b2 = int(m2[8])
        contig_m2 = get_contig(m2[5], ['>', '<'], b2, True)

    #get info of first chunk of read
    #contig_m1 = get_contig(m1[5], ['>', '<'], b1, True)
    #second chunk of read 
    #contig_m2 = get_contig(m2[5], ['>', '<'], b2, False) 
    #return a more legible break point text format
    return [contig_m1[0], contig_m1[1], strand_translation[m1[4]] + strand_translation[m2[4]], contig_m2[0], contig_m2[1], str(min([int(m1[11]), int(m2[11])])),m1[0]]

#function to convert node locations into chromosome/contig locations 
def convert(contig, brk, direct):
    if not contig.startswith(('<', '>')):
        return contig, brk , direct
    locs = contig.split(':')[-1].split('-')
    if contig.startswith('>'):
        return contig.split(':')[0][1:] ,  str(int(locs[0]) + int(brk)), contig.split(':')[0][0]
    else:
        return contig.split(':')[0][1:] ,  str(int(locs[1]) -  int(brk)), contig.split(':')[0][0]

#swap "<<" to ">>" and "<>" to "><" based on symmetries
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
    return original 

#calling all the functions above
#def output_breakpoints_bed(sorted_lines): 
    #prev_line = [0,0]
    #output = [] 
    #go through all sorted lines
    #for line in sorted_lines:
        #get read ids with supp. mappings
       # if line[0] == prev_line[0]:
         #       output.append(adj_graph_paths(breaks))
        #prev_line = line
    #write to stdout for now, but can just pass array straight into merge function
    #sys.stdout.write('\n'.join(output))


def call_breakpoints(read_cluster, min_mapQ=10, min_map_len=2000):
    #sorted_lines = load_gaf_for_breakpoints(gafFile, min_mapQ=10, min_map_len=2000)
    output = []      
    for i in range(1, len(read_cluster)):
        breaks = find_break(read_cluster[i-1], read_cluster[i])
        if breaks:
             output.append(adj_graph_paths(breaks))
    return output 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Identify Break Points from GAF input') 
    parser.add_argument('-i', metavar='<input.gaf>', required=True, help='input GAF file')
    parser.add_argument('-m', metavar= '--min_mapping_quality', required= False, type=int, default = 10, help = 'Minimum mapping quality (integer)') 
    parser.add_argument('-a', metavar= '--min_alignment_length', required= False, type=int, default = 2000, help = 'Minimum alignment length (bps)')

    args = parser.parse_args()
    gafFile = args.i
    min_mapQ = args.m 
    min_map_len = args.a
    call_breakpoints(gafFile, min_mapQ, min_map_len)
