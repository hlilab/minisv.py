import sys
import numpy as np 
import argparse

parser = argparse.ArgumentParser(description='Identify Break Points from GAF input')
parser.add_argument('-n', metavar='<breakpoints.txt>', required=True, help='input normal break points')
parser.add_argument('-t', metavar='<breakpoints.txt>', required=True, help='input tumor break points')
parser.add_argument('-w', metavar= '--merge_window', required= False, type=int, default = 100, help = 'Size of window to merge break points in')
parser.add_argument('-o', metavar='--output', required=True, help='output files header')
parser.add_argument('-regs', metavar= '--centromere_regions', required= False, default = '', help = 'File of poor mapping region locations')

args = parser.parse_args()
normal = args.n
tumor = args.t
margin = args.w
outPref = args.o
cents = args.regs

#read in given difficult mapping regions file
centromeres = {}
if len(cents) > 0:
    with open(cents, 'r') as f:
        for line in f:
            s = line.strip().split('\t')
            c = [int(s[1]),int(s[2])]
            centromeres[s[0]] = [min(c),max(c)]

#combine normal and tumor BNDs
lines = [] 
with open(normal, 'r') as f:
    for line in f:
        lines.append((line.strip() + '\tn').split('\t')) 

with open(tumor, 'r') as g: 
    for line in g: 
        lines.append((line.strip() + '\tt').split('\t'))
#sort the combinded files together by chrs and locations
sorted_lines = sorted(lines, key = lambda x: (x[0],x[3], int(x[1]), int(x[4])))

#for sor in sorted_lines: 
 #   sys.stdout.write('\t'.join(sor) + '\n')

#function to determine distance of break point from given centromeres
def get_dist_to_cen(chrom, loc):
    if chrom not in centromeres.keys():
        return 'N/A'
    cen_locs = centromeres[chrom]
    dist = 0
    if loc < min(cen_locs):
        dist = min(cen_locs) - loc
    elif loc > max(cen_locs):
        dist = loc -max(cen_locs)
    return str(dist) 

#function that takes a cluster of breakpoints within the margin to merge into one line
def get_merges(cluster): 
    if len(cluster) < 5: 
        return None 
    arr = np.array(cluster)
    ns= 0 
    ts = 0 
    for x in arr[:,7]: 
        if x == 'n': 
            ns+=1
        else: 
            ts +=1    
    start_val = round(np.median([int(i) for i in arr[:,1]]))
    start_dist = get_dist_to_cen(arr[0,0], start_val) 
    end_val = round(np.median([int(i) for i in arr[:,4]]))
    end_dist = get_dist_to_cen(arr[0,3], end_val) 
    return '\t'.join([arr[0,0], str(start_val), arr[0,2], arr[0,3], str(end_val), start_dist+','+end_dist , ','.join(arr[:,5]), ','.join(arr[:,6]), str(ns)+ '/' + str(ts)])

#convert the breakpoint text format into vcf format
def vcf_format(original, number): 
    out = [] 
    s = original.strip().split('\t')
    qual = str(round(np.mean([int(x) for x in s[6].split(',')])))
    info = ';NORMAL_COUNTS='+ str(s[8].split('/')[0]) + ';TUMOR_COUNTS='+ str(s[8].split('/')[1])+";READ_IDS=" + s[7]  
    if s[2] == '>>':
        out.append('\t'.join([s[0],s[1],'bnd_' + str(number), 'N', 'N[' + s[3]+':'+ s[4]+'[', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number+1) + ';CEN_DIST='+ s[5].split(',')[0] + info]))
        out.append('\t'.join([s[3],s[4],'bnd_' + str(number+1), 'N', ']' + s[0]+':'+ s[1]+']N', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number) +';CEN_DIST='+ s[5].split(',')[1] + info]))
    else:
        out.append('\t'.join([s[0],s[1],'bnd_' + str(number), 'N', 'N]' + s[3]+':'+ s[4]+']', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number+1) +';CEN_DIST='+ s[5].split(',')[0] + info])) 
        out.append('\t'.join([s[3],s[4],'bnd_' + str(number+1), 'N', '[' + s[0]+':'+ s[1]+'[N', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number) + ';CEN_DIST='+ s[5].split(',')[1] +info]))
    return out

all_merges= [] 
all_vcfs = ['##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']
prev_line = [0,0,0,0,0,0] 
lines_to_merge = [] 
i=1
for l in sorted_lines:
    #cluster BNDs that are within the margin on both sides of the BND
    if l[0] == prev_line[0] and abs(int(prev_line[1]) - int(l[1])) <= margin and l[3] == prev_line[3] and abs(int(prev_line[4]) - int(l[4])) <= margin:
        lines_to_merge.append(l)
    else: 
        out = get_merges(lines_to_merge)
        if out: 
            vcf_out = vcf_format(out, i)
            all_vcfs.append(vcf_out[0])
            all_vcfs.append(vcf_out[1])
            all_merges.append(out)
            i += 2
        lines_to_merge = [l]         
    prev_line = l 
out = get_merges(lines_to_merge)
if out: 
    vcf_out = vcf_format(out, i)
    all_vcfs.append(vcf_out[0])
    all_vcfs.append(vcf_out[1])
    all_merges.append(out) 
    i += 2 
with open(outPref + '.vcf','w') as final_vcf: 
    final_vcf.write('\n'.join(all_vcfs)) 

with open(outPref + '.txt','w') as final_breaks: 
    final_breaks.write('\n'.join(all_merges)) 

#print(all_vcfs)
#sys.stdout.write('\n'.join(all_merges)) 




