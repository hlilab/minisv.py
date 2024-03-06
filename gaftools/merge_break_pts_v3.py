import sys
import numpy as np 
import argparse

parser = argparse.ArgumentParser(description='Identify Break Points from GAF input')
parser.add_argument('-i', metavar='<breakpoints.txt>', required=True, help='input GAF file')
parser.add_argument('-w', metavar= '--merge_window', required= False, type=int, default = 100, help = 'Size of window to merge break points in')
parser.add_argument('-o', metavar='--output', required=True, help='output files header')
parser.add_argument('-regs', metavar= '--centromere_regions', required= False, default = '', help = 'File of poor mapping region locations')

args = parser.parse_args()
brkFile = args.i
margin = args.w
cents = args.regs
outPref =args.o

centromeres = {}
if len(cents) > 0:
     with open(cents, 'r') as f:
         for line in f:
             s = line.strip().split('\t') 
             c = [int(s[1]),int(s[2])] 
             centromeres[s[0]] = [min(c),max(c)] 

lines = [] 
with open(brkFile, 'r') as f:
    for line in f:
        lines.append(line.strip().split('\t')) 
sorted_lines = sorted(lines, key = lambda x: (x[0],x[3], int(x[1]), int(x[4])))

#for sor in sorted_lines: 
 #   sys.stdout.write('\t'.join(sor) + '\n')
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

def get_merges(cluster): 
    if len(cluster) < 5: 
        return None 
    arr = np.array(cluster)
    start_val = round(np.median([int(i) for i in arr[:,1]]))
    start_dist = get_dist_to_cen(arr[0,0], start_val)
    end_val = round(np.median([int(i) for i in arr[:,4]]))
    end_dist = get_dist_to_cen(arr[0,3], end_val)
    return '\t'.join([arr[0,0], str(start_val), arr[0,2], arr[0,3], str(end_val), start_dist+','+end_dist , ','.join(arr[:,5]), ','.join(arr[:,6])]) 
    #return '\t'.join([arr[0,0], str(round(np.median([int(i) for i in arr[:,1]]))), arr[0,2], arr[0,3], str(round(np.median([int(i) for i in arr[:,4]]))), ';'.join(arr[:,5]), ';'.join(arr[:,6])])

def vcf_format(original, number): 
    out = [] 
    s = original.strip().split('\t')
    qual = str(round(np.mean([int(x) for x in s[6].split(',')])))
    read_count = str(len(s[6].split(';')))
    info = ';READ_COUNTS='+ read_count + ';READ_IDS=' + s[7]
    if s[2] == '>>':
        out.append('\t'.join([s[0],s[1],'bnd_' + str(number), 'N', 'N[' + s[3]+':'+ s[4]+'[', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number+1) +  ';CEN_DIST='+ s[5].split(',')[0] + info]))
        out.append('\t'.join([s[3],s[4],'bnd_' + str(number+1), 'N', ']' + s[0]+':'+ s[1]+']N', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number) +  ';CEN_DIST='+ s[5].split(',')[1] + info]))
    else:
        out.append('\t'.join([s[0],s[1],'bnd_' + str(number), 'N', 'N]' + s[3]+':'+ s[4]+']', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number+1) +  ';CEN_DIST='+ s[5].split(',')[0] + info])) 
        out.append('\t'.join([s[3],s[4],'bnd_' + str(number+1), 'N', '[' + s[0]+':'+ s[1]+'[N', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number) +  ';CEN_DIST='+ s[5].split(',')[1] + info]))
    return out
    #make sure to update i when calling!!! and check on IGV that the strandedness is correct too 
    #i +=2

all_merges= [] 
all_vcfs = ['##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']
prev_line = [0,0,0,0,0,0] 
lines_to_merge = [] 
i=1
for l in sorted_lines:
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


#rewrite this to be done while lines are being parsed the first time rather than going through the dict again 
#for temp in all_merges:  
 #   if len(temp) >= 5: 
  #      arr = np.array(temp) 
   #     sys.stdout.write('\t'.join([arr[0,0], str(round(np.median([int(i) for i in arr[:,1]]))), arr[0,2], arr[0,3], str(round(np.median([int(i) for i in arr[:,4]]))), ';'.join(arr[:,5]), ';'.join(arr[:,6])])+ '\n')   

