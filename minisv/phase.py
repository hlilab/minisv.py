import pysam

def extract_phase_HP(bam):
    samfile = pysam.AlignmentFile(bam, "rc")
    for line in samfile:
        hp = 'NA'
        if line.has_tag("HP"):
            hp = line.get_tag("HP")
        print(line.query_name, hp, sep='\t')

def annotate_HP(hpstdin, msv):
    hptagdict = {}
    for line in hpstdin:
        line = line.strip().split()
        if line[-1] == 'NA':
            continue
        hptagdict[line[0]] = line[-1]

    with open(msv) as infile:
        #reads=m84039_230415_005427_s4/170330520/ccs,m84039_230415_005427_s4/21300188/ccs,m84039_230415_005427_s4/58855421/ccs,m84039_230415_005427_s4/71242084/ccs
        for line in infile:
            splitline = line.strip().split()
            reads = splitline[-1].split(';')[-1].replace("reads=", "")
            hp_counts = [0, 0, 0] # non, HP1, HP2
            for read in reads.split(','):
                if read in hptagdict:
                    hp_counts[int(hptagdict[read])] += 1
                else:
                    hp_counts[0] += 1
            print(line.strip(), ','.join(map(str, hp_counts)), sep='\t')
