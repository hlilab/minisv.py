import gzip
import warnings
import re



def parse_readids(tsv):
    severus_id_dict = {}
    with open(tsv) as inf:
        for line in inf:
            if line.startswith("#SV_ID"):
                continue
            line = line.strip().split(',')
            # only query the tumor SVs supported reads from the read id tsv file
            severus_id_dict[line[0]] = line[1].replace("\"", "").split('; ')
    return severus_id_dict


def parse_msvasm(msvasm):
    
    with gzip.open(msvasm) as gzip_file:
        read_ids = []
        for line in gzip_file:
            line = line.strip().split()
            is_bp = False
            if re.match(r"[><]", line[2].decode('utf-8')):
                is_bp = True
            read_col_info = 5 if is_bp else 3
            read_ids.append(line[read_col_info].decode('utf-8'))

    read_ids = set(read_ids)
    return read_ids


def classify_sv_len(svlen=0):
    if svlen == 0:
        sv_len_tag = 'translocation'
    elif svlen >= 1e6:
        sv_len_tag = '>1M'
    elif svlen < 1e6 and svlen >= 1e5:
        sv_len_tag = '100kb-1M'
    elif svlen < 1e5 and svlen >= 20e3:
        sv_len_tag = '20kb-100kb'
    elif svlen < 20e3 and svlen >= 100:
        sv_len_tag = '100bp-20kb'
    else:
        sv_len_tag = '<100bp'
    return sv_len_tag


def call_filterseverus(severusvcf, readidtsv, msvasm, outstat, consensus_sv_ids, asm_count_cutoff=2):
    parsed_id_dict = parse_readids(readidtsv)
    read_ids = parse_msvasm(msvasm)

    consensus_ids = []
    with open(consensus_sv_ids) as inf:
        for line in inf:
            consensus_ids.append(re.sub("_[12]", "", line.strip()))

    outf = open(outstat, 'w')
    outf.write("svlen\tsvlen_range\tsvid\tis_consensus\tDV\tread_name_number\tasm_support\n")
    vcf_svids = []
    inconsistent_read_num = 0
    total = 0
    with open(severusvcf) as inf:
        for line in inf:
            svlen = 0
            if line.startswith("#"):
                print(line.strip())
                continue

            total += 1
            elements = line.strip().split()
            info = elements[-3].split(';')
            for info_field in info:
                 if info_field.startswith('SVLEN='):
                     svlen = int(info_field.replace('SVLEN=', ''))

            sv_len_tag = classify_sv_len(abs(svlen))
            
            ## Severus FORMAT column: GT:GQ:VAF:hVAF:DR:DV
            svid = elements[2]
            svid = re.sub("_[12]", "", svid)
            vcf_svids.append(svid)

            # if there is no read information, do not output the SVs
            assert svid in parsed_id_dict
            form = elements[-1]
            form_list = form.split(':')
            assert int(form_list[-1]) >= len(set(parsed_id_dict[svid])), f"{svid} {form_list[-1]} {len(set(parsed_id_dict[svid]))} {form_list}"
            inconsistent_read_num += int(form_list[-1]) != len(set(parsed_id_dict[svid]))

            if int(form_list[-1]) != len(set(parsed_id_dict[svid])):
                 w_mess = f"{svid} {int(form_list[-1])} {len(set(parsed_id_dict[svid]))} unequal reads between output read names and vcf DV"
                 warnings.warn(w_mess, UserWarning)

            # if there is more than 2 SVs with self-assembly support, output the SVs
            if len(set(parsed_id_dict[svid]) & read_ids) >= asm_count_cutoff:
                print(line.strip())

            if elements[2].endswith('_2'): # BND_2, for paired events such as translocation, output one end for read filtering
                continue

            ## SVLEN, SVID, DV total tumor SV reads, read TSV file tumor SV read num, asm supported SV
            if sv_len_tag != '<100bp':
                outf.write('\t'.join([str(svlen), sv_len_tag, svid, str(svid in consensus_ids), form_list[-1], str(len(set(parsed_id_dict[svid]))), str(len(set(parsed_id_dict[svid]) & read_ids))])+'\n')

    assert set(vcf_svids).issubset(set(parsed_id_dict.keys()))
    #print(total, inconsistent_read_num)
    outf.close()


def call_filtersnf(snfvcfgz, msvasm, outstat, asm_count_cutoff=2):
    # Sniffle2 suport, DV, read names issues:
    #   https://github.com/fritzsedlazeck/Sniffles/issues/496
    #   https://github.com/fritzsedlazeck/Sniffles/issues/382
    read_ids = parse_msvasm(msvasm)

    outf = open(outstat, 'w')
    outf.write("svlen\tsvlen_range\tsupport\tsvid\tread_name_number\tDV\tasm_support\n")
    with gzip.open(snfvcfgz) as inf:
        for line in inf:
            line = line.decode('utf-8')
            svlen = 0
            support = 0
            snf_reads = []

            if line.startswith("#"):
                print(line.strip())
                continue
            elements = line.strip().split()
            svid = elements[2]
            # -2: tumor, -1: normal
            # GT:GQ:DR:DV:ID
            fmt = elements[8].split(':')
            tumor_form  = elements[-2]
            tumor_form_list = tumor_form.split(':')
            normal_form = elements[-1]
            normal_form_list = normal_form.split(':')

            info = elements[-4].split(';')
            for info_field in info:
                 if info_field.startswith('SVLEN='):
                     svlen = int(info_field.replace('SVLEN=', ''))
            sv_len_tag = classify_sv_len(abs(svlen))

            for info_field in info:
                 if info_field.startswith('RNAMES'):
                     snf_reads = info_field.replace('RNAMES=', '').split(',')
                 if info_field.startswith('SUPPORT'):
                     support = info_field.replace('SUPPORT=', '')

            for i in range(len(fmt)):
                if fmt[i] == 'DV':
                   ##assert len(set(snf_reads)) == int(tumor_form_list[-3]) + int(tumor_form_list[-2]) + int(normal_form_list[-2]) + int(normal_form_list[-3]) 
                   if int(normal_form_list[i]) == 0 and int(tumor_form_list[i]) >= asm_count_cutoff:
                       # SVID, support, svid, read name number, DV, asm supported reads
                       if sv_len_tag != '<100bp':
                           outf.write('\t'.join([str(svlen), sv_len_tag, support, svid, str(len(set(snf_reads))), tumor_form_list[-2], str(len(read_ids & set(snf_reads)))])+'\n')
                       # if there is more than 2 SVs with self-assembly support, output the SVs
                       if len(read_ids & set(snf_reads)) >= asm_count_cutoff:
                           print(line.strip())
    outf.close()


def call_filtermsv(msvtg, msvasm, outstat, asm_count_cutoff=2):

    read_ids = parse_msvasm(msvasm)

    outf = open(outstat, 'w')
    outf.write("svlen\tsvlen_range\tread_name_number\tasm_support\n")
    msv_all_reads = []
    with open(msvtg) as inf:
        for line in inf:
            elements = line.strip().split()
            info = elements[-1].split(';')
            svlen = 0
            msv_reads = []

            for info_field in info:
                 if info_field.startswith('reads='):
                     msv_reads = info_field.replace('reads=', '').split(',')
                 if info_field.startswith('SVLEN='):
                     svlen = int(info_field.replace('SVLEN=', ''))
            msv_all_reads += msv_reads

            sv_len_tag = classify_sv_len(abs(svlen))
            if sv_len_tag != '<100bp':
                 # SVID, support, svid, read name number, DV, asm supported reads
                 outf.write('\t'.join([str(svlen), sv_len_tag, str(len(set(msv_reads))), str(len(read_ids & set(msv_reads)))])+'\n')
            if len(read_ids & set(msv_reads)) >= asm_count_cutoff:
                 print(line.strip())
    ##assert set(msv_all_reads).issubset(read_ids)
    outf.close()
