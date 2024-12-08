import gzip
from .regex import re_info
import math
import warnings
import re
from .util import is_gzipped, is_vcf
from .eval import gc_parse_sv, iit_overlap
from .type import get_type, simple_type


def parse_readids(file_path):
    tsv = file_path
    read_id_dict = {}
    caller = 'snf' if is_vcf(tsv) or is_gzipped(file_path) else ''

    if is_gzipped(file_path):
        f = gzip.open(file_path, 'rt')
    else:
        f = open(file_path)

    line_no = 0
    for line in f:
        line_no += 1
        if caller == 'snf':
            if line.startswith("#"):
                continue
        if line.startswith("#SV_ID"):
            caller = 'severus'
            continue
        if line.startswith('VARIANT_ID'):
            caller = 'savana'
            continue
        if line.startswith('chr') and line.strip()[-1] in ['+', '-']:
            caller = 'nanomonsv'
        if 'avg_mapq' in line:
            caller = 'msv'

        if caller == 'severus':
            line = line.strip().split(',')
        elif caller == 'savana':
            line = line.strip().split('\t')
        elif caller == 'nanomonsv':
            line = line.strip().split('\t')
        elif caller == 'snf':
            line = line.strip().split('\t')
        elif caller == 'msv':
            line = line.strip().split('\t')
        else:
            raise Exception("not support yet")

        # only query the tumor SVs supported reads from the read id tsv file
        if caller == 'severus':
            read_id_dict[line[0]] = line[1].replace("\"", "").split('; ')
        elif caller == 'savana':
            read_id_dict[line[0]] = line[1].replace("\"", "").split(',')
        elif caller == 'nanomonsv':
            read_id_dict[line[7]] = read_id_dict.get(line[7], []) + [line[8]]
        elif caller == 'snf':
            info = re_info.findall(line[7])
            snf_reads = []
            for info_field, info_val in info:
                 if info_field.startswith('RNAMES'):
                     snf_reads = info_val.split(',')
                     break
            read_id_dict[line[2]] = snf_reads
        elif caller == 'msv':
            is_bp = False
            if re.match(r"[><]", line[2]):
                is_bp = True
            col_info = 8 if is_bp else 6
            info = line[col_info]
            info = re_info.findall(info)
            msv_reads = []
            for info_field, info_val in info:
                 if info_field.startswith('reads'):
                     msv_reads = info_val.split(',')
                     break
            # msv use line number as sv ids
            read_id_dict[str(line_no)] = msv_reads
        else:
            raise Exception("not support yet")
    f.close()
    return read_id_dict


def parse_msvasm(msvasm):
    with gzip.open(msvasm) as gzip_file:
        read_ids = []
        for line in gzip_file:
            line = line.strip().split()
            is_bp = False
            if re.match(r"[><]", line[2].decode('utf-8')):
                is_bp = True
            col_info = 8 if is_bp else 6
            read_ids.append(line[col_info-3].decode('utf-8'))

    read_ids = set(read_ids)
    return read_ids


def parse_msvasm_withtype(msvasm):

    h = {}
    with gzip.open(msvasm, 'rt') as gzip_file:
        read_ids = []
        for line in gzip_file:
            t = line.strip().split('\t')
            is_bp = False
            if re.match(r"[><]", t[2]):
                is_bp = True

            col_info = 8 if is_bp else 6
            name = t[col_info-3]
            if name not in h:
                h[name] = []
            h[name].append(get_type(t, col_info))
    return h


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


def parse_svid(svid):
    svid = str(svid)
    # remove _1 and _2 for paired SV event
    # severus
    if svid.startswith('severus'):
        query_svid = re.sub("_[12]", "", svid)
    # savana
    elif svid.startswith("ID_"):
        query_svid = re.sub(r"(ID_\d+)_([1,2])", "\\1", svid)
    # nanomonsv
    elif svid.startswith("r_") or svid.startswith("i_") or svid.startswith("d_"):
        query_svid = re.sub(r"(r_\d+)_([0,1])", "\\1", svid)
    # sniffles2
    elif svid.startswith("Sniffles2"):
        query_svid = svid
    else:
        query_svid = svid
    return query_svid
        #raise Exception("not support yet")

def othercaller_filterasm(vcf_file, opt, readidtsv, msvasm, outstat, consensus_sv_ids):
    """
    vcf: input vcf, support Severus, Sniffles, nanomonsv, and savana
    consensus_sv_ids: consensus SV ids from minisv annot
    """
    # NOTE: why some germline SV disappear and no somatic SV
    parsed_id_dict = parse_readids(readidtsv)
    read_ids = parse_msvasm(msvasm)
    read_ids_withtype = parse_msvasm_withtype(msvasm)

    min_read_len = math.floor(opt.min_len * opt.read_len_ratio + 0.499)
    vcf = gc_parse_sv(vcf_file, min_read_len, opt.min_count, opt.ignore_flt, opt.check_gt)

    consensus_ids = []
    if consensus_sv_ids != "":
        with open(consensus_sv_ids) as inf:
            for line in inf:
                consensus_ids.append(re.sub("_[12]", "", line.strip()))

    is_consensus = True
    outf = open(outstat, 'w')
    outf.write("svlen\tsvlen_range\tsvid\tis_consensus\ttype\tcontig1\tpos1\tori\tcontig2\tpos2\tasm_support\tasm_support_onlyreadname\tread_name_number\n")

    # first vcf is to be annotated
    # read id, typing and length
    all_ids = []
    # only read name overlap
    all_ids_onlyname = []
    # msv sv id record
    for i in range(len(vcf)):
 
        t = vcf[i]

        # Not long enough for non-BND type
        # NOTE: here use 100bp as the default length instead of 80
        if t.SVTYPE != "BND" and abs(t.SVLEN) < opt.min_len:
            continue

        if t.SVTYPE == "BND" and t.ctg == t.ctg2 and abs(t.SVLEN) < opt.min_len:
            continue

        if t.count > 0 and t.count < opt.min_count:
            continue

        # filter by VAF
        if t.vaf is not None and t.vaf < opt.min_vaf:
            continue

        if opt.bed is not None:
            if t.ctg not in opt.bed or t.ctg2 not in opt.bed:
                continue
            if len(iit_overlap(opt.bed[t.ctg], t.pos, t.pos + 1)) == 0:
                continue
            if len(iit_overlap(opt.bed[t.ctg2], t.pos2, t.pos2 + 1)) == 0:
                continue

        query_svid = str(parse_svid(t.svid))
        
        if len(consensus_ids) > 0:
            is_consensus = query_svid in consensus_ids
        assert query_svid in parsed_id_dict

        ol_readn = set(parsed_id_dict[query_svid]) & read_ids

        read_num_wt_type = 0
        t_type_flag = simple_type(t.SVTYPE, t.ctg, t.ctg2)
        for read_i in ol_readn:
            assert read_i in read_ids_withtype

            # one read with multiple sv type
            # read_ids_withtype: asm sv results
            for sv_j in read_ids_withtype[read_i]:
                # translocation
                if t_type_flag == 8 and sv_j.flag == 8:
                    read_num_wt_type += 1
                    break

                # 1. >100k(cutoff) DUP/INS/DEL/INV can be translocation in self-assembly
                if t_type_flag in [1, 2, 4] and abs(t.SVLEN) >= 1e5 and sv_j.flag == 8: 
                    read_num_wt_type += 1
                    break

                len_check = (
                    abs(sv_j.len) >= abs(t.SVLEN) * opt.min_len_ratio
                    and abs(t.SVLEN) >= abs(sv_j.len) * opt.min_len_ratio
                ) or abs(abs(sv_j.len) - abs(t.SVLEN)) <= 1000

                # strict sv type and length check
                if (t_type_flag & sv_j.flag) and len_check:
                    read_num_wt_type += 1
                    break

                # patch complex BND such as template insertion
                # or self-assembly intra-chrom BND but nanomonsv reports INV
                if (t_type_flag == 0 or sv_j.flag == 0) and len_check:
                     read_num_wt_type += 1
                     break

                # 2. <1K INDEL could change type INS->DEL, DEL->INS
                #    partly caused by VNTR
                if t_type_flag in [1, 2] and sv_j.flag in [1, 2]:
                    if (t_type_flag == 1 and sv_j.flag == 2) or (t_type_flag == 2 and sv_j.flag == 1):
                        len_check1 = abs(t.SVLEN) + abs(sv_j.len) <= 1000
                        #ignore the SVTYPE for the filtering in this step
                        if len_check1:
                            read_num_wt_type += 1
                            break
                        else:
                            # duplicated check
                            len_check2 = (
                                abs(sv_j.len) >= abs(t.SVLEN) * opt.min_len_ratio
                                and abs(t.SVLEN) >= abs(sv_j.len) * opt.min_len_ratio
                            )
                            if len_check2 and (t_type_flag & sv_j.flag):
                                read_num_wt_type += 1
                                break

                #debug INS
                #if t.svid == 'severus_INS16922':
                #     print(t.svid, read_i, sv_j.flag, sv_j.len, t.SVLEN, t.SVTYPE, len_check, abs(t.SVLEN) * opt.min_len_ratio, abs(sv_j.len) * opt.min_len_ratio, opt.min_len_ratio)
                #if t.svid == 'severus_INS12685':
                #     print(t.svid, read_i, sv_j.flag, sv_j.len, t.SVLEN, t.SVTYPE, len_check, abs(t.SVLEN) * opt.min_len_ratio, abs(sv_j.len) * opt.min_len_ratio, opt.min_len_ratio)
                #if t.svid == 'severus_INS17052':
                #     print(t.svid, read_i, sv_j.flag, sv_j.len, t.SVLEN, t.SVTYPE, len_check, abs(t.SVLEN) * opt.min_len_ratio, abs(sv_j.len) * opt.min_len_ratio, opt.min_len_ratio)
                #if t.svid == 'severus_INS21634':
                #     print(t.svid, read_i, sv_j.flag, sv_j.len, t.SVLEN, sv_j.type, t.SVTYPE, len_check, abs(t.SVLEN) * opt.min_len_ratio, abs(sv_j.len) * opt.min_len_ratio, opt.min_len_ratio)
                #if t.svid == 'severus_BND225_1':
                #     print(t.svid, read_i, sv_j.flag, t_type_flag, sv_j.len, t.SVLEN, sv_j.type, t.SVTYPE, len_check, abs(t.SVLEN) * opt.min_len_ratio, abs(sv_j.len) * opt.min_len_ratio, opt.min_len_ratio)
                #if t.svid == 'r_319_0':
                #     print(t.svid, read_i, sv_j.flag, t_type_flag, sv_j.len, t.SVLEN, sv_j.type, t.SVTYPE, len_check, abs(t.SVLEN) * opt.min_len_ratio, abs(sv_j.len) * opt.min_len_ratio, opt.min_len_ratio, read_num_wt_type)
                #if t.svid == 'severus_INS21583':
                #if t.svid == 'severus_INS21634':
                #if t.svid == 'severus_INS13454':
                #if t.svid == 'severus_INS13602':
                #if t.svid == 'i_1179':
                #     print(t.svid, read_i, sv_j.flag, t_type_flag, sv_j.len, t.SVLEN, sv_j.type, t.SVTYPE, len_check, abs(t.SVLEN) * opt.min_len_ratio, abs(sv_j.len) * opt.min_len_ratio, opt.min_len_ratio, read_num_wt_type)

                #severus_INS10090 m84039_230414_235240_s2/244716513/ccs 8 1 0 1760 BND INS False 1056.0 0.0 0.6 0
                if t.svid == opt.svid:
                    print(t.svid, read_i, sv_j.flag, t_type_flag, sv_j.len, t.SVLEN, sv_j.type, t.SVTYPE, len_check, abs(t.SVLEN) * opt.min_len_ratio, abs(sv_j.len) * opt.min_len_ratio, opt.min_len_ratio, read_num_wt_type)


        filt_cnt = len(ol_readn)

        if opt.print_all:
            outf.write('\t'.join(map(str, [t.SVLEN, classify_sv_len(abs(t.SVLEN)), t.svid, is_consensus, t.SVTYPE, t.ctg, t.pos, t.ori, t.ctg2, t.pos2, 
                       read_num_wt_type, filt_cnt, t.count]))+'\n')

         # output sv ids in vcf or msv
        if read_num_wt_type >= opt.min_count and filt_cnt >= opt.min_count and t.count >= opt.min_count:
            all_ids.append(query_svid)
        if opt.only_readname:
            if filt_cnt >= opt.min_count and t.count >= opt.min_count:
                all_ids_onlyname.append(query_svid)

    outf.close()

    assert set(all_ids).issubset(set(parsed_id_dict.keys()))
    assert set(all_ids_onlyname).issubset(set(parsed_id_dict.keys()))
    ## output vcf format with sv ids above
    if is_gzipped(vcf_file):
        f = gzip.open(vcf_file, 'rt')
    else:
        f = open(vcf_file)
    line_no = 0
    for line in f:
        line_no += 1
        if line[0] == "#":
            print(line.strip())
            continue
        if 'avg_mapq' in line: # MSV
            query_svid = str(parse_svid(line_no))
        else: # other caller
            t = line.strip().split("\t")
            query_svid = str(parse_svid(t[2]))
        if opt.only_readname:
            if query_svid in all_ids_onlyname:
                print(line.strip())
        else:
            if query_svid in all_ids:
                print(line.strip())
    f.close()


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
