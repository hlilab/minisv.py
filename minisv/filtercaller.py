import gzip
import pandas as pd
import csv
import time
import os
import mappy as mp
from .regex import re_info
import math
import warnings
import re
from .util import is_gzipped, is_vcf
from .eval import gc_parse_sv, iit_overlap
from .type import get_type, simple_type
from pathlib import Path
import subprocess
from minisv.read_parser import load_reads
from .ensemble import insilico_truth


# The categories to group by
categories = ["l", "l+s", "l+g+s", "l+g"]
callers = ["severus", "savana", "nanomonsv"]


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
        assert query_svid in parsed_id_dict, f"{query_svid}, {vcf_file}"

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



def hit_to_paf(read_name, read_seq, hit):
    """
    Convert a mappy.Alignment object into a PAF line.
    """

    # Required PAF fields
    qname = read_name
    qlen = len(read_seq)
    qstart = hit.q_st
    qend = hit.q_en

    rname = hit.ctg
    rlen = hit.ctg_len
    rstart = hit.r_st
    rend = hit.r_en

    strand = "+" if hit.strand == 1 else "-"

    # Optional: minimap2-style tags
    tp = "P" if hit.is_primary else "S"       # primary or secondary
    nm = hit.NM if hit.NM is not None else 0
    ms = hit.mlen                            # number of matching bases
    bl = hit.blen                            # alignment block length
    mapq = hit.mapq
    cigar = hit.cigar_str or "*"
    md = hit.MD or ""
    cs = hit.cs or ""
    ds = hit.ds or ""

    # Build PAF line
    fields = [
        qname, qlen, qstart, qend,
        strand,
        rname, rlen, rstart, rend,
        ms, bl, mapq,
        f"NM:i:{nm}",
        f"ms:i:{ms}",
        f"AS:i:{ms - nm}",          # crude score estimate
        f"nn:i:0",
        f"tp:A:{tp}",
        f"cg:Z:{cigar}",
    ]

    if md:
        fields.append(f"MD:Z:{md}")
    if cs:
        fields.append(f"cs:Z:{cs}")
    if ds:
        fields.append(f"ds:Z:{ds}")

    return "\t".join(map(str, fields))


import time
import csv
import psutil
import os
import gc
from pathlib import Path

class Timer:
    """Context manager to track time and max memory usage per step"""
    def __init__(self, name, timings_list):
        self.name = name
        self.timings_list = timings_list
        self.process = psutil.Process(os.getpid())

    def __enter__(self):
        gc.collect()  # Clean up before measuring
        self.start_time = time.time()
        self.start_mem = self.process.memory_info().rss / 1024 / 1024  # MB
        self.max_mem = self.start_mem
        return self

    def __exit__(self, *args):
        gc.collect()
        elapsed = time.time() - self.start_time
        current_mem = self.process.memory_info().rss / 1024 / 1024  # MB

        # Update max memory if current is higher
        self.max_mem = max(self.max_mem, current_mem)

        record = {
            "step": self.name,
            "time_sec": round(elapsed, 3),
            "max_memory_mb": round(self.max_mem, 2)
        }
        self.timings_list.append(record)

        print(f"[TIMING] {self.name}: {elapsed:.2f}s | Max Memory: {self.max_mem:.1f} MB")


def isec(w, gsvs, file_handler=None):
    """Usage: minisv isec base.gsv alt.gsv [...]"""
    from .type import get_type

    g = []
    # file as filter
    for gsv in gsvs[1:]:
        h = {}
        with gzip.open(gsv, 'rt') as gsv_file:
            for line in gsv_file:
                t = line.strip().split("\t")
                if re.match(r"[><]", t[2]):
                    col_info = 8
                else:
                    col_info = 6
                name = t[col_info - 3]
                if name not in h:
                    h[name] = []
                h[name].append(get_type(t, col_info))
        g.append(h)
    # one read may contain somatic and germline sv
    # g: [{readname:[sv_dict]}]
    out = file_handler if file_handler is not None else None

    # file to be filtered
    with gzip.open(gsvs[0], 'rt') as base_gsv:
        for line in base_gsv:
            t = line.strip().split("\t")
            if re.match(r"[><]", t[2]):
                col_info = 8
            else:
                col_info = 6
            name = t[col_info - 3]
            x = get_type(t, col_info)
            n_found = 0

            for h in g:
                if name not in h:
                    break

                a = h[name]
                found = False
                for j in range(len(a)):
                    if x.st - w < a[j].en and a[j].st < x.en + w:
                        # NOTE: why a[i].flag & 8?? what if base sv is not translocation
                        if x.flag == 0 or (x.flag & a[j].flag) or (a[j].flag & 8):
                            found = True
                if not found:
                    break
                n_found += 1

            if n_found == len(g):
                print(line.strip(), file=out)


class MinisvReads:
    def __init__(self, som_vcfs, readid_tsvs, bam_path, ref, hap1_denovo_ref_path, hap2_denovo_ref_path, 
                 graph_ref_path = "/hlilab/hli/minigraph/HPRC-r2/CHM13-464.gfa.gz", work_dir="minisv_work", filtered_readcount_cutoff=2, platform='hifi'):
        self.som_vcfs = som_vcfs
        self.readid_tsvs = readid_tsvs

        self.bam_path = Path(bam_path)

        self.ref = Path(ref)
        self.hap1_denovo_ref_path = Path(hap1_denovo_ref_path)
        self.hap2_denovo_ref_path = Path(hap2_denovo_ref_path)
        self.graph_ref_path = Path(graph_ref_path)

        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.filtered_readcount_cutoff = filtered_readcount_cutoff
        self.platform = platform
        self.timings = []

    def extract_read_ids(self, min_read_len, min_count, ignore_flt, check_gt):
        """ only extract somatic SV read ids """
        with Timer("extract_read_ids", self.timings):
            self.read_ids_file = self.work_dir / "readid.names"
            som_read_ids = set()

            self.read_id_dict = {}
            for (readid_tsv, vcf) in zip(self.readid_tsvs, self.som_vcfs):
                assert os.path.exists(readid_tsv), f"{readid_tsv} not exists"
                assert os.path.exists(vcf), f"{vcf} not exists"

                self.read_id_dict |= parse_readids(readid_tsv)
                vcf = gc_parse_sv(vcf, min_read_len, min_count, ignore_flt, check_gt)

                for i in range(len(vcf)):
                    t = vcf[i]
                    query_svid = str(parse_svid(t.svid))
                    assert query_svid in self.read_id_dict, 'somatic sv not in read id table'
                    som_read_ids.add(query_svid)

            read_names = set()
            for s in list(som_read_ids):
                read_names |= set(self.read_id_dict[s])
       
            with open(self.read_ids_file, "w") as fin:
                for i in sorted(list(read_names)):
                    fin.write(f"{i}\n")

   
    def extract_reads(self):
        """samtools fastq aln.bam | seqtk subseq - read-names.txt > reads.fq"""
        with Timer("extract_reads", self.timings):
            self.fastq_out = self.work_dir / "som_reads.fq.gz"
            if os.path.exists(self.fastq_out):
                return

            cmd = f"samtools fastq {self.bam_path} | seqtk subseq - {self.read_ids_file} | gzip  > {self.fastq_out}"
            subprocess.run(cmd, shell=True, check=True)
            print(f"Reads extracted to {self.fastq_out}")
            return

    def build_pooled_reference(self, out_fa="pooled_reference.fa"):
        with Timer("build_pooled_reference", self.timings):
            self.pooled_ref = self.work_dir / out_fa
            with open(self.pooled_ref, "w") as out:
                for ref in [self.hap1_denovo_ref_path, self.hap2_denovo_ref_path]:
                    with gzip.open(ref, "rt") as f:
                        out.write(f.read())
            return mp.Aligner(str(self.pooled_ref)) # preset??

    def align_reads_to_self(self, paf='denovo_aligned.paf.gz'):
        with Timer("align_reads_to_denovo", self.timings):
            self.paf_out = self.work_dir / paf
            if os.path.exists(self.paf_out):
                return

            #aligner = self.build_pooled_reference()
            #if not aligner:
            #    raise Exception("ERROR: failed to load/build index")
    
            #f = gzip.open(self.paf_out, 'wt')
            #alignments = []
            #for name, seq, qual in mp.fastx_read(str(self.fastq_out)):
            #    for hit in aligner.map(seq, MD=False, cs=False, ds=True):
            #        paf_line = hit_to_paf(name, seq, hit)
            #        f.write(f"{paf_line}\n")
            #        alignments.append(paf_line)
            #f.close()
            #return alignments

            if self.platform == 'hifi':
                cmd = f"../../1a.alignment_sv_tools/minimap2/minimap2 --ds -t 4 -cx map-hifi -s50 <(zcat {self.hap1_denovo_ref_path} {self.hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            else:
                cmd = f"../../1a.alignment_sv_tools/minimap2/minimap2 --ds -t 4 -cx lr:hq <(zcat {self.hap1_denovo_ref_path} {self.hap2_denovo_ref_path}) -I100g --secondary=no {self.fastq_out} | gzip - > {self.paf_out}" 
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"aligned to self assembly")

    def build_reference(self):
        return mp.Aligner(str(self.ref))

    def align_reads_to_grch38(self, paf='grch38_aligned.paf.gz'):
        with Timer("align_reads_to_grch38", self.timings):
            self.grch38_paf_out = self.work_dir / paf
            #aligner = self.build_reference()

            #if os.path.exists(self.grch38_paf_out):
            #    return
            #if not aligner:
            #    raise Exception("ERROR: failed to load/build index")
    
            #f = gzip.open(self.grch38_paf_out, 'wt')
            #alignments = []
            #for name, seq, qual in mp.fastx_read(str(self.fastq_out)):
            #    for hit in aligner.map(seq, MD=False, cs=False, ds=True):
            #        paf_line = hit_to_paf(name, seq, hit)
            #        f.write(f"{paf_line}\n")
            #        alignments.append(paf_line)
            #f.close()
            #return alignments

            if self.platform == 'hifi':
                cmd = f"../../1a.alignment_sv_tools/minimap2/minimap2 --ds -t 4 -cx map-hifi -s50 {self.ref} {self.fastq_out} | gzip - > {self.grch38_paf_out}" 
            else:
                cmd = f"../../1a.alignment_sv_tools/minimap2/minimap2 --ds -t 4 -cx lr:hq {self.ref} {self.fastq_out} | gzip - > {self.grch38_paf_out}" 
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"aligned to grch38")

    def align_reads_to_graph(self, gaf='denovo_aligned.gaf.gz'):
        with Timer("align_reads_to_graph", self.timings):
            self.gaf_out = self.work_dir / gaf
            if os.path.exists(self.gaf_out):
                return
            cmd = f"../../1a.alignment_sv_tools/minigraph/minigraph -cxlr -t 4 {self.graph_ref_path} {self.fastq_out} | gzip - > {self.gaf_out}"
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            print(f"Reads aligned to {self.graph_ref_path} to generate {self.gaf_out}")
            return

    def parse_raw_sv_grch38(self, opt, gsv='grch38.gsv.gz'):
        with Timer("parse_raw_sv_grch38", self.timings):
            self.grch38_gsv_out = self.work_dir / gsv

            f = gzip.open(self.grch38_gsv_out, 'wt')
            load_reads(str(self.grch38_paf_out), opt, f)
            f.close()
            return

    def parse_raw_sv_self(self, opt, gsv='denovo.gsv.gz'):
        with Timer("parse_raw_sv_denovo", self.timings):
            self.denovo_gsv_out = self.work_dir / gsv

            f = gzip.open(self.denovo_gsv_out, 'wt')
            load_reads(str(self.paf_out), opt, f)
            f.close()
            return

    def parse_raw_sv_graph(self, opt, gsv='graph.gsv.gz'):
        with Timer("parse_raw_sv_graph", self.timings):
            self.graph_gsv_out = self.work_dir / gsv

            f = gzip.open(self.graph_gsv_out, 'wt')
            load_reads(str(self.gaf_out), opt, f)
            f.close()
            return

    def isec_g(self, opt, msv='l+g.gsv.gz'):
        # l + g
        with Timer("isec_graph", self.timings):
            self.isec_graph_out = self.work_dir / msv

            f = gzip.open(self.isec_graph_out, 'wt')
            isec(1000, [self.grch38_gsv_out, self.graph_gsv_out], f)
            f.close()
            return

    def isec_s(self, opt, msv='l+s.gsv.gz'):
        # l+s
        with Timer("isec_graph", self.timings):
            self.isec_denovo_out = self.work_dir / msv

            f = gzip.open(self.isec_denovo_out, 'wt')
            isec(1000, [self.grch38_gsv_out, self.denovo_gsv_out], f)
            f.close()
            return

    def isec_gs(self, opt, msv='l+g+s.gsv.gz'):
        # l+s
        with Timer("isec_gs", self.timings):
            self.isec_gs_out = self.work_dir / msv

            f = gzip.open(self.isec_gs_out, 'wt')
            isec(1000, [self.grch38_gsv_out, self.graph_gsv_out, self.denovo_gsv_out], f)
            f.close()
            return

    def parse_ids_from_gsv(self, msv):
        with gzip.open(self.work_dir / msv, 'rt') as fin:
            for line in fin:
                t = line.strip().split('\t')

                is_bp = False
                if re.match(r"[><]", t[2]):
                    is_bp = True

                col_info = 8 if is_bp else 6
                name = t[col_info-3]
                yield name

    def export_filtered_stat(self):
        """
        ##        l+s   l+g   l+g+s   severus savana nanomonsv
        #read1
        #read2
        #read3
        """
        self.filtered_stat_read = self.work_dir / Path("read_stat.tsv")
        self.filtered_stat_sv = self.work_dir / Path("sv_stat.tsv")

        self.filtered_sv_records = [self.work_dir / Path("severus.tsv"), self.work_dir / Path("savana.tsv"), self.work_dir / Path("nanomonsv.tsv")]

        with Timer("export_filtered_stat", self.timings):
            read_names = []
            ls = list(self.parse_ids_from_gsv(Path("l+s.gsv.gz")))
            lgs = list(self.parse_ids_from_gsv(Path("l+g+s.gsv.gz")))
            lg = list(self.parse_ids_from_gsv(Path("l+g.gsv.gz")))

            with open(self.read_ids_file) as inf:
                for line in inf:
                    read_names.append(line.strip())
            df = pd.DataFrame({"read_name": read_names}, index=read_names)

            severus_status = []
            savana_status = []
            nanomonsv_status = []

            for (readid_tsv, vcf, status) in zip(self.readid_tsvs, self.som_vcfs, [severus_status, savana_status, nanomonsv_status]):
                read_to_sv_dict = {}
                parsed_id_dict = parse_readids(readid_tsv)
                for k in parsed_id_dict:
                    for r in parsed_id_dict[k]:
                        # svid -> read name
                        read_to_sv_dict[r] = k
                for r in read_names:
                    status.append(read_to_sv_dict.get(r, 'other_caller'))

            df.loc[:, 'l'] = True
            df.loc[:, 'l+s'] = df.index.isin(ls)
            df.loc[:, 'l+g+s'] = df.index.isin(lgs)
            df.loc[:, 'l+g'] = df.index.isin(lg)
            df.loc[:, 'severus'] = severus_status
            df.loc[:, 'savana'] = savana_status
            df.loc[:, 'nanomonsv'] = nanomonsv_status

        df.to_csv(str(self.filtered_stat_read), sep='\t')

        print(df.head())
        #savana  l+s  l+g+s  l+g
        #0            ID_11080    0      0    0
        #1             ID_1277    0      0    0
        #2            ID_13632    8      8    8
        #3            ID_15396    9      9    9
        #4            ID_15633    0      0    0

        sv_stats = []
        for caller, sv_record in zip(callers, self.filtered_sv_records):
            df_filtered = df.loc[:, categories+[caller]]
            df_filtered = df_filtered.loc[df_filtered[caller]!='other_caller', :]
            df_filtered = df_filtered.groupby(caller).sum()
            df_filtered.to_csv(sv_record, sep='\t')
            df_filtered_stat = (df_filtered > self.filtered_readcount_cutoff).sum(axis=0)
            sv_stats.append(df_filtered_stat)

        sv_stats = pd.concat(sv_stats, axis=1)
        sv_stats.columns = callers
        sv_stats.to_csv(str(self.filtered_stat_sv), sep='\t')
        ##            raw_sv l+s l+g l+g+s
        # severus
        # nanomonsv
        # savana

    def apply_filter_to_vcf(self):
        # l+s
        with Timer("apply_filter_to_vcf", self.timings):
            self.filtered_sv_records = [self.work_dir / Path("severus.tsv"), self.work_dir / Path("savana.tsv"), self.work_dir / Path("nanomonsv.tsv")]
            for (readid_tsv, vcf, filtered_sv_record, caller) in zip(self.readid_tsvs, self.som_vcfs, self.filtered_sv_records, callers):
                df = pd.read_table(filtered_sv_record, index_col=0)
                df = df > self.filtered_readcount_cutoff
                for cat in categories:
                    svids = df.index[df.loc[:, cat].values]
                    ## output vcf format with sv ids above
                    if is_gzipped(vcf):
                        f = gzip.open(vcf, 'rt')
                        suffix = '.vcf.gz'
                    else:
                        f = open(vcf)
                        suffix = '.vcf'

                    fout = open(self.work_dir / Path(f"{caller}_{cat}_filtered.vcf"), 'w')
                    for line in f:
                        if line[0] == "#":
                            print(line.strip(), file=fout)
                            continue
                        t = line.strip().split("\t")
                        query_svid = str(parse_svid(t[2]))
                        if query_svid in svids:
                            print(line.strip(), file=fout)
                    f.close()
                    fout.close()
            return

    def union_filtered_vcf(self, read_min_len, opt):
        from .union import union_sv
        for cat in categories:
            filtered_vcfs = [ str(self.work_dir / Path(f"{caller}_{cat}_filtered.vcf")) for caller in callers ]

            opt.print_sv = True
            f = open(self.work_dir / Path(f"{cat}_union.msv"), 'w')
            union_sv(filtered_vcfs, read_min_len, opt, file_handler=f)
            f.close()

            opt.print_sv = False
            f = open(self.work_dir / Path(f"{cat}_union_stat.msv"), 'w')
            union_sv(filtered_vcfs, read_min_len, opt, file_handler=f)
            f.close()

            f = open(self.work_dir / Path(f"{cat}_union_dedup.msv"), 'w')
            insilico_truth(str(self.work_dir / Path(f"{cat}_union.msv")), f)
            f.close()

    def save_timings(self, tsv_path=None):
        """Save collected timings to a TSV file"""
        if not self.timings:
            print("[WARNING] No timing data collected.")
            return

        if tsv_path is None:
            tsv_path = self.work_dir / "minisv_timings.tsv"

        fieldnames = ["step", "time_sec", "max_memory_mb"]

        with open(tsv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for record in self.timings:
                writer.writerow(record)

        total_time = sum(r["time_sec"] for r in self.timings)
        overall_max_mem = max(r["max_memory_mb"] for r in self.timings)

        print(f"\n[SUMMARY] Timings saved to: {tsv_path}")
        print(f"   Total time: {total_time:.2f} seconds")
        print(f"   Highest memory used: {overall_max_mem:.1f} MB")

## new SV filter from bam -> fastq reads -> gsv

