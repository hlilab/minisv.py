#!/bin/bash -ex


test1() {
#/hlilab/alvin/miniconda3/envs/gafcall/bin/k8 /hlilab/alvin/miniconda3/envs/gafcall/bin/minisv.js annot -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -c 2 -l 100 -p ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf ../../1a.alignment_sv_tools/output/savana12/COLO829_hifi1/grch38/grch38_T_tag.classified.somatic.vcf ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf | awk '$9>=2' | cut -f 8 > colo829.consensusid
#minisv filterasm -m 0.4 -c 2 -i colo829.consensusid -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38_cutoff2_read_ids/read_ids.csv /hlilab/hli/gafcall/pair_v2/COLO829T.self.Q0.gsv.gz test.out ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf

#/hlilab/alvin/miniconda3/envs/gafcall/bin/k8 /hlilab/alvin/miniconda3/envs/gafcall/bin/minisv.js annot -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -c 2 -l 100 -p ../../1a.alignment_sv_tools/output/severus/HCC1395_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf ../../1a.alignment_sv_tools/output/savana12/HCC1395_hifi1/grch38/grch38_T_tag.classified.somatic.vcf ../../1a.alignment_sv_tools/output/nanomonsv/HCC1395_hifi1/grch38_tnpair.vcf | awk '$9>=2' | cut -f 8 > colo829.consensusid
#minisv filterasm -c 2 --svid severus_INS10090 -i colo829.consensusid -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/severus/HCC1395_hifi1/grch38_cutoff2_read_ids/read_ids.csv /hlilab/hli/gafcall/pair_v2/HCC1395T.self.Q0.gsv.gz test.out ../../1a.alignment_sv_tools/output/severus/HCC1395_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf 
#minisv filterasm -r -c 2 --svid severus_INS10090 -i colo829.consensusid -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/severus/HCC1395_hifi1/grch38_cutoff2_read_ids/read_ids.csv /hlilab/hli/gafcall/pair_v2/HCC1395T.self.Q0.gsv.gz test.out ../../1a.alignment_sv_tools/output/severus/HCC1395_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf 
#minisv filterasm -c 2 --svid severus_DEL2309 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../1a.alignment_sv_tools/output/severus_latest/COLO829_hifi1/grch38_cutoff2_read_ids/read_ids.csv /hlilab/hli/gafcall/pair_v2/COLO829T.self.Q0.gsv.gz test.out ../1a.alignment_sv_tools/output/severus_latest/COLO829_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf
echo "------"
#minisv filterasm -c 2 --svid severus_DEL2309 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../1a.alignment_sv_tools/output/severus_latest/COLO829_hifi1/grch38_cutoff2_read_ids/read_ids.csv COLO829_hifi1_all_filtered_whole_preset1_wt_snf2/denovo.gsv.gz test.out2 ../1a.alignment_sv_tools/output/severus_latest/COLO829_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf

minisv filterasm -c 2 --svid Sniffles2.DEL.442MB -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed test_snf2.vcf COLO829_hifi1_all_filtered_whole_preset1_wt_snf2/denovo.gsv.gz test.out2 test_snf2.vcf

#/hlilab/alvin/miniconda3/envs/gafcall/bin/k8 /hlilab/alvin/miniconda3/envs/gafcall/bin/minisv.js annot -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -c 2 -l 100 -p ../../1a.alignment_sv_tools/output/severus/NCI1437_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf ../../1a.alignment_sv_tools/output/savana12/NCI1437_hifi1/grch38/grch38_T_tag.classified.somatic.vcf ../../1a.alignment_sv_tools/output/nanomonsv/NCI1437_hifi1/grch38_tnpair.vcf | awk '$9>=2' | cut -f 8 > colo829.consensusid
#minisv filterasm -c 2 -i colo829.consensusid -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/severus/NCI1437_hifi1/grch38_cutoff2_read_ids/read_ids.csv /hlilab/hli/gafcall/pair_v2/NCI1437T.self.Q0.gsv.gz test.out ../../1a.alignment_sv_tools/output/severus/NCI1437_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf

#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/nanomonsv/NCI1437_hifi1/T/grch38_parse.nanomonsv.supporting_read.txt /hlilab/hli/gafcall/pair_v2/NCI1437T.self.Q0.gsv.gz test.out ../../1a.alignment_sv_tools/output/nanomonsv/NCI1437_hifi1/grch38_tnpair.vcf 

#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/nanomonsv/HCC1395_hifi1/T/grch38_parse.nanomonsv.supporting_read.txt /hlilab/hli/gafcall/pair_v2/HCC1395T.self.Q0.gsv.gz test.out ../../1a.alignment_sv_tools/output/nanomonsv/HCC1395_hifi1/grch38_tnpair.vcf 

#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /hlilab/hli/gafcall/pair_v2/msv/COLO829T.hg38l+tg.pair-c2s1.msv /hlilab/hli/gafcall/pair_v2/COLO829T.self.Q0.gsv.gz test.out /hlilab/hli/gafcall/pair_v2/msv/COLO829T.hg38l+tg.pair-c2s1.msv

##TUMOUR_SUPPORTING_READS
#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/savana12/COLO829_hifi1/grch38/grch38_T_tag.sv_breakpoints_read_support.tsv /hlilab/hli/gafcall/pair_v2/COLO829T.self.Q0.gsv.gz  test.out ../../1a.alignment_sv_tools/output/savana12/COLO829_hifi1/grch38/grch38_T_tag.classified.somatic.vcf 

#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/severus/HCC1937_hifi1/grch38_cutoff2_read_ids/read_ids.csv /hlilab/hli/gafcall/pair_v2/HCC1937T.self.Q0.gsv.gz  test.out ../../1a.alignment_sv_tools/output/severus/HCC1937_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf | wc -l
#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/savana12/HCC1937_hifi1/grch38/grch38_T_tag.sv_breakpoints_read_support.tsv /hlilab/hli/gafcall/pair_v2/HCC1937T.self.Q0.gsv.gz  test.out ../../1a.alignment_sv_tools/output/savana12/HCC1937_hifi1/grch38/grch38_T_tag.classified.somatic.vcf | wc -l

#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/T/grch38_parse.nanomonsv.supporting_read.txt /hlilab/hli/gafcall/pair_v2/COLO829T.self.Q0.gsv.gz test.out ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf

#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/T/grch38_parse.nanomonsv.supporting_read.txt /hlilab/hli/gafcall/pair_v2/COLO829T.self.Q0.gsv.gz test.out ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf

#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/sniffles/COLO829_hifi1/grch38_multi.vcf.gz /hlilab/hli/gafcall/pair_v2/COLO829T.self.Q0.gsv.gz test.out ../../1a.alignment_sv_tools/output/sniffles/COLO829_hifi1/grch38_multi.vcf.gz

#minisv snfpair -t 1 -n 2 ../../1a.alignment_sv_tools/output/sniffles/COLO829_hifi1/grch38_multi.vcf.gz > test1.vcf
#/hlilab/alvin/miniconda3/bin/k8 ../../minisv/minisv.js snfpair -t 1 -n 2 ../../1a.alignment_sv_tools/output/sniffles/COLO829_hifi1/grch38_multi.vcf.gz > test2.vcf
#diff test1 test2 | wc -l
#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed test1.vcf /hlilab/hli/gafcall/pair_v2/COLO829T.self.Q0.gsv.gz test.out test1.vcf > test1.filterasm

#minisv filterasm -c 2 -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed ../../1a.alignment_sv_tools/output/savana12/HCC1937_hifi1/grch38/grch38_T_tag.sv_breakpoints_read_support.tsv /hlilab/hli/gafcall/pair_v2/HCC1937T.self.Q0.gsv.gz  test.out ../../1a.alignment_sv_tools/output/savana12/HCC1937_hifi1/grch38/grch38_T_tag.classified.somatic.vcf | wc -l

}


main() {
test1
}


main
