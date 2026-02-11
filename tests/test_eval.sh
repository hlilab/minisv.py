#!/bin/bash -ex

export PYTHONPATH=/hlilab/alvin/pangenome_sv_benchmarking/minisv.py/minimap2/:${PYTHONPATH}

test1() {
    #minisv eval -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38/somatic_SVs/severus_somatic.vcf  > test1.eval
    #/hlilab/alvin/miniconda3/bin/k8 ../../minisv/minisv.js eval -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38/somatic_SVs/severus_somatic.vcf > test2.eval
    #wc -l test1.eval test2.eval
    #diff test1.eval test2.eval | wc -l

    #minisv eval -s -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/nanomonsv_latest/COLO829_hifi1/grch38_tnpair.vcf > nanomonsv_sv_size_test1.eval
    #minisv eval -s -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/severus_latest/COLO829_hifi1/grch38_cutoff2_read_ids/somatic_SVs/severus_somatic.vcf > severus_sv_size_test1.eval
    #minisv eval -s -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/savana13/COLO829_hifi1/grch38/grch38_T_tag.classified.somatic.vcf > savana_sv_size_test1.eval

    minisv eval -s -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../COLO829_hifi1_all_filtered_whole_preset1_wt_snf2/nanomonsv_l+s_3_filtered.vcf > nanomonsv_sv_size_test1.eval_asm
    minisv eval -s -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../COLO829_hifi1_all_filtered_whole_preset1_wt_snf2/severus_l+s_3_filtered.vcf > severus_sv_size_test1.eval_asm
    minisv eval -s -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../COLO829_hifi1_all_filtered_whole_preset1_wt_snf2/savana_l+s_3_filtered.vcf > savana_sv_size_test1.eval_asm

    #/hlilab/alvin/miniconda3/bin/k8 ../../minisv/minisv.js eval -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf  > test2.eval
    #wc -l test1.eval test2.eval
    #diff test1.eval test2.eval | wc -l

    #minisv eval -M -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38/somatic_SVs/severus_somatic.vcf ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf > test1.eval
    #/hlilab/alvin/miniconda3/bin/k8 ../../minisv/minisv.js eval -M -a -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38/somatic_SVs/severus_somatic.vcf ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf  > test2.eval
    #wc -l test1.eval test2.eval
    #diff test1.eval test2.eval | wc -l
}

test2() {
    minisv eval -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -M -a /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38/somatic_SVs/severus_somatic.vcf ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf ../../1a.alignment_sv_tools/output/savana12/COLO829_hifi1/grch38/grch38_T_tag.classified.somatic.vcf > test1.eval
    /hlilab/alvin/miniconda3/bin/k8 ../../minisv/minisv.js eval -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -M -a /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38/somatic_SVs/severus_somatic.vcf ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf ../../1a.alignment_sv_tools/output/savana12/COLO829_hifi1/grch38/grch38_T_tag.classified.somatic.vcf > test2.eval
    wc -l test1.eval test2.eval
    diff test1.eval test2.eval | wc -l
}

main() {
    test1
    #test2
}

main

