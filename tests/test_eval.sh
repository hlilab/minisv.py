#!/bin/bash -ex


test1() {
minisv eval -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -a /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38/somatic_SVs/severus_somatic.vcf ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf > test1.eval

/hlilab/alvin/miniconda3/bin/k8 ../../minisv/minisv.js eval -b ~/data/pangenome_sv_benchmarking/minisv/data/hs38.reg.bed -a /homes6/hli/hli1/gafcall/COLO829.truth.hs38.vcf ../../1a.alignment_sv_tools/output/severus/COLO829_hifi1/grch38/somatic_SVs/severus_somatic.vcf ../../1a.alignment_sv_tools/output/nanomonsv/COLO829_hifi1/grch38_tnpair.vcf ../../1a.alignment_sv_tools/output/sniffles/COLO829_hifi1/grch38_multi.vcf.gz > test2.eval
wc -l test1.eval test2.eval
diff test1.eval test2.eval | wc -l
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
