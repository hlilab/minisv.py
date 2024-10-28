#!/bin/bash -ex

test1() {
minisv isec -w 1000 pair_v2/COLO829BL.hg38l.Q5.gsv.gz  pair_v2/COLO829BL.hg38g.Q20.gsv.gz > test_isec.gsv
/hlilab/alvin/miniconda3/bin/k8 ../../minisv_js_latest/minisv.js isec -w 1000 pair_v2/COLO829BL.hg38l.Q5.gsv.gz  pair_v2/COLO829BL.hg38g.Q20.gsv.gz > test_isec_orig.gsv

diff test_isec.gsv test_isec_orig.gsv | wc -l
wc -l test_isec_orig.gsv test_isec.gsv
}

test2() {
minisv isec -w 1000 pair_v2/COLO829BL.hg38l.Q5.gsv.gz pair_v2/COLO829BL.chm13g.Q20.gsv.gz pair_v2/COLO829BL.chm13l.Q20.gsv.gz > test_isec.gsv
/hlilab/alvin/miniconda3/bin/k8 ../../minisv_js_latest/minisv.js isec -w 1000 pair_v2/COLO829BL.hg38l.Q5.gsv.gz  pair_v2/COLO829BL.chm13g.Q20.gsv.gz pair_v2/COLO829BL.chm13l.Q20.gsv.gz > test_isec_orig.gsv

diff test_isec.gsv test_isec_orig.gsv | wc -l
wc -l test_isec_orig.gsv test_isec.gsv
}

test3() {
minisv isec -w 1000 pair_v2/COLO829BL.hg38l.Q5.gsv.gz pair_v2/COLO829BL.chm13g.Q20.gsv.gz pair_v2/COLO829BL.chm13l.Q20.gsv.gz /hlilab/hli/gafcall/normal_v2/COLO829BL.self.Q0.gsv.gz > test_isec.gsv
/hlilab/alvin/miniconda3/bin/k8 ../../minisv_js_latest/minisv.js isec -w 1000 pair_v2/COLO829BL.hg38l.Q5.gsv.gz  pair_v2/COLO829BL.chm13g.Q20.gsv.gz pair_v2/COLO829BL.chm13l.Q20.gsv.gz /hlilab/hli/gafcall/normal_v2/COLO829BL.self.Q0.gsv.gz > test_isec_orig.gsv

diff test_isec.gsv test_isec_orig.gsv | wc -l
wc -l test_isec_orig.gsv test_isec.gsv
}

main() {
test1
test2
test3
}

main
