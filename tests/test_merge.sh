#!/bin/bash -ex


test1() {
##minisv isec -w 1000 pair_v2/COLO829BL.hg38l.Q5.gsv.gz pair_v2/COLO829BL.chm13g.Q20.gsv.gz pair_v2/COLO829BL.chm13l.Q20.gsv.gz > test_isec.gsv
cat test_isec.gsv | sort -k1,1 -k2,2n | minisv merge - > test_isec_merge1.gsv
cat test_isec.gsv | sort -k1,1 -k2,2n | /hlilab/alvin/miniconda3/bin/k8 ../../minisv_js_latest/minisv.js merge - > test_isec_merge2.gsv
wc -l test_isec_merge1.gsv test_isec_merge2.gsv
diff test_isec_merge1.gsv test_isec_merge2.gsv | wc -l
}


main() {
test1
}

main
