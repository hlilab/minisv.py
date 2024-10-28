#!/bin/bash 

test1() {
zcat pair_v2/COLO829T.chm13g.paf.gz | head -50000 | gzip > test.gaf.gz
minisv getsv -b ../data/chm13v2.cen-mask.bed -n COLO829T test.gaf.gz > test.gsv
num=$(cat test.gsv | wc -l)
zcat pair_v2/COLO829T.chm13g.def.gsv.gz | head -n ${num} > test2.gsv
diff test.gsv test2.gsv | wc -l
wc -l test.gsv test2.gsv
}


test2() {
zcat pair_v2/HCC1954T.chm13g.paf.gz | head -50000 | gzip > hcc1954_test.gaf.gz
minisv getsv -b ../data/chm13v2.cen-mask.bed -n HCC1954T hcc1954_test.gaf.gz > hcc1954_test.gsv
num=$(cat hcc1954_test.gsv | wc -l)

zcat pair_v2/HCC1954T.chm13g.def.gsv.gz | head -n ${num} > hcc1954_test2.gsv
diff hcc1954_test.gsv hcc1954_test2.gsv | wc -l
wc -l  hcc1954_test.gsv hcc1954_test2.gsv
}


test2_full() {
minisv getsv -b ../data/chm13v2.cen-mask.bed -n HCC1954T pair_v2/HCC1954T.chm13g.paf.gz > hcc1954_test.gsv
num=$(cat hcc1954_test.gsv | wc -l)

zcat pair_v2/HCC1954T.chm13g.def.gsv.gz > hcc1954_test2.gsv
diff hcc1954_test.gsv hcc1954_test2.gsv | wc -l
wc -l hcc1954_test.gsv hcc1954_test2.gsv
}


main() {
test1
test2
#test2_full
}

main

