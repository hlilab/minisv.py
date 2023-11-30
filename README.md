# gaftools


# Installation

``` bash
mamba create -n pangenome python==3.10 poetry
conda activate pangenome
git clone https://github.com/qinqian/gaftools
poetry install
```

# Usage

``` bash
for gaf in *gaf; do
    gaftools parse -m 30 -l 100 --input $gaf -r 3 -p ${gaf/.gaf/}_python_mgutils_mapq30_mlen100_cnt3 &
done
```


