## <a name ="started"></a>Getting started

### <a name="docker"></a>Use docker

```sh
git clone https://github.com/qinqian/gaftools
cd gaftools && docker build -t gaftools .
tumor_gaf=tumor_gaf
normal_gaf=normal_gaf

# tumor-only or normal-only mode
docker run -i -v $(pwd):$(pwd) gaftools gaftools getsv -c 4 -m 5 -l 50 --input tumor.gaf -r 2 -p output_prefix --vntr trf.bed --l1 l1.fasta -a 2000 -s grch38graph --ds

# tumor-normal pair mode
docker run -i -v $(pwd):$(pwd) gaftools gaftools getsv -c 4 -m 5 -l 50 --input tumor.gaf --normal normal.gaf -r 2 -p output_prefix --vntr trf.bed --l1 l1.fasta -a 2000 -s grch38graph --ds

```

### <a name="ubuntu"></a>Use Ubuntu

```sh
sudo apt install build-essential
mamba create -n gaftools python==3.10 poetry cython==3.0.7
conda activate gaftools
git clone https://github.com/qinqian/gaftools
cd gaftools && make
gaf=input_gaf
normalgaf=normal_gaf
# tumor-only or normal-only mode
gaftools getindel-cython -m 30 -l 100 --input $gaf -r 3 -p output_prefix

# tumor-normal pair mode
docker run -i -v $(pwd):$(pwd) gaftools gaftools getindel -c 4 -m 30 -l 100 --input $(pwd)/$gaf --normal $normalgaf -r 3 -p $(pwd)/${gaf/.gaf/}_mapq30_mlen100_cnt3
```

## <a name="intro"></a>Introduction

`gaftools` is a pangenome tool that predicts germline or somatic structural variation from the long read whole genome sequencing data for tumor-normal pairs of samples, tumor-only or normal-only sample. The experiment could be either normal-only, tumor-only or tumor-normal paired sequencing.

## Table of Contents

<!--<img aligh="right" width="278" src="doc/example1.png">-->

- [Getting Started](#started)
  - [Install with docker](#docker)
  - [Install with ubuntu](#ubuntu)
- [Introduction](#intro)
- [Data preprocessing](#process)
- [Usage](#usage)


## <a name="process"></a>Data preprocessing

We provide [all scripts](https://github.com/qinqian/pangenome_sv_benchmarking) for processing the bam, cram file to prepare input data for gaftools.

Get the reference graph and linear genome,

```sh
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
wget -c https://zenodo.org/records/6983934/files/chm13-90c.r518.gfa.gz?download=1 -O chm13-90c.r518.gfa.gz
wget -c https://zenodo.org/records/6983934/files/GRCh38-90c.r518.gfa.gz?download=1 -O GRCh38-90c.r518.gfa.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip chm13v2.0.fa.gz chm13-90c.r518.gfa.gz GRCh38-90c.r518.gfa.gz
```

Use chm13 genome as an example to process the data,

### Minigraph

```sh
#graph genome
samtools fasta -@ 4 input.bam | minigraph -c -t 4 chm13-90c.r518.gfa - > sample_id.gaf
#linear genome
samtools fasta -@ 4 input.bam | minigraph -c -t 4 chm13v2.0.fa - > sample_id.gaf
```

### Minimap2

```sh
mode="map-hifi"
samtools fasta -@ 4 input.bam | minimap2 -c -x $mode -t 4 chm13v2.0.fa - > sample_id.paf
```

## <a name="usage"></a>Usage
