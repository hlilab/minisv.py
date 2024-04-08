## <a name ="started"></a>Getting started

### <a name="docker"></a>Use docker

```sh
sudo apt install build-essential
mamba create -n gaftools python==3.10 poetry cython==3.0.7
conda activate gaftools
git clone https://github.com/qinqian/gaftools
make

tumor_gaf=tumor_gaf
normal_gaf=normal_gaf

# tumor-only or normal-only mode
gaftools sv -b centromere.bed -n tumor tumor_gaf > tumor_sv.bed
sort -k1,1 -k2,2n tumor_sv.bed | gaftools merge - > tumor_mergedsv.bed

# normal-only model
gaftools sv -b centromere.bed -n normal normal_gaf > normal_sv.bed
sort -k1,1 -k2,2n normal_sv.bed | gaftools merge - > normal_mergedsv.bed

# tumor-normal pair mode
gaftools sv -b centromere.bed -n tumor tumor_gaf > tumor_sv.bed
gaftools sv -b centromere.bed -n normal normal_gaf > normal_sv.bed
cat tumor_sv.bed normal_sv.bed | sort -k1,1 -k2,2n - | gaftools merge - > tumor_normal_pair_mergedsv.bed


Or use docker version

docker build -t gaftools .

```

### <a name="ubuntu"></a>Use Ubuntu

```sh
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
