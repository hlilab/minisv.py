def is_vcf(tsv):
    return tsv.endswith("vcf")

def is_gzipped(tsv):
    return tsv.endswith("gz")
