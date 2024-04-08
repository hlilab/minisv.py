import re

path_seg_pattern = re.compile(r"([><])([^><:\s]+):(\d+)-(\d+)")
cigar_pattern = re.compile(r"(\d+)([=XIDM])")
ds_pattern = re.compile(r"([\+\-\*:])([A-Za-z\[\]0-9]+)")
re_tsd = re.compile(r"(\[([A-Za-z]+)\])?([A-Za-z]+)(\[([A-Za-z]+)\])?")
re_info = re.compile(r"([^;\s=]+)=([^;\s=]+)")
sv_region_regex = re.compile(r"(\d+),(\d+)")
