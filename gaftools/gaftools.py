"""An module for parsing gfa alignment from minigraph
"""

import os
import re

CODE2CIGAR = "MIDNSHP=XB"
CIGAR2CODE = dict((y, x) for x, y in enumerate(CODE2CIGAR))
CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")
CHROM_PATTERN = r"^chr[0-9XY]{1,2}$"


class GafParser(object):
    """Gaf file parser 
    """
    def __init__(self, gaf_path: str):
        self.gaf_path = gaf_path

    def parse_gaf(self):
        with open(self.gaf_path) as fin:
            for line in fin:
                read = line.strip().split()
                parsed_read = AlignedRead(read)
                # large_del = parsed_read.get_deletion_blocks()
                # if large_del:
                #     yield large_del
                large_indel = parsed_read.get_indels()
                if large_indel:
                    yield large_indel


class AlignedRead(object):
    """Parse one read from the gaf file
    """
    def __init__(self, read: list):
        self.read = read

    @property
    def querylen(self):
        return int(self.read[1])

    @property
    def strand(self):
        return self.read[4]

    @property
    def path(self):
        return self.read[5]

    @property
    def mapq(self):
        return int(self.read[11])

    @property
    def cigar(self):
        return self.read[-1]

    @property
    def path_start(self):
        """https://github.com/lh3/gfatools/blob/master/doc/rGFA.md"""
        return int(self.read[7])

    @property
    def path_end(self):
        return int(self.read[8])

    def extract_graph_path_exact_coordinate(self, parts:list, min_del_len:int = 50) -> list:
        """ Calculate exact genome coordinate by converting graph path string, start and end 
        into a bed-like coordinate for a given reference

        We focus on simple deletion example with both reads landing on the reference segments:
        path:'>chr1:460900-479179>HG02109#1#JAHEPG010000279.1:461634-461770>chr1:479192-494904'
        start: 795
        end: 25349
        output: chr1\t460900+795+CIGAR_deletion_start\t460900+795+CIGAR_deletion_end

        We ignore complicated case for the future development.
        path:'<HG01891#1#JAGYVO010000066.1:9329914-9343284>HG02055#1#JAHEPK010000072.1:70420115-70420122<HG01891#1#JAGYVO010000066.1:9307976-9325196'
        start:4523    
        end:23231
        output:?
        """
        segs_major_reference_chrom = [seg.split(':')[0] for seg in re.split(r'[><]', self.path) if seg != '']
        # some CHM13 mapped to GRCH38 such as 'GRCh38#0#chr1', we ignore such a case now
        segs_major_reference_chrom = [chrom for chrom in segs_major_reference_chrom if re.match(CHROM_PATTERN, chrom)]
        if len(segs_major_reference_chrom) == 0:
            return []
        assert len(set(segs_major_reference_chrom)) == 1, 'only one major chrom expected after removing non-chrom segments'

        del_blocks = []
        # graph path start and end are relative coordinate
        ref_pos = self.path_start
        for length, op in parts:
            length = int(length)
            if op == 'D' and length >= min_del_len:
                del_blocks.append((self.path, ref_pos, ref_pos + length, self.strand))
            if op in ["M", "D", "N", "=", "X"]:
                ref_pos += length

        # convert major chromosome position to relative positions
        path_segment_list = [seg for seg in re.split(r'[><]', self.path) if seg != '']
        output_exact_coordiantes = []
        path_seg_cumulative_length = 0 
        for seg_interval in path_segment_list:
            seg_interval = seg_interval.split(':')
            seg_chrom = seg_interval[0]
            seg_start, seg_end = map(int, seg_interval[1].split('-'))
            # 0-based coordinate 
            # https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
            path_seg_cumulative_interval = (path_seg_cumulative_length, path_seg_cumulative_length+seg_end-seg_start)
            path_seg_cumulative_length += seg_end - seg_start
            if re.match(CHROM_PATTERN, seg_chrom):
                for del_block in del_blocks:
                    # Both ends in the path segment
                    if path_seg_cumulative_interval[0] <= del_block[1] and del_block[2] <= path_seg_cumulative_interval[1]:
                        # convert them into exact coordinate
                        # by adding segment starting coordinate and deletion block relative coordinate
                        # TODO: confirm with Heng whether we need to consider partial match
                        output_exact_coordiantes.append((seg_chrom, seg_start+del_block[1], seg_start+del_block[2], self.strand))
        return output_exact_coordiantes

    def extract_linear_reference_coordinate(self, parts: list, min_del_len: int = 50) -> list:
        """ extract linear reference coordiante """
        ref_pos = self.path_start
        del_blocks = []
        for length, op in parts:
            length = int(length)
            if op == 'D' and length >= min_del_len:
                # exclude alternative contigs
                if re.match(CHROM_PATTERN, self.path):
                    del_blocks.append((self.path, ref_pos, ref_pos + length, self.strand))
            if op in ["M", "D", "N", "=", "X"]:
                ref_pos += length
        return del_blocks

    def is_skipped(self):
        # 1. low mapping quality then return zero deletion block
        # 2. no major reference contig then return zero deletion block
        #    this only applies to GAF instead of PAF
        return (self.mapq < 30) or ('chr' not in self.path)

    def get_deletion_blocks(self, min_del_len: int =50) -> list:
        """Identify large deletions from the CIGAR in gaf file
        """
        assert len(self.read) > 12, 'incomplete gaf'
        if self.is_skipped():
            return []

        parts = CIGAR_REGEX.findall(self.cigar)
        # Test the path length 
        assert sum([int(t[0]) for t in parts if t[1] in ["M", "D", "N", "=", "X"]]) == self.path_end - self.path_start
        #print(self.read[6])
        #print(self.path_end - self.path_start)
        # linear genome coordinate from minimap2
        # or linear reference segment coordinate from minigraph
        if ('>' not in self.path) and ('<' not in self.path):
            del_blocks = self.extract_linear_reference_coordinate(parts)
        else:
            del_blocks = self.extract_graph_path_exact_coordinate(parts)
        return del_blocks

    def convert_type(self):
        # convert element types
        for i in range(1, 4):
            self.read[i] = int(self.read[i])
        for i in range(6, 12):
            self.read[i] = int(self.read[i])

    def get_indels(self, min_mapq:int = 5, min_len: int =100, max_cnt: int = 5, min_frac = 0.7, dbg: bool = True) -> list:
        """Get indels from one aligned long read in GAF or PAF file
    
        Update script from https://github.com/lh3/minigraph/blob/master/misc/mgutils-es6.js#L232-L401
        """
        cigar_pattern = re.compile(r"(\d+)([=XIDM])")
        path_seg_pattern = re.compile(r"([><])([^><:\s]+):(\d+)-(\d+)")

        # self.read is t variable in the k8 script above
        if len(self.read) < 12:
            return []

        self.convert_type()
        # low mapping quality
        if self.mapq < min_mapq:
            return []

        # mapped fraction is low
        if self.read[3] - self.read[2] < self.read[1] * min_frac:
            return []

        cg = self.cigar[5:]
        cg_segs = cigar_pattern.findall(self.cigar)

        # record cigar indels
        a = []
        x = self.path_start
        for length, op in cg_segs:
            length = int(length)
            if length >= min_len:
                if op == "I":
                    a.append([x-1, x+1, length])
                elif op == "D":
                    a.append([x, x + length, -length])
            if op in ["M", "D", "=", "X"]:
                x += length

        # No indel cigar or too many 
        if len(a) == 0 or len(a) > max_cnt:
            return []
        if dbg:
            print("X0", a)
            print("X0", self.read)

        # path segments
        seg = []
        if re.match(r'[><]', self.path):
            assert self.strand == '+', "reverse strand on path"
            y = 0
            for m in path_seg_pattern.findall(self.path):
                st, en = int(m[2]), int(m[3])
                # [path, absolute path start, absolute path end, path orientation, relative seg start, relative seg end]
                seg.append([m[1], st, en, 1 if m[0] == '>' else -1,
                            y, y + (en - st)])
        else:
            # https://github.com/lh3/miniasm/blob/master/PAF.md
            # https://github.com/lh3/gfatools/blob/master/doc/rGFA.md
            # [path, 0, path length, 1, 0, path length]
            seg.append([self.path, 0, self.read[6], 1, 0, self.read[6]])

        if dbg:
            print("X1", seg)

        # starts and ends points for the indels on the segments
        # also record start and end segments index 
        # these are used to overlap between path segment and indels 
        sts, ens = [], []
        for i in range(len(a)):
            k = 0
            # segment k relative end <= indel relative start
            while k < len(seg) and seg[k][5] <= a[i][0]:
                k += 1
            if k == len(seg):
                raise Exception("failed to find start position")
            # relative distance between indel start and segment start
            start_l = a[i][0] - seg[k][4]
            # graph path > or linear genome
            if seg[k][3] > 0:
                # [seg index, absolute path start + distance to indel]
                sts.append([k, seg[k][1] + start_l])
            # graph path <, ask Heng 
            else:
                # [seg index, absolute path end + distance to indel]
                sts.append([k, seg[k][2] - start_l])

            k = 0
            # segment k relative end <= indel relative end
            while k < len(seg) and seg[k][5] <= a[i][1]:
                k += 1
            if k == len(seg):
                raise Exception("failed to find end position")
            # relative distance between indel end and segment start
            end_l = a[i][1] - seg[k][4]
            # graph path > or linear genome
            if seg[k][3] > 0:
                # [seg index, absolute path start + distance to indel]
                ens.append([k, seg[k][1] + end_l])
            # graph path <, ask Heng 
            else:
                ens.append([k, seg[k][2] - end_l])

        for i in range(len(a)):
            if dbg:
                print("X2", a[i][0], a[i][1], sts[i][0], ens[i][0])
        
        return []
