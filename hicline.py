#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys
import tempfile
from collections import defaultdict
from gzopen import gzopen
from itertools import izip


def trimm_hic_reads(read1_fastq,read2_fastq):
   '''This function trimms each read line at any uncut restriction enzyme site
   (GATC) and conserves the lefmost part. Then it output in fasta format. '''
   
   # Open 2 files to write
   out1 = re.sub(r'.fastq(\.gz)?', 'read1.fasta',fname1)
   out2 = re.sub(r'.fastq(\.gz)?', 'read2.fasta',fname2)
   
   # We cut in enzyme restriction site GATC and make a fasta file
   with gzopen(read1_fastq) as f, gzopen(read2_fastq) as g, \
        open(out1,'w') as y, open(out2,'w') as z:
      for lineno,(line1,line2) in enumerate(izip(f,g)):
         if lineno % 4 != 1: continue
         seq1 = line1.rstrip().split('GATC')[0]
         seq2 = line2.rstrip().split('GATC')[0]
         if len(seq1) > 16 and len(seq2) > 16:
            y.write('>%d\n' % (lineno / 4))
            y.write(seq1 + '\n')
            z.write('>%d\n' % (lineno / 4))
            z.write(seq2 + '\n')
   return(out1,out2)

def call_gem_mapper_on_fasta_file(fname_fasta):
   """This function takes the trimmed sequences extracted from the
   HiC PE sequencing files and calls gem to do the mapping with up to 3
   mismatches and using 4 threads to align."""

   INDEX = '/dmel_gem_index/dm3R5_pT2_unmasked.gem'

   outfname = re.sub('\.fasta$', '', fname_fasta)

   # Skip if file exists.
   if os.path.exists(outfname + '.map'): return outfname + '.map'

   # TODO: specify version info for `gem-mapper`.
   # System call to `gem-mapper` passing the desired arguments.
   subprocess.call([
       'gem-mapper',
       '-I', INDEX ,
       '-i', fname_fasta,
       '-o', outfname,
       '-m3',
       '-T4',
       '--unique-mapping',
   ])
   # gem-mapper adds `.map` to the output file.
   return outfname + '.map'

def extract_hic_pairs(read1_mapped,read2_mapped):
   '''This function takes the PE maped files and constructs a tab
   separated file with coordinate and pair sequences of each contact'''
   paird = defaultdict(int) 
   pairs_fname =  re.sub('\.map$', '', read1_mapped)
   
   with open(sys.argv[1]) as f,open(sys.argv[2]) as g:
      for (line1,line2) in izip(f,g):
         items1 = line1.split()
         items2 = line2.split()
         if items1[-1] == '-' or items2[-1] == '-':
            paird['-'] += 1
            continue
            try:   
               chr1 = items1[-1].split(':')[0]
               pos1 = items1[-1].split(':')[2]
               chr2 = items2[-1].split(':')[0]
               pos2 = items2[-1].split(':')[2]
               pair = chr1 + ':' + pos1 + ';' + chr2 + ':' + pos2
               paird[pair] += 1
            except IndexError:
               sys.stderr.write('Error in this lines %s\n%s\n' % (line1,line2))
               continue
         
   for k in paird:
      npair = paird[k]
      pairs_fname.write('%s\t%s\n' % (k,npair))

   return pairs_fname

def make_HiC_matrix(pairs_fname,chrom_name):
   '''This function constructs the Hi-C matrices at a 2Kb
   resolution for the first 10000 bins of each chromosome from the paired
   file.'''
   chr_tempf = tempfile.NamedTemporaryFile(delete=False)
   with open(pairs_fname) as f:
      
      for line in f:
         try:
            chr1 = line.split()[0].split(';')[0].split(':')[0]
            chr2 = line.split()[0].split(';')[1].split(':')[0]
         except IndexError:
            sys.stderr.write('Abnormal line, %s\n' % line)
            
         if chr1 == chr2 == chrom_name:
            try:
               coor1 = line.split()[0].split(';')[0].split(':')[1]
               coor2 = line.split()[0].split(';')[1].split(':')[1]
               value = line.split()[1]
               out = (coor1,coor2,value)
               chr_tempf.write('%s\t%s\t%s\n' % out)
            except IndexError:
               chr_tempf.write('Shity line, %s\n' % line)


def main(read1_fastq,read2_fastq,*args):
   fnames_fasta = trimm_hic_reads(read1_fastq,read2_fastq)
   fnames_mapped = [call_gem_mapper_on_fasta_file(ffasta) for ffasta in fnames_fasta]
   fname_pairs  = [extract_hic_pairs(fmaped) for fmaped in fnames_mapped]
   fname_mtx = [make_HiC_matrix(fname_pairs,sel_chr) for sel_chr in args]

if __name__ == '__main__':
   main(sys.argv[1], sys.argv[2], *sys.argv[3:])

