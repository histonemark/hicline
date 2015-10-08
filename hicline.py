#!/usr/bin/env python
# -*- coding:utf-8 -*-
import os
import subprocess
import re
import sys
import tempfile
from collections import defaultdict
from gzopen import gzopen
from itertools import izip
import pdb

def trimm_hic_reads(read1_fastq,read2_fastq):
   '''This function trimms each read line at any uncut restriction enzyme site
   (GATC) and conserves the lefmost part. Then it output in fasta format. '''
   
   # Open 2 files to write
   out1 = re.sub(r'.fastq(\.gz)?', 'read1.fasta',read1_fastq)
   out2 = re.sub(r'.fastq(\.gz)?', 'read2.fasta',read2_fastq)
   
   # Continue if files exist
   if os.path.exists(out1) & os.path.exists(out2): return([out1,out2])

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
   print([out1,out2]) 
   return([out1,out2])

def call_gem_mapper_on_fasta_files(ffasta_list):
   """This function takes the trimmed sequences extracted from the
   HiC PE sequencing files and calls gem to do the mapping with up to 3
   mismatches and using 4 threads to align."""

   INDEX = '/dmel_gem_index/dm3R5_pT2_unmasked.gem'
   processed = []

   for ffasta in ffasta_list:
       outfname = re.sub('\.fasta$', '', ffasta)
       if os.path.exists(outfname +  '.map'): 
          processed.append(outfname + '.map')
          continue

       # TODO: specify version info for `gem-mapper`.
       # System call to `gem-mapper` passing the desired arguments.
       subprocess.call([
         'gem-mapper',
         '-I', INDEX ,
         '-i', ffasta,
         '-o', outfname,
         '-m3',
         '-T4',
         '--unique-mapping',
          ])
       processed.append(outfname + '.map')
   # gem-mapper adds `.map` to the output file.
   if len(processed) == 2:
       to_mapf1 = processed[0]
       to_mapf2 = processed[1]	
       return to_mapf1,to_mapf2

def extract_hic_pairs(mapped_pair):
   #pdb.set_trace()
   '''This function takes the PE maped files and constructs a tab
   separated file with coordinate and pair sequences of each contact'''
   paird = defaultdict(int) 
   mapped1 = mapped_pair[0]
   mapped2 = mapped_pair[1]
  
   pairs_fname = re.sub('\.map$', '', mapped1)   
   if os.path.exists(pairs_fname + '.HiC'): return pairs_fname + '.HiC'   

   with open(mapped1) as f,open(mapped2) as g:
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

   with open(pairs_fname + '.HiC','w') as h:      
      for k in paird:
         npair = paird[k]
         h.write('%s\t%s\n' % (k,npair))

   return pairs_fname + '.HiC'

def make_HiC_matrix(pairs_fname,chrom_name):
   '''This function constructs the Hi-C matrices at a 2Kb
   resolution for the first 10000 bins of each selected chromosome from the paired
   file.'''

   N = 10000
   out = [[0]*N for i in range(N)]   

   with open(pairs_fname) as f, \
        open('chr' + chrom_name + '_mtx.tsv','w') as g:
            
      for line in f:
         # Removing PCR duplicates
         if line.split()[-1] != '1': continue
         try:
            A,B = re.sub(r'\t.*', '', line).split(';')
            chrA,locA = A.split(':')
            chrB,locB = B.split(':')
            if chrA != chrom_name or chrB != chrom_name: continue
            #pdb.set_trace()
            locA = int(locA) / 2000
            locB = int(locB) / 2000
            if locA < N and locB < N:
               out[locA][locB] += 1
               out[locB][locA] += 1
         except Exception:
            sys.stderr.write('Error in ' + line)
      
      for i in range(N):
         g.write('\t'.join([str(a) for a in out[i]]) + '\n')
      



def main(read1_fastq,read2_fastq,*args):
   fnames_fasta = trimm_hic_reads(read1_fastq,read2_fastq)
   fnames_mapped = call_gem_mapper_on_fasta_files(fnames_fasta)
   fname_pairs =  extract_hic_pairs(fnames_mapped) 
   fname_mtx = [make_HiC_matrix(fname_pairs,sel_chr) for sel_chr in args]

if __name__ == '__main__':
   main(sys.argv[1], sys.argv[2], *sys.argv[3:])

