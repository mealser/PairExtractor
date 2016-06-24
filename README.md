# PairExtractor
It extracts the reference sequence to which a read maps, given the position in the reference sequence

Compilation
===========

gcc -g -O2 -Wall -o pairExtractor pairExtractor.c -lz -lm


Simulation and evaluation
===========
Program: Pair Extractor  
Description: It extracts the reference sequence to which a read maps, given the position in the reference sequence,)  
Version: 1.0.0  
Contact: Mohammed Alser <mohammedalser@bilkent.edu.tr>  

Usage: ./pairExtractor [options] <human_g1k_v37_prepared.fasta> <in.fastq> <out.fastq> <read length>  

Options: -p <human_g1k_v37.fasta>	reference genome preparation mode  


Example
-------------------------
```
./pairExtractor -p human_g1k_v37.fasta  
./pairExtractor human_g1k_v37_prepared.fasta in.fastq out.fastq 100  
```
