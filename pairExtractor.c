/* The MIT License

   Copyright (c) 2016 Mohammed Alser

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/





#include <zlib.h>  
#include <stdio.h>  
#include "kseq.h"  

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>

#define PACKAGE_VERSION "2.0.0"

// STEP 1: declare the type of file handler and the read() function  
KSEQ_INIT(gzFile, gzread)  

void prepare_genome(const char *fn, const char *fn_new)
{
	FILE *fpnew;
	kseq_t *ks;
	gzFile fp_fa;
	uint64_t tot_len;
	int l, n_ref;
	tot_len = n_ref = 0;

	fpnew = fopen(fn_new, "w");
	fp_fa = gzopen(fn, "r");
	ks = kseq_init(fp_fa);
	if (!fpnew)
		fprintf(stderr, "[pairExtractor] file open error\n");
	
	while ((l = kseq_read(ks)) >= 0) {
		fprintf(fpnew,"%s", ks->seq.s); 
		tot_len += l;
		++n_ref;
	}
	
	if (n_ref==0) 
		fprintf(stderr,"[%s] It seems the reference genome you provided is empty!!\n", __func__);
	else 
		fprintf(stderr,"[%s] %d sequences, total length of the reference sequence: %llu\n", __func__, n_ref, (long long)tot_len);

	kseq_destroy(ks);
	gzclose(fp_fa);
	fclose(fpnew);
}


static int extract_pairs(const char *genome,const char *fpin, const char *fpout, int readlength)
{
	fprintf(stderr,"[%s] Extracting the reference sequence...\n", __func__);
	char mystring [readlength];
	int start,i,j;
	FILE *genomenew, *fpnew, *fp;
	char line[1024];
	char *p;

	genomenew = fopen(genome, "r"); //prepared genome
	fp = fopen(fpin, "r");          //input fastq
	fpnew = fopen(fpout, "w");	 //output fastq
    	if (!genomenew || !fpnew || !fp) {
		fprintf(stderr, "[%s] file open error\n", __func__);
		return 1;
	} 
	/*
	fprintf(fpnew,"Name\t");
	fprintf(fpnew,"Read\t");
	fprintf(fpnew,"SIGAR\t");
	fprintf(fpnew,"Ref\n");
	*/
	for (i = 1; i <= 4000000; i++){
		if ((fgets(line, sizeof(line), fp)!= NULL ) && (i>86)){		
			j=1;

			for (p = strtok(line, "\t"); p != NULL; p = strtok(NULL, "\t")){
				if (j==1)
					fprintf(fpnew,"%s\t",p);
				else if (j==4){
					start = atoi(p);
					fseek ( genomenew , start-1 , SEEK_SET );
					if ( fgets (mystring , readlength+1 , genomenew) != NULL ){
						fprintf(fpnew,"%s\t", mystring);
					}		
				}
				else if (j==6)
					fprintf(fpnew,"%s\t",p);
				else if (j==10)
					fprintf(fpnew,"%s\n",p);
				j=j+1;
			}	
    		}
	}
 
	fclose(fp);
	fclose(fpnew);
	fclose(genomenew);
	return 0;
}


/*static int extract_pairs(const char *genome,const char *fpin, const char *fpout, int readlength)
{
	char mystring [readlength];
    	gzFile fp;  
    	kseq_t *seq;  
	int start, end, l;
    	char *token,*token2,*token3;
	FILE *genomenew, *fpnew;


	
	genomenew = fopen(genome, "r"); //prepared genome
	fp = gzopen(fpin, "r");          //input fastq
	fpnew = fopen(fpout, "w");	 //output fastq
    	if (!genomenew || !fpnew) {
		fprintf(stderr, "[%s] file open error\n", __func__);
		return 1;
	} 
	seq = kseq_init(fp);
fprintf(fpnew,"name\n"); 
	while ((l = kseq_read(seq)) >= 0) {  
        	fprintf(fpnew,"name: %s\n", seq->name.s);  
        	if (seq->comment.l) 
			fprintf(fpnew,"comment: %s\n",seq->comment.s); 
		//split the name 
		token = strtok(seq->name.s, "_");
		token2 = strtok(NULL, "_");
		token3 = strtok(NULL, "_");
		start = atoi(token2);
		end = atoi(token3);
		
		fseek ( genomenew , start+1 , SEEK_SET );
		if ( fgets (mystring , readlength , genomenew) != NULL )
       			puts (mystring);
		fprintf(stderr,"[%s] Extracting the reference sequence...%d\n%d\n%d\n %s\n", __func__,start, end, end-start+1, mystring);
	

		fprintf(fpnew,"%s\n%s\n%s\n", token,token2,token3);
		 
		fprintf(fpnew,"seq: %s\n", seq->seq.s);
		fprintf(fpnew,"ref: %s\n", mystring);  
		if (seq->qual.l) 
			fprintf(fpnew,"qual: %s\n", seq->qual.s);  
	
	}    
	kseq_destroy(seq); // STEP 5: destroy seq  
	gzclose(fp); // STEP 6: close the file handler 
	fclose(fpnew);
	fclose(genomenew);
	return 0;
}

*/
static int extractor_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: Pair Extractor\n");
	fprintf(stderr, "Description: It extracts the reference sequence to which a read maps, given the position in the reference sequence,)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Mohammed Alser <mohammedalser@bilkent.edu.tr>\n\n");
	fprintf(stderr, "Usage: ./pairExtractor [options] <human_g1k_v37_prepared.fasta> <in.fastq> <out.fastq> <read length>\n\n");
	fprintf(stderr, "Options: -p <human_g1k_v37.fasta>	reference genome preparation mode\n");
	fprintf(stderr, "\n");
	return 1;
}


int main(int argc, char *argv[])  
{    	
	int c, prep = 0;
	char *genome=" ";
	

	while ((c = getopt(argc, argv, "p:")) >= 0) {
		switch (c) {
		case 'p': genome = optarg; prep = 1; break;
		}
	}

	//Mode 1: reference genome preparation mode
	if (prep) { 
		if (genome== NULL) 
			return extractor_usage();
		prepare_genome(genome,argv[optind]);
	} 
	//Mode 2: read-reference extraction mode
	else {
		if (argc - optind < 4) 
			return extractor_usage();
		fprintf(stderr, "length %s\n%s\n%s\n%d\n", argv[optind], argv[optind+1], argv[optind+2], atoi(argv[optind+3]));
		extract_pairs(argv[optind], argv[optind+1], argv[optind+2], atoi(argv[optind+3]));
	}
	return 0;  
}  
