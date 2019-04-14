#ifndef M_PARSE_BAM_H
#define M_PARSE_BAM_H

#include "htslib/htslib/sam.h"

long countaln(samFile *in, bam_hdr_t *h, bam1_t *b);
int substitutions(samFile *in, bam_hdr_t *h, bam1_t *b, double f, int outflag);
int n_subs(samFile *in, bam_hdr_t *h, bam1_t *b, int n, int outflag);
void printCtoT(int freq[20][6]);
int t_nuc(char c);
char t_seq(int i);
char t_cig(int i);
void str_rev(char *s);
char t_comp(char c);
void revcomp(char *s);
char n_enc(char c);
char* get_seq(bam1_t *b);
char* get_cig(bam1_t *b);
char* get_md(bam1_t *b);

#endif
