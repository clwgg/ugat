#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

#include "htslib/htslib/khash.h"
#include "htslib/htslib/klist.h"

#include "randl_list.h"
#include "parse_bam.h"
#include "exp_fit.h"

#define nofree(x)
KLIST_INIT(list, char, nofree)

long countaln(samFile *in, bam_hdr_t *h, bam1_t *b)
{

  long i = 0;
  while (sam_read1(in,h,b) >= 0) ++i;
  return i;

}

int substitutions(samFile *in, bam_hdr_t *h, bam1_t *b, double f, int outflag)
{

  int freq[20][6];
  memset(freq, 0, sizeof(freq));

  int r = 0;
  if (f) r = rand();

  while (sam_read1(in,h,b) >= 0) {

    if (b->core.flag & BAM_FUNMAP) continue;

    if (f && f > 0.) {
      int k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ r);
      if ((double)(k&0xffffff) / 0x1000000 >= f) continue;
    }


    // CIGAR
    char *cig = get_cig(b);
    if (cig == NULL) continue;

    // SEQUENCE
    char *bases = get_seq(b);

    // MD
    char *md = get_md(b);
    if ( md == NULL ) continue;

    if (bam_is_rev(b)) {
      str_rev(cig);
      revcomp(bases);
      revcomp(md);
    }

    int i;
    int md_i = 0;
    for (i = 0; i < 20; i++) {
      if (cig[i] == 'S' || cig[i] == 'I') continue;
      else {
        if (md[md_i] == 'N') {
          freq[i][t_nuc(bases[i])]++;
        }
        else if (bases[i] == 'T' && md[md_i] == 'C') {
          freq[i][t_nuc(md[md_i])]++;
          freq[i][5]++;
        }
        else {
          freq[i][t_nuc(md[md_i])]++;
        }
        md_i++;
      }
    }
    free(cig);
    free(bases);
    free(md);

  }

  if (outflag) {
    printCtoT(freq);
  }
  else {
    double pval = run_fit(freq);
    printf("p-value of fit: %.7e\n", pval);
  }

  return 0;

}


int n_subs(samFile *in, bam_hdr_t *h, bam1_t *b, int n, int outflag)
{

  l_list *l = 0;

  while (sam_read1(in,h,b) >= 0) {

    if (b->core.flag & BAM_FUNMAP) continue;

    // CIGAR
    char *cig = get_cig(b);
    if (cig == NULL) continue;

    // SEQUENCE
    char *bases = get_seq(b);

    // MD
    char *md = get_md(b);
    if ( md == NULL ) continue;

    if (bam_is_rev(b)) {
      str_rev(cig);
      revcomp(bases);
      revcomp(md);
    }

    char freq[20];
    memset(freq, 0, sizeof(freq));
    int i;
    int md_i = 0;

    for (i = 0; i < 20; i++) {
      if (cig[i] == 'S' || cig[i] == 'I') continue;
      else {
        if (md[md_i] == 'N') {
          freq[i] = n_enc(bases[i]);
        }
        else if (bases[i] == 'T' && md[md_i] == 'C') {
          freq[i] = n_enc(md[md_i]);
          freq[i] |= 32;
        }
        else {
          freq[i] = n_enc(md[md_i]);
        }
        md_i++;
      }
    }
    free(cig);
    free(bases);
    free(md);

    int r = rand();
    l = sorted_add(l, r, freq);
    if(l->size > n) l = n_shift(l);

  }

  int freq[20][6];
  memset(freq, 0, sizeof(freq));

  node *tmp = l->root;
  while (tmp != 0) {

    int i;
    for (i = 0; i < 20; i++) {
      if      ( tmp->m[i] & 32 ) freq[i][5]++;
      if      ( tmp->m[i] & 1 )  freq[i][0]++;
      else if ( tmp->m[i] & 2 )  freq[i][1]++;
      else if ( tmp->m[i] & 4 )  freq[i][2]++;
      else if ( tmp->m[i] & 8 )  freq[i][3]++;
      else if ( tmp->m[i] & 16 ) freq[i][4]++;
    }

    tmp = tmp->next;
  }

  free_list(l);


  if (outflag) {
    printCtoT(freq);
  }
  else {
    double pval = run_fit(freq);
    printf("p-value of fit: %.7e\n", pval);
  }


  return 0;

}

void printCtoT(int freq[20][6])
{
  int i;
  fprintf(stdout, "CtoT\n");
  for (i = 0; i < 20; i++) {
    if (freq[i][5]) fprintf(stdout, "%f\n", (double)freq[i][5] / freq[i][1]);
    else fprintf(stdout, "0\n");
  }
}

int t_nuc(char c)
{
  int i = -1;
  switch (c) {
    case 'A': i = 0; break;
    case 'C': i = 1; break;
    case 'G': i = 2; break;
    case 'T': i = 3; break;
    case 'N': i = 4; break;
  }
  return i;
}

char t_seq(int i)
{
  char c = 0;
  switch (i) {
    case 1:  c = 'A'; break;
    case 2:  c = 'C'; break;
    case 4:  c = 'G'; break;
    case 8:  c = 'T'; break;
    case 15: c = 'N'; break;
  }
  return c;
}

char t_cig(int i)
{
  char c = 0;
  switch (i) {
    case 0: c = 'M'; break;
    case 1: c = 'I'; break;
    case 2: c = 'D'; break;
    case 3: c = 'N'; break;
    case 4: c = 'S'; break;
    case 5: c = 'H'; break;
    case 6: c = 'P'; break;
    case 7: c = '='; break;
    case 8: c = 'X'; break;
  }
  return c;
}

char t_comp(char c) 
{
  switch (c) {
    case 'A': c = 'T'; break;
    case 'C': c = 'G'; break;
    case 'G': c = 'C'; break;
    case 'T': c = 'A'; break;
    case 'a': c = 't'; break;
    case 'c': c = 'g'; break;
    case 'g': c = 'c'; break;
    case 't': c = 'a'; break;
  }
  return c;
}

char n_enc(char c)
{
  char out = 0;

  switch (c) {
    case 'A': out |= 1;  break;
    case 'C': out |= 2;  break;
    case 'G': out |= 4;  break;
    case 'T': out |= 8;  break;
    case 'N': out |= 16; break;
  }
  return out;
}

void str_rev(char *s)
{
  int len  = strlen(s);
  int last = len - 1;
  int i;

  for (i = 0; i < len/2; i++) {
    char tmp = s[i];
    s[i] = s[last - i];
    s[last - i] = tmp;
  }
}

void revcomp(char *s)
{
  int len  = strlen(s);
  int last = len - 1;
  int i;

  for (i = 0; i < len/2; i++) {
    char tmp = s[i];
    s[i] = t_comp(s[last - i]);
    s[last - i] = t_comp(tmp);
  }
}

char* get_cig(bam1_t *b)
{

  uint32_t *cig = bam_get_cigar(b);
  int i;

  kl_list_t *mylist = kl_init(list);

  for (i = 0; i < b->core.n_cigar; ++i) {
    char op = t_cig( bam_cigar_op(cig[i]) );
    if (op == 'D') continue;
    if (strchr("HNPX=", op)) {
      kl_destroy(list, mylist);
      return NULL;
    }
    int l = bam_cigar_oplen(cig[i]);

    while (l) {
      *kl_pushp(list, mylist) = op; l--;
    }
  }


  char *s = malloc(mylist->size+1);
  memset(s, 0, mylist->size+1);
  i = 0;

  while(kl_shift(list, mylist, &s[i]) == 0) i++; 

  kl_destroy(list, mylist);

  return s;

}

char* get_seq(bam1_t *b)
{

  uint8_t *seq = bam_get_seq(b);
  int i;
  char *bases = malloc(b->core.l_qseq+1);
  memset(bases, 0, b->core.l_qseq+1);

  for (i = 0; i < b->core.l_qseq; ++i) {
    bases[i] = t_seq( bam_seqi(seq, i) );
  }

  return(bases);
}

char* get_md(bam1_t *b)
{

  uint8_t *md = bam_aux_get(b, "MD");
  if ( md == NULL ) return NULL;

  kl_list_t *mylist = kl_init(list);
  char *md_tmp = malloc(strlen((const char *)md)+1);
  memset(md_tmp, 0, strlen((const char *)md)+1);

  do {
    md++;    // start from md[1] since string includes the leading 'Z'
    int len = strlen(md_tmp);

    if (md_tmp[0] == '\0') md_tmp[0] = *md;
    else if ( isalpha(*md) && isalpha(md_tmp[len-1]) ) md_tmp[len] = *md;
    else if ( isalpha(*md) && ispunct(md_tmp[len-1]) ) md_tmp[len] = *md;
    else if ( isdigit(*md) && isdigit(md_tmp[len-1]) ) md_tmp[len] = *md;
    else {
      if ( isdigit(md_tmp[0]) ) {
        int j = atoi(md_tmp);
        while (j) {
          *kl_pushp(list, mylist) = 'N'; j--;
        }
      }
      else if ( isalpha(md_tmp[0]) && md_tmp[0] != '^') {
        int j;
        for(j = 0; j < strlen(md_tmp); j++) {
          *kl_pushp(list, mylist) = md_tmp[j];
        }
      }
      memset(md_tmp, 0, len);
      md_tmp[0] = *md;
    }

  } while (*md != '\0');

  free(md_tmp);

  char *s = malloc(mylist->size+1);
  memset(s, 0, mylist->size+1);
  int i = 0;

  while(kl_shift(list, mylist, &s[i]) == 0) i++;

  kl_destroy(list, mylist);

  return(s);

}

