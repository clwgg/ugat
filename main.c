#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "parse_bam.h"

int usage(char **argv)
{

  printf("\nUsage: %s [options]\n\n", argv[0]);
  printf("Options:\n");
  printf("\t-b\t-\tBAM file\n\n");
  printf("\t-n\t-\tNumber to subsample (in conflict with -f)\n");
  printf("\t-s\t-\tSeed for random subsample\n");
  printf("\t-f\t-\tFraction to subsample (in conflict with -n)\n\n");
  printf("\t-t\t-\tShow C to T fraction at first 20 bases instead of exp-fit p-value\n");
  printf("\t-c\t-\tJust count alignments in BAM file, no subsampling is done (conflict with -f, -n, -t and -s)\n\n");

  return 1;
}


int main(int argc, char **argv) 
{

  samFile *in = NULL;
  bam_hdr_t *h = NULL;
  bam1_t *b = NULL;
  double f = 0;
  int s = 0;
  int n = 0;
  char *bfile = NULL;
  int c = 0;
  int outflag = 0;
  int elem;

  while ((elem = getopt(argc, argv, "b:n:s:f:ct")) >= 0) {
    switch (elem) {
      case 'b': bfile = optarg; break;
      case 'n': n = atoi(optarg); break;
      case 's': s = atoi(optarg); break;
      case 'f': f = atof(optarg); break;
      case 'c': c = 1; break;
      case 't': outflag = 1; break;
    }
  }

  if (!bfile) {
    return usage(argv);
  }

  if (s) srand(s);
  else srand(time(NULL));

  in = sam_open(bfile, "r");
  if (!in) {
    return usage(argv);
  }

  h = sam_hdr_read(in);
  b = bam_init1();
  int ret = 0;

  if (!n && !s && !f && !c) {
    ret = substitutions(in, h, b, 0, outflag);
  }
  else if (c && !n && !s && !f) {
    long count = countaln(in, h, b);
    fprintf(stdout, "%ld\n", count);
  }
  else if (f && !n && !c ) {
    ret = substitutions(in, h, b, f, outflag);
  }
  else if (n && !f && !c) {
    ret = n_subs(in, h, b, n, outflag);
  }
  else {
    ret = usage(argv);
  }

  bam_destroy1(b);
  bam_hdr_destroy(h);
  sam_close(in);

  return ret;

}

