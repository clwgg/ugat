#ifndef L_LIST_H
#define L_LIST_H

#define NROW 20

typedef struct node {

  int x;
  char m[NROW];
  struct node *next;

} node;

typedef struct l_list {

  node *root;
  node *leaf;
  int size;

} l_list;

void n_cp(char a[NROW], char b[NROW]);
l_list* n_init(int a, char arr[NROW]);
l_list* sorted_add(l_list *l, int a, char arr[NROW]);
l_list* n_shift(l_list *l);
void free_list(l_list *l);

#endif
