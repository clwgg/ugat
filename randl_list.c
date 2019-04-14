#include <stdlib.h>

#include "randl_list.h"

void n_cp(char a[NROW], char b[NROW])
{
  int i;
  for (i=0; i<NROW; ++i) {
    b[i] = a[i];
  }
}

l_list* n_init(int a, char arr[NROW]) 
{

  l_list *l = malloc( sizeof(l_list) );
  l->size = 0;
  l->root = malloc( sizeof(node) );
  node *root = l->root;

  root->x = a;
  n_cp( arr, root->m);
  root->next = NULL;

  l->leaf = l->root;
  l->size++;

  return l;

}

l_list* sorted_add(l_list *l, int a, char arr[NROW]) 
{

  if (!l) return n_init(a, arr);

  node *new = malloc( sizeof(node) );
  new->x = a;
  n_cp( arr, new->m );

  if (a < l->root->x) {
    new->next = l->root;
    l->root = new;
  }
  else if (a > l->leaf->x) {
    new->next = NULL;
    l->leaf->next = new;
    l->leaf = new;
  }
  else {
    node *last = l->root;
    while (last != l->leaf) {
      
      if (a >= last->x && a <= last->next->x) {
        new->next = last->next;
        last->next = new;
        break;
      }

      last = last->next;
    }
  }

  l->size++;
  return l;

}

l_list* n_shift(l_list *l)
{
  node *ptr = l->root;
  l->root = l->root->next;
  free (ptr);

  l->size--;

  return l;
}

void free_list(l_list *l) 
{

  node *root = l->root;

  while (root != NULL) {
    node *ptr = root;
    root = root->next;
    free (ptr);
  }

  free(l);

}

