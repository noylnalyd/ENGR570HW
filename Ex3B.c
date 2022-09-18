#include <stdio.h>
#include <stdlib.h>
// #include <stdint.h>
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

int digitCount(int num){
  int runprod = 10;
  int count=1;
  while(runprod<=num){
    count+=1;
    runprod*=10;
  }
  return count;
}

void Ex3B(int N){
    // Part ii
  int i = -1, j = -1, k = -1, ind = -1;

  // Part iii
  int *A = (int *)malloc(N * N * sizeof(int));
  for (i = 0; i < N * N; i++) {
    ind = 0;
    for (j = 0; j < N; j++) {
      if (j & 1) {
        ind += ((i >> j & 1) << (j / 2)) * N;
      } else {
        ind += (i >> j & 1) << (j / 2);
      }
    }
    A[ind] = i + 1;
  }

  // Part iv
  int numSpacesMax = MAX(digitCount(N*N)+1, 2);
  printf("A=[");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < (numSpacesMax - digitCount(A[i*N+j]) - ((i + j) == 0));
           k++)
        printf(" ");
      printf("%d", A[i * N + j]);
    }
    if (i < N - 1)
      printf("\n  ");
    else
      printf("]\n");
  }
  free(A);
}