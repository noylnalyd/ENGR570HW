// #include <math.h>
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

int main(int N) {
  // Part i
  // N is read as an arg to main!

  // Part ii
  int i = -1, j = -1, k = -1, ind = -1;
  int flag = 0;
  for (i = 2; i <= 16; i *= 2) {
    flag += (i == N);
  }
  if (flag == 0) {
    printf("Invalid N.\n");
    return 1;
  }

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
  int numSpacesMax = MAX(digitCount(N*N)+1, 3);
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
  return 0;
}

