#include <stdio.h>
#include <stdlib.h>
int main(int argc, char *argv[]) {
  // Part i
  // N is read as an arg to main!
  int N;
  sscanf(argv[1],"%d",&N);

  // Part ii
  int i = -1;
  int flag = 0;
  for (i = 2; i <= 16; i *= 2) {
    flag += (i == N);
  }
  if (flag == 0) {
    printf("Invalid N.\n");
    exit(1);
  }
  Ex3A(N);
  return 0;
}