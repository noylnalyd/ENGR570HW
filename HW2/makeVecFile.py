import random
import sys

filePrefix = sys.argv[1]
f = open(filePrefix+"IN.txt","w")

n = int(sys.argv[2])
max_val = 1000.0
a = [0.0]*n
random.seed()
for i in range(0,n):
    a[i] = random.uniform(0.0,max_val);
    f.write((str(a[i])+"\n"))

f.close()



    