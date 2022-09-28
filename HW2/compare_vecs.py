import sys

vecFileOut = sys.argv[1]
solution = sys.argv[2]
vo = open(vecFileOut,"r")
so = open(solution,"r")

for x in vo:
    y = so.readline()
    if abs(float(x)-float(y)) > 1e-4:
        print("Not a good match.")
        quit()

print("A good match.")

vo.close()
so.close()
