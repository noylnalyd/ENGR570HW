import sys

vecFileOut = sys.argv[1]
solution = sys.argv[2]
vo = open(vecFileOut,"r")
so = open(solution,"r")

while True:
    x = vo.readline()
    y = so.readline()
    if x=="" or y=="":
        break
    if abs(float(x)-float(y))/max(abs(float(x)),abs(float(y))) > 1e-4:
        print("Vectors are not equal")
        print(x)
        print(y)
        quit()

print("Vectors are equal")

vo.close()
so.close()
