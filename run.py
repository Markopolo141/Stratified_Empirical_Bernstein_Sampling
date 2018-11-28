import os
import sys

A = sys.argv[1]
B = sys.argv[2]
C = int(sys.argv[3])
D = int(sys.argv[4])

for i in range(C+1):
	os.system("python {} {} {} {}".format(A,B.format(i),i,D))
