import os
import sys
a = int(sys.argv[1])
b = int(sys.argv[2])+1
for i in range(a,b):
    os.system(f'qdel {i}')
