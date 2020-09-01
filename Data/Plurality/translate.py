import sys
import io

with open(sys.argv[1]) as f:
    lines = f.readlines()
    print(lines[0],end='')
    print(lines[1],end='')

    cntr = 0
    for i in range(2, len(lines)):
        toks = lines[i].strip().split(',')

        cand = toks[0]
        num = int(toks[1])

        for j in range(num):
            print("1,0_0_{},{}".format(cntr, cand))
            cntr += 1
        
