import sys

txtformat = sys.argv[1]
cid = sys.argv[2]
with open(txtformat, "r") as fid:
    lines = fid.readlines()

    candidates = [x.strip() for x in lines[0].split(',')]

    print(1)
    print("Contest,{},{},{}".format(cid, len(candidates),
        ','.join(candidates)))

    bid = 1
    for l in range(2,len(lines)):
        bline = lines[l];
        bline = bline.replace(")","") 
        bline = bline.replace("(","")

        parts = bline.split(':')
        votes = int(parts[1].strip())

        cands = parts[0].strip()
        clist = [x.strip() for x in parts[0].split(',') if x.strip() != ""] 

        if len(clist) == 0:
            for i in range(votes):
                print("{},{}".format(cid,bid))
                bid += 1
        else:
            for i in range(votes):
                print("{},{},{}".format(cid,bid,','.join(clist)))
                bid += 1 
