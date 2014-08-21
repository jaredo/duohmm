import sys
minp = 0.1


if len(sys.argv)<2:
    print "\nUsage: python mapavg.py file1.rec file2.rec ... fileN.rec\n"
    print "Possible recombination sites will be printed to stdout.\n"


files = sys.argv[1:]
n = 0
nsnp = 0
recmap = {}
pos = []
for i,f in enumerate(files):

    try:
        fin = open(f)
    except IOError:
        print "Problem reading",f
        sys.exit(1)

    pos = fin.next().split()[2:]

    for row in fin:
        child,parent = row[:1000].split()[:2]
        p = [float(val) for val in row.split()[2:]]
        if i==0:
            recmap[(child,parent)] = p   
            nsnp = len(p)
        else:
            for j in range(nsnp):
                recmap[(child,parent)][j] += p[j]
                
nfile = len(files)

print "CHILD\tPARENT\tSTART\tEND\tPROB_RECOMBINATION"
for k,v in recmap.iteritems():
    for i,p in enumerate(v):
        v[i] /= nfile
    i=0
    while(i<nsnp):
        while i<nsnp  and v[i]<minp: i+=1
        start = i
        while i<nsnp and v[i]>=minp: i+=1
        stop = i
        if i<nsnp:
            print "\t".join(map(str,[k[0],k[1],pos[start],pos[stop],v[start]]))
