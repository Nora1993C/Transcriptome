import sys
import os
f=open(sys.argv[1])
h=open(sys.argv[2],'w')
f.readline()
annot=f.readline().split()
for i in range(6,len(annot)): annot[i]=annot[i].split("/")[-1].strip(".Aligned.sortedByCoord.out.bam")
hline=[annot[0]]+annot[6:]
h.write("\t".join(hline)+"\n")
for line in f:
	line=line.split()
	hline=[line[0]]+line[6:]
	h.write("\t".join(hline)+"\n")

f.close()
h.close()

