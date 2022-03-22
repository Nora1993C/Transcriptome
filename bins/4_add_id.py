import os
f=open(sys.argv[1])
h=open(sys.argv[2],'w')
line=f.readline().split()
line=["id"]+line
h.write("\t".join(line)+"\n")
for line in f: h.write(line)
f.close()
h.close()
os.system("mv sys.argv[2] sys.argv[1]")

