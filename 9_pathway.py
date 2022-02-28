#!/usr/bin/env python
# coding=utf-8
def output(outFILE,Info):
    h=open(outFILE,'w')
    Trans={"Control":"-","PAH":"+"}
    h.write("Database\tPathway\tGenes\tPvalue\t#Genes\tScore\tHigh express(+:PAH;-:Control)\n")
    for db in Info.keys():
        for upIn in ["Control","PAH"]:
            for path in Info[db][upIn].keys():
                nGenes=len(Info[db][upIn][path][0])
                p=float(Info[db][upIn][path][1])
                score=-math.log(p)/math.log(2)*nGenes
                hline=[db,path,";".join(Info[db][upIn][path][0]),str(p),str(nGenes),str(score),Trans[upIn]]
                h.write("\t".join(hline)+"\n")
    
    h.close()

def getInfo(file):
    f=open(file)
    for i in range(3): f.readline()
    Info={}
    for line in f:
        line=line.strip("\n").split("\t")
        Info[line[0].strip("\"")]=[[],line[1].strip("\"")]
    f.seek(0)
    for i in range(3): f.readline()
    for line in f:
        line=line.strip("\n").split("\t")
        Info[line[0].strip("\"")][0].append(line[-1].strip("\""))

    f.close()
    return Info

if __name__=='__main__':
    import math,sys
    folder="/project/Transcriptome/Test/8_Enrichment"
    Info={}
    for db in ["KEGG","GO"]:
        #print "*"+db
        Info[db]={}
        for upIn in ["HG","NG"]:
            #print "\t->"+upIn
            Info[db][upIn]=getInfo(folder+"/Control_Up_in_"+upIn+"_"+db+"_genes.txt")
    
    output(folder+"/PathwayEnrichment.txt",Info)

