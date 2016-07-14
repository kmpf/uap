
# Gero Doose gero@bioinf.uni-leipzig.de
# python script for parsing the output of segemehl into a (tophat-like) cufflinks-compatible output.
# Usage: reads the segemehl /remapper/ realigner-output-file (SAM format) either from stdin or from a file given as the first script argument
# Usage: writes the output to stout
# Example: python s2c.py -s mymap.sam -g hg19.fa > mymap_cufflinks_compatible.sam
# Example:  samtools view -h mymap.bam | python s2c.py -s - -g hg19.fa | samtools view -Sb - | samtools sort - mymap_cufflinks_compatible_sorted.bam

from __future__ import division
import sys
import subprocess
import re
import tempfile
import argparse
import os
from Bio import SeqIO
 
parser = argparse.ArgumentParser(description='python script for parsing the output of segemehl (remapper/realigner) into a cufflinks-compatible (tophat-like) output.     Usage: reads the segemehl-output-file (SAM format) either from stdin or from a file given as the first script argument.    Writes the output to stout.   Example1: python s2c.py -s mymap.sam -g ref_genome.fa > mymap_cufflinks_compatible.sam   Example2:  samtools view -h mymap.bam | python s2c.py -s - -g ref_genome.fa | samtools view -Sb - | samtools sort - mymap_cufflinks_compatible_sorted.bam')
parser.add_argument("-s", "--sam", dest='my_sam', required=True,type=argparse.FileType('r'), help="specifies the path to directory where the segemehl.sam input file is located (when piping in with samtools than just use - as argument )")
parser.add_argument("-g", "--genome", dest='my_genome', help="in case the protocol is not strand specific: specifies the path to the genome")
parser.add_argument("-d", "--maxDist", dest='my_maxDist',type=int, help="specifies the maximal distance of a splice junction. junctions with disctance higher than this value are classified as fusions (default is 200.000nt)")
parser.add_argument("-o", "--outdir", dest='my_outdir', help="specifies the path to ouput directory (temp files are also created in this directory) (default is '.')")
args = parser.parse_args()
if(args.my_maxDist):
    myMaxDist=args.my_maxDist
else:
    myMaxDist=200000
if(args.my_outdir):
    outPath=os.path.abspath(args.my_outdir)
else:
    outPath=os.path.abspath(".")
os.chdir(outPath) 
strandSpecific=False
if(args.my_genome):
    refGenome=args.my_genome
else:
    strandSpecific=True
    

myXSunsure=0
normalReads=0
fusiontranscripts=0
fusionsMore=0
transcripts=0
transcriptsMore=0
wrongOrder=0
multiOrder=0
wrongStrand=0
wrongChr=0
wrongDist=0
wrongOverlapp=0
wrongOrderLarge=0
bestOnlyCounter=0

def getGenomicLength(cigar):
    myM=0
    myD=0
    myX=0
    myEqual=0
    myReg = re.compile(r'\D+',re.I)
    endPosAnzahl=-1
    myIter= re.finditer(myReg, cigar)
    item2= None
    for item2 in myIter:
        if(len(cigar[endPosAnzahl+1:])>0):
            amount=int(cigar[endPosAnzahl+1:item2.start()])
            symbol=cigar[item2.start():item2.start()+1]
            endPosAnzahl=item2.start()
            if(symbol=='M'):
                myM+=amount
            else:
                if(symbol=='D'):
                    myD+=amount
                else:
                    if(symbol=='X'):
                        myX+=amount
                    else:
                        if(symbol=='='):
                            myEqual+=amount
    return myM+myD+myX+myEqual

def updateXQs(mapping):
    XX_index=-1
    XQ_index=-1
    for i in range(11,(len(mapping[0]))):
        if "XX:i" in mapping[0][i]:
            XX_index=i
        if "XQ:i" in mapping[0][i]:
            XQ_index=i
    xq_sorted_mapping = sorted(mapping,key=lambda x: int(x[XX_index][5:]))

    for i in range(0,(len(xq_sorted_mapping))):
        newXQ='XQ:i:'+str(i)
        xq_sorted_mapping[i][XQ_index]=newXQ

    pos_sorted_mapping = sorted(xq_sorted_mapping,key=lambda x: int(x[3]))
    return pos_sorted_mapping

def bestOnly(mappingList):
    global bestOnlyCounter
    NM_index=-1
    NH_index=-1
    for i in range(11,(len(mappingList[0]))):
        if "NM:i" in mappingList[0][i]:
            NM_index=i
        if "NH:i" in mappingList[0][i]:
            NH_index=i
    mymappingDist=[]
    for i in range(0,len(mappingList)):
        mymappingDist.append(0)
        for j in range(0,len(mappingList[i])):
            mymappingDist[i]+=int(mappingList[i][j][NM_index][5:])
    minDist=min(mymappingDist)
    new_mappingList=[]
    for i in range(0, len(mymappingDist)):
        if(mymappingDist[i]==minDist):
            new_mappingList.append(mappingList[i])
        else:
            bestOnlyCounter+=1
    new_NH = 'NH:i:'+str(len(new_mappingList))
    for i in range(0, len(new_mappingList)):
        for j in range(0, len(new_mappingList[i])):
            new_mappingList[i][j][NH_index]= new_NH
    return new_mappingList


def getmapping(fragmentList):
    maxXI=0
    for i in range(0,len(fragmentList)):
        XI_index=-1
        for k in range(11,(len(fragmentList[i]))):
            if "XI:i" in fragmentList[i][k]:
                XI_index=k
        if(XI_index==-1):
            sys.stderr.write('\nno XI set \n')
            sys.stderr.write(str(fragmentList[i])+'\n')
        if(maxXI<int(fragmentList[i][XI_index][5:])):
            maxXI=int(fragmentList[i][XI_index][5:])

    mappingList=[]
    for j in range(0,(maxXI+1)):
        mapping=[]
        for i in range(0,len(fragmentList)):
            XI_index=-1
            for k in range(11,(len(fragmentList[i]))):
                if "XI:i" in fragmentList[i][k]:
                    XI_index=k
            if(int(fragmentList[i][XI_index][5:])==j):
                mapping.append(fragmentList[i])
        if mapping:   
            mappingList.append(updateXQs(mapping))
    if(len(mappingList)>1):
        return bestOnly(mappingList)
    else:
        return mappingList


def isFusion(eingabe):
    global multiOrder, wrongOrder, wrongStrand, wrongChr, wrongDist, wrongOverlapp, fusionsMore, myOrderswitchRange, wrongOrderLarge
    abgehts=False
    bWrongStrand=False
    bWrongChr=False
    bWrongDist=False
    bWrongOverlapp=False
    wrongOrderCounterPerRead=0
    fusion=False
    myline=[]
    mylastLine=[]
    XQ_index_last=-1
    for myline in eingabe:
        XQ_index=-1
        for i in range(11,(len(myline))):
            if "XQ:i" in myline[i]:
                XQ_index=i
        if(abgehts):
            #chromosom switch
            if(myline[2]!=mylastLine[2]):
                bWrongChr=True
            #max dist?
            elif(int(myline[3])-(int(mylastLine[3])+getGenomicLength(mylastLine[5])) >= myMaxDist):
                bWrongDist=True
                #strand switch
            elif((int(myline[1]) >= 16 and int(myline[1]) & 16)  != (int(mylastLine[1]) >= 16 and int(mylastLine[1]) & 16)):
                bWrongStrand=True
            #overlap?
            elif(int(mylastLine[3])+ getGenomicLength(mylastLine[5]) > int(myline[3])):
                bWrongOverlapp=True
            # order switch ?
            elif((int(myline[1]) >= 16 and int(myline[1]) & 16)):
                if((int(myline[XQ_index][5:]) != int(mylastLine[XQ_index_last][5:])-1)):
                    wrongOrderCounterPerRead+=1
            else:
                if((int(myline[XQ_index][5:]) != int(mylastLine[XQ_index_last][5:])+1)):
                    wrongOrderCounterPerRead+=1

        abgehts =True
        mylastLine=myline
        XQ_index_last =XQ_index
    if(bWrongChr):
        wrongChr+=1
        fusion=True
    else:
        if(bWrongStrand):
            wrongStrand+=1
            fusion=True
        else:
            if(bWrongDist):
                wrongDist+=1
                fusion=True
            else:
                if(bWrongOverlapp):
                    wrongOverlapp+=1
                    fusion=True
                else:
                    if(wrongOrderCounterPerRead>1):
                        multiOrder+=1
                        fusion=True
                    else:
                        if(wrongOrderCounterPerRead==1):
                            fusion=True
                            wrongOrder+=1
    if(fusion):
        if(len(eingabe)>2):
            fusionsMore+=1
    return fusion


def processLines(eingabe):
    global myXSunsure
    xaSet=False
    xaIsSet=False
    NH=''
    XI=''
    XA=''
    for i in range(11,(len(eingabe[0]))):
        if "NH:i" in eingabe[0][i]:
            NH=str(eingabe[0][i])
        if "XI:i" in eingabe[0][i]:
            XI=str(eingabe[0][i])
        if "XA:Z" in eingabe[0][i]:
            XA=str(eingabe[0][i])
            
    newLine=[]
    laenge=0
    newCigar=''
    startOld=0
    newNucSeq=''
    qualityIsSet=False
    if(eingabe[0][10]!='*'):
        qualityIsSet=True
    newQualityString=''
    NMi=0
    newMD=''
    md1=''
    md2=''
    for i in range (0, len(eingabe)):
        NM_index=-1
        MD_index=-1
        for j in range(11,(len(eingabe[i]))):
            if "NM:i:" in eingabe[i][j]:
                NM_index= j
            else:
                if "MD:i:" in eingabe[i][j]:
                    MD_index= j

        if(i>0 and i<len(eingabe)):
            newCigar += str(int(eingabe[i][3])-(startOld +(laenge)))
            newCigar+='N'
            if not strandSpecific:
                if(not xaIsSet):
                    subseq=str(record_dict[eingabe[0][2]].seq[(startOld+laenge)-1:startOld+laenge+1])
                    subseq+=str(record_dict[eingabe[0][2]].seq[int(eingabe[i][3])-3:int(eingabe[i][3])-1])
                    subseq=subseq.upper()
                    if(subseq=="GTAG"or subseq=="GCAG" or subseq=="ATAC"):
                        xaIsSet=True
                        xaSet=True
                    if(subseq=="CTAC"or subseq=="CTGC" or subseq=="GTAT"):
                        xaIsSet=True
                        xaSet=False
                
        laenge= getGenomicLength(eingabe[i][5])
        newNucSeq+=eingabe[i][9]
        if(qualityIsSet):
            newQualityString+=eingabe[i][10]
        newCigar+=eingabe[i][5]
        startOld=int(eingabe[i][3])
        if(i==0):
            NMi=int(eingabe[i][NM_index][5:])
            md1=eingabe[i][MD_index][5:]
        else:
            NMi+=int(eingabe[i][NM_index][5:])
            md2=eingabe[i][MD_index][5:]
            myReg = re.compile(r'\D+',re.I)
            myIter= re.finditer(myReg, md1)
            item= None
            for item in myIter:
                pass
            myFirstMatch= re.search(myReg, md2)
            if(item==None and myFirstMatch ==None):
                newMD = str(int(md1)+int(md2))
            else:
                if(item!=None and item.end()!= len(md1) and myFirstMatch ==None):
                    newMD = str(md1[:item.end()])+str(int(md1[item.end():])+int(md2))
                else:
                    if(item==None and myFirstMatch!=None and myFirstMatch.start()>0):
                        newMD = str(int(md1)+int(md2[:myFirstMatch.start()]))+str(md2[myFirstMatch.start():])
                    else:
                        if(item!=None and item.end()!= len(md1) and myFirstMatch!=None and myFirstMatch.start()>0):
                            newMD = str(md1[:item.end()])+ str(int(md1[item.end():])+int(md2[:myFirstMatch.start()]))+str(md2[myFirstMatch.start():])
                        else:
                            newMD = str(md1)+str(md2)
            md1=newMD
    NMi='NM:i:'+str(NMi)
    newMD='MD:Z:'+newMD
    for j in range (0,11):
        if(j==5):
            newLine.append(newCigar)
        else:
            if(j==9):
                newLine.append(newNucSeq)
            else:
                if(j==10 and qualityIsSet):
                    newLine.append(newQualityString)
                else:
                    newLine.append(eingabe[0][j])
    newLine.append(NMi)
    newLine.append(newMD)
    newLine.append(NH)
    newLine.append(XI)
    newLine.append(XA)
    if(strandSpecific):
        if( int(eingabe[0][1])>= 16 and int(eingabe[0][1]) & 16):
            newLine.append('XS:A:-')
        else:
            newLine.append('XS:A:+')
    else:
        if(xaIsSet):
            if(xaSet):
                newLine.append('XS:A:+')
            else:
                newLine.append('XS:A:-')
        else:
            myXSunsure+=1
            newLine.append('XS:A:+')
    return newLine


if strandSpecific:
    sys.stderr.write('since you did not specified a path to a reference genome it is assumed that the protocol is strand specific!  \n')
else:
    sys.stderr.write('since you specified a path to a reference genome it is assumed that the protocol is not strand specific!  \n')
temp = tempfile.NamedTemporaryFile(dir=outPath)
mytempName = temp.name
temp_mate = tempfile.NamedTemporaryFile(dir=outPath)
mytemp_mateName = temp_mate.name
sys.stderr.write('collect all split reads.  \n')
for line in args.my_sam:
    splitRead=False
    columns = line.strip().split('\t')
    if(columns[0][:1]=="@"):
        for fields in columns:
            sys.stdout.write(fields)
            if(not fields==columns[len(columns)-1]):
                sys.stdout.write("\t")
        sys.stdout.write("\n")
    else:
        for i in range(11,(len(columns))):
            if "XQ:i" in columns[i]:
                splitRead=True
        if(splitRead):
            if(int(columns[1]) >= 128 and int(columns[1]) & 128):
                for fields in columns:
                    temp_mate.write(fields)
                    if(not fields==columns[len(columns)-1]):
                        temp_mate.write("\t")
                temp_mate.write("\n")
            else:
                for fields in columns:
                    temp.write(fields)
                    if(not fields==columns[len(columns)-1]):
                        temp.write("\t")
                temp.write("\n")
        else:
            normalReads+=1
            for fields in columns:
                sys.stdout.write(fields)
                if(not fields==columns[len(columns)-1]):
                    sys.stdout.write("\t")
            sys.stdout.write("\n")

temp.flush()
temp.seek(0)
temp2 = tempfile.TemporaryFile(dir=outPath)
p=subprocess.Popen(['sort', '-t','\t','-k', '1,1', '-k', '4,4n',mytempName],stdin=subprocess.PIPE, stdout=subprocess.PIPE)
temp2.write(p.communicate()[0])
temp.close()
temp2.flush()
temp2.seek(0)

listToMerge = []
lastLine = []
fusion= False
losgehts = False
if not strandSpecific:
    sys.stderr.write('read fasta records into memory (for strand unspecific protocols only).  \n')
    handle = open(refGenome, "rU")
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
sys.stderr.write('process the sorted split reads.  \n')
for line in temp2.readlines():
    columns = line.strip().split('\t')
    if(losgehts):
        if(columns[0]==lastLine[0]):
            listToMerge.append(columns)
        else:
            listofmappings = []
            listofmappings=getmapping(listToMerge)
            for i in range(0,len(listofmappings)):
                transcripts+=1
                if(len(listofmappings[i])>2):
                    transcriptsMore+=1
                fusion=isFusion(listofmappings[i])
                if(not fusion):
                    processedLine=processLines(listofmappings[i])
                    for fields in processedLine:
                        sys.stdout.write(fields)
                        if(not fields==processedLine[len(processedLine)-1]):
                            sys.stdout.write("\t")
                    sys.stdout.write("\n")
                else:
                    fusiontranscripts+=1
            listToMerge = []
            listToMerge.append(columns)
            fusion=False
    else:
        listToMerge = []
        listToMerge.append(columns)
        fusion=False
    lastLine = columns
    losgehts=True
if(len(listToMerge)>0):
    listofmappings = []
    listofmappings=getmapping(listToMerge)
    for i in range(0,len(listofmappings)):
        transcripts+=1
        if(len(listofmappings[i])>2):
            transcriptsMore+=1
        fusion=isFusion(listofmappings[i])
        if(not fusion):
            processedLine=processLines(listofmappings[i])
            for fields in processedLine:
                sys.stdout.write(fields)
                if(not fields==processedLine[len(processedLine)-1]):
                    sys.stdout.write("\t")
            sys.stdout.write("\n")
        else:
            fusiontranscripts+=1
temp2.close()

temp_mate.flush()
temp_mate.seek(0)
temp2_mate = tempfile.TemporaryFile(dir=outPath)
p_mate=subprocess.Popen(
    ['sort', '-t','\t','-k', '1,1', '-k', '4,4n', '-T', outPath, mytemp_mateName],
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE)
temp2_mate.write(p_mate.communicate()[0])
temp_mate.close()
temp2_mate.flush()
temp2_mate.seek(0)
listToMerge = []
lastLine = []
fusion= False
losgehts = False
for line in temp2_mate.readlines():
    columns = line.strip().split('\t')
    if(losgehts):
        if(columns[0]==lastLine[0]):
            listToMerge.append(columns)
        else:
            listofmappings = []
            listofmappings=getmapping(listToMerge)
            for i in range(0,len(listofmappings)):
                transcripts+=1
                if(len(listofmappings[i])>2):
                    transcriptsMore+=1
                fusion=isFusion(listofmappings[i])
                if(not fusion):
                    processedLine=processLines(listofmappings[i])
                    for fields in processedLine:
                        sys.stdout.write(fields)
                        if(not fields==processedLine[len(processedLine)-1]):
                            sys.stdout.write("\t")
                    sys.stdout.write("\n")
                else:
                    fusiontranscripts+=1
            listToMerge = []
            listToMerge.append(columns)
            fusion=False
    else:
        listToMerge = []
        listToMerge.append(columns)
        fusion=False
    lastLine = columns
    losgehts=True
if(len(listToMerge)>0):
    listofmappings = []
    listofmappings=getmapping(listToMerge)
    for i in range(0,len(listofmappings)):
        transcripts+=1
        if(len(listofmappings[i])>2):
            transcriptsMore+=1
        fusion=isFusion(listofmappings[i])
        if(not fusion):
            processedLine=processLines(listofmappings[i])
            for fields in processedLine:
                sys.stdout.write(fields)
                if(not fields==processedLine[len(processedLine)-1]):
                    sys.stdout.write("\t")
            sys.stdout.write("\n")
        else:
            fusiontranscripts+=1
temp2_mate.close()

sys.stderr.write('all done.  \n')
sys.stderr.write('removed mappings to retain bestOnly philosophy: '+str(bestOnlyCounter)+'\n')
sys.stderr.write('read mappings (total)            :\t'+str(transcripts+normalReads)+'\n')
sys.stderr.write('read mappings (one fragment)     :\t'+str(normalReads)+'\n')
sys.stderr.write('split reads (total)              :\t'+str(transcripts)+'\tratio to all reads\t'+str(transcripts/(transcripts+normalReads))+'\n')
if(transcripts>0):
    if not strandSpecific:
        sys.stderr.write('split reads not at canonical motives:\t'+str(myXSunsure)+'\tratio to all split reads: '+str(myXSunsure/transcripts)+'\n')
    sys.stderr.write('split reads >2 fragments         :\t'+str(transcriptsMore)+'\tratio to all reads\t'+str(transcriptsMore/(transcripts+normalReads))+'\tratio to all split reads\t'+str(transcriptsMore/(transcripts))+'\n')
    sys.stderr.write('fusion split reads               :\t'+str(fusiontranscripts)+'\tratio to all reads\t'+str(fusiontranscripts/(transcripts+normalReads))+'\tratio to all split reads\t'+str(fusiontranscripts/(transcripts))+'\n')
    if(fusiontranscripts>0):
        sys.stderr.write('fusion split reads > 2 fragments :\t'+str(fusionsMore)+'\tratio to all reads\t'+str(fusionsMore/(transcripts+normalReads))+'\tratio to all split reads\t'+str(fusionsMore/(transcripts))+'\tratio to all fusions\t'+str(fusionsMore/(fusiontranscripts))+'\n')
        sys.stderr.write('  chromosom switch               :\t'+str(wrongChr)+'\tratio to all split reads\t'+str(wrongChr/transcripts)+'\tratio to all fusion reads\t'+str(wrongChr/fusiontranscripts)+'\n')
        sys.stderr.write('  fragment dist > '+str(myMaxDist)+'         :\t'+str(wrongDist)+'\tratio to all split reads\t'+str(wrongDist/transcripts)+'\tratio to all fusion reads\t'+str(wrongDist/fusiontranscripts)+'\n')
        sys.stderr.write('  strand switch                  :\t'+str(wrongStrand)+'\tratio to all split reads\t'+str(wrongStrand/transcripts)+'\tratio to all fusion reads\t'+str(wrongStrand/fusiontranscripts)+'\n')
        sys.stderr.write('  fragment overlapp              :\t'+str(wrongOverlapp)+'\tratio to all split reads\t'+str(wrongOverlapp/transcripts)+'\tratio to all fusion reads\t'+str(wrongOverlapp/fusiontranscripts)+'\n')
        sys.stderr.write('  circular order switch          :\t'+str(wrongOrder)+'\tratio to all split reads\t'+str(wrongOrder/transcripts)+'\tratio to all fusion reads\t'+str(wrongOrder/fusiontranscripts)+'\n')
        sys.stderr.write('  non-circular order switch      :\t'+str(multiOrder)+'\tratio to all split reads\t'+str(multiOrder/transcripts)+'\tratio to all fusion reads\t'+str(multiOrder/fusiontranscripts)+'\n')
sys.stdout.flush()
sys.stdout.close()


# Gero Doose gero@bioinf.uni-leipzig.de
# python script for parsing the output of segemehl into a cufflinks-compatible output.
# Usage: reads the segemehl /remapper/ realigner-output-file (SAM format) either from stdin or from a file given as the first script argument
# Usage: writes the output to stout
# Example: python s2c.py -s mymap.sam -g hg19.fa > mymap_cufflinks_compatible.sam
# Example:  samtools view -h mymap.bam | python s2c.py -s - -g hg19.fa | samtools view -Sb - | samtools sort - mymap_cufflinks_compatible_sorted.bam
