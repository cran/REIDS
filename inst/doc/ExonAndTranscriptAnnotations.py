import pandas as pd
import argparse


parser = argparse.ArgumentParser(description="Retrieving the manufacturer's exon and transcript information")
parser.add_argument('input', metavar='input', type=str, nargs=1, help="The name of the information file")
parser.add_argument('output', metavar='output', type=str, nargs=1, help="The name of the output file")
args=parser.parse_args()


def main(argv):
    input1=argv.input[0]
    output=argv.output[0]
    with open(input1,"r") as file:
        with open(output,"w") as outfile:
            for line in file:    
                if line[1]!="J" and line[1]!="P":
                    continue            
                line=line.split(",")
                PSR=line[0].strip('\"')
                Strand=line[2].strip('\"')
                TC=line[6].strip('\"')
                TC=TC.split("///")
                TCs=list(set(TC))
                Tc=TC[0]
                E=line[8].strip('\"')
                TR=line[14].strip('\"')
                Trstemp=TR.split("///")
                Trs=list(Trstemp)
                for tr in range(len(Trs)):
                    if "to" in E:
                        #print(tr)
                        Estemp=E.split("///")
                        Es=list(Estemp)
                        #print(Es)
                #for e in Es:
                #junction: write lines for each with association next to it
                        Ese=Es[tr].split("_")
                        row1=[Tc,PSR,Ese[0].strip('\"'),3,Trs[tr].strip('\"'),Strand]
                        outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(row1[0].split('.')[0],row1[1].split('.')[0],row1[2].split('.')[0],row1[3],row1[4].split('.')[0],row1[5]))        
                        row2=[Tc,PSR,Ese[2].strip('\"'),5,Trs[tr].strip('\"'),Strand]
                        outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(row2[0].split('.')[0],row2[1].split('.')[0],row2[2].split('.')[0],row2[3],row2[4].split('.')[0],row2[5]))    
                    else:
                        Estemp=E.split("///")
                        Es=list(Estemp)
                        e=Es[tr]
                        row=[Tc,PSR,e.strip('\"'),"",Trs[tr].strip('\"'),Strand]
                        outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(row[0].split('.')[0],row[1].split('.')[0],row[2].split('.')[0],row[3],row[4].split('.')[0],row[5]))                    
        outfile.close()
    file.close()    
    ExonAnnot=pd.read_table(output,header=None)
    ExonAnnot=ExonAnnot.sort_values(by=[0,1])
    ExonAnnot.to_csv(output,sep="\t",header=False,index=False)
	
    with open(output,"r") as file1:
        with open("LineIndexing_EAnnot_TrAnnot.txt","w") as outfile1:
            line=file1.readline()
            line=line.split("\t")
            TC=line[0]
            begin=1
            linenr=1
            length=1
            for line in file1:
                line=line.split("\t")
                linenr+=1
                TCtemp=line[0]
                if TC==TCtemp:
                    length+=1
                else:
                    row=[TC.split('.')[0],begin,length]
                    outfile1.write("{0}\t{1}\t{2}\n".format(row[0],row[1],row[2]))
                    begin=linenr
                    length=1
                    TC=TCtemp
            row=[TC.split('.')[0],begin,length]        
            outfile1.write("{0}\t{1}\t{2}\n".format(row[0],row[1],row[2]))
        outfile1.close()            
    file1.close()

if __name__ == "__main__":
    main(args)                    
















