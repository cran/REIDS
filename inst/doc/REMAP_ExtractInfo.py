'''
Created on 1-mei-2018

@author: lucp8409
'''

#This function extract info from the REMAP output
#The EAnnot File
#The TrAnnot File
#The JAnnot File
#The split up file to be used in the PivotTransformation function

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Retrieving the REMAP exon and transcript information")
parser.add_argument('input1', metavar='input1', type=str, nargs=1, help="The location of the GeneAllInfo.txt file")
parser.add_argument('input2', metavar='input2', type=str, nargs=1, help="The location of the GeneInfo.csv file")
parser.add_argument('output', metavar='output', type=str, nargs=1, help="The name of the output file")
args=parser.parse_args()

#PSRAnnotations=pd.read_table("Output/GeneAllInfo.txt",sep="\t",header=0)
#PSRAnnotationsOrdered=PSRAnnotations.sort_values(by=['TCID_GeneID', 'Probeset'])
#PSRAnnotationsOrdered.to_csv("Output/GeneAllInfo.txt",sep='\t', index=None)

def main(argv):
    input1=argv.input1[0]
    input2=argv.input2[0]
    output=argv.output[0]
    i=0
    with open(input1,"r") as REMAPResult:
        header=REMAPResult.readline()
        line1=REMAPResult.readline()
        print(line1)
        line=line1.split("\t")
        print(line)
        TC0=line[0].split("_")[0]
        TC_Gene0=line[0]
        ProbeSet=line[1].split("_")[0]
        GenePSRs=list()
        GenePSRs.append(ProbeSet)
        EAnnot=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
        LineIndexE_Start=1
        ELines=1
        #LineIndexSplit_Start=1
        #SplitLines=1
        #print(ProbeSet)
        if ProbeSet[0]=="P":
            Exons=line[4].split("|")
            if len(Exons)>1:
                Es = pd.Series(Exons)
                PSR = pd.Series(np.tile([ProbeSet], len(Exons)) )
                TCG = pd.Series(np.tile([TC_Gene0], len(Exons)) )
                T = pd.Series(np.tile([TC0], len(Exons)) )
                Types=pd.Series(np.tile([""], len(Exons)) )
          
               # EAnnot["TC"]=T
                EAnnot["TC_Gene"]=TCG
                EAnnot["PSR_ID"]=PSR
                EAnnot["EAnnot"]=Es
                EAnnot["Type"]=Types
            
            elif len(Exons)==1:
                EAnnot=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
                EAnnot.loc[0] = [TC_Gene0,ProbeSet,Exons[0],""]
            else:
                EAnnot=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
                EAnnot.loc[0] = [TC_Gene0,ProbeSet,"-",""] 
        
        if ProbeSet[0]=="J":
            Exons=line[4].split("//")
            Exons=[x for x in Exons if "|" in x]
            if len(Exons)>1: 
                Exons=[x.split("|") for x in Exons]
                Exons=[item for sublist in Exons for item in sublist]
                Es = pd.Series(Exons)
                PSR = pd.Series(np.tile([ProbeSet], len(Exons)) )
                TCG = pd.Series(np.tile([TC_Gene0], len(Exons)) )
                T = pd.Series(np.tile([TC0], len(Exons)) )
                Types=pd.Series(np.tile(["3","5"], len(Exons)) )
          
                #EAnnot["TC"]=T
                EAnnot["TC_Gene"]=TCG
                EAnnot["PSR_ID"]=PSR
                EAnnot["EAnnot"]=Es
                EAnnot["Type"]=Types
            
            elif len(Exons)==1:
                EAnnot=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
                EAnnot.loc[0] = [TC_Gene0,ProbeSet,Exons[0],""]
            else:
                EAnnot=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
                EAnnot.loc[0] = [TC_Gene0,ProbeSet,"-",""] 
            
        ELines=len(Exons)     
        for line1 in REMAPResult:
            line=line1.split("\t")
            i+=1
            print(i)
            TC=line[0].split("_")[0]
            TC_Gene=line[0]
            ProbeSet=line[1].split("_")[0]  
            result=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
            #if TC_Gene==TC:
            #print("91")
            if TC_Gene==TC_Gene0:
                GenePSRs.append(ProbeSet)
            else:
                Row=[TC0,TC_Gene0,"|".join(GenePSRs)]
                with open('TC_Gene_SplitFile.txt', 'a') as out:
                     out.write("{0}\t{1}\t{2}\n".format(Row[0],Row[1],Row[2]))
                    #SplitLines+=1
 
                RowLines=[TC_Gene0,LineIndexE_Start,ELines] 
                with open('TC_Gene_REMAP_Indices.txt', 'a') as out2:
                    out2.write("{0}\t{1}\t{2}\n".format(RowLines[0],RowLines[1],RowLines[2]))
                
                LineIndexE_Start=LineIndexE_Start+ELines
                ELines=0
                TC_Gene0=TC_Gene
                GenePSRs=list()
                GenePSRs.append(ProbeSet)
					
            Temp=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
            if ProbeSet[0]=="P":
                Exons=line[4].split("|")
                if len(Exons)>1:
                    Temp=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
                    Es = pd.Series(Exons)
                    PSR = pd.Series(np.tile([ProbeSet], len(Exons)) )
                    TCG = pd.Series(np.tile([TC_Gene], len(Exons)) )
                    T = pd.Series(np.tile([TC], len(Exons)) )
                    Types=pd.Series(np.tile([""], len(Exons)) )
                  
                        #Temp["TC"]=T
                    Temp["TC_Gene"]=TCG
                    Temp["PSR_ID"]=PSR
                    Temp["EAnnot"]=Es
                    Temp["Type"]=Types
                    
                elif len(Exons)==1:                   
                    Temp=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
                    Temp.loc[0] = [TC_Gene,ProbeSet,Exons[0],""]
                else:
                    Temp=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
                    Temp.loc[0] = [TC_Gene,ProbeSet,"-",""]
                    Exons="-" 
                
            elif ProbeSet[0]=="J":
                Exons=line[4].split("//")
                Exons=[x for x in Exons if "|" in x]
                
                if len(Exons)>0: 
                    Temp=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
                    Exons=[x.split("|") for x in Exons]
                    Exons=[item for sublist in Exons for item in sublist]
                    Es = pd.Series(Exons)
                    PSR = pd.Series(np.tile([ProbeSet], len(Exons)) )
                    TCG = pd.Series(np.tile([TC_Gene], len(Exons)) )
                    T = pd.Series(np.tile([TC], len(Exons)) )
                    Types=pd.Series(np.tile(["3","5"], len(Exons)) )
                  
                        #Temp["TC"]=T
                    Temp["TC_Gene"]=TCG
                    Temp["PSR_ID"]=PSR
                    Temp["EAnnot"]=Es
                    Temp["Type"]=Types
                    
                #elif len(Exons)==1:                   
                #    Temp=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
					
                #    Temp.loc[0] = [TC_Gene,ProbeSet,Exons[0],""]
                else:
                    Temp=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
                    Temp.loc[0] = [TC_Gene,ProbeSet,"-",""]
                    Exons="-" 

            EAnnot=EAnnot.append(Temp,ignore_index=True)
            ELines=ELines+len(Exons)

            EAnnot.reset_index()
           
            EAnnot.to_csv('temp.txt',mode="a",sep='\t', index=None, header=None)
 			
            TC0=TC
            TC_Gene0=TC_Gene
            EAnnot=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type"])
			
        Row=[TC0,TC_Gene0,"|".join(GenePSRs)]
        with open('TC_Gene_SplitFile.txt', 'a') as out:
            out.write("{0}\t{1}\t{2}\n".format(Row[0],Row[1],Row[2]))
                    #SplitLines+=1
 
        RowLines=[TC_Gene0,LineIndexE_Start,ELines] 
        with open('TC_Gene_REMAP_Indices.txt', 'a') as out2:
            out2.write("{0}\t{1}\t{2}\n".format(RowLines[0],RowLines[1],RowLines[2]))		
            
        
    LineIndices=pd.read_table("TC_Gene_REMAP_Indices.txt",sep="\t",header=None)     

    with open(input2,"r") as trinput:
        TR=list()
        line1=trinput.readline()
        line=line1.split(",")
        del line[1]
        del line[-1]
        TR.append(line)
        TC0=line[0].split("_")[0]
        TC_Gene0=line[0]   
        print(TC_Gene0)
        StartTr=1
        LenTr=0
        i=0
        for line1 in trinput:
            line=line1.split(",")
            del line[1]
            del line[-1]
            i+=1
            print(i)
            TC=line[0].split("_")[0]
            TC_Gene=line[0]
            print(TC_Gene0)
            if TC==TC0:
                TR.append(line)
            else:
                TRDF=pd.DataFrame(columns=["TC_Gene","Exons","Tr","Strand","ENSG","GeneName"])
                for l in TR:
                    TCG=l[0]
                    Tr=l[2]
                    Strand=np.sign(int(l[6]))
                    ENSG=l[7]
                    Gene=l[8]
                    Exons=l[1].split("|")
                
                    if len(Exons)>1:
                        Temp=pd.DataFrame(columns=["TC_Gene","Exons","Tr","Strand","ENSG","GeneName"])
                        Es = pd.Series(Exons)
                        TCG = pd.Series(np.tile([TCG], len(Exons)) )
                        Tr = pd.Series(np.tile([Tr], len(Exons)) )
                        Strand = pd.Series(np.tile([Strand], len(Exons)) )
                        ENSG = pd.Series(np.tile([ENSG], len(Exons)) )
                        Gene = pd.Series(np.tile([Gene], len(Exons)) )

                        Temp["TC_Gene"]=TCG
                        Temp["Exons"]=Es
                        Temp["Tr"]=Tr
                        Temp["Strand"]=Strand
                        Temp["ENSG"]=ENSG
                        Temp["Gene"]=Gene

                    elif len(Exons)==1:
                        Temp=pd.DataFrame(columns=["TC_Gene","Exons","Tr","Strand","ENSG","GeneName"])
                        Temp.loc[0] = [TCG,Exons[0],Tr,Strand,ENSG,Gene]
                    else:
                        Temp=pd.DataFrame(columns=["TC_Gene","Exons","Tr","Strand","ENSG","GeneName"])
                        Temp.loc[0] = [TCG,"",Tr,Strand,ENSG,Gene]
                
                    TRDF=TRDF.append(Temp,ignore_index=True)    
                TRDF=TRDF[["TC_Gene","Exons","Tr","Strand","ENSG","GeneName"]] 
            
                #try:
                #    Lines=LineIndices.ix[LineIndices.ix[:,0]==TC_Gene0,[1,2]].values.tolist()[0]
                #    OtherIndexing=False
                #except IndexError: 
                Lines=LineIndices.ix[LineIndices.ix[:,0].str.contains(TC0),[1,2]].values.tolist()
                Skip=Lines[0][0]
                NRows=sum([x[1] for x in Lines])
				
                EAnnot=pd.read_table("temp.txt",skiprows=(Skip-1),nrows=NRows,header=None,sep="\t")
                EAnnot.reset_index()
                EAnnotTr=pd.DataFrame()

                for index, row in EAnnot.iterrows():
                    TrEAnnot=TRDF.ix[TRDF.ix[:,1]==str(row[2]).strip("*"),[2,3,4]]
                    TrEAnnot.reset_index()
                    T=0
      
                    row.reset_index()
         
                    for index1, row1 in TrEAnnot.iterrows():
                        row1.reset_index()
                  
                        l1=row.values.tolist()
                        l2=row1.values.tolist()
                        L=l1+l2
                  
                        R=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type","Tr","Strand","ENSG"])
                        R.loc[0] = L
                
                        EAnnotTr=EAnnotTr.append(R,ignore_index=True)
                        T=T+1
                    LenTr=LenTr+len(TrEAnnot.index)
      
            
                EAnnotTr.to_csv(output,mode="a",sep='\t', index=None, header=None)    

                TC0=TC
                TC_Gene0=TC_Gene
                TR=list()
                TR.append(line)
    
        TRDF=pd.DataFrame(columns=["TC_Gene","Exons","Tr","Strand","ENSG","GeneName"])
        for l in TR:
            TCG=l[0]
            Tr=l[2]
            Strand=np.sign(int(l[6]))
            ENSG=l[7]
            Gene=l[8]
            Exons=l[1].split("|")
        
            if len(Exons)>1:
                Temp=pd.DataFrame(columns=["TC_Gene","Exons","Tr","Strand","ENSG","GeneName"])
                Es = pd.Series(Exons)
                TCG = pd.Series(np.tile([TCG], len(Exons)) )
                Tr = pd.Series(np.tile([Tr], len(Exons)) )
                Strand = pd.Series(np.tile([Strand], len(Exons)) )
                ENSG = pd.Series(np.tile([ENSG], len(Exons)) )
                Gene = pd.Series(np.tile([Gene], len(Exons)) )

                Temp["TC_Gene"]=TCG
                Temp["Exons"]=Es
                Temp["Tr"]=Tr
                Temp["Strand"]=Strand
                Temp["ENSG"]=ENSG
                Temp["Gene"]=Gene

            elif len(Exons)==1:
                Temp=pd.DataFrame(columns=["TC_Gene","Exons","Tr","Strand","ENSG","GeneName"])
                Temp.loc[0] = [TCG,Exons[0],Tr,Strand,ENSG,Gene]
            else:
                Temp=pd.DataFrame(columns=["TC_Gene","Exons","Tr","Strand","ENSG","GeneName"])
                Temp.loc[0] = [TCG,"",Tr,Strand,ENSG,Gene]
        
            TRDF=TRDF.append(Temp,ignore_index=True)    
            TRDF=TRDF[["TC_Gene","Exons","Tr","Strand","ENSG","GeneName"]] 
    
        Lines=LineIndices.ix[LineIndices.ix[:,0].str.contains(TC0),[1,2]].values.tolist()
        Skip=Lines[0][0]
        NRows=sum([x[1] for x in Lines])
        EAnnot=pd.read_table("temp.txt",skiprows=(Skip-1),nrows=NRows,header=None,sep="\t")
        EAnnot.reset_index()
        EAnnotTr=pd.DataFrame()
        for index, row in EAnnot.iterrows():
            TrEAnnot=TRDF.ix[TRDF.ix[:,1]==str(row[2]).strip("*"),[2,3,4]]
            TrEAnnot.reset_index()
            T=0
    
            row.reset_index()
       
            for index1, row1 in TrEAnnot.iterrows():
                row1.reset_index()
            
                l1=row.values.tolist()
                l2=row1.values.tolist()
                L=l1+l2
        
                R=pd.DataFrame(columns=["TC_Gene","PSR_ID","EAnnot","Type","Tr","Strand","ENSG"])
                R.loc[0] = L
       
                EAnnotTr=EAnnotTr.append(R,ignore_index=True)
                T=T+1
            LenTr=LenTr+len(TrEAnnot.index)

        EAnnotTr.to_csv(output,mode="a",sep='\t', index=None, header=None)    

        # RowLines=[TC_Gene0,StartTr,LenTr] 
        # with open('REMAP_Indices.txt', 'a') as out2:
            # out2.write("{0}\t{1}\t{2}\n".format(RowLines[0],RowLines[1],RowLines[2]))
    
    with open(output,'r') as out1file:
        with open("REMAP_Indices.txt",'w') as out2file:
            out2file.write('{0}\t{1}\t{2}\n'.format('TC_Gene','Begin','NRows'))
            header=out1file.readline()
            line1=out1file.readline()
            row=line1.split("\t")
            print(row)
            ID=row[0]   
            prevID=ID                
            count=1
            begin=1
            for line in out1file:
                row=line.split()
                ID=row[0]  
                if ID==prevID :
                    count +=1
                    continue
                else:
                    out2file.write('{0}\t{1}\t{2}\n'.format(prevID,begin,count))
                    begin+=count
                    count=1
                    prevID=ID
            out2file.write('{0}\t{1}\t{2}\n'.format(ID,begin,count))   	
        out2file.close()
    out1file.close()       
	
if __name__ == "__main__":
    main(args)					
	
            
            
            
            
  
    
    
    
    
                