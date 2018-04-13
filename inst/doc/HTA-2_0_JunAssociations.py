import pandas as pd
import numpy as np 

ExonAnnot=pd.read_table("HTA-2_0_ExonAnnotations.txt",header=None)

D_PSR = ExonAnnot[ExonAnnot.ix[:, 1].str.match("PSR")]
D_JUN = ExonAnnot[ExonAnnot.ix[:, 1].str.match("JUC")]
D_PSR=D_PSR.reset_index(drop=True)
D_JUN=D_JUN.reset_index(drop=True)

with open("HTA-2_0_JunAssociations.txt","a") as outfile:
	for j in D_JUN[1].unique():
		
		PSR3_Final=""
		PSR5_Final=""
		
		#Side 3 annots
		SubJ=D_JUN.loc[(D_JUN.ix[:, 1]==j)&(D_JUN.ix[:, 4]==3), :]
		TC=SubJ[0].iloc[0]
		Strand=SubJ[2].iloc[0]
		
		Es=list(set(SubJ.ix[:,3]))
		PSRs3=D_PSR.loc[D_PSR.ix[:, 3].isin(Es), :]
		

		Count=PSRs3[1].value_counts()
		PSR3=list(Count[Count==len(Es)].index)
		PSRs3temp=list(D_PSR.loc[D_PSR.ix[:, 1].isin(PSR3), 1])

		if Strand=="-" and len(PSRs3temp)>0:
			PSR3_Final=PSRs3temp[0]
		elif Strand=="+" and len(PSRs3temp)>0:
			PSR3_Final=PSRs3temp[-1]
		outfile.write("{0}\t{1}\t{2}\t{3}\n".format(TC,PSR3_Final,j,"3"))	
		
		#Side 5 annots
		SubJ=D_JUN.loc[(D_JUN.ix[:, 1]==j)&(D_JUN.ix[:, 4]==5), :]
		Strand=SubJ[2].iloc[0]
		
		Es=list(set(SubJ.ix[:,3]))
		PSRs5=D_PSR.loc[D_PSR.ix[:, 3].isin(Es), :]
		
		Count=PSRs5[1].value_counts()
		PSR5=list(Count[Count==len(Es)].index)
		PSRs5temp=list(D_PSR.loc[D_PSR.ix[:, 1].isin(PSR5), 1])

		if Strand=="-"and len(PSRs5temp)>0:
			PSR5_Final=PSRs5temp[-1]
		elif Strand=="+"and len(PSRs5temp)>0:
			PSR5_Final=PSRs5temp[0]
		outfile.write("{0}\t{1}\t{2}\t{3}\n".format(TC,PSR5_Final,j,"5"))	
			
		
		#exclusions
		if PSR3_Final!="" and PSR5_Final!="" and Strand=="+":
			I1=D_PSR.index[D_PSR[1]==PSR3_Final].tolist()[-1]
			I2=D_PSR.index[D_PSR[1]==PSR5_Final].tolist()[0]
		elif PSR3_Final!="" and PSR5_Final!="" and Strand=="-":
			I1=D_PSR.index[D_PSR[1]==PSR5_Final].tolist()[-1]
			I2=D_PSR.index[D_PSR[1]==PSR3_Final].tolist()[0]
	
		if I1!=(I2+1):
			R=list(range(I1+1,I2))
			ExclPSR_temp=D_PSR.iloc[R]
			Exrows=ExclPSR_temp[0].count()
			if Exrows!=0:
				ExclPSR=list(set(ExclPSR_temp.ix[:,1]))
				for e in ExclPSR:
					if e!=PSR5_Final and e!=PSR3_Final:
						outfile.write("{0}\t{1}\t{2}\t{3}\n".format(TC,e,j,"exclusion"))
outfile.close()

JunAnnot=pd.read_table("HTA-2_0_JunAssociations.txt",header=None)
print("here")
JunAnnot=JunAnnot.drop_duplicates()
JunAnnot=JunAnnot.sort_values(by=[0,1])
JunAnnot.to_csv("HTA-2_0_JunAssociations.txt",sep="\t",header=False,index=False)


with open("HTA-2_0_JunAssociations.txt","r") as file3:
	with open("LineIndexing_JAnnot.txt","w") as outfile3:
		line=file3.readline()
		line=line.split("\t")
		TC=line[0]
		begin=0
		linenr=0
		length=1
		for line in file3:
			line=line.split("\t")
			linenr+=1
			TCtemp=line[0]
			if TC==TCtemp:
				length+=1
			else:
				row=[TC,begin,length]
				outfile3.write("{0}\t{1}\t{2}\n".format(row[0],row[1],row[2]))
				begin=linenr
				length=1
				TC=TCtemp
		row=[TC,begin,length]
		outfile3.write("{0}\t{1}\t{2}\n".format(row[0],row[1],row[2]))		
	outfile3.close()		
file3.close()



















