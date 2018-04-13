import pandas as pd

with open("HTA-2_0.na35.hg19.probeset.csv","r") as file:
	with open("HTA-2_0_ExonAnnotations.txt","a") as outfile:
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
			if "to" in E:
				Estemp=E.split("///")
				Es=list(set(Estemp))
				for e in Es:
				#junction: write lines for each with association next to it
					Ese=e.split("_")
					row1=[Tc,PSR,Strand, Ese[0].strip('\"'),3]
					outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(row1[0].split('.')[0],row1[1].split('.')[0],row1[2],row1[3].split('.')[0],row1[4]))		
					row2=[Tc,PSR,Strand, Ese[2].strip('\"'),5]
					outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(row2[0].split('.')[0],row2[1].split('.')[0],row2[2],row2[3].split('.')[0],row2[4]))	
			else:
				Estemp=E.split("///")
				Es=list(set(Estemp))
				for e in Es:
					row=[Tc,PSR,Strand, e.strip('\"'),'']
					outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(row[0].split('.')[0],row[1].split('.')[0],row[2],row[3].split('.')[0],''))	
	outfile.close()
file.close()	

ExonAnnot=pd.read_table("HTA-2_0_ExonAnnotations.txt",header=None)
ExonAnnot=ExonAnnot.sort_values(by=[0,1])
ExonAnnot.to_csv("HTA-2_0_ExonAnnotations.txt",sep="\t",header=False,index=False)

with open("HTA-2_0_ExonAnnotations.txt","r") as file1:
	with open("LineIndexing_ExonAnnot.txt","w") as outfile1:
		line=file1.readline()
		line=line.split("\t")
		TC=line[0]
		begin=0
		linenr=0
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
	outfile1.close()			
file1.close()

					
with open("HTA-2_0.na35.hg19.probeset.csv","r") as file:					
	with open("HTA-2_0_TranscriptAnnotations.txt","a") as outfile2:
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
					Estemp=E.split("///")
					Es=list(Estemp)
					#for e in Es:
						#junction: write lines for each with association next to it
					Ese=Es[tr].split("_")
					row=[Tc,PSR,Strand, Ese[0].strip('\"'),Trs[tr].strip('\"')]
					outfile2.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(row[0].split('.')[0],row[1].split('.')[0],row[2],row[3].split('.')[0],row[4].split('.')[0]))			
					row=[Tc,PSR,Strand, Ese[2].strip('\"'),Trs[tr].strip('\"')]
					outfile2.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(row[0].split('.')[0],row[1].split('.')[0],row[2],row[3].split('.')[0],row[4].split('.')[0]))	
			
				else:
					Estemp=E.split("///")
					Es=list(Estemp)
					#for e in Es:
					e=Es[tr]
					row=[Tc, PSR, Strand,e.strip('\"'),Trs[tr].strip('\"')]
					outfile2.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(row[0].split('.')[0],row[1].split('.')[0],row[2],row[3].split('.')[0],row[4].split('.')[0]))	
		
	outfile2.close()
file.close()		

TrAnnot=pd.read_table("HTA-2_0_TranscriptAnnotations.txt",header=None)
TrAnnot=TrAnnot.sort_values(by=[0,1])
TrAnnot.to_csv("HTA-2_0_TranscriptAnnotations.txt",sep="\t",header=False,index=False)


with open("HTA-2_0_TranscriptAnnotations.txt","r") as file2:
	with open("LineIndexing_TrAnnot.txt","w") as outfile2:
		line=file2.readline()
		line=line.split("\t")
		TC=line[0]
		begin=0
		linenr=0
		length=1
		for line in file2:
			line=line.split("\t")
			linenr+=1
			TCtemp=line[0]
			if TC==TCtemp:
				length+=1
			else:
				row=[TC.split('.')[0],begin,length]
				outfile2.write("{0}\t{1}\t{2}\n".format(row[0],row[1],row[2]))
				begin=linenr
				length=1
				TC=TCtemp
	outfile2.close()					
file2.close()



















