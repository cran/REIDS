with open("HTA-2_0.na35.hg19.probeset.csv","r") as file:
	with open("HTA-2_0_PRobesetID_Name_TC.txt","a") as outfile:
		for line in file:	
			if line[1]!="J" and line[1]!="P":
				continue			
			line=line.split(",")
			PSR=line[0].strip('\"')
			TC=line[6].strip('\"')
			TC=TC.split("///")
			TCs=list(set(TC))
			Tc=TC[0]
			row1=[PSR,Tc]
			outfile.write("{0}\t{1}\n".format(row1[0],row1[1]))
	outfile.close()
file.close()


import pandas as pd
D1=pd.read_csv("HTA-2_0.r3.psr",sep="\t")
D2=pd.read_csv("HTA-2_0_PRobesetID_Name_TC.txt",sep="\t",header=None,names=["probesetName","TCID"])

D=pd.merge(D1,D2,how="inner",on="probesetName")

D=D[["PSID","probesetName","TCID"]]

D.to_csv('PSID_Name_TCID.txt', index=None, sep='\t', mode='a')

