
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description="Mapping the probe set ID, probe set name and transcript ID")
parser.add_argument('input1', metavar='input1', type=str, nargs=1, help="The name of the probe set annotation file")
parser.add_argument('input2', metavar='input2', type=str, nargs=1, help="The name of the .psr file")
parser.add_argument('output', metavar='output', type=str, nargs=1, help="The name of the output file")
args=parser.parse_args()


def main(argv):
	input1=argv.input1[0]
	input2=argv.input2[0]
    output=argv.output[0]

	with open(input1,"r") as file:
		with open("ProbesetID_Name_TC.txt","a") as outfile:
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

	D1=pd.read_csv(input2,sep="\t")
	D2=pd.read_csv("ProbesetID_Name_TC.txt",sep="\t",header=None,names=["probesetName","TCID"])

	D=pd.merge(D1,D2,how="inner",on="probesetName")

	D=D[["PSID","probesetName","TCID"]]

	D.to_csv(output, index=None, sep='\t', mode='a')
if __name__ == "__main__":
    main(args)	
