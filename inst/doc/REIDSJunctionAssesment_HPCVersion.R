# Project: Detection of Splicing Variants
# 
# Author: lucp8409
###############################################################################

# Necessary Pacackages
library(REIDS)

#Necessary Data
#load the AS probe sets
#load DABG Samples
#load Line indexes of annotation files
#load transcript data
#load position data

# Reading the data in form the line indexing file created with the cluster function
args <- commandArgs(TRUE)
file_name <- as.character(args[1])
file_pos <- as.numeric(args[2])
line_length <- as.numeric(args[3])
conn <- file(file_name, 'rb')
current.pos <- seek(conn, where = file_pos, origin = 'start')
data <- readBin(conn, 'raw', n = line_length)
s <- paste(sapply(data, rawToChar), collapse='')

t=strsplit(s,"\",\"")

d=unlist(t)

d[1]=substr(d[1],start=2,stop=nchar(d[1]))
d[5]=substr(d[5],start=1,stop=nchar(d[5])-1)


if(d[1]!="GeneID"){
	geneID=as.character(d[1])
	if(length(strsplit(geneID,"\n")[[1]])==2){
		geneID=strsplit(geneID,"[\n\"]")[[1]][3]
	}
	if(length(strsplit(geneID,"\"")[[1]])==2){
		geneID=strsplit(geneID,"[\"]")[[1]][2]
	}
	exonID=as.character(unlist(strsplit(d[2],",")))
	lengthexons=as.integer(unlist(strsplit(d[3],",")))
	
	npersample=sum(lengthexons)
	
	allsamples=d[4]
	samples=as.numeric(unlist(strsplit(allsamples,",")))
	samplenames=as.character(unlist(strsplit(d[5],",")))
	nsamples=length(samplenames)
	
	splitsamples<-function(x,samples,npersample){
		start=1+npersample*(x-1)
		end=npersample*x
		values=samples[start:end]
		return(values)
	}
	
	samplevalues=lapply(c(1:nsamples),function(i) splitsamples(i,samples,npersample) )
	TempData=rbindlist(list(samplevalues))
	data.table::setnames(TempData,colnames(TempData),samplenames)
	
	geneData=data.frame(geneID=rep(geneID,npersample),exonID=rep(exonID,lengthexons))
	geneData=data.frame(lapply(geneData, as.character), stringsAsFactors=FALSE)
	DataS=cbind(geneData,TempData)

	REIDS_IsoformInfo=REIDSJunctionAssesment_HPCVersion(geneID=geneID,DataS=DataS,ASPSR=ASPSR,Juninfo="User",JAnnotI,JAnnot=NULL,EandTrAnnotI=NULL,EandTrAnnot=NULL,positionData=NULL,transcriptData=NULL,Groups=list(),Low_AllSamples=c(),Low_GSamples=c(),Plot=FALSE,Location=NULL,Name="")
	save(REIDS_IsoformInfo, file=paste("....","/REIDS_IsoformInfo_", geneID, ".RData", sep="")) ## FILL IN  A LOCATION TO SAVE THE DATA
	
}


