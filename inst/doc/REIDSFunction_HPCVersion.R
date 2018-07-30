
# Necessary Pacackages
library(REIDS)

#Necessary Data
#load the DABG probe set vectors

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


if(d[1]!="geneID"){
	
	geneID=as.character(d[1])
	geneID=as.character(d[1])
	if(length(strsplit(geneID,"\n")[[1]])==2){
		geneID=strsplit(geneID,"[\n\"]")[[1]][3]
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
	setnames(TempData,colnames(TempData),samplenames)
	
	geneData=data.frame(geneID=rep(geneID,npersample),exonID=rep(exonID,lengthexons))
	geneData=data.frame(lapply(geneData, as.character), stringsAsFactors=FALSE)
	geneData=cbind(geneData,TempData)
	
	REIDS_Gene=REIDSFunction_HPCVersion(geneID=geneID,geneData=geneData,ASPSR=c(),nsim=1000,informativeCalls=FALSE,Summarize=c("EqualAll","WeightedAll"),rho=0.5,Low_AllSamples=c())
	save(REIDS_Gene, file=paste("....","/REIDS_Gene_", geneID, ".RData", sep="")) ## FILL IN  A LOCATION TO SAVE THE DATA
	
}