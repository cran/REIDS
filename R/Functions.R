## THE REIDS PACKAGE ##


## IMPORTS ##

#' @import aroma.core
#' @import GenomeGraphs
#' @import biomaRt
#' @importFrom MCMCpack riwish
#' @importFrom data.table rbindlist 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom lmtest lrtest
#' @importFrom methods hasArg
#' @importFrom aroma.affymetrix AffymetrixCdfFile AffymetrixCelSet RmaBackgroundCorrection writeCdf QuantileNormalization getCellIndices getData getCdf setCdf ExonRmaPlm getChipEffectSet extractDataFrame FirmaModel getFirmaScores

##### DATA PROCESSING FUNCTIONS #####

## DataProcessing ##
## A wrapper function of some of the functions in the aroma.affymetrix package to make sure the data is processed correctly


#' @title DataProcessing
#' 
#' @description The DataProcessing function processes raw .CEL files to probe intensities values with the help of functions of the aroma.affymetrix package. It returns a data frame and saves it as an .RData file.
#' 
#' @export
#' @param chipType The name of the chip type of the array data.
#' @param tags Tags that is added to the chipType.
#' @param ExonSummarization Logical. Should the data be summarized at the exon level?
#' @param GeneSummarization Logical. Should the data be summarized at the gene level?
#' @param FIRMA Logical. Should the FIRMA model be performed on the data?
#' @param Location The location where the .rda file is to be stored. If NULL, a list containing the requested objects is returned to the user.
#' @param Name A string indicating the prefix for the names of the outputs to be saved at the Location. Defaults to "". 
#' @param verbose Logical. If TRUE, messages are printed during the data processing.
#' @return An .rda file that is saved at the specified location.
#' @details The DataProcessing function is a wrapper of several functions of the aroma.affymetrix package. To obtain the data to perform the REIDS model on the raw .CEL files are 
#' background corrected with the rma background correction and normalization is performed with the quantile normalization. In order for the function to run properly, a chipType and
#' its possible tags need to be specified. It is also important to have the same folder structure as required by the aroma.affymetrix package. This implies the following:
#' a rawData folder with therein a folder with the "Name" parameter. This "Name" folder should contain a folder with the chipType name and herein the .CEL files should be placed.
#' Also a folder annotationData should be present. Herein a folder chipTypes should be make which contains folders for type of chips with the respective names. In the folder of
#' each chiptype the corresponding .cdf file should be saved. If specified, the processed data will be saved at a specific location as a data frame with the first colum the gene IDs and the second column the exon IDs. All
#' other columns contain the sample values. Further the object also contains a vector of the unique gene ID and a vector of the  unique exon IDs. If requested, exon and gene level summarization are performed and saved as data frames at the specified location. Further,the option is provided to perform the FIRMA model on the data as well.
#' @examples
#' \dontrun{
#' DataProcessing(chipType="HTA-2_0",tags="*,r",
#' ExonSummarization=TRUE,GeneSummarization=TRUE,FIRMA=TRUE,
#' Location="HTAData",Name="HTAData",verbose=TRUE)
#' }
DataProcessing<-function(chipType="HuEx-1_0-st-v2",tags="coreR3,A20071112,EP",ExonSummarization=TRUE,GeneSummarization=TRUE,FIRMA=TRUE,Location=NULL,Name="",verbose=TRUE){
	
	if(verbose==TRUE){
		verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
	}
	
	#Setting up the CDF file
	if(is.null(tags)){
		cdf <- AffymetrixCdfFile$byChipType(chipType)
	}
	else{	
		cdf <- AffymetrixCdfFile$byChipType(chipType, tags=tags)
	}
	print(cdf)
	
	#Setting up CELset
	cs <- AffymetrixCelSet$byName(Name, cdf=cdf)
	print(cs)
	
	setCdf(cs, cdf)
	
	
	#Backgroundcorrection
	bc <- RmaBackgroundCorrection(cs,tags=tags)
	csBC <- process(bc,verbose=verbose)
	
	
	#(Quantile) Normalization
	qn <- QuantileNormalization(csBC, typesToUpdate="pm")
	csN <- process(qn, verbose=verbose)
	
	
	if(chipType=="HTA-2_0"){
		flattenCellIndices <- function(cells, ...) {
			# Flatten cell data
			cells <- unlist(cells);
			
			# Do some tricks to clean up the names
			names(cells) <- gsub("[.](groups|indices)", "", names(cells));
			
			# Extract the vector of unit names
			unitNames <- gsub("[.].*", "", names(cells))
			
			# Extract the unit group names
			temp <- gsub("[T].{12}[.]", "",names(cells))
			groupNames<-gsub("[.].*", "", temp)
			
			# Merge data
			data.frame(unitNames=unitNames, groupNames=groupNames, cell=cells)
		} 
	}
	else{
		flattenCellIndices <- function(cells, ...) {
			# Flatten cell data
			cells <- unlist(cells);
			
			# Do some tricks to clean up the names
			names(cells) <- gsub("[.](groups|indices)", "", names(cells));
			
			# Extract the vector of unit names
			unitNames <- gsub("[.].*", "", names(cells))
			
			# Extract the unit group names
			groupNames <- gsub(".*[.]", "", names(cells))
			
			# Merge data
			data.frame(unitNames=unitNames, groupNames=groupNames, cell=cells)
		} 
	}
	
	cells1 <-getCellIndices(cdf, units=1:nbrOfUnits(cdf))
	cells2 <-flattenCellIndices(cells1)

	probeintensities <- getData(csN, indices=cells2$cell, fields=c("intensities"))$intensities
	probeintensities <- log2(probeintensities)
	colnames(probeintensities) <- cs$names
	
	probeintensities=as.data.frame(probeintensities)
	
	GeneID <- as.character(cells2$unitNames)
	ExonID <- as.character(cells2$groupNames)
	#ProbeID<-as.character(cells2$cell)
	
	UniqueGeneID <- unique(GeneID)
	UniqueExonID <- unique(ExonID)
	
	Data=cbind(GeneID,ExonID,probeintensities)
	Data=as.data.frame(Data)
	for(i in c(4:ncol(Data))){
		Data[,i]=round(as.numeric(as.character(Data[,i])),4)
	}
	
	
	
	if(!(is.null(Location))){
		assign(Name,Data,envir=environment())
		eval(parse(text=paste("save(",Name,",UniqueGeneID,UniqueExonID,file=\"",Location,"/", Name, ".RData\")", sep=""))) 			
	}
	
	
	if(ExonSummarization==TRUE){
		
		getCdf(csN)
		plmEx <- ExonRmaPlm(csN, mergeGroups=FALSE)
		fit(plmEx, verbose=verbose)
		
		cesEx <- getChipEffectSet(plmEx)
		exFit <- extractDataFrame(cesEx, units=NULL, addNames=TRUE)
		
		ExonLevelSummarized_rma=exFit[,-c(3,4,5)]
		colnames(ExonLevelSummarized_rma)[c(1,2)]=c("GeneID","ExonID")
		ExonLevelSummarized_rma[,-c(1,2)]=log2(ExonLevelSummarized_rma[,-c(1,2)])
		
		if(!(is.null(Location))){
			assign(paste(Name,"ExonLevelSummarized",sep="_"),ExonLevelSummarized_rma,envir=environment())
			eval(parse(text=paste("save(",Name, "_ExonLevelSummarized, file=\"",Location,"/",Name,"", "_ExonLevelSummarized.RData\")", sep=""))) 
		}
		
		

	}
	
	if(GeneSummarization==TRUE){
		
		plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)
		fit(plmTr, verbose=verbose)
		
		cesTr <- getChipEffectSet(plmTr)
		TrFit <- extractDataFrame(cesTr, units=NULL, addNames=TRUE)
		
		GeneLevelSummarized_rma=TrFit[,-c(2,3,4,5)]
		colnames(GeneLevelSummarized_rma)[1]=c("GeneID")
		GeneLevelSummarized_rma[,-c(1)]=log2(GeneLevelSummarized_rma[,-c(1)])
		
		if(!(is.null(Location))){
			assign(paste(Name,"GeneLevelSummarized",sep="_"),GeneLevelSummarized_rma,envir=environment())
			eval(parse(text=paste("save(",Name, "_GeneLevelSummarized, file=\"",Location,"/",Name,"", "_GeneLevelSummarized.RData\")", sep=""))) 
		}
		
		
	}
	
	if(FIRMA==TRUE){
		plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)
		fit(plmTr,verbose=verbose)

		firma <- FirmaModel(plmTr)
		fit(firma, verbose=verbose)
		
		fs <- getFirmaScores(firma)
		exFit <- extractDataFrame(fs, units=NULL, addNames=TRUE)
		exFit=exFit[,-c(3:5)]
		colnames(exFit)[c(1,2)]=c("GeneID","ExonID")

		FIRMA_Output=exFit
		
		if(!(is.null(Location))){
			assign(paste(Name,"FIRMA_Output",sep="_"),FIRMA_Output,envir=environment())
			eval(parse(text=paste("save(",Name, "_FIRMA_Output, file=\"",Location,"/",Name,"", "_FIRMA_Output.RData\")", sep=""))) 
		}
	}

}


#' "REMAP_SplitProbeSets"
#' 
#' The REMAP_SplitProbeSets function converts the original data frame to a data frame based on the retrieved REMAP annotation.
#' @export
#' @param Data The data frame to be transformed.
#' @param REMAPSplitFile The name of the file with the REMAP information regarding the split of the probe sets if the TC ID is annotated to mutiple genes.
#' @param NotAnnotated Logical. Should the probe sets which are not annotated to a gene still be included? If FALSE, these are excluded. If TRUE, these are included. Default is FALSE.
#' @param Location The location where the file should be saved. If NULL, the object is returned to the user. Otherwise, a file with the specified name is created.
#' @param Name The name of the output file. Defaults to "REMAP".
#' @return A data frame with similar information as the data frame specified in Data with an altered gene ID (first) column in order to match the REMAP annotation.
#' @examples
#' \dontrun{
#' data(TC1500264)
#' 
#' REMAPTest=REMAP_SplitProbeSets(Data=TC1500264,REMAPSplitFile=
#' "TC1500264_Gene_Split.txt")
#' }
REMAP_SplitProbeSets<-function(Data,REMAPSplitFile,NotAnnotated=FALSE,Location=NULL,Name="REMAP"){
	
	SplitFile=utils::read.table(REMAPSplitFile,header=FALSE,stringsAsFactors=FALSE)
	
	Data[,1]=as.character(Data[,1])
	
	Changed=rep(0,nrow(Data))
	for(i in 1:nrow(SplitFile)){
		set=strsplit(SplitFile[i,3],"[|]")[[1]]
		if(NotAnnotated){
			Data[which(as.character(Data[,2])%in%set),1]=SplitFile[i,2]
			Changed[which(as.character(Data[,2])%in%set)]=1	
		}
		else{
			if(grepl("_",SplitFile[i,2])){
				Data[which(as.character(Data[,2])%in%set),1]=SplitFile[i,2]
				Changed[which(as.character(Data[,2])%in%set)]=1	
			}
		}
	}
	
	if(any(Changed==0)){
		Data=Data[-c(which(Changed==0)),]
	}
	
	Data=Data[order(Data[,1]),]

	if(!(is.null(Location))){
		save(Data,file=paste(Location,"/",Name,".RData",sep=""))	
	}
	else{
		return(Data)
	}
}

## Pivot Transformation
## Data format should be an data.frame with an GeneID colum, ExonID column and a column per array ##

#' "PivotTransformation"
#' 
#' The PivotTransformation function converts a data frame with multiple rows per gene into a .csv file with one row per gene. This is the first step in data transformation to apply the REIDS function on a HPC Cluster.
#' @export
#' @param Data The data frame to be transformed.
#' @param GeneID A character vector of the the gene IDs that correspond to the rows of the data frame. Necessary if no GeneID column is present in the data frame
#' @param ExonID A character vector of the the gene IDs that correspond to the rows of the data frame. Necessary if no ExonID column is present in the data frame
#' @param REMAPSplitFile The name of the file with the REMAP information regarding the split of the probe sets if the TC ID is annotated to mutiple genes.
#' @param NotAnnotated Logical. Should the probe sets which are not annotated to a gene still be included? If FALSE, these are excluded. If TRUE, these are included. Default is FALSE.
#' @param Location The location where the file should be saved. If NULL, the object is returned to the user. Otherwise, a file with the specified name is created.
#' @param Name The name of the output file. Defaults to "Pivot".
#' @return A data frame with one row per gene. This row contains the values for each exon per sample and is convenient for processing on a HPC cluster. Futher also a data frame with a column of the gene ID's is returned.
#' @details All information concerning one gene is gathered. The first column of the returned data frame is the gene ID, the second column contains the exon IDs of all exons of that gene. The third colum indicates the number of probes per exon, the fourth contains the values of thos probes per sample and the last column contains the sample names.This way a .csv file is created for processing on a HPC cluster.
#' @examples
#' data(TC12000010)
#' 
#' PivotTest=PivotTransformData(Data=TC12000010,GeneID=NULL,ExonID=NULL,
#' Location=NULL)
#' 
#' \dontrun{
#' data(TC1500264)
#' 
#' PivotTransformData(Data=TC1500264,GeneID=NULL,ExonID=NULL,
#' REMAPSplitFile="TC1500264_Gene_SplitFile.txt",Location=
#' "Output",Name="TC1500264_Pivot")
#' }
PivotTransformData<-function(Data,GeneID=NULL,ExonID=NULL,REMAPSplitFile=NULL,NotAnnotated=FALSE,Location=NULL,Name="Pivot"){
	Data=as.data.frame(Data,stringsAsFactor=FALSE)
	
	if(!is.null(REMAPSplitFile)){
		SplitFile=utils::read.table(REMAPSplitFile,header=FALSE,stringsAsFactors=FALSE)
	}
	else{
		SplitFile=NULL
	}
	
	
	if(is.null(Data$GeneID)&is.null(Data$ExonID)){
		DataTemp=data.frame(GeneID=GeneID,ExonID=ExonID)
		DataTemp=cbind(DataTemp,Data)
		Data=DataTemp
	}else{
		GeneID=as.character(Data$GeneID)
		ExonID=as.character(Data$ExonID)
	}
	
	## From this point we assume that the first column of the data is a Gene ID and the second column is the Exon ID
	
	Transformation<-function(gID,Data,SplitFile,NotAnnotated){
		Subset=Data[which(Data$GeneID==gID),]
		if(!is.null(SplitFile)){
			Split=SplitFile[which(SplitFile[,1]==gID),,drop=FALSE]
			if(nrow(Split)!=1){
				#multiple genes annotated to the same TC ID
				gIDs=SplitFile[,2]
				generow=c()
				for(g in 1:length(gIDs)){
					
					if(!NotAnnotated){
						if(!grepl("_",gIDs[g])){
							next
						}
					}
					
					gIDTC=as.character(gIDs[g])
					
					eIDs=unique(Subset$ExonID)
					Set=strsplit(SplitFile[g,3],"[|]")[[1]]
					eID=eIDs[which(eIDs%in%Set)]
					lengthe=sapply(eID,function(x) return(length(which(Subset$ExonID==x))))
					
					eIDall=paste(as.character(eID),sep="",collapse=",")
					lengthe=paste(as.character(lengthe),sep="",collapse=",")
					
					#get all samples at once: apply on columns of the Subset data and discard the geneID and exonID columns
					samples=apply(Subset[which(Subset[,2]%in%eID),-c(1,2)],2,function(x) paste(round(as.numeric(as.character(x)),4),sep="",collapse=","))
					allsamples=paste(samples,sep="",collapse=",")
					samplenames=paste(colnames(Data)[-c(1,2)],sep="",collapse=",")
					#put everything intro c()
					generow=rbind(generow,c(gIDTC,eIDall,lengthe,allsamples,samplenames))

				}
			}
			else{
				gID=Split[1,2]
				generow=c()
				gID=as.character(gID)
				
				eID=unique(Subset$ExonID)
				lengthe=sapply(eID,function(x) return(length(which(Subset$ExonID==x))))
				
				eID=paste(as.character(eID),sep="",collapse=",")
				lengthe=paste(as.character(lengthe),sep="",collapse=",")
				
				#get all samples at once: apply on columns of the Subset data and discard the geneID and exonID columns
				samples=apply(Subset[,-c(1,2)],2,function(x) paste(as.character(x),sep="",collapse=","))
				allsamples=paste(samples,sep="",collapse=",")
				samplenames=paste(colnames(Data)[-c(1,2)],sep="",collapse=",")
				#put everything intro c()
				generow=c(gID,eID,lengthe,allsamples,samplenames)
			}
		}
		else{
			generow=c()
			gID=as.character(gID)
			
			eID=unique(Subset$ExonID)
			lengthe=sapply(eID,function(x) return(length(which(Subset$ExonID==x))))
			
			eID=paste(as.character(eID),sep="",collapse=",")
			lengthe=paste(as.character(lengthe),sep="",collapse=",")
			
			#get all samples at once: apply on columns of the Subset data and discard the geneID and exonID columns
			samples=apply(Subset[,-c(1,2)],2,function(x) paste(as.character(x),sep="",collapse=","))
			allsamples=paste(samples,sep="",collapse=",")
			samplenames=paste(colnames(Data)[-c(1,2)],sep="",collapse=",")
			#put everything intro c()
			generow=c(gID,eID,lengthe,allsamples,samplenames)
		}	
		return(generow)
	}
	#use rbindlist to get a full data file
	DataPivot=lapply(unique(GeneID),Transformation,Data,SplitFile,NotAnnotated)
	DataBind=do.call("rbind",DataPivot)
	colnames(DataBind)=c("GeneID","ExonID","lengthexons","allsamples","samplenames")
	rownames(DataBind)=unique(DataBind[,1])
	
	DataBind=as.data.frame(DataBind)
	DataBind$GeneID=as.character(DataBind$GeneID)
	DataBind$ExonID=as.character(DataBind$ExonID)
	DataBind$lengthexons=as.character(DataBind$lengthexons)
	DataBind$allsamples=as.character(DataBind$allsamples)
	DataBind$samplenames=as.character(DataBind$samplenames)
	
	GeneID=unique(as.character(DataBind$GeneID))
	GeneTable=data.frame("GeneID"=GeneID)
	
	
	if(!(is.null(Location))){
		utils::write.table(DataBind,file=paste(Location,"/",Name,".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",",quote=TRUE,qmethod="double")
		save(GeneTable,file=paste(Location,"/",Name,"_GeneTable.RData",sep=""))	
	}
	else{
		return(list(DataBind,GeneTable))
	}
	
}

##### REIDS FUNCTIONS ######


# I/NI Calls model

#' "iniREIDS"
#' 
#' The function is part of the larger REIDS function and performs the I/NI calls filtering model on the data for a single gene to select only the informative probesets. The method is performed if Informative=TRUE in the REIDS function before
#' before applying the REIDS model itself. The function is written to fit in the flow of the REIDS model
#' @export
#' @param SubgeneData A subset of the data. Particularly, a subset of the data corresponding to one gene.
#' @param nsim The number of iterations to perform.
#' @return A list with an item per gene. Per gene it is mentioned wich probesets are informative and which are not.
iniREIDS<- function(SubgeneData, nsim=1000) {
	
	# Preparation 
	nprobes <- nrow(SubgeneData)
	nsamples<- ncol(SubgeneData)
	nexon <- G <- length(unique(rownames(SubgeneData)))
	totalObs <- nprobes*nsamples
	geneExp<- as.vector(as.matrix(SubgeneData))
	probes <- as.factor(as.vector(matrix(c(1:nprobes),nrow=nprobes,ncol=nsamples)))
	grp <- as.vector(matrix(rownames(SubgeneData),nrow=nprobes,ncol=nsamples))
	id <-   as.factor(as.vector(t(matrix(c(1:nsamples),nrow=nsamples,ncol=nprobes))))
	
	# Fit of model
	if(length(levels(probes))==1){
		geneRes=rep(0,length(probes))
	}
	else{
		ft <- stats::lm(geneExp~probes-1,x=TRUE) # Why is the model fitted? To get resiudal for estimation of random intercepts?
		
		beta<- as.vector(summary(ft)$coefficient[,1]) #not used?
		fixedDesignMatrix <- as.matrix(as.data.frame(ft$x)) ## create designsamples matrix: not used?
		geneRes <- as.vector(summary(ft)$residuals)  #get residuals ==> for estimation of the ramdom intercepts??
	}
	geneResMat <- matrix(geneRes,nprobes,nsamples) #filled in by column: ok, one column per sample, calculation happens per sample
	
	zid <- matrix(c(1:nsamples),nrow=nsamples,ncol=G)
	uz <- unique(grp)
	zcls2<- sapply(uz,function(x)ifelse(grp==x,1,0))
	
	idv <- as.numeric(id)
	Z<-matrix(0,totalObs ,G*nsamples)
	for (ii in 1:nsamples) Z[idv==ii,(G*(ii-1)+1):(G*ii)]<-zcls2[idv==ii,]
	
	###########
	# Priors  #
	###########
	d0<-g0<-.001		# alpha and beta	
	tau <- 1 #sigma^-2 of measurement error
	
	#################
	# Store Results #
	#################
	nused <- floor(nsim/2)
	nburnin <-  nsim-nused
	outputbik <- matrix(0,nused,nsamples*G)
	outputerror <- NULL
	outputsigma <- matrix(0,nused,G)	# Fixed Effects
	
	
	###############################
	# Fixed Posterior Hyperparameter #
	#    for tau and taub		#
	###############################
	d<-d0+totalObs/2  #alpah+0.5*N
	
	###################
	# GIBBS SAMPLER   #
	###################
	
	nu0 <- G
	c0<-diag(G)			# Scale matrix for Wishart Prior onsamples Sigma.b
	nu<-nu0+nsamples   #gamma+n 
	Taub<-diag(G)			# Ransamplesdom Effects Prec Matrix
	z <- grp
	Taub <- diag(1,G,G) #not needed
	
	for (i in 1:nsim){
		
		set.seed(123*i+i)
		
		# Update b
		cZZ <- diag(colSums(Z))
		#vb <- diag(nsamples)%x%Taub+tau*cZZ #This is Var =D^-1+ sigma^(-2)*ni)
		vb <- solve(diag(nsamples)%x%Taub+tau*cZZ) # This is Var^(-1) =(D^-1+ sigma^(-2)*ni)^-1
		#vbInv=solve(vb)
		mb2<-(tau*crossprod(Z,geneRes ))  #Theta
		mb <- apply(vb,1,function(x) sum(x*mb2) ) #Var^-1*Theta
		vb1 <- diag(vb)  #only need variances here
		b<-sapply(c(1:(nsamples*G)),function(x) stats::rnorm(1,mean=mb[x],sd=sqrt(vb1[x])))  #aren't these the variances^-1 instead of the variances?
		bmat<-matrix(b,ncol=G,byrow=T)  		# Put insamples nsamples x G matrix form for updatinsamplesg taub
		
		
		# Update tau
		zb <- rowSums(sapply(c(1:G),function(x) t(matrix(bmat[,x],nsamples,nprobes))*matrix(zcls2[,x],nprobes,nsamples)))	
		geneReszb <-geneRes-zb	
		g<-g0+crossprod(geneReszb,geneReszb)/2
		tau<-stats::rgamma(1,d,g)
		
		
		# Update Taub 
		
		Sigma.b<-riwish(nu,c0+crossprod(bmat,bmat))
		Taub<-solve(Sigma.b)
		
		if(i > nburnin ){ 
			outputbik[(i-nburnin),] <- b
			outputerror[(i-nburnin)] <- 1/tau
			outputsigma[(i-nburnin),] <- as.vector(diag(Sigma.b))
			
		}
	}
	
	tp1 <- apply(outputsigma,2,function(x) mean(x))
	tp2 <- mean(outputerror)
	icc <- tp1/(tp1+tp2)
	
	Output <- icc
	return(Output)
}


# REIDS model

#' "REIDSmodel_intern"
#' 
#' The function is part of the larger REIDS function and performs the REIDS model on the data for a single gene.
#' @export
#' @param SubgeneData A subset of the data. Particularly, a subset of the data corresponding to one gene.
#' @param nsim The number of iterations to perform. Defaults to 1000.
#' @return A list with 2 items per gene. The first item is the exon scores of the corresponding probesets and the second contains a data frame with the array scores of the exons across the samples. If the iniREIDS model was performed. The items will be added to the previously made list.
REIDSmodel_intern<- function(SubgeneData, nsim=1000) {
	
	nprobes <- nrow(SubgeneData)
	nsamples<- ncol(SubgeneData)
	nexon <- G <- length(unique(rownames(SubgeneData)))
	totalObservations <- nprobes*nsamples
	geneExp<- as.vector(as.matrix(SubgeneData))
	probes <- as.factor(as.vector(matrix(c(1:nprobes),nrow=nprobes,ncol=nsamples)))
	grp <- as.vector(matrix(rownames(SubgeneData),nrow=nprobes,ncol=nsamples))
	id <-   as.factor(as.vector(t(matrix(c(1:nsamples),nrow=nsamples,ncol=nprobes))))
	
	if(length(levels(probes))==1){
		geneRes=rep(0,length(probes))
		Probes<-geneExp
		ID<-rep(0,length(unique(id)))
	}
	else{
		ft <- stats::lm(geneExp~probes+id-1,x=TRUE)  #Difference from inigds: id(samples) is involved in the calculation	
		Probes<- as.vector(summary(ft)$coefficient[1:nprobes,1])
		ID<-as.vector(summary(ft)$coefficient[(nprobes+1):(nprobes+length(unique(id))-1),1])
		ID<-c(0,ID)
		fixedDesignMatrix <- as.matrix(as.data.frame(ft$x)) ## create designsamples matrix
		geneRes <- as.vector(summary(ft)$residuals)
	}
	geneResMat <- matrix(geneRes,nprobes,nsamples)
	zid <- matrix(c(1:nsamples),nrow=nsamples,ncol=G)
	uz <- unique(grp)
	zcls2<- sapply(uz,function(x)ifelse(grp==x,1,0))
	idv <- as.numeric(id)
	Z<-matrix(0,totalObservations ,G*nsamples)
	for (ii in 1:nsamples) Z[idv==ii,(G*(ii-1)+1):(G*ii)]<-zcls2[idv==ii,]
	
	###########
	# Priors  #
	###########
	d0<-g0<-.001			
	tau <- 1
	
	#################
	# Store Results #
	#################
	nused <- floor(nsim/2)
	nburnin <-  nsim-nused
	outputbik <- matrix(0,nused,nsamples*G)
	outputerror <- NULL
	outputsigma <- matrix(0,nused,G)	# Fixed Effects
	###############################
	# Fixed Posterior Hyperparameter #
	#    for tau and taub		#
	###############################
	d<-d0+totalObservations/2
	
	###################
	# GIBBS SAMPLER   #
	###################
	
	nu0 <- G
	c0<-diag(G)			# Scale matrix for Wishart Prior onsamples Sigma.b
	nu<-nu0+nsamples
	Taub<-diag(G)			# Ransamplesdom Effects Prec Matrix
	z <- grp
	Taub <- diag(1,G,G)
	for (i in 1:nsim){
#		if(i==1){
#			ptm<-proc.time()
#		}	
#		set.seed(123*i+i)
		
		# Update b
		cZZ <- diag(colSums(Z))
		#vb <- diag(nsamples)%x%Taub+tau*cZZ #This is Var =D^-1+ sigma^(-2)*ni)
		vb <- solve(diag(nsamples)%x%Taub+tau*cZZ) # This is Var^(-1) =(D^-1+ sigma^(-2)*ni)^-1
		#vbInv=solve(vb)
		mb2<-(tau*crossprod(Z,geneRes))
		mb <- apply(vb,1,function(x) sum(x*mb2))
		vb1 <- diag(vb)
		b<-sapply(c(1:(nsamples*G)),function(x) stats::rnorm(1,mean=mb[x],sd=sqrt(vb1[x])))
		bmat<-matrix(b,ncol=G,byrow=TRUE)  		# Put insamples nsamples x G matrix form for updatinsamplesg taub
		
		
		# Update tau
		zb <- rowSums(sapply(c(1:G),function(x) t(matrix(bmat[,x],nsamples,nprobes))*matrix(zcls2[,x],nprobes,nsamples)))	
		geneReszb <-geneRes-zb	
		g<-g0+crossprod(geneReszb,geneReszb)/2
		tau<-stats::rgamma(1,d,g)
		
		
		# Update Taub 
		
		Sigma.b<-riwish(nu,c0+crossprod(bmat,bmat))
		Taub<-solve(Sigma.b)
		
		if(i > nburnin ){ 
			outputbik[(i-nburnin),] <- b
			outputerror[(i-nburnin)] <- 1/tau
			outputsigma[(i-nburnin),] <- as.vector(diag(Sigma.b))
			
		}
		
#		if(i==1){
#			time<-proc.time()-ptm
#			if(time[1]>2){
#				FailedItems=read.table("FailedGenes.csv",header=FALSE)
#				FailedItems[length(FailedItems)+1]=geneID
#				write.table(FailedItems,file="FailedGenes.csv",row.names=FALSE)
#				stop("Computational time of gds exceeded 2 seconds per iteration")
#			}
#		}
	}
	
	exon_variances=apply(outputsigma,2,function(x) stats::quantile(x, probs = c(0.5)))
	epsilon_error=stats::quantile(outputerror, probs = c(0.5))
	icc <- sapply(exon_variances,function(x) x/(x+epsilon_error))
	icc <- data.frame(exon=uz,type="icc",icc)
	nik <- as.vector(t(sapply(uz,function(x) rep(x,nsamples))))
	exon_effects<- apply(outputbik,2,function(x) stats::quantile(x, probs = c(0.5)))
	exon_effects<- data.frame(exon=nik,type="bik",exon_effects)
	exonscore=data.frame(icc=icc[,3])
	rownames(exonscore)=icc$exon
	
	arrayscore=exon_effects[,3]
	dim(arrayscore)<- c(G,nsamples)
	rownames(arrayscore)<-exon_effects$exon[1:G]
	colnames(arrayscore) <- colnames(SubgeneData)
	
	Output <- list(exonScores=exonscore,arrayScores=arrayscore,errorVar=epsilon_error,exonVar=exon_variances,probeEffects=Probes,sampleEffects=ID)
	return(Output)
}


# REIDS Function - Cluster version
# This function is accompagnied with a .pbs file in the documentation folder. It is advised to run this model an a HPC cluster and not on a regular laptop as it will
# consume time and memory

#' "REIDS_HPCVersion"
#' 
#' The REIDS_ClusterVersion performs the REIDS model and was adapted for use on a HPC cluster. This function should be used with the REIDS_HPCVersion.R file and REIDS_HPCVersion.pbs script in the documentation folder of the package.
#' After running this function on the cluster, the output files should be binded together with the CreateOutput function.
#' @export
#' @param geneID The gene ID
#' @param geneData The data with as rows the probesets and as columns the samples. Note that the first column should contain the gene IDs and the second column the exon IDs
#' @param ASPSR A vector with alternatively spliced probe sets which are taken out of the analysis and summarization. This is useful if Summarize is "WeightedConst" and/or "EqualConst".
#' @param nsim The number of iterations to perform. Defaults to 1000.
#' @param informativeCalls Logical. Should the I/NI calls method be perform before applying the REIDS model?
#' @param Summarize A character vector specifying wich summarization method is to be performed. The choices are "EqualAll", "WeightedAll", "EqualConst" and "WeightedConst". The former two use all probe sets while the latter use only the constituitive probe sets. Summarization on the constistuitive probe sets will only be performed if ASPSR is specified.
#' @param rho The threshold for filtering in the I/NI calls method. Probesets with scores higher than rho are kept.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @return A .RData file will be saved for each gene with the elements returned by the iniREIDS and REIDS functions. The outputs can be bound together by CreateOutput.
REIDSFunction_HPCVersion<- function(geneID,geneData,ASPSR=c(),nsim=1000,informativeCalls=TRUE,Summarize=FALSE,rho=0.5,Low_AllSamples=c()){
	
	if(length(ASPSR)>0){
		AS=which(geneData[,2]%in%ASPSR)
		if(length(AS)>0){
			geneData=geneData[-c(which(geneData[,2]%in%AS)),,drop=FALSE]
		}
	}
	
	DABGs=which(geneData[,2]%in%Low_AllSamples)
	if(length(DABGs)>0){
		geneData=geneData[-c(which(geneData[,2]%in%Low_AllSamples)),,drop=FALSE]
	}
	
	Juncs=which(sapply(geneData[,2],function(x) substr(x,1,3))=="JUC")
	if(length(Juncs)>0){
		geneData=geneData[-c(which(sapply(geneData[,2],function(x) substr(x,1,3))=="JUC")),,drop=FALSE]
	}
	
	exonScore <- arrayScore <- informativeData<- NULL
	
	output=list()
	output[[1]]=list()
	names(output)[1]=geneID
	lcmmData <- geneData[,-c(1,2)]
	lcmmData=as.matrix(lcmmData)
	enames <- geneData$exonID
	names(enames)=NULL
	
	rownames(lcmmData)<- enames
	
	##informative calls
	i=1
	if(informativeCalls){			
		fit <- iniREIDS(SubgeneData=lcmmData, nsim) 
		fit2 <- data.frame(exonNames=unique(enames),Score=fit,informative=fit>rho) # iniREIDS returns one value per exon: filtering on exon level,no replicates
		output[[1]][[i]]=fit2
		names(output[[1]])[i]="Informative"
		iniData <- lcmmData[which(rownames(lcmmData)%in%fit2$exonNames[fit2$informative]),] # Of those that pass filtering step, retrieve the replicates and samples
		lcmmData <-  iniData   
		i=i+1
	}
	
	
	if(!is.null(lcmmData)&length(unique(rownames(lcmmData)))>=3){  
		fit <- REIDSmodel_intern(SubgeneData=lcmmData, nsim) 
		
		exonScore <-fit$exonScores
		arrayScore <- fit$arrayScores
		output[[1]][[i]]=exonScore
		names(output[[1]])[i]="exonScore"
		output[[1]][[i+1]]=arrayScore
		names(output[[1]])[i+1]="arrayScore"
		i=i+2
		
		if(any(c("EqualAll","WeightedAll")%in%Summarize)){
			exonVariances<-fit$exonVar	
			probeEffects<-fit$probeEffects
			sampleEffects<-fit$sampleEffects
			estimatedvalues<-fit$arrayScores
			
			if("WeightedAll"%in%Summarize){## with weights
				TotalExonVar=sum(1/exonVariances)
				Weights=(1/exonVariances)/TotalExonVar
				names(Weights)=unique(rownames(lcmmData))
				
				arrayLevels=apply(estimatedvalues,2,function(j) sum(Weights*j))
				
				WeightedGeneLevelEstimates=mean(probeEffects)+sampleEffects+arrayLevels
				
				output[[1]][[i]]=Weights
				names(output[[1]])[i]="WeightsAllProbesets"
				
				output[[1]][[i+1]]=WeightedGeneLevelEstimates
				names(output[[1]])[i+1]="WeightedAll"
			}
			if("EqualAll"%in%Summarize){
				
				## without weights		
				arrayLevels=apply(estimatedvalues,2,mean)
				GeneLevelEstimates=mean(probeEffects)+sampleEffects+arrayLevels
				
				output[[1]][[i+2]]=GeneLevelEstimates
				names(output[[1]])[i+2]="EqualAll"
			}
			i=i+3
			
			#ExonLevelValues
			
			ProbeLevel=matrix(probeEffects)
			rownames(ProbeLevel)=geneData[,2]
			PSR=unique(geneData[,2])
			ExonLevelEstimates=matrix(0,nrow=length(PSR),ncol=length(sampleEffects))
			for(e in 1:length(PSR)){
				E=mean(ProbeLevel[which(rownames(ProbeLevel)==PSR[e]),])+sampleEffects+estimatedvalues[which(rownames(estimatedvalues)==PSR[e]),]
				ExonLevelEstimates[e,]=E
			}
			rownames(ExonLevelEstimates)=PSR
			
			output[[1]][[i]]=ExonLevelEstimates
			names(output[[1]])[i]="ExonLevel"
			i=i+1
		}
		if(any(c("EqualConst","WeightedConst")%in%Summarize)){
			exonVariances<-fit$exonVar	
			names(exonVariances)=unique(rownames(lcmmData))
			probeEffects<-fit$probeEffects
			sampleEffects<-fit$sampleEffects
			estimatedvalues<-fit$arrayScores
			
			if("WeightedConst"%in%Summarize){## with weights
				TotalExonVar=sum(1/exonVariances)
				Weights=(1/exonVariances)/TotalExonVar
				names(Weights)=unique(rownames(lcmmData))
				
				arrayLevels=apply(estimatedvalues,2,function(j) sum(Weights*j))
				
				WeightedGeneLevelEstimates=mean(probeEffects)+sampleEffects+arrayLevels
				
				output[[1]][[i]]=Weights
				names(output[[1]])[i]="WeightsConstProbesets"
				
				output[[1]][[i+1]]=WeightedGeneLevelEstimates
				names(output[[1]])[i+1]="WeightedConst"
			}
			if("EqualConst"%in%Summarize){
				
				## without weights		
				arrayLevels=apply(estimatedvalues,2,mean)
				GeneLevelEstimates=mean(probeEffects)+sampleEffects+arrayLevels
				
				output[[1]][[i+2]]=GeneLevelEstimates
				names(output[[1]])[i+2]="EqualConst"
			}
			i=i+3
			
			#ExonLevelValues
			
			ProbeLevel=matrix(probeEffects)
			rownames(ProbeLevel)=geneData[,2]
			PSR=unique(geneData[,2])
			ExonLevelEstimates=matrix(0,nrow=length(PSR),ncol=length(sampleEffects))
			for(e in 1:length(PSR)){
				E=mean(ProbeLevel[which(rownames(ProbeLevel)==PSR[e]),])+sampleEffects+estimatedvalues[which(rownames(estimatedvalues)==PSR[e]),]
				ExonLevelEstimates[e,]=E
			}
			rownames(ExonLevelEstimates)=PSR
			
			output[[1]][[i]]=ExonLevelEstimates
			if(length(AS)==0){
				names(output[[1]])[i]="ExonLevel"
				i=i+1
			}
			else{
				names(output[[1]])[i]="ExonLevel_Const"
				i=i+1
			}
		}
	}

	return(output)
}


# CreateOutput

#' "CreateOutput"
#' 
#' The CreateOutput functions writes the .RData files returned by the REIDS_HPCVersion to .txt files: "Name_INICalls.txt", Name_ExonScores.txt", "Name_ArrayScores.txt" and the summarized values distributed across "Name_WeightedAll.txt", "Name_EqualAll.txt", "Name_WeightedConst.txt" and "Name_EqualConst.txt". 
#' The function can also be used for the returned files of REIDSIsoformAssesment_HPCVersion for which it will create the files: "Name_ASInfo.txt" "Name_Compositions.txt","Name_GroupTranscripts.txt" and "Name_NovelConnections".
#' The function is advised to be used with the CreateOutput.R and CreateOutput.pbs file in the documentation folder.
#' @export
#' @param ID A data frame with a "geneID" column.
#' @param Groups A list with elements specifying the columns of the data in each group.
#' @param Name A name for the returned list.
#' @param Location The location where the file should be saved.
#' @return .txt files with the information of the REIDSFunction and REIDSJunctionAssesment.
CreateOutput<-function(ID,Groups,Location="",Name){
	ID=as.vector(as.matrix(ID))
	REIDSOutputFiles=list.files(path=Location,pattern="^REIDS_Gene_")
	if(length(REIDSOutputFiles)>0){
		Output=list()
		for(i in as.character(ID)){
			Data=try(get(load(paste(Location,"/REIDS_Gene_",as.character(i),".RData",sep=""))),silent=TRUE)
			if(class(Data)!="try-error"){
				Output[length(Output)+1]=Data
				names(Output)[length(Output)]=i
			}
		}
		
		assign(paste(Name,"_REIDS_Output",sep=""),Output,envir=environment())
		eval(parse(text=paste("save(",Name, "_REIDS_Output, file=\"",Location,"/",Name, "_REIDS_Output.RData\")", sep=""))) 
		
		for(l in 1:length(Output)){
			utils::write.table(t(c("TC_ID","PSR_ID","Informative")), file = paste(Location,"/",Name,"_INICalls.txtt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			utils::write.table(t(c("TC_ID","PSR_ID","ExonScore")), file = paste(Location,"/",Name,"_ExonScores.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			utils::write.table(t(c("TC_ID","PSR_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_ArrayScores.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			utils::write.table(t(c("TC_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_WeightedAll.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			utils::write.table(t(c("TC_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_EqualAll.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			utils::write.table(t(c("TC_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_WeightedConst.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			utils::write.table(t(c("TC_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_EqualConst.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			utils::write.table(t(c("TC_ID","PSR_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_ExonLevel.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			
			#if gene and exon level summ values are available: collect these.
			if("Informative"%in%names(Output[[l]])){
				
				utils::write.table(t(Output[[l]]$"Informative"), paste(Location,"/",Name,"_INICalls.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	

			}
			if("exonScore"%in%names(Output[[l]])){
				utils::write.table(t(Output[[l]]$"exonScore"), paste(Location,"/",Name,"_ExonScores.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				
			}
			if("arrayScore"%in%names(Output[[l]])){
				utils::write.table(t(Output[[l]]$"arrayScore"), paste(Location,"/",Name,"_ArrayScores.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				
			}
			if("WeightedAll"%in%names(Output[[l]])){
				utils::write.table(t(Output[[l]]$"WeightedAll"), paste(Location,"/",Name,"_WeightedAll.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
#				
#				WeightedAll=do.call("rbind",lapply(Output,function(x) x$WeightedGeneLevelEstimatesAllProbesets))	
#				assign(paste("EventPointer_REIDS_Output_WeightedAll",sep=""),WeightedAll,envir=environment())
#				save(EventPointer_REIDS_Output_WeightedAll, file="EventPointer_REIDS_Output_WeightedAll.RData")	
			}
			if("WeightedAll"%in%names(Output[[l]])){
				utils::write.table(t(Output[[l]]$"WeightedAll"), paste(Location,"/",Name,"_WeightedAll.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
#				
#				WeightedAll=do.call("rbind",lapply(Output,function(x) x$WeightedGeneLevelEstimatesAllProbesets))	
#				assign(paste("EventPointer_REIDS_Output_WeightedAll",sep=""),WeightedAll,envir=environment())
#				save(EventPointer_REIDS_Output_WeightedAll, file="EventPointer_REIDS_Output_WeightedAll.RData")	
			}
			if("EqualAll"%in%names(Output[[l]])){
				utils::write.table(t(Output[[l]]$"EqualAll"), paste(Location,"/",Name,"_EqualAll.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
#				
#				EqualAll=do.call("rbind",lapply(Output,function(x) x$GeneLevelEstimatesAllProbesets))	
#				assign(paste("EventPointer_REIDS_Output_EqualAll",sep=""),EqualAll,envir=environment())
#				save(EventPointer_REIDS_Output_WeightedAll, file="EventPointer_REIDS_Output_EqualAll.RData")
			}
			if("ExonLevel"%in%names(Output[[l]])){
				utils::write.table(Output[[l]]$"ExonLevel", paste(Location,"/",Name,"_ExonLevel.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				
#				ExonLevel=do.call("rbind",lapply(Output,function(x) x$ExonLevelEstimates))	
#				assign(paste("EventPointer_REIDS_Output_ExonLevel",sep=""),ExonLevel,envir=environment())
#				save(EventPointer_REIDS_Output_ExonLevel, file="EventPointer_REIDS_Output_ExonLevel.RData")
			}
			if("WeightedConst"%in%names(Output[[l]])){
				utils::write.table(t(Output[[l]]$"WeightedConst"), paste(Location,"/",Name,"_WeightedConst.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)				
#				WeightedConst=do.call("rbind",lapply(Output,function(x) x$WeightedGeneLevelEstimatesConstProbesets))	
#				assign(paste("EventPointer_REIDS_Output_WeightedConst",sep=""),WeightedConst,envir=environment())
#				save(EventPointer_REIDS_Output_WeightedConst, file="EventPointer_REIDS_Output_WeightedConst.RData")
			}
			if("EqualConst"%in%names(Output[[l]])){
				utils::write.table(t(Output[[l]]$"EqualConst"), paste(Location,"/",Name,"_EqualConst.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)					
#				EqualCons=do.call("rbind",lapply(Output,function(x) x$GeneLevelEstimatesConstProbesets))	
#				assign(paste("EventPointer_REIDS_Output_EqualConst",sep=""),EqualCons,envir=environment())
#				save(EventPointer_REIDS_Output_WeightedConst, file="EventPointer_REIDS_Output_EqualConst.RData")
			}
		}

		file.remove(paste(Location,REIDSOutputFiles,sep="/"))
		
	}
	
	REIDSIsoFiles=list.files(path=Location,pattern="^REIDS_IsoformInfo_")
	if(length(REIDSIsoFiles)>0){
		Output=list()
		for(i in as.character(ID)){
			Data=try(get(load(paste(Location,"/REIDS_IsoformInfo_",as.character(i),".RData",sep=""))),silent=TRUE)
			if(class(Data)!="try-error"){
				Output[length(Output)+1]=Data
				names(Output)[length(Output)]=i
			}
		}
		
		assign(paste(Name,"_REIDS_IsoformInfo",sep=""),Output,envir=environment())
		eval(parse(text=paste("save(",Name, "_REIDS_IsoformInfo, file=\"",Location,"/",Name, "_REIDS_IsoformInfo.RData\")", sep=""))) 
		
		for(l in 1:length(Output)){
			utils::write.table(t(c("TC_ID","PSR_ID","Type","Unreliable Junctions","Linking Junctions","Exclusion Junctions","Supported by","Identified by","Fold Change","Exons")), file = paste(Location,"/",Name,"_ASInfo.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			utils::write.table(t(c("TC_ID","TranscriptName","Junctions","ProbeSets")), file = paste(Location,"/",Name,"_Compositions.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			utils::write.table(t(c("TC_ID","TranscriptName",paste("Present in Group",c(1:length(Groups))))), file = paste(Location,"/",Name,"_GroupTranscripts.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
			
			if(class(Output[[l]]$"TranscriptInformation")!="character"&!is.null(Output[[l]]$"TranscriptInformation")){
				
				if(!is.null(Output[[l]]$"TranscriptInformation"[[1]])){
					utils::write.table(Output[[l]]$"TranscriptInformation"[[1]], paste(Location,"/",Name, "_ASInfo.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				}
				if(!is.null(Output[[l]]$"TranscriptInformation"[[2]])){
					utils::write.table(Output[[l]]$"TranscriptInformation"[[2]], paste(Location,"/",Name, "_Compositions.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)		
					
				}
				if(!is.null(Output[[l]]$"TranscriptInformation"[[3]])){
					utils::write.table(Output[[l]]$"TranscriptInformation"[[3]], paste(Location,"/",Name, "_GroupTranscripts.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				}
				if(!is.null(Output[[l]]$"TranscriptInformation"[[4]])){
					utils::write.table(t(c(names(Output)[l],Output[[l]]$"TranscriptInformation"[[4]])), paste(Location,"/",Name, "_NovelConnections.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				}
#				if(!is.null(Output[[l]]$"IsoformIndication")){
#					utils::write.table(cbind(names(Output)[l],Output[[l]]$"IsoformIndication"),file = paste(Location,"/",Name, "_IsoformIndication.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
#				}
#				if(!is.null(Output[[l]]$"ExonTesting")){
#					utils::write.table(cbind(names(Output)[l],Output[[l]]$"ExonTesting"),file = paste(Location,"/",Name, "_ExonDE.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
#					
#				}
#				if(!is.null(Output[[l]]$"Possible DE Isoforms")){
#					utils::write.table(cbind(names(Output)[l],Output[[l]]$"Possible DE Isoforms"[[2]]),file = paste(Location,"/",Name, "_PossibleDEIsoforms.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)				
#				}
			}
			else{
				if(is.null(Output[[l]]$"TranscriptInformation")){
					utils::write.table(t(c(names(Output)[l],"Transcript not assessed")), paste(Location,"/",Name, "_NotAssessedTranscripts.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)		
					
				}
				else{
					utils::write.table(t(c(names(Output)[l],Output[[l]]$"TranscriptInformation")), paste(Location,"/",Name, "_NotAssessedTranscripts.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)		
				}
			}
		}
	}
}


# REIDS Function - Regular version
# It is advised to run this model on a HPC cluster and not on a regular laptop as it will consume time and memory

#' "reidsfunction_genebygene"
#' 
#' The reidsfunction_genebygene performs the REIDS model and is an internal function of the REIDSFunction. The function calls on the pivot transformed .csv file and transforms the read lines into a data frame on which the REIDS model is performed. 
#' @export
#' @param file_name The name of the pivot transformed .csv file.
#' @param file_pos The position in the file where to start reading.
#' @param line_length The length of the line to read.
#' @param ASPSR A vector with alternatively spliced probe sets which are taken out of the analysis and summarization. This is useful if Summarize is "WeightedConst" and/or "EqualConst".
#' @param Summarize A character vector specifying wich summarization method is to be performed. The choices are "EqualAll", "WeightedAll", "EqualConst" and "WeightedConst". The former two use all probe sets while the latter use only the constituitive probe sets. Summarization on the constistuitive probe sets will only be performed if ASPSR is specified.
#' @param nsim The number of iterations to perform. Defaults to 1000.
#' @param informativeCalls Logical. Should the I/NI calls method be perform before applying the REIDS model?
#' @param rho The threshold for filtering in the I/NI calls method. Probesets with scores higher than rho are kept.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @param Location A character string indication the place where the output should be saved.
#' @param Name A name for the output to be saved at Location. Defaults to "REIDS". 
#' @return The functions writes the obtained information to .txt files: "Name_INICalls.txt", Name_ExonScores.txt", "Name_ArrayScores.txt" and the summarized values distributed across "Name_WeightedAll.txt", "Name_EqualAll.txt", "Name_WeightedConst.txt" and "Name_EqualConst.txt". 
reidsfunction_genebygene <- function(file_name,file_pos,line_length,ASPSR=c(),nsim=1000,informativeCalls=TRUE,Summarize=FALSE,rho=0.5,Low_AllSamples=c(),Location=NULL,Name="REIDS"){
	file_pos=as.numeric(as.matrix(file_pos))
	line_length=as.numeric(as.matrix(line_length))
	
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
		geneData=cbind(geneData,TempData)

		if(length(ASPSR)>0){
			AS=which(geneData[,2]%in%ASPSR)
			if(length(AS)>0){
				geneData=geneData[-c(AS),,drop=FALSE]
			}
		}
		
		DABGs=which(geneData[,2]%in%Low_AllSamples)
		if(length(DABGs)>0){
			geneData=geneData[-c(which(geneData[,2]%in%Low_AllSamples)),,drop=FALSE]
		}
		
		Juncs=which(sapply(geneData[,2],function(x) substr(x,1,3))=="JUC")
		if(length(Juncs)>0){
			geneData=geneData[-c(which(sapply(geneData[,2],function(x) substr(x,1,3))=="JUC")),,drop=FALSE]
		}
		
		exonScore <- arrayScore <- informativeData<- NULL
		
		output=list()
		output[[1]]=list()
		names(output)[1]=geneID
		lcmmData <- geneData[,-c(1,2)]
		lcmmData=as.matrix(lcmmData)
		enames <- geneData$exonID
		names(enames)=NULL
		
		rownames(lcmmData)<- enames
		
		##informative calls
		i=1
		if(informativeCalls&length(ASPSR)==0){			
			fit <- iniREIDS(SubgeneData=lcmmData, nsim) 
			fit2 <- data.frame(geneID=geneID,exonNames=unique(enames),Score=fit,informative=fit>rho) # iniREIDS returns one value per exon: filtering on exon level,no replicates
			output[[1]][[i]]=fit2
			names(output[[1]])[i]="Informative"
			iniData <- lcmmData[which(rownames(lcmmData)%in%fit2$exonNames[fit2$informative]),] # Of those that pass filtering step, retrieve the replicates and samples
			lcmmData <-  iniData   
			i=i+1
			
			if(!is.null(Location)){
				utils::write.table(fit2, paste(Location,"/",Name,"_INICalls.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
			}
		}
		
			
		if(!is.null(lcmmData)&length(unique(rownames(lcmmData)))>=3){  
			fit <- REIDSmodel_intern(SubgeneData=lcmmData, nsim) 
			
			exonScore <-fit$exonScores
			arrayScore <- fit$arrayScores
			
			if(length(ASPSR)==0){
				output[[1]][[i]]=cbind(geneID,rownames(arrayScore),exonScore)
				names(output[[1]])[i]="exonScore"
				utils::write.table(cbind(geneID,rownames(arrayScore),exonScore), paste(Location,"/",Name,"_ExonScores.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
			
				output[[1]][[i+1]]=cbind(geneID,rownames(arrayScore),arrayScore)
				names(output[[1]])[i+1]="arrayScore"
				i=i+2
				if(!is.null(Location)){
					utils::write.table(cbind(geneID,rownames(arrayScore),arrayScore), paste(Location,"/",Name,"_ArrayScores.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				}
			}
			
			if(any(c("EqualAll","WeightedAll")%in%Summarize)&length(ASPSR)==0){
				exonVariances<-fit$exonVar	
				probeEffects<-fit$probeEffects
				sampleEffects<-fit$sampleEffects
				estimatedvalues<-fit$arrayScores
				
				if("WeightedAll"%in%Summarize){## with weights
					TotalExonVar=sum(1/exonVariances)
					Weights=(1/exonVariances)/TotalExonVar
					names(Weights)=unique(rownames(lcmmData))
					
					arrayLevels=apply(estimatedvalues,2,function(j) sum(Weights*j))
					
					WeightedGeneLevelEstimates=mean(probeEffects)+sampleEffects+arrayLevels
					
					output[[1]][[i]]=Weights
					names(output[[1]])[i]="WeightsAllProbesets"
					
					output[[1]][[i+1]]=c(geneID,WeightedGeneLevelEstimates)
					names(output[[1]])[i+1]="WeightedAll"
					if(!is.null(Location)){
						utils::write.table(t(c(geneID,WeightedGeneLevelEstimates)), paste(Location,"/",Name,"_WeightedAll.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
					}
				}
				if("EqualAll"%in%Summarize){
					
					## without weights		
					arrayLevels=apply(estimatedvalues,2,mean)
					GeneLevelEstimates=mean(probeEffects)+sampleEffects+arrayLevels
					
					output[[1]][[i+2]]=c(geneID,GeneLevelEstimates)
					names(output[[1]])[i+2]="EqualAll"
					if(!is.null(Location)){
						utils::write.table(t(c(geneID,GeneLevelEstimates)), paste(Location,"/",Name,"_EqualAll.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
					}
				}
				i=i+3
				
				#ExonLevelValues
				
				ProbeLevel=matrix(probeEffects)
				rownames(ProbeLevel)=geneData[,2]
				PSR=unique(geneData[,2])
				ExonLevelEstimates=matrix(0,nrow=length(PSR),ncol=length(sampleEffects))
				for(e in 1:length(PSR)){
					E=mean(ProbeLevel[which(rownames(ProbeLevel)==PSR[e]),])+sampleEffects+estimatedvalues[which(rownames(estimatedvalues)==PSR[e]),]
					ExonLevelEstimates[e,]=E
				}
				rownames(ExonLevelEstimates)=PSR
				
				output[[1]][[i]]=cbind(geneID,PSR,ExonLevelEstimates)
				names(output[[1]])[i]="ExonLevel"
				i=i+1
				if(!is.null(Location)){
					utils::write.table(cbind(geneID,PSR,ExonLevelEstimates), paste(Location,"/",Name,"_ExonLevel.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				}
				
			}
			if(any(c("EqualConst","WeightedConst")%in%Summarize)&length(ASPSR)>0){
				exonVariances<-fit$exonVar	
				names(exonVariances)=unique(rownames(lcmmData))
				probeEffects<-fit$probeEffects
				sampleEffects<-fit$sampleEffects
				estimatedvalues<-fit$arrayScores
				
				if("WeightedConst"%in%Summarize){## with weights
					TotalExonVar=sum(1/exonVariances)
					Weights=(1/exonVariances)/TotalExonVar
					names(Weights)=unique(rownames(lcmmData))
					
					arrayLevels=apply(estimatedvalues,2,function(j) sum(Weights*j))
					
					WeightedGeneLevelEstimates=mean(probeEffects)+sampleEffects+arrayLevels
					
					output[[1]][[i]]=Weights
					names(output[[1]])[i]="WeightsConstProbesets"
					
					output[[1]][[i+1]]=c(geneID,WeightedGeneLevelEstimates)
					names(output[[1]])[i+1]="WeightedConst"
					if(!is.null(Location)){
						utils::write.table(t(c(geneID,WeightedGeneLevelEstimates)), paste(Location,"/",Name,"_WeightedConst.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
					}
				}
				if("EqualConst"%in%Summarize){
					
					## without weights		
					arrayLevels=apply(estimatedvalues,2,mean)
					GeneLevelEstimates=mean(probeEffects)+sampleEffects+arrayLevels
					
					output[[1]][[i+2]]=c(geneID,GeneLevelEstimates)
					names(output[[1]])[i+2]="EqualConst"
					if(!is.null(Location)){
						utils::write.table(t(c(geneID,GeneLevelEstimates)), paste(Location,"/",Name,"_EqualConst.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
					}
				}
				i=i+3
				
				#ExonLevelValues
				
				ProbeLevel=matrix(probeEffects)
				rownames(ProbeLevel)=geneData[,2]
				PSR=unique(geneData[,2])
				ExonLevelEstimates=matrix(0,nrow=length(PSR),ncol=length(sampleEffects))
				for(e in 1:length(PSR)){
					E=mean(ProbeLevel[which(rownames(ProbeLevel)==PSR[e]),])+sampleEffects+estimatedvalues[which(rownames(estimatedvalues)==PSR[e]),]
					ExonLevelEstimates[e,]=E
				}
				rownames(ExonLevelEstimates)=PSR
				
				output[[1]][[i]]=cbind(geneID,PSR,ExonLevelEstimates)
				if(length(AS)==0){
					names(output[[1]])[i]="ExonLevel"
					i=i+1
				}
				else{
					names(output[[1]])[i]="ExonLevel_Const"
					i=i+1
				}				
			}
		}
		
		#save(output, file=paste(Location,"/REIDS_Gene_", geneID, ".RData", sep=""))
		#rm(output)
		#gc()
	
		close(conn)
		if(!is.null(Location)){
			return(output)
		}
	}
	else{
		close(conn)
	}
}


#' "REIDSFunction"
#' 
#' The REIDSFunction performs the REIDS model on the pivot transformed data by calling on the line indexed file. The REIDS model is performed gene by gene and the returned outputs are knitted together.
#' @export
#' @param ASPSR A vector with alternatively spliced probe sets which are taken out of the analysis and summarization. This is useful if Summarize is "WeightedConst" and/or "EqualConst".
#' @param Indices The .csv file created by Line_Indexer.py which contains indices for every gene in geneIDs.
#' @param DataFile The .csv file created by PivotTransformation. 
#' @param nsim The number of iterations to perform. Defaults to 1000.
#' @param informativeCalls Logical. Should the I/NI calls method be perform before applying the REIDS model?
#' @param Summarize A character vector specifying wich summarization method is to be performed. The choices are "EqualAll", "WeightedAll", "EqualConst" and "WeightedConst". The former two use all probe sets while the latter use only the constituitive probe sets. Summarization on the constistuitive probe sets will only be performed if ASPSR is specified.
#' @param rho The threshold for filtering in the I/NI calls method. Probesets with scores higher than rho are kept.
#' @param Groups A list with elements specifying the columns of the data in each group.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @param Location A character string indication the place where the outputs are saved.
#' @param Name A name for the output to be saved at Location. Defaults to "REIDS". 
#' @return The functions writes the obtained information to .txt files: "Name_INICalls.txt", Name_ExonScores.txt", "Name_ArrayScores.txt" and the summarized values distributed across "Name_WeightedAll.txt", "Name_EqualAll.txt", "Name_WeightedConst.txt" and "Name_EqualConst.txt". 
#' @examples
#' \dontrun{
#' data(TC1500264)
#' PivotTransformData(Data=TC1500264,GeneID=NULL,ExonID=NULL,
#' REMAPSplitFile="TC1500264_Gene_SplitFile.txt",Location="Output/",Name="TC1500264_Pivot")
#' 
#' REIDSFunction(ASPSR=c(), Indices="Output/TC1500264_LineIndex.csv",
#' DataFile="Output/TC1500264_Pivot.csv",nsim=50,informativeCalls=FALSE,
#' Summarize=c("WeightedAll","EqualAll"),
#' rho=0.5,Low_AllSamples=c(),Groups=list(c(1:3),c(4:6)),Location="Output",Name="TC1500264")
#' }
REIDSFunction<-function(ASPSR=c(),Indices,DataFile,nsim=1000,informativeCalls=TRUE,Summarize=FALSE,rho=0.5,Low_AllSamples=c(),Groups,Location=NULL,Name="REIDS"){
	Lines=utils::read.table(Indices,header=TRUE,sep=",",stringsAsFactors=FALSE)
	Lines=data.frame(Lines)	
	Lines[,1]=as.numeric(Lines[,1])
	Lines[,2]=as.numeric(Lines[,2])

	if(!is.null(Location)){
		Files=list.files(path=Location,pattern=paste("^",Name,sep=""))	
		if(!paste(Name,"_INICalls.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID","PSR_ID","Informative")), file = paste(Location,"/",Name,"_INICalls.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_ExonScores.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID","PSR_ID","ExonScore")), file = paste(Location,"/",Name,"_ExonScores.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_ArrayScores.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID","PSR_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_ArrayScores.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_WeightedAll.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_WeightedAll.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_EqualAll.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_EqualAll.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_WeightedConst.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_WeightedConst.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_EqualConst.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_EqualConst.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_ExonLevel.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID","PSR_ID",paste("Sample",c(1:length(unlist(Groups)))))), file = paste(Location,"/",Name,"_ExonLevel.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
	}
	REIDSOut=apply(Lines,1, function(x) reidsfunction_genebygene(file_name=DataFile,file_pos=x[1],line_length=x[2],ASPSR,nsim,informativeCalls,Summarize,rho,Low_AllSamples,Location,Name))
	
	#CreateOutput(ID=geneIDs,Name,Location)
}

##### OUTPUT ANALYSIS FUNCTIONS #####


#' "ExonTesting"
#' 
#' The ExonTesting function performs a t-test (2 groups) or F-test (more than 2 groups) between the array scores of predefined groups. If specified, probe sets are filtered out on exon scores and significance level. The function is the internal function of ASExons.
#' @export
#' @param ExonScores The path to the file with the exon scores of the probe sets.
#' @param ArrayScores The path to the file with the array scores of the probe sets.
#' @param Exonthreshold The exon score threshold to be maintained. If not NULL, probe sets with an exon score lower than this value are not considered further and the p-values will be adjusted for multiplicity after testing. If NULL, all probesets are considered and a multiplicity correction is not performed.
#' @param Groups A list with elements specifying the columns of the data in each group.
#' @param paired Logical. Are the groups paired? If TRUE the mean paired differences are calculated and tested whether these are significantly different from zero or not.
#' @param significancelevel The significance level to be maintained on the p-values. The filtering on the significance is conducted only if an Exonthreshold is specified and the p-value are adjusted for multiplicity.
#' @return A data frame with one line per exon. The columns contain the gene ID, the exon ID, the exon score, the test statistic, a p-value and an adjusted p-value. If the groups are paired also the mean paired difference is given. The p-values are adjusted for multiplicity and filtered on significance if significancelevel is not NULL.
ExonTesting <- function(ExonScores,ArrayScores,Exonthreshold=NULL,Groups=list(),paired=FALSE,significancelevel=NULL){

	if(is.null(Exonthreshold)){
		Exonthreshold=0
	}
	
	## Step 1 : filtering on the Exon value. If 0, no filtering occurs
	SubsetExonScores=ExonScores[which(round(ExonScores[,3],2)>Exonthreshold),,drop=FALSE]
	Exons=unique(SubsetExonScores[,2])
	
	## Step 2 : Testing of the Array Scores -- paired or not paired
	ArrayScoreTTest=matrix(0,nrow=length(Exons),ncol=3)
	colnames(ArrayScoreTTest)=c("statistic","p.value","adj.p.value")
	
	ttest<-function(data,groups,pairs){
		if(length(groups)==2){
			C=c(data[,groups[[1]]],data[,groups[[2]]])
			l1=C>0&C<0.5
			l2=C<0&C>-0.5
			l=l1+l2
			if(length(which(l==1))>=(0.80*length(l))){
				out2=c(0,1)
			}
			else{
				#out1=ks.test(x=as.numeric(data[,groups[[1]]]),y=as.numeric(data[,groups[[2]]]))
				out1=stats::t.test(x=data[,groups[[1]]],y=data[,groups[[2]]],paired=pairs)
				out2=cbind(out1$statistic,out1$p.value)
				return(out2)	
			}
		}
		else{
			C=c(data[,groups[[1]]])
			l=C>0&C<0.5
			if(length(which(l==1))>=(0.80*length(l))){
				out2=c(0,1)
			}
			else{
				out1=stats::t.test(x=data[,groups[[1]]])
				out2=cbind(out1$statistic,out1$p.value)
				return(out2)
			}
		}
	}
	
	anovaFtest<-function(data,Groups){
		fit<-stats::aov(as.vector(data)~Groups)
		out1=cbind(summary(fit)[[1]][1,4],summary(fit)[[1]][1,5])	
		return(out1)
		
	}
	
	
	if(paired==FALSE){  # Test between two groups of Array Scores
		
		names(Groups)=c(1:length(Groups))

		if(length(Groups)<=2){
			ArrayScoreTTest[,c(1,2)] = t(sapply(Exons,function(i) {ttest(data=ArrayScores[which(ArrayScores[,2]==i),-c(1,2)],groups=Groups, pairs = paired)}))	
		}
		else{
			groups=rep(0,length(unlist(Groups)))
			for(j in 1:length(Groups)){
				positions=Groups[[j]]
				groups[positions]=names(Groups)[j]
			}
			groups=as.numeric(groups)
			groups=factor(groups)
			
			ArrayScoreTTest[,c(1,2)] = t(sapply(Exons,function(i) {anovaFtest(data=ArrayScores[which(ArrayScores[,2]==i),-c(1,2),drop=FALSE],Groups=groups)}))	
			
		}
		p_vals=as.vector(as.matrix((ArrayScoreTTest[,grep("^(p.value)",colnames(ArrayScoreTTest))])))
		adj_p_vals=matrix(stats::p.adjust(p_vals,"BH"),nrow=nrow(ArrayScoreTTest),ncol=length(grep("^(p.value)",colnames(ArrayScoreTTest))))
		ArrayScoreTTest[,which(seq(1,ncol(ArrayScoreTTest))%%3==0)]=adj_p_vals
		
		ArrayScoreTTest=as.data.frame(ArrayScoreTTest)
		Out=cbind(SubsetExonScores,ArrayScoreTTest)
		
	}
	
	else if(paired==TRUE){
		if(!(is.null(groups[[1]])) & length(Exons)!=0){
			
			ArrayScore_group1=ArrayScores[which(ArrayScores[,2]%in%Exons),groups[[1]]+2,drop=FALSE]
			ArrayScore_group2=ArrayScores[which(ArrayScores[,2]%in%Exons),groups[[2]]+2,drop=FALSE]
			
			mean_paired_diff<-function(g1,g2){
				Paired_Diff=g1-g2
				Mean_Diff=mean(Paired_Diff)
				
				out1=stats::t.test(x=Paired_Diff)
				out2=cbind(out1$statistic,out1$p.value)
				out3=cbind(Mean_Diff,out2)
				
				return(out3)
			}
			
			ArrayScoreTest = t(sapply(c(1:length(Exons)),function(i) mean_paired_diff(g1=ArrayScore_group1[i,],g2=ArrayScore_group2[i,]) ))
			ArrayScoreTest=data.frame("Mean_Diff"=ArrayScoreTest[,1],"t-statistic"=ArrayScoreTest[,2],p.value=ArrayScoreTest[,3])
			rownames(ArrayScoreTest)=rownames(ArrayScore_group1)
			Out=cbind(SubsetExonScores,ArrayScoreTTest)
			
		}
		else{
			Out= NULL
		}
	}
	
	if(is.null(Out)){
		Out=SubsetExonScores
	}

	colnames(Out)[1]="GeneID"
	colnames(Out)[2]="PSR_ID"
	colnames(Out)[3]="ExonScore"
	
	if(nrow(Out)==0){
		Out=NULL
		return(Out)
	}
	
	if(!is.null(significancelevel)&!is.null(Out)){
		Out=Out[which(Out$adj.p.value<=significancelevel),,drop=FALSE]
	}
	if(nrow(Out)>0){
		rownames(Out)=seq(1:nrow(Out))
	}
	return(Out)
}


## Identification of the AS exons

#' "ASExons"
#' 
#' The ASExons functions alternatively spliced exons from the exon scores and array scores. It filters probesets on their exon scores, adjusts p-values for multiplicity and only keeps the significant probesets.
#' @export
#' @param ExonScores The path to the file with the exon scores of the probe sets.
#' @param ArrayScores The path to the file with the array scores of the probe sets.
#' @param Exonthreshold The exon score threshold to be maintained. If not NULL, probe sets with an exon score lower than this value are not considered further and the p-values will be adjusted for multiplicity after testing. If NULL, all probesets are considered and a multiplicity correction is not performed.
#' @param Groups A list with elements specifying the columns of the data in each group.
#' @param paired Logical. Are the groups paired? If TRUE the mean paired differences are calculated and tested whether these are significantly different from zero or not.
#' @param significancelevel The significance level to be maintained on the p-values. The filtering on the significance is conducted only if an Exonthreshold is specified and the p-value are adjusted for multiplicity.
#' @param Location A character string indication the place where the outputs are saved.
#' @param Name A character string with the name of the ouput file. Defaults to "REIDSAS".
#' @return A data frame with one line per exon. The columns contain the gene ID, the exon ID, the exon score the test statistic, a p-value and an adjusted p-value. If the groups are paired also the mean paired difference is given. Only the probesets with high enough exon scores and a significant test are kept in the data frame.
#' @examples 
#' \dontrun{
#' data(TC1500264)
#' 
#' PivotTransformData(Data=TC1500264,GeneID=NULL,ExonID=NULL,
#' REMAPSplitFile="TC1500264_Gene_SplitFile.txt",Location="Output/",Name="TC1500264_Pivot")
#' 
#' REIDSFunction(ASPSR=c(), Indices="Output/TC1500264_LineIndex.csv",
#' DataFile="Output/TC1500264_Pivot.csv",nsim=50,informativeCalls=FALSE,Summarize=
#' c("WeightedAll","EqualAll"),rho=0.5,Low_AllSamples=c(),Groups=list(c(1:3),c(4:6)),
#' Location="Output",Name="TC1500264")
#' 
#' TC1500264_1vs2=ASExons(ExonScores="Output/TC1500264_ExonScores.txt",ArrayScores=
#' "Output/TC1500264_ArrayScores.txt",Exonthreshold=0.5,Groups=list(c(1:3),c(4:6)),
#' paired=FALSE,significancelevel=0.05)
#' }
ASExons<-function(ExonScores,ArrayScores,Exonthreshold=0.5,Groups=list(group1=NULL,group2=NULL),paired=FALSE,significancelevel=0.05,Location=NULL,Name="REIDSAS"){
	
	message("The used threshold for the exon scores is 0.5")
	message("The used significance level for the p-values is 0.05")
	
	ExonScores=utils::read.table(ExonScores,sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses=c("character","character","numeric"))
	
	ArrayScores=utils::read.table(ArrayScores,header=TRUE,sep="\t",stringsAsFactors=FALSE,colClasses=c("character","character",rep("numeric",length(unlist(Groups)))))
		
	TestedData=ExonTesting(ExonScores,ArrayScores,Exonthreshold=Exonthreshold,Groups=Groups,paired=paired,significancelevel=significancelevel)
		
	if(nrow(TestedData)!=0){
		message("Ordering data in from high to low significance")
		Data_Sign_Ordered=TestedData[order(TestedData[,ncol(TestedData)]),]
		rownames(Data_Sign_Ordered)=c(1:nrow(Data_Sign_Ordered))
	}
	else{
		Data_Sign_Ordered=NULL
	}

	if(!is.null(Location)){
		utils::write.table(t(c("TC_ID","PSR_ID","ExonScore","Statistic","Pvalue","Adj.PValue")), file = paste(Location,"/",Name,"_ASTesting.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)		
		utils::write.table(Data_Sign_Ordered,file=paste(Location,"/",Name,"_ASTesting.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
	}
	else{
		return(Data_Sign_Ordered)
	}
	
}

#' JunInfo
#' 
#' JunInfo functions asses the junction information for a single gene
#' @export
#' @param file_name The name of the pivot transformed .csv file.
#' @param file_pos The position in the file where to start reading.
#' @param line_length The length of the line to read.
#' @param ASPSR The AS probe sets as identified by ASExons.
#' @param Juninfo A parameter specifying wether the annotations are user of Ensembl defined. If JunInfo is "User" (default) the annotations provided in EandTrAnnot are used. If JunInfo is "Ensembl" the annotations in EandTrAnnot are used to set up tje junction associations but the gene name and position in transcriptData and positionData are used to connect with the Ensembl data base and retrieve corresponding information. 
#' @param JAnnotI The file name with line indices for the junction associations.
#' @param JAnnot The file name with the junction associations.
#' @param EandTrAnnotI The file name with line indices for the exon and isoform annotations.
#' @param EandTrAnnot The file name with the exon and isoform annotations.
#' @param PartiallyAnnotated Logical. Should the exon annotations with partially annotated probe sets still be included? If FALSE, these are excluded. If TRUE, these are included. Default is FALSE.
#' @param positionData The file with the chromosome start and ends for the probe sets. Only needed in JunInfo=Ensembl.
#' @param transcriptData The file with gene name of the transcripts. Only needed in JunInfo=Ensembl.
#' @param Groups A list with  elements speficifing the columns of the data in each group.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @param Low_GSamples A list with a  character vector per group containing the probe sets which are not DABG in that group.
#' @param Plot Should a plot of the gene model be made?
#' @param Location A character string indication the place where the outputs are saved.
#' @param Name A character string with the name of the ouput file.
#' @details The plot is produced by the arcplot function of the arcdiagram package (https://github.com/gastonstat/arcdiagram)
#' @return The function returns four files. The first file has name "Name_ASInfo.txt" and contains a line per probe set. It shows the reached decision regarding the probe set (Const/AS/not DABG),its linking and exclusion junctions, the fold change, the AS type and its annotated exons. The second file, "Name_Compositions.txt", is a list of all found transcripts for a particular TC ID. The third file,"Name_GroupTranscripts.txt" indicates whether a specific transcript is present or absent in a group. The fourth file "Name_NovelConnections.txt" contains junctions which are showing an undocumented connection between probe sets.
JunInfo<-function(file_name,file_pos,line_length,ASPSR=c(),Juninfo="User",JAnnotI,JAnnot=NULL,EandTrAnnotI=NULL,EandTrAnnot=NULL,PartiallyAnnotated=FALSE,positionData=NULL,transcriptData=NULL,Groups=list(),Low_AllSamples=c(),Low_GSamples=c(),Plot=FALSE,Location=NULL,Name=""){

	file_pos=as.vector(as.matrix(file_pos))
	line_length=as.vector(as.matrix(line_length))
	
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
	
	
		if(!is.null(EandTrAnnotI)){
			ETrI=utils::read.table(EandTrAnnotI,header=FALSE,stringsAsFactors=FALSE)
			Lines=ETrI[which(ETrI[,1]==geneID),]
			ETrAnnot=utils::read.table(EandTrAnnot,header=FALSE,sep="\t",nrows=as.numeric(Lines[3]),skip=as.numeric(Lines[2]),stringsAsFactors=FALSE)
			TrAnnot=ETrAnnot[,c(1,2,8,4,7)]
			TrAnnot=TrAnnot[order(TrAnnot[,2]),]
			colnames(TrAnnot)=c("TC_ID","PSR_ID","strand","EAnnot","TrAnnot")
			
			if(!PartiallyAnnotated){
				Incl=unlist(sapply(TrAnnot[,4],function(x) !grepl("[*]",x)))
				TrAnnot=TrAnnot[Incl,]
			}
			
			
			if(!is.null(JAnnotI)){
				JI=utils::read.table(JI,header=FALSE)
				Lines=JI[which(JI[,1]==geneID),]
				JAnnot=utils::read.table(EandTrAnnot,header=FALSE,sep="\t",nrows=as.numeric(Lines[3]),skip=as.numeric(Lines[2]),stringsAsFactors=FALSE)		
			}
			else{
				
				JEAnnot=ETrAnnot[,-c(5,6,7)]
				JEAnnot=JEAnnot[!duplicated(JEAnnot),]
				JEAnnot=JEAnnot[order(JEAnnot[,2]),]
				JUCs=which(sapply(JEAnnot[,2],function(x) substr(x,1,1)=="J"))
				JUC=JEAnnot[JUCs,]
				PSR=JEAnnot[-c(JUCs),]
				
				
				Incl=unlist(sapply(JUC[,4],function(x) !grepl("[*]",x)))
				JUC=JUC[Incl,]
				
				if(!PartiallyAnnotated){
					Incl=unlist(sapply(PSR[,4],function(x) !grepl("[*]",x)))
					PSR=PSR[Incl,]
				}
				
				
				J=unique(JUC[,2])
				JAnnot=c()
				for(j in J){
					
					PSR3_Final=""
					PSR5_Final=""
					
					#Side 3 annots
					SubJ=JUC[which(JUC[,2]==j&JUC[,3]==3), ]
					TC=SubJ[1,1]
					Strand=SubJ[1,5]
					
					Es=unique(SubJ[,4])
					PSRs3=PSR[which(PSR[,4]%in%Es), ]
					
					Count=table(PSRs3[,2])
					PSR3=names(Count)[which.min(table(PSRs3[,2])-length(Es))]
					PSRs3temp=PSR[which(PSR[,2]%in%PSR3), 2]
					
					if((Strand=="-"|Strand==-1|Strand=="-1.0") & length(PSRs3temp)>0){
						PSR3_Final=PSRs3temp[1]
					}			
					else if((Strand=="+"|Strand==1|Strand=="1.0") & length(PSRs3temp)>0){
						PSR3_Final=PSRs3temp[length(PSRs3temp)]
					}	
					Row=c(TC,PSR3_Final,j,"3")
					JAnnot=rbind(JAnnot,Row)
					
					#Side 5 annots
					SubJ=JUC[which(JUC[,2]==j&JUC[,3]==5), ]
					TC=SubJ[1,1]
					Strand=SubJ[1,5]
					
					Es=unique(SubJ[,4])
					PSRs5=PSR[which(PSR[,4]%in%Es), ]
					
					Count=table(PSRs5[,2])				
					PSR5=names(Count)[which.min(table(PSRs5[,2])-length(Es))]
					PSRs5temp=PSR[which(PSR[,2]%in%PSR5), 2]
					
					if((Strand=="-"|Strand==-1|Strand=="-1.0") & length(PSRs5temp)>0){
						PSR5_Final=PSRs5temp[length(PSRs5temp)]
					}			
					else if((Strand=="+"|Strand==1|Strand=="1.0") & length(PSRs5temp)>0){
						PSR5_Final=PSRs5temp[1]
					}	
					Row=c(TC,PSR5_Final,j,"5")
					JAnnot=rbind(JAnnot,Row)
					
					#exclusions
					I1=0
					I2=0
					if(PSR3_Final!="" & PSR5_Final!="" & (Strand=="+"|Strand==1|Strand=="1.0")){
						I1=which(PSR[,2]==PSR3_Final)[length(which(PSR[,2]==PSR3_Final))]
						I2=which(PSR[,2]==PSR5_Final)[1]
					}
					else if(PSR3_Final!="" & PSR5_Final!="" & (Strand=="-"|Strand==-1|Strand=="-1.0")){
						I1=which(PSR[,2]==PSR5_Final)[length(which(PSR[,2]==PSR5_Final))]
						I2=which(PSR[,2]==PSR3_Final)[1]
					}
					if(I1!=0&I2!=0&I1!=(I2+1)&I2!=(I1+1)){
						if(I1<I2){
							R=seq(I1+1,I2,1)
						}
						else if(I1>I2){
							R=seq(I1-1,I2,-1)
						}
						ExclPSR_temp=PSR[R,]
						Exrows=nrow(ExclPSR_temp)
						
						if(Exrows!=0){
							ExclPSR=unique(ExclPSR_temp[,2])
							if(any(is.na(ExclPSR))){
								ExclPSR=ExclPSR[-c(which(is.na(ExclPSR)))]
							}
							for(e in ExclPSR){
								if(e!=PSR5_Final & e!=PSR3_Final){
									Row=c(TC,e,j,"exclusion")
									JAnnot=rbind(JAnnot,Row)
								}
								
							}
						}
					}
					
				}
				if(!is.null(JAnnot)){
					JAnnot=JAnnot[order(JAnnot[,2]),]
					colnames(JAnnot)=c("TC_ID","PSR_ID","JUC_ID","as_type")			
					JAnnot[,1]=as.character(JAnnot[,1])
					JAnnot[,2]=as.character(JAnnot[,2])
					JAnnot[,3]=as.character(JAnnot[,3])
					JAnnot[,4]=as.character(JAnnot[,4])
				}
				else{
					JAnnot=c()
				}
			}
			if(Juninfo=="Ensemble"){
				ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
				attributes.region<- c("chromosome_name", "start_position", "end_position", "ensembl_gene_id")							
				filter.symbol<- "hgnc_symbol"
				
				TrAnnot=c()
				
				trim <- function(s, ...) {
					s <- as.character(s);
					s <- sub("^[\t\n\f\r[:punct:]  ]*", "", s);
					s <- sub("[\t\n\f\r ]*$", "", s);
					s;
				} 
				
				trans=paste(geneID,".hg",sep="")
				exons=positionData[positionData$transcript_cluster_id == trans, "probeset_id"]
				
				PSR=sapply(exons, function(x) substr(x,1,3)!="JUC")
				exons=exons[PSR]
				trans=paste(trans,".1",sep="")
				symbol.to.annotate <- AnnotateGenes(trans,transcriptData)$symbol # get HBC identity
				if(all(symbol.to.annotate!="---")){
					ensembl.output     <- AnnotateGeneSymbol(symbol.to.annotate) # get region
					if(nrow(ensembl.output) > 0){
						gene  = makeGene(id = ensembl.output$ensembl_gene_id, biomart = ensembl)
						
						strand <- transcriptData[transcriptData[,1] == trans,][2]
						gene.positions   <- transcriptData[transcriptData$probeset_id %in% exons,]
						gene.positions$start=as.numeric(gene.positions$start)
						gene.positions$stop=as.numeric(gene.positions$stop)
						gene.positions<-gene.positions[order(gene.positions$start,decreasing=FALSE),]
						gene.positions[,1]=sapply(gene.positions[,1],function(x) strsplit(x,"[.]")[[1]][1])
						PSR=gene.positions[,1]
						G=gene@ens
						for(p in PSR){
							Start=gene.positions[which(gene.positions[,1]==p),3]
							Stop=gene.positions[which(gene.positions[,1]==p),4]
							for(r in 1:nrow(G)){
								if(G[r,4]<=Start&G[r,5]>=Stop){
									Row=c(geneID,p,strand,G[r,3],G[r,2])
									TrAnnot=rbind(TrAnnot,Row)
								}
							}
						}
						
						if(is.null(TrAnnot)){
							TrAnnot=c()
						}
						else{
							colnames(TrAnnot)=c("TC_ID","PSR_ID","strand","EAnnot","TrAnnot")
							
							TrAnnot=as.data.frame(TrAnnot)
							TrAnnot[,1]=as.character(TrAnnot[,1])
							TrAnnot[,2]=as.character(TrAnnot[,2])
							TrAnnot[,3]=as.character(TrAnnot[,3])
							TrAnnot[,4]=as.character(TrAnnot[,4])
							TrAnnot[,5]=as.character(TrAnnot[,5])
							TrAnnot=as.matrix(TrAnnot)
						}
						
					}
					else{
						TrAnnot=c()
					}
				}
				else{
					TrAnnot=c()
				}
			}
			
		}
		else{
			ETrAnnot=c()
			JAnnot=c()
		}
		
		if(length(JAnnot)==0|length(TrAnnot)==0){
			
			output=paste(geneID,"No valuable isoform composition information available",sep="\t")
			if(!is.null(Location)){
				utils::write.table(output, paste(Location,"/",Name,"_NotAssessedTranscripts.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				#save(output, file=paste(Location,"/REIDS_IsoformInfo_", geneID, ".RData", sep=""))
				#rm(output)
				#gc()
				close(conn)
				return(paste(geneID," Completed",sep=""))
			}
			else{
				close(conn)
				return(output)
			}
		}
		
		
		colnames(DataS)[2]="ExonID"
		if(any(is.na(JAnnot[,2]))){
			D=which(is.na(JAnnot[,2]))
			JAnnot=JAnnot[-c(D),]
		}
		
		
		OutRange=matrix(0,nrow=length(unique(DataS$ExonID)),ncol=length(Groups))
		for(i in 1:length(unique(DataS$ExonID))){
			for(j in 1:length(Groups)){
				Subset=as.vector(as.matrix(DataS[which(DataS$ExonID==unique(DataS$ExonID)[i]),Groups[[j]]+2]))
				Range=range(Subset)
				OutRange[i,j]=Range[2]-Range[1]
			}	
		}
		rownames(OutRange)=unique(DataS$ExonID)
		
		PSRsExons=unique(DataS$ExonID)[which(substr(unique(DataS$ExonID),1,3)=="PSR")]
		Continue=TRUE
		TempRange=as.vector(as.matrix(OutRange))
		
		while(Continue){
			M=stats::median(TempRange)
			SD=stats::sd(TempRange)
			
			RangeCutOff=M+1.5*SD
			
			if(length(which(TempRange>RangeCutOff))==0){
				Continue=FALSE
			}	
			else{
				TempRange=TempRange[-c(which(TempRange>RangeCutOff))]
			}	
			
			if(RangeCutOff<2){
				RangeCutOff=2
			}
		}
		
		GroupData=list()
		DelJ=c()
		Rem=c()
		for(g in 1:length(Groups)){
			
			GData=DataS[,c(1,2,Groups[[g]]+2)]
			
			for(e in unique(GData$ExonID)){
				Subset=GData[which(GData$ExonID==e),-c(1,2)]
				Temp=as.vector(as.matrix(GData[which(GData$ExonID==e),-c(1,2)]))
				
				StopDel=c()
				if(length(Temp)<=12){
					StopDel=0
					next
				}
				else if(length(Temp)<18){
					StopDel=2
				}
				else{
					StopDel=3
				}
				
				if(round(OutRange[e,g],2)>RangeCutOff){
					Continue=TRUE	
					CutOff=c()
					Flagged=c()
					while(Continue){
						M=mean(Temp)
						SD=stats::sd(Temp)	
						Dist=abs(Temp-M)
						Test=Temp[which.max(Dist)]
						if(Test>=M){
							P=stats::pnorm(Test,m=M,sd=SD,lower.tail=FALSE)
						}
						else{
							P=stats::pnorm(Test,m=M,sd=SD,lower.tail=TRUE)
						}
						if(round(P,2)>=0.05){
							CutOff=c(CutOff)
							Continue=FALSE
						}	
						else{
							CutOff=c(CutOff,Temp[which.max(Dist)])
							Temp=Temp[-c(which.max(Dist))]
							R=range(Temp)
							RR=R[2]-R[1]
							if(round(RR,2)<=RangeCutOff){
								Continue=FALSE								
							}		
						}
					}
					Temp=Temp[-c(which.max(Dist))]
					R=range(Temp)
					RR=R[2]-R[1]
					if(round(RR,2)>2*RangeCutOff&substr(e,1,3)=="JUC"){
						DelJ=c(DelJ,e)
					}
					
					if(ncol(Subset)<=6){
						AtLeast=ncol(Subset)
					}
					else{
						AtLeast=ncol(Subset)-1
					}
					
					if(length(CutOff)>0){
						F=c()
						P=c()
						for(c in CutOff){
							if(c>=M){
								P=c(P,stats::pnorm(c,m=M,sd=SD,lower.tail=FALSE))
							}
							else{
								P=c(P,stats::pnorm(c,m=M,sd=SD,lower.tail=TRUE))
							}
							
							F=c(F,rownames(which(Subset==c,arr.ind=TRUE))[1])
						}
						FlagP=c()
						names(P)=F
						for(f in unique(F)){
							FlagP=c(FlagP,mean(P[f]))
						}
						names(FlagP)=unique(F)
						Flagged=names(which(table(F)>=AtLeast))
						FlagP=FlagP[Flagged]
					}
					
					if(length(Flagged)>0){
						#print(e)
						#print(length(Flagged))
						if(length(Flagged)>length(StopDel)){
							Flagged=names(sort(FlagP)[1:StopDel])
							Rem=c(Rem,Flagged)
						}
						#GData=GData[-c(which(rownames(GData)%in%Flagged)),]
					}
				}
			}
			
			
		}
		
		if(length(Rem)>0){
			if(any(is.na(Rem))){
				Rem=Rem[-c(which(is.na(Rem)))]
			}
			if(length(Rem)>0){
				for(g in 1:length(Groups)){
					GData=DataS[,c(1,2,Groups[[g]]+2)]
					GData=GData[-c(which(rownames(GData)%in%Rem)),]
					GroupData[[g]]=GData	
				}
			}
			else{
				for(g in 1:length(Groups)){
					GData=DataS[,c(1,2,Groups[[g]]+2)]
					GroupData[[g]]=GData	
				}
			}
		}
		else{
			for(g in 1:length(Groups)){
				GData=DataS[,c(1,2,Groups[[g]]+2)]
				GroupData[[g]]=GData	
			}
		}
		
		print(geneID)	
		
		DelJuncs=names(which(table(DelJ)==length(Groups)))
		Jucs=unique(JAnnot[,3])[which(!is.na(unique(JAnnot[,3])))]
		
		
		DABGPSR=unique(PSRsExons)[which(unique(PSRsExons)%in%Low_AllSamples)]
		if(length(DABGPSR)>0){
			PSRsExons=PSRsExons[-c(which(PSRsExons%in%DABGPSR))]
		}
		
		Const=unique(PSRsExons)[which(!unique(PSRsExons)%in%ASPSR)]
		AS=unique(PSRsExons)[which(unique(PSRsExons)%in%ASPSR)]
		JAsses=data.frame(matrix(0,nrow=length(Jucs),ncol=(3+length(Groups))))
		colnames(JAsses)=c("Pattern","Flat","Low",paste("Low-",c(1:length(Groups)),sep=""))
		Hold=1
		Connections=data.frame(matrix(0,ncol=(2+length(Groups)),nrow=(nrow(JAsses)*length(Groups))))
		Place=1
		#ExonAssesment=data.frame(matrix(0,ncol=2,nrow=length(unique(JAnnot$PSR_ID))))
		#rownames(ExonAssesment)=unique(JAnnot$PSR_ID)
		ASJ=c()
		ExclDef=c()
		
		# Assesment of junctions
		if(length(Jucs)>0){
			for(k in 1:length(Jucs)){
				j=Jucs[k]
				
				if(is.na(j)){
					next
				}
				
				JValues=list()
				JList=list()
				JLengths=list()
				JAssesP=c()
				JAssesV=c()
				JFlat=c()
				JLow=c()
				
				#valid probes?
				if(nrow(DataS[which(DataS$ExonID==j),-c(1,2)])<=3){
					JAssesP=c(JAssesP,"Junction has too few valid probes",rep("-",(2+length(Groups))))
					JAsses[k,]=JAssesP
					Hold=Hold+1
					next				
				}
				
				#Low_AllSamples
				DABG=TRUE
				if(j%in%Low_AllSamples){
					JAssesP=c(JAssesP,"Junction is not DABG",c("-",TRUE,rep("-",length(Groups))))
					#DABG=FALSE	
				}
				
				if(j%in%DelJuncs){
					JAssesP=c(JAssesP,"Junction has a too large spread of values",rep("-",(2+length(Groups))))
					JAsses[k,]=JAssesP
					Hold=Hold+1
					next
				}
				
				if(j%in%unique(DataS$ExonID)){
					JD=list()
					for(g in 1:length(Groups)){
						JD[[g]]=as.vector(as.matrix(GroupData[[g]][which(GroupData[[g]]$ExonID==j),-c(1,2)]))
					}	
					#JUC_Ranks=sort(JD,index.return=TRUE)$ix
					JUC_Ranks=rank(unlist(JD),ties.method="random")
					JList[[length(JList)+1]]=JUC_Ranks
					names(JList)[length(JList)]=j
					JLengths=c(JLengths,length(JUC_Ranks))
					JValues=list(JD)
					names(JValues)=j
				}
				
				Set=sort(JAnnot[which(JAnnot[,3]==j&(JAnnot[,4]!="exclusion")),2])
				L=list()
				D=list()
				if(!all(Set%in%DataS$ExonID)){
					Set=Set[-c(which(!Set%in%DataS$ExonID))]
				}
				if(length(Set)==2&length(unique(Set))==1){
					print(j)
				}
				if(length(Set)>1){
					for(s in Set){	
						PSR=list()
						for(g in 1:length(Groups)){
							PSR[[g]]=as.vector(as.matrix(GroupData[[g]][which(GroupData[[g]]$ExonID==s),-c(1,2)]))
						}
						D[[s]]=PSR
						Ranks=list(rank(unlist(PSR),ties.method="random"))
						L=c(L,Ranks)
					}
					names(L)=Set	
				}
				else{
					JAssesP=c(JAssesP,"Junction does not have 2 anchor points",rep("-",(2+length(Groups))))
					JAsses[k,]=JAssesP
					Hold=Hold+1
					next			
				}
				L=c(L,JList)
				D=c(D,JValues)
				
				LL=sapply(unlist(D,recursive=FALSE),length)
				if(length(unique(LL))>1){
					MinLength=min(LL)
					Index=which(LL>MinLength)
					for(sj in 1:length(Index)){
						si=Index[sj]
						t=MinLength/(length(Groups[[1]]))
						Dnew=c()
						temp=c()
						g=as.numeric(substr(names(si),nchar(names(si)),nchar(names(si))))
						p=substr(names(si),1,nchar(names(si))-1)
						temp=GroupData[[g]][which(GroupData[[g]]$ExonID==p),-c(1,2)]
						M=apply(as.matrix(temp),2,stats::median)
						for(c1 in 1:ncol(temp)){
							Dist=abs(temp[[c1]]-M[c1])
							I=rank(Dist)
							Select=which(I<=t)
							if(length(Select)<t){
								Add=t-length(Select)
								SelectAdd=which((I>t)&(I<=(t+Add)))[c(1:Add)]
								Select=c(Select,SelectAdd)
							}
							else if(length(Select)>t){
								Del=length(Select)-t
								SelectDel=which(Select>=t)[c(1:Del)]
								Select=Select[-c(SelectDel)]
							}
							NewC=temp[Select,c1]
							Dnew=cbind(Dnew,NewC)
						}
						D[[p]][[g]]=as.vector(as.matrix(Dnew))
						PSR=unlist(D[[p]])
						#Ranks=list(sort(PSR,index.return=TRUE)$ix)
						Ranks=rank(PSR)
						L[[p]]=Ranks	
						
					}
				}
				
				if(length(Set)>1){
					L[[1]]=(L[[1]]+L[[2]])/2
					L=L[-c(2)]
					L[[1]]=sort(L[[1]],index.return=TRUE)$ix
				}	
				Ranks=do.call("cbind",L)
				
				Y=as.vector(as.matrix(Ranks))
				exon=c()
				for(c in 1:ncol(Ranks)){
					exon=c(exon,rep(c,nrow(Ranks)))
				}
				tissuetemp=c()
				for(g in 1:length(Groups)){
					tissuetemp=c(tissuetemp,rep(g,nrow(Ranks)/(length(Groups))))
				}
				tissue=rep(tissuetemp,ncol(Ranks))
				
				exon=as.factor(exon)
				tissue=as.factor(tissue)
				ft1<-stats::lm(Y~exon*tissue-1)
				
				Inter=summary(ft1)$coefficients[((length(L)+2):nrow(summary(ft1)$coefficients)),4]
				
				if(all(round(Inter,2)>=0.05)){
					# Junction is product of its supporting exons
					JAssesP=c(JAssesP,"Pattern Supported")						
				}
				else{
#					if(DABG==FALSE){
#						JAssesP=c(JAssesP,"Junction is not DABG",c("-",TRUE,rep("-",length(Groups))))
#						JAsses[k,]=JAssesP
#						Hold=Hold+1
#						next
#					}
					JAssesP=c(JAssesP,"Pattern Not Supported")
				}
				
				#Junction Value Assessment
				Set=sort(JAnnot[which(JAnnot[,3]==j&(JAnnot[,4]!="exclusion")),2])
#				
				#Flat line
				tissuetemp=c()
				for(g in 1:length(Groups)){
					tissuetemp=c(tissuetemp,rep(g,length(D[[1]][[1]])))
				}
				Flattest=stats::lm(unlist(D[[length(D)]])~tissuetemp-1)
				Flat=summary(Flattest)$coefficients[1,4]
				if(round(Flat,2)>0.05){
					JFlat=c(JFlat,TRUE)
				}
				else{
					JFlat=c(JFlat,FALSE)
					
				}	
				
				#Low??
				Low=rep(FALSE,length(Groups))
				JLow=FALSE
				if(j%in%unlist(Low_GSamples)){
					R=lapply(Low_GSamples,function(x) j%in%x)
					R=unlist(R)
					if(all(R)){
						if(JFlat){
							Low[which(R)]=TRUE
						}	
					}
					else if(any(R)){
						Low[which(R)]=TRUE
					}
				}					
				Row=c(JAssesP,JFlat,JLow,Low)
				JAsses[k,]=Row
				
			}
			JAsses=cbind(Jucs,JAsses)
			
			
			# PSR reflection on junctions
			
			for(r in 1:nrow(JAsses)){
				R=JAsses[r,]
				J=as.character(JAsses[r,1])
				PSRs=unique(JAnnot[,2][which(JAnnot[,3]==J)])
				supp=sort(JAnnot[,2][which(JAnnot[,2]%in%PSRs&JAnnot[,3]==J&JAnnot[,4]!="exclusion")])
				if(!all(supp%in%DataS$ExonID)){
					supp=supp[-c(which(!supp%in%DataS$ExonID))]
				}
				excl=sort(JAnnot[,2][which(JAnnot[,2]%in%PSRs&JAnnot[,3]==J&JAnnot[,4]=="exclusion")])
				if(!all(excl%in%DataS$ExonID)){
					excl=excl[-c(which(!excl%in%DataS$ExonID))]
				}
				if(R[2]%in%c("Junction has too few valid probes","Junction has a too large spread of values","Junction does not have 2 anchor points")){
					next
				}
				#Do the support points connect? Check Low, 1-Low and 2-Low
				for(g in 1:length(Groups)){
					
					GData=GroupData[[g]]
					
					if(R[2]=="Junction is not DABG"){
						Connections[Place,c(1,2,g+2)]=c(J,paste(supp,collapse="-"),"never")
						Place=Place+1
					}
					
					else if(as.logical(R[5+g-1])){
						Connections[Place,c(1,2,g+2)]=c(J,paste(supp,collapse="-"),"never")
						Place=Place+1
						#junction is low	
					}	
					else{
						Connections[Place,c(1,2,g+2)]=c(J,paste(supp,collapse="-"),"present")
						Place=Place+1
						#junction is present ==> both linking points are present and connection is present
						#if excl junction: exon is indeed excluded => AS
						
						if(length(excl)>0){					
							for(e in excl){
								ExclDef=c(ExclDef,e)
								
								if(JAsses[r,2]=="Pattern Not Supported"){
									
									JAsses[r,2]="Pattern Supported"
								}		
							}
						}		
					}	
				}	
			}		
		}
		
		
		if(any(Connections[,1]==0)){
			Connections=Connections[-c(which(Connections[,1]==0)),]
		}
		if(any(duplicated(Connections))){
			Connections=Connections[!duplicated(Connections),]
		}
		Links=c()		
		UL=unique(Connections[,c(1,2)])
		for(c in 1:nrow(UL)){
			rowLinks=c(cbind(UL[c,1],UL[c,2]))
			Set=which(Connections[,1]==UL[c,1]&Connections[,2]==UL[c,2])
			for(g in 1:length(Groups)){
				GAsses=Connections[Set,g+2]
				Get=GAsses[which(GAsses!=0)]
				if(all(Get=="never")){
					rowLinks=c(rowLinks,"never")
				}
				else{
					rowLinks=c(rowLinks,"present")
				}
			}
			Links=rbind(Links,rowLinks)
			
		}	
		
		Edges=c()
		#Exons first
		Exons=sort(unique(c(unique(JAnnot[,2]),c(unique(PSRsExons),unique(DABGPSR)))))
		Exons=Exons[which(Exons%in%DataS[,2])]
		if(length(Exons)==0){
			output=paste(geneID,"No PSR's present",sep="\t")
			if(!is.null(Location)){
				utils::write.table(output, paste(Location,"/",Name,"_NotAssessedTranscripts.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				#save(output, file=paste(Location,"/REIDS_IsoformInfo_", geneID, ".RData", sep=""))
				rm(output)
				gc()
				close(conn)
				return(paste(geneID," Completed",sep=""))
			}
			else{
				close(conn)
				return(output)
			}
		}
		for(i in 1:(length(Exons)-1)){
			R=c(Exons[i],Exons[i+1])
			Edges=rbind(Edges,R)
		}
		col1=rep("grey",nrow(Edges))
		width1=rep(2,nrow(Edges))
		G=list()
		Cols=list()
		Widths=list()
		for(g in 1:(length(Groups)+1)){
			G[[g]]=Edges
			Cols[[g]]=col1	
			Widths[[g]]=width1
		}
		#Links
		if(length(Jucs)>0&length(Links)>0){
			for(l in 1:nrow(Links)){
				E=strsplit(Links[l,2],"-")[[1]]
				if(length(E)>1){
					for(g in 1:(length(Groups)+1)){
						G[[g]]=rbind(G[[g]],E)
						if(g==1){
							Cols[[1]]=c(Cols[[1]],"blue")
							Widths[[1]]=c(Widths[[1]],2)
						}
						else{
							L=Links[l,g+1]
							if(L=="present"){
								Cols[[g]]=c(Cols[[g]],"blue")
								Widths[[g]]=c(Widths[[g]],2)
							}
							else{
								Cols[[g]]=c(Cols[[g]],"red")
								Widths[[g]]=c(Widths[[g]],2)
							}
						}
					}			
				}
			}
			#rownames(Edges)=1:nrow(Edges)
		}
		
		if(length(Groups)==2){
			FC=c()
			for(e in Exons){				
				A=mean(as.vector(as.matrix(GroupData[[1]][which(GroupData[[1]]$ExonID==e),-c(1,2)])))
				B=mean(as.vector(as.matrix(GroupData[[2]][which(GroupData[[2]]$ExonID==e),-c(1,2)])))
				FC=c(FC,A-B)
			}	
			
			names(FC)=Exons
			
			
			if(length(Const)==0){
				MedFC=stats::median(FC)
				SDFC=0.5
			}
			else{
				ConstFC=FC[Const]
				MedFC=stats::median(ConstFC)
				SDFC=stats::sd(ConstFC)
				if(is.na(SDFC)){
					SDFC=0.5
				}
			}
			ASFC=FC[AS]
			ASPvals=c()
			for(a in FC){
				if(a>MedFC){
					ASPvals=c(ASPvals,stats::pnorm(a,MedFC,SDFC,lower.tail=FALSE))
				}
				else{
					ASPvals=c(ASPvals,stats::pnorm(a,MedFC,SDFC,lower.tail=TRUE))
				}	
			}
			names(ASPvals)=names(FC)
			ASPvals=stats::p.adjust(ASPvals,"fdr")
			ASfctemp=names(which(ASPvals<0.05))
			ASfc=ASfctemp[which(!ASfctemp%in%c(AS,ASJ))]
		}
		else{
			FC=rep("-",length(Exons))
			ASfc=NULL
		}
		colNodes=rep("grey",length(Exons))
		names(colNodes)=Exons
		ColN=list()
		for(g in 1:(length(Groups)+1)){
			ColN[[g]]=colNodes			
		}
		for(c in Exons){
			if(c%in%DABGPSR){
				for(g in c(1:length(Groups)+1)){
					ColN[[g]][c]="grey"
				}	
			}
			else if(c%in%Const){
				for(g in c(1:length(Groups)+1)){
					ColN[[g]][c]="black"
				}	
			}
			else{
				if(length(Groups)==2){
					#if(round(ASPvals[c],2)<0.05){
					if(FC[c]<0){
						ColN[[2]][c]="red"
						ColN[[3]][c]="green"
					}
					else{
						ColN[[2]][c]="green"
						ColN[[3]][c]="red"
					}
					#}
				}
			}	
		}
		
		
		DelJ=unique(c(DelJ,as.character(JAsses[,1][which(JAsses[,2]%in%c("Junction has too few valid probes","Junction has a too large spread of values","Junction does not have 2 anchor points"))])))
		OutList=data.frame("TC_ID"=rep(DataS[1,1],length(Exons)),"PSR_ID"=Exons,"Type"=rep(0,length(Exons)),"Unreliable Junctions"=rep(0,length(Exons)),"Linking Junctions"=rep(0,length(Exons)),"Exclusion Junctions"=rep(0,length(Exons)),"Supported by"=rep(0,length(Exons))
				,"Identified by"=rep(0,length(Exons)),"Fold Change"=rep(0,length(Exons)),"Exons"=rep(0,length(Exons)))
		
		#Reduce transcripts in all samples:
		NeverPresent=which(apply(Links[,3:ncol(Links),drop=FALSE],1,function(x) all(x==rep("never",length(Groups)))))
		#NeverLinkedPSR=list()
		NeverLinkedTr=list()
		if(length(NeverPresent)>0){
			RLinks=Links[NeverPresent,,drop=FALSE]		
			for(r in 1:nrow(RLinks)){
				PSR=RLinks[r,2]
				PSRs=strsplit(PSR,"-")[[1]]
				set=c()
				for(t in unique(TrAnnot[,5])){
					PSRset=TrAnnot[which(TrAnnot[,5]==t),2]
					if(all(PSRs%in%PSRset)){
						set=c(set,t)
					}
				}
				NeverLinkedTr[[r]]=set
#			set=TrAnnot[which(TrAnnot[,2]==RLinks[r,1]),5]
#			#set=strsplit(RLinks[r,2],"-")[[1]]		
#			#NeverLinkedPSR[[r]]=set
#			NeverLinkedTr[[r]]=set
			}
			NeverLinkedTr=unique(unlist(NeverLinkedTr))
		}
		
		
		#Event Type
		TC=DataS[1,1]
		Strand=TrAnnot[1,3]
		TRS=list()
		N=c()
		for(tr in unique(TrAnnot[,5])){
			if(tr%in%NeverLinkedTr){
				next
			}
			Set=TrAnnot[which(TrAnnot[,5]==tr),2]
			PSet=Set[which(substr(Set,1,1)=="P")]
			ESete=TrAnnot[which(TrAnnot[,5]==tr),4]
			ESet=ESete[which(substr(Set,1,1)=="P")]
#			Remove=FALSE
#			if(length(NeverLinkedPSR)>0){
#				for(nl in NeverLinkedPSR){
#					Pos1=which(PSet==nl[1])
#					Pos2=which(PSet==nl[2])
#					if(length(Pos1)>0&length(Pos2)>0){
#						if(Pos1+1==Pos2){
#							#linked in transcript but the link is not present thus transcript can be removed from the collection
#							Remove=TRUE
#						}
#					}
#				}
#			}
			if(length(DABGPSR)>0){
				Stop=FALSE
				for(d in DABGPSR){
					if(any(PSet==d)){
						Stop=TRUE
						break					
					}
				}
				if(Stop){
					next
				}
			}
			
			if(Strand=="-"|Strand==-1){
				PSet=rev(sort(PSet))
				ESet=rev(sort(ESet))
			}
			names(ESet)=PSet
			TRS[[length(TRS)+1]]=ESet	
			N=c(N,tr)
		}
		names(TRS)=N
		if(length(TRS)==0){
			
			output=list("Output"=OutList)
			if(!is.null(Location)){
				if(!is.null(output[[1]])){
					utils::write.table(output[[1]], paste(Location,"/",Name,"_ASInfo.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
				}
				#save(output, file=paste(Location,"/REIDS_IsoformInfo_", geneID, ".RData", sep=""))
				rm(output)
				gc()
				close(conn)
				return(paste(geneID," Completed",sep=""))	
			}	
			else{
				close(conn)
				return(output)
			}	
		}
		
		AllExons=unique(unlist(TRS))		
		for(e in Exons){
			EAnnotp=TrAnnot[which(TrAnnot[,2]==e),4]
			EAnnotp=EAnnotp[which(EAnnotp%in%AllExons)]
			EAnnotp=paste(EAnnotp,collapse="|")
			s=JAnnot[which(JAnnot[,2]==e),c(3,4),drop=FALSE]
			Dels=s[which(s[,1]%in%DelJ),1]
			if(length(Dels)>0){
				s=s[which(!s[,1]%in%DelJ),,drop=FALSE]
			}
			else{
				Dels="-"
			}
			if(e%in%DABGPSR){
				r=c("not DABG",rep("-",6),EAnnotp)
			}
			
			else if(e%in%Const){			
				if(!(all(is.na(s)))){
					Jss=c()
					Jsse=c()
					Js=s[which(s[,2]%in%c(3,5)),1]
					if(length(Js)>0){
						Jss=as.character(JAsses[which(JAsses[,1]%in%Js&JAsses[,2]=="Pattern Supported"&JAsses[,4]==FALSE),1])
						if(!(all(Js%in%Jss))){
							Dels=c(Dels,Js[which(!Js%in%Jss)])
							if(length(Jss)==0){
								Jss="-"
							}
						}
					}
					else{
						Jss="-"
					}
					Jse=s[which(s[,2]=="exclusion"),1]
					if(length(Jse)>0){
						Jsse=as.character(JAsses[which(JAsses[,1]%in%Jse&JAsses[,4]==FALSE),1])
					}
					else{
						Jsse="-"
					}
					SupportedBy=c(Jss,Jsse)
					if(all(SupportedBy=="-")){
						SupportedBy="-"
					}
					else if(all(SupportedBy%in%s[,1])){
						SupportedBy="All"
					}	
					else if(!any(SupportedBy%in%s[,1])){
						SupportedBy="None"
					}
					else{
						if(any(SupportedBy=="-")){
							SupportedBy=SupportedBy[-c(which(SupportedBy=="-"))]
						}
					}
					if(length(Dels)>1&any(Dels=="-")){
						Dels=Dels[-c(which(Dels=="-"))]
					}
					
					r=c("Const",paste(Dels,collapse="|"),paste(Jss,collapse="|"),paste(Jsse,collapse="|"),paste(SupportedBy,collapse="|"),"-",FC[e],EAnnotp)
				}	
				else{
					r=c("Const",rep("-",5),FC[e],EAnnotp)
				}
			}
			else if(e%in%AS){
				
				if(!(all(is.na(s)))){
					Jss=c()
					Jsse=c()
					Js=s[which(s[,2]%in%c(3,5)),1]
					if(length(Js)>0){
						Jss=JAsses[which(JAsses[,1]%in%Js&JAsses[,2]=="Pattern Supported"&JAsses[,4]==FALSE),1]
						if(!(all(Js%in%Jss))){
							Dels=c(Dels,Js[which(!Js%in%Jss)])
							if(length(Jss)==0){
								Jss="-"
							}
						}
					}
					else{
						Jss="-"
					}
					Jse=s[which(s[,2]=="exclusion"),1]
					if(length(Jse)>0){
						Jsse=JAsses[which(JAsses[,1]%in%Jse&JAsses[,4]==FALSE),1]
					}
					else{
						Jsse="-"
					}
					SupportedBy=c(as.character(Jss),as.character(Jsse))
					if(all(SupportedBy=="-")){
						SupportedBy="-"
					}
					else if(all(SupportedBy%in%s[,1])){
						SupportedBy="All"
					}	
					else if(!any(SupportedBy%in%s[,1])){
						SupportedBy="None"
					}	
					else{
						if(any(SupportedBy=="-")){
							SupportedBy=SupportedBy[-c(which(SupportedBy=="-"))]
						}
					}
					if(length(Dels)>1&any(Dels=="-")){
						Dels=Dels[-c(which(Dels=="-"))]
					}
					r=c("AS",paste(Dels,collapse="|"),paste(Jss,collapse="|"),paste(Jsse,collapse="|"),paste(SupportedBy,collapse="|"),"REIDS",FC[e],EAnnotp)
				}	
				else{
					r=c("AS",rep("-",4),"REIDS",FC[e],EAnnotp)
				}
			}
			else if(e%in%ASJ){
				if(!(all(is.na(s)))){
					Jss=c()
					Jsse=c()
					Js=s[which(s[,2]%in%c(3,5)),1]
					if(length(Js)>0){
						Jss=JAsses[which(JAsses[,1]%in%Js&JAsses[,2]=="Pattern Supported"&JAsses[,4]==FALSE),1]
						if(!(all(Js%in%Jss))){
							Dels=c(Dels,Js[which(!Js%in%Jss)])
							if(length(Jss)==0){
								Jss="-"
							}
						}
					}
					else{
						Jss="-"
					}
					Jse=s[which(s[,2]=="exclusion"),1]
					if(length(Jse)>0){
						Jsse=JAsses[which(JAsses[,1]%in%Jse&JAsses[,4]==FALSE),1]
					}
					else{
						Jsse="-"
					}
					SupportedBy=c(as.character(Jss),as.character(Jsse))
					if(all(SupportedBy=="-")){
						SupportedBy="-"
					}
					else if(all(SupportedBy%in%s[,1])){
						SupportedBy="All"
					}	
					else if(!any(SupportedBy%in%s[,1])){
						SupportedBy="None"
					}
					else{
						if(any(SupportedBy=="-")){
							SupportedBy=SupportedBy[-c(which(SupportedBy=="-"))]
						}
					}
					if(length(Dels)>1&any(Dels=="-")){
						Dels=Dels[-c(which(Dels=="-"))]
					}
					r=c("AS",paste(Dels,collapse="|"),paste(Jss,collapse="|"),paste(Jsse,collapse="|"),paste(SupportedBy,collapse="|"),"Junction Support",FC[e],EAnnotp)
				}
				else{
					r=c("AS",rep("-",4),"Junction Support",FC[e],EAnnotp)
				}	
			}
			else if(e%in%ASfc){
				
				if(!(all(is.na(s)))){
					Jss=c()
					Jsse=c()
					Js=s[which(s[,2]%in%c(3,5)),1]
					if(length(Js)>0){
						Jss=JAsses[which(JAsses[,1]%in%Js&JAsses[,2]=="Pattern Supported"&JAsses[,4]==FALSE),1]
						if(!(all(Js%in%Jss))){
							Dels=c(Dels,Js[which(!Js%in%Jss)])
							if(length(Jss)==0){
								Jss="-"
							}
						}
					}
					else{
						Jss="-"
					}
					Jse=s[which(s[,2]=="exclusion"),1]
					if(length(Jse)>0){
						Jsse=JAsses[which(JAsses[,1]%in%Jse&JAsses[,4]==FALSE),1]
					}else{
						Jsse="-"
					}
					SupportedBy=c(as.character(Jss),as.character(Jsse))
					if(all(SupportedBy=="-")){
						SupportedBy="-"
					}
					else if(all(SupportedBy%in%s[,1])){
						SupportedBy="All"
					}	
					else if(!any(SupportedBy%in%s[,1])){
						SupportedBy="None"
					}	
					else{
						if(any(SupportedBy=="-")){
							SupportedBy=SupportedBy[-c(which(SupportedBy=="-"))]
						}
					}
					if(length(Dels)>1&any(Dels=="-")){
						Dels=Dels[-c(which(Dels=="-"))]
					}
					r=c("AS",paste(Dels,collapse="|"),paste(Jss,collapse="|"),paste(Jsse,collapse="|"),paste(SupportedBy,collapse="|"),"Fold Change",FC[e],EAnnotp)
				}	
				else{
					r=c("AS",rep("-",4),"Fold Change",FC[e],EAnnotp)
				}
				
			}
			else{
#				if(e%in%DataS[,2]&(!e%in%JAnnot[,2])){
#					r=c("INI Filtered",rep("-",6))
#				}
				#else{
				print(e)
				#}
			}
			OutList[OutList[,2]==e,c(3:10)]=r
		}
		
		L=matrix(0,ncol=3,nrow=nrow(Links))
		for(i in 1:nrow(Links)){
			r=Links[i,c(2)]
			PSRs=strsplit(r,"-")[[1]]
			L[i,]=c(Links[i,1],PSRs[1],PSRs[2])
		}
		
		Table1=matrix(0,nrow=length(TRS)*2,ncol=4)
		l=c()
		for(i in 1:length(TRS)){
			N=names(TRS)[i]	
			l=c(l,L[which(L[,2]%in%names(TRS[[i]])&L[,3]%in%names(TRS[[i]])),1])
			J=paste(L[which(L[,2]%in%names(TRS[[i]])&L[,3]%in%names(TRS[[i]])),1],collapse="|")
			R1=paste(names(TRS[[i]]),collapse="|")
			R2=paste(TRS[[i]],collapse="|")
			Table1[i*2-1,]=c(TC,N,J,R1)
			Table1[i*2,]=c(TC,N,J,R2)
		}
		l=unique(l)
		if(any(!Links[,1]%in%l)){
			NovelConnections=Links[which(!Links[,1]%in%l),1]
			Novel=cbind(NovelConnections,Links[which(Links[,1]%in%NovelConnections),])
			NovelConnections=cbind(TC,NovelConnections)
		}	
		else{
			NovelConnections=c()
		}
		
		
		#utils::write.table(Table1, "Transcripts.txt", sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		#Transcripts per Groups
		
		GroupTranscripts=list()
		for(g in 1:length(Groups)){
			N=c()
			GroupTranscripts[[g]]=list()
			NeverPresent=which(sapply(Links[,g+2],function(x) all(x=="never")))
			RLinks=Links[NeverPresent,,drop=FALSE]
			#NeverLinkedPSR=list()
			NeverLinkedTr=list()
			if(nrow(RLinks)>0){
				for(r in 1:nrow(RLinks)){
					set=c()
					for(t in unique(TrAnnot[,5])){
						PSRset=TrAnnot[which(TrAnnot[,5]==t),2]
						if(all(PSRs%in%PSRset)){
							set=c(set,t)
						}
					}
					#set=Transcripts[which(Transcripts[,2]==RLinks[r,1]),5]
					NeverLinkedTr[[r]]=set
					#set=strsplit(RLinks[r,1],"-")[[1]]
					#NeverLinkedPSR[[r]]=set
				}
				NeverLinkedTr=unique(unlist(NeverLinkedTr))
			}	
			for(tr in unique(TrAnnot[,5])){
				if(tr%in%NeverLinkedTr){
					next
				}
				Set=TrAnnot[which(TrAnnot[,5]==tr),2]
				PSet=Set[which(substr(Set,1,1)=="P")]
				ESete=TrAnnot[which(TrAnnot[,5]==tr),4]
				ESet=ESete[which(substr(Set,1,1)=="P")]
				
				
				if(length(Low_GSamples[[g]])>0){
					if(any(PSet%in%Low_GSamples[[g]])){
						next
					}
				}
				
				if(Strand=="-"|Strand==-1){
					PSet=rev(sort(PSet))
					ESet=rev(sort(ESet))
				}
				
				GroupTranscripts[[g]][[length(GroupTranscripts[[g]])+1]]=ESet	
				N=c(N,tr)
			}
			names(GroupTranscripts[[g]])=N
			
		}
		
		Table2=matrix(0,nrow=length(TRS),ncol=(2+length(Groups)))
		for(i in 1:length(TRS)){
			N=names(TRS)[i]
			Present=rep("no",length(Groups))
			for(g in 1:length(Groups)){
				if(N%in%names(GroupTranscripts[[g]])){
					Present[g]="Yes"
				}
			}
			Table2[i,]=c(TC,N,Present)
		}
		#utils::write.table(Table2, "Transcripts_Groups.txt", sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		
		
		if(Strand=="-"|Strand==-1){
			OutList=OutList[nrow(OutList):1,]
		}
		
		Decision=OutList[,3]
		names(Decision)=OutList[,2]
		C=names(Decision)[which(Decision=="Const")]
		
		
		Category=rep("",nrow(OutList))
		
		for(r in 1:nrow(OutList)){
			if(Decision[r]=="Const"){
				Category[r]="-"
				next
			}
			else if(Decision[r]=="not DABG"){
				Category[r]="-"
				next
			}
			PSR=as.character(OutList[r,2])
			EAnnotp=TrAnnot[which(TrAnnot[,2]==PSR),,drop=FALSE]
			if(nrow(EAnnotp)==0){
				Category[r]="Intron Retention"	
				next
			}	
			else{
				Temp=c()
				for(ea in 1:nrow(EAnnotp)){
					OtherPSR=as.character(TrAnnot[which(TrAnnot[,4]==as.character(EAnnotp[ea,4])),2])
					OtherPSR=OtherPSR[which(OtherPSR!=PSR)]
					if(length(OtherPSR)>0){					
						S=substr(OtherPSR,1,3)
						Temp=c(Temp,OtherPSR[which(S=="PSR")])
					}
					else{
						Temp=c(Temp,"")
					}
					
				}
				if(any(Temp!="")){
					#exon has more than 1 probe set annotated to itself
					Temp=unique(Temp)
					if(any(Temp=="")){
						Temp=Temp[-c(which(Temp==""))]
					}
					if(all(which(OutList[,2]%in%Temp)>r)&r==1){
						Category[r]="Alternative First'"
						next
					}
					else if(all(which(OutList[,2]%in%Temp)>r)){
						Category[r]="Alternative 5'"
						for(tra in 1:length(TRS)){
							if(names(TRS[[tra]])[1]==PSR){
								Category[r]="Alternative First"
								break
							}
						}
						next
					}
					else if(all(which(OutList[,2]%in%Temp)<r)&r==nrow(OutList)){
						Category[r]="Alternative Last"
						next
					}
					else if(all(which(OutList[,2]%in%Temp)<r)){
						Category[r]="Alternative 3'"
						for(tra in c(1:length(TRS))){
							if(names(TRS[[tra]])[length(TRS[[tra]])]==PSR){
								Category[r]="Alternative Last"
								break
							}
						}
						next
					}
					else{
						Category[r]="Complex Event"
						next
					}
					
				}
				else{
					NeighboursIn=c()
					NeighboursOut=c()
					Positions=c()
					for(tra in 1:length(TRS)){
						if(length(TRS[[tra]])==0){
							next
						}
						
						if(names(TRS[[tra]])[1]==PSR){
							Category[r]="Alternative First"
							break
						}
						else if(names(TRS[[tra]])[length(TRS[[tra]])]==PSR){
							Category[r]="Alternative Last"
							break
						}
						else{
							Pos=which(names(TRS[[tra]])==PSR)
							if(length(Pos)>0){
								Left=names(TRS[[tra]])[Pos-1]
								Right=names(TRS[[tra]])[Pos+1]
								NeighboursIn=c(NeighboursIn,Left,Right)	
							}
							else if(names(TRS[[tra]])[1]%in%OutList[,2]&names(TRS[[tra]])[length(TRS[[tra]])]%in%OutList[,2]){
								if(which(OutList[,2]==names(TRS[[tra]])[1])<r&which(OutList[,2]==names(TRS[[tra]])[length(TRS[[tra]])])>r){				
									F=c(names(TRS[[tra]]),PSR)
									F=sort(F)
									Pos=which(F==PSR)
									Left=F[Pos-1]
									Right=F[Pos+1]
									NeighboursOut=c(NeighboursOut,Left,Right)	
								}
							}
							
						}
						
					}
					
					if(Category[r]==""){
						if(length(NeighboursIn)>0&length(NeighboursOut)>0){
							NIN=sort(unique(NeighboursIn))
							NOUT=sort(unique(NeighboursOut))
							
							if(length(NIN)==2&length(NOUT)==2){
								if(all(NIN%in%NOUT)){
									if(Decision[OutList[which(OutList[,2]==NIN[1]),2]]=="Const"&Decision[OutList[which(OutList[,2]==NIN[2]),2]]=="Const"){
										Category[r]="Cassette Exon"
									}
									else{
										Category[r]="Complex Event"
									}
								}
								
								else if(NIN[1]==NOUT[1]&NOUT[2]==OutList[r+1,2]&NIN[2]!=OutList[r+1,2]){
									if(r<(length(Decision)-2)){
										if(Decision[OutList[r+2,2]]=="Const"){
											Category[r]="Mutually Exclusive Exon"
										}
										else{
											Category[r]="Complex Event"
										}
									}
									else{
										Category[r]="Complex Event"
									}
								}
								
								else if(NIN[2]==NOUT[2]&NOUT[1]==OutList[r-1,2]&NIN[1]!=OutList[r-1,2]){
									if(r>2){
										if(Decision[OutList[r-2,2]]=="Const"){
											Category[r]="Mutually Exclusive Exon"
										}	
										else{
											Category[r]="Complex Event"
										}
									}
									else{
										Category[r]="Complex Event"
									}
									
									
								}
								else{
									Category[r]="Complex Event"
								}
							}
							else{
								Category[r]="Complex Event"
							}
						}
						else if(length(NeighboursIn)>0){
							NIN=sort(unique(NeighboursIn))
							if(length(NIN)==2){
								if(NIN[1]%in%OutList[,2]&NIN[2]%in%OutList[,2]){
									if(Decision[OutList[which(OutList[,2]==NIN[1]),2]]=="Const"&Decision[OutList[which(OutList[,2]==NIN[2]),2]]=="Const"){
										Category[r]="Cassette Exon"
									}
									else{
										Category[r]="Complex Event"
									}
								}
								else{
									Category[r]="Complex Event"
								}
							}
							else{
								Category[r]="Complex Event"
							}
						}	
						else{
							Category[r]="Complex Event"
						}
					}
					else{
						next
					}
				}
			}
		}
		
		if(Plot){		
			#pdf("Figures/Tr_TC17001298.pdf",width=18,height=10)
			graphics::par(mfrow=c(2,1))
			graphics::par(mar=c(12,1,1,1))
			for(g in 2:3){
				if(Strand=="-"|Strand==-1){
					R=G[[g]][which(rownames(G[[g]])=="R"),]
					R=R[nrow(R):1,c(2,1)]
					E=G[[g]][which(rownames(G[[g]])=="E"),]
					E=E[nrow(E):1,c(2,1)]
					G1=rbind(R,E)
					
					C1=rev(Cols[[g]][which(rownames(G[[g]])=="R")])
					C2=rev(Cols[[g]][which(rownames(G[[g]])=="E")])
					C=c(C1,C2)
					
					arcplot(G1,lwd.arcs=rev(Widths[[g]])+1,above=NULL,col.arcs=C,ljoin = 2,lend=2,pch.nodes=15,cex.nodes=2,cex.labels=1.2,col.nodes=rev(ColN[[g]]))
				}
				else{
					arcplot(G[[g]],lwd.arcs=Widths[[g]],above=NULL,col.arcs=Cols[[g]],ljoin = 2,lend=2,pch.nodes=15,cex.nodes=2,cex.labels=1.5,col.nodes=ColN[[g]])
				}
			}
			#dev.off()
		}
		#warnings()
		OutList=cbind(OutList,Category)
		#utils::write.table(OutList, paste(Name,".txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		#rm(OutList)
		#gc() 
		output=list("Output"=OutList,"Compositions"=Table1,"GroupTranscripts"=Table2,"NovelConnections"=NovelConnections)
		if(!is.null(Location)){
	
			#save(output, file=paste(Location,"/REIDS_IsoformInfo_", geneID, ".RData", sep=""))
			if(!is.null(output[[1]])){
				utils::write.table(output[[1]], paste(Location,"/",Name,"_ASInfo.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
			}
			if(!is.null(output[[2]])){
				utils::write.table(output[[2]], paste(Location,"/",Name,"_Compositions.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)		
			
			}
			if(!is.null(output[[3]])){
				utils::write.table(output[[3]], paste(Location,"/",Name,"_GroupTranscripts.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
			}
			if(!is.null(output[[4]])){
				utils::write.table(output[[4]], paste(Location,"/",Name,"_NovelConnections.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)	
			}
			rm(output)
			gc()
			close(conn)
			return(paste(geneID," Completed"))
		}
		else{
			close(conn)
			return(output)
		}

	}	
	close(conn)
}
	

#' REIDSJunctionAssessment_HPCVersion
#' 
#' REIDSJunctionAssessment_HPCVersion is the HPC version of REIDSJunctionAssessment. This function should be used with the REIDSJunctionAssessment_HPCVersion.R file and REIDSJunctionAssessment_HPCVersion.pbs script in the documentation folder of the package.
#' After running this function on the cluster, the output files should be binded together with the CreateOutput function.
#' @export
#' @param geneID The gene ID.
#' @param DataS The data with as rows the probesets and as columns the samples. Note that the first column should contain the gene IDs and the second column the exon IDs
#' @param ASPSR The AS probe sets as identified by ASExons.
#' @param Juninfo A parameter specifying wether the annotations are user of Ensembl defined. If JunInfo is "User" (default) the annotations provided in EandTrAnnot are used. If JunInfo is "Ensembl" the annotations in EandTrAnnot are used to set up tje junction associations but the gene name and position in transcriptData and positionData are used to connect with the Ensembl data base and retrieve corresponding information. 
#' @param JAnnotI The file name with line indices for the junction associations.
#' @param JAnnot The file name with the junction associations.
#' @param EandTrAnnotI The file name with line indices for the exon and isoform annotations.
#' @param EandTrAnnot The file name with the exon and isoform annotations.
#' @param PartiallyAnnotated Logical. Should the exon annotations with partially annotated probe sets still be included? If FALSE, these are excluded. If TRUE, these are included. Default is FALSE.
#' @param positionData The file with the chromosome start and ends for the probe sets. Only needed in JunInfo=Ensembl.
#' @param transcriptData The file with gene name of the transcripts. Only needed in JunInfo=Ensembl.
#' @param Groups A list with  elements speficifing the columns of the data in each group.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @param Low_GSamples A list with a  character vector per group containing the probe sets which are not DABG in that group.
#' @param Plot Should a plot of the gene model be made?
#' @details The plot is produced by the arcplot function of the arcdiagram package (https://github.com/gastonstat/arcdiagram)
#' @return A .RData file will be saved for each gene with the four elements returned REIDSJunctionAssessment function. The outputs can be bound together by CreateOutput.
REIDSJunctionAssessment_HPCVersion<-function(geneID,DataS=DataS,ASPSR=ASPSR,Juninfo="User",JAnnotI,JAnnot=NULL,EandTrAnnotI=NULL,EandTrAnnot=NULL,PartiallyAnnotated,positionData=NULL,transcriptData=NULL,Groups=list(),Low_AllSamples=c(),Low_GSamples=c(),Plot=FALSE){

	if(!is.null(EandTrAnnotI)){
		ETrI=utils::read.table(EandTrAnnotI,header=FALSE,stringsAsFactors=FALSE)
		Lines=ETrI[which(ETrI[,1]==geneID),]
		ETrAnnot=utils::read.table(EandTrAnnot,header=FALSE,sep="\t",nrows=as.numeric(Lines[3]),skip=as.numeric(Lines[2]),stringsAsFactors=FALSE)
		TrAnnot=ETrAnnot[,c(1,2,8,4,7)]
		TrAnnot=TrAnnot[order(TrAnnot[,2]),]
		colnames(TrAnnot)=c("TC_ID","PSR_ID","strand","EAnnot","TrAnnot")
		
		if(!PartiallyAnnotated){
			Incl=unlist(sapply(TrAnnot[,4],function(x) !grepl("[*]",x)))
			TrAnnot=TrAnnot[Incl,]
		}
		
		
		if(!is.null(JAnnotI)){
			JI=utils::read.table(JI,header=FALSE)
			Lines=JI[which(JI[,1]==geneID),]
			JAnnot=utils::read.table(EandTrAnnot,header=FALSE,sep="\t",nrows=as.numeric(Lines[3]),skip=as.numeric(Lines[2]),stringsAsFactors=FALSE)		
		}
		else{
			
			JEAnnot=ETrAnnot[,-c(5,6,7)]
			JEAnnot=JEAnnot[!duplicated(JEAnnot),]
			JEAnnot=JEAnnot[order(JEAnnot[,2]),]
			JUCs=which(sapply(JEAnnot[,2],function(x) substr(x,1,1)=="J"))
			JUC=JEAnnot[JUCs,]
			PSR=JEAnnot[-c(JUCs),]
			
			
			Incl=unlist(sapply(JUC[,4],function(x) !grepl("[*]",x)))
			JUC=JUC[Incl,]
			
			if(!PartiallyAnnotated){
				Incl=unlist(sapply(PSR[,4],function(x) !grepl("[*]",x)))
				PSR=PSR[Incl,]
			}
			
			J=unique(JUC[,2])
			JAnnot=c()
			for(j in J){
					
					PSR3_Final=""
					PSR5_Final=""
					
					#Side 3 annots
					SubJ=JUC[which(JUC[,2]==j&JUC[,3]==3), ]
					TC=SubJ[1,1]
					Strand=SubJ[1,5]
					
					Es=unique(SubJ[,4])
					PSRs3=PSR[which(PSR[,4]%in%Es), ]
					
					Count=table(PSRs3[,2])
					PSR3=names(Count)[which.min(table(PSRs3[,2])-length(Es))]
					PSRs3temp=PSR[which(PSR[,2]%in%PSR3), 2]
					
					if((Strand=="-"|Strand==-1|Strand=="-1.0") & length(PSRs3temp)>0){
						PSR3_Final=PSRs3temp[1]
					}			
					else if((Strand=="+"|Strand==1|Strand=="1.0") & length(PSRs3temp)>0){
						PSR3_Final=PSRs3temp[length(PSRs3temp)]
					}	
					Row=c(TC,PSR3_Final,j,"3")
					JAnnot=rbind(JAnnot,Row)
					
					#Side 5 annots
					SubJ=JUC[which(JUC[,2]==j&JUC[,3]==5), ]
					TC=SubJ[1,1]
					Strand=SubJ[1,5]
					
					Es=unique(SubJ[,4])
					PSRs5=PSR[which(PSR[,4]%in%Es), ]
					
					Count=table(PSRs5[,2])				
					PSR5=names(Count)[which.min(table(PSRs5[,2])-length(Es))]
					PSRs5temp=PSR[which(PSR[,2]%in%PSR5), 2]
					
					if((Strand=="-"|Strand==-1|Strand=="-1.0") & length(PSRs5temp)>0){
						PSR5_Final=PSRs5temp[length(PSRs5temp)]
					}			
					else if((Strand=="+"|Strand==1|Strand=="1.0") & length(PSRs5temp)>0){
						PSR5_Final=PSRs5temp[1]
					}	
					Row=c(TC,PSR5_Final,j,"5")
					JAnnot=rbind(JAnnot,Row)
					
					#exclusions
					I1=0
					I2=0
					if(PSR3_Final!="" & PSR5_Final!="" & (Strand=="+"|Strand==1|Strand=="1.0")){
						I1=which(PSR[,2]==PSR3_Final)[length(which(PSR[,2]==PSR3_Final))]
						I2=which(PSR[,2]==PSR5_Final)[1]
					}
					else if(PSR3_Final!="" & PSR5_Final!="" & (Strand=="-"|Strand==-1|Strand=="-1.0")){
						I1=which(PSR[,2]==PSR5_Final)[length(which(PSR[,2]==PSR5_Final))]
						I2=which(PSR[,2]==PSR3_Final)[1]
					}
					if(I1!=0&I2!=0&I1!=(I2+1)&I2!=(I1+1)){
						if(I1<I2){
							R=seq(I1+1,I2,1)
						}
						else if(I1>I2){
							R=seq(I1-1,I2,-1)
						}
						ExclPSR_temp=PSR[R,]
						Exrows=nrow(ExclPSR_temp)
						
						if(Exrows!=0){
							ExclPSR=unique(ExclPSR_temp[,2])
							if(any(is.na(ExclPSR))){
								ExclPSR=ExclPSR[-c(which(is.na(ExclPSR)))]
							}
							for(e in ExclPSR){
								if(e!=PSR5_Final & e!=PSR3_Final){
									Row=c(TC,e,j,"exclusion")
									JAnnot=rbind(JAnnot,Row)
								}
								
							}
						}
					}
					
				}
				if(!is.null(JAnnot)){
				JAnnot=JAnnot[order(JAnnot[,2]),]
				colnames(JAnnot)=c("TC_ID","PSR_ID","JUC_ID","as_type")			
				JAnnot[,1]=as.character(JAnnot[,1])
				JAnnot[,2]=as.character(JAnnot[,2])
				JAnnot[,3]=as.character(JAnnot[,3])
				JAnnot[,4]=as.character(JAnnot[,4])
			}
			else{
				JAnnot=c()
			}
		}
		if(Juninfo=="Ensemble"){
			ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
			attributes.region<- c("chromosome_name", "start_position", "end_position", "ensembl_gene_id")							
			filter.symbol<- "hgnc_symbol"
			
			TrAnnot=c()
			
			trim <- function(s, ...) {
				s <- as.character(s);
				s <- sub("^[\t\n\f\r[:punct:]  ]*", "", s);
				s <- sub("[\t\n\f\r ]*$", "", s);
				s;
			} 
			
			trans=paste(geneID,".hg",sep="")
			exons=positionData[positionData$transcript_cluster_id == trans, "probeset_id"]
			
			PSR=sapply(exons, function(x) substr(x,1,3)!="JUC")
			exons=exons[PSR]
			trans=paste(trans,".1",sep="")
			symbol.to.annotate <- AnnotateGenes(trans,transcriptData)$symbol # get HBC identity
			if(all(symbol.to.annotate!="---")){
				ensembl.output     <- AnnotateGeneSymbol(symbol.to.annotate) # get region
				if(nrow(ensembl.output) > 0){
					gene  = makeGene(id = ensembl.output$ensembl_gene_id, biomart = ensembl)
					
					strand <- transcriptData[transcriptData[,1] == trans,][2]
					gene.positions   <- transcriptData[transcriptData$probeset_id %in% exons,]
					gene.positions$start=as.numeric(gene.positions$start)
					gene.positions$stop=as.numeric(gene.positions$stop)
					gene.positions<-gene.positions[order(gene.positions$start,decreasing=FALSE),]
					gene.positions[,1]=sapply(gene.positions[,1],function(x) strsplit(x,"[.]")[[1]][1])
					PSR=gene.positions[,1]
					G=gene@ens
					for(p in PSR){
						Start=gene.positions[which(gene.positions[,1]==p),3]
						Stop=gene.positions[which(gene.positions[,1]==p),4]
						for(r in 1:nrow(G)){
							if(G[r,4]<=Start&G[r,5]>=Stop){
								Row=c(geneID,p,strand,G[r,3],G[r,2])
								TrAnnot=rbind(TrAnnot,Row)
							}
						}
					}
					
					if(is.null(TrAnnot)){
						TrAnnot=c()
					}
					else{
						colnames(TrAnnot)=c("TC_ID","PSR_ID","strand","EAnnot","TrAnnot")
						
						TrAnnot=as.data.frame(TrAnnot)
						TrAnnot[,1]=as.character(TrAnnot[,1])
						TrAnnot[,2]=as.character(TrAnnot[,2])
						TrAnnot[,3]=as.character(TrAnnot[,3])
						TrAnnot[,4]=as.character(TrAnnot[,4])
						TrAnnot[,5]=as.character(TrAnnot[,5])
						TrAnnot=as.matrix(TrAnnot)
					}
					
				}
				else{
					TrAnnot=c()
				}
			}
			else{
				TrAnnot=c()
			}
		}
		
	}
	else{
		ETrAnnot=c()
		JAnnot=c()
	}
	
	if(length(JAnnot)==0|length(TrAnnot)==0){
		
		output=paste(geneID,"No valuable isoform composition information available",sep="\t")
		
		return(output)
		
	}
	
	
	colnames(DataS)[2]="ExonID"
	if(any(is.na(JAnnot[,2]))){
		D=which(is.na(JAnnot[,2]))
		JAnnot=JAnnot[-c(D),]
	}
	
	
	OutRange=matrix(0,nrow=length(unique(DataS$ExonID)),ncol=length(Groups))
	for(i in 1:length(unique(DataS$ExonID))){
		for(j in 1:length(Groups)){
			Subset=as.vector(as.matrix(DataS[which(DataS$ExonID==unique(DataS$ExonID)[i]),Groups[[j]]+2]))
			Range=range(Subset)
			OutRange[i,j]=Range[2]-Range[1]
		}	
	}
	rownames(OutRange)=unique(DataS$ExonID)
	
	PSRsExons=unique(DataS$ExonID)[which(substr(unique(DataS$ExonID),1,3)=="PSR")]
	Continue=TRUE
	TempRange=as.vector(as.matrix(OutRange))
	
	while(Continue){
		M=stats::median(TempRange)
		SD=stats::sd(TempRange)
		
		RangeCutOff=M+1.5*SD
		
		if(length(which(TempRange>RangeCutOff))==0){
			Continue=FALSE
		}	
		else{
			TempRange=TempRange[-c(which(TempRange>RangeCutOff))]
		}	
		
		if(RangeCutOff<2){
			RangeCutOff=2
		}
	}
	
	GroupData=list()
	DelJ=c()
	Rem=c()
	for(g in 1:length(Groups)){
		
		GData=DataS[,c(1,2,Groups[[g]]+2)]
		
		for(e in unique(GData$ExonID)){
			Subset=GData[which(GData$ExonID==e),-c(1,2)]
			Temp=as.vector(as.matrix(GData[which(GData$ExonID==e),-c(1,2)]))
			
			StopDel=c()
			if(length(Temp)<=12){
				StopDel=0
				next
			}
			else if(length(Temp)<18){
				StopDel=2
			}
			else{
				StopDel=3
			}
			
			if(round(OutRange[e,g],2)>RangeCutOff){
				Continue=TRUE	
				CutOff=c()
				Flagged=c()
				while(Continue){
					M=mean(Temp)
					SD=stats::sd(Temp)	
					Dist=abs(Temp-M)
					Test=Temp[which.max(Dist)]
					if(Test>=M){
						P=stats::pnorm(Test,m=M,sd=SD,lower.tail=FALSE)
					}
					else{
						P=stats::pnorm(Test,m=M,sd=SD,lower.tail=TRUE)
					}
					if(round(P,2)>=0.05){
						CutOff=c(CutOff)
						Continue=FALSE
					}	
					else{
						CutOff=c(CutOff,Temp[which.max(Dist)])
						Temp=Temp[-c(which.max(Dist))]
						R=range(Temp)
						RR=R[2]-R[1]
						if(round(RR,2)<=RangeCutOff){
							Continue=FALSE								
						}		
					}
				}
				Temp=Temp[-c(which.max(Dist))]
				R=range(Temp)
				RR=R[2]-R[1]
				if(round(RR,2)>2*RangeCutOff&substr(e,1,3)=="JUC"){
					DelJ=c(DelJ,e)
				}
				
				if(ncol(Subset)<=6){
					AtLeast=ncol(Subset)
				}
				else{
					AtLeast=ncol(Subset)-1
				}
				
				if(length(CutOff)>0){
					F=c()
					P=c()
					for(c in CutOff){
						if(c>=M){
							P=c(P,stats::pnorm(c,m=M,sd=SD,lower.tail=FALSE))
						}
						else{
							P=c(P,stats::pnorm(c,m=M,sd=SD,lower.tail=TRUE))
						}
						
						F=c(F,rownames(which(Subset==c,arr.ind=TRUE))[1])
					}
					FlagP=c()
					names(P)=F
					for(f in unique(F)){
						FlagP=c(FlagP,mean(P[f]))
					}
					names(FlagP)=unique(F)
					Flagged=names(which(table(F)>=AtLeast))
					FlagP=FlagP[Flagged]
				}
				
				if(length(Flagged)>0){
					#print(e)
					#print(length(Flagged))
					if(length(Flagged)>length(StopDel)){
						Flagged=names(sort(FlagP)[1:StopDel])
						Rem=c(Rem,Flagged)
					}
					#GData=GData[-c(which(rownames(GData)%in%Flagged)),]
				}
			}
		}
		
		
	}
	
	if(length(Rem)>0){
		if(any(is.na(Rem))){
			Rem=Rem[-c(which(is.na(Rem)))]
		}
		if(length(Rem)>0){
			for(g in 1:length(Groups)){
				GData=DataS[,c(1,2,Groups[[g]]+2)]
				GData=GData[-c(which(rownames(GData)%in%Rem)),]
				GroupData[[g]]=GData	
			}
		}
		else{
			for(g in 1:length(Groups)){
				GData=DataS[,c(1,2,Groups[[g]]+2)]
				GroupData[[g]]=GData	
			}
		}
	}
	else{
		for(g in 1:length(Groups)){
			GData=DataS[,c(1,2,Groups[[g]]+2)]
			GroupData[[g]]=GData	
		}
	}
	
	print(geneID)	
	
	DelJuncs=names(which(table(DelJ)==length(Groups)))
	Jucs=unique(JAnnot[,3])[which(!is.na(unique(JAnnot[,3])))]
	
	
	DABGPSR=unique(PSRsExons)[which(unique(PSRsExons)%in%Low_AllSamples)]
	if(length(DABGPSR)>0){
		PSRsExons=PSRsExons[-c(which(PSRsExons%in%DABGPSR))]
	}
	
	Const=unique(PSRsExons)[which(!unique(PSRsExons)%in%ASPSR)]
	AS=unique(PSRsExons)[which(unique(PSRsExons)%in%ASPSR)]
	JAsses=data.frame(matrix(0,nrow=length(Jucs),ncol=(3+length(Groups))))
	colnames(JAsses)=c("Pattern","Flat","Low",paste("Low-",c(1:length(Groups)),sep=""))
	Hold=1
	Connections=data.frame(matrix(0,ncol=(2+length(Groups)),nrow=(nrow(JAsses)*length(Groups))))
	Place=1
	#ExonAssesment=data.frame(matrix(0,ncol=2,nrow=length(unique(JAnnot$PSR_ID))))
	#rownames(ExonAssesment)=unique(JAnnot$PSR_ID)
	ASJ=c()
	ExclDef=c()
	
	# Assesment of junctions
	if(length(Jucs)>0){
		for(k in 1:length(Jucs)){
			j=Jucs[k]
			
			if(is.na(j)){
				next
			}
			
			JValues=list()
			JList=list()
			JLengths=list()
			JAssesP=c()
			JAssesV=c()
			JFlat=c()
			JLow=c()
			
			#valid probes?
			if(nrow(DataS[which(DataS$ExonID==j),-c(1,2)])<=3){
				JAssesP=c(JAssesP,"Junction has too few valid probes",rep("-",(2+length(Groups))))
				JAsses[k,]=JAssesP
				Hold=Hold+1
				next				
			}
			
			#Low_AllSamples
			DABG=TRUE
			if(j%in%Low_AllSamples){
				JAssesP=c(JAssesP,"Junction is not DABG",c("-",TRUE,rep("-",length(Groups))))
				#DABG=FALSE	
			}
			
			if(j%in%DelJuncs){
				JAssesP=c(JAssesP,"Junction has a too large spread of values",rep("-",(2+length(Groups))))
				JAsses[k,]=JAssesP
				Hold=Hold+1
				next
			}
			
			if(j%in%unique(DataS$ExonID)){
				JD=list()
				for(g in 1:length(Groups)){
					JD[[g]]=as.vector(as.matrix(GroupData[[g]][which(GroupData[[g]]$ExonID==j),-c(1,2)]))
				}	
				#JUC_Ranks=sort(JD,index.return=TRUE)$ix
				JUC_Ranks=rank(unlist(JD),ties.method="random")
				JList[[length(JList)+1]]=JUC_Ranks
				names(JList)[length(JList)]=j
				JLengths=c(JLengths,length(JUC_Ranks))
				JValues=list(JD)
				names(JValues)=j
			}
			
			Set=sort(JAnnot[which(JAnnot[,3]==j&(JAnnot[,4]!="exclusion")),2])
			L=list()
			D=list()
			if(!all(Set%in%DataS$ExonID)){
				Set=Set[-c(which(!Set%in%DataS$ExonID))]
			}
			if(length(Set)==2&length(unique(Set))==1){
				print(j)
			}
			if(length(Set)>1){
				for(s in Set){	
					PSR=list()
					for(g in 1:length(Groups)){
						PSR[[g]]=as.vector(as.matrix(GroupData[[g]][which(GroupData[[g]]$ExonID==s),-c(1,2)]))
					}
					D[[s]]=PSR
					Ranks=list(rank(unlist(PSR),ties.method="random"))
					L=c(L,Ranks)
				}
				names(L)=Set	
			}
			else{
				JAssesP=c(JAssesP,"Junction does not have 2 anchor points",rep("-",(2+length(Groups))))
				JAsses[k,]=JAssesP
				Hold=Hold+1
				next			
			}
			L=c(L,JList)
			D=c(D,JValues)
			
			LL=sapply(unlist(D,recursive=FALSE),length)
			if(length(unique(LL))>1){
				MinLength=min(LL)
				Index=which(LL>MinLength)
				for(sj in 1:length(Index)){
					si=Index[sj]
					t=MinLength/(length(Groups[[1]]))
					Dnew=c()
					temp=c()
					g=as.numeric(substr(names(si),nchar(names(si)),nchar(names(si))))
					p=substr(names(si),1,nchar(names(si))-1)
					temp=GroupData[[g]][which(GroupData[[g]]$ExonID==p),-c(1,2)]
					M=apply(as.matrix(temp),2,stats::median)
					for(c1 in 1:ncol(temp)){
						Dist=abs(temp[[c1]]-M[c1])
						I=rank(Dist)
						Select=which(I<=t)
						if(length(Select)<t){
							Add=t-length(Select)
							SelectAdd=which((I>t)&(I<=(t+Add)))[c(1:Add)]
							Select=c(Select,SelectAdd)
						}
						else if(length(Select)>t){
							Del=length(Select)-t
							SelectDel=which(Select>=t)[c(1:Del)]
							Select=Select[-c(SelectDel)]
						}
						NewC=temp[Select,c1]
						Dnew=cbind(Dnew,NewC)
					}
					D[[p]][[g]]=as.vector(as.matrix(Dnew))
					PSR=unlist(D[[p]])
					#Ranks=list(sort(PSR,index.return=TRUE)$ix)
					Ranks=rank(PSR)
					L[[p]]=Ranks	
					
				}
			}
			
			if(length(Set)>1){
				L[[1]]=(L[[1]]+L[[2]])/2
				L=L[-c(2)]
				L[[1]]=sort(L[[1]],index.return=TRUE)$ix
			}	
			Ranks=do.call("cbind",L)
			
			Y=as.vector(as.matrix(Ranks))
			exon=c()
			for(c in 1:ncol(Ranks)){
				exon=c(exon,rep(c,nrow(Ranks)))
			}
			tissuetemp=c()
			for(g in 1:length(Groups)){
				tissuetemp=c(tissuetemp,rep(g,nrow(Ranks)/(length(Groups))))
			}
			tissue=rep(tissuetemp,ncol(Ranks))
			
			exon=as.factor(exon)
			tissue=as.factor(tissue)
			ft1<-stats::lm(Y~exon*tissue-1)
			
			Inter=summary(ft1)$coefficients[((length(L)+2):nrow(summary(ft1)$coefficients)),4]
			
			if(all(round(Inter,2)>=0.05)){
				# Junction is product of its supporting exons
				JAssesP=c(JAssesP,"Pattern Supported")						
			}
			else{
#					if(DABG==FALSE){
#						JAssesP=c(JAssesP,"Junction is not DABG",c("-",TRUE,rep("-",length(Groups))))
#						JAsses[k,]=JAssesP
#						Hold=Hold+1
#						next
#					}
				JAssesP=c(JAssesP,"Pattern Not Supported")
			}
			
			#Junction Value Assessment
			Set=sort(JAnnot[which(JAnnot[,3]==j&(JAnnot[,4]!="exclusion")),2])
#				
			#Flat line
			tissuetemp=c()
			for(g in 1:length(Groups)){
				tissuetemp=c(tissuetemp,rep(g,length(D[[1]][[1]])))
			}
			Flattest=stats::lm(unlist(D[[length(D)]])~tissuetemp-1)
			Flat=summary(Flattest)$coefficients[1,4]
			if(round(Flat,2)>0.05){
				JFlat=c(JFlat,TRUE)
			}
			else{
				JFlat=c(JFlat,FALSE)
				
			}	
			
			#Low??
			Low=rep(FALSE,length(Groups))
			JLow=FALSE
			if(j%in%unlist(Low_GSamples)){
				R=lapply(Low_GSamples,function(x) j%in%x)
				R=unlist(R)
				if(all(R)){
					if(JFlat){
						Low[which(R)]=TRUE
					}	
				}
				else if(any(R)){
					Low[which(R)]=TRUE
				}
			}					
			Row=c(JAssesP,JFlat,JLow,Low)
			JAsses[k,]=Row
			
		}
		JAsses=cbind(Jucs,JAsses)
		
		
		# PSR reflection on junctions
		
		for(r in 1:nrow(JAsses)){
			R=JAsses[r,]
			J=as.character(JAsses[r,1])
			PSRs=unique(JAnnot[,2][which(JAnnot[,3]==J)])
			supp=sort(JAnnot[,2][which(JAnnot[,2]%in%PSRs&JAnnot[,3]==J&JAnnot[,4]!="exclusion")])
			if(!all(supp%in%DataS$ExonID)){
				supp=supp[-c(which(!supp%in%DataS$ExonID))]
			}
			excl=sort(JAnnot[,2][which(JAnnot[,2]%in%PSRs&JAnnot[,3]==J&JAnnot[,4]=="exclusion")])
			if(!all(excl%in%DataS$ExonID)){
				excl=excl[-c(which(!excl%in%DataS$ExonID))]
			}
			if(R[2]%in%c("Junction has too few valid probes","Junction has a too large spread of values","Junction does not have 2 anchor points")){
				next
			}
			#Do the support points connect? Check Low, 1-Low and 2-Low
			for(g in 1:length(Groups)){
				
				GData=GroupData[[g]]
				
				if(R[2]=="Junction is not DABG"){
					Connections[Place,c(1,2,g+2)]=c(J,paste(supp,collapse="-"),"never")
					Place=Place+1
				}
				
				else if(as.logical(R[5+g-1])){
					Connections[Place,c(1,2,g+2)]=c(J,paste(supp,collapse="-"),"never")
					Place=Place+1
					#junction is low	
				}	
				else{
					Connections[Place,c(1,2,g+2)]=c(J,paste(supp,collapse="-"),"present")
					Place=Place+1
					#junction is present ==> both linking points are present and connection is present
					#if excl junction: exon is indeed excluded => AS
					
					if(length(excl)>0){					
						for(e in excl){
							ExclDef=c(ExclDef,e)
							
							if(JAsses[r,2]=="Pattern Not Supported"){
								
								JAsses[r,2]="Pattern Supported"
							}		
						}
					}		
				}	
			}	
		}		
	}
	
	
	if(any(Connections[,1]==0)){
		Connections=Connections[-c(which(Connections[,1]==0)),]
	}
	if(any(duplicated(Connections))){
		Connections=Connections[!duplicated(Connections),]
	}
	Links=c()		
	UL=unique(Connections[,c(1,2)])
	for(c in 1:nrow(UL)){
		rowLinks=c(cbind(UL[c,1],UL[c,2]))
		Set=which(Connections[,1]==UL[c,1]&Connections[,2]==UL[c,2])
		for(g in 1:length(Groups)){
			GAsses=Connections[Set,g+2]
			Get=GAsses[which(GAsses!=0)]
			if(all(Get=="never")){
				rowLinks=c(rowLinks,"never")
			}
			else{
				rowLinks=c(rowLinks,"present")
			}
		}
		Links=rbind(Links,rowLinks)
		
	}	
	
	Edges=c()
	#Exons first
	Exons=sort(unique(c(unique(JAnnot[,2]),c(unique(PSRsExons),unique(DABGPSR)))))
	Exons=Exons[which(Exons%in%DataS[,2])]
	if(length(Exons)==0){
		output=paste(geneID,"No PSR's present",sep="\t")
		
		return(output)
		
	}
	for(i in 1:(length(Exons)-1)){
		R=c(Exons[i],Exons[i+1])
		Edges=rbind(Edges,R)
	}
	col1=rep("grey",nrow(Edges))
	width1=rep(2,nrow(Edges))
	G=list()
	Cols=list()
	Widths=list()
	for(g in 1:(length(Groups)+1)){
		G[[g]]=Edges
		Cols[[g]]=col1	
		Widths[[g]]=width1
	}
	#Links
	if(length(Jucs)>0&length(Links)>0){
		for(l in 1:nrow(Links)){
			E=strsplit(Links[l,2],"-")[[1]]
			if(length(E)>1){
				for(g in 1:(length(Groups)+1)){
					G[[g]]=rbind(G[[g]],E)
					if(g==1){
						Cols[[1]]=c(Cols[[1]],"blue")
						Widths[[1]]=c(Widths[[1]],2)
					}
					else{
						L=Links[l,g+1]
						if(L=="present"){
							Cols[[g]]=c(Cols[[g]],"blue")
							Widths[[g]]=c(Widths[[g]],2)
						}
						else{
							Cols[[g]]=c(Cols[[g]],"red")
							Widths[[g]]=c(Widths[[g]],2)
						}
					}
				}			
			}
		}
		#rownames(Edges)=1:nrow(Edges)
	}
	
	if(length(Groups)==2){
		FC=c()
		for(e in Exons){				
			A=mean(as.vector(as.matrix(GroupData[[1]][which(GroupData[[1]]$ExonID==e),-c(1,2)])))
			B=mean(as.vector(as.matrix(GroupData[[2]][which(GroupData[[2]]$ExonID==e),-c(1,2)])))
			FC=c(FC,A-B)
		}	
		
		names(FC)=Exons
		
		
		if(length(Const)==0){
			MedFC=stats::median(FC)
			SDFC=0.5
		}
		else{
			ConstFC=FC[Const]
			MedFC=stats::median(ConstFC)
			SDFC=stats::sd(ConstFC)
			if(is.na(SDFC)){
				SDFC=0.5
			}
		}
		ASFC=FC[AS]
		ASPvals=c()
		for(a in FC){
			if(a>MedFC){
				ASPvals=c(ASPvals,stats::pnorm(a,MedFC,SDFC,lower.tail=FALSE))
			}
			else{
				ASPvals=c(ASPvals,stats::pnorm(a,MedFC,SDFC,lower.tail=TRUE))
			}	
		}
		names(ASPvals)=names(FC)
		ASPvals=stats::p.adjust(ASPvals,"fdr")
		ASfctemp=names(which(ASPvals<0.05))
		ASfc=ASfctemp[which(!ASfctemp%in%c(AS,ASJ))]
	}
	else{
		FC=rep("-",length(Exons))
		ASfc=NULL
	}
	colNodes=rep("grey",length(Exons))
	names(colNodes)=Exons
	ColN=list()
	for(g in 1:(length(Groups)+1)){
		ColN[[g]]=colNodes			
	}
	for(c in Exons){
		if(c%in%DABGPSR){
			for(g in c(1:length(Groups)+1)){
				ColN[[g]][c]="grey"
			}	
		}
		else if(c%in%Const){
			for(g in c(1:length(Groups)+1)){
				ColN[[g]][c]="black"
			}	
		}
		else{
			if(length(Groups)==2){
				#if(round(ASPvals[c],2)<0.05){
				if(FC[c]<0){
					ColN[[2]][c]="red"
					ColN[[3]][c]="green"
				}
				else{
					ColN[[2]][c]="green"
					ColN[[3]][c]="red"
				}
				#}
			}
		}	
	}
	
	
	DelJ=unique(c(DelJ,as.character(JAsses[,1][which(JAsses[,2]%in%c("Junction has too few valid probes","Junction has a too large spread of values","Junction does not have 2 anchor points"))])))
	OutList=data.frame("TC_ID"=rep(DataS[1,1],length(Exons)),"PSR_ID"=Exons,"Type"=rep(0,length(Exons)),"Unreliable Junctions"=rep(0,length(Exons)),"Linking Junctions"=rep(0,length(Exons)),"Exclusion Junctions"=rep(0,length(Exons)),"Supported by"=rep(0,length(Exons))
			,"Identified by"=rep(0,length(Exons)),"Fold Change"=rep(0,length(Exons)),"Exons"=rep(0,length(Exons)))
	
	#Reduce transcripts in all samples:
	NeverPresent=which(apply(Links[,3:ncol(Links),drop=FALSE],1,function(x) all(x==rep("never",length(Groups)))))
	#NeverLinkedPSR=list()
	NeverLinkedTr=list()
	if(length(NeverPresent)>0){
		RLinks=Links[NeverPresent,,drop=FALSE]		
		for(r in 1:nrow(RLinks)){
			PSR=RLinks[r,2]
			PSRs=strsplit(PSR,"-")[[1]]
			set=c()
			for(t in unique(TrAnnot[,5])){
				PSRset=TrAnnot[which(TrAnnot[,5]==t),2]
				if(all(PSRs%in%PSRset)){
					set=c(set,t)
				}
			}
			NeverLinkedTr[[r]]=set
#			set=TrAnnot[which(TrAnnot[,2]==RLinks[r,1]),5]
#			#set=strsplit(RLinks[r,2],"-")[[1]]		
#			#NeverLinkedPSR[[r]]=set
#			NeverLinkedTr[[r]]=set
		}
		NeverLinkedTr=unique(unlist(NeverLinkedTr))
	}
	
	
	#Event Type
	TC=DataS[1,1]
	Strand=TrAnnot[1,3]
	TRS=list()
	N=c()
	for(tr in unique(TrAnnot[,5])){
		if(tr%in%NeverLinkedTr){
			next
		}
		Set=TrAnnot[which(TrAnnot[,5]==tr),2]
		PSet=Set[which(substr(Set,1,1)=="P")]
		ESete=TrAnnot[which(TrAnnot[,5]==tr),4]
		ESet=ESete[which(substr(Set,1,1)=="P")]
#			Remove=FALSE
#			if(length(NeverLinkedPSR)>0){
#				for(nl in NeverLinkedPSR){
#					Pos1=which(PSet==nl[1])
#					Pos2=which(PSet==nl[2])
#					if(length(Pos1)>0&length(Pos2)>0){
#						if(Pos1+1==Pos2){
#							#linked in transcript but the link is not present thus transcript can be removed from the collection
#							Remove=TRUE
#						}
#					}
#				}
#			}
		if(length(DABGPSR)>0){
			Stop=FALSE
			for(d in DABGPSR){
				if(any(PSet==d)){
					Stop=TRUE
					break					
				}
			}
			if(Stop){
				next
			}
		}
		
		if(Strand=="-"|Strand==-1){
			PSet=rev(sort(PSet))
			ESet=rev(sort(ESet))
		}
		names(ESet)=PSet
		TRS[[length(TRS)+1]]=ESet	
		N=c(N,tr)
	}
	names(TRS)=N
	if(length(TRS)==0){
		
		output=list("Output"=OutList)
		
		return(output)
	
	}
	
	AllExons=unique(unlist(TRS))		
	for(e in Exons){
		EAnnotp=TrAnnot[which(TrAnnot[,2]==e),4]
		EAnnotp=EAnnotp[which(EAnnotp%in%AllExons)]
		EAnnotp=paste(EAnnotp,collapse="|")
		s=JAnnot[which(JAnnot[,2]==e),c(3,4),drop=FALSE]
		Dels=s[which(s[,1]%in%DelJ),1]
		if(length(Dels)>0){
			s=s[which(!s[,1]%in%DelJ),,drop=FALSE]
		}
		else{
			Dels="-"
		}
		if(e%in%DABGPSR){
			r=c("not DABG",rep("-",6),EAnnotp)
		}
		
		else if(e%in%Const){			
			if(!(all(is.na(s)))){
				Jss=c()
				Jsse=c()
				Js=s[which(s[,2]%in%c(3,5)),1]
				if(length(Js)>0){
					Jss=as.character(JAsses[which(JAsses[,1]%in%Js&JAsses[,2]=="Pattern Supported"&JAsses[,4]==FALSE),1])
					if(!(all(Js%in%Jss))){
						Dels=c(Dels,Js[which(!Js%in%Jss)])
						if(length(Jss)==0){
							Jss="-"
						}
					}
				}
				else{
					Jss="-"
				}
				Jse=s[which(s[,2]=="exclusion"),1]
				if(length(Jse)>0){
					Jsse=as.character(JAsses[which(JAsses[,1]%in%Jse&JAsses[,4]==FALSE),1])
				}
				else{
					Jsse="-"
				}
				SupportedBy=c(Jss,Jsse)
				if(all(SupportedBy=="-")){
					SupportedBy="-"
				}
				else if(all(SupportedBy%in%s[,1])){
					SupportedBy="All"
				}	
				else if(!any(SupportedBy%in%s[,1])){
					SupportedBy="None"
				}
				else{
					if(any(SupportedBy=="-")){
						SupportedBy=SupportedBy[-c(which(SupportedBy=="-"))]
					}
				}
				if(length(Dels)>1&any(Dels=="-")){
					Dels=Dels[-c(which(Dels=="-"))]
				}
				
				r=c("Const",paste(Dels,collapse="|"),paste(Jss,collapse="|"),paste(Jsse,collapse="|"),paste(SupportedBy,collapse="|"),"-",FC[e],EAnnotp)
			}	
			else{
				r=c("Const",rep("-",5),FC[e],EAnnotp)
			}
		}
		else if(e%in%AS){
			
			if(!(all(is.na(s)))){
				Jss=c()
				Jsse=c()
				Js=s[which(s[,2]%in%c(3,5)),1]
				if(length(Js)>0){
					Jss=JAsses[which(JAsses[,1]%in%Js&JAsses[,2]=="Pattern Supported"&JAsses[,4]==FALSE),1]
					if(!(all(Js%in%Jss))){
						Dels=c(Dels,Js[which(!Js%in%Jss)])
						if(length(Jss)==0){
							Jss="-"
						}
					}
				}
				else{
					Jss="-"
				}
				Jse=s[which(s[,2]=="exclusion"),1]
				if(length(Jse)>0){
					Jsse=JAsses[which(JAsses[,1]%in%Jse&JAsses[,4]==FALSE),1]
				}
				else{
					Jsse="-"
				}
				SupportedBy=c(as.character(Jss),as.character(Jsse))
				if(all(SupportedBy=="-")){
					SupportedBy="-"
				}
				else if(all(SupportedBy%in%s[,1])){
					SupportedBy="All"
				}	
				else if(!any(SupportedBy%in%s[,1])){
					SupportedBy="None"
				}	
				else{
					if(any(SupportedBy=="-")){
						SupportedBy=SupportedBy[-c(which(SupportedBy=="-"))]
					}
				}
				if(length(Dels)>1&any(Dels=="-")){
					Dels=Dels[-c(which(Dels=="-"))]
				}
				r=c("AS",paste(Dels,collapse="|"),paste(Jss,collapse="|"),paste(Jsse,collapse="|"),paste(SupportedBy,collapse="|"),"REIDS",FC[e],EAnnotp)
			}	
			else{
				r=c("AS",rep("-",4),"REIDS",FC[e],EAnnotp)
			}
		}
		else if(e%in%ASJ){
			if(!(all(is.na(s)))){
				Jss=c()
				Jsse=c()
				Js=s[which(s[,2]%in%c(3,5)),1]
				if(length(Js)>0){
					Jss=JAsses[which(JAsses[,1]%in%Js&JAsses[,2]=="Pattern Supported"&JAsses[,4]==FALSE),1]
					if(!(all(Js%in%Jss))){
						Dels=c(Dels,Js[which(!Js%in%Jss)])
						if(length(Jss)==0){
							Jss="-"
						}
					}
				}
				else{
					Jss="-"
				}
				Jse=s[which(s[,2]=="exclusion"),1]
				if(length(Jse)>0){
					Jsse=JAsses[which(JAsses[,1]%in%Jse&JAsses[,4]==FALSE),1]
				}
				else{
					Jsse="-"
				}
				SupportedBy=c(as.character(Jss),as.character(Jsse))
				if(all(SupportedBy=="-")){
					SupportedBy="-"
				}
				else if(all(SupportedBy%in%s[,1])){
					SupportedBy="All"
				}	
				else if(!any(SupportedBy%in%s[,1])){
					SupportedBy="None"
				}
				else{
					if(any(SupportedBy=="-")){
						SupportedBy=SupportedBy[-c(which(SupportedBy=="-"))]
					}
				}
				if(length(Dels)>1&any(Dels=="-")){
					Dels=Dels[-c(which(Dels=="-"))]
				}
				r=c("AS",paste(Dels,collapse="|"),paste(Jss,collapse="|"),paste(Jsse,collapse="|"),paste(SupportedBy,collapse="|"),"Junction Support",FC[e],EAnnotp)
			}
			else{
				r=c("AS",rep("-",4),"Junction Support",FC[e],EAnnotp)
			}	
		}
		else if(e%in%ASfc){
			
			if(!(all(is.na(s)))){
				Jss=c()
				Jsse=c()
				Js=s[which(s[,2]%in%c(3,5)),1]
				if(length(Js)>0){
					Jss=JAsses[which(JAsses[,1]%in%Js&JAsses[,2]=="Pattern Supported"&JAsses[,4]==FALSE),1]
					if(!(all(Js%in%Jss))){
						Dels=c(Dels,Js[which(!Js%in%Jss)])
						if(length(Jss)==0){
							Jss="-"
						}
					}
				}
				else{
					Jss="-"
				}
				Jse=s[which(s[,2]=="exclusion"),1]
				if(length(Jse)>0){
					Jsse=JAsses[which(JAsses[,1]%in%Jse&JAsses[,4]==FALSE),1]
				}else{
					Jsse="-"
				}
				SupportedBy=c(as.character(Jss),as.character(Jsse))
				if(all(SupportedBy=="-")){
					SupportedBy="-"
				}
				else if(all(SupportedBy%in%s[,1])){
					SupportedBy="All"
				}	
				else if(!any(SupportedBy%in%s[,1])){
					SupportedBy="None"
				}	
				else{
					if(any(SupportedBy=="-")){
						SupportedBy=SupportedBy[-c(which(SupportedBy=="-"))]
					}
				}
				if(length(Dels)>1&any(Dels=="-")){
					Dels=Dels[-c(which(Dels=="-"))]
				}
				r=c("AS",paste(Dels,collapse="|"),paste(Jss,collapse="|"),paste(Jsse,collapse="|"),paste(SupportedBy,collapse="|"),"Fold Change",FC[e],EAnnotp)
			}	
			else{
				r=c("AS",rep("-",4),"Fold Change",FC[e],EAnnotp)
			}
			
		}
		else{
#				if(e%in%DataS[,2]&(!e%in%JAnnot[,2])){
#					r=c("INI Filtered",rep("-",6))
#				}
			#else{
			print(e)
			#}
		}
		OutList[OutList[,2]==e,c(3:10)]=r
	}
	
	L=matrix(0,ncol=3,nrow=nrow(Links))
	for(i in 1:nrow(Links)){
		r=Links[i,c(2)]
		PSRs=strsplit(r,"-")[[1]]
		L[i,]=c(Links[i,1],PSRs[1],PSRs[2])
	}
	
	Table1=matrix(0,nrow=length(TRS)*2,ncol=4)
	l=c()
	for(i in 1:length(TRS)){
		N=names(TRS)[i]	
		l=c(l,L[which(L[,2]%in%names(TRS[[i]])&L[,3]%in%names(TRS[[i]])),1])
		J=paste(L[which(L[,2]%in%names(TRS[[i]])&L[,3]%in%names(TRS[[i]])),1],collapse="|")
		R1=paste(names(TRS[[i]]),collapse="|")
		R2=paste(TRS[[i]],collapse="|")
		Table1[i*2-1,]=c(TC,N,J,R1)
		Table1[i*2,]=c(TC,N,J,R2)
	}
	l=unique(l)
	if(any(!Links[,1]%in%l)){
		NovelConnections=Links[which(!Links[,1]%in%l),1]
		Novel=cbind(NovelConnections,Links[which(Links[,1]%in%NovelConnections),])
		NovelConnections=cbind(TC,NovelConnections)
	}	
	else{
		NovelConnections=c()
	}
	
	
	#utils::write.table(Table1, "Transcripts.txt", sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
	#Transcripts per Groups
	
	GroupTranscripts=list()
	for(g in 1:length(Groups)){
		N=c()
		GroupTranscripts[[g]]=list()
		NeverPresent=which(sapply(Links[,g+2],function(x) all(x=="never")))
		RLinks=Links[NeverPresent,,drop=FALSE]
		#NeverLinkedPSR=list()
		NeverLinkedTr=list()
		if(nrow(RLinks)>0){
			for(r in 1:nrow(RLinks)){
				set=c()
				for(t in unique(TrAnnot[,5])){
					PSRset=TrAnnot[which(TrAnnot[,5]==t),2]
					if(all(PSRs%in%PSRset)){
						set=c(set,t)
					}
				}
				#set=Transcripts[which(Transcripts[,2]==RLinks[r,1]),5]
				NeverLinkedTr[[r]]=set
				#set=strsplit(RLinks[r,1],"-")[[1]]
				#NeverLinkedPSR[[r]]=set
			}
			NeverLinkedTr=unique(unlist(NeverLinkedTr))
		}	
		for(tr in unique(TrAnnot[,5])){
			if(tr%in%NeverLinkedTr){
				next
			}
			Set=TrAnnot[which(TrAnnot[,5]==tr),2]
			PSet=Set[which(substr(Set,1,1)=="P")]
			ESete=TrAnnot[which(TrAnnot[,5]==tr),4]
			ESet=ESete[which(substr(Set,1,1)=="P")]
			
			
			if(length(Low_GSamples[[g]])>0){
				if(any(PSet%in%Low_GSamples[[g]])){
					next
				}
			}
			
			if(Strand=="-"|Strand==-1){
				PSet=rev(sort(PSet))
				ESet=rev(sort(ESet))
			}
			
			GroupTranscripts[[g]][[length(GroupTranscripts[[g]])+1]]=ESet	
			N=c(N,tr)
		}
		names(GroupTranscripts[[g]])=N
		
	}
	
	Table2=matrix(0,nrow=length(TRS),ncol=(2+length(Groups)))
	for(i in 1:length(TRS)){
		N=names(TRS)[i]
		Present=rep("no",length(Groups))
		for(g in 1:length(Groups)){
			if(N%in%names(GroupTranscripts[[g]])){
				Present[g]="Yes"
			}
		}
		Table2[i,]=c(TC,N,Present)
	}
	#utils::write.table(Table2, "Transcripts_Groups.txt", sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
	
	
	if(Strand=="-"|Strand==-1){
		OutList=OutList[nrow(OutList):1,]
	}
	
	Decision=OutList[,3]
	names(Decision)=OutList[,2]
	C=names(Decision)[which(Decision=="Const")]
	
	
	Category=rep("",nrow(OutList))
	
	for(r in 1:nrow(OutList)){
		if(Decision[r]=="Const"){
			Category[r]="-"
			next
		}
		else if(Decision[r]=="not DABG"){
			Category[r]="-"
			next
		}
		PSR=as.character(OutList[r,2])
		EAnnotp=TrAnnot[which(TrAnnot[,2]==PSR),,drop=FALSE]
		if(nrow(EAnnotp)==0){
			Category[r]="Intron Retention"	
			next
		}	
		else{
			Temp=c()
			for(ea in 1:nrow(EAnnotp)){
				OtherPSR=as.character(TrAnnot[which(TrAnnot[,4]==as.character(EAnnotp[ea,4])),2])
				OtherPSR=OtherPSR[which(OtherPSR!=PSR)]
				if(length(OtherPSR)>0){					
					S=substr(OtherPSR,1,3)
					Temp=c(Temp,OtherPSR[which(S=="PSR")])
				}
				else{
					Temp=c(Temp,"")
				}
				
			}
			if(any(Temp!="")){
				#exon has more than 1 probe set annotated to itself
				Temp=unique(Temp)
				if(any(Temp=="")){
					Temp=Temp[-c(which(Temp==""))]
				}
				if(all(which(OutList[,2]%in%Temp)>r)&r==1){
					Category[r]="Alternative First'"
					next
				}
				else if(all(which(OutList[,2]%in%Temp)>r)){
					Category[r]="Alternative 5'"
					for(tra in 1:length(TRS)){
						if(names(TRS[[tra]])[1]==PSR){
							Category[r]="Alternative First"
							break
						}
					}
					next
				}
				else if(all(which(OutList[,2]%in%Temp)<r)&r==nrow(OutList)){
					Category[r]="Alternative Last"
					next
				}
				else if(all(which(OutList[,2]%in%Temp)<r)){
					Category[r]="Alternative 3'"
					for(tra in c(1:length(TRS))){
						if(names(TRS[[tra]])[length(TRS[[tra]])]==PSR){
							Category[r]="Alternative Last"
							break
						}
					}
					next
				}
				else{
					Category[r]="Complex Event"
					next
				}
				
			}
			else{
				NeighboursIn=c()
				NeighboursOut=c()
				Positions=c()
				for(tra in 1:length(TRS)){
					if(length(TRS[[tra]])==0){
						next
					}
					
					if(names(TRS[[tra]])[1]==PSR){
						Category[r]="Alternative First"
						break
					}
					else if(names(TRS[[tra]])[length(TRS[[tra]])]==PSR){
						Category[r]="Alternative Last"
						break
					}
					else{
						Pos=which(names(TRS[[tra]])==PSR)
						if(length(Pos)>0){
							Left=names(TRS[[tra]])[Pos-1]
							Right=names(TRS[[tra]])[Pos+1]
							NeighboursIn=c(NeighboursIn,Left,Right)	
						}
						else if(names(TRS[[tra]])[1]%in%OutList[,2]&names(TRS[[tra]])[length(TRS[[tra]])]%in%OutList[,2]){
							if(which(OutList[,2]==names(TRS[[tra]])[1])<r&which(OutList[,2]==names(TRS[[tra]])[length(TRS[[tra]])])>r){				
								F=c(names(TRS[[tra]]),PSR)
								F=sort(F)
								Pos=which(F==PSR)
								Left=F[Pos-1]
								Right=F[Pos+1]
								NeighboursOut=c(NeighboursOut,Left,Right)	
							}
						}
						
					}
					
				}
				
				if(Category[r]==""){
					if(length(NeighboursIn)>0&length(NeighboursOut)>0){
						NIN=sort(unique(NeighboursIn))
						NOUT=sort(unique(NeighboursOut))
						
						if(length(NIN)==2&length(NOUT)==2){
							if(all(NIN%in%NOUT)){
								if(Decision[OutList[which(OutList[,2]==NIN[1]),2]]=="Const"&Decision[OutList[which(OutList[,2]==NIN[2]),2]]=="Const"){
									Category[r]="Cassette Exon"
								}
								else{
									Category[r]="Complex Event"
								}
							}
							
							else if(NIN[1]==NOUT[1]&NOUT[2]==OutList[r+1,2]&NIN[2]!=OutList[r+1,2]){
								if(r<(length(Decision)-2)){
									if(Decision[OutList[r+2,2]]=="Const"){
										Category[r]="Mutually Exclusive Exon"
									}
									else{
										Category[r]="Complex Event"
									}
								}
								else{
									Category[r]="Complex Event"
								}
							}
							
							else if(NIN[2]==NOUT[2]&NOUT[1]==OutList[r-1,2]&NIN[1]!=OutList[r-1,2]){
								if(r>2){
									if(Decision[OutList[r-2,2]]=="Const"){
										Category[r]="Mutually Exclusive Exon"
									}	
									else{
										Category[r]="Complex Event"
									}
								}
								else{
									Category[r]="Complex Event"
								}
								
								
							}
							else{
								Category[r]="Complex Event"
							}
						}
						else{
							Category[r]="Complex Event"
						}
					}
					else if(length(NeighboursIn)>0){
						NIN=sort(unique(NeighboursIn))
						if(length(NIN)==2){
							if(NIN[1]%in%OutList[,2]&NIN[2]%in%OutList[,2]){
								if(Decision[OutList[which(OutList[,2]==NIN[1]),2]]=="Const"&Decision[OutList[which(OutList[,2]==NIN[2]),2]]=="Const"){
									Category[r]="Cassette Exon"
								}
								else{
									Category[r]="Complex Event"
								}
							}
							else{
								Category[r]="Complex Event"
							}
						}
						else{
							Category[r]="Complex Event"
						}
					}	
					else{
						Category[r]="Complex Event"
					}
				}
				else{
					next
				}
			}
		}
	}
	
	if(Plot){		
		#pdf("Figures/Tr_TC17001298.pdf",width=18,height=10)
		graphics::par(mfrow=c(2,1))
		graphics::par(mar=c(12,1,1,1))
		for(g in 2:3){
			if(Strand=="-"|Strand==-1){
				R=G[[g]][which(rownames(G[[g]])=="R"),]
				R=R[nrow(R):1,c(2,1)]
				E=G[[g]][which(rownames(G[[g]])=="E"),]
				E=E[nrow(E):1,c(2,1)]
				G1=rbind(R,E)
				
				C1=rev(Cols[[g]][which(rownames(G[[g]])=="R")])
				C2=rev(Cols[[g]][which(rownames(G[[g]])=="E")])
				C=c(C1,C2)
				
				arcplot(G1,lwd.arcs=rev(Widths[[g]])+1,above=NULL,col.arcs=C,ljoin = 2,lend=2,pch.nodes=15,cex.nodes=2,cex.labels=1.2,col.nodes=rev(ColN[[g]]))
			}
			else{
				arcplot(G[[g]],lwd.arcs=Widths[[g]],above=NULL,col.arcs=Cols[[g]],ljoin = 2,lend=2,pch.nodes=15,cex.nodes=2,cex.labels=1.5,col.nodes=ColN[[g]])
			}
		}
		#dev.off()
	}
	#warnings()
	OutList=cbind(OutList,Category)
	#utils::write.table(OutList, paste(Name,".txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
	#rm(OutList)
	#gc() 
	output=list("ASInfo"=OutList,"Compositions"=Table1,"GroupTranscripts"=Table2,"NovelConnections"=NovelConnections)
	
	return(output)
	
}	
	

#' REIDS_JunctionAssessment
#' 
#' The REIDS_JunctionAssessment functions assess identified AS exons based on their 5'end and 3'end and exclusion junction support.
#' @export
#' @param Indices The .csv file created by Line_Indexer.py which contains indices for every gene in geneIDs.
#' @param DataFile The .csv file created by PivotTransformation. 
#' @param ASProbeSets The AS probe sets as identified by ASExons.
#' @param Juninfo A parameter specifying wether the annotations are user of Ensembl defined. If JunInfo is "User" the annotations provided in EandTrAnnot are used. If JunInfo is "Ensembl" the annotations in EandTrAnnot are used to set up tje junction associations but the gene name and position in transcriptData and positionData are used to connect with the Ensembl data base and retrieve corresponding information. 
#' @param JAnnotI The file name with line indices for the junction associations.
#' @param JAnnot The file name with the junction associations.
#' @param EandTrAnnotI The file name with line indices for the exon and isoform annotations.
#' @param EandTrAnnot The file name with the exon and isoform annotations.
#' @param PartiallyAnnotated Logical. Should the exon annotations with partially annotated probe sets still be included? If FALSE, these are excluded. If TRUE, these are included. Default is FALSE.
#' @param positionData The file with the chromosome start and ends for the probe sets. Only needed in JunInfo=Ensembl.
#' @param transcriptData The file with gene name of the transcripts. Only needed in JunInfo=Ensembl.
#' @param Groups A list with  elements speficifing the columns of the data in each group.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @param Low_GSamples A list with a  character vector per group containing the probe sets which are not DABG in that group.
#' @param Location A character string indication the place where the outputs are saved.
#' @param Name A character string with the name of the ouput file. Defaults to "REIDS_Jun".
#' @return The function returns four files. The first file has name "Name_ASInfo.txt" and contains a line per probe set. It shows the reached decisionregarding the probe set (Const/AS/not DABG),its linking and exclusion junctions, the fold change, the AS type and its annotated exons. The secondfile, "Name_Compositions.txt", is a list of all found transcripts for a particular TC ID. The third file,"Name_GroupTranscripts.txt" indicates whether a specific transcript is present or absent in a group. The fourth file "Name_NovelConnections.txt" contains junctions which are showing an undocumented connection between probe sets.
#' @examples
#' \dontrun{
#' data(TC1500264)
#' PivotTransformData(Data=TC1500264,GeneID=NULL,ExonID=NULL,
#' REMAPSplitFile="TC1500264_Gene_SplitFile.txt",Location="Output/",Name="TC1500264_Pivot")
#' REIDSFunction(ASPSR=c(), Indices="Output/TC1500264_LineIndex.csv",
#' DataFile="Output/TC1500264_Pivot.csv",nsim=50,informativeCalls=FALSE,Summarize=
#' c("WeightedAll","EqualAll"),rho=0.5,Low_AllSamples=c(),Groups=list(c(1:3),c(4:6)),
#' Location="Output",Name="TC1500264")
#' 
#' TC1500264_1vs2=ASExons(ExonScores="Output/TC1500264_ExonScores.txt",ArrayScores=
#' "Output/TC1500264_ArrayScores.txt",Exonthreshold=0.5,Groups=list(c(1:3),c(4:6)),paired=FALSE,
#' significancelevel=0.05)
#' 
#' ASPSR_PSR=TC1500264_1vs2[which(round(TC1500264_1vs2$ExonScore,2)>=0.50&
#' TC1500264_1vs2$adj.p.value<0.05),2]
#' 
#' REIDS_JunctionAssesment(Indices="Output/TC1500264_LineIndex.csv",DataFile=
#' "Output/TC1500264_Pivot.csv",ASProbeSets=ASPSR_PSR,Juninfo="User",JAnnotI=NULL,
#' JAnnot=NULL,EandTrAnnotI="Output/REMAP_Indices.txt",EandTrAnnot=
#' "Output/HJAY_REMAP.txt",positionData=NULL,transcriptData=NULL,Groups=
#' list(c(1:3),c(4:6)),Low_AllSamples=c(),Low_GSamples=c(),Location="Output",
#' Name="TC1500264")
#' }
REIDS_JunctionAssesment<-function(Indices,DataFile,ASProbeSets=c(),Juninfo="User",JAnnotI=NULL,JAnnot=NULL,EandTrAnnotI=NULL,EandTrAnnot=NULL,PartiallyAnnotated=FALSE,positionData=NULL,transcriptData=NULL,Groups=list(c(3,4,5),c(6,7,8)),Low_AllSamples=c(),Low_GSamples=c(),Location=NULL,Name="REIDS_Jun"){
	Lines=utils::read.table(Indices,header=TRUE,sep=",",stringsAsFactors=FALSE)
	Lines=data.frame(Lines)	
	Lines[,1]=as.numeric(Lines[,1])
	Lines[,2]=as.numeric(Lines[,2])
	
	if(!is.null(Location)){
		Files=list.files(path=Location,pattern=paste("^",Name,sep=""))	
		if(!paste(Name,"_ASInfo.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID","PSR_ID","Type","Unreliable Junctions","Linking Junctions","Exclusion Junctions","Supported by","Identified by","Fold Change","Exons")), file = paste(Location,"/",Name,"_ASInfo.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_Compositions.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID","TranscriptName","Junctions","ProbeSets")), file = paste(Location,"/",Name,"_Compositions.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_GroupTranscripts.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID","TranscriptName",paste("Present in Group",c(1:length(Groups))))), file = paste(Location,"/",Name,"_GroupTranscripts.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
	}
	ScoresOutput=apply(Lines,1, function(x) JunInfo(file_name=DataFile,file_pos=x[1],line_length=x[2],ASPSR=ASProbeSets,Juninfo,JAnnotI,JAnnot,EandTrAnnotI,EandTrAnnot,PartiallyAnnotated,positionData,transcriptData,Groups=Groups,Low_AllSamples,Low_GSamples,Plot=FALSE,Location,Name))
}

#' REIDS_IsoformAssesment
#' 
#' The REIDS_IsoformAssesment is an experimental function to analyze the isoform information based on the exon level values and the isoform composition.
#' @export
#' @param geneIDs A vector with the geneIDs to analyze.
#' @param IsoformInfo The path to the Composition file created by REIDS_JunctionAssesment or CreateOutput.
#' @param ExonLevel The path to the ExonLevel.txt file
#' @param Groups A list with  elements speficifing the columns of the data in each group.
#' @param paired Logical. Are the groups paired samples?
#' @param Location A character string indication the place where the outputs are saved.
#' @param Name A character string with the name of the ouput file. Defaults to "REIDSIsoforms".
#' @return The function returns three files. The first file has name "Name_IsoformIndication.txt" and contains an assesment of the relative expression levels of the isoforms. A second file,"Name_ExonTesting.txt", shows information regarding the differential expression of exon. A final file,"Name_PossibelDEIsoforms" lists isoforms which might be differentially expressed between the groups.
REIDS_IsoformAssesment<-function(geneIDs,IsoformInfo,ExonLevel,Groups,paired,Location=NULL,Name="REIDSIsoforms"){
	
	Info=utils::read.table(IsoformInfo,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	
	ExonLevel=utils::read.table(ExonLevel,sep="\t",header=TRUE,stringsAsFactors=FALSE,colClasses=c("character","character",rep("numeric",length(unlist(Groups)))))
	
	if(!is.null(Location)){
		Files=list.files(path=Location,pattern=paste("^",Name,sep=""))	
		if(!paste(Name,"_IsoformIndication.txt",sep="")%in%Files){
			mainname=c()
			for(g in 1:length(Groups)){
				mainname=c(paste("Probe Set Group ",g,sep=""),paste("Expression Group ",g,sep=""),mainname)
			}
			utils::write.table(t(c("TC_ID","Transcript",mainname)), file = paste(Location,"/",Name,"_IsoformIndication.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_ExonTesting.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID","Probe Set","Statistic","PValue","Adj.PValue")), file = paste(Location,"/",Name,"_ExonTesting.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
		if(!paste(Name,"_PossibleDEIsoforms.txt",sep="")%in%Files){
			utils::write.table(t(c("TC_ID","Percentage Isoforms DE","Isoforms")), file = paste(Location,"/",Name,"_PossibleDEIsoforms.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
		}
	}
	Isoformanalysis<-function(geneID,IsoformI,ExonData,groups,paired,Location){
		output=list()
		Comp=matrix(0,nrow=nrow(ExonData),ncol=nrow(IsoformI)/2)
		rownames(Comp)=unique(ExonData[,2])
		SQ=seq(1,nrow(IsoformI),2)
		colnames(Comp)=IsoformI[SQ,2]
		for(tr in 1:length(SQ)){
			S=strsplit(IsoformI[SQ[tr],4],"[|]")
			Comp[which(rownames(Comp)%in%S[[1]]),tr]=1
			S=strsplit(IsoformI[SQ[tr],3],"[|]")
			Comp[which(rownames(Comp)%in%S[[1]]),tr]=1
		}
		if(any(colSums(Comp)==0)){
			Comp=Comp[,-c(which(colSums(Comp)==0))]
		}
		
		if(length(unique(ExonData[,2]))==1){
			Comp=t(as.matrix(Comp))
			rownames(Comp)=unique(ExonData[,2])
		}
		if(ncol(Comp)==1){
			Set=which(Comp[,1]==0)
			if(length(Set)>0){
				if(any(rownames(ExonData)%in%names(Set))){
					V=ExonData[-c(which(rownames(ExonData)%in%names(Set))),,drop=FALSE]
				}	
				else{
					V=ExonData
				}
				Comp=Comp[-c(which(Comp[,1]==0)),,drop=FALSE]
			}
			else{
				V=ExonData
			}		
		}
		else{
			if(any(rowSums(Comp)==0)){
				Set=which(rowSums(Comp)==0)
				if(any(ExonData[,2]%in%names(Set))){
					V=ExonData[-c(which(ExonData[,2]%in%names(Set))),,drop=FALSE]
				}
				else{
					V=ExonData
				}
				Comp=Comp[-c(which(rowSums(Comp)==0)),,drop=FALSE]				
			}
			else{
				V=ExonData
			}
		}

		#unique probe sets
		if(ncol(Comp)==1){
			UPandJ=rownames(Comp)[which(Comp[,1]==1)]
			UP=UPandJ[which(substr(UPandJ,1,1)=="P")]
		}
		else{
			UPandJ=which(rowSums(Comp)==1)
			UP=names(UPandJ)[which(substr(names(UPandJ),1,1)=="P")]
		}
		
		
		Isoform=c()
		if(length(UP)>0){
			Trs=colnames(Comp)[sapply(UP,function(x) which(Comp[x,]==1))]
			IsLevel=c()
			for(g in 1:length(groups)){
				ILevel=sapply(UP,function(x) round(mean(as.vector(as.matrix(V[which(V[,2]==x),groups[[g]]+2]))),2))
				IsLevel=cbind(IsLevel,UP,ILevel)
			}
			Isoform=cbind(Trs,IsLevel)
			T=c()
			for(t in unique(Trs)){
				temp=c()
				for(g in 1:length(Groups)){
					temp=c(temp,paste(Isoform[which(Isoform[,1]==t),2],collapse="|"),mean(as.numeric(Isoform[which(Isoform[,1]==t),2*g+1])))
				}
				T=rbind(T,c(t,temp))
			}
			Isoform=T
		}
		else{
			Trs=c()
		}
		
		#tr lowest levels
		IsLevel=c()
		TrOthers=colnames(Comp)[which(!colnames(Comp)%in%Trs)]
		if(length(TrOthers)>0){	
			for(f in 1:length(TrOthers)){
				T=c()		
				JandPSR_Present=rownames(Comp)[which(Comp[,TrOthers[f]]==1)]
				PSR_Present=JandPSR_Present[which(substr(JandPSR_Present,1,1)=="P")]
				for(p in PSR_Present){
					M=c()
					for(g in 1:length(groups)){
						S=V[which(V[,2]==p),c(groups[[g]]+2)]
						M=c(M,p,round(mean(as.vector(as.matrix(S))),2))
					}	
					T=rbind(T,M)
					
				}
				
				T=as.matrix(T)
				rownames(T)=c(1:nrow(T))
				T=data.frame(T,stringsAsFactors=FALSE)
				T[,2]=as.numeric(T[,2])
				T[,4]=as.numeric(T[,4])
				
				R=c()
				for(g in 1:length(groups)){
					Select=which.min(T[,2*g])
					R=c(R,T[Select,2*g-1],T[Select,2*g])
				}
				IsLevel=rbind(IsLevel,c(TrOthers[f],R))
			}
			
			IsoformIndicationLevels=rbind(Isoform,IsLevel)
			
		}
		
		IsoformIndicationLevels=rbind(Isoform,IsLevel)
		if(class(IsoformIndicationLevels)!="matrix"){
			IsoformIndicationLevels=matrix(IsoformIndicationLevels,nrow=1)
		}
		rownames(IsoformIndicationLevels)=c(1:nrow(IsoformIndicationLevels))
		IsoformIndicationLevels=data.frame(IsoformIndicationLevels,stringsAsFactors=FALSE)
		for(g in 1:length(groups)){
			IsoformIndicationLevels[,2*g+1]=as.numeric(IsoformIndicationLevels[,2*g+1])
			IsoformIndicationLevels[,2*g+1]=as.numeric(IsoformIndicationLevels[,2*g+1])
		}
		colnames(IsoformIndicationLevels)=c("Transcript",paste(rep(c("Lowest Expr Probe Set","Group"),length(groups)),sort(rep(c(1:length(groups)),length(groups)))))
		IsoformIndicationLevels=cbind("GeneID"=geneID,IsoformIndicationLevels)
		#Indication not estimate
		
		ttest<-function(x,y,pairs){	
			out1=stats::t.test(x,y,paired=pairs)
			out2=cbind(out1$statistic,out1$p.value)
			return(out2)	
		}
		
		el=1
		output[[el]]=IsoformIndicationLevels
		names(output)[el]="IsoformIndication"
		el=el+1
		
		# DE assesment
		IsoformDETest=matrix(0,nrow=nrow(V),ncol=3)
		colnames(IsoformDETest)=c("statistic","p-value","adj.p-value")
		
		anovaFtest1<-function(data,Groups){
			fit<-stats::aov(as.vector(data)~Groups)
			out1=cbind(summary(fit)[[1]][1,4],summary(fit)[[1]][1,5])	
			return(out1)
			
		}
		
		if(length(groups)<=2){
			IsoformDETest[,c(1,2)] = t(sapply(as.vector(V[,2]),function(i) {ttest(x=V[which(V[,2]==i),groups[[1]]+2],y=V[which(V[,2]==i),groups[[2]]+2], pairs = paired)}))	
		}
		else{
			Groups=rep(0,length(unlist(groups)))
			for(j in 1:length(groups)){
				positions=groups[[j]]
				Groups[positions]=names(groups)[j]
			}
			Groups=as.numeric(Groups)
			Groups=factor(Groups)
			
			IsoformDETest[,c(1,2)]  = t(sapply(as.vector(V[,2]),function(i) {anovaFtest1(data=V[which(V[,2]==i),,drop=FALSE],Groups=Groups)}))				
		}
		
		p_vals=as.vector(as.matrix((IsoformDETest[,grep("^(p-value)",colnames(IsoformDETest))])))
		adj_p_vals=stats::p.adjust(p_vals,"BH")
		IsoformDETest[,3]=adj_p_vals
		
		IsoformDETest=as.data.frame(IsoformDETest)
		IsoformDE=cbind(geneID,V[,2],IsoformDETest)
		
		PossibleDETr=c()
		for(r in 1:nrow(IsoformDE)){
			if(IsoformDE[r,4]>=0.01){
				next
			}
			else{
				p=as.character(IsoformDE[r,2])
				PossibleDETr=c(PossibleDETr,unique(names(which(Comp[p,]==1))))
			}
		}
		PossibleDETr=unique(PossibleDETr)
		Percentage=length(PossibleDETr)/ncol(Comp)
		
		
		output[[el]]=IsoformDE
		names(output)[el]="ExonTesting"
		output[[el+1]]=c(geneID,Percentage,paste(PossibleDETr,collapse="|"))
		names(output)[el+1]="Possible DE Isoforms"
		
		
		if(!is.null(Location)){
			#save(output, file=paste(Location,"/REIDS_IsoformAssesment_", geneID, ".RData", sep=""))
			if(!is.null(output[[1]])){
				utils::write.table(output[[1]],file = paste(Location,"/",Name,"_IsoformIndication.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
			}
			if(!is.null(output[[2]])){
				utils::write.table(output[[2]],file = paste(Location,"/",Name,"_ExonTesting.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
				
			}
			if(!is.null(output[[3]])){
				utils::write.table(t(output[[3]]),file = paste(Location,"/",Name,"_PossibleDEIsoforms.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)				
			}
			
		}
		else{
			return(output)
		}
	
	}
	
	G=unique(Info[,1])
	
	invisible(lapply(G,function(x) Isoformanalysis(geneID=x,IsoformI=Info[which(Info[,1]==x),],ExonData=ExonLevel[which(ExonLevel[,1]==x),],groups=Groups,paired=paired,Location=Location)))

}


#' REIDS_Analysis
#' 
#' The REIDS_Analysis is a wrapper function for the REIDSFunction, the ASExon function, the REIDS_JunctionAssesment function and the REIDS_IsoformAssesment function.
#' @export
#' @param geneIDs A vector with the geneIDs to analyze.
#' @param Indices The .csv file created by Line_Indexer.py which contains indices for every gene.
#' @param DataFile The .csv file created by PivotTransformation. 
#' @param nsim The number of iterations to perform. Defaults to 1000.
#' @param informativeCalls Logical. Should the I/NI calls method be perform before applying the REIDS model?
#' @param Summarize A character vector specifying the wich summarization method to be performed. The choices are using "EqualAll", "WeightedAll", "EqualConst", "WeightedConst". The former two use all probe sets while the latter to use only the consituitive probe sets. Summarization on the consistuitive probe sets will only be performed if ASPSR is specified.
#' @param rho The threshold for filtering in the I/NI calls method. Probesets with scores higher than rho are kept.
#' @param Exonthreshold The exon score threshold to be maintained. If not NULL, probesets with an exon score lower than this value are not considered further and the p-values will be adjusted for multiplicity after testing. If NULL, all probesets are considered and a multiplicity correction is not performed.
#' @param Groups A list with elements specifying the columns of the data in each group.
#' @param paired Logical. Are the groups paired? If TRUE the mean paired differences are calculated and tested whether these are significantly different from zero or not.
#' @param significancelevel The significance level to be maintained on the p-values. The filtering on the significance is conducted only if an Exonthreshold is specified and the p-value are adjusted for multiplicity.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @param Low_GSamples A list with a  character vector per group containing the probe sets which are not DABG in that group.
#' @param Juninfo A parameter specifying wether the annotations are user of Ensembl defined. If JunInfo is "User" (default) the annotations provided in EandTrAnnot are used. If JunInfo is "Ensembl" the annotations in EandTrAnnot are used to set up tje junction associations but the gene name and position in transcriptData and positionData are used to connect with the Ensembl data base and retrieve corresponding information. 
#' @param JAnnotI The file name with line indices for the junction associations.
#' @param JAnnot The file name with the junction associations.
#' @param EandTrAnnotI The file name with line indices for the exon and isoform annotations.
#' @param EandTrAnnot The file name with the exon and isoform annotations.
#' @param PartiallyAnnotated Logical. Should the exon annotations with partially annotated probe sets still be included? If FALSE, these are excluded. If TRUE, these are included. Default is FALSE.
#' @param positionData The file with the chromosome start and ends for the probe sets. Only needed in JunInfo=Ensembl.
#' @param transcriptData The file with gene name of the transcripts. Only needed in JunInfo=Ensembl.
#' @param Location A character string indication the place where the outputs are saved. Defaults to Output.
#' @param Name A character string with the name of the ouput file. Defaults to "REIDSAnalysis".
#' @return The output will be written to each of the corresponding .txt files of the called upon functions.
#' @examples
#' \dontrun{
#' data(TC1500264)
#' PivotTransformData(Data=TC1500264,GeneID=NULL,ExonID=NULL,
#' REMAPSplitFile="TC1500264_Gene_SplitFile.txt",Location="Output/",Name="TC1500264_Pivot")
#'
#' REIDS_Analysis(Indices="Output/TC1500264_LineIndex.csv",DataFile="Output/TC1500264_Pivot.csv",
#' nsim=100,informativeCalls=FALSE,Summarize=c("WeightedAll","EqualAll","WeightedConst","EqualConst"),
#' rho=0.5,Exonthreshold=0.5,significancelevel=0.05,Groups=Groups,paired=FALSE,Low_AllSamples=c()
#' ,Low_GSamples=c(),Juninfo="User",JAnnotI=NULL,JAnnot=NULL,EandTrAnnotI="Output/REMAP_Indices.txt",
#' EandTrAnnot="Output/HJAY_REMAP.txt",positionData=NULL,transcriptData=NULL,
#' Location="OutputREIDSAnalysis",Name="TC1500264")
#' }
REIDS_Analysis<-function(geneIDs,Indices,DataFile,nsim=5000,informativeCalls=FALSE,Summarize=TRUE,rho=0.5,Exonthreshold=0.5,significancelevel=0.05,Groups,paired=FALSE,Low_AllSamples=c(),Low_GSamples=c(),
		Juninfo="User",JAnnotI=NULL,JAnnot=NULL,EandTrAnnotI=NULL,EandTrAnnot=NULL,PartiallyAnnotated=FALSE,positionData=NULL,transcriptData=NULL,Location="Output",Name="REIDSAnalysis"){
	
	#1) Perform the REIDS function
	print("Performing the REIDS function")
	REIDSFunction(ASPSR=c(),Indices,DataFile,nsim,informativeCalls,Summarize=Summarize[which(Summarize%in%c("WeightedAll","EqualAll"))],rho,Low_AllSamples,Groups,Location,Name)

	#2) AS Identification
	print("Performing the AS identification")
	ExonS=paste(Location,"/",Name, "_ExonScores.txt", sep="")
	ArrayS=paste(Location,"/",Name, "_ArrayScores.txt", sep="")
	ASTest=ASExons(ExonS,ArrayS,Exonthreshold=0,Groups,paired,significancelevel=NULL,Location=NULL)

	utils::write.table(t(c("TC_ID","PSR_ID","ExonScore","Statistic","Pvalue","Adj.PValue")), file = paste(Location,"/",Name,"_ASTesting.txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)		
	utils::write.table(ASTest,file=paste(Location,"/",Name, "_ASTesting.txt", sep=""), sep = "\t", row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
	ASPSR=unique(ASTest[which(round(ASTest[,3],2)>Exonthreshold&round(ASTest$adj.p.value,2)<significancelevel),2])
	
	print(paste(length(ASPSR)," found as AS probe sets",sep=""))
	
	#3) Summarize on Const probe sets
	if(any(Summarize%in%c("WeightedConst","EqualConst"))){
		print("Performing summarization on the constituitive probe sets")
		REIDSFunction(ASPSR,Indices,DataFile,nsim,informativeCalls,Summarize=Summarize[which(Summarize%in%c("WeightedConst","EqualConst"))],rho,Low_AllSamples,Groups,Location,Name)
	}
	#3) Junction information (if available)
	if(!is.null(EandTrAnnotI)){
		print("Performing junction and isoform assessment")
		REIDS_JunctionAssesment(Indices,DataFile,ASProbeSets=ASPSR,Juninfo,JAnnotI,JAnnot,EandTrAnnotI,EandTrAnnot,PartiallyAnnotated,positionData,transcriptData,Groups,Low_AllSamples,Low_GSamples,Location,Name)
		ExonD=paste(Location,"/",Name, "_ExonLevel.txt", sep="")
		IsinfoName=paste(Location,"/",Name,"_Compositions.txt",sep="")
		
		REIDS_IsoformAssesment(geneIDs,IsoformInfo=IsinfoName,ExonLevel=ExonD,Groups,paired,Location,Name)
	}	
}


## FIRMA Model Analysis

#' "FIRMAScores"
#' 
#' The FIRMAScores function performs an analyais on the FIRMA scores of the sanples. If the groups are not paired, the Firma all sample score will be the minimum value of group 1. A test statistic is perforrmed on the grouping of interest to see if there is a significant difference between them. If the data is paired, the mean paired differences are obtained and tested versus zero.
#' @export
#' @param Data The output of the FIRMA model
#' @param InformativeExons A character vector of exon IDs. As for the REIDS model probesets are filtered out by I/NI calls model and later on exon score, the remaining exons can be specified here. Only these shall be considered in the FIRMA analysis to make the results between REIDS and FIRMA more comparable
#' @param groups A list with two elements speficifing the columns of the data of group 1 in group1 and those of group 2 in group2. Here is the possibility to specifify multiple subgroups in group 1 for the calculation of the all FIRMA sample score.
#' @param paired Logical. Are the groups paired? If TRUE the mean paired differences are calculated and tested whether these are significantly different from zero or not.
#' @param significancelevel If specified, filtering is conducted on the p-values.
#' @return A data frame with one line per exon. The columns contain the gene ID, the exon ID, the test statistic, a p-value and an adjusted p-value. If the groups are paired also the mean paired difference is given. Only the probesets with high enough exon scores and a significant test are kept in the data frame.
#' @details The input to this function is the output of the FIRMA model. On the FIRMA scores, a t-test is performed. This is either between groups if the data is not paired and on the mean paired differences if the data is paired. If no pairing is present, an all sample FIRMA score is computed as the maximum of the minimum values of the 
#' of the subgroups in group 1. The returned p-value is adjusted for multiplicity. If a vector is given in InformativeExons, only these exon IDS are considered in the calculations.
#' @examples 
#' data(ExampleFirmaOutput)
#' FIRMATest=FIRMAScores(Data=ExampleFirmaOutput,InformativeExons=NULL,groups=list(group1=
#' list(group1a=c(1,2,3),group1b=NULL),group2=c(4,5,6)),paired=FALSE,significancelevel=NULL) 
FIRMAScores <- function(Data,InformativeExons=NULL,groups=list(group1=list(group1a=NULL,group1b=NULL),group2=NULL),paired=FALSE,significancelevel=NULL){
#	if(is.null(groups$group1)){
#		message("A test on the array scores will NOT be performed")
#	}
	
	if(!is.null(InformativeExons)){
		Data=Data[which(Data$ExonID%in%as.character(InformativeExons)),]
	}
	
	data=Data[,-c(1,2,3,4,5)]
	#Data=Data[-which(is.na(data)),]
	#data=data[-which(is.na(data)),]
	
	if(paired==FALSE){
		
		# 1) Compute the FIRMA scores as Sj=max(min(Fij_group1a),min(Fij_group1b))
		subgroupsgroup1=groups$group1
		if(is.null(subgroupsgroup1$group1b)){
			minofeachgroup=sapply(subgroupsgroup1[1],function(x) apply(data[,x],1,min))
		}
		else{
			minofeachgroup=sapply(subgroupsgroup1,function(x) apply(data[,x],1,min))
		}
	
		AllSample_FIRMAScore=apply(minofeachgroup,1,max)	
		
	
		Out1=data.frame('geneID'=Data[,1],'exonID'=Data[,2],'AllSample_FIRMAScore'=AllSample_FIRMAScore)
	
	
		# 2) Perform a t-test on the FIRMA scores
	
		data_group1=data[,unlist(groups$group1)]
		data_group2=data[,groups$group2]
	
		filterASExons1<-function(i,group1,group2){
			out1=stats::t.test(x=group1,y=group2)
			out2=cbind(out1$statistic,out1$p.value)
			return(out2)
		
		}
	
		Out_ttest=t(sapply(1:nrow(data),function(i) filterASExons1(i,group1=data_group1[i,],group2=data_group2[i,])))
	
		Out2=cbind(Out1,Out_ttest)	
		colnames(Out2)=c("GeneID","ExonID","AllSample_FIRMAScore","t.statistic","p.value")
		Out2$adj.p.value=stats::p.adjust(Out2$p.value,"fdr")
		Out2$GeneID=as.character(Out2$GeneID)
		Out2$ExonID=as.character(Out2$ExonID)
	}
	
	else if(paired==TRUE){
		
		filterASExons2<-function(Subset){
			out1=stats::t.test(x=Subset)
			out2=cbind(out1$statistic,out1$p.value)
			return(out2)
			
		}
		data_group1=data[,groups$group1]
		data_group2=data[,groups$group2]
		
		Paired_Diff=as.matrix(data_group1-data_group2)
		
		Mean_Diff=apply(Paired_Diff,1,mean)
		
		Out=data.frame('geneID'=Data[,1],'exonID'=Data[,2],'Mean_Diff'=Mean_Diff)
		Paired_Diff=Paired_Diff[-which(is.na(Out$Mean_Diff)),]
		Out=Out[-which(is.na(Out$Mean_Diff)),]
		
		Out_ttest=t(sapply(1:nrow(Paired_Diff),function(i) filterASExons2(Subset=Paired_Diff[i,])))
		
		Out2=cbind(Out,Out_ttest)
		
		colnames(Out2)=c("GeneID","ExonID","Mean_Diff","t.statistic","p.value")
		Out2$adj.p.value=stats::p.adjust(Out2$p.value,"fdr")
		Out2$GeneID=as.character(Out$GeneID)
		Out2$ExonID=as.character(Out$ExonID)

	}
	
	
	if(!is.null(significancelevel)){
		Out2$adj.p.value=stats::p.adjust(Out2$p.value,"fdr")
		Out2=Out2[which(Out2$adj.p.value<0.05),]
	}
	
	return(Out2)
	
}


#PLOTTING FUNCTION

#' "ExpressionLevelPlot"
#' 
#' The ExpressionLevelPlot produces a plot of the expression levels of a specific exon and its corresponding gene.
#' @export
#' @param GeneID The gene ID of the gene of interest.
#' @param ExonID The exon ID of the exon of interest.
#' @param Data The processed data as returned by DataProcessing. This is were the observed probe intensities will be retrieved.
#' @param GeneLevelData The gene level summarized data to retrieve the gene level values.
#' @param ExonLevelData The exon level summarized data to retrieve the exon level values.
#' @param Groups The groups of interest in the data.
#' @param ylabel The label for the y-axis.
#' @param title A title for the plot.
#' @examples
#' \dontrun{
#' data(TC12000010)
#' data(TC12000010_ExonLevel)
#' data(TC12000010_GeneLevel)
#' ExpressionLevelPlot(GeneID="TC12000010",ExonID="PSR12000150",
#' Data=TC12000010,GeneLevelData=TC12000010_GeneLevel,ExonLevelData
#' =TC12000010_ExonLevel,Groups=list(c(1:9),c(10:18)),ylabel="",
#' title="PSR12000150")
#' }

ExpressionLevelPlot<-function(GeneID=NULL,ExonID=NULL,Data,GeneLevelData=NULL,ExonLevelData=NULL,Groups,ylabel=NULL,title=""){
	message("The gene and exon level data will be log2 transformed")
	
	if(is.null(GeneID) | is.null(ExonID)){
		stop("no GeneID and/or ExonID specified")
	}
	
	groups=list()
	l=0
	for(g in 1:length(Groups)){
		groups[[g]]=l+seq_along(Groups[[g]])
		l=l+length((Groups[[g]]))
	}
	
	Exon_ExonID=ExonLevelData[which(ExonLevelData[,2]==ExonID),]
	Exon_ExonID=Exon_ExonID[,-c(1,2)]
	#Exon_ExonID=log2(Exon_ExonID)
	Exon_ExonID=Exon_ExonID[,unlist(Groups)]
	
	
	# Gene Level Data
	Gene_GeneID=GeneLevelData[which(GeneLevelData[,1]==GeneID),]
	#Gene_GeneID[,-c(1)]=log2(Gene_GeneID[,-c(1)])
	Gene_GeneID=Gene_GeneID[,-c(1)]
	Gene_GeneID=Gene_GeneID[,unlist(Groups)]

	
	ProbeIntensities=Data[which(Data[,2]==ExonID),]
	ProbeIntensities_ExonID=ProbeIntensities[,-c(1,2)]
	ProbeIntensities_ExonID=ProbeIntensities_ExonID[,unlist(Groups)]
	
	ProbeIntensities_ExonID=apply(ProbeIntensities_ExonID,2,as.character)
	ProbeIntensities_ExonID=apply(ProbeIntensities_ExonID,2,as.numeric)
	#plot 1 : gene level + exon level + probe intensities
	
	max1=max(Gene_GeneID,Exon_ExonID)
	max_ylim=max(as.numeric(ProbeIntensities_ExonID),max1)

	if(is.null(ylabel)){
		ylabel=paste("Transcript",GeneID," - PSR",ExonID,sep=" ")
	}
	
	gg=groups
	graphics::plot(0,0,typ="n",xlab="",ylab=ylabel,ylim=c(0,max_ylim+2),xlim=c(1,ncol(ProbeIntensities_ExonID)),xaxt="n",frame=TRUE,,cex.axis=2,cex.lab=2)
	for(g in 1:length(Groups)){
		graphics::lines(x=gg[[g]],y=Gene_GeneID[gg[[g]]],col="black",lwd=2) #Exon Level data
	}
	if(length(ExonLevelData)>0){
		for(g in 1:length(Groups)){
			graphics::lines(x=gg[[g]],y=Exon_ExonID[gg[[g]]],col="#f0d130",lwd=2) #Exon Level data
		}	
	}
	for(i in 1:nrow(ProbeIntensities_ExonID)){
		graphics::points(x=c(1:ncol(ProbeIntensities_ExonID)),y=as.matrix(ProbeIntensities_ExonID[i,]),pch=19,col="#f0d130")
	}
	graphics::axis(1,labels=colnames(Gene_GeneID),at=c(1:ncol(ProbeIntensities_ExonID)),las=2,cex.axis=2)	
	title(main = title,cex.main=2)
	
}

#' "trim"
#' 
#' The trim function trims a long character string on gene annotation to a short single gene annotation.
#' @param s Character vector to trim.
trim <- function(s) {
	s <- as.character(s);
	s <- sub("^[\t\n\f\r[:punct:]  ]*", "", s);
	s <- sub("[\t\n\f\r ]*$", "", s);
	s;
} 

#' "AnnotateGenes"
#' 
#' The AnnotateGenes function annotates the TC ID to a HGNC Gene symbol.
#' @param transcript.IDs The transcript IDs to annotate.
#' @param trinfo The transcript data frame from which to retrieve the annotation symbol.
AnnotateGenes <- function(transcript.IDs,trinfo){
	transcripts <- as.character(trinfo[trinfo[,1] %in% transcript.IDs,][1])
	gene.assignment <- strsplit(as.character(trinfo[trinfo[,1] %in% transcript.IDs,][,3]) , "//")
	gene.symbols <- unlist(lapply(gene.assignment, function(assignment){
						length <- length(assignment)
						if(length > 2){
							HGNC <- seq(from = 2, to = length, by = 5)
							assignment <- unique(trim(assignment[HGNC]))
							#assignment <- paste(assignment, collapse = "/")
						} else{
							assignment <- "---"} 
					}))
	if(length(gene.symbols)>0){
		gene.symbols[is.na(gene.symbols)]	<- "---"
		annotate      <- data.frame(transcript = transcripts, symbol = gene.symbols)
	}else{
		annotate=data.frame(transcript = transcripts, symbol = "---")
	}	
	annotate
}

#' "AnnotateGeneSymbol"
#' 
#' The AnnotateGeneSymbol function annotates HGNC Gene symbol to an ensemble region.
#' @param symbol.to.annotate Gene symbol to annotate to an ensembl region.
AnnotateGeneSymbol	<- function(symbol.to.annotate){
	ensembl = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
	attributes.region	<- c("chromosome_name", "start_position", "end_position", "ensembl_gene_id")							
	filter.symbol		  <- "hgnc_symbol"
	
	annotated.genes<-biomaRt::getBM(attributes = attributes.region, 
			filters = filter.symbol,
			values = symbol.to.annotate, 
			mart = ensembl)
	return(annotated.genes)
}

#' "GetIntensities"
#' 
#' The GetIntensities function retrieves the exon expression level in the data.
#' @param  trans The transcript ID. 
#' @param Data The data frame with the exon level expression of the transcript.
#' @param Groups A list with the groups (columns) of interest
GetIntensities <- function(trans, Data,Groups){

	Subset=Data[Data[,1]==trans,]
	if(any(is.na(Subset[,2]))){
		Subset=Subset[-c(which(is.na(Subset[,2]))),]
	}	
#	for(i in 3:ncol(Subset)){
#		Subset[,i]=log2(Subset[,i])			
#	}
	
	M_Subset=Subset[,-c(1,2)]
	
	mean.intensities <- lapply(Groups, function(group){
				if(length(group) > 1){
					apply(M_Subset[,group], 1, mean)
				} else {
					M_Subset[, group]
				}
			})
	mean.intensities <- t(do.call(rbind, mean.intensities))
	gene.intensities <- as.data.frame(cbind(Subset, mean.intensities))
	
	return(gene.intensities)
}


#' "TranscriptsPlot "
#' 
#' The TranscriptsPlot function plots the known Ensemble transcript isoform composition plots.
#' @export 
#' @param trans The TC ID of the transcript to be plotted.
#' @param positions A table with the start and stop positions of the probe sets.
#' @param transcriptinfo A table with the transcript information of the TC ID.
#' @param display.probesets Logical. Should the probe sets be shown?
#' @param Data The exon level summarized data.
#' @param Groups A list with the groups (columns) of interest in the data.
#' @param Start Specify a specific start point on in the genome.
#' @param Stop Specify a specific stop point on in the genome.
#' @param Highlight A character string specifying a probe set to be highlighted in the transcript composition.
#' @examples 
#' \dontrun{
#' data(positions_36)
#' data(transcript.clusters.NetAffx.36)
#' data(TC12000010_ExonLevel)
#' TranscriptsPlot(trans="TC12000010", display.probesets = TRUE,Data=TC12000010_ExonLevel,
#' Groups=list(c(10:18),c(19:27)),Start=NULL,Stop=NULL,Highlight="PSR12000150")
#' } 
TranscriptsPlot <- function(trans, positions, transcriptinfo, display.probesets = TRUE,Data,Groups=list(),Start=NULL,Stop=NULL,Highlight=NULL){
	ensembl = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
	attributes.region	<- c("chromosome_name", "start_position", "end_position", "ensembl_gene_id")							
	filter.symbol		  <- "hgnc_symbol"
	
	trans1=paste(trans,".hg",sep="")
	exons=positions[positions$transcript_cluster_id == trans1, "probeset_id"]
	
	PSR=sapply(exons, function(x) substr(x,1,3)!="JUC")
	exons=exons[PSR]
	trans2=paste(trans1,".1",sep="")
	symbol.to.annotate <- AnnotateGenes(transcript.IDs=trans2,trinfo=transcriptinfo)$symbol # get HBC identity
	ensembl.output     <- AnnotateGeneSymbol(symbol.to.annotate) # get region
	print(paste("Gene position: ",paste(ensembl.output,collapse=" "),sep=""))
	if(nrow(ensembl.output) > 0){
		strand <- transcriptinfo[transcriptinfo[,1] == trans2,][2]
		title = GenomeGraphs::makeTitle(text = paste(ensembl.output[4], 
						" (", strand,")", sep = ""), 
				color ='darkred')
		gene  = GenomeGraphs::makeGene(id = ensembl.output$ensembl_gene_id, biomart = ensembl)
		transcript = GenomeGraphs::makeTranscript(id = ensembl.output$ensembl_gene_id, biomart = ensembl)
		
		gene.positions   <- positions[positions$probeset_id %in% exons,]
		gene.positions$start=as.numeric(gene.positions$start)
		gene.positions$stop=as.numeric(gene.positions$stop)
		gene.positions<-gene.positions[order(gene.positions$start,decreasing=FALSE),]
		gene.positions[,1]=sapply(gene.positions[,1],function(x) strsplit(x,"[.]")[[1]][1])
		print(paste("Min Exon position: ",min(as.numeric(gene.positions$start)), " ; Max Exon position ", max(as.numeric(gene.positions$stop)),sep=""))
		gene.n.probes	 <- rep(1, nrow(gene.positions))
		
		gene.intensities <- GetIntensities(trans, Data,Groups)
		
		if(nrow(gene.intensities)>nrow(gene.positions)){
			gene.intensities<-gene.intensities[which(gene.intensities[,2]%in%gene.positions[,1]),]
		}
		else{
			gene.positions<-gene.positions[which(gene.positions[,1]%in%gene.intensities[,2]),]
		}
		gene.intensities=gene.intensities[match(gene.positions[,1],gene.intensities[,2]),]
		gene.intensities=as.matrix(gene.intensities[,-c(1,2)])
		#** user-defined colour
		if(length(Groups) > 2){
			palette <- grDevices::rainbow(length(Groups))
		} else {
			palette <- c("red", "blue")
		}
		colour<- NULL
		for(i in 1:length(Groups)){
			colour[Groups[[i]]]   <- palette[i]
		}
		
		lwd<- rep(0.05, ncol(gene.intensities)-2)
		
		
		#** user-defined colour
		colour[(length(unlist(Groups)) + 1):(length(unlist(Groups)) + length(Groups))] <- palette
		
		#** user-defined line width
		lwd[(length(Groups) + 1):(length(Groups) + length(Groups))]  <- 3
		
		Names=sapply(gene.positions[,1],function(x) strsplit(x,"[.]")[[1]][1])
		

		exon  = GenomeGraphs::makeExonArray(intensity = gene.intensities, 
				probeStart = as.numeric(gene.positions[,3]), 
				probeEnd   = as.numeric(gene.positions[,4]),
				probeId    = Names, 
				nProbes    = gene.n.probes,
				dp         = DisplayPars(color = colour, lwd = lwd, 
						mapColor ='dodgerblue2'), 
				displayProbesets =display.probesets )

		
		if(!is.null(Highlight)){
			R=list()
			for(p in 1:length(Highlight)){
				gene.pos<-gene.positions[gene.positions[,1]==Highlight[p],]
				print(paste("Highlighted Region: ",paste(gene.pos,collapse=" ")))
				rOverlay <- makeRectangleOverlay(start = (as.numeric(gene.pos$start)-500), end = (as.numeric(gene.pos$stop)+500),region=c(2,3),dp = DisplayPars(alpha = .5, fill = "pink"))
				R[[p]]=rOverlay
			}
			if(is.null(Start)&is.null(Stop)){
				gdPlot(list(exon,gene,transcript),overlays = R)
			}
			else{
				gdPlot(list(exon,gene,transcript),Start,Stop,overlays = R)
			}
		}
		else{
			
			if(is.null(Start)&is.null(Stop)){
				gdPlot(list(exon,gene,transcript))
			}
			else{
				gdPlot(list(exon,gene,transcript),Start,Stop)
			}
		}
	}
}

# SI Index
#' "SpliceIndex"
#' 
#' The SpliceIndex function computes the ratio of the splice indices of defined groups. Further, it performs the SI algorithm if the length of the groups is two and the MiDAS algorithm is more groups are specified. Both algorithms are implemented as defined by Affymetrix.
#' @export
#' @param GeneData The microarray data summarized at gene level.
#' @param ExonData The microarray data summarized at exon level.
#' @param InformativeExons A character vector of exon IDs. As for the REIDS model probesets are filtered out by I/NI calls model and later on exon score, the remaining exons can be specified here. Only these shall be considered in the FIRMA analysis to make the results between REIDS and FIRMA more comparable
#' @param Groups A list with the groups (columns) of interest in the data.
#' @param paired Logical. Are the groups paired? only used if two groups are present.
#' @param significancelevel If specified, filtering is conducted on the p-values.
#' @return A data frame wiith one line per exon. The columns conatin the gene ID, the exon ID, the ratio of the splice indices if two groups are present, a t- or F-statitic, a p-value and an adjusted p-value.
#' @details Given the gene level and exon level summarized data, the splice index method for the detection of alternative splicing is performed. The first step is to normalize the exon
#' data by taking the ratio with the gene level data. These values are referred to as the splice indices. If only two groups are specified, the ratio of their splice indices is taken as a measure for alternative splicing. The more the ratio deviates from zero, the more there is an indication of alternative splicing.
#' A t-test is conducted on the splice indices of the two groups to test their difference. If more than two groups are specified, an ANOVA model is fitted on the splice indices to discover with an F-test whether there is a difference between the groups somewhere. If a vector of informative exons is given 
#' to the function, only these are considered for the analysis. Finally, the p-values are adjusted for multiplicity and if a significance level is specified only the significant p-valuesare kept in the data frame.
#' @examples
#' data(TC12000010_ExonLevel)
#' data(TC12000010_GeneLevel)
#' SI_Test=SpliceIndex(GeneData=TC12000010_GeneLevel,ExonData=TC12000010_ExonLevel
#' ,InformativeExons=NULL,Groups=list(group1=c(1:9),group2=c(10:18)),
#' paired=FALSE,significancelevel=NULL)
SpliceIndex<-function(GeneData,ExonData,InformativeExons=NULL,Groups=list(group1=NULL,group2=NULL),paired=FALSE,significancelevel=NULL){

	
	if(!is.null(InformativeExons)){
		ExonData=ExonData[which(ExonData$ExonID%in%as.character(InformativeExons)),]
	}

	GeneID=intersect(GeneData[,1],ExonData[,1])  #special measure for the HTA Data : ExonLevel has Gene IDs from the mapping while GeneLevel has those of the ENSg cdf: different genes
	
	
	si<-function(ID,DG,DE,groups,paired){
		DataG=as.matrix(DG[,-c(1)])  #Gene Expression Data Level
		DataE=as.matrix(DE[,-c(1,2)]) #Exon Expression Data Level (exons coming from the same gene)
		
		DataE=DataE[!duplicated(rownames(DataE)), ,drop=FALSE]
		
		grouping=c(groups[[1]],groups[[2]])
		
		DataE=DataE[,grouping,drop=FALSE]
		DataG=DataG[,grouping,drop=FALSE]
		
		group1=seq_along(groups[[1]])
		group2=length(group1)+seq_along(groups[[2]])
		groups=list(group1=group1,group2=group2)
		
		NormData=sweep(DataE,2,DataG,"-")  #The Normalized Expression Data ( "-" because  log2 is already taken)
		
		SI_group1=apply(NormData[,group1,drop=FALSE],1,mean)
		
		SI_group2=apply(NormData[,group2,drop=FALSE],1,mean)
		
		SI=as.matrix(SI_group1-SI_group2,drop=FALSE)
		colnames(SI)="Splice Index"
		rownames(SI)=DE$ExonID
		
		stats<-function(i,Data,groups,paired){
			out1=stats::t.test(Data[i,groups$group1],y=Data[i,groups$group2],paired=paired)
			out2=cbind(out1$statistic,out1$p.value)
			
			return(out2)
		}
		
		tstats=t(sapply(c(1:nrow(NormData)),stats, Data=NormData, groups=groups,paired=paired)) 
		colnames(tstats)=c("t-statistic","p.value")
		
		out_temp=cbind(SI,tstats)
		
		return(out_temp)
	}
	

	midas<-function(ID,DG,DE,groups){
		DataG=as.matrix(DG[,-c(1)])  #Gene Expression Data Level
		DataE=as.matrix(DE[,-c(1,2)]) #Exon Expression Data Level (exons coming from the same gene)
		
		DataE=DataE[!duplicated(rownames(DataE)), ,drop=FALSE]
		
		grouping=as.vector(unlist(groups))
		
		DataE=DataE[,grouping,drop=FALSE]
		DataG=DataG[,grouping,drop=FALSE]
	
		
		groupstemp=list()
		for(i in 1:length(groups)){
			groupstemp[[i]]=rep(i,length(groups[[i]]))		
		}
		names(groupstemp)=names(groups) #in the assumptution that the names of groups are group1, group2, group3,...
		groups=as.factor(unlist(groupstemp))
		
		#normalized intensities : responses of the ANOVA model
		NormData=sweep(DataE,2,DataG,"-") 
	
		
		#anova model
		stats=t(sapply(1:nrow(DataE),function(i) summary(stats::aov(NormData[i,]~groups))[[1]][1,c(4,5)]))
		colnames(stats)=c("F-statistic","p.value")
		
		out_temp=stats
		return(out_temp)
		
	}
	
	
	if(length(Groups)==2){
		#2 groups ; perform splice index analysis with t-test
		out=lapply(GeneID,function(x) si(ID=x,DG=GeneData[which(GeneData[,1]==x),],DE=ExonData[which(ExonData[,1]==x),],groups=Groups,paired=paired)) #for every gene
		names(out)=GeneID
	}
	
	if(length(Groups)>2){
		out=lapply(GeneID,function(x) midas(ID=x,DG=GeneData[which(GeneData[,1]==x),],DE=ExonData[which(ExonData[,1]==x),],groups=Groups)) #for every gene
		names(out)=GeneID
	}
	
	out1<-do.call(rbind.data.frame, out)
	out1$adj.p.value=stats::p.adjust(out1$p.value,"fdr")
	
	fac=as.factor(sapply(c(1:nrow(out1)),function(x) strsplit(rownames(out1)[x],"[.]")[[1]][1]))
	out2<-split(out1,fac)
	
	replacerownames<-function(ID,names,Data){
		if(any(is.na(names))){
			Data=Data[-c(which(is.na(names))),]
			names=names[-c(which(is.na(names)))]		
		}
		rownames(Data)=names
		return(Data)
	}
	
	output=lapply(1:length(GeneID),function(x) replacerownames(ID=x,names=ExonData[which(ExonData[,1]==GeneID[x]),2],Data=out2[[which(names(out2)==GeneID[x])]])) 
	names(output)=GeneID
	
	output2=do.call(rbind.data.frame, output)
	
	Names=strsplit(rownames(output2),"[.]")
	Names=do.call(rbind.data.frame, output)
	colnames(Names)=c("GeneID","ExonID")
	Names$GeneID=as.character(Names[,1])
	Names$ExonD=as.character(Names[,1])
	
	SI=cbind(GeneID=Names$GeneID,ExonID=Names$ExonID,output2)
	rownames(SI)=seq(1:nrow(SI))
	
	if(!is.null(significancelevel)){
		SI$adj.p.value=stats::p.adjust(SI$p.value,"fdr")
		SI=SI[which(SI$adj.p.value<0.05),]
	}

	return(SI)
	
}

#' @title Graph Information
#' 
#' @description Gets graph information
#' 
#' @param edgelist basically a two-column matrix with edges 
#' @param vertices optional vector of vertex names corresponding with 
#' those in the edgelist
#' @param sorted logical to indicate if nodes should be sorted 
#' (default \code{FALSE})
#' @param decreasing logical to indicate type of sorting 
#' (used only when \code{sorted=TRUE})
#' @param ordering optional numeric or string vector providing the 
#' ordering of nodes. When provided, this parameter overrides 
#' \code{sorted=TRUE}). See the details section for more information.
#' @param labels optional string vector with labels for the nodes
#' @keywords internal
graph_info <- 
		function(edgelist, vertices, sorted = FALSE, decreasing = FALSE, 
				ordering = NULL, labels = NULL)
{
	# ======================================================
	# Checking arguments
	# ======================================================
	# edgelist as a two-column matrix
	if (!is.matrix(edgelist) || ncol(edgelist) != 2)
		stop("\nSorry, 'edgelist' must be a two column matrix")
	
	num_edges = nrow(edgelist)
	# get nodes (this could be numeric or character)
	if(hasArg(vertices)){
		#to deal with singleton nodes
		nodes = vertices 
	}else{
		nodes = unique(as.vector(t(edgelist)))  
	}
	num_nodes = length(nodes)
	# check labels (i.e. node names)
	if (!is.null(labels))
	{
		if (length(labels) != num_nodes)
			stop("\nLength of 'labels' differs from number of nodes")
	} else {
		labels = nodes
	}
	
	# auxiliar order (this may change if sorted or ordering required)
	aux_ord = 1:num_nodes  
	
	# If sorted is required, ennumerate nodes
	if (sorted) {
		ordered_nodes = order(nodes, decreasing = decreasing)
		nodes = nodes[ordered_nodes]
		labels = labels[ordered_nodes]
		# auxiliar order
		aux_ord = ordered_nodes
	}
	
	# If ordering is provided, re-ennumerate nodes
	if (!is.null(ordering)) 
	{
		if (length(ordering) != num_nodes) {
			stop("\nLength of 'ordering' differs from number of nodes")      
		}
		
		if (is.character(ordering)) {
			# make sure labels contains elements in ordering
			unmatched_ordering <- !(ordering %in% labels)
			if (any(unmatched_ordering)) {
				undetected = ordering[unmatched_ordering]
				stop(sprintf("\nUnrecognized values in ordering: '%s'", undetected))
			}
			ordering = match(ordering, labels)
		}
		
		nodes = nodes[ordering]
		labels = labels[ordering]
		# auxiliar order
		aux_ord = ordering
	}
	
	## output
	list(
			nodes = nodes,
			labels = labels,
			num_nodes = num_nodes,
			num_edges = num_edges,
			aux_ord = aux_ord
	)
}


#' @title X or Y coordinates of node locations
#' 
#' @description
#' Gives axis locations of each node
#'  
#' @param num_nodes number of nodes
#' @param aux_ord vector with the index number for ordering the nodes
#' @param labels optional string vector with labels for the nodes
xynodes <- function(num_nodes, aux_ord, labels)
{
	# ======================================================
	# Coordinates of nodes (i.e. vertices)
	# ======================================================
	# node labels at equal distances from each other
	nf = rep(1 / num_nodes, num_nodes)
	# center coordinates of node labels
	fin = cumsum(nf)
	ini = c(0, cumsum(nf)[-num_nodes])
	centers = (ini + fin) / 2
	names(centers) = labels[aux_ord]
	
	# output
	centers
}


#' @title Arc Radius Locations
#' 
#' @description Computes the location and radius of each arc
#' 
#' @param edgelist 2-column matrix
#' @param nodes vector of nodes
#' @param centers vector with xy-positions of nodes
#' @return a list with locations and radios
#' @return \item{locs}{locations}
#' @return \item{radios}{radius values}
#' @keywords internal
arc_radius_locs <- function(edgelist, nodes, centers)
{
	# ======================================================
	# Coordinates of arcs (i.e. edges)
	# ======================================================
	# handy matrix with numeric indices '1:FROM' , '2:TO'
	edges_from_to = matrix(0, nrow(edgelist), 2)
	for (i in 1L:nrow(edgelist))
	{
		edges_from_to[i,1] = centers[which(nodes == edgelist[i,1])]
		edges_from_to[i,2] = centers[which(nodes == edgelist[i,2])]
	}
	
	# maximum radius of arcs 
	radios = abs(edges_from_to[,1] - edges_from_to[,2]) / 2
	max_radios = which(radios == max(radios))
	max_rad = unique(radios[max_radios] / 2)    
	
	# arc locations
	locs = rowSums(edges_from_to) / 2
	
	# output
	list(locs = locs, radios = radios)
}


#' @title Above or Below
#' 
#' @description Determines how arcs should be displayed
#' @details 
#' If \code{horizontal = TRUE} then arcs can be plotted above or 
#' below the horizontal axis \cr
#' If \code{horizontal = FALSE} then arcs can be plotted to the right or 
#' left of the vertical axis
#' @param edgelist two-column matrix
#' @param above optional numeric or logical vector indicating what edges 
#' (arcs) should be plotted above (or to the right of) of chosen axis
#' If \code{above = NULL} then all arcs are plotted above (or to the right) 
#' If \code{above} is numeric, it cannot contain both positive and negative
#' indices.
#' If \code{above} is logical, its length must equal the number of rows in
#' \code{edgelist}
#' @return a logical vector indicating how arcs should be displayed
#' @keywords internal
above_below <- function(edgelist, above)
{
	# ======================================================
	# Coordinates of arcs (i.e. edges) below the axis
	# ======================================================
	# check above
	if (is.null(above)) {
		above = rep(TRUE, nrow(edgelist))     
	} else {
		if (length(above) > nrow(edgelist))
			stop("\nlength of 'above' exceeds number of rows in 'edgelist'")
		# check numeric above and convert to logical
		if (is.numeric(above)) {
			above_positive <- any(above > 0)
			above_negative <- any(above < 0)
			if (above_positive & above_negative)
				stop("\n'above' cannot contain both negative and positive indices")
			# convert to logical
			if (all(above > 0)) {
				above = 1:nrow(edgelist) %in% above
			}
			if (all(above < 0)) {
				above <- !(-(1:nrow(edgelist)) %in% above)
			}
			if (all(above == 0)) {
				above = rep(FALSE, nrow(edgelist))          
			}
		}
		# check logical above
		if (is.logical(above)) {
			if (length(above) != nrow(edgelist))
				stop("\nlength of 'above' must equal number of rows in 'edgelist'")
		}
	}
	
	# output
	above
}


#' @title Minimum and Maximum Margin Limits
#' @description Computes the minimum and maximum margin limits of 
#' plotting region
#' @param radios vector of arc radius
#' @param above logical vectors indicating whether arcs should be displayed
#' @return list with minimum and maximum margin limits
#' @keywords internal
min_max_margin <- function(radios, above)
{
	# determine maximum radius
	max_radios = which(radios == max(radios))
	
	# minimum and maximum margin limits
	lim_min = 0
	lim_max = 0
	
	above_radios = radios[above]
	if (length(above_radios > 0)) {
		max_above_radios = which(above_radios == max(above_radios))[1]
		lim_max = above_radios[max_above_radios]
	}
	
	below_radios = radios[!above]
	if (length(below_radios > 0)) {
		max_below_radios = which(below_radios == max(below_radios))[1]
		lim_min = -1 * below_radios[max_below_radios]
	}
	
	# margin limits
	list(min = lim_min, max = lim_max)
}




#' @title Arc Diagram Plot
#' 
#' @description
#' Give me an edgelist and I'll help you plot a pretty damn arc diagram
#' 
#' @details
#' The arcs are scaled such that they fit in a plot region with its
#' x-axis ranging from zero to one. Node symbols and labels can be
#' optionally displayed. Node symbols are displayed through
#' the function \code{points}. In turn, node labels are displayed
#' through the function \code{mtext}.
#' 
#' When \code{ordering} is provided in numeric format and node labels are 
#' strings, the labels are alphabetically ordered first, and then nodes are 
#' sorted according to the provided \code{ordering}.
#' 
#' If \code{ordering} is provided in string format, the node labels must be 
#' strings as well. The nodes will be sorted according to \code{ordering}.
#' 
#' @param edgelist basically a two-column matrix with edges 
#' @param vertices optional vector of vertex names corresponding with 
#' those in the edgelist
#' @param sorted logical to indicate if nodes should be sorted 
#' (default \code{FALSE})
#' @param decreasing logical to indicate type of sorting 
#' (used only when \code{sorted=TRUE})
#' @param ordering optional numeric or string vector providing the 
#' ordering of nodes. When provided, this parameter overrides 
#' \code{sorted=TRUE}). See the details section for more information.
#' @param labels optional string vector with labels for the nodes
#' @param horizontal logical indicating whether to plot 
#' in horizontal orientation
#' @param above optional vector indicating which arcs should be displayed
#' above (or to the right) and below (or to the left) of the axis
#' @param col.arcs color for the arcs (default \code{"gray50"})
#' @param lwd.arcs line width for the arcs (default 1)
#' @param lty.arcs line type for the arcs (see \code{\link{par}})
#' @param lend the line end style for the arcs (see \code{\link{par}})
#' @param ljoin the line join style for the arcs (see \code{\link{par}})
#' @param lmitre the line mitre limit for the arcs (see \code{\link{par}})
#' @param show.nodes logical indicating whether to show node symbols
#' @param pch.nodes plotting 'character', i.e. symbol to use when
#' plotting nodes (\code{pch.nodes=0:25})
#' @param cex.nodes expansion of the node symbols (default 1)
#' @param col.nodes color of the node symbols (default \code{"gray50"})
#' @param bg.nodes background (fill) color for the node symbols 
#' given by \code{pch.nodes=21:25}
#' @param lwd.nodes line width for drawing node symbols 
#' (see \code{\link{points}})
#' @param show.labels logical indicating whether to show node labels
#' @param col.labels color of the node labels (default \code{"gray50"})
#' @param cex.labels expansion of node labels (default \code{"gray50"})
#' @param las numeric in {0,1,2,3}; the style of axis labels 
#' (see \code{\link{par}})
#' @param font font used for node labels (see \code{\link{par}})
#' @param line on which margin line the node labels are displayed, 
#' starting at 0 counting outwards (see \code{\link{mtext}})
#' @param outer use outer margins, if available, to plot node labels
#' (see \code{\link{mtext}})
#' @param adj adjustment for each string in reading direction 
#' (see \code{\link{mtext}})
#' @param padj adjustment for each string perpendicular to 
#' the reading direction (see \code{\link{mtext}})
#' @param axes logical indicating whether to plot the axes 
#' (default \code{FALSE})
#' @param ... further graphical parameters (see \code{\link{par}}), including
#' \code{family}, \code{xpd}, \code{main}, \code{asp}, etc.
#' @author Gaston Sanchez
arcplot <- function(
		edgelist, vertices, sorted = FALSE, decreasing = FALSE, ordering = NULL, 
		labels = NULL, horizontal = TRUE, above = NULL, 
		col.arcs = "#5998ff77", lwd.arcs = 1.8, lty.arcs = 1, 
		lend = 1, ljoin = 2, lmitre = 1, show.nodes = TRUE, pch.nodes = 19, 
		cex.nodes = 1, col.nodes = "gray80", bg.nodes = "gray80", lwd.nodes = 1,
		show.labels = TRUE, col.labels = "gray55",
		cex.labels = 0.9, las = 2, font = 1, line = 0, 
		outer = FALSE, adj = NA, padj = NA, axes = FALSE, ...)
{
	# Get graph information
	if (hasArg(vertices)) { 
		nodes_edges = graph_info(edgelist, vertices = vertices, sorted = sorted, 
				decreasing = decreasing, 
				ordering = ordering, labels = labels)
	} else {
		nodes_edges = graph_info(edgelist, sorted = sorted, decreasing = decreasing, 
				ordering = ordering, labels = labels)
	}
	nodes = nodes_edges$nodes
	num_nodes = nodes_edges$num_nodes
	num_edges = nodes_edges$num_edges
	aux_ord = nodes_edges$aux_ord
	labels = nodes_edges$labels
	
	# x-y node coordinates
	centers = xynodes(num_nodes, aux_ord, labels)
	
	# determine above or below display of arcs
	above = above_below(edgelist, above)
	
	# arc radius and locations
	radios_locs = arc_radius_locs(edgelist, nodes, centers)
	radios = radios_locs$radios
	locs = radios_locs$locs
	
	# ======================================================
	# Graphical parameters for Arcs
	# ======================================================
	# color of arcs
	if (length(col.arcs) != num_edges) 
		col.arcs = rep(col.arcs, length=num_edges)
	# line width of arcs
	if (length(lwd.arcs) != num_edges) 
		lwd.arcs = rep(lwd.arcs, length=num_edges)
	# line type of arcs
	if (length(lty.arcs) != num_edges) 
		lty.arcs = rep(lty.arcs, length=num_edges)
	
	# ======================================================
	# Graphical parameters for Nodes
	# ======================================================
	# pch symbol of nodes
	if (length(pch.nodes) != num_nodes) {
		pch.nodes = rep(pch.nodes, length = num_nodes)    
	}
	pch.nodes = pch.nodes[aux_ord]
	# cex of nodes
	if (length(cex.nodes) != num_nodes) {
		cex.nodes = rep(cex.nodes, length = num_nodes)    
	}
	cex.nodes = cex.nodes[aux_ord]
	# color of nodes
	if (length(col.nodes) != num_nodes) {
		col.nodes = rep(col.nodes, length = num_nodes)    
	}
	col.nodes = col.nodes[aux_ord]
	# bg of nodes
	if (length(bg.nodes) != num_nodes) {
		bg.nodes = rep(bg.nodes, length = num_nodes)    
	}
	bg.nodes = bg.nodes[aux_ord]
	# line widths of nodes
	if (length(lwd.nodes) != num_nodes) {
		lwd.nodes = rep(lwd.nodes, length = num_nodes)    
	}
	lwd.nodes = lwd.nodes[aux_ord]
	
	# ======================================================
	# Graphical parameters for Node Labels
	# ======================================================
	# color of labels
	if (length(col.labels) != num_nodes) {
		col.labels = rep(col.labels, length = num_nodes)    
	} 
	col.labels = col.labels[aux_ord]
	# cex of labels
	if (length(cex.labels) != num_nodes) {
		cex.labels = rep(cex.labels, length = num_nodes)    
	}
	cex.labels = cex.labels[aux_ord]
	
	# ======================================================
	# Plot arc diagram (horizontally or vertically)
	# ======================================================
	# auxiliar vector for plotting arcs
	z = seq(0, pi, length.out = 100)
	
	if (horizontal) {
		side = 1 
	} else {
		side = 2
	}
	xlim=NULL
	ylim=NULL
	if (is.null(xlim)) {
		if (horizontal) {
			xlim = c(-0.015, 1.015)
			x_nodes = centers
		} else {
			xlims = min_max_margin(radios, above)
			xlim = c(xlims$min, xlims$max)
			x_nodes = rep(0, num_nodes)
		}
	} else {
		if (horizontal) {
			x_nodes = centers
		} else {
			x_nodes = rep(0, num_nodes)
		}
	}
	
	if (is.null(ylim)) {
		if (horizontal) {
			ylims = min_max_margin(radios, above)
			ylim = c(ylims$min, ylims$max)
			y_nodes = rep(0, num_nodes) 
		} else {
			ylim = c(-0.015, 1.015)
			y_nodes = centers
		}
	} else {
		if (horizontal) {
			y_nodes = rep(0, num_nodes) 
		} else {
			y_nodes = centers
		}
	}
	
	# open empty plot window
	graphics::plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", 
			xlab = "", ylab = "", axes = axes, ...)
	# add each edge
	for (i in 1L:num_edges)
	{
		# get radius length
		radio = radios[i]
		if (horizontal) {
			# x-y coords of each arc
			x_arc = locs[i] + radio * cos(z)
			if (above[i]) { # above axis
				y_arc = radio * sin(z)
			} else {  # below axis
				y_arc = radio * sin(-z)
			}
		} else {
			# x-y coords of each arc
			y_arc = locs[i] + radio * cos(z)
			if (above[i]) { # above axis
				x_arc = radio * sin(z)
			} else {  # below axis
				x_arc = radio * sin(-z)
			}      
		}
		
		# plot arc connecting nodes
		graphics::lines(x_arc, y_arc, col=col.arcs[i], lwd=lwd.arcs[i], lty=lty.arcs[i],
				lend=lend, ljoin=ljoin, lmitre=lmitre)
		# add node symbols with points
		if (show.nodes) {
			graphics::points(x=x_nodes, y=y_nodes, pch=pch.nodes, 
					col=col.nodes, bg=bg.nodes, cex=cex.nodes, lwd=lwd.nodes)    
		}
		# add node labels with mtext
		if (show.labels) {
			graphics::mtext(labels, side=side, line=line, at=centers, cex=cex.labels, outer=outer,
					col=col.labels, las=las, font=font, adj=adj, padj=padj, ...)    
		}
	}
}


#' @title Node Coordinates
#' 
#' @description
#' Computes axis locations of each node. This function can be helpful when 
#' you want to separately plot the node labels using the function mtext.
#'
#' @param edgelist basically a two-column matrix with edges 
#' @param sorted logical to indicate if nodes should be sorted
#' @param decreasing logical to indicate type of sorting 
#' @param ordering optional numeric vector providing the ordering of nodes
#' @param labels character vector with labels for the nodes
#' @return a vector with the location of nodes in the x-axis
#' @author Gaston Sanchez
node_coords <- function(
		edgelist, sorted = FALSE, decreasing = FALSE, ordering = NULL, 
		labels = NULL)
{
	# Get graph information
	nodes_edges = graph_info(edgelist, sorted = sorted, decreasing = decreasing, 
			ordering = ordering, labels = labels)
	
	nodes = nodes_edges$nodes
	num_nodes = nodes_edges$num_nodes
	num_edges = nodes_edges$num_edges
	aux_ord = nodes_edges$aux_ord
	labels = nodes_edges$labels
	
	# x-y node coordinates
	centers = xynodes(num_nodes, aux_ord, labels)
}
if(getRversion() >= "2.15.1"){
	globalVariables(c("Arguments"))
}
