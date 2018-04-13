## THE REIDS PACKAGE ##


## IMPORTS ##

#' @import aroma.affymetrix
#' @import aroma.core
#' @import GenomeGraphs
#' @import biomaRt
#' @importFrom MCMCpack riwish
#' @importFrom data.table rbindlist 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom lmtest lrtest
#' @importFrom methods hasArg

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
#' @param Name The name of the data.
#' @param ExonSummarization Logical. Should the data be summarized at the exon level?
#' @param GeneSummarization Logical. Should the data be summarized at the gene level?
#' @param FIRMA Logical. Should the FIRMA model be performed on the data?
#' @param location The location where the .rda file is to be stored. If NULL, a list containing the requested objects is returned to the user.
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
#' DataProcessing(chipType="HTA-2_0",tags="*,r",Name="HTAData",
#' ExonSummarization=TRUE,GeneSummarization=TRUE,FIRMA=TRUE,
#' location="HTAData",verbose=TRUE)
#' }
DataProcessing<-function(chipType="HuEx-1_0-st-v2",tags="coreR3,A20071112,EP",Name="ColonCancer",ExonSummarization=TRUE,GeneSummarization=TRUE,FIRMA=TRUE,location=NULL,verbose=TRUE){
	
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
	
	
	
	if(!(is.null(location))){
		assign(Name,Data,envir=environment())
		eval(parse(text=paste("save(",Name,",UniqueGeneID,UniqueExonID,file=\"",location,"/", Name, ".RData\")", sep=""))) 			
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
		
		if(!(is.null(location))){
			assign(paste(Name,"ExonLevelSummarized",sep="_"),ExonLevelSummarized_rma,envir=environment())
			eval(parse(text=paste("save(",Name, "_ExonLevelSummarized, file=\"",location,"/",Name,"", "_ExonLevelSummarized.RData\")", sep=""))) 
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
		
		if(!(is.null(location))){
			assign(paste(Name,"GeneLevelSummarized",sep="_"),GeneLevelSummarized_rma,envir=environment())
			eval(parse(text=paste("save(",Name, "_GeneLevelSummarized, file=\"",location,"/",Name,"", "_GeneLevelSummarized.RData\")", sep=""))) 
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
		
		if(!(is.null(location))){
			assign(paste(Name,"FIRMA_Output",sep="_"),FIRMA_Output,envir=environment())
			eval(parse(text=paste("save(",Name, "_FIRMA_Output, file=\"",location,"/",Name,"", "_FIRMA_Output.RData\")", sep=""))) 
		}
	}
	
	if(is.null(location)){
		Out=list(Data,ExonLevelSummarized_rma,GeneLevelSummarized_rma,FIRMA_Output)
		names(Out)=c("Data","ExonLevelSummarized","GeneLevelSummarized","FIRMA")
		return(Out)
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
#' @param Location The location and name where the file should be saved. If NULL, the object is returned to the user. Otherwise, a file with the specified name is created.
#' @return A data frame with one row per gene. This row contains the values for each exon per sample and is convenient for processing on a HPC cluster.
#' @details All information concerning one gene is gathered. The first column of the returned data frame is the gene ID, the second column contains the exon IDs of all exons of that gene. The third colum indicates the number of probes per exon, the fourth contains the values of thos probes per sample and the last column contains the sample names.This way a .csv file is created for processing on a HPC cluster.
#' @examples
#' data(TC12000010)
#' 
#' PivotTest=PivotTransformData(Data=TC12000010,GeneID=NULL,ExonID=NULL,
#' Location=NULL)
PivotTransformData<-function(Data, GeneID=NULL,ExonID=NULL,Location=NULL){
	Data=as.data.frame(Data)
	
	if(is.null(Data$GeneID)&is.null(Data$ExonID)){
		DataTemp=data.frame(GeneID=GeneID,ExonID=ExonID)
		DataTemp=cbind(DataTemp,Data)
		Data=DataTemp
	}else if(is.null(GeneID)){
		GeneID=Data$GeneID
	}else if(is.null(ExonID)){
		ExonID=Data$ExonID
	}
	
	## From this point we assume that the first column of the data is a Gene ID and the second column is the Exon ID
	
	Transformation<-function(gID,Data){
		Subset=Data[which(Data$GeneID==gID),]
		
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
	#use rbindlist to get a full data file
	DataPivot=lapply(unique(GeneID),Transformation,Data)
	DataBind=t(rbindlist(list(DataPivot)))
	colnames(DataBind)=c("GeneID","ExonID","lengthexons","allsamples","samplenames")
	rownames(DataBind)=unique(GeneID)
	
	DataBind=as.data.frame(DataBind)
	DataBind$GeneID=as.character(DataBind$GeneID)
	DataBind$ExonID=as.character(DataBind$ExonID)
	DataBind$lengthexons=as.character(DataBind$lengthexons)
	DataBind$allsamples=as.character(DataBind$allsamples)
	DataBind$samplenames=as.character(DataBind$samplenames)
	
	if(!(is.null(Location))){
		utils::write.table(DataBind,file=Location,row.names=FALSE,col.names=TRUE,sep=",",quote=TRUE,qmethod="double")
	}
	else{
		return(DataBind)
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
#' @param nsim The number of iterations to perform.
#' @return A list with 2 items per gene. The first item is the exon scores of the corresponding probesets and the second contains a data frame with the array scores of the exons across the samples.
#' If the iniREIDS model was performed. The items will be added to the previously made list.
REIDSmodel_intern<- function(SubgeneData, nsim=100) {
	
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
	}
	else{
		ft <- stats::lm(geneExp~probes+id-1,x=TRUE)  #Difference from iniREIDS: id(samples) is involved in the calculation
		
		beta<- as.vector(summary(ft)$coefficient[,1])
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
		
		set.seed(123*i+i)
		
		# Update b
		cZZ <- diag(colSums(Z))
		#vb <- diag(nsamples)%x%Taub+tau*cZZ #This is Var =D^-1+ sigma^(-2)*ni)
		vb <- solve(diag(nsamples)%x%Taub+tau*cZZ) # This is Var^(-1) =(D^-1+ sigma^(-2)*ni)^-1
		#vbInv=solve(vb)
		mb2<-(tau*crossprod(Z,geneRes))
		mb <- apply(vb,1,function(x) sum(x*mb2))
		vb1 <- diag(vb)
		b<-sapply(c(1:(nsamples*G)),function(x) stats::rnorm(1,mean=mb[x],sd=sqrt(vb1[x])))
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
	
	tp1 <- apply(outputsigma,2,function(x) stats::quantile(x, probs = c(0.025, 0.5, 0.95)))
	tp2 <- stats::quantile(outputerror, probs = c(0.025, 0.5, 0.95))
	icc <- t(apply(tp1,2,function(x) x/(x+tp2)))
	icc <- data.frame(exon=uz,type="icc",icc )
	nik <- as.vector(t(sapply(uz,function(x) rep(x,nsamples))))
	ikEffect2  <- t(apply(outputbik,2,function(x) stats::quantile(x, probs = c(0.025, 0.5, 0.95))))
	ikEffect<- data.frame(exon=nik,type="bik",ikEffect2)
	exonscore <- icc[,c(1,4)]
	arrayscore <- ikEffect[,4] 
	
	dim(arrayscore)<- c(G,nsamples)
	rownames(arrayscore)<-ikEffect$exon[1:G] 
	colnames(arrayscore) <- colnames(SubgeneData)
	Output <- list(exonScores=exonscore,arrayScores=arrayscore)
	return(Output)
}


# REIDS Function - Cluster version
# This function is accompagnied with a .pbs file in the documentation folder. It is advised to run this model an a HPC cluster and not on a regular laptop as it will
# consume time and memory

#' "REIDS_ClusterVersion"
#' 
#' The REIDS_ClusterVersion performs the REIDS model and was adapted for use on a HPC cluster. This function should be used with the REIDS_ClusterVersion.R file and REIDS_ClusterVersion.pbs script in the documentation folder of the package.
#' After running this function on the cluster, the output files should be binded together with the CreateOutput function.
#' @export
#' @param geneData The data with as rows the probesets and as columns the samples. Note that the first column should contain the gene IDs and the second column the exon IDs
#' @param nsim The number of iterations to perform.
#' @param geneID A vector of the gene IDs.
#' @param informativeCalls Logical. Should the I/NI calls method be perform before applying the REIDS model?
#' @param rho The threshold for filtering in the I/NI calls method. Probesets with scores higher than rho are kept.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @return A .RData file will be saved for each gene with the elements returned by the iniREIDS and REIDS functions.
REIDSFunction_ClusterVersion<- function(geneData,nsim=1000,geneID,informativeCalls=TRUE,rho=0.5,Low_AllSamples){
	
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
		fit2 <- data.frame(exonNames=unique(enames),Score=fit,informative=fit>rho) # iniREIDS returns one value per exon: filtering on exon level, no replicates
		output[[1]][[i]]=fit2
		names(output[[1]])[i]="Informative"
		iniData <- lcmmData[which(rownames(lcmmData)%in%fit2$exonNames[fit2$informative]),] # Of those that pass filtering step, retrieve the replicates and samples
		lcmmData <-  iniData   
		i=i+1
	}
	
	
	if(!is.null(lcmmData)&length(unique(rownames(lcmmData)))>1){  
		fit <- REIDSmodel_intern(SubgeneData=lcmmData, nsim) 
		exonScore <-fit$exonScores
		arrayScore <- fit$arrayScores
		output[[1]][[i]]=exonScore
		names(output[[1]])[i]="exonScore"
		output[[1]][[i+1]]=arrayScore
		names(output[[1]])[i+1]="arrayScore"
	}
	
	return(output)
}


# CreateOutput

#' "CreateOutput"
#' 
#' The CreateOutput functions bind the .RData files returned by the REIDS_ClusterVersion together into one list with an element per gene. The function is advised to be used with the CreateOutput.R and CreateOutput.pbs file in the documentation folder.
#' @export
#' @param ID A data frame with a "geneID" column.
#' @param Name A name for the returned list.
#' @param Location The location where the file should be saved.
#' @return A list with the output of the REIDS_ClusterVersion binded together for all genes.
CreateOutput<-function(ID,Name,Location=""){
	Output=list()
	for(i in as.character(ID$geneID)){
		Data=get(load(paste(Location,"/REIDS_Gene_",as.character(i),".RData",sep="")))
		
		Output[length(Output)+1]=Data
		names(Output)[length(Output)]=i
	}
	
	assign(paste(Name,"_REIDS_Output",sep=""),Output,envir=environment())
	eval(parse(text=paste("save(",Name, "_REIDS_Output, file=\"",Location,"/",Name, "_REIDS_Output.RData\")", sep=""))) 
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
#' @param nsim The number of iterations to perform.
#' @param informativeCalls Logical. Should the I/NI calls method be perform before applying the REIDS model?
#' @param rho The threshold for filtering in the I/NI calls method. Probesets with scores higher than rho are kept.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @param Location A character string indication the place where the output should be saved.
#' @return An .RData file for each gene with the values returned by the iniREIDS and REIDS functions.
reidsfunction_genebygene <- function(file_name,file_pos,line_length,nsim,informativeCalls=TRUE,rho=0.5,Low_AllSamples,Location){
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
			
			
			if(!is.null(lcmmData)&length(unique(rownames(lcmmData)))>1){  
				fit <- REIDSmodel_intern(SubgeneData=lcmmData, nsim) 
				exonScore <-fit$exonScores
				arrayScore <- fit$arrayScores
				output[[1]][[i]]=exonScore
				names(output[[1]])[i]="exonScore"
				output[[1]][[i+1]]=arrayScore
				names(output[[1]])[i+1]="arrayScore"
			}
			
		
		save(output, file=paste(Location,"/REIDS_Gene_", geneID, ".RData", sep=""))
		rm(output)
		gc()
	}	
	
	close(conn)
}


#' "REIDSFunction"
#' 
#' The REIDSFunction performs the REIDS model on the pivot transformed data by calling on the line indexed file. The REIDS model is performed gene by gene and the returned outputs are knitted together.
#' @export
#' @param geneIDs A data frame with a "geneID" column.
#' @param Name A name for the returned list.
#' @param Indices The .csv file created by Line_Indexer.py which contains indices for every gene.
#' @param DataFile The .csv file created by PivotTransformation. 
#' @param nsim The number of iterations to perform.
#' @param informativeCalls Logical. Should the I/NI calls method be perform before applying the REIDS model?
#' @param rho The threshold for filtering in the I/NI calls method. Probesets with scores higher than rho are kept.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @param Location A character string indication the place where the outputs are saved.
#' @return A list with an element for each gene with per gene the values returned by the iniREIDS and REIDS functions.
REIDSFunction<-function(geneIDs,Name,Indices,DataFile,nsim=5000,informativeCalls=TRUE,rho=0.5,Low_AllSamples,Location){
	Lines=utils::read.table(Indices,header=TRUE,sep=",",stringsAsFactors=FALSE)
	Lines=data.frame(Lines)	
	Lines[,1]=as.numeric(Lines[,1])
	Lines[,2]=as.numeric(Lines[,2])
	
	REIDSOut=apply(Lines,2, function(x) reidsfunction_genebygene(file_name=DataFile,file_pos=x[1],line_length=x[2],nsim,informativeCalls,rho,Low_AllSamples,Location))
	
	CreateOutput(ID=geneIDs,Name,Location)
}

##### OUTPUT ANALYSIS FUNCTIONS #####

#FILTERING FUNCTION

#' "FilterInformativeGenes"
#'
#' The FilterInformativeGenes function works on the list returned by the REIDS functions. It returns a list that only retains the genes that have informative probesets.
#' @export 
#' @param Data The output of the REIDS function.
#' @return A subset of the list provided in Data.
#' @examples
#' data(TC12000010_REIDS_Output)
#' Test_F=FilterInformativeGenes(TC12000010_REIDS_Output)
FilterInformativeGenes <- function(Data){
	Out1=list()
	filter <- function(Subset){
		if("exonScore"%in% names(Subset)){
			return(Subset)
		}
	}
	Out1=lapply(Data,filter)
	
	Out2=Out1[vapply(Out1, Negate(is.null), NA)]
	
	return(Out2)
}


#SEARCH FUNCTION

#' "Search"
#' 
#' The Search function investigates whether specific genes and/or exons are present in the data. It returns the elements of the found elements and a list of those that were not found.
#' @export 
#' @param WhatToLookFor A data frame with a GeneID and ExonID column. Only the ExonID column is necessary, the corresponding gene ID will be sought in the data
#' @param Data The data in which the given gene and exon ID should be sought.
#' @param AggregateResults Logical. Should the results be aggregated on gene level? This results in a list with one item per gene
#' @param NotFound Not be specified by the user.
#' @return The returned value is a list with two elements. The first element is SearchResults which contains a list of the found instances in the data. The list contains an element per found instance. If AggregateResults is TRUE, the list is reduced to an element per gene. The secod is a data frame calles NotFound which contain the gene and exon ID which were not found in the data. If only exon ID were specified, the gene ID are NA in this data frame.
#' @examples
#' data(TC12000010_REIDS_Output)
#' Test_F=FilterInformativeGenes(TC12000010_REIDS_Output)
#' Test_S=Search(WhatToLookFor=data.frame(ExonID=c("PSR12000150")),Data=Test_F,
#' AggregateResults=FALSE,NotFound=NULL)
Search <- function(WhatToLookFor=data.frame(GeneID=NULL,ExonID=NULL), Data, AggregateResults=FALSE,NotFound=NULL){
	if(!(is.null(WhatToLookFor$GeneID))){
		
		if(is.null(NotFound)){
			NotFound<-data.frame(GeneID=integer(),ExonID=integer())
		}
		
		if(class(Data[[1]])=="list"){
			retrieve<-function(GID,EID,Data){
				
				listnumber=which(names(Data)==as.character(GID))
				
				if(length(listnumber)==0){
					NotFound=rbind(NotFound,data.frame(geneID=GID,exonID=EID))
					#assign("NotFound", notfound, envir = .GlobalEnv)
					
					ExonInfo=NULL
					#return(ExonInfo)
				}
				else{
					SelectGeneData=Data[[listnumber]]
					

					if(class(SelectGeneData)=="list"){ #here we assume that the data is the direct output of the REIDSfunction
						
						if(names(SelectGeneData)[1]=="exonScore"){
							ExonScore=SelectGeneData$exonScore[which(SelectGeneData$exonScore$exon==EID),,drop=FALSE]
							ArrayScore=SelectGeneData$arrayScore[which(rownames(SelectGeneData$arrayScore)==EID),,drop=FALSE]	
							
							ExonInfo=list(exonScore=ExonScore,arrayScore=ArrayScore)
						}
						else{
							InformativeExon=SelectGeneData$Informative[which(SelectGeneData$Informative$exonNames==EID),,drop=FALSE]
						
							if(InformativeExon$informative==TRUE & length(which(SelectGeneData$Informative$informative==TRUE))>1){
								ExonScore=SelectGeneData$exonScore[which(SelectGeneData$exonScore$exon==EID),,drop=FALSE]
								ArrayScore=SelectGeneData$arrayScore[which(rownames(SelectGeneData$arrayScore)==EID),,drop=FALSE]
						
								ExonInfo=list(Informative=InformativeExon,exonScore=ExonScore,arrayScore=ArrayScore)
							}
							else{
								ExonInfo=list(Informative=InformativeExon)
							}
						}
					
					}
					else if(class(SelectGeneData)=="data.frame"){
						ExonInfo=SelectGeneData[which(SelectGeneData$ExonID==EID),,drop=FALSE]
					
						if(nrow(ExonInfo)==0){
							ExonInfo=NULL
						
							NotFound=rbind(NotFound,data.frame(geneID=GID,exonID=EID))
							#assign("NotFound", notfound, envir = .GlobalEnv)
						}
					
					}
				}	
				
				return(ExonInfo)
				
			}
			
		}
		
		else if(class(Data)=="data.frame"){
			retrieve<-function(GID,EID,Data){
				
				Rownumbers=which(Data$GeneID==GID)
				if(length(Rownumbers)==0){
					NotFound=rbind(NotFound,data.frame(geneID=GID,exonID=EID))
					
					ExonInfo=NULL

				}
				else{
					
					SelectGeneData=Data[Rownumbers,]
					
					ExonInfo=SelectGeneData[which(as.character(SelectGeneData$ExonID)==as.character(EID)),,drop=FALSE]
					
					if(nrow(ExonInfo)==0){
						ExonInfo=NULL
						
						NotFound=rbind(NotFound,data.frame(geneID=GID,exonID=EID))
						#assign("NotFound", notfound, envir = .GlobalEnv)
					}			
				}
			
				return(ExonInfo)		
			}
		}
		
		Out1=lapply(1:nrow(WhatToLookFor),function(x) retrieve(GID=as.character(WhatToLookFor[x,1]),EID=as.character(WhatToLookFor[x,2]),Data = Data))
		names(Out1)=WhatToLookFor$GeneID
		Out2=Out1[vapply(Out1, Negate(is.null), NA)]
		
		SearchResults=Out2
		
		if(AggregateResults==TRUE){
			message("Results are aggregated on the gene level...")
			
			genes=unique(names(Out2))
			
			aggregate<-function(gene,Data){
				if(class(Data[[1]])=="data.frame"){  #assume a data.frame per gene
					Out3=rbindlist(Data)
					Out3=as.data.frame(Out3)
				}
				else if(class(Data[[1]])=="list"){ #assume a list of data frames per gene
					
					ntables=max(sapply(Data, length))
					
					aggregatetables<-function(data){
						data=lapply(data,as.data.frame)
						Out4=rbindlist(data)
						Out4=as.data.frame(Out4)
						return(Out4)					 
					}
					
					Out3=lapply(1:ntables,function(i) aggregatetables(data=lapply(Data,function(x) ifelse(length(x)>=i,x[i],x[NA]))))	
					names(Out3)=names(Data[[ which(sapply(Data,length)==max(sapply(Data, length)))[1]]])
					
				}
				
				return(Out3)
				
			}
			
			Out5=lapply(1:length(genes),function(x) aggregate(gene = genes[x], Data=Out2[which(names(Out2)==genes[x])]))
			names(Out5)=genes
			
			SearchResults=Out5
			
		}
		
		
	}
	
	else if(is.null(WhatToLookFor$GeneID) & !(is.null(WhatToLookFor$ExonID))){
		
		if(class(Data[[1]])=="list"){
			
			if(names(Data[[1]])[1]=="exonScore"){
				ExonNames=lapply(1:length(Data),function(x) as.character(Data[[x]]$exonScore$exon))
			}
			else{
				ExonNames=lapply(1:length(Data),function(x) as.character(lapply(Data[x],"[[",1)[[1]]$exonNames))
			}
			
			names(ExonNames)=names(Data)		
			GeneNames=sapply(1:length(ExonNames),function(i) rep(names(ExonNames)[i],length(ExonNames[[i]])))			
			GenesAndExons=data.frame(GeneID=as.character(unlist(GeneNames)),ExonID=as.character(unlist(ExonNames)))
			
			exons=WhatToLookFor$ExonID
			
			LookFor=lapply(exons, function(x) GenesAndExons[which(as.character(GenesAndExons$ExonID)==x),])			
			LookFor=as.data.frame(rbindlist(LookFor))
			
			exons_notfound=exons[which(!exons%in%LookFor$ExonID)]
			
			if(length(exons_notfound)==0){
				NF=NULL
			}
			else{
				NF=data.frame(GeneID=rep(NA,length(exons_notfound)),ExonID=exons_notfound)
			}
			
		}
		else if(class(Data)=="data.frame"){
			
			#ExonNames=lapply(1:length(Data),function(x) Data[x][[1]]$ExonID)
			ExonNames=Data$ExonID
			GeneNames=Data$GeneID
		
			
			#GeneNames=sapply(1:length(ExonNames),function(i) rep(names(ExonNames)[i],length(ExonNames[[i]])))			
			#GenesAndExons=data.frame(geneID=as.integer(unlist(GeneNames)),exonID=as.integer(unlist(ExonNames)))
			GenesAndExons=data.frame(GeneID=GeneNames,ExonID=ExonNames)
			exons=WhatToLookFor$ExonID
			
			LookFor=lapply(exons, function(x) ifelse(nrow(GenesAndExons[which(GenesAndExons$ExonID==x),])!=0,return(GenesAndExons[which(GenesAndExons$ExonID==x),]),return(data.frame(GeneID="NA",ExonID=x))))
			LookFor=as.data.frame(rbindlist(LookFor))
			LookFor$GeneID=suppressWarnings(as.integer(as.character(LookFor$GeneID)))
			exons_notfound=LookFor[which(is.na(LookFor$GeneID)),]
			LookFor=LookFor[-which(is.na(LookFor$GeneID)),]
			
			if(nrow(exons_notfound)==0){
				NF=NULL
			}
			else{
				NF=exons_notfound
			}
			
		}
		
		return(Search(WhatToLookFor=LookFor, Data=Data, AggregateResults=AggregateResults,NotFound=NF))
	}
	
	NMatched=length(Out2)
	message(paste(NMatched," instances found in the Data..."))
	return(list("SearchResults"=SearchResults,"NotFound"=NotFound))
}


#' "ExonTesting"
#' 
#' The ExonTesting function performs a t-test (2 groups) or F-test (more than 2 groups) on the array score of predefined groups. If specified, probesets are filtered out on exon scores and test significance.
#' @export
#' @param Data The Data on which testing of the array scores should be conducted. This is preferably output of the REIDS function
#' @param Exonthreshold The exon score threshold to be maintained. If not NULL, probesets with an exon score lower than this value are not considered further and the p-values will be adjusted for multiplicity after testing. If NULL, all probesets are considered and a multiplicity correction is not performed.
#' @param Groups A list with elements specifying the columns of the data in each group.
#' @param paired Logical. Are the groups paired? If TRUE the mean paired differences are calculated and tested whether these are significantly different from zero or not.
#' @param significancelevel The significance level to be maintained on the p-values. The filtering on the significance is conducted only if an Exonthreshold is specified and the p-value are adjusted for multiplicity.
#' @return A data frame with one line per exon. The columns contain the gene ID, the exon ID, the test statistic, a p-value and an adjusted p-value. If the groups are paired also the mean paired difference is given. The p-values are adjusted for multiplicity and filtered on significance if significancelevel is not NULL.
#' @examples
#' data(TC12000010_REIDS_Output)
#' ExonTest=ExonTesting(Data=TC12000010_REIDS_Output,Exonthreshold=NULL,Groups=list(c(1:9),c(10:18)),
#' paired=FALSE,significancelevel=NULL)
ExonTesting <- function(Data, Exonthreshold=NULL,Groups=list(),paired=FALSE,significancelevel=NULL){
	
	if(is.null(Exonthreshold)){
		Exonthreshold=0
	}
	
	if(is.null(Groups[[1]])){
		message("A test on the array scores will NOT be performed")
	}
	
	message("Data is filtered to only contain genes with informative exons")
	Data_Filtered1=FilterInformativeGenes(Data)  #only those genes with informative exons remain
	
	
	filterASExons<-function(i,geneID,Subset,Exon,groups,pairing){
		
		## Step 1 : filtering on the Exon value. If 0, no filtering occurs
		ExonsExon = Subset$exonScore
		ArrayScores = Subset$arrayScore
		
		SelectExonsExon = ExonsExon[which(ExonsExon$X50>Exon),]
		ExonsPassedExon = SelectExonsExon$exon
		
		SelectRowsArray = which(rownames(ArrayScores)%in%ExonsPassedExon)
		
		## Step 2 : Testing of the Array Scores -- paired or not paired
		if(pairing==FALSE){  # Test between two groups of Array Scores
			
			names(groups)=c(1:length(groups))
			
			ArrayScoreTTest=matrix(0,nrow=length(SelectRowsArray),ncol=3)
			colnames(ArrayScoreTTest)=c("statistic","p.value","adj.p.value")
			rownames(ArrayScoreTTest)=rownames(ArrayScores[SelectRowsArray,])
			
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
			
			
			if(length(groups)<=2){
				ArrayScoreTTest[,c(1,2)] = t(sapply(SelectRowsArray,function(i) {ttest(data=ArrayScores[i,,drop=FALSE],groups=groups, pairs = paired)}))	
			}
			else{
				Groups=rep(0,length(unlist(groups)))
				for(j in 1:length(groups)){
					positions=groups[[j]]
					Groups[positions]=names(groups)[j]
				}
				Groups=as.numeric(Groups)
				Groups=factor(Groups)
				
				ArrayScoreTTest[,c(1,2)] = t(sapply(SelectRowsArray,function(i) {anovaFtest(data=ArrayScores[i,,drop=FALSE],Groups=Groups)}))	
				
			}
			p_vals=as.vector(as.matrix((ArrayScoreTTest[,grep("^(p.value)",colnames(ArrayScoreTTest))])))
			adj_p_vals=matrix(stats::p.adjust(p_vals,"BH"),nrow=nrow(ArrayScoreTTest),ncol=length(grep("^(p.value)",colnames(ArrayScoreTTest))))
			ArrayScoreTTest[,which(seq(1,ncol(ArrayScoreTTest))%%3==0)]=adj_p_vals
			
			ArrayScoreTTest=as.data.frame(ArrayScoreTTest)
			
		}
		else if(pairing==TRUE){ # Test the mean paired difference against zero
			
			
			if(!(is.null(groups[[1]])) & length(SelectRowsArray)!=0){
				
				ArrayScore_group1=ArrayScores[SelectRowsArray,groups[[1]],drop=FALSE]
				ArrayScore_group2=ArrayScores[SelectRowsArray,groups[[2]],drop=FALSE]
				
				mean_paired_diff<-function(g1,g2){
					Paired_Diff=g1-g2
					Mean_Diff=mean(Paired_Diff)
					
					out1=stats::t.test(x=Paired_Diff)
					out2=cbind(out1$statistic,out1$p.value)
					out3=cbind(Mean_Diff,out2)
					
					return(out3)
				}
				
				ArrayScoreTest = t(sapply(c(1:length(ExonsPassedExon)),function(i) mean_paired_diff(g1=ArrayScore_group1[i,],g2=ArrayScore_group2[i,]) ))
				ArrayScoreTest=data.frame("Mean_Diff"=ArrayScoreTest[,1],"t-statistic"=ArrayScoreTest[,2],p.value=ArrayScoreTest[,3])
				rownames(ArrayScoreTest)=rownames(ArrayScore_group1)
				
				
			}
			else{
				ArrayScoreTest = NULL
			}
		}	
		
		if(!is.null(ArrayScoreTTest)){
			Output=cbind("geneID"=rep(geneID,length(ExonsPassedExon)), SelectExonsExon,ArrayScoreTTest)	
		}
		else{
			Output=cbind("geneID"=rep(geneID,length(ExonsPassedExon)),SelectExonsExon)
		}
		
		colnames(Output)[2]="ExonID"
		
		if(nrow(Output)==0){
			Output=NULL
		}
		
		return(Output)
		
	}
	
	Out1=lapply(1:length(Data_Filtered1), function(i) filterASExons(i,geneID = names(Data_Filtered1[i]), Subset = Data_Filtered1[[i]], Exon = Exonthreshold , groups=Groups, pairing = paired))
	names(Out1)=names(Data_Filtered1)
	
	Out2=Out1[vapply(Out1, Negate(is.null), NA)]
	
	Out3<-do.call(rbind.data.frame, Out2)
	
	
	
	## Step 3 : if Exon threshold specified: we adjust for multiplicity ;  filter on significance level 
	if(nrow(Out3)!=0){
		message("The p-value are adjusted")
		if(!is.null(Out3$p.value)){
			Out3$adj.p.value=stats::p.adjust(Out3$p.value,"fdr")
		}	
		if(!is.null(significancelevel)){
			Out3=Out3[which(Out3$adj.p.value<=significancelevel),]
		}
		rownames(Out3)=seq(1:nrow(Out3))
	}
	
	return(Out3)
	
}


## Identification of the AS exons

#' "ASExons"
#' 
#' The ASExons functions can be performed either on the output of the REIDS model or on the ExonTesting model and identifies the alternatively spliced exons. It filters probesets on their exon scores, adjusts p-values for multiplicity and only keeps the significant probesets.
#' @export
#' @param Data The Data on which testing of the array scores should be conducted. This can be either the output of the REIDS model or the ExonTesting function.
#' @param Exonthreshold The exon score threshold to be maintained. If not NULL, probesets with an exon score lower than this value are not considered further and the p-values will be adjusted for multiplicity after testing. If NULL, all probesets are considered and a multiplicity correction is not performed.
#' @param Groups A list with elements specifying the columns of the data in each group.
#' @param paired Logical. Are the groups paired? If TRUE the mean paired differences are calculated and tested whether these are significantly different from zero or not.
#' @param significancelevel The significance level to be maintained on the p-values. The filtering on the significance is conducted only if an Exonthreshold is specified and the p-value are adjusted for multiplicity.
#' @return A data frame with one line per exon. The columns contain the gene ID, the exon ID, the test statistic, a p-value and an adjusted p-value. If the groups are paired also the mean paired difference is given. Only the probesets with high enough exon scores and a significant test are kept in the data frame.
#' @examples 
#' data(TC12000010_REIDS_Output)
#' ASTest=ASExons(Data=TC12000010_REIDS_Output,Exonthreshold=0.5,Groups=list(c(1:9),c(10:18)),
#' paired=FALSE,significancelevel=0.05)
ASExons<-function(Data,Exonthreshold=0.5,Groups=list(group1=NULL,group2=NULL),paired=FALSE,significancelevel=0.05){
	
	if(class(Data)=="list"){
		message("The data is assumed to be output of the REIDS model. Filtering of the probesets and testing of the array scores will be performed")
		message("The used threshold for the exon scores is 0.5")
		message("The used significance level for the p-values is 0.05")
		
		TestedData=ExonTesting(Data=Data,Exonthreshold=Exonthreshold,Groups=Groups,paired=paired,significancelevel=significancelevel)
		
		if(nrow(TestedData)!=0){
			message("Ordering data in from high tolow exon scores")
			Data_Sign_Ordered=TestedData[order(-TestedData$X50.),]
			rownames(Data_Sign_Ordered)=c(1:nrow(Data_Sign_Ordered))
		}
		else{
			Data_Sign_Ordered=NULL
		}
		
	}
	else if(class(Data)=="data.frame"){
		
		message("In using this function please make sure that the data has not been filtered yet and still has unadjusted p-values.")
		
		message(paste("Keep probesets with exon score greater than",Exonthreshold,sep=" "))
		Data_Filt1=Data[which(Data$X50.>Exonthreshold),]
		
		if(nrow(Data_Filt1)!=0){
			message("Adjusting p-values for multiplicity")	
			Data_Filt1$adj.p.value=stats::p.adjust(Data_Filt1$p.value,"fdr")
		
			message(paste("Keep probesets with a p-value lower than",significancelevel,sep=" "))
			Data_Sign=Data_Filt1[which(Data_Filt1$adj.p.value<significancelevel),]

			if(nrow(Data_Sign)!=0){
				message("Ordering data in from high to low exon scores")
				Data_Sign_Ordered=Data_Sign[order(-Data_Sign$X50.),]
				rownames(Data_Sign_Ordered)=c(1:nrow(Data_Sign_Ordered))
			}
			else{
				Data_Sign_Ordered=NULL
			}

		}
		else{
			Data_Sign_Ordered=NULL
		}
	}
	return(Data_Sign_Ordered)
	
}

#' JunInfo
#' 
#' JunInfo functions asses the junction information for a single gene
#' @export
#' @param x The TC ID for which to retrieve and asses the junction information .
#' @param ASPSR The AS probe sets as identified by ASExons.
#' @param JLines The lines which contain information on the TC ID in the junction association file.
#' @param TrLines The lines which contain information on the TC ID in the transcript annotation file.
#' @param ELines The lines which contain information on the TC ID in the exon annotation file.
#' @param DataS The TC ID subset of the probe level data.
#' @param Groups A list with  elements speficifing the columns of the data in each group.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @param Low_GSamples A list with a  character vector per group containing the probe sets which are not DABG in that group.
#' @param Plot Should a plot of the gene model be made?
#' @param Name A character string with the name of the ouput file.
#' @details The plot is produced by the arcplot function of the arcdiagram package (https://github.com/gastonstat/arcdiagram)
#' @return The function returns three files. The first file has name "Name.txt" and contains a line per probe set. It shows the reached decision
#' regarding the probe set (const/AS/not DABG),its linking and exclusion junctions, the fold change, the AS type and its annotated exons. The second
#' file is a list of all found transcripts for a particular TC I. The third file indicates whether a specific transcript is present or absent in a group.
JunInfo<-function(x,ASPSR,JLines,TrLines,ELines,DataS,Groups,Low_AllSamples=c(),Low_GSamples=c(),Plot,Name){
	
	EAnnot=utils::read.table("HTA-2_0_ExonAnnotations.txt",header=FALSE,sep="\t",nrows=as.numeric(ELines[3]),skip=as.numeric(ELines[2]),stringsAsFactors=FALSE)
	JAnnot=utils::read.table("HTA-2_0_JunAssociations.txt",header=FALSE,sep="\t",nrows=as.numeric(JLines[3]),skip=as.numeric(JLines[2]),stringsAsFactors=FALSE)
	colnames(JAnnot)=c("TC_ID","PSR_ID","JUC_ID","as_type")		
	Transcripts=utils::read.table("HTA-2_0_TranscriptAnnotations.txt",header=FALSE,sep="\t",nrows=as.numeric(TrLines[3]),skip=as.numeric(TrLines[2]),stringsAsFactors=FALSE)
	
	
	colnames(DataS)[2]="ExonID"
	if(any(is.na(JAnnot$PSR_ID))){
		D=which(is.na(JAnnot$PSR_ID))
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
					}
					GData=GData[-c(which(rownames(GData)%in%Flagged)),]
				}
			}
		}
		
		GroupData[[g]]=GData			
	}
	
	
	print(x)	
	
	DelJuncs=names(which(table(DelJ)==length(Groups)))
	Jucs=unique(JAnnot$JUC_ID)[which(!is.na(unique(JAnnot$JUC_ID)))]
	
	
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
			
			Set=sort(JAnnot[which(JAnnot$JUC_ID==j&(JAnnot$as_type!="exclusion")),2])
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
							SelectAdd=which((I>t)&(I<=(t+Add)))[Add]
							Select=c(Select,SelectAdd)
						}
						else if(length(Select)>t){
							Del=length(Select)-t
							SelectDel=which(Select>=t)[Del]
							Select=Select[-c(Del)]
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
			Set=sort(JAnnot[which(JAnnot$JUC_ID==j&(JAnnot$as_type!="exclusion")),2])
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
			PSRs=unique(JAnnot$PSR_ID[which(JAnnot$JUC_ID==J)])
			supp=sort(JAnnot$PSR_ID[which(JAnnot$PSR_ID%in%PSRs&JAnnot$JUC_ID==J&JAnnot$as_type!="exclusion")])
			if(!all(supp%in%DataS$ExonID)){
				supp=supp[-c(which(!supp%in%DataS$ExonID))]
			}
			excl=sort(JAnnot$PSR_ID[which(JAnnot$PSR_ID%in%PSRs&JAnnot$JUC_ID==J&JAnnot$as_type=="exclusion")])
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
	Exons=sort(unique(c(unique(JAnnot$PSR_ID),c(unique(PSRsExons),unique(DABGPSR)))))
	Exons=Exons[which(Exons%in%DataS[,2])]
	if(length(Exons)==0){
		return(c("No PSR's present"))
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
		if(c%in%DABG){
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
			set=Transcripts[which(Transcripts[,2]==RLinks[r,1]),5]
			#set=strsplit(RLinks[r,2],"-")[[1]]		
			#NeverLinkedPSR[[r]]=set
			NeverLinkedTr[[r]]=set
		}
		NeverLinkedTr=unique(unlist(NeverLinkedTr))
	}
	
	
	#Event Type
	TC=DataS[1,1]
	Strand=Transcripts[1,3]
	TRS=list()
	N=c()
	for(tr in unique(Transcripts[,5])){
		if(tr%in%NeverLinkedTr){
			next
		}
		Set=Transcripts[which(Transcripts[,5]==tr),2]
		PSet=Set[which(substr(Set,1,1)=="P")]
		ESete=Transcripts[which(Transcripts[,5]==tr),4]
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
			for(d in DABGPSR){
				if(any(PSet==d)){
					next
				}
			}	
		}
		
		if(Strand=="-"){
			PSet=rev(sort(PSet))
			ESet=rev(sort(ESet))
		}
		names(ESet)=PSet
		TRS[[length(TRS)+1]]=ESet	
		N=c(N,tr)
	}
	names(TRS)=N
	
	
	AllExons=unique(unlist(TRS))		
	for(e in Exons){
		EAnnotp=EAnnot[which(EAnnot[,2]==e),4]
		EAnnotp=EAnnotp[which(EAnnotp%in%AllExons)]
		EAnnotp=paste(EAnnotp,collapse="|")
		s=JAnnot[which(JAnnot$PSR_ID==e),c(3,4)]
		Dels=s[which(s[,1]%in%DelJ),1]
		if(length(Dels)>0){
			s=s[which(!s[,1]%in%DelJ),]
		}
		else{
			Dels="-"
		}
		if(e%in%DABG){
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
	
	OutList[,9]=as.numeric(OutList[,9])
	
	
	Table1=matrix(0,nrow=length(TRS)*2,ncol=3)
	for(i in 1:length(TRS)){
		N=names(TRS)[i]
		R1=paste(names(TRS[[i]]),collapse="|")
		R2=paste(TRS[[i]],collapse="|")
		Table1[i*2-1,]=c(TC,N,R1)
		Table1[i*2,]=c(TC,N,R2)
	}
	utils::write.table(Table1, "Transcripts.txt", sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
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
				set=Transcripts[which(Transcripts[,2]==RLinks[r,1]),5]
				NeverLinkedTr[[r]]=set
				#set=strsplit(RLinks[r,1],"-")[[1]]
				#NeverLinkedPSR[[r]]=set
			}
			NeverLinkedTr=unique(unlist(NeverLinkedTr))
		}	
		for(tr in unique(Transcripts[,5])){
			if(tr%in%NeverLinkedTr){
				next
			}
			Set=Transcripts[which(Transcripts[,5]==tr),2]
			PSet=Set[which(substr(Set,1,1)=="P")]
			ESete=Transcripts[which(Transcripts[,5]==tr),4]
			ESet=ESete[which(substr(Set,1,1)=="P")]
			
			
			if(length(Low_GSamples[[g]])>0){
				if(any(PSet%in%Low_GSamples[[g]])){
					next
				}
			}
			
			if(Strand=="-"){
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
	utils::write.table(Table2, "Transcripts_Groups.txt", sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
	
	
	if(Strand=="-"){
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
		EAnnotp=EAnnot[which(EAnnot[,2]==PSR),]
		if(nrow(EAnnotp)==0){
			Category[r]="Intron Retention"	
			next
		}	
		else{
			Temp=c()
			for(ea in 1:nrow(EAnnotp)){
				OtherPSR=as.character(EAnnot[which(EAnnot[,4]==as.character(EAnnotp[ea,4])),2])
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
			if(Strand=="-"){
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
	utils::write.table(OutList, paste(Name,".txt",sep=""), sep = "\t", col.names = FALSE, row.names=FALSE,append = TRUE,quote=FALSE)
	rm(OutList)
	gc() 
	#return(list("Output"=OutList,"Transcripts"=GroupTranscripts,"PlotInfo"=list(G,Widths,Cols)))
}


#' REIDS_JunctionAssessment
#' 
#' The REIDS_JunctionAssessment functions assess identified AS exons based on their 5'end and 3'end and exclusion junction support.
#' @export
#' @param ASProbeSets The AS probe sets as identified by ASExons.
#' @param JAnnotI A table with the line indexing for the juncion associations.
#' @param EAnnotI A table with the line indexing for the exon annotations.
#' @param TrAnnotI A table with the line indexing for the transcript annotations.
#' @param Data The probe level data.
#' @param Groups A list with  elements speficifing the columns of the data in each group.
#' @param Low_AllSamples A character vector containing the probe sets which are not DABG in all samples.
#' @param Low_GSamples A list with a  character vector per group containing the probe sets which are not DABG in that group.
#' @param Name A character string with the name of the ouput file.
#' @return The function returns three files. The first file has name "Name.txt" and contains a line per probe set. It shows the reached decision
#' regarding the probe set (const/AS/not DABG),its linking and exclusion junctions, the fold change, the AS type and its annotated exons. The second
#' file is a list of all found transcripts for a particular TC ID. The third file indicates whether a specific transcript is present or absent in a group.
REIDS_JunctionAssesment<-function(ASProbeSets=c(),JAnnotI,EAnnotI,TrAnnotI,Data,Groups=list(c(3,4,5),c(6,7,8)),Low_AllSamples=c(),Low_GSamples=c(),Name){
	
	Genes=unique(Data$GeneID)
	writeLines(c("TC_ID\tPSR_ID\tType\tUnreliable Junctions\tLinking Junctions\tExclusion Junctions\tSupported by\tIdentified by\tFold Change\tExons"),paste(Name,".txt",sep=""),sep="\n")
	ScoresOutput=lapply(Genes, function(x) JunInfo(x,ASPSR=ASProbeSets,JLines=JAnnotI[which(JAnnotI[,1]==x),],TrLines=TrAnnotI[which(TrAnnotI[,1]==x),],ELines=EAnnotI[which(EAnnotI[,1]==x),],DataS=Data[which(as.character(Data[,1])==x),],Groups=Groups,Low_AllSamples,Low_GSamples,Plot=FALSE,Name))
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
#' =TC12000010_ExonLevel,Groups=list(c(10:18),c(19:27)),ylabel="",
#' title="PSR12000150")
#' }

ExpressionLevelPlot<-function(GeneID=NULL,ExonID=NULL,Data,GeneLevelData=NULL,ExonLevelData=NULL,Groups,ylabel=NULL,title=""){
	message("The gene and exon level data will be log2 transformed")
	
	if(is.null(GeneID) | is.null(ExonID)){
		stop("no GeneID and/or ExonID specified")
	}
	
	
	Exon_ExonID=ExonLevelData[which(ExonLevelData$ExonID==ExonID),]
	Exon_ExonID=Exon_ExonID[,-c(1,2)]
	#Exon_ExonID=log2(Exon_ExonID)
	Exon_ExonID=Exon_ExonID[,unlist(Groups)]
	
	
	# Gene Level Data
	Gene_GeneID=GeneLevelData[which(GeneLevelData$GeneID==GeneID),]
	#Gene_GeneID[,-c(1)]=log2(Gene_GeneID[,-c(1)])
	Gene_GeneID=Gene_GeneID[,-c(1)]
	Gene_GeneID=Gene_GeneID[,unlist(Groups)]
	

	# Observed Probe intensities
	ProbeIntensities=Data[which(Data$ExonID==ExonID),]
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
	
	graphics::plot(0,0,typ="n",xlab="",ylab=ylabel,ylim=c(0,max_ylim+2),xlim=c(1,ncol(ProbeIntensities_ExonID)),xaxt="n",frame=TRUE,,cex.axis=1.5)
	graphics::lines(x=c(1:ncol(ProbeIntensities_ExonID)),y=Gene_GeneID) #Gene level data...
	graphics::lines(x=c(1:ncol(ProbeIntensities_ExonID)),y=Exon_ExonID,col="blue") #Exon Level data
	for(i in 1:nrow(ProbeIntensities_ExonID)){
		graphics::points(x=c(1:ncol(ProbeIntensities_ExonID)),y=as.matrix(ProbeIntensities_ExonID[i,]),pch=19,col="blue")
	}
	graphics::axis(1,labels=colnames(Gene_GeneID),at=c(1:ncol(ProbeIntensities_ExonID)),las=2,cex.axis=1.5)	
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
							assignment <- "--"} 
					}))
	gene.symbols[is.na(gene.symbols)]	<- "--"
	annotate      <- data.frame(transcript = transcripts, symbol = gene.symbols)
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
			gene.pos<-gene.positions[gene.positions[,1]==paste(Highlight),]
			print(paste("Highlighted Region: ",paste(gene.pos,collapse=" ")))
			rOverlay <- makeRectangleOverlay(start = (as.numeric(gene.pos$start)-500), end = (as.numeric(gene.pos$stop)+500),region=c(2,3),dp = DisplayPars(alpha = .5, fill = "pink"))
			if(is.null(Start)&is.null(Stop)){
				gdPlot(list(exon,gene,transcript),overlays = rOverlay)
			}
			else{
				gdPlot(list(exon,gene,transcript),Start,Stop,overlays = rOverlay)
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
