% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{PivotTransformData}
\alias{PivotTransformData}
\title{"PivotTransformation"}
\usage{
PivotTransformData(Data, GeneID = NULL, ExonID = NULL,
  REMAPSplitFile = NULL, NotAnnotated = FALSE, Location = NULL,
  Name = "Pivot")
}
\arguments{
\item{Data}{The data frame to be transformed.}

\item{GeneID}{A character vector of the the gene IDs that correspond to the rows of the data frame. Necessary if no GeneID column is present in the data frame}

\item{ExonID}{A character vector of the the gene IDs that correspond to the rows of the data frame. Necessary if no ExonID column is present in the data frame}

\item{REMAPSplitFile}{The name of the file with the REMAP information regarding the split of the probe sets if the TC ID is annotated to mutiple genes.}

\item{NotAnnotated}{Logical. Should the probe sets which are not annotated to a gene still be included? If FALSE, these are excluded. If TRUE, these are included. Default is FALSE.}

\item{Location}{The location where the file should be saved. If NULL, the object is returned to the user. Otherwise, a file with the specified name is created.}

\item{Name}{The name of the output file. Defaults to "Pivot".}
}
\value{
A data frame with one row per gene. This row contains the values for each exon per sample and is convenient for processing on a HPC cluster. Futher also a data frame with a column of the gene ID's is returned.
}
\description{
The PivotTransformation function converts a data frame with multiple rows per gene into a .csv file with one row per gene. This is the first step in data transformation to apply the REIDS function on a HPC Cluster.
}
\details{
All information concerning one gene is gathered. The first column of the returned data frame is the gene ID, the second column contains the exon IDs of all exons of that gene. The third colum indicates the number of probes per exon, the fourth contains the values of thos probes per sample and the last column contains the sample names.This way a .csv file is created for processing on a HPC cluster.
}
\examples{
data(TC12000010)

PivotTest=PivotTransformData(Data=TC12000010,GeneID=NULL,ExonID=NULL,
Location=NULL)

\dontrun{
data(TC1500264)

PivotTransformData(Data=TC1500264,GeneID=NULL,ExonID=NULL,
REMAPSplitFile="TC1500264_Gene_SplitFile.txt",Location=
"Output",Name="TC1500264_Pivot")
}
}
