
############################################################
##### PREPROCESSING OF THE AFFYMETRIX ARRAY INFORMATION ####
############################################################

In order for the REIDS function to work an exon annotation file, a transcript annotation file and a junction annotation file are needed.
These can be created by using the "HTA-2_0.na35.hg19.probeset.csv" file and the python scripts included in the package. Run the script either in the folder
which contains "HTA-2_0.na35.hg19.probeset.csv" or alter the python script to include the path to the file.

The first python script will create exon and transcript annotations. The output files are respectively "HTA-2_0_ExonAnnotations.txt", "HTA-2_0_TranscriptAnnotations.txt", "LineIndexing_ExonAnnot.txt" and "LineIndexing_TrAnnot.txt".
The former two files contain the annotations, the latter two files are indexing files which are used in the R code to prevent the loading of the entire files. The script is executed as:

$ python HTA-2_0_ExonAndTranscriptAnnotations.py   

The second python script will create junction associations based on the exon annotations. The output files are respectively "HTA-2_0_JunAnnotations.txt" and "LineIndexing_JAnnot.txt". The script is executed as:

$ python HTA-2_0_JunAssociations.py

The first script will not take long, the second script will (a day). Please feel free to alter the code to speed things up or retrieve junction associations from the exon annotations in the R code directly
(changes are needed in the JunctionAssesment function). 

###############################
###### OTHER MAPPING FILES ####
###############################

Other mappings than Affymetrix' can be used as well if the files have the structure. 

The "HTA-2_0_ExonAnnotations.txt" is:

	TC01000001	JUC01000001	+	EX01051517	3.0
	TC01000001	JUC01000001	+	EX01051518	5.0
	TC01000001	JUC01000002	+	EX01051516	3.0
	TC01000001	JUC01000002	+	EX01051515	5.0
	....
	TC01000001	PSR01000002	+	EX01022607	
	TC01000001	PSR01000002	+	EX01055113	
	TC01000001	PSR01000004	+	EX01022608	
	TC01000001	PSR01000004	+	EX01051517	
	
The "HTA-2_0_TranscriptAnnotations.txt" is:	

	TC01000001	JUC01000001	+	EX01051517	TR01023991
	TC01000001	JUC01000001	+	EX01051518	TR01023991
	TC01000001	JUC01000002	+	EX01051516	TR01023992
	TC01000001	JUC01000002	+	EX01051515	TR01023992
	...
	TC01000001	PSR01000002	+	EX01055113	TR01009585
	TC01000001	PSR01000002	+	EX01022607	TR01021558
	TC01000001	PSR01000002	+	EX01022607	TR01023991
	TC01000001	PSR01000002	+	EX01022607	TR01015438
	TC01000001	PSR01000002	+	EX01022607	TR01023992
	
The "HTA-2_0_JuntAnnotations.txt" is:	

	TC01000205	PSR01003404	JUC01001839	3
	TC01000205	PSR01003405	JUC01001839	5
	TC01000205	PSR01003405	JUC01001844	3
	TC01000205	PSR01003407	JUC01001844	exclusion
	TC01000205	PSR01003413	JUC01001829	3
	TC01000205	PSR01003413	JUC01001833	3
	TC01000205	PSR01003413	JUC01001835	3
	TC01000205	PSR01003413	JUC01001840	3
	TC01000205	PSR01003413	JUC01001842	3
	TC01000205	PSR01003413	JUC01001843	3
	
The exon and transcript annotations can thus easily be replaced by, for example, Ensemble annotations.


###########################
## CREATING THE CDF FILE ##
###########################

The .CDF file necessary to process the array with the R package aroma.affymetrix is created by combining the .clf en .pgf file of the Affymetrix library.
Credit goes to the arome.affymetrix maintainers and writes of the necessary function which I bundle here.

The following steps are required to create the .CDF file:

1) Combine the .pgf and .clf information with "combineProbeInfo.pl":
	
	$ perl combineProbeInfo HTA-2_0.r3
	
	The output consists of two files: "HTA-2_0.r3.probeflat" and "HTA-2_0.r3.psr".

	
2) 	Run "ID_Name_TC.py" to combine probe set ID, probe set name and TC ID. This function retrieves information of "HTA-2_0.na35.hg19.probeset.csv" and the "HTA-2_0.r3.psr" file.

	$ python ID_Name_TC.py
	
	The output is one file: "PSID_Name_TCID.txt".
	
	
3) 	Create an empty file .flat file and run "addGeneId_2.pl" on the "HTA-2_0.r3.probeflat" file and "PSID_Name_TCID.txt" with the empty file as third argument.

	$ perl addGeneId_2 HTA-2_0.r3.probeflat PSID_Name_TCID.txt HTA-2_0.r3.flat
	
	The output will be written to "HTA-2_0.r3.flat".
	
	
4) Run "Flat2CDF.R" in R on the created outfile to obtain the .CDF file.

	The output file is "HTA-2_0,r3.cdf".


(Note: I had to rename "HTA-2_0,r3.cdf" to "HTA-2_0,r.cdf" in order for aroma.affymetrix to recognize the .cdf file and process the .CEL files)

############################
#### PREPROCESSING DONE ####
############################








































	
	
	