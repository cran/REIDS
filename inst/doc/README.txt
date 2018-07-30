
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

	$ python ID_Name_TC.py HTA-2_0.na35.hg19.probeset.csv HTA-2_0.r3.psr PSID_Name_TCID.txt
	
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








































	
	
	