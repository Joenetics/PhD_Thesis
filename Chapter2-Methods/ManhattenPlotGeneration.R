##This code is for generating Manhatten Plots. Gets trait data (etc) and joins it together with SNP data.


##I think I only used these packages...
install.packages("qqman")
install.packages("MASS")
install.packages("lme4")
install.packages("car")
install.packages("ggplot2")
install.packages("ggplot")
install.packages("data.table")
install.packages("backports")
library(data.table)
library (qqman)
library (MASS)
library (lme4)
library (car)
library (ggplot2)

##Format these two lines, then just run the code below to input files
Plates <- c(1,2,3,4,5,6,7)          ##E.g c(1,2,3,4)
PlatesName <- "1-2-3-4-5-6-7"       ##E.g "1-2-3-4"
#PartialOrBinary <- "Partial" 
PartialOrBinary <- "Binary"
#Plates <- c("ZYGO")          ##E.g c(1,2,3,4)
#PlatesName <- "ZYGO"       ##E.g "1-2-3-4"
#Species <- "ZYGO" #not implemented correctly.
Species <- "SacchOnly" #only SC strains
IncludingMasked <- "NoMasked"  # dont put masked SNPs in final output
#IncludingMasked <-  ""  # put masked in final output
#S28COnly <- "S28C_Only"  # this is only activated for kinship calculations, but not for actual linear regressions
S28COnly <- ""
#SCset <- "Full"
#SCset <- "Reduced"
SCset <- "DoubleReduced"
cigar_stuff <- "_CIGAR" # using CIGAR stuff
#cigar_stuff <- ""  # normal

#CoreOnly = "CoreOnly"
CoreOnly = ""

##For all strains, YNB metabolites
Plates <- c(1,2,3,4,5,6,7)          ##E.g c(1,2,3,4)
PlatesName <- "1-2-3-4-5-6-7"       ##E.g "1-2-3-4"
#PartialOrBinary <- "Partial" 
PartialOrBinary <- "Binary"
#Plates <- c("ZYGO")          ##E.g c(1,2,3,4)
#PlatesName <- "ZYGO"       ##E.g "1-2-3-4"
#Species <- "ZYGO" #not implemented correctly.
#Species <- "SacchOnly" #only SC strains
Species <- "" # all strains
IncludingMasked <- "NoMasked"  # dont put masked SNPs in final output
#IncludingMasked <-  ""  # put masked in final output
S28COnly <- "S28C_Only"  # this is only activated for kinship calculations, but not for actual linear regressions
#S28COnly <- ""
#SCset <- "Full"
#SCset <- "Reduced"
#SCset <- "DoubleReduced"
SCset <- ""
cigar_stuff <- "_CIGAR" # using CIGAR stuff
#cigar_stuff <- ""  # normal

#CoreOnly = "CoreOnly"
CoreOnly = ""


##
# Use these if you want to use SC with DE appended to the end.
Plates <- c(1,2,3,4,5,6,7, "DE")          ##old= c("DE")
PlatesName <- "1-2-3-4-5-6-7-DE"       ##old= "DE"
#PartialOrBinary <- "Partial" 
PartialOrBinary <- "Binary"
Species <- "" #if nothing. Or 'All'
IncludingMasked <- "NoMasked"  # dont put masked SNPs in final output
#IncludingMasked <-  ""  # put masked in final output
S28COnly <- ""
SCset <- ""
#CoreOnly = "CoreOnly"
CoreOnly = ""
cigar_stuff <- "_CIGAR" # using CIGAR stuff
#cigar_stuff <- ""  # normal



base_directory = "C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\"


if (PlatesName == "1-2-3-4-5-6-7-DE"){
  print("All plates and DE...")
  SNPs <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joined_MAF_files\\Linear_Regression_", PartialOrBinary,"_SNPs_MappedToReference_Plates_", PlatesName, Species, IncludingMasked, S28COnly,SCset, cigar_stuff, ".txt", sep =''), sep = ' ', header= FALSE)
  print("DE SNPs done...")
  Letter_SNPs <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joined_MAF_files\\ATGC_File_AllStrains_Plates_", PlatesName, Species, IncludingMasked, S28COnly,SCset, cigar_stuff, ".csv", sep= ''), header=FALSE)
  print("DE Letter SNPs done...")
  #TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\NCYC_Excel\\PythonStrainFinder\\WekaPredictedResistanceLogGrowthSacchOnlyDoubleReduced666.csv", header= TRUE)  # plates 1-9, not DE
  TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\NCYC_Excel\\PythonStrainFinder\\WekaPredictedResistanceLogGrowthDE666.csv", header= TRUE)  # DE growth!!
  #TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\NCYC_Excel\\PythonStrainFinder\\JoinedSlopeResultsLogGrowthBinaryResistance_8OrGreater.csv", header= TRUE)
  #TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Josephh_NMR_TraitData.csv", header= TRUE, row.names = 1)
  print("DE Trait Data done...")
  GeneNames <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joined_MAF_files\\MAFfile_HighQuality_SNPsOnly_Plates_", PlatesName, Species, IncludingMasked,S28COnly, SCset, cigar_stuff, ".csv", sep = ''), header= FALSE)
  print("DE Gene Names done...")
  SelectListOfStrains <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joined_MAF_files\\OnlyPlates_", PlatesName,"_OnlySpecies_", Species,".csv", sep =''), header = TRUE)
  print("DE Select List Of Strains done...")
  Qmatrix <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_CoreOnly", Species,".Q", sep =''), sep = " ", header = FALSE, row.names = 1)
  print("DE QMatrix done...")
} else if (length(Plates) ==1){
  print("Single plate analysis")
  SNPs <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Plate", PlatesName, "VCF", cigar_stuff, "_AfterR\\", CoreOnly, "Linear_Regression_",PartialOrBinary,"_SNPs_MappedToReference_Plate", PlatesName, IncludingMasked, S28COnly, SCset, cigar_stuff, ".txt", sep = ''), sep = ' ', header= FALSE)
  print("SNPs done...")
  Letter_SNPs <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Plate", PlatesName, "VCF", cigar_stuff, "_AfterR\\", CoreOnly, "ATGC_File_AllStrains_Plate", PlatesName, IncludingMasked, S28COnly,SCset, cigar_stuff, ".csv", sep = ''), header=FALSE)
  print("Letter SNPs done...")
  #TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\NCYC Excel\\PythonStrainFinder\\All9PlatesForWeka.csv", header= TRUE)
  TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\NCYC_Excel\\PythonStrainFinder\\WekaPredictedResistanceLogGrowthDE666.csv", header= TRUE)  # DE growth!!
  print("TraitData done...")
  GeneNames <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Plate", PlatesName, "VCF", cigar_stuff, "_AfterR\\", CoreOnly, "MAFfile_HighQuality_SNPsOnly_Plate", PlatesName, IncludingMasked, S28COnly, SCset, cigar_stuff, ".csv", sep = ''), header= FALSE)
  print("Gene Names done...")
  #Qmatrix <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Plate", PlatesName, "VCF_AfterR\\QMatrixForPlates_", PlatesName,".Q", sep =''), sep = " ", row.names = 1, header = FALSE)
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Josephh_SNP.Q", sep = " ", header = FALSE)
  Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_NoMasked_DEonly.Q", sep = " ", header = FALSE, row.names = 1) # for DE?
  #Qmatrix <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_CoreOnly", Species,".Q", sep =''), sep = " ", header = FALSE, row.names = 1)
  print("Single Q-Matrix read in...")
  }else {
  print("Starting Multi-plate...")
  SNPs <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joined_MAF_files\\Linear_Regression_", PartialOrBinary,"_SNPs_MappedToReference_Plates_", PlatesName, Species, IncludingMasked, S28COnly,SCset, cigar_stuff, ".txt", sep =''), sep = ' ', header= FALSE)
  print("Multi-plate SNPs done...")
  Letter_SNPs <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joined_MAF_files\\ATGC_File_AllStrains_Plates_", PlatesName, Species, IncludingMasked, S28COnly,SCset, cigar_stuff, ".csv", sep= ''), header=FALSE)
  print("Multi-plate Letter SNPs done...")
  TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\NCYC_Excel\\PythonStrainFinder\\WekaPredictedResistanceLogGrowthSacchOnlyDoubleReduced666.csv", header= TRUE)  # plates 1-9, not DE. BEST FURFURAL ONE
  #CoreOnlyATGC_File_AllStrains_Plates_1-2-3-4-5-6-7SacchOnlyNoMaskedDoubleReduced_ATGCClusters_JCD5_3branches.q
  #TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\RStudio_YNB_Metabolites_Updated.csv", header= TRUE, row.names = 1)  # plates 1-9 for all strains. YNB NMR metabolite data. NOT FURFURAL
  #TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\RStudio_Malt_Metabolites.csv", header= TRUE, row.names = 1)  # plates 1-9 for all strains. Malt NMR metabolite data. NOT FURFURAL
  #TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\NCYC_Excel\\PythonStrainFinder\\JoinedSlopeResultsLogGrowthBinaryResistance_8OrGreater.csv", header= TRUE)
  #TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Josephh_NMR_TraitData.csv", header= TRUE, row.names = 1)
  print("Multi-plate Trait Data done...")
  GeneNames <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joined_MAF_files\\MAFfile_HighQuality_SNPsOnly_Plates_", PlatesName, Species, IncludingMasked,S28COnly, SCset, cigar_stuff, ".csv", sep = ''), header= FALSE)
  print("Multi-plate Gene Names done...")
  #SelectListOfStrains <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joined_MAF_files\\OnlyPlates_", PlatesName,"_OnlySpecies_", Species,".csv", sep =''), header = TRUE)
  print("Multi-plate Select List Of Strains done...")
  #Qmatrix <- read.csv(paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_CoreOnly", Species,".Q", sep =''), sep = " ", header = FALSE, row.names = 1)
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\SANE_QMatrix_Cigar_CoreOnly_SacchOnly_TamD.q", sep = "\t", header = FALSE, row.names = 1)
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_NoMasked_SacchOnly.Q", sep = " ", header = FALSE, row.names = 1) #
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_CoreOnly_NoMasked_SacchOnly.Q", sep = " ", header = FALSE, row.names = 1) #.
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_NoMasked_SacchOnly.Q", sep = " ", header = FALSE, row.names = 1) #
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_CoreOnly_NoMasked_SacchOnly.Q", sep = " ", header = FALSE, row.names = 1) # DECENT 17/09/2021
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_CoreOnly_SacchOnly.Q", sep = " ", header = FALSE, row.names = 1) # DECENT 17/09/2021
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_NoMasked_SacchOnly.Q", sep = " ", header = FALSE, row.names = 1)# DECENT 17/09/2021
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_NoMasked_SacchOnly.Q", sep = " ", header = FALSE, row.names = 1)# DECENT 17/09/2021
  Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_NoMasked.Q", sep = " ", header = FALSE, row.names = 1)# All strains! 17/09/2021
  Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\SANE_QMatrix_Cigar_NoMasked_K2P.Q", sep = "\t", header = FALSE, row.names = 1)# All strains! 17/09/2021
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_SacchOnly.Q", sep = " ", header = FALSE, row.names = 1)# DECENT 17/09/2021
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\PSIKO_QMatrix_Cigar_CoreOnly_NoMasked_SacchOnly.Q", sep = " ", header = FALSE, row.names = 1)  # DECENT 18/09/2021
  #Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\All_QMatrices\\SANE_QMatrix_Cigar_CoreOnly_NoMasked_SacchOnly_TamD.Q", sep = "\t", header = FALSE, row.names = 1)  # own TAMD distance (best distribution on mantel test figure.) 23/09/2021
  
  print("QMatrix done...")
  
}



kai_data <- read.csv("Kai_Metabolites.csv", header = TRUE)
row.names(kai_data) <- kai_data$Strain

Strain_list <- SelectListOfStrains$NCYC.Number
for (i in Strain_list){
  print(i)
}


#Run this first, always.

  Genes <- data.frame(unique(GeneNames$V3))
  row.names(SelectListOfStrains) <- SelectListOfStrains$NCYC.Number
  #row.names(TraitData) <- TraitData$StrainNumber
  row.names(TraitData) <- TraitData$Strain.Number  # do not use this line when doing NMR data 
  GeneNames <- data.frame(GeneNames)
  #rownames(Qmatrix) <- Qmatrix$V1    # should no longer be necessary, as import rownames from csv
  #Qmatrix <- Qmatrix[,-1]   #remove strain names.
  
  #run this second, always
  ##This makes a DF for the SNPs and removes strain data, making them column names.
  SNPs <- t(SNPs)
  dfSNPs <- data.frame(SNPs)
  colnames(dfSNPs) <- dfSNPs[1,]
  dfSNPs <- dfSNPs[-1,]
  row.names(dfSNPs) <- (GeneNames$V1)


#Run this third, always. (then skip to ###ONE or ###TWO to start.)
##This is for letter SNPs, to match them later to p-value things
Letter_SNPs <- t(Letter_SNPs)
df_Letter_SNPs <- data.frame(Letter_SNPs)
colnames(df_Letter_SNPs) <- colnames(dfSNPs)
df_Letter_SNPs <- df_Letter_SNPs[-1,]
row.names(df_Letter_SNPs) <- (GeneNames$V1)


#length(row.names(df_Letter_SNPs))
#length(GeneNames$V1)
##Trait data part. NEED to run all of this part.
#for (i in c(611,2435,2474,2516,2628,2677)){##Remove strains that failed genome sequencing in plate 7 
#  TraitData <- TraitData[row.names(TraitData) != i,]  
#}

#Only next line if doing Kai malt data stuff
Kai_data1 <- kai_data[row.names(TraitData) == "DeleteAll"]

TraitData1 <- TraitData[row.names(TraitData) == "DeleteAll",]
dfSNPs2 <- dfSNPs[,colnames(dfSNPs) == "DeleteAll"]
df_Letter_SNPs2 <- df_Letter_SNPs[,colnames(df_Letter_SNPs) == "DeleteAll"]
Qmatrix2 <- Qmatrix[row.names(Qmatrix) == "DeleteAll",]
SelectListOfStrains2 <- SelectListOfStrains[row.names(SelectListOfStrains) == "DeleteAll",]


for (j in c(colnames(dfSNPs))){             ##Include strains in new df that are present in dfSNPs
    #print(j)
    if ( j %in% c(colnames(dfSNPs)) && j %in% row.names(TraitData) && j %in% row.names(Qmatrix)){
    #if ( j %in% c(colnames(dfSNPs)) && j %in% row.names(kai_data) && j %in% row.names(Qmatrix)){
    #Only next line if doing Kai data , AND line above.
    #Kai_data1 <- rbind(Kai_data1, kai_data[row.names(kai_data) == j,] )
    TraitData1 <- rbind(TraitData1, TraitData[row.names(TraitData) == j,] )
    dfSNPs2 <- cbind(dfSNPs2, dfSNPs[,colnames(dfSNPs) ==  j] )
    df_Letter_SNPs2 <- cbind(df_Letter_SNPs2, df_Letter_SNPs[,colnames(df_Letter_SNPs) ==  j] )
    Qmatrix2 <- rbind(Qmatrix2, Qmatrix[row.names(Qmatrix) == j,] )
    SelectListOfStrains2 <- rbind(SelectListOfStrains2, SelectListOfStrains[row.names(SelectListOfStrains) == j,] )
    
  }


}
colnames(dfSNPs2) <- rownames(TraitData1)##bring back names to columns of dfSNPs2
colnames(df_Letter_SNPs2) <- rownames(TraitData1)##bring back names to columns of dfSNPs2


dfSNPs2 <- lapply(dfSNPs2, as.integer)  # for some reason, DE becomes character. need it as int.
dfSNPs2 <- as.data.frame(dfSNPs2)  # for some reason, DE becomes character. need it as int.
colnames(dfSNPs2) <- rownames(TraitData1)##bring back names to columns of dfSNPs2


#This is for Furfural MaxOD data.
Resistances <- TraitData1$Resistance.Points..1.10.  # Resistances are o longer MaxOD; here, they are resistance from 1-15.
Resistances <- TraitData1$Inflection.Cluster.Score  # ]
Resistances <- TraitData1$Inflection.Cluster.Score
Resistances <- TraitData1$Resistance.Points..1.10.
Resistances <- TraitData1$Slope.Cluster.Score + TraitData1$Slope.Point.Cluster.Score # this one for inflection point 
Resistances <- TraitData1$Inflection.Cluster.Score + TraitData1$Slope.Cluster.Score + TraitData1$Slope.Point.Cluster.Score
Resistances <- TraitData1$Inflection.Cluster.Score + TraitData1$Slope.Point.Cluster.Score
Resistances <- TraitData1$Inflection.Cluster.Score+ TraitData1$Slope.Cluster.Score
Resistances <- TraitData1$Inflection.Cluster.Score + TraitData1$MaxOD.Cluster.Score + TraitData1$C_value.Cluster.Score + TraitData1$Slope.Point.Cluster.Score



Resistances <- TraitData1$Inflection.Cluster.Score + TraitData1$Slope.Cluster.Score  # one at -6
Resistances <- TraitData1$Inflection.Cluster.Score + TraitData1$C_value.Cluster.Score + TraitData1$Slope.Cluster.Score + TraitData1$MaxOD.Cluster.Score + TraitData1$Slope.Point.Cluster.Score


Resistances <- TraitData1$Rough.Inflection.Cluster.Score + TraitData1$MaxOD.Cluster.Score + TraitData1$Slope.Cluster.Score  
Resistances <- TraitData1$Inflection.Cluster.Score + TraitData1$MaxOD.Cluster.Score + TraitData1$Slope.Point.Cluster.Score # best!

Resistances <- TraitData1$Slope.Point.Cluster.Score
Resistances <- TraitData1$Succinate

colnames(TraitData1)
Resistances <- TraitData1$Ethanol
Resistances <- TraitData1$Succinate
Resistances <- TraitData1$Methanol
Resistances <- TraitData1$Glycerol
Resistances <- TraitData1$Glucose
Resistances <- TraitData1$Acetoin

Resistances <- as.numeric(as.character(Kai_data1$Ethanol))
Resistances <- as.numeric(as.character(Kai_data1$Butyrate))
Resistances <- as.numeric(as.character(Kai_data1$Acetoin))
Resistances <- as.numeric(as.character(Kai_data1$Succinate))



###ONE 
#This is code for all strains in dataset

###This is code for making GWAS data frame!!! Each SNP gains a P-value, for later addition into manhatten plot.
PValues <- array(0,dim(dfSNPs2)[1])  #
LmerPValues <- array (0, dim(dfSNPs2)[1])
LmerPValues_positivenegative <- array (0, dim(dfSNPs2)[1])
logLmerPValues <- array (0, dim(dfSNPs2)[1])
Positive_Negative_Correlation <- array (0, dim(dfSNPs2)[1])



#Only Linear- no Qmatrix LMM
for (i in c(1:dim(dfSNPs2)[1])){
  if (0 %in% unlist(dfSNPs2[i,]) && 1 %in% unlist(dfSNPs2[i,])){  # test there are both 0 and 1 in dataset... in small ones, you can have solely mutant allele and lm fails.
   
  Linears <- lm(Resistances~ unlist(dfSNPs2[i,]))   ##Remember, resistances MUST == strain numbers! (multiples of 96 for wells.) Then, LR from 0/1 to OD is better!
  if (is.na(summary(Linears)$coefficients[8])){
    PValues[i] <- 1
  }
  else {
    PValues[i] <- summary(Linears)$coefficients[8]
  }
  if (is.na((Linears)$coefficients[2])){
    Positive_Negative_Correlation[i] <- "NA"
  }
  else if (as.double((Linears)$coefficients[2]) >= 0){
    Positive_Negative_Correlation[i] <- "+"
  }
  else {
    Positive_Negative_Correlation[i] <- "-"
  }
  }
  else{
    PValues[i] <- 1  # only alt allele in dataset. cannot do LM.
    Positive_Negative_Correlation[i] <- "NA"
  }

}


for (i in c(1:dim(dfSNPs2)[1])){
 #= logLmerians <- lmer(logg ~ unlist(dfSNPs2[i,])+ (1 | Qmatrix2[,1]) + (1 | Qmatrix2[,2]) + (1 | Qmatrix2[,3]) , REML = FALSE, options(warn = -1))
  if (0 %in% unlist(dfSNPs2[i,]) && 1 %in% unlist(dfSNPs2[i,])){
  Lmerians <- lmer(Resistances ~ unlist(dfSNPs2[i,])+ (1 | Qmatrix2[,1]) + (1 | Qmatrix2[,2]) + (1 | Qmatrix2[,3]) , REML = FALSE, options(warn = -1))
  Linears <- lm(Resistances~ unlist(dfSNPs2[i,]))   ##Remember, resistances MUST == strain numbers! (multiples of 96 for wells.) Then, LR from 0/1 to OD is better!
  if (is.na(summary(Linears)$coefficients[8])){
    PValues[i] <- 1
  }
  else {
    PValues[i] <- summary(Linears)$coefficients[8]
  }
  if (is.na(Anova(Lmerians)$Pr[1])){
  LmerPValues[i]<- 1
  }
  else {
  LmerPValues[i] <- (Anova(Lmerians)$Pr)[1]
  }
  if (is.na((Linears)$coefficients[2])){
    Positive_Negative_Correlation[i] <- "NA"
  }
  else if (as.double((Linears)$coefficients[2]) >= 0){
    Positive_Negative_Correlation[i] <- "+"
  }
  else {
    Positive_Negative_Correlation[i] <- "-"
  }
  #if (is.na(Anova(logLmerians)$Pr[1])){
  #  logLmerPValues[i]<- 1
  #}
  #else {
  #  logLmerPValues[i] <- (Anova(Lmerians)$Pr)[1]
  #}
  }
  else{  # cannot make correlation if all strains are 0 or 1
    PValues[i] <- 1
    LmerPValues[i]<- 1
    Positive_Negative_Correlation[i] <- "NA"
  }
}


#Initise MyORFs (just numbers, use ORF numbers as fake chromosomes)
   #genomewideline = -log10(5e-08)
manhattan(MyLmerGWAS, col=c("blue4","red3", "green2", "orange3"), chrlabs= ORFnames,  suggestiveline = FALSE, genomewideline = FALSE) #genomewideline = -log10(5e-08)



##Without QMATRIX, With dfSNPs2
for (i in c(1:dim(dfSNPs2)[1])){
  Linears <- lm(Succinate~ unlist(dfSNPs2[i,]))   ##Remember, resistances MUST == strain numbers! (multiples of 96 for wells.) Then, LR from 0/1 to OD is better!
  if (is.na(summary(Linears)$coefficients[8])){
    PValues[i] <- 1
  }
  else {
    PValues[i] <- summary(Linears)$coefficients[8]
  }
  
}


#printing names of genes with selected p-values in list

for (i in c(1:length(PValues))){
  if (grepl( "YKL071W", GeneNames[i,1], fixed = TRUE)){
    print(paste("Name of Gene:", GeneNames[i,1], "PValue:", PValues[i]))
  }
  if (PValues[i] < 0.0001 ){
    ACount = 0
    GCount = 0
    CCount = 0
    TCount = 0
    
    print(paste("Name of Gene:", GeneNames[i,1], "PValue:", PValues[i]))
    #print(Resistances)
    #print(df_Letter_SNPs[i,],  row.names = FALSE)
    for (f in c(1:length(df_Letter_SNPs2[i,]))){
      if (df_Letter_SNPs2[i,][f] == "A"){
        ACount <- ACount +1
      }
      else if (df_Letter_SNPs2[i,][f] == "G"){
        GCount <- GCount +1
      }
      else if (df_Letter_SNPs2[i,][f] == "C"){
        CCount <- CCount +1
      }
      else if (df_Letter_SNPs2[i,][f] == "T"){
        TCount <- TCount +1
      }
      else{
        print(paste("what: ", df_Letter_SNPs2[i,][f]))
      }
    }
    print(paste("A: ", toString(ACount), "T: ", toString(TCount), "G: ", toString(GCount), "C: ", toString(CCount)), sep ='')
  }
}


#printing names of genes with selected lmer-p-values in list
for (i in c(1:length(LmerPValues))){
  
  if (LmerPValues[i] < 0.00001 ){
    ACount = 0
    GCount = 0
    CCount = 0
    TCount = 0
    print(paste("Name of Gene:", GeneNames[i,1], "PValue:", LmerPValues[i]))
    #print(Resistances)
    #print(df_Letter_SNPs[i,],  row.names = FALSE)
    for (f in c(1:length(df_Letter_SNPs2[i,]))){
      if (df_Letter_SNPs2[i,][f] == "A"){
        ACount <- ACount +1
      }
      else if (df_Letter_SNPs2[i,][f] == "G"){
        GCount <- GCount +1
      }
      else if (df_Letter_SNPs2[i,][f] == "C"){
        CCount <- CCount +1
      }
      else if (df_Letter_SNPs2[i,][f] == "T"){
        TCount <- TCount +1
      }
      else{
        print(paste("what: ", df_Letter_SNPs2[i,][f]))
      }
    }
    print(paste("A: ", toString(ACount), "T: ", toString(TCount), "G: ", toString(GCount), "C: ", toString(CCount)), sep ='')
  }
}



###FOR DE- get top SNP hits from non-DE file. 
#here : JOINED_MAF_Files\lmerPValues_Plates_1-2-3-4-5-6-7_SacchOnlyNoMaskedDoubleReduced666infpoints_SlopePoint_maxOD_withQ.csv
df_Letter_SNPs2["7825-Q0105_NumOfGenes_4:6255_C/", ]  # example
SNPs_of_interest <- data.frame(read.csv("Top1000Genes.csv", header = FALSE))
comparison_table <- data.frame(comparison_table, stringsAsFactors = TRUE)

temp_vector_snp <- c(1:length(SNPs_of_interest[,1]))
  for (i in c(1:length(SNPs_of_interest[,1]))){
  
   current_SNP <- SNPs_of_interest[i,]
  
    ACount = 0
    GCount = 0
    CCount = 0
    TCount = 0
    
    snp_list <- df_Letter_SNPs2[current_SNP,]
    snp_list <- paste(shQuote(df_Letter_SNPs2[current_SNP,]), collapse=",")
    snp_list <- gsub("\"", "", snp_list)
    #print(paste("Name of Gene:", current_SNP))
    temp_vector_snp[i] <- snp_list
    ##print(paste("SNP List:", as.vector(snp_list)))
    #for (f in c(1:length(df_Letter_SNPs2[i,]))){
    #  if (df_Letter_SNPs2[i,][f] == "A"){
    #    ACount <- ACount +1
    #  }
    #  else if (df_Letter_SNPs2[i,][f] == "G"){
    #    GCount <- GCount +1
    #  }
    #  else if (df_Letter_SNPs2[i,][f] == "C"){
    #    CCount <- CCount +1
    #  }
    #  else if (df_Letter_SNPs2[i,][f] == "T"){
    #    TCount <- TCount +1
    #  }
    #  else{
    #    print(paste("what: ", df_Letter_SNPs2[i,][f]))
    #  }
    #}
    #print(paste("A: ", toString(ACount), "T: ", toString(TCount), "G: ", toString(GCount), "C: ", toString(CCount)), sep ='')
  
}
comparison_table <- data.frame(temp_vector_snp)
row.names(comparison_table) <- SNPs_of_interest[,1]
stains_list <- paste(shQuote(row.names(TraitData1)), collapse=",")
stains_list <- gsub("\"", "", stains_list)
colnames(comparison_table) <- stains_list

write.table(comparison_table, paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joined_MAF_files\\DE_SNPs_For_Top_1001Hits.csv", sep = ''), sep = ",", col.names = TRUE, row.names = TRUE)


###FOR DE


for (i in c(1:length(Resistances))){
  print(Resistances[i], row.names = FALSE)
}

#Initise MyORFs (just numbers, use ORF numbers as fake chromosomes)
MyORFs <- c(1:length(GeneNames$V2))
#Make unique list of ORF names, for chrlabel bit later
ORFnames <- unique(GeneNames$V3)
#This vector of whole gene name must be a character for Manhatten plot.
GeneNames$V3 <- as.character(GeneNames$V3)
#Make DF for manhatten plot
MyGWAS <- data.frame(SNP= GeneNames$V1, CHR= MyORFs, BP= c(1:length(GeneNames$V1)), P= PValues)
MyLmerGWAS <- data.frame(SNP= GeneNames$V1, CHR= MyORFs, BP= c(1:length(GeneNames$V1)), P= LmerPValues)
MylogLmerGWAS <- data.frame(SNP= GeneNames$V1, CHR= MyORFs, BP= c(1:length(GeneNames$V1)), P= logLmerPValues)
#Run Manhatten Plot
manhattan(MyGWAS, col=c("blue4","red3", "green2", "orange3"), chrlabs= ORFnames,  suggestiveline = FALSE, genomewideline = FALSE) #genomewideline = -log10(5e-08)
manhattan(MyLmerGWAS, col=c("blue4","red3", "green2", "orange3"), chrlabs= ORFnames,  suggestiveline = FALSE, genomewideline = FALSE) #genomewideline = -log10(5e-08)
manhattan(MylogLmerGWAS, col=c("blue4","red3", "green2", "orange3"), chrlabs= ORFnames,  suggestiveline = FALSE, genomewideline = FALSE) #genomewideline = -log10(5e-08)


manhattan(MyGWAS, col=c("blue4"), chrlabs= ORFnames,  suggestiveline = FALSE, genomewideline = FALSE) #genomewideline = -log10(5e-08)
manhattan(MyLmerGWAS, col=c("red3"), chrlabs= ORFnames,  suggestiveline = FALSE, genomewideline = FALSE) #genomewideline = -log10(5e-08)



#Random error, one value was 0 in pvalues, so -log10(p) was infinity... this fixes that.
logLmerPValues[logLmerPValues == 0] <- 0.00001
for (i in c(1: length(LmerPValues))){
  if (LmerPValues[i] == 0.0003){
    print(LmerPValues[i])
  }
}

##Printing Stuff to 'Studio_Hits' file. Then, use RStudio_Hits_mapped_To_Yeast_genome.py to map hits to genome.
CounterForHits <- 0
PValueThreshold <- 0.001         #Modify this for different threshold on hits file.
for (i in c(1:length(PValues))){
  
  if (PValues[i] < PValueThreshold ){
    
    CounterForHits <- CounterForHits + 1
  }
  
}
OutputCSV <- as.data.frame(cbind(c(1:CounterForHits), c(1:CounterForHits)))
colnames(OutputCSV) <- c("GeneNames", "Pvalues")
OutputNumber <- 1
for (i in c(1:length(PValues))){
  
  if (PValues[i] < PValueThreshold ){
    
    OutputCSV[OutputNumber,1] <- GeneNames[i,3]
    OutputCSV[OutputNumber,2] <- PValues[i]
    OutputNumber <- OutputNumber +1
    
  }
  
}
write.csv(OutputCSV, file = paste("Succinate_96Sacch_NOQ_RegularRegression", "_RStudioHits.csv"))





################################################################
###TWO
#This is code for reduced dataset (e.g, only sacch, etc)

##List Strains
colnames(SelectListOfStrains2) <- c("Index", "NCYC", "Strain", "Resistant?", "Binary Resistance")
SelectListOfStrainsNames <- SelectListOfStrains2$NCYC ##Only NCYC number list
SelectListOfStrainsIndex <- SelectListOfStrains2$Index ##Only index of strains

#Resistances becomes the final OD of each strain (in furfural tests). For NMR data, skip this bit. 
ResistancesForSelectStrains <- TraitData1[c(SelectListOfStrainsIndex),7]
##Select List of strains SNP df
SelectSNPs2 <- dfSNPs[,SelectListOfStrainsIndex]
Select_Letter_SNPs <- df_Letter_SNPs[, SelectListOfStrainsIndex]

##put here a way to mix SC & NMR columns only. Start by creating empty DFs
SelectTraitData <- TraitData[row.names(TraitData) == "DeleteAll",]
SelectSNPs <- dfSNPs[,colnames(dfSNPs) == "DeleteAll"]
Select_Letter_SNPs <- df_Letter_SNPs[,colnames(dfSNPs) == "DeleteAll"]

for (j in c(colnames(SelectSNPs2))){             ##Include strains in new df that are present in dfSNPs
  
  if ( j %in% c(colnames(SelectSNPs2)) && j %in% row.names(TraitData1)){  ##TraitData1 has successful NMR strains; SelectSNPs2 has SC (or other) strainss.
    SelectTraitData <- rbind(SelectTraitData, TraitData[row.names(TraitData) == j,] )
    SelectSNPs <- cbind(SelectSNPs, dfSNPs[,colnames(dfSNPs) ==  j] )
    Select_Letter_SNPs <- cbind(Select_Letter_SNPs, df_Letter_SNPs[,colnames(df_Letter_SNPs) ==  j] )
    
  }
  
}
colnames(SelectSNPs) <- rownames(SelectTraitData)##bring back names to columns of dfSNPs2
colnames(Select_Letter_SNPs) <- rownames(SelectTraitData)##bring back names to columns of dfSNPs2


###This is code for making Sacch-only GWAS data frame!!! Each SNP gains a P-value, for later addition into manhatten plot.
PValuesListed <- array(0,dim(SelectSNPs)[1])
for (i in c(1:dim(SelectSNPs)[1])){
  Linears <- lm(SelectTraitData$Succinate~ unlist(SelectSNPs[i,]))
  if (is.na(summary(Linears)$coefficients[8])){FileForPythonThenWeka <- data.frame(GeneNames$V1, LmerPValues )

FileForPythonThenWeka <- data.frame(GeneNames$V1, LmerPValues )
colnames(FileForPythonThenWeka) <- c("GeneNames", "Pvalue_of_Gene")
write.table(FileForPythonThenWeka, paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\PValue_Files\\lmerPValues_Plates_", PlatesName,"_",Species, IncludingMasked,S28COnly,SCset, cigar_stuff,  "_MaltEthanol.csv", sep = ''), sep = ",", col.names = FALSE, row.names = FALSE)
FileForPythonThenWeka <- data.frame(GeneNames$V1, PValues)
colnames(FileForPythonThenWeka) <- c("GeneNames", "Pvalue_of_Gene")
write.table(FileForPythonThenWeka, paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\PValue_Files\\PValues_Plates_", PlatesName,"_", PartialOrBinary, Species, IncludingMasked, S28COnly,SCset, cigar_stuff, "_MaltEthanol.csv", sep = ''), sep = ",", col.names = FALSE, row.names = FALSE)

FileForPythonThenWeka <- data.frame(GeneNames$V1, LmerPValues, Positive_Negative_Correlation)
colnames(FileForPythonThenWeka) <- c("GeneNames", "Pvalue_of_Gene", "Positive_Or_Negative_Correlation")
write.table(FileForPythonThenWeka, paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\PValue_Files\\lmerPValues_Plates_", PlatesName,"_",Species, IncludingMasked,S28COnly,SCset, cigar_stuff,  "Glycerol_with_PSIKO_Cigar_NoMasked.csv", sep = ''), sep = ",", col.names = FALSE, row.names = FALSE)

write.table(FileForPythonThenWeka, paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\PValue_Files\\lmerPValues_Plates_", PlatesName,"_",Species, IncludingMasked,S28COnly,SCset, cigar_stuff,  "666infpoints_SlopePoint_maxOD_with_PSIKO_DE_NoMasked.csv", sep = ''), sep = ",", col.names = FALSE, row.names = FALSE)
FileForPythonThenWeka <- data.frame(GeneNames$V1, PValues, Positive_Negative_Correlation)
colnames(FileForPythonThenWeka) <- c("GeneNames", "Pvalue_of_Gene", "Positive_Or_Negative_Correlation")
write.table(FileForPythonThenWeka, paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\PValue_Files\\PValues_Plates_", PlatesName,"_", PartialOrBinary, Species, IncludingMasked, S28COnly,SCset, cigar_stuff, "666infpoints_SlopePoint_maxOD_NoQ(DE).csv", sep = ''), sep = ",", col.names = FALSE, row.names = FALSE)

write.table(FileForPythonThenWeka, paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\PValue_Files\\PValues_Plates_", PlatesName,"_", PartialOrBinary, Species, IncludingMasked, S28COnly,SCset, cigar_stuff, "666infpoints_SlopePoint_maxOD_NoQ.csv", sep = ''), sep = ",", col.names = FALSE, row.names = FALSE)
#below is temp for kai
w#rite.table(FileForPythonThenWeka, paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\PValue_Files\\PValues_Plates_Succinate_Kai.csv", sep = ''), sep = ",", col.names = FALSE, row.names = FALSE)

    PValuesListed[i] <- 1
  }
  else {
    PValuesListed[i] <- summary(Linears)$coefficients[8]
  }
}


#This array will be printed to csv, then used by PValuesForML.py to make .arff file


SameAsTopForOD <- data.frame(row.names(TraitData1), TraitData1$Binary.Resistance..1..resistant.)
colnames(SameAsTopForOD) <- c("Strain", "Strain_Resistance")
write.table(SameAsTopForOD, paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joined_MAF_files\\Strain_and_Resistance_", PlatesName, "_", Species, "12.csv", sep = ''), sep = ",", col.names = FALSE, row.names = FALSE)


#printing names of genes with selected p-values in list
for (i in c(1:length(PValuesListed))){
  
  if (PValuesListed[i] < 0.00001 ){
    ACount = 0
    GCount = 0
    CCount = 0
    TCount = 0
    print(paste("Name of Gene:", GeneNames[i,1], "PValue:", PValuesListed[i]))
    #print(Resistances)
    #print(df_Letter_SNPs[i,],  row.names = FALSE)
    for (f in c(1:length(Select_Letter_SNPs[i,]))){
      if (Select_Letter_SNPs[i,][f] == "A"){
        ACount <- ACount +1
      }
      else if (Select_Letter_SNPs[i,][f] == "G"){
        GCount <- GCount +1
      }
      else if (Select_Letter_SNPs[i,][f] == "C"){
        CCount <- CCount +1
      }
      else if (Select_Letter_SNPs[i,][f] == "T"){
        TCount <- TCount +1
      }
      else{
        print(paste("what: ", Select_Letter_SNPs[i,][f]))
      }
    }
    print(paste("A: ", toString(ACount), "G: ", toString(GCount), "C: ", toString(CCount), "T: ", toString(TCount)), sep ='')
  }
}

for (i in c(1:length(SelectTraitData$Succinate))){
  print(SelectTraitData$Succinate[i], row.names = FALSE)
}

#Initise MyORFs (just numbers, use ORF numbers as fake chromosomes)
MyORFs <- GeneNames$V2
#Make unique list of ORF names, for chrlabel bit later
ORFnames <- unique(GeneNames$V3)
#This vector of whole gene name must be a character for Manhatten plot.
GeneNames$V3 <- as.character(GeneNames$V3)
#Make DF for manhatten plot
MyGWAS <- data.frame(SNP= GeneNames$V1, CHR= MyORFs, BP= c(1:length(GeneNames$V1)), P= PValuesListed)
#Run Manhatten Plot
manhattan(MyGWAS, col=c("blue4","red3", "green2", "orange3"), chrlabs= ORFnames,  suggestiveline = FALSE, genomewideline = -log10(1e-11))



##small histogram images. 'Characteristic frequency by K-means cluster size'
input_kmeans_data <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\K_Means_Changes_To_Characteristics.csv", header = TRUE)
input_kmeans_data <- data.frame(input_kmeans_data)

titles <- c("Inflection Point", "MaxOD", "Timepoint of ??")



x11()
par(mfrow=c(length(input_kmeans_data)/3,1)) 
par(mar=c(3,4,4,1))
for (column_n in c(1:(length(input_kmeans_data)/3))){

  little_DF <- as.matrix(t(input_kmeans_data[,c(column_n, column_n + 3, column_n + 6)]))
  barplot(little_DF, beside = TRUE, ylim = c(0,50),
          main = titles[column_n],
          xlab = "Cluster K-Number",
          ylab = "Number of Strains",
          col = c("gray1","gray50", "gray80")
  )
  axis(1,seq(1,28,4),c(1,2,3,4,5,6,7))
  title(line=1.8, xlab="Unique numbered cluster",cex.lab=1)
  
  legend("topleft",
         c("K = 5","K = 6", "K = 7"),
         fill = c("gray1","gray50", "gray80")
  
  )
  #cluster_mean <- mean(little_DF[cluster,])
  #cluster_sd <- sd(little_DF[cluster,])
  #little_DF <- cbind(little_DF, cluster_mean)  
  #little_DF <- cbind(little_DF, cluster_sd)  
  #arrows(little_DF, cluster_mean-cluster_sd,
  #      little_DF, cluster_mean+cluster_sd,angle=90,code=1)
}

###
#'Histogram highlighting strain resistance from PCA' on thesis
#Highly Sensitive	19
#Sensitive	54
#Resistant	71
#Highly Resistant	24  
score_table <- c(17, 16, 16, 16, 16, 16, 16, 16, 16, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 3)
summed_table <- data.frame(cbind( c(3:18), c(1, 5, 4, 9, 9, 11, 15, 19, 13, 19, 22, 17, 15, 8, 1, 0)))
colnames(summed_table) <- c("Resistance Score", "Number of Strains")
as.matrix(t(summed_table[,2]))

colours <- c(rep("DarkRed", 4), rep("Red", 4), rep("Blue", 4), rep("Cyan", 4))
little_DF <- as.matrix(t(summed_table[,c(2)]))
row.names(little_DF) <- c("Resistance Score")
par(mai =  c(0.8, 1.1, 0.2, 0.2))
# par(mai =  c(bottom, left, top, right))
barplot(little_DF, col = colours, ylim = c(0, 25), beside = TRUE, cex.axis = 2, ylab = "Frequency", xlab = "Resistance score", cex.lab =3)
axis(1,seq(1,32, 2), c(3:18), cex.axis =2)





##myu_slope <- TraitData1$Highest.Slope.furfural
C_Value <- (TraitData1$C.Value)
inflection_point <- (TraitData1$Rough.Inflection.Point)/2
myu_point <- (TraitData1$Timepoint.furfural)/2
max_od <- TraitData1$Max.OD.furfural
ratio_myu <- TraitData1$Slope.Ratio.furfural.control
ratio_inflection <- TraitData1$Inflection.Ratio.furfural.control

df_hist <- data.frame(myu_point, C_Value, max_od, inflection_point, ratio_myu, ratio_inflection)
colnames(df_hist) <- c("Timepoint ?? (A)", "C value (B)", "Maximal OD (C)", "Inflection Point (D)", "Control/Furfural ?? Ratio (E)", "Control/Furfural Inflection Ratio (F)")
units <- c("Time (Hrs)", "OD at Time 0", "Log OD value", "Time (Hrs)", "Ratio Value", "Ratio Value")
print("ENSURE MYU WORKS AFTER LOADING SAVE!")

x11()
par(mfrow=c(3,2)) 
par(mar=c(3,4,4,1))
for (plot in 1:length(df_hist[1,])){
  
  #if (plot%%2){
  if (1==1){
  hist(df_hist[,plot], main = colnames(df_hist)[plot], ylim = c(1,110), ylab = "Number of Strains", xlab = "", cex.lab =1.3)
  title(xlab=units[plot], line=2, cex.lab=1.2)
    #print(df_hist[,plot])#
  }
  else{
  hist(df_hist[,plot], main = colnames(df_hist)[plot],ylim = c(1,110), yaxt = "n", ylab = "", xlab = "", cex.lab =1.3)
  title(xlab=units[plot], line=2, cex.lab=1.2)
    }
}

## PCA plots I USED THIS FIRST ONE

library(devtools)
install_github("htmltools")
install.packages("rlang")
install.packages("rgl")
install.packages("magick")
install.packages("htmltools")
library(htmltools)
library(rgl)
library(magick)
TraitData1[,9]
forPCA <- TraitData1[,c(7,8,9)]
colnames(forPCA) <- c("Slope Point", "Max OD", "Inflection")
pc <- princomp(forPCA, cor=TRUE, scores=TRUE, )
#pc <- princomp(TraitData1[,c(5,7,8)], cor=TRUE, scores=TRUE)
summary(pc)
pc$loadings
traitvector <- TraitData1$Resistance.Points..1.18.
traitvector[traitvector == "3"] <- "darkred"
traitvector[traitvector == "4"] <- "darkred"
traitvector[traitvector == "5"] <- "darkred"
traitvector[traitvector == "6"] <- "darkred"
traitvector[traitvector == "7"] <- "red"
traitvector[traitvector == "8"] <- "red"
traitvector[traitvector == "9"] <- "red"
traitvector[traitvector == "10"] <- "red"
traitvector[traitvector == "11"] <- "blue"
traitvector[traitvector == "12"] <- "blue"
traitvector[traitvector == "13"] <- "blue"
traitvector[traitvector == "14"] <- "blue"
traitvector[traitvector == "15"] <- "cyan"
traitvector[traitvector == "16"] <- "cyan"
traitvector[traitvector == "17"] <- "cyan"
traitvector[traitvector == "18"] <- "cyan"
#traitvector[traitvector == "19"] <- "cyan"
#plot3d(pc$scores[,1:3], col=traitvector, xlab = "X", ylab = "Y", zlab = "Z", box = TRUE ,axes=FALSE)
plot3d(pc$scores[,1:3], col=traitvector, xlab = "X", ylab = "Y", zlab = "Z")

#movie3d(spin3d(axis = c(0, 0, 1)), 5, dev = rgl.cur(), dir = tempdir("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\GIF\\"), convert = TRUE)
movie3d(spin3d(axis = c(0, 1, 0)), 12, dev = rgl.cur(), dir = getwd())



#THIS SECOND ONE ISNT AS GOOD. NOT IN THESIS
library(devtools)
install_github("vqv/ggbiplot")
install.packages("ggbiplot")
library(ggbiplot)

write.csv((TraitData1), file = "C:\\Users\\Joseph\\Desktop\\Test_stuff\\TOM_Stuff\\TraitData.csv")
#>11 is resistant on 18-point scale of 3
  plot(TraitData1[,c(6)], TraitData1[,c(9)])
traitvector <- TraitData1$Slope.Ratio.furfural.control.score
traitvector <- TraitData1$Inflection.Ratio.furfural.control.score
traitvector <- TraitData1$Slope.Cluster.Score + TraitData1$Inflection.Cluster.Score + TraitData1$Max.OD.furfural
traitvector[traitvector < 6] <- "Non Resistant"
traitvector[traitvector == 6] <- "Resistant"
traitvector[traitvector == 7] <- "Resistant"
traitvector[traitvector > 8] <- "Resistant"
twochar.pca <- prcomp(TraitData1[,c(7,9)], center = TRUE, scale. = TRUE)

summary(twochar.pca)
ggbiplot(twochar.pca,  var.axes=FALSE, circle = TRUE, ellipse = TRUE,  labels = row.names(TraitData1), groups = traitvector)



all_possible_colnames <- colnames(TraitData1)
Sane_or_Psiko <- "SANE"
for (column_name in all_possible_colnames){
  
Resistances <- unlist(TraitData1[column_name])
print(column_name)
##same, but automated for each metabolite
PValues <- array(0,dim(dfSNPs2)[1])  #
LmerPValues <- array (0, dim(dfSNPs2)[1])
LmerPValues_positivenegative <- array (0, dim(dfSNPs2)[1])
logLmerPValues <- array (0, dim(dfSNPs2)[1])
Positive_Negative_Correlation <- array (0, dim(dfSNPs2)[1])


#do both normal and LMER
for (i in c(1:dim(dfSNPs2)[1])){
  #= logLmerians <- lmer(logg ~ unlist(dfSNPs2[i,])+ (1 | Qmatrix2[,1]) + (1 | Qmatrix2[,2]) + (1 | Qmatrix2[,3]) , REML = FALSE, options(warn = -1))
  if (0 %in% unlist(dfSNPs2[i,]) && 1 %in% unlist(dfSNPs2[i,])){
    Lmerians <- lmer(Resistances ~ unlist(dfSNPs2[i,])+ (1 | Qmatrix2[,1]) + (1 | Qmatrix2[,2]) + (1 | Qmatrix2[,3]) , REML = FALSE, options(warn = -1))
    Linears <- lm(Resistances~ unlist(dfSNPs2[i,]))   ##Remember, resistances MUST == strain numbers! (multiples of 96 for wells.) Then, LR from 0/1 to OD is better!
    if (is.na(summary(Linears)$coefficients[8])){
      PValues[i] <- 1
    }
    else {
      PValues[i] <- summary(Linears)$coefficients[8]
    }
    if (is.na(Anova(Lmerians)$Pr[1])){
      LmerPValues[i]<- 1
    }
    else {
      LmerPValues[i] <- (Anova(Lmerians)$Pr)[1]
    }
    if (is.na((Linears)$coefficients[2])){
      Positive_Negative_Correlation[i] <- "NA"
    }
    else if (as.double((Linears)$coefficients[2]) >= 0){
      Positive_Negative_Correlation[i] <- "+"
    }
    else {
      Positive_Negative_Correlation[i] <- "-"
    }
    #if (is.na(Anova(logLmerians)$Pr[1])){
    #  logLmerPValues[i]<- 1
    #}
    #else {
    #  logLmerPValues[i] <- (Anova(Lmerians)$Pr)[1]
    #}
  }
  else{  # cannot make correlation if all strains are 0 or 1
    PValues[i] <- 1
    LmerPValues[i]<- 1
    Positive_Negative_Correlation[i] <- "NA"
  }
}
FileForPythonThenWeka <- data.frame(GeneNames$V1, LmerPValues, Positive_Negative_Correlation)
colnames(FileForPythonThenWeka) <- c("GeneNames", "Pvalue_of_Gene", "Positive_Or_Negative_Correlation")
write.table(FileForPythonThenWeka, paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\PValue_Files\\lmerPValues_",cigar_stuff, IncludingMasked, Species, S28COnly, SCset,  column_name, "_with_", Sane_or_Psiko, "_Cigar_NoMasked.csv", sep = ''), sep = ",", col.names = FALSE, row.names = FALSE)
}
#1) 
