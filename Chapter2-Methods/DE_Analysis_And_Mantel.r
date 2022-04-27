#This script is based on PlateReader.R , and moved here for clarity.
install.packages("reshape2")
install.packages("ggplot2")
install.packages("Hmisc", dependencies=T)
install.packages('inflection')
library('inflection')
library("Hmisc")
library("reshape2")
library("ggplot2")



##### HEATMAP of DE Resistance scores.
file_for_series <- read.csv("GrwothFiles_Folder\\WekaPredictedResistanceDEOnly.csv", header = TRUE)

row.names(file_for_series) <- file_for_series$Well  # set row names
file_for_series <- file_for_series[,-1]             # remove 'well' column with rownames.

length(file_for_series[,1])/96
plate_1_data = c()
plate_2_data = c()
all_plate_data = c()
plate_1_matrix <- data.frame(matrix(, nrow = 8, ncol = 96))
plate_2_matrix <- data.frame(matrix(, nrow = 8, ncol = 96))
all_plates_matrix <- data.frame(matrix(, nrow = 8, ncol = 192))
phenotype_column <- 7  # Max OD
phenotype_column <- 13  # resistance score

for (week in c(1:length(file_for_series[,1]))){

  if ((week == 1 || (week %% 192 == 1)) && week != length(file_for_series[,1])){
    print("********")
    print(paste("Beginning start is ", week))
    print(paste("Beginning end is ", week + 95))
    print(paste("Halfway start is", week + 96))
    print(paste("Halfway end is", week + 191))
    
    
    #print(((week -1)/192) + 1)
    plate_1_data <- append(plate_1_data, file_for_series[week:(week + 95), phenotype_column])
    plate_1_matrix[((week -1)/192) + 1,] <- file_for_series[week:(week + 95), phenotype_column]
    
    plate_2_data <- append(plate_2_data, file_for_series[(week + 96): (week + 191), phenotype_column])
    plate_2_matrix[((week -1)/192) + 1,] <- file_for_series[(week + 96): (week + 191), phenotype_column]
    
    all_plate_data <- append(all_plate_data, file_for_series[week:(week + 191), phenotype_column])
    all_plates_matrix[((week -1)/192) + 1,] <-  file_for_series[week:(week + 191), phenotype_column]
    
  }
}


heatmap(as.matrix(plate_1_matrix), Rowv=NA, Colv = NA)
heatmap(as.matrix(plate_2_matrix), Rowv=NA, Colv = NA)
heatmap(as.matrix(all_plates_matrix), Rowv=NA, Colv = NA)




#MANTEL TEST to compare own custom Q-Matrix to PSIKO one.
install.packages("ape")
install.packages("ade4")
library("ape")
library(ade4)

distance_method <- c("Percentage", "TamD", "JCD", "K2P", "TND")

for (d_method in distance_method){
  #psiko_Q_data <- read.table("All_QMatrices\\PSIKO_QMatrix_Cigar_CoreOnly.Q", sep = " ", header = FALSE, row.names = 1)
  #sane_Q_data <- read.table(paste("All_QMatrices\\SANE_QMatrix_Cigar_CoreOnly_",d_method, ".Q", sep = ""), sep = "\t", header = FALSE, row.names = 1)
  psiko_Q_data <- read.table("All_QMatrices\\PSIKO_QMatrix_Cigar_CoreOnly_NoMasked_SacchOnly.Q", sep = " ", header = FALSE, row.names = 1)
  sane_Q_data <- read.table(paste("All_QMatrices\\SANE_QMatrix_Cigar_CoreOnly_NoMasked_SacchOnly_",d_method, ".Q", sep = ""), sep = "\t", header = FALSE, row.names = 1)
  
  colnames(psiko_Q_data) <- c("Cluster1", "Cluster2", "Cluster3")
  colnames(sane_Q_data) <- c("Cluster1", "Cluster2", "Cluster3")
  psiko_rownames <- row.names(psiko_Q_data)
  new_sane <- array()
  for (each_strain in psiko_rownames){
  new_sane <- rbind(new_sane,sane_Q_data[each_strain,])
  }
  new_sane <- new_sane[-1,]
  
  psiko_dist <- (dist(cbind(psiko_Q_data$Cluster1, psiko_Q_data$Cluster2, psiko_Q_data$Cluster3)))
  custom_dist <- (dist(cbind(new_sane$Cluster1, new_sane$Cluster2, new_sane$Cluster3)))
  
  test_1 <- mantel.rtest(m1 = psiko_dist, m2 = custom_dist, nrepet = 999)
  #test_1$pvalue
  #test_1$expvar[3]
  print(paste("--Variance: ", test_1$expvar[3], "--P-value: ", test_1$pvalue, "Distance method: ", d_method))
  mantel.test(m1 = as.matrix(psiko_dist), m2 = as.matrix(custom_dist), nperm = 999, graph = TRUE, alternative = "two.sided", main = paste (d_method, " Genetic Distance"))
  
  }


#old system

psiko_Q_data <- read.csv("PSIKO_Q_Matrix.csv", row.names = 1, header = FALSE)
custom_Q_data <- read.csv("Custom_Q_Matrix.csv", row.names = 1, header = FALSE)
colnames(psiko_Q_data) <- c("Cluster1", "Cluster2", "Cluster3")
colnames(custom_Q_data) <- c("Cluster1", "Cluster2", "Cluster3")
#psiko_Q_data <- as.matrix(psiko_Q_data)
#custom_Q_data <- as.matrix(custom_Q_data)

psiko_dist <- (dist(cbind(psiko_Q_data$Cluster1, psiko_Q_data$Cluster2, psiko_Q_data$Cluster3)))
custom_dist <- (dist(cbind(custom_Q_data$Cluster1, custom_Q_data$Cluster2, custom_Q_data$Cluster3)))

mantel.rtest(m1 = psiko_dist, m2 = custom_dist, nrepet = 999)
mantel.test(m1 = as.matrix(psiko_dist), m2 = as.matrix(custom_dist), nperm = 999, graph = TRUE, alternative = "two.sided")




###Directed evolution stuff

FileName1 <- paste("DE_Week14_Plate1.csv", sep = '')
FileName2 <- paste("DE_Week14_Plate2.csv", sep = '')
FileName3 <- paste("DE_Week15_Plate1.csv", sep = '')
FileName4 <- paste("DE_Week15_Plate2.csv", sep = '')
FileName5 <- paste("DE_Week16_Plate1.csv", sep = '')
FileName6 <- paste("DE_Week16_Plate2.csv", sep = '')
FileName7 <- paste("DE_Week17_Plate1.csv", sep = '')
FileName8 <- paste("DE_Week17_Plate2.csv", sep = '')


MainFile1 <- read.csv(paste ("Directed_Evolution\\", FileName1, sep = ''))
MainFile2 <- read.csv(paste ("Directed_Evolution\\", FileName2, sep = ''))
MainFile3 <- read.csv(paste ("Directed_Evolution\\", FileName3, sep = ''))
MainFile4 <- read.csv(paste ("Directed_Evolution\\", FileName4, sep = ''))
MainFile5 <- read.csv(paste ("Directed_Evolution\\", FileName5, sep = ''))
MainFile6 <- read.csv(paste ("Directed_Evolution\\", FileName6, sep = ''))
MainFile7 <- read.csv(paste ("Directed_Evolution\\", FileName7, sep = ''))
MainFile8 <- read.csv(paste ("Directed_Evolution\\", FileName8, sep = ''))

# blabla[row,column]
JustData1 <- MainFile1[3:98,4:100]            #Take all the datapoints (every half hour) for all wells.
JustData2 <- MainFile2[3:98,4:100]            #Take all the datapoints (every half hour) for all wells.
JustData3 <- MainFile3[3:98,4:100]            #Take all the datapoints (every half hour) for all wells.
JustData4 <- MainFile4[3:98,4:100]            #Take all the datapoints (every half hour) for all wells.
JustData5 <- MainFile5[3:98,4:100]            #Take all the datapoints (every half hour) for all wells.
JustData6 <- MainFile6[3:98,4:100]            #Take all the datapoints (every half hour) for all wells.
JustData7 <- MainFile7[3:98,4:100]            #Take all the datapoints (every half hour) for all wells.
JustData8 <- MainFile8[3:98,4:100]            #Take all the datapoints (every half hour) for all wells.

####8-colour plot of quadruple-replicate for DE =
#FourPlot( Do you want it as 12*8 matrix? Then TRUE, do you want an output pdf? then also TRUE)
FourPlotDE(TRUE, FALSE)

FourPlotDE <- function(plotted, yesno, platenumber){
  
  if (plotted == TRUE) {
    x11()
    par(mfrow=c(8,12))
    
    par(mar=c(1,1,1,1))
  }
  
  WellNames <- c("A01" ,"A02" ,"A03" ,"A04" ,"A05" ,"A06" ,"A07" ,"A08" ,"A09" ,"A10" ,"A11" ,"A12" ,"B01" ,"B02" ,"B03" ,"B04" ,"B05" ,"B06" ,"B07" ,"B08" ,"B09" ,"B10" ,"B11" ,"B12" ,"C01" ,"C02" ,"C03" ,"C04" ,"C05" ,"C06" ,"C07" ,"C08" ,"C09" ,"C10" ,"C11" ,"C12" ,"D01" ,"D02" ,"D03" ,"D04" ,"D05" ,"D06" ,"D07" ,"D08" ,"D09" ,"D10" ,"D11" ,"D12" ,"E01" ,"E02" ,"E03" ,"E04" ,"E05" ,"E06" ,"E07" ,"E08" ,"E09" ,"E10" ,"E11" ,"E12" ,"F01" ,"F02" ,"F03" ,"F04" ,"F05" ,"F06" ,"F07" ,"F08" ,"F09" ,"F10" ,"F11" ,"F12" ,"G01" ,"G02" ,"G03" ,"G04" ,"G05" ,"G06" ,"G07" ,"G08" ,"G09" ,"G10" ,"G11" ,"G12" ,"H01" ,"H02" ,"H03" ,"H04" ,"H05" ,"H06" ,"H07" ,"H08" ,"H09" ,"H10" ,"H11" ,"H12")
  for (i in 1:96){
    #Tab <- cbind(JustData1[i,], JustData3[i,], JustData5[i,], JustData7[i,])
    Tab <- cbind(JustData2[i,], JustData4[i,], JustData6[i,], JustData8[i,])
    
    Tab <- t(Tab)
    Tab <- data.frame(Tab)
    Tab[,2] <- rep(c(1:97),4)
    Tab[,3] <- c(rep("Yellow", 97),rep("Orange", 97), rep("Pink", 97), rep("Red", 97))
    colnames(Tab) <- c("OpticalDensity", "Cycle", "Group")
    plot(Tab$Cycle, Tab$OpticalDensity, col = Tab$Group, pch = 16, yaxt = 'n',xaxt = 'n',  cex = 1, ylim=c(0,3), xlab = "Cycle(30Min)", ylab = "OpticalDensity", main = paste("Replicate Test Well", WellNames[i], sep = ' '))  
    
    axis(1,at=c(1:nrow(Tab)),labels=rownames(Tab))
  } 
  cat("Ignore warnings. They're for row names. Still need to fix that")
  if (yesno == TRUE){
    print("potate")
    dev.copy2pdf(file = paste("DE_", "Plate", platenumber,".pdf", sep = ''))
  }
}


#Function ends here. DE Stuff
#Finds Scores for all files (or DE files specifically.)


all_files = c(rep("No File", 16) )
for (filename in c(11:18)){
  
  all_files[((filename - 11) * 2) + 1] <- paste("DE_Week", filename,"_Plate1.csv", sep = "")
  all_files[((filename - 11) * 2) + 2] <- paste("DE_Week", filename,"_Plate2.csv", sep = "")
}


number_of_files <- length(all_files)
well_names <- rep("S", 96)
letters <- c("A", "B", "C", "D", "E", "F", "G", "H")
for (i in c(1:8)){
  well_names[((i-1) *12) + 1] <- paste(letters[i], "01", sep = "")
  well_names[((i-1) *12) + 2] <- paste(letters[i],"02",  sep = "")
  well_names[((i-1) *12) + 3] <- paste(letters[i],"03",  sep = "")
  well_names[((i-1) *12) + 4] <- paste(letters[i],"04",  sep = "")
  well_names[((i-1) *12) + 5] <- paste(letters[i],"05",  sep = "")
  well_names[((i-1) *12) + 6] <- paste(letters[i],"06",  sep = "")
  well_names[((i-1) *12) + 7] <- paste(letters[i],"07",  sep = "")
  well_names[((i-1) *12) + 8] <- paste(letters[i],"08",  sep = "")
  well_names[((i-1) *12) + 9] <- paste(letters[i],"09",  sep = "")
  well_names[((i-1) *12) + 10] <- paste(letters[i],"10",  sep = "")
  well_names[((i-1) *12) + 11] <- paste(letters[i],"11",  sep = "")
  well_names[((i-1) *12) + 12] <- paste(letters[i],"12",  sep = "")
  
  }



DE_df <- data.frame(replicate(number_of_files * 10, c(rep(0,96))))
row.names(DE_df) <- well_names
#Highest Slope Control
#Timepoint control
#Max OD Control
# Rough Inflection Control
#Highest slope furfural
# Timepoint Furfural
#MaxOD Furfural
#Rough INflection furfural
#C value
#Strain Number
pre_names <- rep(c("Highest Slope Control", "Timepoint control", "Max OD Control", "Rough Inflection Control","Highest Slope Furfural", "Timepoint Furfural", "MaxOD Furfural", "Rough Inflection Furfural", "C Value", "Strain Number"), number_of_files)

for (name in c(1: (number_of_files/2))){
  
  pre_names[((name-1)*20) + 1] <- paste(pre_names[((name-1)*20) + 1], "_Plate1_", name)
  pre_names[((name-1)*20) + 2] <- paste(pre_names[((name-1)*20) + 2], "_Plate1_", name)
  pre_names[((name-1)*20) + 3] <- paste(pre_names[((name-1)*20) + 3], "_Plate1_", name)
  pre_names[((name-1)*20) + 4] <- paste(pre_names[((name-1)*20) + 4], "_Plate1_", name)
  pre_names[((name-1)*20) + 5] <- paste(pre_names[((name-1)*20) + 5], "_Plate1_", name)
  pre_names[((name-1)*20) + 6] <- paste(pre_names[((name-1)*20) + 6], "_Plate1_", name)
  pre_names[((name-1)*20) + 7] <- paste(pre_names[((name-1)*20) + 7], "_Plate1_", name)
  pre_names[((name-1)*20) + 8] <- paste(pre_names[((name-1)*20) + 8], "_Plate1_", name)
  pre_names[((name-1)*20) + 9] <- paste(pre_names[((name-1)*20) + 9], "_Plate1_", name)
  pre_names[((name-1)*20) + 10] <- paste(pre_names[((name-1)*20) + 10], "_Plate1_", name)
  
  
  pre_names[((name-1)*20) + 11] <- paste(pre_names[((name-1)*20) + 11], "_Plate2_", name)
  pre_names[((name-1)*20) + 12] <- paste(pre_names[((name-1)*20) + 12], "_Plate2_", name)
  pre_names[((name-1)*20) + 13] <- paste(pre_names[((name-1)*20) + 13], "_Plate2_", name)
  pre_names[((name-1)*20) + 14] <- paste(pre_names[((name-1)*20) + 14], "_Plate2_", name)
  pre_names[((name-1)*20) + 15] <- paste(pre_names[((name-1)*20) + 15], "_Plate2_", name)
  pre_names[((name-1)*20) + 16] <- paste(pre_names[((name-1)*20) + 16], "_Plate2_", name)
  pre_names[((name-1)*20) + 17] <- paste(pre_names[((name-1)*20) + 17], "_Plate2_", name)
  pre_names[((name-1)*20) + 18] <- paste(pre_names[((name-1)*20) + 18], "_Plate2_", name)
  pre_names[((name-1)*20) + 19] <- paste(pre_names[((name-1)*20) + 19], "_Plate2_", name)
  pre_names[((name-1)*20) + 20] <- paste(pre_names[((name-1)*20) + 20], "_Plate2_", name)
  
  
  
}
colnames(DE_df) <- pre_names

#
#DE_df <- data.frame(replicate(number_of_files * 4, c(rep(0,96))))
#row.names(DE_df) <- well_names
#pre_names <- rep(c("MaxOD", "Inflection Point", "Myu Point", "Resistance Score"), number_of_files)#
#
#for (name in c(1: (number_of_files/2))){
#  
#  pre_names[((name-1)*8) + 1] <- paste(pre_names[((name-1)*8) + 1], "_Plate1_", name)
#  pre_names[((name-1)*8) + 2] <- paste(pre_names[((name-1)*8) + 2], "_Plate1_", name)
#  pre_names[((name-1)*8) + 3] <- paste(pre_names[((name-1)*8) + 3], "_Plate1_", name)
#  pre_names[((name-1)*8) + 4] <- paste(pre_names[((name-1)*8) + 4], "_Plate1_", name)
  
#  pre_names[((name-1)*8) + 5] <- paste(pre_names[((name-1)*8) + 5], "_Plate2_", name)
#  pre_names[((name-1)*8) + 6] <- paste(pre_names[((name-1)*8) + 6], "_Plate2_", name)
#  pre_names[((name-1)*8) + 7] <- paste(pre_names[((name-1)*8) + 7], "_Plate2_", name)
#  pre_names[((name-1)*8) + 8] <- paste(pre_names[((name-1)*8) + 8], "_Plate2_", name)
#
#}
#colnames(DE_df) <- pre_names


for (file_number in c(1: length(all_files))){
  file_name <- all_files[file_number]
  FileName <- paste(file_name, sep = '')
  #print(FileName)
  MainFile <- read.csv(paste("Directed_Evolution\\", FileName, sep = ''))
  JustData1 <- MainFile[3:98,4:100]            #Take all the datapoints (every half hour) for all wells.
  
  std <- function(x) sd(x)/sqrt(length(x))
  SteepestSlope1 <- c(1:96)
  Match1 <- c(1:96)
  MaxOD1 <- c(1:96)
  Inlfection_Point_Angle <- rep(NA, 96)
  Inlfection_Point_Angle_control <- rep(NA, 96)
  C_scores <- rep(NA, 96)
  Slope_Control <- rep(0, 96)
  timepoint_control <- rep(0,96)
  maxod_control <- rep(0, 96)
  inflection_control <- rep(0,96)
  highest_slope_furfural <- rep(0,96)
  c_value <- rep(0, 96)
  strains <- well_names
  #Strain Number
  
  WellNames <- c("A01" ,"A02" ,"A03" ,"A04" ,"A05" ,"A06" ,"A07" ,"A08" ,"A09" ,"A10" ,"A11" ,"A12" ,"B01" ,"B02" ,"B03" ,"B04" ,"B05" ,"B06" ,"B07" ,"B08" ,"B09" ,"B10" ,"B11" ,"B12" ,"C01" ,"C02" ,"C03" ,"C04" ,"C05" ,"C06" ,"C07" ,"C08" ,"C09" ,"C10" ,"C11" ,"C12" ,"D01" ,"D02" ,"D03" ,"D04" ,"D05" ,"D06" ,"D07" ,"D08" ,"D09" ,"D10" ,"D11" ,"D12" ,"E01" ,"E02" ,"E03" ,"E04" ,"E05" ,"E06" ,"E07" ,"E08" ,"E09" ,"E10" ,"E11" ,"E12" ,"F01" ,"F02" ,"F03" ,"F04" ,"F05" ,"F06" ,"F07" ,"F08" ,"F09" ,"F10" ,"F11" ,"F12" ,"G01" ,"G02" ,"G03" ,"G04" ,"G05" ,"G06" ,"G07" ,"G08" ,"G09" ,"G10" ,"G11" ,"G12" ,"H01" ,"H02" ,"H03" ,"H04" ,"H05" ,"H06" ,"H07" ,"H08" ,"H09" ,"H10" ,"H11" ,"H12")
  for (i in 1:96){
    Tab <- log10(cbind(JustData1[i,]))  # log values
    #Tab <- cbind(JustData1[i,])        # non-logged.
    
    Tab <- t(Tab)
    Tab <- data.frame(Tab)
    Tab[,2] <- rep(c(1:97),1)
    Tab[,3] <- c(rep("FireBrick1", 97))
    colnames(Tab) <- c("OpticalDensity", "Cycle", "Group")
    
    ##Extra bit for slope calc.
    Tab6 <- rbind(log10(JustData1[i,]))  # log of data
    #Tab6 <- rbind(JustData1[i,])        # non-log of data
    
    Tab6 <- t(Tab6)
    ##this loop makes average columns (7&8)

    
    
    ##This does 3-point average (46+4 =49)
    Slopes1 <- rep(0, 96)
    Slopes2 <- rep(0, 96)
    Slopes3 <- rep(0, 96)
    c_values <- rep(0, 96)
    for (p in 6:93){
      
      #print("*****")
      y_values <- c(Tab6[p-4,1], Tab6[p-3,1], Tab6[p-2,1], Tab6[p-1,1], Tab6[p,1], Tab6[p+1,1], Tab6[p+2,1], Tab6[p+3,1], Tab6[p+4,1])
      differences <- sort(c(abs(y_values[1] -y_values[2]), abs(y_values[2] -y_values[3]), abs(y_values[3] -y_values[4]), abs(y_values[4] -y_values[5]), abs(y_values[5] -y_values[6]), abs(y_values[6] -y_values[7]), abs(y_values[7] -y_values[8]), abs(y_values[8] -y_values[9])))
      if  (round(differences[8], 9) == (round(differences[7],9))){
        differences[7] <- differences[7] * 2
      }
      if (differences[8] > (5*differences[7]) || differences[7] > (5*differences[6])){
        
        Slopes1[p] <- 0
        Slopes2[p] <- 0
        Slopes3[p] <- 0
        c_values[p] <- 0
      }
      else{
        summed <- sum(y_values)
        average_y <- summed/9
        x_values <- c(p-4, p-3, p-2, p-1, p, p+1, p+2, p+3, p+4)
        average_x <- p
        x_change <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
        y_change <- c(y_values[1] - average_y, y_values[2] - average_y, y_values[3] - average_y, y_values[4] - average_y, y_values[5] - average_y, y_values[6] - average_y, y_values[7] - average_y, y_values[8] - average_y, y_values[9] - average_y)
        multiplied_xy <- c(y_change[1] * x_change[1], y_change[2] * x_change[2], y_change[3] * x_change[3], y_change[4] * x_change[4], y_change[5] * x_change[5], y_change[6] * x_change[6], y_change[7] * x_change[7], y_change[8] * x_change[8], y_change[9] * x_change[9] )
        summed_XY <- sum(multiplied_xy)
        multiplied_xx <- c(16, 9, 4, 1, 0, 1, 4, 9, 16)
        summed_xx <- sum(multiplied_xx)
        Slopes2[p] <- (summed_XY/summed_xx) * 2  # multiplied by two as this is per half-hour, we want per-hour
        Slopes3[p] <- (summed_XY/summed_xx) * 2
        c_values[p] <- average_y - (average_x * Slopes1[p])
      }

    }
    
    #plot(1:48, Slopes3, main = "Slopes")
    start_OD_furfural <- (Tab6[2,1] + Tab6[3,1] + Tab6[4,1])/3 # first point is often bad; ignore it
    
    angles <- rep(0, 88)   
    for (slope in 4:88){  # 88 + 6 safety +next 2 points check == 48
      if (Tab6[slope+1, 1] > Tab6[slope,1] && Tab6[slope+2,1]> Tab6[slope+1,1] && Tab6[slope-2, 1] < Tab6[slope+2,1]){  # old version (check next two points for increase)
        
        
        angles[slope] <- abs(atan(abs((Slopes3[slope-2] - Slopes3[slope + 2])/(1 + (Slopes3[slope-2] * Slopes3[slope + 2])))))
        #angles[slope] <- abs(atan(abs((((Slopes3[slope - 1] + Slopes3[slope - 2])/2) - Slopes3[slope +1])/(1 + (Slopes3[slope+1] * ((Slopes3[slope - 1] + Slopes3[slope - 2])/2))))))
        
        
      }
    }
    angles_a <- angles
    for (pop in c(4: 87)){
      angles_a[pop] <- (angles[pop] + angles[pop+1])/2
    }
    
    largest_angle = sort(angles_a)[88]
    Inlfection_Point_Angle[i] <- match(largest_angle, angles_a)
    #plot(1:46, angles, main = "Gradient swap")
    ## Original!!
    #angles1 <- sort(angles)
    #Inlfection_Point_Angle[i] <- match(c(angles1[46]),angles)
    #if (angles1[46] == 0){
    #  Inlfection_Point_Angle[i] <- 48
    #}
    
    ## TOMS:
    
    
    point_of_highest_slope <- match(c(sort(Slopes2)[96]),Slopes2)
    per_half_hour <- (sort(Slopes2)[96])/2  # dont convert
    #per_half_hour <- (sort(Slopes2)[48])/2  # convert back to half-hour
    Inlfection_Point_Angle[i] <- point_of_highest_slope - abs(round((Tab6[point_of_highest_slope,1] - start_OD_furfural)/(per_half_hour)))
    if (Inlfection_Point_Angle[i] < 1){
      Inlfection_Point_Angle[i] <- 1
    }
    
    point_of_highest_slope <- match(c(sort(Slopes1)[96]),Slopes1)
    per_half_hour <- (sort(Slopes1)[96])/2  # dont convert


    C_scores[i] <- c_values[point_of_highest_slope]
    
    
    MaxOD1[i] <- Tab6[97,1]
    Slopes4 <- sort(Slopes2)
    
    Match1[i] <- match(c(Slopes4[96]),Slopes2) ##this is the literal point of steepest slope (middle line)
    

    SteepestSlope1[i] <- Slopes4[96]
    
    
    
  DE_df[,((file_number-1) * 10) + 1] <- Slope_Control
  DE_df[,((file_number-1) * 10) + 2] <- timepoint_control
  DE_df[,((file_number-1) * 10) + 3] <- maxod_control
  
  DE_df[,((file_number-1) * 10) + 4] <- inflection_control
  DE_df[,((file_number-1) * 10) + 5] <- highest_slope_furfural
  DE_df[,((file_number-1) * 10) + 6] <- Match1
  DE_df[,((file_number-1) * 10) + 7] <- MaxOD1
  DE_df[,((file_number-1) * 10) + 8] <- Inlfection_Point_Angle
  DE_df[,((file_number-1) * 10) + 9] <- C_scores
  DE_df[,((file_number-1) * 10) + 10] <- well_names
  
  
}


for_csv <- data.frame(rbind(rep(0,10)))
colnames(for_csv) <- c("Highest Slope Control", "Timepoint control", "Max OD Control", "Rough Inflection Control","Highest Slope Furfural", "Timepoint Furfural", "MaxOD Furfural", "Rough Inflection Furfural", "C Value", "Well Number")
for (bloc in c(0:(length(DE_df) - 1))){
  if (!bloc %% 10){
    #print("****")
    #print(bloc)
    week_number <- round((((bloc/20) + 1) - 0.1),0)
    version <- ((bloc/20) + 1)
    if (week_number< version){
      version <- 2
    }
    else{
      version <- 1
    }
    #print(week_number)
    #print(version)
    row_name_init <- paste("Week", week_number, "_Rep", version, "_", sep = "")
    rows_names <- as.vector(rep(row_name_init, 96))
    for (each_well in c(1: length(well_names))){
    rows_names[each_well] <- paste(rows_names[each_well], well_names[each_well], sep = "")
    }
    mini_addition <- DE_df[, (bloc + 1): (bloc + 10)]
    colnames(for_csv) <- colnames(mini_addition)
    row.names(mini_addition) <- rows_names
    #print(colnames(DE_df[, (bloc + 1): (bloc + 10)]))
    #for_csv <- rbind(DE_df[, (bloc + 1): (bloc + 10)], setNames(rev(for_csv), names(for_csv)))
    #for_csv <- rbind(for_csv, colnames(DE_df[, (bloc + 1): (bloc + 10)]), DE_df[, (bloc + 1): (bloc + 10)]) 
    
    for_csv <- rbind(for_csv, mini_addition) 
    #print(DE_df[, (bloc + 1): (bloc + 10)])
  }
}
colnames(for_csv) <- c("Highest Slope Control", "Timepoint control", "Max OD Control", "Rough Inflection Control","Highest Slope Furfural", "Timepoint Furfural", "MaxOD Furfural", "Rough Inflection Furfural", "C Value", "Well Number")

for_csv <- for_csv[-1,]
write.table(for_csv, paste("Directed_Evolution\\", "DE_Data.csv", sep = ''), sep = ",", col.names = TRUE, row.names = TRUE)


##########END HERE

for (well in c(1:96)){
  maxod_vector <- rep(0, length(DE_df[1,])/4)
  inflection_vector <- rep(0, length(DE_df[1,])/4)
  myu_point_vector <- rep(0, length(DE_df[1,])/4)
  for (trait in c(1: length(DE_df[1,]))){
    if ((trait-1)%%4){
      #maxod
      print(trait)
      maxod_vector[trait/4] <- DE_df[well, trait]
    }
    if ((trait-2)%%4){
      #inflection
      inflection_vector[trait/4] <- DE_df[well, trait]
    }
    if ((trait-3)%%4){
      #myu point
      myu_point_vector[trait/4] <- DE_df[well, trait]
    }
    
    
  }
}

trait <- 2
if (trait%%4){
  #maxod
  print(trait)
}
0%%4











#stop here if  just making the DE Resistance Scores for each strain
#########################################################################################################


#test

Tab6 <- data.frame(matrix(NA, nrow = 48, ncol = 97))
for (i in 1:95){
  if (!i %% 12){ #row, column
    val <- (i/2)
    Tab6[val + 1, ] <- JustData1[i+1,] + JustData1[i+7,]
    Tab6[val + 2, ] <- JustData1[i+2,] + JustData1[i+8,]
    Tab6[val + 3, ] <- JustData1[i+3,] + JustData1[i+9,]
    Tab6[val + 4, ] <- JustData1[i+4,] + JustData1[i+10,]
    Tab6[val + 5, ] <- JustData1[i+5,] + JustData1[i+11,]
    Tab6[val + 6, ] <- JustData1[i+6,] + JustData1[i+12,]
  }
  else if(i==1){
    Tab6[1, ] <- JustData1[1,] + JustData1[7,]
    Tab6[2, ] <- JustData1[2,] + JustData1[8,]
    Tab6[3, ] <- JustData1[3,] + JustData1[9,]
    Tab6[4, ] <- JustData1[4,] + JustData1[10,]
    Tab6[5, ] <- JustData1[5,] + JustData1[11,]
    Tab6[6, ] <- JustData1[6,] + JustData1[12,]
  }
}
#test

# DE sixliner
sixlinerDE(TRUE,FALSE, TRUE,FALSE, JustData1, "13", "1")


sixlinerDE <- function(LogOfValues, matrix, pdfout, printall,justdata, replicate, platenumber){  # two replicates per plate. 
  std <- function(x) sd(x)/sqrt(length(x))
  if (matrix == TRUE){
    x11()
    par(mfrow=c(8,12))
    par(mar=c(1,1,1,1))
    
  }
  SteepestSlope1 <- c(1:48)
  SteepestSlope2 <- c(1:48)
  Match1 <- c(1:48)
  Match2 <- c(1:48)
  MaxOD1 <- c(1:48)
  MaxOD2 <- c(1:48)
  Inlfection_Point_Angle <- rep(NA, 48)
  Inlfection_Point_Angle_control <- rep(NA, 48)
  C_scores <- rep(NA, 48)
  
  prep_table <- data.frame(matrix(NA, nrow = 48, ncol = 97))
  for (i in 1:95){
    if (!i %% 12){ #row, column
      val <- (i/2)
      prep_table[val + 1, ] <- (justdata[i+1,] + justdata[i+7,])
      prep_table[val + 2, ] <- (justdata[i+2,] + justdata[i+8,])
      prep_table[val + 3, ] <- (justdata[i+3,] + justdata[i+9,])
      prep_table[val + 4, ] <- (justdata[i+4,] + justdata[i+10,])
      prep_table[val + 5, ] <- (justdata[i+5,] + justdata[i+11,])
      prep_table[val + 6, ] <- (justdata[i+6,] + justdata[i+12,])
    }
    else if(i==1){
      prep_table[1, ] <- (justdata[1,] + justdata[7,])
      prep_table[2, ] <- (justdata[2,] + justdata[8,])
      prep_table[3, ] <- (justdata[3,] + justdata[9,])
      prep_table[4, ] <- (justdata[4,] + justdata[10,])
      prep_table[5, ] <- (justdata[5,] + justdata[11,])
      prep_table[6, ] <- (justdata[6,] + justdata[12,])
    }
  }
  
  WellNames <- c("A01" ,"A02" ,"A03" ,"A04" ,"A05" ,"A06","B01" ,"B02" ,"B03" ,"B04" ,"B05" ,"B06" ,"C01" ,"C02" ,"C03" ,"C04" ,"C05" ,"C06","D01" ,"D02" ,"D03" ,"D04" ,"D05" ,"D06" ,"E01" ,"E02" ,"E03" ,"E04" ,"E05" ,"E06" ,"F01" ,"F02" ,"F03" ,"F04" ,"F05" ,"F06" ,"G01" ,"G02" ,"G03" ,"G04" ,"G05" ,"G06","H01" ,"H02" ,"H03" ,"H04" ,"H05" ,"H06" )
  for (i in 1:48){
    if (LogOfValues == TRUE){
      Tab <- log10(cbind(prep_table[i,]))
      
    }
    else{
      Tab <- cbind(prep_table[i,])
      
    }
    Tab <- t(Tab)
    Tab <- data.frame(Tab)
    Tab[,2] <- c(1:97)
    Tab[,3] <- c(rep("FireBrick1", 97))
    colnames(Tab) <- c("OpticalDensity", "Cycle", "Group")
    
    ##Extra bit for slope calc.
    
    if (LogOfValues == TRUE){
      Tab6 <- rbind(log10(prep_table[i,]))
      Tab6 <- t(Tab6)
      lo <- loess(exp(Tab6[,1])~ Tab6[,1])
      xl <- seq(min(Tab6[,1]),max(Tab6[,1]), (max(Tab6[,1]) - min(Tab6[,1]))/1000)
      out = predict(lo,xl)
      infl <- c(FALSE, diff(diff(out)>0)!=0)
      
    }
    else{
      Tab6 <- rbind(prep_table[i,])
      Tab6 <- t(Tab6)
      lo <- loess(Tab6[,1]~log10(Tab6[,1]))
      xl <- seq(min(Tab6[,1]),max(Tab6[,1]), (max(Tab6[,1]) - min(Tab6[,1]))/1000)
      out = predict(lo,xl)
      infl <- c(FALSE, diff(diff(out)>0)!=0)
    } 
    
    
    ##This does 3-point average (46+4 =49)
    Slopes1 <- rep(0, 96)
    Slopes2 <- rep(0, 96)
    Slopes3 <- rep(0, 96)
    c_values <- rep(0, 96)
    for (p in 6:93){
      ##THIS BIT IS NON-FURFURAL
      #print(Tab6[,1])
      y_values <- c(Tab6[p-4,1], Tab6[p-3,1], Tab6[p-2,1], Tab6[p-1,1], Tab6[p,1], Tab6[p+1,1], Tab6[p+2,1], Tab6[p+3,1], Tab6[p+4,1])
      differences <- sort(c(abs(y_values[1] -y_values[2]), abs(y_values[2] -y_values[3]), abs(y_values[3] -y_values[4]), abs(y_values[4] -y_values[5]), abs(y_values[5] -y_values[6]), abs(y_values[6] -y_values[7]), abs(y_values[7] -y_values[8]), abs(y_values[8] -y_values[9])))
      #print(y_values)
      #print(differences)
      #print(round(differences[8], 9))
      #print(round(differences[7], 9))
      if  (round(differences[8], 9) == (round(differences[7],9))){
        differences[7] <- differences[7] * 2
      }
      if (differences[8] > (5*differences[7]) || differences[7] > (5*differences[6])){
        
        Slopes1[p] <- 0
        Slopes2[p] <- 0
        Slopes3[p] <- 0
        c_values[p] <- 0
      }
      else{
        summed <- sum(y_values)
        average_y <- summed/9
        x_values <- c(p-4, p-3, p-2, p-1, p, p+1, p+2, p+3, p+4)
        average_x <- p
        x_change <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
        y_change <- c(y_values[1] - average_y, y_values[2] - average_y, y_values[3] - average_y, y_values[4] - average_y, y_values[5] - average_y, y_values[6] - average_y, y_values[7] - average_y, y_values[8] - average_y, y_values[9] - average_y)
        multiplied_xy <- c(y_change[1] * x_change[1], y_change[2] * x_change[2], y_change[3] * x_change[3], y_change[4] * x_change[4], y_change[5] * x_change[5], y_change[6] * x_change[6], y_change[7] * x_change[7], y_change[8] * x_change[8], y_change[9] * x_change[9] )
        summed_XY <- sum(multiplied_xy)
        multiplied_xx <- c(16, 9, 4, 1, 0, 1, 4, 9, 16)
        summed_xx <- sum(multiplied_xx)
        Slopes1[p] <- (summed_XY/summed_xx) * 2  # multiplied by two as this is per half-hour, we want per-hour
        
        #c_values[p] <- average_y - (average_x * Slopes1[p])
      }
      
      
      #print("*****")
      y_values <- c(Tab6[p-4,1], Tab6[p-3,1], Tab6[p-2,1], Tab6[p-1,1], Tab6[p,1], Tab6[p+1,1], Tab6[p+2,1], Tab6[p+3,1], Tab6[p+4,1])
      differences <- sort(c(abs(y_values[1] -y_values[2]), abs(y_values[2] -y_values[3]), abs(y_values[3] -y_values[4]), abs(y_values[4] -y_values[5]), abs(y_values[5] -y_values[6]), abs(y_values[6] -y_values[7]), abs(y_values[7] -y_values[8]), abs(y_values[8] -y_values[9])))
      if  (round(differences[8], 9) == (round(differences[7],9))){
        differences[7] <- differences[7] * 2
      }
      if (differences[8] > (5*differences[7]) || differences[7] > (5*differences[6])){
        
        Slopes1[p] <- 0
        Slopes2[p] <- 0
        Slopes3[p] <- 0
        c_values[p] <- 0
      }
      else{
        summed <- sum(y_values)
        average_y <- summed/9
        x_values <- c(p-4, p-3, p-2, p-1, p, p+1, p+2, p+3, p+4)
        average_x <- p
        x_change <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
        y_change <- c(y_values[1] - average_y, y_values[2] - average_y, y_values[3] - average_y, y_values[4] - average_y, y_values[5] - average_y, y_values[6] - average_y, y_values[7] - average_y, y_values[8] - average_y, y_values[9] - average_y)
        multiplied_xy <- c(y_change[1] * x_change[1], y_change[2] * x_change[2], y_change[3] * x_change[3], y_change[4] * x_change[4], y_change[5] * x_change[5], y_change[6] * x_change[6], y_change[7] * x_change[7], y_change[8] * x_change[8], y_change[9] * x_change[9] )
        summed_XY <- sum(multiplied_xy)
        multiplied_xx <- c(16, 9, 4, 1, 0, 1, 4, 9, 16)
        summed_xx <- sum(multiplied_xx)
        Slopes2[p] <- (summed_XY/summed_xx) * 2  # multiplied by two as this is per half-hour, we want per-hour
        Slopes3[p] <- (summed_XY/summed_xx) * 2
        c_values[p] <- average_y - (average_x * Slopes1[p])
      }
      #print("slope")
      
    }
    
    #plot(1:48, Slopes3, main = "Slopes")
    start_OD_control <- (Tab6[2,1] + Tab6[3,1] + Tab6[4,1])/3 # first point is often bad; ignore it
    
    angles <- rep(0, 88)   
    for (slope in 4:88){  # 46+next 2 points check == 48
      if (Tab6[slope+1, 1] > Tab6[slope,1] && Tab6[slope+2,1]> Tab6[slope+1,1] && Tab6[slope-2, 1] < Tab6[slope+2,1]){  # old version (check next two points for increase)
        
        
        angles[slope] <- abs(atan(abs((Slopes3[slope-2] - Slopes3[slope + 2])/(1 + (Slopes3[slope-2] * Slopes3[slope + 2])))))
        #angles[slope] <- abs(atan(abs((((Slopes3[slope - 1] + Slopes3[slope - 2])/2) - Slopes3[slope +1])/(1 + (Slopes3[slope+1] * ((Slopes3[slope - 1] + Slopes3[slope - 2])/2))))))
        
        
      }
    }
    angles_a <- angles
    for (pop in c(4: 87)){
      angles_a[pop] <- (angles[pop] + angles[pop+1])/2
    }
    
    largest_angle = sort(angles_a)[88]
    Inlfection_Point_Angle[i] <- match(largest_angle, angles_a)
    
    ## TOMS:
    
    
    
    point_of_highest_slope <- match(c(sort(Slopes1)[96]),Slopes1)
    per_half_hour <- (sort(Slopes1)[96])/2  # dont convert
    #per_half_hour <- (sort(Slopes2)[48])/2  # convert back to half-hour
    Inlfection_Point_Angle_control[i] <- point_of_highest_slope - abs(round((Tab6[point_of_highest_slope,1] - start_OD_control)/(per_half_hour)))
    if (Inlfection_Point_Angle_control[i] < 1){
      Inlfection_Point_Angle_control[i] <- 1
    }
    
    C_scores[i] <- c_values[point_of_highest_slope]
    
    
    MaxOD1[i] <- Tab6[96,1]
    Slopes3 <- sort(Slopes1)
    Slopes4 <- sort(Slopes2)
    Match1[i] <- match(c(Slopes3[96]),Slopes1)
    
    #print(Slopes3)
    SteepestSlope1[i] <- Slopes3[96]
    SteepestSlope2[i] <- Slopes4[96]
    
    #####
    if (LogOfValues == TRUE && printall == TRUE){
      plot(Tab$Cycle, Tab$OpticalDensity, col = Tab$Group ,pch = 16, yaxt = 'n',xaxt = 'n',  cex = 1, ylim=c(-1,1) , ylab = "OpticalDensity", main = paste("Replicate Test Well", WellNames[i], sep = ' '))  
      legend("topleft", legend=c(signif(Slopes3[96], digits = 3), signif(Slopes4[96], digits = 3)), col=c("red", "blue"), cex = 1, x.intersp = 0.1,y.intersp = 0.7)
      
      ##axis(1,at=c(1:nrow(Tab)),labels=rownames(Tab))
      axis(2, 0:3)
      points(xl[infl ], out[infl ], col="black")
      lines(xl, out, col='red', lwd=2)
    }
    else{
      plot(Tab$Cycle, Tab$OpticalDensity, col = Tab$Group ,pch = 16, yaxt = 'n',xaxt = 'n',  cex = 1, ylim=c(-3,3) , ylab = "OpticalDensity", main = paste("Replicate Test Well", WellNames[i], sep = ' '))  
      legend("topleft", legend=c(signif(Slopes3[96], digits = 3), signif(Slopes4[96], digits = 3)), col=c("red", "blue"), cex = 1, x.intersp = 0.1,y.intersp = 0.7)
      
      ##axis(1,at=c(1:nrow(Tab)),labels=rownames(Tab))
      axis(2, 0:3)
      points(xl[infl ], out[infl ], col="black")
      lines(xl, out, col='red', lwd=2)  
    }
    
  } 
  #MainFrame <- data.frame(SteepestSlope1, Match1,MaxOD1, SteepestSlope2, Match2, MaxOD2, InfPoints)
  MainFrame <- data.frame(SteepestSlope1, Match1,MaxOD1,Inlfection_Point_Angle_control,  C_scores)
  colnames(MainFrame) <- c("Highest Slope", "Timepoint", "Max OD control", "Rough Inflection Control", "C_Value")
  rownames(MainFrame) <- WellNames
  if (LogOfValues == TRUE){
    write.csv(MainFrame, file = paste("DE_Plate", platenumber, "_Replicate", replicate, ".csv", sep = ''))
  }
  else{
    write.csv(MainFrame, file = paste("DE_Plate", platenumber, "_Replicate", replicate, ".csv", sep = ''))
  }
  sum(SteepestSlope1)/48
  sum(SteepestSlope2)/48
  cat("Ignore warnings. They're for row names. Still need to fix that")
  if (pdfout == TRUE){
    dev.copy2pdf(file = paste("C:\\Users\\Joe\\Desktop\\PhD_stuff\\DE_Plate", platenumber, "_Replicate", replicate, ".pdf", sep=''))
  }
}
sixlinerDE(TRUE,FALSE, FALSE,FALSE, JustData8, "17", "2")
#function ends here



###
###DE end