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
install.packages("gplots")
library("gplots")
library(data.table)
library (qqman)
library (MASS)
library (lme4)
library (car)
library (ggplot2)



install.packages('Ternary')
library('Ternary')
citation("Ternary")

Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joseph_SNPs_PSIKO.Q", sep = " ", header = FALSE, row.names = 1)  # Using PSIKO Q matrix from Jo & Popescu
#Qmatrix <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Joseph_SNPs_SANE.Q", sep = " ", header = FALSE, row.names = 1)  # Using SANE Q Matrix I designed.


TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Malt_concentrations_Updated.csv", header= TRUE, row.names = 1)  # plates 1-9 for all strains. YNB NMR metabolite data. NOT FURFURAL
TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\RStudio_YNB_Metabolites_Updated.csv", header= TRUE, row.names = 1)  # plates 1-9 for all strains. Malt NMR metabolite data. NOT FURFURAL


##data prep, match QMatrix to NMR data
TraitData1 <- TraitData[row.names(TraitData) == "DeleteAll",]
Qmatrix1 <- Qmatrix[row.names(Qmatrix) == "DeleteAll",]

for (j in c(row.names(Qmatrix))){             ##Include strains in new df that are present in dfSNPs
  if (j %in% row.names(TraitData)){
    TraitData1 <- rbind(TraitData1, TraitData[row.names(TraitData) == j,] )
    Qmatrix1 <- rbind(Qmatrix1, Qmatrix[row.names(Qmatrix) == j,] )

  }
}


for (column_data in c(1: length(TraitData1))){
#For selected trait, create percentages out of previous values
current_trait <- TraitData1[, column_data] # specific column
highest_value <- max(current_trait)
current_metabolite <- colnames(TraitData1)[column_data]

for (i in c(1: length(current_trait))){
  current_trait[i] <- (current_trait[i]/highest_value)*100
  #print(current_trait[i])
}

# create plot
library("gplots")
scale_fill_continuous <- colorpanel(1000, "red", "blue") # red is low, blue is high?
#Generate data

mapping <- ceiling(current_trait* 10) #round to index in colour scale
coloured_data <- scale_fill_continuous[mapping] #map to colour sale

TernaryPlot(point = "up", atip = '1', btip = '2', ctip = '3', alab = 'Founder 1', blab = 'Founder 2', clab = 'Founder 3', main = current_metabolite)  # base
TernaryText(list(A = c(10, 1, 1), B = c(1, 10, 1), C = c(1, 1, 10)),col = cbPalette8[4], font = 2)  # green numbers 

#add dots with QMatrix founder percentage as co-ordinates
data_points <- list()  
for(i in 1:nrow(Qmatrix1)) {             # Using for-loop to add columns to list
  data_points <- c(data_points, i = list(c(Qmatrix1[ i, 1], Qmatrix1[ i, 2], Qmatrix1[ i, 3]))) 
}
 
names(data_points) <- row.names(Qmatrix1)

#AddToTernary(points, data_points, pch = 21, cex = 1, bg = vapply(current_trait, function (x) rgb((x[1]) * 2.55, 0, (100-x[1]) * 2.55, 128, maxColorValue = 255), character(1)))
#legend('topright', 
#       legend = c('High', 'Medium', 'Low', 'None'),
#       cex = 0.8, bty = 'n', pch = 21, pt.cex = 1.8,
#       pt.bg = c(rgb(255,   0,   0, 128, NULL, 255), 
#                 rgb(127, 0, 63, 128, NULL, 255),
#                 rgb(63, 0, 127, 128, NULL, 255),
#                 rgb(0, 0, 255, 128, NULL, 255)),)


AddToTernary(points, data_points, pch = 21, cex = 1, bg = coloured_data)

legend('topright', 
      legend = c('High', 'Medium', 'Low', 'None'),
       cex = 0.8, bty = 'n', pch = 21, pt.cex = 1.8,
       pt.bg = c(scale_fill_continuous[1000], 
                 scale_fill_continuous[750],
                 scale_fill_continuous[500],
                 scale_fill_continuous[1]),)



}


#Heatmap of metabolite expression


##### HEATMAP of DE Resistance scores.


TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Malt_concentrations_Updated.csv", header= TRUE, row.names = 1)  # plates 1-9 for all strains. YNB NMR metabolite data. NOT FURFURAL



x11()
par(mfrow=c(1,2)) 
par(mar=c(3,4,4,1))
plot(1~1)
for (i in c(1)){
#Malt
  
TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Malt_concentrations_Updated.csv", header= TRUE, row.names = 1)  # plates 1-9 for all strains. YNB NMR metabolite data. NOT FURFURAL
  
big_metabolites <- data.frame(TraitData$Glucose, TraitData$Ethanol, TraitData$Maltose, TraitData$Glycerol, TraitData$Acetate, TraitData$Methanol, TraitData$Malate)
colnames(big_metabolites) <- c("Glucose", "Ethanol", "Maltose", "Glycerol", "Acetate", "Methanol", "Malate")
big_metabolites <- as.matrix(big_metabolites)
row.names(big_metabolites) <- row.names(TraitData)
heatmap((big_metabolites), Rowv=NA, Colv = NA, margins = c(7,0), labRow = row.names(big_metabolites))

TraitData_WithoutEthanol <-  TraitData
TraitData_WithoutEthanol$Ethanol <- NULL  # remove ethanol.
TraitData_WithoutEthanol$Glucose <- NULL  # remove Glucose
TraitData_WithoutEthanol$Glycerol <- NULL  # remove Glycerol
TraitData_WithoutEthanol$Acetate <- NULL  # remove Acetate
TraitData_WithoutEthanol$Methanol <- NULL  # remove Methanol- 
TraitData_WithoutEthanol$Malate <- NULL  # remove Malate
TraitData_WithoutEthanol$Maltose <- NULL  # remove Maltose
TraitData_WithoutEthanol <- as.matrix(TraitData_WithoutEthanol)
row.names(TraitData_WithoutEthanol) <- row.names(TraitData)

heatmap((TraitData_WithoutEthanol), Rowv=NA, Colv = NA, margins = c(11,0), labRow = row.names(TraitData_WithoutEthanol))
}

#end malt here
plot((TraitData), col="navy", main="Matrix Scatterplot")



x11()
par(mfrow=c(2,1)) 
par(mar=c(3,4,4,1))

TraitData <- read.csv("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\RStudio_YNB_Metabolites_Updated.csv", header= TRUE, row.names = 1)  # plates 1-9 for all strains. Malt NMR metabolite data. NOT FURFURAL

#YNB
big_metabolites <- data.frame(TraitData$Glucose, TraitData$Ethanol, TraitData$Glycerol, TraitData$Acetate, TraitData$Methanol)
colnames(big_metabolites) <- c("Glucose", "Ethanol", "Glycerol", "Acetate", "Methanol")
big_metabolites <- as.matrix(big_metabolites)
row.names(big_metabolites) <- row.names(TraitData)

heatmap(as.matrix(big_metabolites), Rowv=NA, Colv = NA, margins = c(7,0), labRow = row.names(big_metabolites))


TraitData_WithoutEthanol <-  TraitData
TraitData_WithoutEthanol$Ethanol <- NULL  # remove ethanol.
TraitData_WithoutEthanol$Glucose <- NULL  # remove Glucose
TraitData_WithoutEthanol$Glycerol <- NULL  # remove Glycerol
TraitData_WithoutEthanol$Acetate <- NULL  # remove Acetate
TraitData_WithoutEthanol$Methanol <- NULL  # remove Methanol- 
TraitData_WithoutEthanol <- as.matrix(TraitData_WithoutEthanol)
row.names(TraitData_WithoutEthanol) <- row.names(TraitData)
heatmap(as.matrix(TraitData_WithoutEthanol), Rowv=NA, Colv = NA, margins = c(11,0), labRow = row.names(TraitData_WithoutEthanol))

#end YNB here

plot((TraitData), col="navy", main="Matrix Scatterplot")






























##### correlations between metabolite quantities

for (first_column in c(1:length(TraitData))){
  
  if (first_column != length(TraitData)){  # only do loop if it isnt the last column; cant correlate to non-existant next one...
    #print(paste("Column number: ", first_column))
  for (second_column in c ((first_column + 1): length(TraitData))){
    #print(paste("Column 2 number: ", second_column))
    ##correlate
    col_1 <- TraitData[first_column]
    col_2 <- TraitData[second_column]
    #Linears <- lm(col1~ unlist(col1))
    Linears <- lm(as.matrix(col_1)~ unlist(col_2))
    #print(summary(Linears))
    #print(summary(Linears)$coefficients[8])
    if (is.na(summary(Linears)$coefficients[8])){
      # dont want NA
    }
    else if (summary(Linears)$coefficients[8] < 0.1 && summary(Linears)$coefficients[8] > 0){
      print(paste("Correlation between ",  as.character(colnames(col_1)),  " and ", as.character(colnames(col_2)), " is ",  as.character(summary(Linears)$coefficients[8]), sep = ""))
    }
  }}
  
}


plot(TraitData$Glucose ~ TraitData$Glycerol)
