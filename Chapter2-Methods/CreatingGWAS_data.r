##This code creates datasets for manhatten plots (etc)
##Run CreatingGWAS_data.r top bit, then SNP_and_MAF_finding.py, then CreatingGWAS_data.r second bit, then MatrixMaker.py


#source("https://bioconductor.org/biocLite.R")
install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "snpStats"))
#biocLite("VariantAnnotation")
#biocLite("snpStats")
install.packages("rlang")
library("snpStats")
library(VariantAnnotation)
library(GenomicFeatures)
#Playe1
Plate1 <- c(1,2,4,6,8,9,10,16,17,18,20,21,22,23,26,31,36,39,40,43,44,45,46,49,51,52,54,55,56,57,58,59,60,61,62,63,64,65,68,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,17,388,416,505,523,568,731,768,777,1006,1026,1187,1228,1245,1681,2433,2489,2572,2577,2578,2600,2629,2701,2791,2804,2888,2889,2890,2904)
#Plate2
Plate2 <- c(166,371,372,373,745,746,747,783,794,894, 2321, 2322, 2395, 2396, 2480, 2486, 2491, 2492, 2521, 2529, 2580, 2753, 3047, 3096, 3129, 3398, 3400, 3719, 3792, 3870, 140, 141, 147, 161, 408, 411, 492, 524, 559, 566, 582, 585, 608, 677, 678, 696, 820, 2473, 2568 ,2741, 2885, 3024, 3141, 3239, 3255, 3506, 3877, 62, 135, 138, 142, 154, 155, 158, 159, 162, 195, 377, 502, 539, 541, 758, 759, 796, 797, 844, 845, 930, 931, 974, 1401, 1645, 1646, 1647, 1648, 1649, 1650, 1651, 1659, 1660, 2439, 2440, 2581, 2599, 2605, 2666)
#Plate3
Plate3 <- c(2752, 2864, 2873, 2972, 3056, 3057, 3072, 3120, 3401, 3411, 3444, 3504, 3536, 3721, 3722, 3725, 3735, 3772, 3775, 3816, 3817, 3820, 3821, 3832, 3833, 3834, 3835, 3836, 3837, 3838, 3867, 3872, 2826, 3284, 3290, 3312, 3452, 3468, 99, 104, 107, 108, 109, 110, 113, 118, 121, 122, 124, 125, 126, 167, 176, 177, 181, 182, 183, 185, 186, 187, 190, 192, 196, 197, 198, 199, 200, 201, 202, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231)
#Plate4
Plate4 <- c(392, 971, 2449, 2560, 2729, 2703, 2991, 2827, 1417, 2483, 2702, 814, 3853, 2878, 2450, 3108, 2754, 2403, 1416, 1495, 2999, 2789, 2513, 2875, 2508, 2644, 543, 546, 538, 3788, 2739, 2976, 17, 1063, 1064, 1151, 1337, 2904, 1449, 2908, 3612, 3630, 3662, 826, 232, 754, 667, 975, 956, 695, 739, 235, 3264, 3265, 3266, 3311, 3313, 3314, 3315, 3318, 3319, 3445, 3447, 3448, 3449, 3451, 3453, 3454, 3455, 3456, 3457, 3458, 3460, 3461, 3462, 3466, 3467, 3469, 3470, 3471, 3472, 3486, 3487, 360, 431, 2592, 2945, 361, 4000, 1444, 1603,1606, 2397, 2733, 2737, 3406)
#Plate5
Plate5 <- c(356, 357, 358, 430, 463, 478, 479, 482, 619, 620, 621, 671, 672, 684, 816, 1406, 1407, 1408, 1409, 1410, 1411, 1412, 1413, 1414, 1415, 1431, 2401, 2402, 2517, 2587, 2683, 2688, 2695, 2804, 2808, 2809, 2855, 2947, 2948, 489, 490, 491, 525, 694, 995, 996, 1529, 1530, 1765, 3020, 3021, 3022, 238, 241, 341, 911, 1384, 1510, 3267, 3431, 100, 111, 143, 151, 152, 179, 188, 243, 244, 350, 426, 587, 744, 827, 851, 906, 970, 1424, 1425, 1426, 1429, 1441, 2265, 2559, 2597, 2675, 2886, 2887, 2907, 3104, 3303, 3344, 3396, 3502, 3519, 3537)
#Plate6
Plate6 <- c(128, 171, 385, 417, 464, 469, 548, 551, 563, 570, 571, 573, 575, 580, 609, 752, 776, 807, 929, 935, 1315, 1368, 1400, 1427, 1474, 1477, 1496, 1497, 1498, 1515, 1520, 1521, 1548, 1553, 1554, 1555, 1556, 1557, 1558, 1572, 1573, 1591, 1592, 1656, 1673, 1766, 2258, 2307, 2693, 2742, 2775, 2790, 2797, 2897, 2898, 2927, 2931, 2932, 2933, 2934, 2935, 2956, 2980, 2981, 2995, 3000, 3001, 3034, 3041, 3053, 3090, 3091, 3092, 3134, 3146, 3298, 3302, 3307, 3354, 3358, 3378, 3379, 3407, 3410, 3414, 3557, 3562, 3716, 3724, 3776, 3961, 3963, 3964, 3966, 3968, 4020)
#Plate7
Plate7 <- c(2582,407,444,476,576,597,601,610,611,854,925,951,1363,1369,1393,1466,1467,1468,1469,1470,1471,1472,1473,2423,2432,2435,2457,2458,2474,2479,2515,2516,2579,2628,2670,2677,2726,2745,2746,2748,2776,2777,2778,2779,2780,2786,2798,2833,2866,2913,2965,2966,2967,2974,2979,3025,3026,3027,3028,3029,3030,3031,3032,3033,3035,3036,3037,3038,3039,3048,3051,3052,3076,3077,3078,3080,3114,3115,3121,3122,3123,3124,3125,3126,3127,234,3133,3324,3325,3326,3331,3332,3333,3334,3338,3339)
#Plate8 
Plate8 <- c(3340,3341,3342,3343,3373,3392,3397,3402,3403,3464,3465,3491,3492,3493,3497,3498,3499,3500,3501,3510,3511,3512,3513,3514,3515,3516,3520,3521,3522,3523,3527,3528,3529,3778,3779,120,329,332,336,337,338,347,872,993,1398,2658,2698,2853,3242,3256,3309,3393,3740,3745,3751,3752,4001,4002,4015,4017,4034,47,50,2471,2499,2602,3138,3272,3759,3777,3879,119,191,436,437,438,442,558,1653,233,236,239,324,343,353,363,367,368,400,401,622,963,1001,1004,1007,1010)
#Plate9
Plate9 <- c(1013,1017,1020,1023,1030,1033,1037,1040,1044,1049,1052,1055,1060,1066,1069,1072,1076,1079,1082,1085,1089,1093,1097,1102,1106,1111,1114,1118,1122,1126,1129,1132,1138,1141,1147,1156,1159,1163,1167,1171,1175,1179,1183,1186,1190,1193,1196,1199,1202,1205,1208,1211,1215,1218,1221,1225,1231,1235,1240,1243,1246,1249,1254,1257,1260,1264,1270,1274,1277,1280,1283,1286,1289,1292,1298,1304,1308,1311,1314,1318,1321,1333,1336,1339,2732,2736,3306,3546,3549,3552,3997,4045,4051,4063,4068,4081)
#ZYGO plate
Plate <- c(1416, 128, 385, 417, 464, 563, 573, 580, 1400, 1427, 1520, 1521, 1553, 1554, 1556, 1557, 1558, 1572, 1573, 1591, 1592, 1766, 2790, 2927, 2931, 2932, 2933, 2934, 2995, 3090, 3091, 3092, 3146, 3302, 3307, 3378, 3379, 3407, 3410, 3414, 3724)
#DE plate
# Plate <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15")  # old
Plate <- c("1A", "1B", "1C", "1D", "1E", "1F", "1G", "1H", "2A", "2B", "2C", "2D", "2E", "2F", "2G")  # uncertain, as sequencing man isnt certain...

plate_strains <- data.frame(Plate1, Plate2, Plate3, Plate4, Plate5, Plate6, Plate7, Plate8, Plate9)
plate_strains[3, 1] # row, column

#Put here name of plate, for purposes of directories
PlateNumber <- "1"
IncludeMasked <- "No"  # if =="no", then remove masked SNPs (predicted)
#PlateNumber <- "ZYGO"
#PlateNumber <- "DE"
cigar_stuff <- "_CIGAR" # using CIGAR stuff
#cigar_stuff <- ""  # normal
#FOR PYTHON CODE

###Use this (each plate...), then use SNP_and_MAF_finding.py.  Then MatrixMaker.py
plates_to_do <- c(1,2,3,4,5,6,7)
plate_strains <- data.frame(Plate1, Plate2, Plate3, Plate4, Plate5, Plate6, Plate7, Plate8, Plate9)

for (ikk in plates_to_do){
  print(paste("Now on plate: ", ikk))
Plate <- plate_strains[,ikk]
for (i in Plate){
  Variablename <- paste("All_VCF\\Plate", ikk, "VCF", cigar_stuff, "\\NCYC",i,"_freebayes_SNP_genome.vcf", sep = '')
  vcf <- readVcf(Variablename)
  vcf1 <- vcf[lapply(info(vcf)$TYPE,length)==1]     #Only single-nucleotide changes
  vcf2 <- vcf1[as.character(info(vcf1)$TYPE)=="snp"]#Only those classed as SNPs
  vcf3 <- vcf2[rowData(vcf2)$QUAL > 30]             #Only those of quality >30
  PercentReference <- (data.frame(info(vcf3)$AO, info(vcf3)$RO))[c(3,4)] #AO = number alternate alleles, RO= number reference alleles
  ToPrint <- paste("Counter (1-96): ", match(i, Plate)," NCYC Number: ", i, " ORF matches: ", dim(vcf3)[1])
  print(ToPrint)
  output <- paste(rownames(vcf3), PercentReference$value, PercentReference$info.vcf)
  write(output, file = paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Plate", ikk, "VCF", cigar_stuff, "_AfterR\\NCYC",i,"_freebayes_SNP_genome_R.txt", sep = ''))
  
}
}


## same, but for ZYGO 'plate'
for (i in Plate){
  Variablename <- paste("All_VCF\\Zygomyces_VCFs\\NCYC",i,"_freebayes_SNP_genome.vcf", sep = '')
  vcf <- readVcf(Variablename)
  vcf1 <- vcf[lapply(info(vcf)$TYPE,length)==1]     #Only single-nucleotide changes
  vcf2 <- vcf1[as.character(info(vcf1)$TYPE)=="snp"]#Only those classed as SNPs
  vcf3 <- vcf2[rowData(vcf2)$QUAL > 30]             #Only those of quality >30
  PercentReference <- (data.frame(info(vcf3)$AO, info(vcf3)$RO))[c(3,4)] #AO = number alternate alleles, RO= number reference alleles
  ToPrint <- paste("Counter (1-96): ", match(i, Plate)," NCYC Number: ", i, " ORF matches: ", dim(vcf3)[1])
  print(ToPrint)
  output <- paste(rownames(vcf3), PercentReference$value, PercentReference$info.vcf)
  write(output, file = paste("C:\\Users\\Joe\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\Zygomyces_VCFs_AfterR\\NCYC",i,"_freebayes_SNP_genome_R.txt", sep = ''))
  
}
## same, but for DE 'plate'
for (i in Plate){
  Variablename <- paste("All_VCF\\PlateDEVCF", cigar_stuff,"\\",i,"_freebayes_SNP_genome.vcf", sep = '')
  vcf <- readVcf(Variablename)
  vcf1 <- vcf[lapply(info(vcf)$TYPE,length)==1]     #Only single-nucleotide changes
  vcf2 <- vcf1[as.character(info(vcf1)$TYPE)=="snp"]#Only those classed as SNPs
  vcf3 <- vcf2[rowData(vcf2)$QUAL > 30]             #Only those of quality >30
  PercentReference <- (data.frame(info(vcf3)$AO, info(vcf3)$RO))[c(3,4)] #AO = number alternate alleles, RO= number reference alleles
  ToPrint <- paste("Counter (1-15): ", match(i, Plate)," NCYC Number: ", i, " ORF matches: ", dim(vcf3)[1])
  print(ToPrint)
  output <- paste(rownames(vcf3), PercentReference$value, PercentReference$info.vcf)
  write(output, file = paste("C:\\Users\\Joseph\\Desktop\\PhD_stuff\\Jo_Stuff\\Programs_for_joShare\\PlateDEVCF", cigar_stuff, "_AfterR\\NCYC", i,"_freebayes_SNP_genome_R.txt", sep = ''))
  
}



#END HERE



#Find partial SNP values? Must be a way to extract... Look at documentation?





count12 <- 0
for (i in elementNROWS(allele2)){
  if (res$map[["ignore"]][i] == FALSE){
    #print("These alleles dont.")
    #print(allele2[[i]])
    #print(res$map[["ignore"]][i])
  }
  else {
    count12 <- count12 + 1
    print("These alleles count.")
  }
  print(count12)
  }


for (i in Plate1){
  Variablename <- paste("All_VCF\\NCYC",i,"_freebayes_SNP_genome.vcf", sep = '')
  VCF <- readVcf(Variablename)
  MatrixVCF <- genotypeToSnpMatrix(VCF)
  head(MatrixVCF)
  allele2 <- MaxtrixVCF$map[["allele.2"]]
  
  for (k in length(elementNROWS(allele2))){
   for (l in c(1,96) ){
     #Check each VCF file's list of SNPs, if the SNP (maybe even just match gene?) is found in all genome matrices, save it.
     #Later use this SNP list to create SNP matrix files (A -> 1000, G -> 0100, ...)
     #In the end, file with 96 genomes and their list of SNPs. save a file with SNP names (or header?).
   }
     
  }
  
}


vcf1 <- vcf[lapply(info(vcf)$TYPE,length)==1]     #Only single-nucleotide changes
vcf2 <- vcf1[as.character(info(vcf1)$TYPE)=="snp"]#Only those classed as SNPs
vcf3 <- vcf2[rowData(vcf2)$QUAL > 30]             #Only those of quality >30

dim(vcf3)                                         #number of SNPs?
median(info(vcf3)$DP)                             #DP is read-depth
sd(info(vcf3)$DP)                                 #calculate SD of read-depth
length(info(vcf3)$AB[as.numeric(info(vcf3)$AB)==0]) #Number of SNPs


for (jk in c(1:10)){
  print(geno(vcf3)$GT[jk])
  rownames(vcf3)
}
print(vcf3)
info(vcf3)
potato <- (vcf3)


# test
vcf <- readVcf("All_VCF\\NCYC1_freebayes_SNP_genome.vcf")
header(vcf)
res <- genotypeToSnpMatrix(vcf)
head(res)
allele2 <- res$map[["allele.2"]]
## number of alternate alleles per variant
unique(elementNROWS(allele2))

head(allele2)

res$map

print(length(elementNROWS(allele2)))

"7-EC1118_1F14_0089g" %in% res$map[["snp.names"]]






