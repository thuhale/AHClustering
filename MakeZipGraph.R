rm(list = ls())
library(data.table)
library(readr)
library(igraph)
setwd("/Users/lehathu/Desktop/AHC")
library(Matrix)
library(sparseAHC)
library(MCL)
library(dplyr)
library(zipcode)
data(zipcode); 
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]

compare1 = read.csv("Data/National_Downloadable_File.csv", stringsAsFactors = F, sep = ",", row.names = NULL)
compare1 = compare1[,c(1:27)]
compare1 = compare1[nchar(compare1$Zip.Code) >= 5,] #remove 1 obs

compare2 = read.csv("Data/Refresh_Data_Archive_November_2015/National_Downloadable_File.csv", stringsAsFactors = F, sep = ",", row.names = NULL, colClasses = "character")
compare2 = compare2[,c(1:27)]


#None patient-facing specialities
spec.excl = c("MASS IMMMUNIZATION ROSTER BILLER", 
              "ANESTHESIOLOGY", "PHYSICAL THERAPY", "RADIOLOGIST", 
"CERTIFIED NURSE MIDWIFE", "CERTIFIED REGISTERED NURSE ANESTHETIST", 
"ANESTHESIOLOGY ASSISTANT", "CHIROPRACTIC", "CLINICAL NURSE SPECIALIST", 
"CLINICAL PSYCHOLOGIST", "CLINICAL SOCIAL WORKER", 
"OCCUPATIONAL THERAPY", "NURSE PRACTITIONER", "PATHOLOGY", 
"PSYCHOLOGIST BILLING INDEPENDENTLY", "REGISTERED DIETITIAN OR
NUTRITION PROFESSIONAL", "SPEECH LANGUAGE PATHOLOGIST")

compare1 = compare1[(!(compare1$Primary.specialty %in% spec.excl))
& (!(compare1$Secondary.specialty.1 %in% spec.excl)), ]

compare2 = compare2[(!(compare2$Primary.specialty %in% spec.excl))
& (!(compare2$Secondary.specialty.1 %in% spec.excl)), ]

## remove organization
ppes = fread("Data/npidata_20050523-20160911.csv")
ppes = as.data.frame(ppes)
col = c(1:2,25)
ppes = ppes[,col]
colnames(ppes) = c("NPI", "EntityType", "zip2")
ppes = ppes[ppes$EntityType == 1,]


compare1 = compare1[compare1$NPI %in% ppes$NPI, ]
compare2 = compare2[compare2$NPI %in% ppes$NPI, ]

compare1$zip = substr(compare1$Zip.Code, 1,5)
compare2$zip = substring(compare2$Zip.Code, 1,5)


short1 = compare1[,c("NPI", "zip")]
short2 = compare2[,c("NPI", "zip")]
removeDuplicates = function(zipp){
  setnames(zipp, 1:2, c("NPI", "zip"))
  setkey(zipp)  # setkey to both columns.  this ensures "unique" looks at all columns.
  zipp = unique(zipp)
  # remove any row that has NPI repeated **including first**
  # thanks!  http://stackoverflow.com/questions/16265808/identify-duplicates-and-mark-first-occurrence-and-all-others
  remt = duplicated(zipp,by="NPI") | duplicated(zipp,by="NPI", fromLast=TRUE)
  zipp = zipp[!remt,]
  if(anyDuplicated(zipp, by = "NPI")>0) print("removeDuplicates didn't work")  # should be zero! 
  return(zipp)
}
short = removeDuplicates(rbindlist(list(short1,short2)))

short = short[short$zip %in% zipcode$zip, ]
write.csv(short, "NPI_Zip_2014_2015.csv")

ref = read_csv(file = "Data/physician-shared-patient-patterns-2014-days180.txt", col_names = F,col_types = "cciii")

colnames(ref) = c("NPI", "NPI2", "Ties", "Unique", "Sameday")
ref = ref[ref$NPI %in% short$NPI & ref$NPI2 %in% short$NPI,]

ref = merge(ref, short, by = "NPI", all.x = T)
ref = merge(ref, short, by.x = "NPI2", by.y = "NPI", all.x = T)
colnames(ref)[6:7] = c("zip", "zip2")
ref = ref[, c(2,1,3,4,5,6,7)]

#Make the graph
g = graph.edgelist(as.matrix(ref[,6:7]), directed = F)
E(g)$weight = ref$Unique
g = simplify(g)

W = get.adjacency(g, attr = 'weight', sparse = T)
## Regularize W
O = rowSums(W)
tau = rep(mean(O), length(O))

L_reg = Diagonal(length(O), 1/sqrt(O+tau))%*%W%*%Diagonal(length(O), 1/sqrt(O+tau))

g = graph.adjacency(L_reg, weighted = T, diag = F)
save(g, file = "Data/Zip_Referral_graph.RData")
###




