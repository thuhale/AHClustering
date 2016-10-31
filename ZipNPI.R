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

compare1 = read.csv("Data/National_Downloadable_File.csv", stringsAsFactors = F, sep = ",", row.names = NULL)
compare1 = compare1[,c(1:27)]
compare1 = compare1[nchar(compare1$Zip.Code) >= 5,]

compare2 = read.csv("Data/Refresh_Data_Archive_November_2015/National_Downloadable_File.csv", stringsAsFactors = F, sep = ",", row.names = NULL)
compare2 = compare2[,c(1:27)]
compare2 = compare2[nchar(compare2$Zip.Code) >= 5,]


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

zipp = rbind(short1, short2)
zipp = unique(zipp)
df1 <- aggregate(zip ~ NPI, data = zipp, FUN = length)
doc = df1$NPI[df1$zip==1]
zipp = zipp[zipp$NPI %in% doc,]


data(zipcode); 
zipcode = as.data.table(zipcode); 
setkey(zipcode, zip)
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH")
zipcode = zipcode[!(zipcode$state %in% peri),]

zipp = zipp[zipp$zip %in% zipcode$zip, ]
zipp = merge(zipp, zipcode, by = "zip", all.x = T)
write.csv(zipp, "NPIZip.csv", row.names = F)


