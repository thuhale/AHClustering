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

compare = read.csv("Data/National_Downloadable_File.csv", stringsAsFactors = F, sep = ",", row.names = NULL)  ##2014 data
compare = compare[,c(1:27)]
compare = compare[nchar(compare$Zip.Code) >= 5,]


#None patient-facing specialities
spec.excl = c("MASS IMMMUNIZATION ROSTER BILLER", 
              "ANESTHESIOLOGY", "PHYSICAL THERAPY", "RADIOLOGIST", 
"CERTIFIED NURSE MIDWIFE", "CERTIFIED REGISTERED NURSE ANESTHETIST", 
"ANESTHESIOLOGY ASSISTANT", "CHIROPRACTIC", "CLINICAL NURSE SPECIALIST", 
"CLINICAL PSYCHOLOGIST", "CLINICAL SOCIAL WORKER", 
"OCCUPATIONAL THERAPY", "NURSE PRACTITIONER", "PATHOLOGY", 
"PSYCHOLOGIST BILLING INDEPENDENTLY", "REGISTERED DIETITIAN OR
NUTRITION PROFESSIONAL", "SPEECH LANGUAGE PATHOLOGIST")


compare = compare[(!(compare$Primary.specialty %in% spec.excl))
& (!(compare$Secondary.specialty.1 %in% spec.excl)), ]

compare$zip = substr(compare$Zip.Code,1,5)

## remove organization
ppes = fread("Data/npidata_20050523-20160911.csv")
ppes = as.data.frame(ppes)
col = c(1:2,25)
ppes = ppes[,col]
colnames(ppes) = c("NPI", "EntityType", "zip2")
ppes = ppes[ppes$EntityType == 1,]


compare = compare[compare$NPI %in% ppes$NPI, ]
rm(ppes);
short = compare[, c("NPI", "Primary.specialty", "zip")]
rm(compare)
#Takeout PCP
pcp_spec = c("FAMILY PRACTICE", "GENERAL PRACTICE",
             "GERIATRIC MEDICINE", "OSTEOPATHIC MANIPULATIVE MEDICINE",
             "PEDIATRIC MEDICINE", "PREVENTATIVE MEDICINE", "HOSPICE/PALLIATIVE CARE")
PCP = short[short$Primary.specialty %in% pcp_spec,]
write.csv(PCP, "PCP-2014.csv", row.names = F)

#read_in referral data
ref = read_csv(file = "Data/physician-shared-patient-patterns-2014-days180.txt", col_names = F,col_types = "cciii")

colnames(ref) = c("NPI", "NPI2", "Ties", "Unique", "Sameday")
ref = ref[ref$NPI %in% short$NPI & ref$NPI2 %in% short$NPI,]
g = graph.edgelist(as.matrix(ref[,1:2]))
E(g)$weight = ref$Unique

# Create PCP Network
W = get.adjacency(g, attr = 'weight', sparse = T)

B = W[rownames(W) %in% PCP$NPI,]
B = B[, !(colnames(W) %in% PCP$NPI)] ##B is the bi-partite matrix
A = B%*%t(B) # A is PCP-referral adjacency matrix
diag(A) = 0 ## remove self-loop

## Regularize A
O = rowSums(A)
tau = rep(mean(O), length(O))

L_reg = Diagonal(length(O), 1/sqrt(O+tau))%*%A%*%Diagonal(length(O), 1/sqrt(O+tau))

#Make PCP referral graph
gref = graph.adjacency(L_reg, weighted = T, diag = F)
core = coreness(gref)

gsm = induced_subgraph(gref,core>300)

Els = as_edgelist(gsm)
Els = data.frame(Els, E(gsm)$weight,stringsAsFactors = F)
colnames(Els) = c("NPI", "NPI2", "Weight")
write.csv(Els, "PCP_Referral_List.csv", row.names = F)


