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
data(zipcode)
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]

compare = read.csv("Data/National_Downloadable_File.csv", stringsAsFactors = F, sep = ",", row.names = NULL)  ##2014 data
compare = compare[,c(1:27)]
compare = compare[nchar(compare$Zip.Code) >= 5,]
compare$zip = substr(compare$Zip.Code,1,5)

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


## remove organization
ppes = fread("Data/npidata_20050523-20160911.csv")
ppes = as.data.frame(ppes)
col = c(1:2,25)
ppes = ppes[,col]
colnames(ppes) = c("NPI", "EntityType", "zip2")
ppes = ppes[ppes$EntityType == 1,]

compare = compare[compare$NPI %in% ppes$NPI, ]

## Only consider NPIs with zipcode within US continent
compare = compare[compare$zip %in% zipcode$zip,]


#Make a table of doctors with zipcode (unique)
short = compare[, c("NPI", "zip")]
short = short[short$zip %in% zipcode$zip,]
short = unique(short)

tab = as.data.frame(table(short$NPI))
colnames(tab) = c("NPI", "Freq")
short1 = short[short$NPI %in% tab$NPI[tab$Freq == 1],] ## short1 is doctors with only one zipcode
short2 = short[short$NPI %in% tab$NPI[tab$Freq == 2],] ## doctor with 2 zipcodes
short3 = short[short$NPI %in% tab$NPI[tab$Freq >2],] ## doctor with 3 zipcodes

getFirst = function(v){return(v[1])}

getMode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
short2 = aggregate(short2$zip, by = list(short2$NPI),FUN = getFirst )
colnames(short2) = c("NPI", "zip")

short3 = aggregate(short3$zip, by = list(short3$NPI),FUN = getMode )
colnames(short3) = c("NPI", "zip")

short = rbind(short1, short2, short3)
write.table(short, "PCP_zip.csv", quote = T, row.names = F, sep = ",")
#For the doctors that have zip information, we have to consolidate zip of the doctor
pcp_spec = c("FAMILY PRACTICE", "GENERAL PRACTICE",
             "GERIATRIC MEDICINE", "OSTEOPATHIC MANIPULATIVE MEDICINE",
             "PEDIATRIC MEDICINE", "PREVENTATIVE MEDICINE", "HOSPICE/PALLIATIVE CARE")
PCP = compare[compare$Primary.specialty %in% pcp_spec,]

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
hist(core)
gsm = induced_subgraph(gref,core>300)
M = get.adjacency(gsm,attr = "weight", sparse = T)

hclust = sparseAHC(M, linkage = "single")

x = sample(V(gsm), 1000)
y = induced_subgraph(gsm, x)


Els = as_edgelist(gsm)
Els = data.frame(Els, E(y)$weight,stringsAsFactors = F)
colnames(Els) = c("NPI", "NPI2", "Weight")
write.csv(Els, "PCP_Referral_List_sample.csv", row.names = F)

#create a merge matrix
order = hclust$order
merge = hclust$merge
sorted_order = data.frame(colnames(M), 1:length(colnames(M)), order)
colnames(sorted_order) = c("NPI", "index", "order")
sorted_order = merge(sorted_order, short, by = "NPI", all.x = T)
sorted_order = merge(sorted_order, zipcode, by = "zip", all.x = T)
sorted_order = sorted_order[,1:6]
sorted_order = sorted_order[order(sorted_order$index),]
write.csv(sorted_order, "sorted_order.csv", row.names = F)
write.csv(merge, "merge.csv", row.names = F)

