rm(list = ls())
library(zipcode)
data(zipcode)
library(geosphere)


## SIDE DATA
peri = c("AA", "AE", "AK", "AP","AS", "FM", "GU", "HI", "MP", "PR", "PW", "VI", "MH") ## remove non-us continential 
zipcode = zipcode[!(zipcode$state %in% peri),]
short = read.csv("Public/NPI_zip.csv", colClasses = "character")

#DATA CREATED FROM "MakeNPIGraph.R"
load(file = "Data/M_2014.RData")
npi = data.frame(colnames(M))
colnames(npi) = c("NPI")
npi = merge(npi, short, by ="NPI", all.x = T)
coreZip = data.frame(unique(npi$zip),stringsAsFactors = F)
colnames(coreZip) = "zip"
coreZip = merge(coreZip, zipcode, by = "zip", all.x = T)

# Assign the non-core zipcode to 5 closest zipcode in geodesic distance. 

excl = zipcode[!(zipcode$zip %in% npi$zip), ]
excl$zip1 = '00000'
excl$zip2 = '00000'
excl$zip3 = '00000'
excl$zip4 = '00000'
excl$zip5 = '00000'
print("Process started")
for(i in 1:dim(excl)[1]){
	tmp = excl[1,]
	loc = c(tmp$longitude, tmp$latitude)
	m = apply(cbind(coreZip$longitude, coreZip$latitude), 1, function(p, loc) distCosine(p,loc), loc = loc)
	l = order(m)
	excl$zip1[i] = coreZip$zip[l[1]]
	excl$zip2[i] = coreZip$zip[l[2]]
	excl$zip3[i] = coreZip$zip[l[3]]
	excl$zip4[i] = coreZip$zip[l[4]]
	excl$zip5[i] = coreZip$zip[l[5]]
	if(i%%1000==0){print(i)}
}
write.csv(excl, "Public/KnnZip.csv", row.names = F)



