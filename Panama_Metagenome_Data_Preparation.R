###################################
## R script for data import and analysis
## Panama metagenome article
##

# G. Lajoie
# S. Kembel

### optional - uncomment to run - install required libraries
#install.packages("ade4", dependencies = TRUE)
#install.packages("reshape2", dependencies = TRUE)
#install.packages("vegan", dependencies = TRUE)
#install.packages("splitstackshape", dependencies = TRUE)
#install.packages("phytools", dependencies = TRUE)
#install.packages("plyr", dependencies = TRUE)
#install.packages("dplyr", dependencies = TRUE)
#install.packages("candisc", dependencies = TRUE)
#install.packages("ape", dependencies = TRUE)
#install.packages("adespatial", dependencies = TRUE)
#install.packages("entropart", dependencies = TRUE)
#install.packages("betapart", dependencies = TRUE)


### load required libraries

library(ade4) # Load ade4 before vegan
library(reshape2)
library(vegan)
library(splitstackshape)
library(phytools)
library(plyr)
library(dplyr)
library(candisc)
library(ape)
library(adespatial)
library(entropart)
library(zoo)
library(stringr)

############################
## CREATE METADATA OBJECT ##
############################

# Sample names
samplenames <- as.character(read.table("/data/shared/panama_metagenome/samplenames.txt", header=FALSE)[,1])
## taxonomy metadata
metadata.taxo <- read.csv("/data/shared/panama_metagenome/metadata/sample-taxonomy-metadata.csv")
rownames(metadata.taxo) <- metadata.taxo$sequencing.samplename
## host traits
metadata.traits <- read.csv("/data/shared/panama_metagenome/metadata/hosttraits.sub.csv")

# Create single metadata file
metadata <- merge(metadata.taxo, metadata.traits, by.x="panama.bact.sampleID", by.y="X", all.x=TRUE)
rownames(metadata) <- metadata$sequencing.samplename

# FINAL OBJECT
metadata <- metadata[rownames(metadata.taxo),]

##################################
## IMPORT TAXONOMIC DATA ##
##################################

reads<-c('R1','R2')

# Import taxonomic annotations
res.t <- NULL
for (i in samplenames) {
  for (j in reads){
    # for each sample
    # make filename fn by pasting sample name i20nto path
    # change details to load different taxonomic levels or methods
    fn <- paste("/data/shared/panama_metagenome/kaiju3/kaiju/bin/", i, "-greedy5-kaiju-names-", j, ".out", sep="")
    # read file fn (kaiju taxon summary file) into data.frame named samp
    samp <- read.table(fn, sep="\t", quote="", header=F, fill=T, col.names=paste("V", 1:7, sep="")) # Specifying colnames string to set the number of columns to 7.
    # add colnames as headers
    colnames(samp)[2] <- c("read.num")
    # create column named samp in object samp and fill it with the name of this sample
    samp$samp <- i
    samp$read <- j
    # add current sample to result object
    res.t <- rbind(res.t, samp)
  }
}

res.t[res.t==""]<-NA # Fill empty cells with NAs

taxan<-NULL

# Simple merging among R1 and R2
for (i in samplenames){
  u1<-res.t[which(res.t$samp==i&res.t$read=='R1'),]
  u2<-res.t[which(res.t$samp==i&res.t$read=='R2'),]
  u<-merge(u1,u2, by=c('read.num'))
  taxan<-rbind(taxan,u)
}

taxo.merge<-taxan
# Change "NA" by NA
taxo.merge[taxo.merge=="NA"]<-NA


###############################################
## IMPORT METAGENOMIC DATA ##
###############################################

# Import functionally annotated sequences

# Trimmed

res <- NULL
for (i in samplenames) {
  for (j in reads){
    # for each sample
    # make filename fn by pasting sample name into path
    # change details to load different taxonomic levels or methods
    #fn <- paste("/data/shared/panama_metagenome/cognizer_trimmed/fixedName/", i, "-q10-trimmed_", j, "CogName.txt", sep="")
    # fn <- paste("/data/shared/panama_metagenome/panama_2018/function/pair/pear_assembled/", i, "-q10-trimmed.assembled_cog2/assignments.txt", sep="")
    fn <- paste("/data/shared/panama_metagenome/panama_2018/function/trim/", i, "-q10-trimmed_", j, "_cog2/assignments.txt", sep="")
    # read file fn (kaiju taxon summary file) into data.frame named samp
    samp <- read.table(fn, skip=2, sep="\t", quote="", header=F, fill=T)
    # add colnames as headers
    colnames(samp) <- c("read.num",	"cog.group", "COG.id",	"Kegg.id",	"PF.id",	"FIG.id",	"GO.id",	"cog.gene",	"kegg.gene",	"pf.gene",	"fig.gene",	"GO.gene")
    # create column named samp in object samp and fill it with the name of this sample
    samp$samp <- i
    samp$read <- j
    # add current sample to result object
    res <- rbind(res, samp)
  }
}

res$read.num<-sapply(str_split(res$read.num," "),'[',1)# Keeping only the portion of the name before the space. The final portion of initial name is only indicative of the read number.

# Change "NA" by NA
res[res=="NA"]<-NA

# Keep only annotated sequences
taxan<-NULL

# Simple merging of R1 and R2
for (i in samplenames){
  u1<-res[which(res$samp==i&res$read=='R1'),]
  u2<-res[which(res$samp==i&res$read=='R2'),]
  u<-merge(u1,u2, by=c('read.num'))
  taxan<-rbind(taxan,u)
}

# Change "-" for NA
function.merge<-taxan
function.merge[function.merge=="-"]<-NA

###############################
## MERGE AND CLEAN DATASETS ##
##############################

taxo.fct<-merge(taxo.merge, function.merge, by=c('read.num'), all=T)

# Reduced file, keeping only a reduced set of columns
tf<-taxo.fct[, which(colnames(taxo.fct)%in%c('read.num','V1.x','V7.x','read.x.x','V1.y','V7.y','read.y.x','samp.x.x','Kegg.id.x','read.x.y','Kegg.id.y','read.y.y'))]

#####
# Clean taxonomic annotations
#####

# Split V7.x (taxonomic annotations of read 1) by taxonomic level
tf<- tf %>% separate(V7.x, c("V7_1.x","V7_2.x","V7_3.x","V7_4.x","V7_5.x","V7_6.x","V7_7.x"), sep=';')

# Split V7.y (taxonomic annotations of read 2) by taxonomic level
tf<- tf %>% separate(V7.y, c("V7_1.y","V7_2.y","V7_3.y","V7_4.y","V7_5.y","V7_6.y","V7_7.y"), sep=';')

# Objects
colnum.x<-which(colnames(tf)%in%c("V7_1.x","V7_2.x","V7_3.x","V7_4.x","V7_5.x","V7_6.x","V7_7.x"))
colnum.y<-which(colnames(tf)%in%c("V7_1.y","V7_2.y","V7_3.y","V7_4.y","V7_5.y","V7_6.y","V7_7.y"))

# Change "NA" by NA
tf[tf==" NA"]<-NA

# Remove NA lines to reduce loop calculations (from 7.3 million to 2.1 million lines)  (keep only sequences that have at least one taxonomically annotated read)
nonna<-as.data.frame(cbind(is.na(tf[,colnum.x[1]]),is.na(tf[,colnum.y[1]])))
nonna.sub<-which(nonna$V1==F&nonna$V2==F)
tf<-tf[nonna.sub,]

# True/False Match Matrix, turn it into characters to sort F/T before Nas
code<-as.data.frame(tf[,colnum.x]==tf[,colnum.y])
code[] <- lapply(code, as.character)
code[is.na(code)]<-'NA'

# Add output columns
tf<-cbind(tf, as.data.frame(matrix(NA,nrow(tf),7)))
colnum.out<-c((ncol(tf)-6):ncol(tf))

# Import "True" matches into matrix
true.code<-which(code=='TRUE',arr.in=TRUE) # obtain index of "TRUE"s in row and column form rather than vector
tf[,colnum.out][true.code]<-tf[,colnum.x][true.code] # Add taxonomy to matching scenarios

# Create matrix with NA matches (1 is identified and the other not, or both are not identified)
x.code<-which(is.na(tf[,colnum.x])==F&is.na(tf[,colnum.y])==T, arr.in=T) # missing data in y
y.code<-which(is.na(tf[,colnum.x])==T&is.na(tf[,colnum.y])==F, arr.in=T) # missing data in x

# Import NA matches into matrix
tf[,colnum.out][x.code]<-tf[,colnum.x][x.code]
tf[,colnum.out][y.code]<-tf[,colnum.y][y.code]

# Manage taxonomy for F scenarios
dd<-which(code == 'FALSE', arr.in=T) # index 'FALSE's
ff<-aggregate(col~row,dd,FUN=function(x) c(min(x):7)) # minimum FALSE per row
reps<-aggregate(col~row,dd,FUN=function(x) length(min(x):7)) # length of replacement vector
rep.vecs<-cbind(rep(ff[,1],reps[,2]),unlist(ff[,2])) # expand matrix of replacement cells

# Import F replacements into matrix
tf[,colnum.out][rep.vecs]<-NA

#####
# Clean functional annotations
#####

# Approach 1: Keeping unique functional genes from R1 and R2 reads for each sample (under the variable name $Kegg.uni) # Attention, très long!!
tf$Kegg.uni<-rep(NA,nrow(tf))

for (j in 1:nrow(tf)){
  
  if (is.na(tf$Kegg.id.x[j])==F&is.na(tf$Kegg.id.y[j])==F){
    tf$Kegg.uni[j]<-paste(unique(unlist(list(str_split(tf$Kegg.id.x[j],', '), str_split(tf$Kegg.id.y[j], ', ')))), collapse=", ") # Choosing to collapse on the Kegg.id, could choose another functional annotation
    next
  }
  if (is.na(tf$Kegg.id.x[j])==F&is.na(tf$Kegg.id.y[j])==T){
    tf$Kegg.uni[j]<-as.character(tf$Kegg.id.x[j])
    next
  }
  if (is.na(tf$Kegg.id.x[j])==T&is.na(tf$Kegg.id.y[j])==F){
    tf$Kegg.uni[j]<-as.character(tf$Kegg.id.y[j])
    next
  }
}

#write.csv(tf,'tf.csv') # Change path

######################
## RAREFYING FUNCTIONAL DATASET ##
####################

funtax.bact<-tf # Written to file tf.csv in code_data
# Keep only simplified dataset
funtax.bact<-funtax.bact[,which(colnames(funtax.bact)%in%c('read.num','samp.x.x','V1','V2','V3','V4','V5','V6','V7','Kegg.uni'))]

# Removing Eukaryota/Archaea/Virus/Unannotated sequences
funtax.bact<-funtax.bact[which(funtax.bact$V1%in%c('Bacteria')),] # See which taxonomic "type" to filter on (trimmed or paired, or both)

### Keeping only bacterial sequences with functional annotations
funtax.rare<-funtax.bact
funtax.rare<-funtax.rare[which(is.na(funtax.rare$Kegg.uni)==F),]
minseq<-min(aggregate(read.num~samp.x.x, data=funtax.rare, length)[,2]) # Min of sequences: 10864 # Max: 58510 

# Round to the lower hundred (doesn't work, do it by hand)
minseq<-10800

# Rarefy at minseq
rareseq<-NULL

for (i in samplenames){
  sub<-sort(sample(rownames(funtax.rare[which(funtax.rare$samp==i),]), minseq, replace=F))
  rareseq<-rbind(rareseq,funtax.rare[which(rownames(funtax.rare)%in%sub),])
}

rareseqfun<-rareseq # File added to the code_data folder
# write.csv(rareseqfun, 'rareseqfun.csv') # Change path

##############################
## CREATE TAX4FUN DATATABLE ##
##############################

# Tax4Fun reads already rarefied (see Tax4Fun script)

# Importing Tax4Fun annotations
tax4<-read.csv("/data/shared/panama_metagenome/panama_2018/16S_data/tax4fun-testing/tax4fun.gene.abunds.corr.rare.csv", header=T, row.names=1)

# Colnames represented only by the Kegg numbers
colnames(tax4)<-substr(colnames(tax4),1,6)

# Keeping samples also represented in the metagenomics dataset
tax4.sub<-tax4[rownames(tax4)%in%metadata$panama.bact.sampleID,]

# Renaming with metagenomics samplenames
tax4.sub<-tax4.sub[match(metadata$panama.bact.sampleID,rownames(tax4.sub)),]
rownames(tax4.sub)<-metadata$sequencing.samplename

# Numbers in the table represent the proportion of total trimmed sequences that have that given function.
# Corrected files represent relative abundances which have also been adjusted to reflect the copy number of the genes 
# (eg in theory they represent the relative abundance of the genomes of the organisms)

# FINAL OBJECT: relative abundance table
tax4.rel<-decostand(tax4.sub, method='total')


##################################
## CREATE METAGENOMIC DATATABLE ##
##################################

# Start with rareseqfun for rarefied sequences
# Importing Metagen annotations
# rareseqfun<-read.csv("/data/shared/panama_metagenome/code_data/rareseqfun.csv", header=T, row.names=1)

# Making community function table, merging functions from across sequences 
comag<-cSplit(rareseqfun, "Kegg.uni", sep=', ')
kegg.uni.names<-paste("Kegg.uni", sprintf('%0.2d', 1:17), sep='_') # 17 different columns in the split of kegg.uni
zz<-which(colnames(comag)%in%kegg.uni.names)
comag$weight<-apply(comag[,min(zz):max(zz)], 1, FUN=function(x) 1/length(which(is.na(x)==F))) # weight to give each kegg in each read

# Sample x function table
comag.m<-melt(comag, id.vars=c('samp.x.x','weight'), measure.vars=kegg.uni.names, na.rm=T)
comag.mat<-acast(comag.m, samp.x.x~value, value.var='weight', sum) # raw counts of occurrences of functions in each sample

# FINAL OBJECT: relative abundance table
meta.rel<-decostand(comag.mat, method='total')
meta.rel<-meta.rel[match(rownames(tax4.rel),rownames(meta.rel)),] # reorder meta.rel as tax4.rel


#################################
## CREATE TAXONOMIC DATATABLES ##
#################################

### From dataset rarefied from the pool of all sequences with taxonomic annotations

# Rarefy on taxonomically annotated sequences (funtax.rare)
funtax.rare<-funtax.bact
funtax.rare<-funtax.rare[which(is.na(funtax.rare$V2)==F),] # Keep only sequences for which there is at least a taxonomic annotation at the Phylum level
minseq<-min(aggregate(read.num~samp.x.x, data=funtax.rare, length)[,2]) # Min of sequences: 27166 (59_S9) Max: 159654 ## Second-to-last and third-to-last do not have duplicates, let's rarefy here.

# Round to the lower hundred
minseq<-27100

rareseq<-NULL

for (i in samplenames){
  sub<-sort(sample(rownames(funtax.rare[which(funtax.rare$samp==i),]), minseq, replace=F))
  rareseq<-rbind(rareseq,funtax.rare[which(rownames(funtax.rare)%in%sub),])
}

rareseqtax<-rareseq # File added to the code_data folder


### From dataset rarefied from the pool of sequences with both taxonomic AND functional annotation
# Use object rareseqfun
# rareseqfun<-read.csv("/data/shared/panama_metagenome/code_data/rareseqfun.csv", header=T, row.names=1)

########################################################
## ASSIGN HIERARCHICAL FUNCTIONAL CATEGORIES TO KEGGS ##
########################################################

# KO pre-cleaned only involves splitting rows into columns, starting from ko0000 hierarchy file downloaded from Kegg
hier<-read.csv('/data/shared/panama_metagenome/code_data/KO_precleaned.csv', header=F, sep=",")
hier<-hier[complete.cases(hier),]

res<-as.data.frame(matrix(NA,nrow(hier),5))
colnames(res)<-c('A','B','C','Kegg','Kegg.def')

for (i in 1:nrow(hier)){
  if (hier[i,1]=='A'){
    res[i,1]<-as.character(hier[i,3])
    next
  }
  if (hier[i,1]=='B'){
    res[i,2]<-as.character(hier[i,3])
    next
  }
  if (hier[i,1]=='C'){
    res[i,3]<-as.character(hier[i,3])
    next
  }
  if (hier[i,1]=='D'){
    res[i,4]<-as.character(hier[i,2])
    res[i,5]<-as.character(hier[i,3])
    next
  }
}

res[,c(1:3)]<-na.locf(res[,c(1:3)], na.rm=F)
res<-res[complete.cases(res),]

# Second step of cleaning

kegg.cat<-res

st<-duplicated(kegg.cat$Kegg)
stt<-duplicated(kegg.cat$Kegg, fromLast=T)
dc<-kegg.cat[which(st==T|stt==T),]
dc<-dc[order(dc$Kegg),]

subsub<-NULL

for (i in unique(dc$Kegg)){
  ff<-dc[which(dc$Kegg==i),]
  st<-duplicated(ff$B)
  stt<-duplicated(ff$B, fromLast=T)
  sub<-ff[which(st==T|stt==T),]
  subsub<-rbind(subsub,sub)
}

subsub2<-substr(subsub$C,1,10)
subsub<-cbind(subsub,subsub2)

tubtub<-NULL

for (i in unique(subsub$Kegg)){
  ff<-subsub[which(subsub$Kegg==i),]
  st<-duplicated(ff$subsub2)
  stt<-duplicated(ff$subsub2, fromLast=T)
  sub<-ff[which(st==T|stt==T),]
  tubtub<-rbind(tubtub,sub)
}

tubtub2<-tubtub[grep('BR:',tubtub$C),]

tubtub3<-tubtub[which(tubtub$Kegg%in%unique(tubtub2$Kegg)),]

# Keep |PATH| designations when both BR and PATH are present for a same function.

grep.rm<-grep('BR:',tubtub3$C)

tob<-rownames(tubtub3[grep.rm,])

kegg.cat2<-kegg.cat

kegg.cat2<-kegg.cat2[-which(rownames(kegg.cat2)%in%tob),]

# Deal with other BR designations: assign them to the nearer PATH

grep2<-grep('BR:', kegg.cat2$C)
length(unique(kegg.cat2$C[grep2]))

# Replace with alternatives PATH_BR

kegg.cat2$C[which(kegg.cat2$C=='Photosynthesis proteins [BR:ko00194]')]<-"Photosynthesis [PATH:ko00195]"
kegg.cat2$C[which(kegg.cat2$C=='Glycosylphosphatidylinositol (GPI)-anchored proteins [BR:ko00537]')]<-"Glycosylphosphatidylinositol(GPI)-anchor biosynthesis [PATH:ko00563]"
kegg.cat2$C[which(kegg.cat2$C=='Lipopolysaccharide biosynthesis proteins [BR:ko01005]')]<-"Lipopolysaccharide biosynthesis [PATH:ko00540]"
kegg.cat2$C[which(kegg.cat2$C=='Peptidoglycan biosynthesis and degradation proteins [BR:ko01011]')]<-"Peptidoglycan biosynthesis [PATH:ko00550]"
kegg.cat2$C[which(kegg.cat2$C=='Polyketide sugar unit biosynthesis [PATH:ko00523]')]<-"Polyketide biosynthesis proteins [BR:ko01008]"
kegg.cat2$C[which(kegg.cat2$C=='Spliceosome [BR:ko03041]')]<-"Spliceosome [PATH:ko03040]"
kegg.cat2$C[which(kegg.cat2$C=='Ribosome [BR:ko03011]')]<-"Ribosome [PATH:ko03010]"
kegg.cat2$C[which(kegg.cat2$C=='Ubiquitin mediated proteolysis [PATH:ko04120]')]<-"Ubiquitin system [BR:ko04121]"
kegg.cat2$C[which(kegg.cat2$C=='Proteasome [BR:ko03051]')]<-"Proteasome [PATH:ko03050]"
kegg.cat2$C[which(kegg.cat2$C=='DNA replication proteins [BR:ko03032]')]<-"DNA replication [PATH:ko03030]"
kegg.cat2$C[which(kegg.cat2$C=='ABC transporters [PATH:ko02010]')]<-"Transporters [BR:ko02000]"
kegg.cat2$C[which(kegg.cat2$C=='Bacterial secretion system [PATH:ko03070]')]<-"Secretion system [BR:ko02044]"
kegg.cat2$C[which(kegg.cat2$C=='Two-component system [BR:ko02022]')]<-"Two-component system [PATH:ko02020]"
kegg.cat2$C[which(kegg.cat2$C=='Cell adhesion molecules and their ligands [BR:ko04516]')]<-"Cell adhesion molecules (CAMs) [PATH:ko04514]"
kegg.cat2$C[which(kegg.cat2$C=='ABC transporters [PATH:ko02010]')]<-"Transporters [BR:ko02000]"
kegg.cat2$C[which(kegg.cat2$C=='Bacterial chemotaxis [PATH:ko02030]')]<-"Bacterial motility proteins [BR:ko02035]"

unique(kegg.cat2$C)

kegg.cat2<-kegg.cat2[-which(kegg.cat2$A%in%c('Human Diseases','Organismal Systems')),]
kegg.cat2<-droplevels(kegg.cat2)

# Strip extra space in kegg names
kegg.cat2$Kegg<-trimws(kegg.cat2$Kegg, which='right')
kegg.cat<-kegg.cat2

#############################################
## CREATE FUNCTIONAL CATEGORIES DATATABLES ##
#############################################

### Metagenomics datatables
data.cat <- meta.rel

# Melt kegg file
data.ks<-melt(as.matrix(data.cat))
data.ks<-as.data.frame(data.ks[complete.cases(data.ks),])
data.ks$Var2<-as.character(data.ks$Var2)

# Merge files by Kegg names
data.funcat<-merge(x=data.ks, y=kegg.cat, by.x=c('Var2'), by.y=c('Kegg'), all.x=T)
colnames(data.funcat)[1:2]<-c('Kegg.id','samp')

# Aggregating by Kegg categories
A.agg<-aggregate(value~Kegg.id+A+samp, data=data.funcat, FUN=unique)
A.agg<-aggregate(value~A+samp, data=A.agg, FUN=sum)
A.mat<-acast(samp~A, data=A.agg, fill=0)
A.mat.rel<-decostand(A.mat, 'total', 1) # Relative abundances

B.agg<-aggregate(value~Kegg.id+B+samp, data=data.funcat, FUN=unique)
B.agg<-aggregate(value~B+samp, data=B.agg, FUN=sum)
B.mat<-acast(samp~B, data=B.agg, fill=0)
B.mat.rel<-decostand(B.mat, 'total', 1) # Relative abundances

C.agg<-aggregate(value~Kegg.id+C+samp, data=data.funcat, FUN=unique)                             #Node~C+samp, data=funcat[-which(funcat$Kegg.id%in%c('K03781','K00249', 'K02298','K00356','K00341','K00336','K00382','K02160','K11263','K01615','K00208','K00123','K00257','K01209','K01869','K01958','K01995','K01998','K11754','K02890','K08348','K03781','K06044','K11263', 'K06861')),], FUN=flen) # 2 functions overly abundant in a few samples : K03781, K00249
C.agg<-aggregate(value~C+samp, data=C.agg, FUN=sum)
C.mat<-acast(samp~C, data=C.agg, fill=0)
# Remove traits with little overall abundance if wanted, e.g. C.mat<-C.mat[,which(colSums(C.mat)>=100)]
C.mat.rel<-decostand(C.mat, 'total', 1) # Relative abundances: USE!

# FINAL OBJECTS
# Metagenomics
meta.catA<-A.mat.rel
meta.catB<-B.mat.rel
meta.catC<-C.mat.rel

### tax4fun datatables
data.cat<-tax4.rel

# Melt kegg file
data.ks<-melt(as.matrix(data.cat))
data.ks<-as.data.frame(data.ks[complete.cases(data.ks),])

# Merge files by Kegg names
data.funcat<-merge(x=data.ks, y=kegg.cat, by.x=c('Var2'), by.y=c('Kegg'), all.x=T)
colnames(data.funcat)[1:2]<-c('Kegg.id','samp')

# Aggregating by Kegg categories
A.agg<-aggregate(value~Kegg.id+A+samp, data=data.funcat, FUN=unique)
A.agg<-aggregate(value~A+samp, data=A.agg, FUN=sum)
A.mat<-acast(samp~A, data=A.agg, fill=0)
A.mat.rel<-decostand(A.mat, 'total', 1) # Relative abundances

B.agg<-aggregate(value~Kegg.id+B+samp, data=data.funcat, FUN=unique)
B.agg<-aggregate(value~B+samp, data=B.agg, FUN=sum)
B.mat<-acast(samp~B, data=B.agg, fill=0)
B.mat.rel<-decostand(B.mat, 'total', 1) # Relative abundances

C.agg<-aggregate(value~Kegg.id+C+samp, data=data.funcat, FUN=unique)                             #Node~C+samp, data=funcat[-which(funcat$Kegg.id%in%c('K03781','K00249', 'K02298','K00356','K00341','K00336','K00382','K02160','K11263','K01615','K00208','K00123','K00257','K01209','K01869','K01958','K01995','K01998','K11754','K02890','K08348','K03781','K06044','K11263', 'K06861')),], FUN=flen) # 2 functions overly abundant in a few samples : K03781, K00249
C.agg<-aggregate(value~C+samp, data=C.agg, FUN=sum)
C.mat<-acast(samp~C, data=C.agg, fill=0)
# Remove traits with little overall abundance if wanted, e.g. C.mat<-C.mat[,which(colSums(C.mat)>=100)]
C.mat.rel<-decostand(C.mat, 'total', 1) # Relative abundances: USE!

# FINAL OBJECTS
# Tax4Fun
tax4.catA<-A.mat.rel
tax4.catB<-B.mat.rel
tax4.catC<-C.mat.rel


###### Rarefied datasets

### From dataset rarefied from the pool of all sequences with taxonomic annotations.
rareseqtax <- read.csv("/data/shared/panama_metagenome/code_data/rareseqtax.csv", header=T, row.names=1)
### From dataset rarefied from the pool of sequences with both taxonomic AND functional annotation
# Use rareseqfun
rareseqfun <- read.csv("/data/shared/panama_metagenome/code_data/rareseqfun.csv", header=T, row.names=1)

## TODO update this section to be reproducible - load data into two separate objects

## TODO for now we use rareseqtax as the basis for meta.tax and tax objects
meta.tax <- rareseqtax # either rareseqtax or rareseqfun
# Keep taxonomic annotations
meta.tax<-meta.tax[,which(colnames(meta.tax)%in%c('samp.x.x','V1.1','V2','V3','V4','V5','V6','V7'))]
# Make list to hold data tables of taxon abundances at different taxonomic levels
# (index 1 for phylum, to 6 for species)
tax<-list()
lev<-paste('V',c(2:7), sep='')
for (i in 1:length(lev)){
  prep<-aggregate(rownames(meta.tax)~samp.x.x+meta.tax[,i+2], data=meta.tax, FUN=function(x) length(unique(x)))
  colnames(prep)<-c('A','B','C')
  ee<-acast(prep, A~B, value.var='C')
  ee[is.na(ee)]<-0
  ee<-decostand(ee, method='total')
  tax[[i]]<-ee
}
meta.tax <- tax
