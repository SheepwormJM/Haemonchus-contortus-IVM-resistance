getwd()
setwd("C:/Users/jenni/Documents/Working from home Mar 2020/BUG/From Ali/")

# L1 = a, b, c
# ShL3 = d,e,f

a<-read.table("Galaxy62-[L1_ERS209879_count].tabular", header=F)
b<-read.table("Galaxy64-[L1_ERS209880_count].tabular", header=F)
c<-read.table("Galaxy66-[L1_ERS209881_count].tabular", header=F)

d<-read.table("Galaxy21-[Sheathed_L3_ERS092633_count].tabular", header = F)
e<-read.table("Galaxy22-[Sheathed_L3_ERS092634_count].tabular", header = F)
f<-read.table("Galaxy23-[Sheathed_L3_ERS092635_count].tabular", header = F)

summary(a)
summary(b)
summary(c)
summary(d)
summary(e)
summary(f)

# All loaded ok

# To merge L1 data and then L3 data by gene name:

tmp<-merge(a,b, by="V1")
L1<-merge(tmp, c, by="V1")
names(L1)<-c("Gene", "L1_ERS209879_count", "L1_ERS209880_count", "L1_ERS209881_count")
summary(L1)
write.table(L1, "L1_featureCount_WBPS14.txt", sep="\t", quote=FALSE)

tmp<-merge(d,e, by="V1")
L3<-merge(tmp, f, by="V1")
names(L3)<-c("Gene", "Sheathed_L3_ERS092633_count", "Sheathed_L3_ERS092634_count", "Sheathed_L3_ERS092635_count")
summary(L3)
write.table(L3, "L3_featureCount_WBPS14.txt", sep="\t", quote=FALSE)

# To merge RNAi genes and L1/L3 counts:

getwd()
setwd("C:/Users/jenni/Documents/Working from home Mar 2020/BUG/From Ali/")

a<-read.table("L1_featureCount_WBPS14.txt", sep="\t", header = T)
b<-read.table("L3_featureCount_WBPS14.txt", sep="\t", header = T)

summary(a)
summary(b)

tmp<-merge(a,b, by="Gene")
summary(tmp)

# Load list of Hcon homologues of Cele RNAi genes:

setwd("C:/Users/jenni/Documents/Working from home Mar 2020/BUG/RNAi (excl. qPCRs)/")

c<-read.table("RNAi genes Cele Hcon.txt", sep="\t", header = T)
summary(c)

RNAi<-merge(c,tmp, by="Gene")
summary(RNAi)

# Great. Notice that one Hcon gene is in 3 times, and another in twice. 
write.table(RNAi, "RNAi Cele genes Hcon homologues feature counts from Ali.txt", sep="\t", quote=FALSE)

library(ggplot2)
library(reshape2)
library(viridis)
new_RNAi<-melt(RNAi)
## Using Gene, Cele_gene, Present_in_Hcon_protein.fa. as id variables
x<-ggplot(data=new_RNAi)+
geom_col(aes(x=Cele_gene, y=value, fill=variable), position = "dodge")+
theme (axis.text.x = element_text (angle = 45))+
scale_fill_viridis(direction=1, option="E", alpha=0.6, discrete=T)+
labs(x = "C. elegans gene name", y = "H. contortus RNA-Seq Feature Count", caption = "Data source: Feature count from Ali", tag = "A", fill = "Lifecycle stage", subtitle="*ego-1, rrf-1 and rrf-3 all have same Hcon homologue\n*alg-1 and alg-2 have same Hcon homologue")
x


# To get the Feature counts of the RNA-Seq data as aligned to the WBPS15 version of Hcon.
# Their metadata (find here: )
# ERR229531	ERR229531	L1, inbred ISE	L1		inbred ISE
ERR229532	ERR229532	L1, inbred ISE	L1		inbred ISE
ERR229533	ERR229533	L1, inbred ISE	L1

# Search on ENA with the ERR codes - provides secondary sample accessions of: ERS209879, ERS209880 and ERS209881
# These correspond to Ali's previous L1. 

#ERR229546	ERR229546	sheathed L3, inbred ISE	sheathed L3		Inbred ISE
ERR229547	ERR229547	sheathed L3, inbred ISE	sheathed L3		Inbred ISE
ERR229548	ERR229548	sheathed L3, inbred ISE	sheathed L3		Inbred ISE
# SEarch on ENA with the ERR codes - provides secondary sample accesssions of: ERS092633, ERS092634 and  ERS092635	
# These correspond to Ali's previous shL3.

# Then get from WBPS15 the 'Counts of aligned reads per run (FeatureCounts)'

# Will want to extract just the columns you want in R:

a<-read.table("ERP002173.counts_per_run.tsv", sep="\t", header = T)
summary(a)

# Extract just the L1 and L3 data: 
cols_to_select = c("gene_id", "ERR229531", "ERR229532", "ERR229533", "ERR229546", "ERR229547", "ERR229548") # Select several columns using a variable
new_a<-a[,cols_to_select]
summary(new_a)
names(new_a) = c("Gene","L1_ERR229531", "L1_ERR229532", "L1_ERR229533", "shL3_ERR229546", "shL3_ERR229547", "shL3_ERR229548") # supply a new list of names
head(new_a)
# Perfect

c<-read.table("RNAi genes Cele Hcon.txt", sep="\t", header = T)
summary(c)

RNAi<-merge(c,new_a, by="Gene")
summary(RNAi)

# Great. Notice that one Hcon gene is in 3 times, and another in twice. 
write.table(RNAi, "RNAi Cele genes Hcon homologues feature counts from WBPS15.txt", sep="\t", quote=FALSE)

library(ggplot2)
library(reshape2)
library(viridis)
new_RNAi<-melt(RNAi)
## Using Gene, Cele_gene_name, Present_in_Hcon_protein.fa as id variables
x<-ggplot(data=new_RNAi)+
geom_col(aes(x=Cele_gene_name, y=value, fill=variable), position = "dodge")+
theme (axis.text.x = element_text (angle = 45))+
scale_fill_viridis(direction=1, option="E", discrete=T)+
labs(x = "C. elegans gene name", y = "H. contortus RNA-Seq Feature Count", caption = "Data source: Feature count from WBPS15", tag = "A", fill = "Lifecycle stage", subtitle="*ego-1, rrf-1 and rrf-3 all have same Hcon homologue\n*alg-1 and alg-2 have same Hcon homologue")
x

#theme(axis.text.x=element_text(face="italic"))+


