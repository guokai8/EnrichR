## ------------------------------------------------------------------------
set.seed(1234)   
library(EnrichR)   
# To check if your the current species if supported !!!
showData()   
showensemble()  
showplant()
# Make the GO and KEGG Pathway data for your analysis
# find suitable species name by using showensemble()    
hsa_go<-makeGOdat(species="human",keytype="SYMBOL")
hsa_ko<-makeKOdat(species = "human",keytype="SYMBOL",builtin = F)
# find suitable species name supported by reactome by using showAvailableRO()
# hsa_ro<-makeROdata(species = "Homo_sapiens") 

## ------------------------------------------------------------------------
# rice_go<-makeplantann(species="Oryza sativa Japonica",ann_type = "GO")   #check the species name by using showplant()
# rice_ko<-makeplantann(species="Oryza sativa Japonica",ann_type = "KEGG") 
# rice_pfam<-makeplantann(species="Oryza sativa Japonica",ann_type = "PFAM")
# rice_inter<-makeplantann(species="Oryza sativa Japonica",ann_type = "InterPro")
# rice_ro<-makePlantROdat(species = "Oryza_sativa") #check the species name by using showAvailablePlants()
## MSU version GO and KEGG infromation also supported named ricego,riceko   
## Zea may V2 GO and KEGG annotation data also supported named zm_v2_go and zm_v2_ko
# we also collect Reactome database for plant, you can just use makePlantROdat function to get RO data.   

## ----fig.height=6,fig.width=6,fig.align="center",dpi=300-----------------
df<-data.frame(gene=sample(unique(hsa_go$SYMBOL),2000),padj=abs(rnorm(2000,0,0.01)))
rownames(df)<-df$gene
res<-GE(df,GO_FILE = hsa_go,gene.cutoff = 0.01)
head(res)

## ----fig.height=6,fig.width=6,fig.align="center"-------------------------
GE.plot(resultFis =res,top=20,usePadj=F,pvalue.cutoff=0.05)
resk<-KE(df,KO_FILE = hsa_ko,gene.cutoff = 0.05,builtin = F)
head(resk)
KE.plot(resultFis = resk,top=10,pvalue.cutoff = 0.05)

## ----fig.height=6,fig.width=6,fig.align="center",dpi=300-----------------
richplot(res,usePadj=F,top=20)
## ----fig.height=6,fig.width=6,fig.align="center",dpi=300-----------------
###df could also be a vector of the genes you used for enrichment analysis
netmap(df=df,rhs=res,top=20,pvalue.cutoff = 0.05,weightcut = 0.01,visNet = T,nodeselect=T)
## ----fig.height=6,fig.width=6,fig.align="center",dpi=300-----------------
gnet(df=df,rhs=res,top=20,pvalue.cutoff = 0.05,weightcut = 0.01,vertex.label.cex=4)
## ----fig.height=6,fig.width=6,fig.align="center",dpi=300-----------------
mnetmap(df=df,gores=res[1:30,],kores=resk,pvalue.cutoff = 0.05,top=50)
## ------------------------------------------------------------------------
resgo<-getdetail(res,df)
head(resgo,6)
resko<-getdetail(resk,df)
head(resko,6)
