library(PharmacoGx)
library(Biobase)
library(scales)
library(piano)
library(GSA)



source("helpers.R")


# retrieve expression profiles of breast cancer cell lines
rnaseq_BC <- readRDS("Data/BC_rnaSeq.rda")
rnaseq_BC <- t(exprs(rnaseq_BC))

# retrieve area-above-the-curve (AACs) represeting the response of cell lines to UNC0642
AACs <- readRDS("Data/UNC0642_AAC.rda")
AACs <- rescale(AACs,to = c(0,1),from = c(0,25)) 

# mappings between different gene IDs such as ENSEMBLEID, Symbol, ENTREZ and so on
geneMapping <- readRDS("Data/geneMapping.rda")



#########################################
#########################################

# Finding the level of associations between each gene expression level and AAC. 
# This will serve as the ranking list fed to GSEA

drug_assoc <- PharmacoGx:::rankGeneDrugSensitivity(data = rnaseq_BC,drugpheno = AACs,single.type = T,nthread = 4)
drug_assoc_drug <- as.matrix(drug_assoc[[1]])
class(drug_assoc_drug) <- "numeric"

drug_assoc_drug <- as.data.frame(drug_assoc_drug,stringsAsFactors=F)
drug_assoc_drug <- drug_assoc_drug[
  with(drug_assoc_drug, order(drug_assoc_drug[,"fdr"], -1*drug_assoc_drug[,"estimate"])),
  ]



# Loading the MYC regulation signatures as gene sets
gSets <- GSA.read.gmt("Data/MYC_regulation_signatures.gmt")
dfgSets <- as.data.frame(cbind(unlist(gSets$genesets))) ## genes
dfgSNames <- as.data.frame(cbind(unlist(gSets$geneset.names)))
listS <- lapply(gSets$genesets,length)
gTogs_Final <- data.frame(dfgSets$V1,rep(dfgSNames$V1,listS))
names(gTogs_Final) <- c("V1","V2")
gTogs_Final <- gTogs_Final[gTogs_Final$V1!="",]

ibx <- match(gTogs_Final$V1,geneMapping[,"Symbol"])
idx <- which(is.na(ibx))
gTogs_Final <- gTogs_Final[-idx,]
ibx <- match(gTogs_Final$V1,geneMapping[,"Symbol"])
gTogs_Final$V1 <- rownames(geneMapping)[ibx]

gsc1 <- loadGSC(gTogs_Final)

drug_assoc_drug_NA <- drug_assoc_drug
drug_assoc_drug_NA[which(is.na(drug_assoc_drug_NA$pvalue)),c("estimate","tstat","pvalue")] <- c(0,0,1)
drug_assoc_drug_NA$fdr <- p.adjust(drug_assoc_drug_NA$pvalue,method = "fdr")
genelevelstats <- drug_assoc_drug_NA[,"tstat"]

names(genelevelstats) <- rownames(drug_assoc_drug_NA)
genelevelstats <- sort(genelevelstats,decreasing = T)

## running GSEA on 10 different MYC regulations signatures using rakning of the genes based on their univariate associations with UNC0642 response.
gsea_out <- piano::runGSA(geneLevelStats=genelevelstats, geneSetStat="fgsea", gsc=gsc1,  nPerm=1000000, ncpus=4, adjMethod="none", verbose=FALSE)
#saveRDS(gsea_out,"gsea_out.rda")

gseares <- try(piano::GSAsummaryTable(gsea_out))

gseares$pval <- apply(gseares[,c("p (dist.dir.up)"),drop=F],1,min,na.rm=T)
gseares$fdr <- p.adjust(gseares$pval,method = "fdr")

gseares <- gseares[,c("Name","Stat (dist.dir)","Genes (tot)","Genes (up)","Genes (down)" ,"pval","fdr")]

gseares <- gseares[order(gseares$pval),]



# The top four MYC signatures associated with UNC0642 response
sigPathways <- c("MYC_UP.V1","HALLMARK_MYC_TARGETS_V1","DANG_BOUND_BY_MYC","FERNANDEZ_BOUND_BY_MYC")


library(fgsea)

pdf("Plots/MYC_UP.V1_PathwayEnrichment.pdf",height = 5,width = 7)
x <- "MYC_UP.V1"
plotEnrichment(gsc1[[1]][[x]],
               genelevelstats) + labs(title=x,  subtitle=paste("Enrichment: " ,sprintf("%.3g", gseares[which(gseares$Name==x),"Stat (dist.dir)"]), ", FDR: ",sprintf("%.1E",gseares[which(gseares$Name==x),"fdr"]),sep = "" )) +
  theme( panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,face = "bold"), plot.subtitle  = element_text(hjust = 0.5,face = "bold") 
         ,axis.title=element_text(face="bold"),axis.text=element_text(face="bold") )
dev.off()
pdf("Plots/PathwayEnrichmentAnalysis.pdf",height = 15,width = 20)

plots <- lapply(gseares$Name, function(x){
  p <- plotEnrichment(gsc1[[1]][[x]],
                      genelevelstats) + labs(title=x,  subtitle=paste("Enrichment: " ,sprintf("%.3g", gseares[which(gseares$Name==x),"Stat (dist.dir)"]), ", FDR: ",sprintf("%.1E",gseares[which(gseares$Name==x),"fdr"]),sep = "" )) +
    theme( panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,face = "bold"), plot.subtitle  = element_text(hjust = 0.5,face = "bold") 
           ,axis.title=element_text(face="bold"),axis.text=element_text(face="bold") )
  return(p)
})

multiplot(plots[[1]],plots[[5]],plots[[8]],
          plots[[3]],plots[[6]],plots[[7]],
          plots[[4]],plots[[9]],plots[[10]],cols=3)


dev.off()


pdf("Plots/metagene_MYC_UP.V1.pdf",width = 7,height = 7)
lapply(sigPathways[1],function(x){
  
  genesForPathway <- gTogs_Final[which(gTogs_Final[,2]==x),1]
  leadingEdgeGenes <- genelevelstats[genesForPathway]
  leadingEdgeGenes <- sort(leadingEdgeGenes,decreasing = T)
  
  #1 Pearson-Mean
  testEnrichPA <- lapply(1:length(leadingEdgeGenes),function(x){
    metaGene <- apply(rnaseq_BC[,names(leadingEdgeGenes)[1:x],drop=F],1,mean)
    A <- cor.test(metaGene,AACs[names(metaGene)],method = "p")
    return(c("cor"=A$estimate,"pval"=A$p.value))
  })
  
  corrsPA <- unlist(lapply(testEnrichPA,"[[",1))
  pvalsPA <- unlist(lapply(testEnrichPA,"[[",2))
  fdrPA <- p.adjust(pvalsPA, method = "fdr")

  l <- which(corrsPA==max(corrsPA))
  ChosenGenes <- names(leadingEdgeGenes)[1:l]
  names(ChosenGenes) <- geneMapping[ChosenGenes,"Symbol"]
  
  
  
  metaGene <- apply(rnaseq_BC[,ChosenGenes,drop=F],1,mean)
  A <- cor.test(metaGene,AACs[names(metaGene)],method = "p")
  results <- c("cor"=A$estimate,"pval"=A$p.value)
  results
  
  plot(metaGene,AACs[names(metaGene)],main = paste("Metagene based on ",x," signature [top ",l," leading genes]\n and its association with G9A inhibitor")
       ,ylab = "AAC",xlab = "Metagene expression level",pch=16)

  reg1 <- glm(formula =AACs[names(metaGene)]~metaGene)
  abline(reg1,col="blue")
  
  legendLabel <- c(paste("Pearson cor.: " ,sprintf("%.3g", corrsPA[l]),sep = ""), paste("FDR: ",sprintf("%.2E",fdrPA[l]),sep = "" )) 
  
  if(corrsPA[l]>0){
    legend("topleft",legend = legendLabel,bty = "n")
  }else{
    legend("topright",legend = legendLabel,bty = "n")
  }
})
  
dev.off()





pdf("Plots/metagene_supp.pdf",width = 21,height = 7)
par(mfrow=c(1,3))
lapply(sigPathways[2:4],function(x){
  
  genesForPathway <- gTogs_Final[which(gTogs_Final[,2]==x),1]
  leadingEdgeGenes <- genelevelstats[genesForPathway]
  leadingEdgeGenes <- sort(leadingEdgeGenes,decreasing = T)
  
  #1 Pearson-Mean
  testEnrichPA <- lapply(1:length(leadingEdgeGenes),function(x){
    metaGene <- apply(rnaseq_BC[,names(leadingEdgeGenes)[1:x],drop=F],1,mean)
    A <- cor.test(metaGene,AACs[names(metaGene)],method = "p")
    return(c("cor"=A$estimate,"pval"=A$p.value))
  })
  
  corrsPA <- unlist(lapply(testEnrichPA,"[[",1))
  pvalsPA <- unlist(lapply(testEnrichPA,"[[",2))
  fdrPA <- p.adjust(pvalsPA, method = "fdr")
  
  l <- which(corrsPA==max(corrsPA))
  ChosenGenes <- names(leadingEdgeGenes)[1:l]
  names(ChosenGenes) <- geneMapping[ChosenGenes,"Symbol"]
  
  
  
  metaGene <- apply(rnaseq_BC[,ChosenGenes,drop=F],1,mean)
  A <- cor.test(metaGene,AACs[names(metaGene)],method = "p")
  results <- c("cor"=A$estimate,"pval"=A$p.value)
  results
  
  plot(metaGene,AACs[names(metaGene)],main = paste("Metagene based on ",x," signature [top ",l," leading genes]\n and its association with G9A inhibitor")
       ,ylab = "AAC",xlab = "Metagene expression level",pch=16)
  
  #  text(x = metaGene,y = AACs[names(metaGene)],labels = names(metaGene))
  reg1 <- glm(formula =AACs[names(metaGene)]~metaGene)
  abline(reg1,col="blue")
  
  
  legendLabel <- c(paste("Pearson cor.: " ,sprintf("%.3g", corrsPA[l]),sep = ""), paste("FDR: ",sprintf("%.2E",fdrPA[l]),sep = "" )) 
  
  if(corrsPA[l]>0){
    legend("topleft",legend = legendLabel,bty = "n")
  }else{
    legend("topright",legend = legendLabel,bty = "n")
  }
})

dev.off()






