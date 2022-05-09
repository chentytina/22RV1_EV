### ++++++++++++++++++++++++++++++++++++
### load library & function
### ++++++++++++++++++++++++++++++++++++

libs<-c("tidyr","dplyr","ggplot2","ggrepel","RColorBrewer","ggpubr","ggbeeswarm",
        "tibble","rmarkdown","pheatmap","colorspace","colormap","reshape2","plyr",
        "data.table","rpart","scales","fitdistrplus","plyranges","tidyverse","viridis",
        "doParallel","factoextra","org.Hs.eg.db","GO.db","GSEABase",
        "limma","edgeR","AnnotationDbi","enrichR","GSVA","annotables")

lapply(libs, require, character.only = TRUE)
#registerDoParallel(makeCluster(6)) # Repeat if cores get reasigned and only if varPar is used.
rm(libs)
#functions
my_filter_list_by_row <- function(x){
  ii <- NULL
  for (i in 1:length(x)) { ii[[i]]<- nrow(x[[i]]) > 0}
  ii <- unlist(ii)
  return(x[ii])
}

### ++++++++++++++++++++++++++++++++++++
### RNASEQ : 22RV1
### ++++++++++++++++++++++++++++++++++++
###read raw counts/phenotype
fcounts <- read.csv("22RV1_RNAseq_raw_counts.csv", row.names = 1)
phenotype <- read.csv("22RV1_RNAseq_metadata.csv")

###filter counts 
x <- fcounts
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
keep.rows <- rowSums(cpm>10)>=1

###Create a DGEList object
dcpmx <- DGEList(counts=fcounts[keep.rows,],group=phenotype$condition, genes = my_labels[keep.rows,])
dcpmx <- calcNormFactors(dcpmx, method = "TMM") #method TMM
#model
design = model.matrix( ~ 0 + condition + concentration + libprep, data=phenotype)
#dispersion
dcpmx = estimateCommonDisp(dcpmx, verbose=TRUE)
dcpmx = estimateTagwiseDisp(dcpmx)

### Voom normalization
v_m3n1 <- voom(dcpmx, design, plot=TRUE, save.plot = TRUE)
fit_imodel <- lmFit(v_m3n1, design)
efit_imodel <- eBayes(fit_imodel)
summary(decideTests(efit_imodel))
contr.matrix <- makeContrasts(
  "10B_vs_19A" = condition10B - condition19A, # EV vs Cell
  levels = colnames(design))
fit_imodel <- lmFit(v_m3n1, design)
efit_imodel <- contrasts.fit(fit_imodel, contrasts=contr.matrix)
efit_imodel <- eBayes(efit_imodel)
summary(decideTests(efit_imodel))

### Accounting for replicate
```{r echo=TRUE, message=FALSE, warning=FALSE}
#adjutment per replicate
phenotype3 <- phenotype[phenotype$full_sample %in% colnames(v_m3n1),]
#block containing the replicates 
block <- as.numeric(as.factor(paste(phenotype3$cc2,phenotype3$rep,sep="")))
#function estimates the correlation between repeated observations
dupcor <- duplicateCorrelation(v_m3n1, design, block=block)
#model accounting for replicates
fit_imodel <- lmFit(v_m3n1, design, block=block, correlation=dupcor$consensus)
efit_imodel <- contrasts.fit(fit_imodel, contrasts=contr.matrix)
efit_imodel <- eBayes(efit_imodel)
summary(decideTests(efit_imodel))
#obtain the results
results_DEGs_10B_vs_19A_3ng <- topTable(efit_imodel,n=Inf,coef=1)
colnames(results_DEGs_10B_vs_19A_3ng)[1:6]<- c("geneid","chr","start","end","strand","length")
results_DEGs_10B_vs_19A_3ng$ensembl <- as.character(gsub("\\.[0-9]*$", "", results_DEGs_10B_vs_19A_3ng$geneid))
results_DEGs_10B_vs_19A_3ng$geneSymbols <- mapIds(org.Hs.eg.db, keys = as.character(gsub("\\.[0-9]*$", "", results_DEGs_10B_vs_19A_3ng$geneid)), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
results_DEGs_10B_vs_19A_3ng$entrez <- mapIds(org.Hs.eg.db, keys = as.character(gsub("\\.[0-9]*$", "", results_DEGs_10B_vs_19A_3ng$geneid)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")

### ++++++++++++++++++++++++++++++++++++
### FIGURE 4B: PCA
### ++++++++++++++++++++++++++++++++++++
#get principal component
mydata <- prcomp(t(v_m3n1$E) , scale=TRUE ) 

#scree plot
Scree.plot <- fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw() + ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=rel(4.8)))
eigen <- get_eigenvalue(mydata)
pca1 <- format(round(eigen$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(eigen$variance.percent[2], 2), nsmall = 2)
Scree.plot

#plot PCA
pca_mts<-as.data.frame(mydata$x)
pca_mts$full_sample <- rownames(pca_mts)
pca_mts <- merge(pca_mts, phenotype, by="full_sample")
pca_mts$lib.size <- dcpmx$samples$lib.size
pca_mts$norm.factors <- dcpmx$samples$norm.factors

PCA <- ggplot(pca_mts, aes(PC1, PC2, color = condition)) + #,shape=libprep)) + 
  #geom_point(aes(size=lib.size)) + geom_label_repel(aes(label=condition),show.legend = FALSE) +
  geom_point(size = 3) + theme_bw() + ggtitle('') + 
  labs(shape = "Library\nPreparation") + guides(color = F) +
  scale_color_manual(values=c("royalblue4", "red3")) +
  xlab("PC1 42.5%") + ylab("PC2 7.5%") +      
  #geom_label_repel(aes(label = cc)) +
  annotate(geom="text", x=50, y=-25, label='atop(bold("EVs"))',color="royalblue4",size=10, parse = T) +
  annotate(geom="text", x=-50, y=-25, label='atop(bold("Cell"))',color="red3",size=10, parse = T) +
  theme(text = element_text(size=rel(4.8)))
PCA

### ++++++++++++++++++++++++++++++++++++
### FIGURE 4C: VOLCANO PLOT OF DEG
### ++++++++++++++++++++++++++++++++++++
#identify FDR < 0.05 genes
sig.volcano <- results_DEGs_10B_vs_19A_3ng[results_DEGs_10B_vs_19A_3ng$adj.P.Val<0.05,]
pos_sig.genes <-sig.volcano[sig.volcano$logFC>5,]
pos.sig.genes <- head(arrange(pos_sig.genes, adj.P.Val, desc(logFC)),5)
neg_sig.genes <- sig.volcano[sig.volcano$logFC< -5,]
neg.sig.genes <- head(arrange(neg_sig.genes,adj.P.Val, logFC),5)
i <- read.csv("i.csv")
pos.sig.genes$geneSymbols <- i$geneSymbols[1:5]
neg.sig.genes$geneSymbols <- i$geneSymbols[11:15]

#plot volcano plot
volcano.plot<- ggplot(data=results_DEGs_10B_vs_19A_3ng,aes(x=logFC,y=-log10(adj.P.Val))) + 
  geom_point(data=results_DEGs_10B_vs_19A_3ng[results_DEGs_10B_vs_19A_3ng$adj.P.Val>0.05,],color="slategray",alpha=0.7, size = 3) + 
  geom_point(data=sig.volcano[sig.volcano$logFC > 0 ,],color="lightskyblue2",size = 3) +      
  geom_point(data=sig.volcano[sig.volcano$logFC < 0 ,],color="rosybrown2",size = 3) + 
  geom_point(data= pos.sig.genes, color = "royalblue4",size = 3) +
  geom_point(data= neg.sig.genes, color = "red3",size = 3) +
  geom_label_repel(data=pos.sig.genes, aes(label = geneSymbols), max.overlaps = 50, color = "royalblue4", size = 5) +
  geom_label_repel(data=neg.sig.genes, aes(label = geneSymbols), max.overlaps = 50, color = "red3", size = 5) +
  geom_vline(xintercept = 0, linetype = 2, color = 'black') +
  geom_hline(yintercept = 1.3, linetype = 2, color = 'black') +
  xlab("Log(Fold Change)") + ylab("-Log(Adjusted P Value)") +
  annotate(geom="text", x=-9, y=2, label='atop(bold("Cells"))',color="red3",size=10, parse = T) +
  annotate(geom="text", x=9, y=2, label='atop(bold("EVs"))',color="royalblue4",size=10, parse = T) +
  theme_bw() + #ggtitle("Volcano Plot: Exosome vs Cell") + 
  theme(text = element_text(size=rel(4.8))) 
volcano.plot

### ++++++++++++++++++++++++++++++++++++
### FIGURE 4D: GSVA HEATMAP FOR HALLMARK GENES
### ++++++++++++++++++++++++++++++++++++
#Perform GSVA   
#download hallmark genes
gene.set <- getGmt("/Users/chentytina/Documents/Tina Chen/ICAHN MSSM/Dogra/22RV1/Analysis/h.all.v6.2.symbols.gmt")
gene.matrix <- as.matrix(v_m3n1$E)
rownames(gene.matrix) <- v_m3n1$genes$symbol
GSVA_matrix <- gsva(gene.matrix, gene.set,method = "ssgsea", 
    kcdf = "Poisson", # Gaussian for microarrays, Poisson for RNAseq 
    verbose=TRUE, ssgsea.norm = TRUE)

GSVA_hallmark_EVs_vs_Cells_score_matrix <- GSVA_matrix

GSVA_score <- GSVA_hallmark_EVs_vs_Cells_score_matrix
colnames(GSVA_score) <- phenotype3$condition

scoring_stat <- matrix(ncol=2,nrow=nrow(GSVA_score))
for (i in 1:nrow(GSVA_score)){
scoring_stat[i,1] <- rownames(GSVA_score)[i]
scoring_stat[i,2] <- wilcox.test(GSVA_score[i,phenotype3$condition %in% "10B"], GSVA_score[i,phenotype3$condition %in% "19A"])$p.value
}
colnames(scoring_stat) <- c("Pathway","p.value")
scoring_stat <- as.data.frame(scoring_stat)
scoring_stat$p.value <- as.numeric(as.character(scoring_stat$p.value))
scoring_stat$FDR <- p.adjust(scoring_stat$p.value,method="BH") 
scoring_stat <- scoring_stat[scoring_stat$FDR<0.05,]
scoring_stat <- scoring_stat[order(scoring_stat$FDR), ]
ii <- which(rownames(GSVA_score) %in% scoring_stat[,1])

generic_var <- data.frame(
  Normal = rowMedians(GSVA_score[ii,colnames(GSVA_score) %in% "10B"]),
  Resistant = rowMedians(GSVA_score[ii,colnames(GSVA_score) %in% "19A"])
)

rownames(generic_var) <- rownames(GSVA_score)[ii]
generic_var$Difference <- abs(generic_var$Normal - generic_var$Resistant)
generic_var <- generic_var[order(generic_var$Difference), ]
generic_var <- generic_var[(generic_var$Difference) > 0.1, ]

col <- data.frame(Condition =  phenotype3$condition  )
colnames(GSVA_score) <- paste( 1:12, colnames(GSVA_score),sep="_")
rownames(col) <- colnames(GSVA_score)

col$Condition[ col$Condition %in% "10B"] <- "Exosomes"
col$Condition[ col$Condition %in% "19A"] <- "Cells"

ii <- which(rownames(GSVA_score) %in% rownames(generic_var))

rownames(GSVA_score) <- gsub("HALLMARK_","",rownames(GSVA_score))

annotation_colors <- list (Condition = c( Exosomes="royalblue4", Cells="red3") )

pheatmap(GSVA_score[ii,], color = colorRampPalette(c("green4", "yellow","red"))(255),  
         cluster_cols = T, cluster_rows=T, 
         show_colnames = F,
         show_rownames = T, 
         main="ssGSVA Scoring", # treeheight_row = 0,  treeheight_col = 0, #hides trees
         annotation_col = col,
         annotation_colors = annotation_colors,
         border_color="black",
         clustering_distance_rows ="euclidean", scale="none")


### ++++++++++++++++++++++++++++++++++++
### FIGURE 4E - 4G: SPEARMAN CORRELATIONS
### ++++++++++++++++++++++++++++++++++++
#labels for top genes for each condition
sig.results_rnaseq_22RV1 <- results_DEGs_10B_vs_19A_3ng[results_DEGs_10B_vs_19A_3ng$adj.P.Val<0.05,]
sig_50_logFC <- rbind(head(arrange(sig.results_rnaseq_22RV1,desc(logFC)),50), 
                arrange(head(arrange(sig.results_rnaseq_22RV1,logFC),50), desc(logFC)))
sig_50_logFC$status <- ifelse(sig_50_logFC$logFC > 0, "EVs", "Cells")
sig_5_logFC <- rbind(head(arrange(sig_50_logFC,desc(logFC)),5), 
                     arrange(head(arrange(sig_50_logFC,logFC),5), desc(logFC)))

sig_5_list <- sig_5_logFC[,c(7,15,17)]

generic_var <- data.frame(v1=v_m3n1$E[,4], v2=v_m3n1$E[,3])
generic_var$ensembl_gene_id <- rownames(generic_var)
generic_var <- merge(generic_var, sig_5_list, by = "ensembl_gene_id", all = T)
exo.exo <- ggplot(generic_var, aes(x = v1, y = v2)) + geom_point() +
  #geom_point(data = generic_var[generic_var$ensembl_gene_id %in% 
  #    sig_5_list$ensembl_gene_id[sig_5_list$status == "EVs"],], color="royalblue1",alpha=0.7, size = 2) +
  xlab("log2(EVs CPM)") + ylab("Log2(EVs CPM)")+ theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=rel(6))) +  
  stat_cor(data = generic_var, method = "spearman", label.x = 1, label.y = 15, size = 10) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 2) + 
  stat_smooth(method = 'lm', color = 'purple') +
  stat_density_2d(show.legend = F)
  #geom_label_repel(data=generic_var[generic_var$ensembl_gene_id %in% 
  #    sig_5_list$ensembl_gene_id[sig_5_list$status == "EVs"],],
  #    aes(label = geneSymbols), color = "royalblue1", size = 5, ylim = c(NA, 7.5), xlim = c(7.5, NA))
exo.exo

generic_var <- data.frame(v1=v_m3n1$E[,5], v2=v_m3n1$E[,10])
#rho <- paste("Rho = ", format(as.numeric(cor.test(x = generic_var$v2, y = generic_var$v1, method = "spearman", exact = F)$estimate), digits = 3))
cell30.cell30 <- ggplot(generic_var, aes(x = v1, y = v2)) + geom_point() +
  xlab("Log2(Cells CPM)") + ylab("Log2(Cells CPM)")+ theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=rel(6))) +  
  stat_cor(data = generic_var, method = "spearman", label.x = 1, label.y = 15, size = 10, r.digits = 2,
  #aes(label = paste("Rho = ",..r.label.., sep = " "))
  ) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 2) + 
  stat_smooth(method = 'lm', color = 'purple') +
  stat_density_2d(show.legend = F)
cell30.cell30

generic_var <- data.frame(v1=v_m3n1$E[,3], v2=v_m3n1$E[,10])
generic_var$ensembl_gene_id <- rownames(generic_var)
generic_var <- merge(generic_var, sig_25_list, by = "ensembl_gene_id", all = T)
cell30.exo <- ggplot(generic_var, aes(x = v1, y = v2)) + geom_point() +
  geom_point(data = head(generic_var[generic_var$ensembl_gene_id %in% 
                     sig_25_list$ensembl_gene_id & generic_var$v2 < 1 & generic_var$v2 > 0,],5), color="royalblue1", size = 2) +
  geom_point(data = head(generic_var[generic_var$ensembl_gene_id %in% 
                     sig_25_list$ensembl_gene_id & generic_var$v1 < 1 & generic_var$v1 > 0,],5), color="red", size = 2) +
  ylab("Log2(Cells CPM)") + xlab("Log2(EVs CPM)") + theme_bw() + ylim(0,17.5) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=rel(6))) +  
  stat_cor(data = generic_var, method = "spearman", label.x = 1, label.y = 17, size = 10) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 2) + 
  stat_smooth(method = 'lm', color = 'purple') +
  stat_density_2d(show.legend = F) +
  geom_label_repel(data = head(generic_var[generic_var$ensembl_gene_id %in% 
                     sig_25_list$ensembl_gene_id & generic_var$v2 < 1 & generic_var$v2 > 0,],5),
      aes(label = geneSymbols), color = "royalblue1", size = 7, ylim = c(NA, 7.5), xlim = c(7.5, NA)) + 
  geom_label_repel(data = generic_var[generic_var$ensembl_gene_id %in% 
                     sig_25_list$ensembl_gene_id & generic_var$v1 < 1 & generic_var$v1 > 0,],
      aes(label = geneSymbols), color = "red", size = 7,ylim = c(7.5,NA), xlim = c(1, NA))
cell30.exo

### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### RNASEQ : External Dataset (Hinger et al., Herrara et al.)
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### EXTERNAL DATASET : Herrara et al.
MOESM3_raw <- read.csv("MOESM3.csv")
#make phenotype
MOESM3_phenotype <- data.frame(full_sample = colnames(MOESM3_raw)[-1])
MOESM3_phenotype$condition <- ifelse(grepl("EXO", MOESM3_phenotype$full_sample), "EVs", "Cells")
MOESM3_phenotype$cell_line <- ifelse(grepl("NF", MOESM3_phenotype$full_sample), "NF", "CAF")
MOESM3_phenotype$sample_id <- paste(MOESM3_phenotype$cell_line, substr(MOESM3_phenotype$full_sample, 2, 2), sep = "_")
MOESM3_phenotype$cell_line_type <- paste(MOESM3_phenotype$cell_line, MOESM3_phenotype$condition, sep = "_")

#extract gene symbol/biotype/ensembl id
library(EnsDb.Hsapiens.v75)
MOESM3_raw$transcript <- substring(MOESM3_raw$NAME, 7,21)
MOESM3_raw$biotype <- substr(MOESM3_raw$NAME, 23, nchar(MOESM3_raw$NAME))
MOESM3_raw$symbol <- mapIds(EnsDb.Hsapiens.v75, keys = MOESM3_raw$transcript, column="SYMBOL", keytype="TXID", multiVals="first")
MOESM3_raw$ensembl_gene_id <- mapIds(EnsDb.Hsapiens.v75, keys = MOESM3_raw$transcript, column="GENEID", keytype="TXID", multiVals="first")

#filtering low expression genes
MOESM3_data <- MOESM3_raw[,colnames(MOESM3_raw) %in% MOESM3_phenotype$full_sample]
rownames(MOESM3_data) <- MOESM3_raw$transcript
table(rowSums(cpm(MOESM3_data)>0)>=1) 
keep <- rowSums(cpm(MOESM3_data)>0)>=1
dcpmx <- DGEList(counts=MOESM3_data[keep, ], group=MOESM3_phenotype$condition, genes = rownames(MOESM3_data))
dcpmx <- calcNormFactors(dcpmx, method = "TMM") #method TMM

identical(MOESM3_phenotype$full_sample, colnames(dcpmx))

#dream normalization
library(BiocParallel)
library(variancePartition)
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)

form <- ~ 0 + cell_line_type + (1|sample_id)

MOESM3_dream = voomWithDreamWeights(dcpmx, form, MOESM3_phenotype)

L =  makeContrastsDream(form, MOESM3_phenotype, 
                         contrasts = c("cell_line_typeNF_Cells - cell_line_typeNF_EVs", 
                                       "cell_line_typeCAF_Cells - cell_line_typeCAF_EVs"))

fitmm = dream(MOESM3_dream, form, MOESM3_phenotype, L)

#saving result
lmfreq_results <- list() ; for ( i in 1:2 ) { lmfreq_results[[i]] <- topTable(fitmm, coef=i,n=Inf,adjust.method="fdr") 
lmfreq_results[[i]]$Comparison <- colnames(fitmm$coefficients)[i]
lmfreq_results[[i]]$Marker <- rownames(topTable(fitmm, coef=i,n=Inf,adjust.method="fdr") ) }
external_dataset_result <- do.call(rbind,lmfreq_results)
rownames(external_dataset_result) <- NULL
external_dataset_result$nLogFDR <- -log10(external_dataset_result$adj.P.Val)
external_dataset_result$Comparison <- gsub('cell_line_type','',external_dataset_result$Comparison)
results_MOESM3 <- external_dataset_result
results_MOESM3$symbol <- MOESM3_raw$symbol[match(results_MOESM3$genes, MOESM3_raw$transcript)]
results_MOESM3$ensembl_gene_id <- MOESM3_raw$ensembl_gene_id[match(results_MOESM3$genes, MOESM3_raw$transcript)]

#organizing dream object
MOESM3_dream$genes$symbol <- MOESM3_raw$symbol
MOESM3_dream$genes$ensembl_gene_id <- MOESM3_raw$ensembl_gene_id



### EXTERNAL DATASET : Hinger et al.
Hinger_raw <- read.csv("Hinger et al RNA-Seq Read Counts.csv")
Hinger_raw <- Hinger_raw[,-grep("X",colnames(Hinger_raw))]
#make phenotype
Hinger_phenotype <- data.frame(full_sample = colnames(Hinger_raw)[grep("cell|exo",colnames(Hinger_raw))])
Hinger_phenotype$condition <- ifelse(grepl("exo", Hinger_phenotype$full_sample), "EVs", "Cells")
Hinger_phenotype$cell_line <- substr(Hinger_phenotype$full_sample, 1, 5)
Hinger_phenotype$sample_id <- paste(Hinger_phenotype$cell_line, 
		substr(Hinger_phenotype$full_sample, nchar(Hinger_phenotype$full_sample), 
			nchar(Hinger_phenotype$full_sample)), sep = "_")
Hinger_phenotype$cell_line_type <- paste(Hinger_phenotype$cell_line, Hinger_phenotype$condition, sep = "_")

#filtering genes with no expression
Hinger_data <- Hinger_raw[,colnames(Hinger_raw) %in% Hinger_phenotype$full_sample]
rownames(Hinger_data) <- Hinger_raw$gene
table(rowSums(cpm(Hinger_data)>0)>=1)
keep <- rowSums(cpm(Hinger_data)>0)>=1
dcpmx <- DGEList(counts=Hinger_data[keep, ], group=Hinger_phenotype$condition, genes = rownames(Hinger_data)[keep])
dcpmx <- calcNormFactors(dcpmx, method = "TMM") #method TMM

identical(Hinger_phenotype$full_sample, colnames(dcpmx))

#dream normalization
form <- ~ 0 + cell_line_type + (1|sample_id)

Hinger_dream = voomWithDreamWeights(dcpmx, form, Hinger_phenotype)

L =  makeContrastsDream(form, Hinger_phenotype, 
                         contrasts = c("cell_line_typeDKs.8_Cells - cell_line_typeDKs.8_EVs",
                         			   "cell_line_typeDKO.1_Cells - cell_line_typeDKO.1_EVs", 
                                       "cell_line_typeDLD.1_Cells - cell_line_typeDLD.1_EVs"
                                       ))

fitmm = dream(Hinger_dream, form, Hinger_phenotype, L)

#saving result
Hinger_lmfreq_results <- list() ; for ( i in 1:3 ) { Hinger_lmfreq_results[[i]] <- topTable(fitmm, coef=i,n=Inf,adjust.method="fdr") 
Hinger_lmfreq_results[[i]]$Comparison <- colnames(fitmm$coefficients)[i]
Hinger_lmfreq_results[[i]]$Marker <- rownames(topTable(fitmm, coef=i,n=Inf,adjust.method="fdr") ) }
Hinger_external_dataset_result <- do.call(rbind,Hinger_lmfreq_results)
rownames(Hinger_external_dataset_result) <- NULL
Hinger_external_dataset_result$nLogFDR <- -log10(Hinger_external_dataset_result$adj.P.Val)
Hinger_external_dataset_result$Comparison <- gsub('cell_line_type','',Hinger_external_dataset_result$Comparison)
results_Hinger <- Hinger_external_dataset_result
results_Hinger$symbol <- Hinger_raw$gene.symbol[match(results_Hinger$genes, Hinger_raw$gene)]

#organizing dream object
Hinger_dream$genes$symbol <- Hinger_raw$gene.symbol[match(Hinger_dream$genes$genes, Hinger_raw$gene)]
Hinger_dream$genes$ensembl_gene_id <- as.character(gsub("\\.[0-9]*$", "", Hinger_dream$genes$genes))

### Combine counts
raw_dcount <- raw.voom$E
rownames(raw_dcount) <- raw.voom$genes$ensembl_gene_id
MOESM3_dcount <- MOESM3_dream$E[!is.na(MOESM3_dream$genes$ensembl_gene_id),]
rownames(MOESM3_dcount) <- MOESM3_dream$genes$ensembl_gene_id[!is.na(MOESM3_dream$genes$ensembl_gene_id)]
Hinger_dcount <- Hinger_dream$E[!is.na(Hinger_dream$genes$ensembl_gene_id),]
rownames(Hinger_dcount) <- Hinger_dream$genes$ensembl_gene_id[!is.na(Hinger_dream$genes$ensembl_gene_id)]

dcount <- merge(raw_dcount, Hinger_dcount, by = 0, all.x = T)
rownames(dcount) <- dcount$Row.names
dcount <- dcount[,-1]

dcount <- merge(dcount, MOESM3_dcount, by = 0, all.x = T)
rownames(dcount) <- dcount$Row.names
dcount <- dcount[,-1]
dcount[is.na(dcount)] = 0  
dcount$symbol <- raw.voom$genes$symbol[match(rownames(dcount), raw.voom$genes$ensembl_gene_id)]

#combine phenotype
dcount_phenotype <- data.frame(full_sample = colnames(dcount)[-ncol(dcount)])
dcount_phenotype$condition[dcount_phenotype$full_sample %in% 
				ds22RV1_phenotype$name] <- ds22RV1_phenotype$condition
dcount_phenotype$condition[dcount_phenotype$full_sample %in% 
				Hinger_phenotype$full_sample] <- Hinger_phenotype$condition
dcount_phenotype$condition[dcount_phenotype$full_sample %in% 
				MOESM3_phenotype$full_sample] <- MOESM3_phenotype$condition
dcount_phenotype$condition <- gsub("19A", "Cells", dcount_phenotype$condition)
dcount_phenotype$condition <- gsub("10B", "EVs", dcount_phenotype$condition)
dcount_phenotype$dataset[dcount_phenotype$full_sample %in% 
				ds22RV1_phenotype$name] <- "22RV1"
dcount_phenotype$dataset[dcount_phenotype$full_sample %in% 
				Hinger_phenotype$full_sample] <- "Hinger et al"
dcount_phenotype$dataset[dcount_phenotype$full_sample %in% 
				MOESM3_phenotype$full_sample] <- "MOESM3"
dcount_phenotype$cell_line[dcount_phenotype$dataset == "22RV1"] <- "22RV1"
dcount_phenotype$cell_line[dcount_phenotype$dataset == "Hinger et al"] = 
				substr(dcount_phenotype$full_sample[dcount_phenotype$dataset == "Hinger et al"], 1, 5)
dcount_phenotype$cell_line[dcount_phenotype$dataset == "MOESM3"] <- 
	ifelse(dcount_phenotype$condition[dcount_phenotype$dataset == "MOESM3"] == "Cells", 
		substr(dcount_phenotype$full_sample[dcount_phenotype$dataset == "MOESM3"], 4,
			nchar(dcount_phenotype$full_sample[dcount_phenotype$dataset == "MOESM3"])-5),
		substr(dcount_phenotype$full_sample[dcount_phenotype$dataset == "MOESM3"], 4,
			nchar(dcount_phenotype$full_sample[dcount_phenotype$dataset == "MOESM3"])-4))


### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### FIGURE 4H: SPEARMAN CORRELATIONS OF ALL NORMALIZED COUNTS
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#find optimal samples for scatter plots via spearman correlation
corr_matrix <- cor(dcount, method = "spearman")
scatter_list <- NULL
opt <- melt(corr_matrix)
for (i in 1:length(unique(dcount_phenotype$cell_line))) {
    cell.line <- unique(dcount_phenotype$cell_line)[i]
    exo <- dcount_phenotype$full_sample[dcount_phenotype$cell_line %in% cell.line & 
    		dcount_phenotype$condition == "EVs"]
    cell <- dcount_phenotype$full_sample[dcount_phenotype$cell_line %in% cell.line & 
    		dcount_phenotype$condition == "Cells"]
    exo_df <- opt[opt$Var1 %in% exo,]
    exo_df <- exo_df[exo_df$Var2 %in% cell,]
    exo_df <- exo_df[which.min(exo_df$value),]
    scatter_list <- rbind(scatter_list, exo_df)
}
colnames(scatter_list) <- c("x", "y", "corr")
scatter_list$cell_line <- dcount_phenotype$cell_line[dcount_phenotype$full_sample 
								%in% scatter_list$x]

#normalize counts across datasets
norm_count <- dcount[,grep(paste(c(scatter_list$x, scatter_list$y), collapse = "|"), colnames(dcount))]
norm_count$symbol <- dcount$symbol
norm_count$ensembl <- rownames(dcount)
lst <- c("10B|19A", "DKO", "DLD", "DKs", "NF", "CAF")

df_list <- NULL
for (i in 1:length(lst)) {
    df <- norm_count[,grep(lst[i], colnames(norm_count))]
    min.x <- min(df[,grep("EVs|exo|EXO",colnames(df))])
    min.y <- min(df[,grep("Cell|CELL|cell",colnames(df))])
    max.x <- max(df[,grep("EVs|exo|EXO",colnames(df))])
    max.y <- max(df[,grep("Cell|CELL|cell",colnames(df))])    
    diff.x <- max.x - min.x
    diff.y <- max.y - min.y

    df_list <- rbind(df_list, 
		    data.frame(y = (df[,grep("Cell|CELL|cell",colnames(df))] - min.y)/ diff.y, 
		    	x = (df[,grep("EVs|exo|EXO",colnames(df))] - min.x)/ diff.x, 
		    	cell_line = scatter_list$cell_line[i], symbol = norm_count$symbol, 
		    	ensembl = norm_count$ensembl)
    	)
}

#colors for plotting
col <- c("22RV1" = "springgreen4", "CAF" = "steelblue4", "DKO.1" = "goldenrod3", 
	"DKs.8" = "firebrick4", "DLD.1" = "mediumpurple3", "NF" = "lightcoral")

#plotting 
ggplot(df_list) + theme_bw() + 
  	#ylim(min(df$y)-2, max(df$y)+2) + xlim(min(df$x)-4, max(df$x)+4) + 
    xlab(paste0("Log2(EVs CPM)")) + 
    ylab(paste0("Log2(Cells CPM)")) + 
    geom_point(data = df_list[df_list$cell_line == "DKs.8",], 
    	aes(x = x, y = y, color = cell_line), size = 2, alpha = 0.5) +    
    geom_point(data = df_list[df_list$cell_line == "22RV1",], 
    	aes(x = x, y = y, color = cell_line), size = 2, alpha = 0.5) + 
    geom_point(data = df_list[df_list$cell_line == "DKO.1",], 
    	aes(x = x, y = y, color = cell_line), size = 2, alpha = 0.5) + 
    geom_point(data = df_list[df_list$cell_line == "DLD.1",], 
    	aes(x = x, y = y, color = cell_line), size = 2, alpha = 0.5) + 
    geom_point(data = df_list[df_list$cell_line == "CAF",], 
    	aes(x = x, y = y, color = cell_line), size = 2, alpha = 0.5) +
    geom_point(data = df_list[df_list$cell_line == "NF",], 
    	aes(x = x, y = y, color = cell_line), size = 2, alpha = 0.5) +    
    scale_color_manual(values=col) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    theme(plot.title = element_text(hjust = 0.5,size=rel(6)), 
    	text = element_text(size=rel(6)),legend.text=element_text(size=rel(6))) +  
	labs(color = "Cell Lines") + 
	guides(colour = guide_legend(override.aes = list(size=8)))



### ++++++++++++++++++++++++++++++++++++
### EV MASSPEC : 22RV1/22RV1DR
### ++++++++++++++++++++++++++++++++++++
mass_spec_results <- read.csv("EV_Mass.Spec_raw_counts.csv", row.names = 1)
#Separate values from main table
mass_spec_data <- mass_spec_results
mass_spec_genes <- data.frame(uniprot_id = rownames(mass_spec_results))
#add gene symbol
mass_spec_genes$gene_symbol <- as.character(mapIds(org.Hs.eg.db, keys = as.character(mass_spec_genes$uniprot_id), column="SYMBOL", keytype="UNIPROT", multiVals="first")) 

#identifying sensitive vs resistant signatures
mass_spec_results$resistant <-rowMeans(mass_spec_results[,grep("DR", colnames(mass_spec_results))])
mass_spec_results$sensitive <- rowMeans(mass_spec_results[,grep("RV1_", colnames(mass_spec_results))])

masspec_exo_sensitive <- data.frame(mass_spec_results[mass_spec_results$sensitive > 0, c("uniprot_id","gene_symbol","sensitive")])
colnames(masspec_exo_sensitive) <- c("uniprot_id", "gene_symbol", "counts")
masspec_exo_sensitive$label_color <- "normal" 
masspec_exo_resist <- data.frame(mass_spec_results[mass_spec_results$resistant > 0, c("uniprot_id","gene_symbol","resistant")])
colnames(masspec_exo_resist) <- c("uniprot_id", "gene_symbol", "counts")
masspec_exo_resist$label_color <- "resistant" 
masspec_exo <- rbind(masspec_exo_sensitive,masspec_exo_resist)
#rowSums(mass_spec_data>0)



### ++++++++++++++++++++++++++++++++++++
### CELL MASSPEC : 22RV1/22RV1DR
### ++++++++++++++++++++++++++++++++++++
masspec_lysate <- read.csv("Cell_Mass.Spec_raw_counts.csv", row.names = 1)
masspec_lysate$uniprot_id <- rownames(masspec_lysate)
masspec_lysate[,grep("X", colnames(masspec_lysate))] <- 10^masspec_lysate[,grep("X", colnames(masspec_lysate))]
masspec_lysate[, grep("X", colnames(masspec_lysate))][is.na(masspec_lysate[, grep("X", colnames(masspec_lysate))])] <- 0
#add gene symbol
masspec_lysate$gene_symbol <- as.character(mapIds(org.Hs.eg.db, keys = as.character(masspec_lysate$uniprot_id), column="SYMBOL", keytype="UNIPROT", multiVals="first")) 

#identifying sensitive vs resistant signatures
masspec_lysate$resistant <-rowMeans(masspec_lysate[,grep("DR", colnames(masspec_lysate))])
masspec_lysate$sensitive <- rowMeans(masspec_lysate[,grep("RV1_", colnames(masspec_lysate))])

masspec_lysate_sensitive <- data.frame(masspec_lysate[masspec_lysate$sensitive > 0, c("uniprot_id","gene_symbol","sensitive")])
colnames(masspec_lysate_sensitive) <- c("uniprot_id", "gene_symbol", "counts")
masspec_lysate_sensitive$label_color <- "normal" 
masspec_lysate_resist <- data.frame(masspec_lysate[masspec_lysate$resistant > 0, c("uniprot_id","gene_symbol","resistant")])
colnames(masspec_lysate_resist) <- c("uniprot_id", "gene_symbol", "counts")
masspec_lysate_resist$label_color <- "resistant" 
masspec_lysate_all <- rbind(masspec_lysate_sensitive,masspec_lysate_resist)



### ++++++++++++++++++++++++++++++++++++
### MASSPEC : External Dataset (Hurwitz et al.)
### ++++++++++++++++++++++++++++++++++++
#load cell data
breast_cell_dataset <- read.csv("breast_cell_dataset.csv")
colnames(breast_cell_dataset) <- ifelse(breast_cell_dataset[2,] == "", breast_cell_dataset[1,], breast_cell_dataset[2,])
colnames(breast_cell_dataset)[5:34] <- paste(colnames(breast_cell_dataset)[5:34], 
                    str_sub(breast_cell_dataset[1,5:34], start = -1,end = -1), sep = "_")
breast_cell_dataset <- breast_cell_dataset[-c(1:2),]
#extract gene names
breast_cell_dataset$`Gene names` <- sub(';.*$','', breast_cell_dataset$`Gene names`)
breast_cell_dataset <- breast_cell_dataset[!duplicated(breast_cell_dataset$`Gene names`),]
#create count/protein name dataframe
breast_cell_data <- breast_cell_dataset[,grep("1|2|3",colnames(breast_cell_dataset))]
rownames(breast_cell_data) <- breast_cell_dataset$`Gene names`
breast_cell_protein <- breast_cell_dataset[,grep("protein|Protein|names",colnames(breast_cell_dataset))]

#load ev data
breast_ev_dataset <- read.csv("breast_ev_dataset.csv")
colnames(breast_ev_dataset) <- ifelse(breast_ev_dataset[2,] == "", breast_ev_dataset[1,], breast_ev_dataset[2,])
colnames(breast_ev_dataset)[5:44] <- paste(colnames(breast_ev_dataset)[5:44], 
                    str_sub(breast_ev_dataset[1,5:44], start = -1,end = -1), sep = "_")
breast_ev_dataset <- breast_ev_dataset[-c(1:3),]
#extract gene names
breast_ev_dataset$`Gene names` <- sub(';.*$','', breast_ev_dataset$`Gene names`)
breast_ev_dataset <- breast_ev_dataset[!duplicated(breast_ev_dataset$`Gene names`),]
#create count/protein name dataframe
breast_ev_data <- breast_ev_dataset[,grep("1|2|3|4",colnames(breast_ev_dataset))]
rownames(breast_ev_data) <- breast_ev_dataset$`Gene names`
breast_ev_protein <- breast_ev_dataset[,1:4]

#find overlapping cell lines between cell and ev
breast_ev_data <- breast_ev_data[,colnames(breast_ev_data) %in% colnames(breast_cell_data)]
breast_cell_data <- breast_cell_data[,colnames(breast_cell_data) %in% colnames(breast_ev_data)]
colnames(breast_cell_data) <- paste(colnames(breast_cell_data), "Cell", sep = "_")
colnames(breast_ev_data) <- paste(colnames(breast_ev_data), "EV", sep = "_")

#create phenotype
phenotype_breast <- data.frame(samples = colnames(breast_cell_data), type = "Cell")
phenotype_breast <- rbind(phenotype_breast, data.frame(samples = colnames(breast_ev_data), type = "EV"))
phenotype_breast$cell_line <-sub("_.*", "", phenotype_breast$samples)
phenotype_breast$cell_line_type <- paste(phenotype_breast$cell_line, phenotype_breast$type, sep = "_")

#merge cell/ev dataframe
i <- breast_ev_data[rownames(breast_ev_data) %in% rownames(breast_cell_data),]
i$genes <- rownames(i)
ii <- breast_cell_data[rownames(breast_cell_data) %in% rownames(breast_ev_data),]
ii$genes <- rownames(ii)
breast_merged_data <- merge(i, ii, by = "genes")
breast_merged_data <- breast_merged_data[-c(1:8),]
rownames(breast_merged_data) <- breast_merged_data$genes
breast_merged_data <- mutate_all(breast_merged_data[,-1], function(x) as.numeric(as.character(x)))
#keep overlapping cell lines/ proteins for phenotype
rownames(phenotype_breast) = phenotype_breast$samples
phenotype_breast <- phenotype_breast[match(colnames(breast_merged_data), rownames(phenotype_breast)),]
phenotype_breast$sample_id <- phenotype_breast$samples
phenotype_breast$sample_id <-  gsub("_EV","",phenotype_breast$sample_id)
phenotype_breast$sample_id <-  gsub("_Cell","",phenotype_breast$sample_id)

#dream noramlization to find cell vs ev for overlapping cell lines
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)

form <- ~ 0 + cell_line_type + (1|sample_id)
L =  makeContrastsDream(form, phenotype_breast, 
                         contrasts = c("cell_line_typeBT549_Cell - cell_line_typeBT549_EV", 
                                       "cell_line_typeHCC1419_Cell - cell_line_typeHCC1419_EV",
                                       "cell_line_typeHCC1954_Cell - cell_line_typeHCC1954_EV",
                                       "cell_line_typeHs578T_Cell - cell_line_typeHs578T_EV",
                                       "cell_line_typeJIMT1_Cell - cell_line_typeJIMT1_EV",
                                       "cell_line_typeLM2_Cell - cell_line_typeLM2_EV",
                                       "cell_line_typeMCF10A_Cell - cell_line_typeMCF10A_EV",
                                       "cell_line_typeMCF7_Cell - cell_line_typeMCF7_EV",
                                       "cell_line_typeMDAMB231_Cell - cell_line_typeMDAMB231_EV",
                                       "cell_line_typeSKBR3_Cell - cell_line_typeSKBR3_EV"
                                       ))
fitmm = dream( breast_merged_data, form, phenotype_breast, L)
#save result
lmfreq_results <- list() ; for ( i in 1:10 ) { lmfreq_results[[i]] <- topTable(fitmm, coef=i,n=Inf,adjust.method="fdr") 
lmfreq_results[[i]]$Comparison <- colnames(fitmm$coefficients)[i]
lmfreq_results[[i]]$Marker <- rownames( topTable(fitmm, coef=i,n=Inf,adjust.method="fdr") ) }
external_dataset_result <- do.call(rbind,lmfreq_results)
rownames(external_dataset_result) <- NULL
external_dataset_result$nLogFDR <- -log10(external_dataset_result$adj.P.Val)
external_dataset_result$Comparison <- gsub('cell_line_type','',external_dataset_result$Comparison)
sig <- external_dataset_result[external_dataset_result$adj.P.Val < 0.05, ]
table(sig$Comparison)
sig_cell <- data.frame(sig %>% group_by(Comparison) %>% arrange(desc(logFC)) %>% slice(1:100))
sig_ev <- data.frame(sig %>% group_by(Comparison) %>% arrange(logFC) %>% slice(1:100))


### ++++++++++++++++++++++++++++++++++++
### FIGURE 3A : MASSPEC PIE CHART
### ++++++++++++++++++++++++++++++++++++
#select enrichR database
dbs <- listEnrichrDbs()
go_dbs <- dbs$libraryName[grep(paste("GO_Cellular_Component_2021",collapse="|"),dbs$libraryName)]

### Resistant EV pathways
resist.exo_enrichr <- NULL
generic_var_x <- NULL

#create genes list for enrichR input
gene_list <- list(as.character(masspec_exo_resist$gene_symbol))
my_genes<-NULL
for (i in 1:length(gene_list)) {
  my_genes <- gene_list[[i]]  
  generic_var_x <- enrichr( genes=my_genes, databases = go_dbs)
  generic_var_x <- my_filter_list_by_row(generic_var_x)
  for (x in 1:length(generic_var_x)) {  generic_var_x[[x]]$database <- names(generic_var_x)[x] }
  resist.exo_enrichr[[i]] <- do.call(rbind,generic_var_x) } 
#identify relatively significant pathways & create log scale score
sig_resist.exo_enrichr<-NULL
for( i in 1:length(resist.exo_enrichr)) { 
  sig_resist.exo_enrichr[[i]] <-  resist.exo_enrichr[[i]][ resist.exo_enrichr[[i]]$Adjusted.P.value < 1,] }
names(sig_resist.exo_enrichr) <- names(gene_list)
resist.exo_enrichr_masspec <- sig_resist.exo_enrichr[[1]]
resist.exo_enrichr_masspec$Score<- -log10(resist.exo_enrichr_masspec$Adjusted.P.value)
resist.exo_enrichr_masspec$input.val <- as.numeric(sub("\\/.*", "", resist.exo_enrichr_masspec$Overlap))
resist.exo_enrichr_masspec$Status <- "Resistant_EV"



### Sensitive EV pathways
normal.exo_enrichr <- NULL
generic_var_x <- NULL

#create genes list for enrichR input
gene_list <- list( as.character(masspec_exo_sensitive$gene_symbol))
for (i in 1:length(gene_list)) {
  my_genes <- gene_list[[i]]  
  generic_var_x <- enrichr( genes=my_genes, databases = go_dbs)
  generic_var_x <- my_filter_list_by_row(generic_var_x)
  for (x in 1:length(generic_var_x)) { generic_var_x[[x]]$database <- names(generic_var_x)[x] }
  normal.exo_enrichr[[i]] <- do.call(rbind,generic_var_x) } 
#identify relatively significant pathways & create log scale score
sig_normal.exo_enrichr<-NULL
for( i in 1:length(normal.exo_enrichr)) {
  sig_normal.exo_enrichr[[i]] <-  normal.exo_enrichr[[i]][ normal.exo_enrichr[[i]]$Adjusted.P.value < 1,] }
names(sig_normal.exo_enrichr) <- names(gene_list)
normal.exo_enrichr_masspec <- sig_normal.exo_enrichr[[1]]
normal.exo_enrichr_masspec$Score<- -(-log10(normal.exo_enrichr_masspec$Adjusted.P.value))
normal.exo_enrichr_masspec$input.val <- as.numeric(sub("\\/.*", "", normal.exo_enrichr_masspec$Overlap))
normal.exo_enrichr_masspec$Status <- "Normal_EV"



### Resistant Cell pathways
resist.cell_enrichr <- NULL
generic_var_x <- NULL

#create genes list for enrichR input
gene_list <- list(as.character(masspec_lysate_resist$gene_symbol))
my_genes<-NULL
for (i in 1:length(gene_list)) {
  my_genes <- gene_list[[i]]  
  generic_var_x <- enrichr( genes=my_genes, databases = go_dbs)
  generic_var_x <- my_filter_list_by_row(generic_var_x)
  for (x in 1:length(generic_var_x)) {  generic_var_x[[x]]$database <- names(generic_var_x)[x] }
  resist.cell_enrichr[[i]] <- do.call(rbind,generic_var_x) } 
#identify relatively significant pathways & create log scale score
sig_resist.cell_enrichr<-NULL
for( i in 1:length(resist.cell_enrichr)) { sig_resist.cell_enrichr[[i]] <-  
  resist.cell_enrichr[[i]][ resist.cell_enrichr[[i]]$Adjusted.P.value < 1,] }
names(sig_resist.cell_enrichr) <- names(gene_list)
resist.cell_enrichr_masspec <- sig_resist.cell_enrichr[[1]]
resist.cell_enrichr_masspec$Score<- -log10(resist.cell_enrichr_masspec$Adjusted.P.value)
resist.cell_enrichr_masspec$input.val <- as.numeric(sub("\\/.*", "", resist.cell_enrichr_masspec$Overlap))
resist.cell_enrichr_masspec$Status <- "Resistant_cell"



### Sensitive Cell pathways
normal.cell_enrichr <- NULL
generic_var_x <- NULL

#create genes list for enrichR input
gene_list <- list( as.character(masspec_lysate_sensitive$gene_symbol) )
for (i in 1:length(gene_list)) {
  my_genes <- gene_list[[i]]  
  generic_var_x <- enrichr( genes=my_genes, databases = go_dbs)
  generic_var_x <- my_filter_list_by_row(generic_var_x)
  for (x in 1:length(generic_var_x)) { generic_var_x[[x]]$database <- names(generic_var_x)[x] }
  normal.cell_enrichr[[i]] <- do.call(rbind,generic_var_x) } 
#identify relatively significant pathways & create log scale score
sig_normal.cell_enrichr<-NULL
for( i in 1:length(normal.cell_enrichr)) { sig_normal.cell_enrichr[[i]] <-  
	normal.cell_enrichr[[i]][ normal.cell_enrichr[[i]]$Adjusted.P.value < 1,] }
names(sig_normal.cell_enrichr) <- names(gene_list)
normal.cell_enrichr_masspec <- sig_normal.cell_enrichr[[1]]
normal.cell_enrichr_masspec$Score<- -(-log10(normal.cell_enrichr_masspec$Adjusted.P.value))
normal.cell_enrichr_masspec$input.val <- as.numeric(sub("\\/.*", "", normal.cell_enrichr_masspec$Overlap))
normal.cell_enrichr_masspec$Status <- "Normal_cell"



### Sensitive Cell Pathways
normal.cell_enrichr <- NULL
generic_var_x <- NULL

#create genes list for enrichR input
gene_list <- list( as.character(masspec_lysate_sensitive$gene_symbol) )
for (i in 1:length(gene_list)) {
  my_genes <- gene_list[[i]]  
  generic_var_x <- enrichr( genes=my_genes, databases = go_dbs)
  generic_var_x <- my_filter_list_by_row(generic_var_x)
  for (x in 1:length(generic_var_x)) { generic_var_x[[x]]$database <- names(generic_var_x)[x] }
  normal.cell_enrichr[[i]] <- do.call(rbind,generic_var_x) } 
#identify relatively significant pathways & create log scale score
sig_normal.cell_enrichr<-NULL
for( i in 1:length(normal.cell_enrichr)) { sig_normal.cell_enrichr[[i]] <-  
	normal.cell_enrichr[[i]][ normal.cell_enrichr[[i]]$Adjusted.P.value < 1,] }
names(sig_normal.cell_enrichr) <- names(gene_list)
normal.cell_enrichr_masspec <- sig_normal.cell_enrichr[[1]]
normal.cell_enrichr_masspec$Score<- -(-log10(normal.cell_enrichr_masspec$Adjusted.P.value))
normal.cell_enrichr_masspec$input.val <- as.numeric(sub("\\/.*", "", normal.cell_enrichr_masspec$Overlap))
normal.cell_enrichr_masspec$Status <- "Normal_cell"

#merge all pathway files
cc_all.list <- rbind(normal.exo_enrichr_masspec,resist.exo_enrichr_masspec,resist.cell_enrichr_masspec,normal.cell_enrichr_masspec)
#separate GO number
cc_all.list$GO.num <- gsub("\\(([^()]*)\\)|.", "\\1", x = cc_all.list$Term, perl = T)



### External Dataset Pathways; pathway files extracted via online EnrichR (https://maayanlab.cloud/Enrichr)
#load text files containing enrichR result for cellular component 2021
df_list = list.files(path = "breast_enrichr",pattern = ".*.txt", full.names = T)
datalist = lapply(df_list, FUN=read.delim, header=TRUE)
#extract cell line and condition (cell/ev) and add status column
Status <- substring(df_list, 16, nchar(df_list)-4)
for (i in 1:length(datalist)) {
    datalist[[i]]$Status <- Status[i]
}
#combine all pathway files
cc.breast = do.call("rbind", datalist) 
#calculate log score
cc.breast$Score<- -log10(cc.breast$Adjusted.P.value)
cc.breast$input.val <- as.numeric(sub("\\/.*", "", cc.breast$Overlap))
cc.breast$GO.num <- gsub("\\(([^()]*)\\)|.", "\\1", x = cc.breast$Term, perl = T)



### Classification for different cellular componenets; referencing QuickGO (https://www.ebi.ac.uk/QuickGO/annotations)
nucleus <- c("nuclear","nucleolus", "Cdc73", "mRNA", "U12", "U2", "spliceosomal complex", 
  "nucleolar", "histone","Cajal","0001650", "cyclin-dependent", "THO", "MLL", "cataly", "GO:0071564",
  "cyclin/CDK", "ESC", "prespliceosome", "chromosome", "chromosomal", "chromatin", "Sin3", "GO:0071565",
  "DNA","meiotic", "nucleopl", "anaphase","PML","polymerase", "Swr","type complex", "0005634",
  "bodies","elong","MLH1","fork","PRC1", "TERF2","NURF","Mut","U2AF2","telo","TAF7","TAF9","GO:0016514",
  "Set1", "DMC1","TRRAP", "Ino80","MED8", "MED17","RBBP","SMARCC","TFIIH","PHF10","0048500", "GO:0033503",
  "0000152","0005721","strand", "0000243", "GO:0005669", "GO:0071162")
cytoplasm <- c("cytoplasm", "granule lumen", "clathrin-sculpted", "clathrin-coated vesicle", "deuterosome",
  "0030118", "0030136", "0030125", "0030130", "clathrin vesicle", "peroxi","0061645", "GO:0043231", "GO:0010369",
  "granule","vacuol", "microbody", "secretory","0005764","0005766", "0043292", "cortex", "GO:0043232", "GO:0002080",
  "insul", "coated vesicle", "sarcoplasm", "0030120", "COP", "cytosolic", "0140275", "GO:0097418", "GO:0001669",
  "chaperonin", "endosoma","lipid","mela","oid","lamellar","-body", "IkappaB","BLOC-1","0034709", "GO:0002080",
  "IPAF","NLRP3","AIM2","N-ter","oxog", "0034774", "phago","0060293", "0045009", "GO:0031095", "GO:0031094",
  "GO:0070013", "GO:0098588", "GO:0065010", "GO:0031970", "GO:0031968", "GO:0031301", "GO:0019866", "GO:0097708")
ribonucleoprotein <- c("preribosome","polysom", "snoRNP", "ribonuclease P", "GO:0031428",
  "micro-ribonucleoprotein", "snRNP","1990124", "RISC","RNAi")
ribosome <- c("0005840", "large ribosomal", "small ribosomal","unit")
pmembrane <- c("raft", "actin cortical", "ruffle mem", "axon", "channel", "apical","glycoprotein", "0016323", "GO:0005890", "GO:0090533", "GO:0033180",
  "flotillin", "recept", "of plasma membrane", "catenin","0043190","0033176","brush","caveola", "GO:0016011", 
  "0042383", "GO:0005834", "GO:0099060", "GO:0060170", "GO:0005922", "0032591", "GO:0098839", "GO:0099634", "GO:0016012",
   "MHC", "0032590","axolemma","attack","0089717","death", "gamma-secretase", "intrinsic", "exocytic", 
   "GO:0045259", "GO:0099056", "GO:0099699", "GO:0016327", 
   "0030425", "0036128", "0043020", "0001401", "0001533", "GO:0045271", "GO:0030132", "GO:0097025", "GO:0099061")
golgi <- c("Golgi")
cytoskeleton <- c("cytoskeleton", "myo", "cytoskeletal", "mitotic", "microtubule", "spindle", "0032059", "0030139","0030666", 
  "centrosome", "centri","filament", "0001725","kinesin", "HAUS", "gamma-tubulin", "0036157", "0002102")
endosome <- c("endosome", "endolysosom", "vesicular", "0045334", "0071682", "0030128", 
  "0030669", "adaptor", "0031982", "extracellular vesicle", "GO:0030666", "GO:0030658", "GO:0030285", "GO:0012506", "GO:0030672")
mitochondria <- c("mitochon", "MICOS")
lysosome <- c("lysosomal", "primary lysosome", "secondary", "PIK3C3;SQSTM1;FTL", "0005764", "0044754")
extracellular <- c("projection","filo", "podium", "villus", "bicellular","focal", "0044291", "0042627", "GO:0036056",
  "disc","zonu","node", "junction", "0062023", "cilium", "GO:0036057", "GO:0001527", "GO:0060076", "GO:0005587",
  "0030057", "0030056", "0036126", "GO:0098978", "GO:0032279", "GO:0014069", "GO:0005604", "GO:0097060")
ER <- c("endoplasmic", "sarcoplasmic", "junctional", "0071953", "GO:0005784","GO:0017059")
other <- c("Cul","ubiquitin","CUL1","CUL4A","CUL1","poprotein", "PRKRA", "serine/threonine", "GO:0043527",
  "aggre","FACIT","TFIIIC", "1902911","1902555", "GO:0099512", "GO:0072357", "GO:0000159", "GO:0061695",
  "GO:0031931", "GO:0031932", "GO:0098984", "0098985", "0060077", "GO:0098982")

### add corresponding categories
#for in house data
cc_all.list$Category <- NA
cc_all.list$Category[grep(paste(cytoskeleton,collapse="|"),cc_all.list$Term)] <- "cytoskeleton" 
cc_all.list$Category[with(cc_all.list, grepl(paste(nucleus,collapse="|"), 
  paste(Term, Genes)))] <- "nucleus"
cc_all.list$Category[grep(paste(cytoplasm,collapse="|"),cc_all.list$Term)] <- "cytoplasm"
cc_all.list$Category[grep(paste(ribonucleoprotein,collapse="|"),cc_all.list$Term)] <- "ribonucleoprotein"
cc_all.list$Category[grep(paste(ribosome,collapse="|"),cc_all.list$Term)] <- "ribosome"
cc_all.list$Category[grep(paste(pmembrane,collapse="|"),cc_all.list$Term)] <- "pmembrane"
cc_all.list$Category[grep(paste(ER,collapse="|"),cc_all.list$Term)] <- "ER"
cc_all.list$Category[grep(paste(golgi,collapse="|"),cc_all.list$Term)] <- "golgi"
cc_all.list$Category[grep(paste(endosome,collapse="|"),cc_all.list$Term)] <- "endosome"
cc_all.list$Category[grep(paste(mitochondria,collapse="|"),cc_all.list$Term)] <- "mitochondria"
cc_all.list$Category[with(cc_all.list, grepl(paste(lysosome,collapse="|"), 
  paste(Term, Genes)))] <- "lysosome"
cc_all.list$Category[grep(paste(extracellular,collapse="|"),cc_all.list$Term)] <- "extracellular"
cc_all.list$Category[with(cc_all.list, grepl(paste(other,collapse="|"), 
  paste(Term, Genes)))] <- "other"

#for external data
cc.breast$Category <- NA
cc.breast$Category[grep(paste(cytoskeleton,collapse="|"),cc.breast$Term)] <- "cytoskeleton" 
cc.breast$Category[with(cc.breast, grepl(paste(nucleus,collapse="|"), 
  paste(Term, Genes)))] <- "nucleus"
cc.breast$Category[grep(paste(cytoplasm,collapse="|"),cc.breast$Term)] <- "cytoplasm"
cc.breast$Category[grep(paste(ribonucleoprotein,collapse="|"),cc.breast$Term)] <- "ribonucleoprotein"
cc.breast$Category[grep(paste(ribosome,collapse="|"),cc.breast$Term)] <- "ribosome"
cc.breast$Category[grep(paste(pmembrane,collapse="|"),cc.breast$Term)] <- "pmembrane"
cc.breast$Category[grep(paste(ER,collapse="|"),cc.breast$Term)] <- "ER"
cc.breast$Category[grep(paste(golgi,collapse="|"),cc.breast$Term)] <- "golgi"
cc.breast$Category[grep(paste(endosome,collapse="|"),cc.breast$Term)] <- "endosome"
cc.breast$Category[grep(paste(mitochondria,collapse="|"),cc.breast$Term)] <- "mitochondria"
cc.breast$Category[with(cc.breast, grepl(paste(lysosome,collapse="|"), 
  paste(Term, Genes)))] <- "lysosome"
cc.breast$Category[grep(paste(extracellular,collapse="|"),cc.breast$Term)] <- "extracellular"
cc.breast$Category[with(cc.breast, grepl(paste(other,collapse="|"), 
  paste(Term, Genes)))] <- "other"



### Create percentage/position for pie chart
#22RV1DR ; resistant
cc.resist <- cc_all.list[grep("Resistant", cc_all.list$Status),]
exo.rsum = sum(cc.resist$input.val[cc.resist$Status == "Resistant_EV"])
cell.rsum = sum(cc.resist$input.val[cc.resist$Status == "Resistant_cell"])
cc.resist$percent <- ifelse(cc.resist$Status =="Resistant_EV", (cc.resist$input.val/exo.rsum)*100, 
  (cc.resist$input.val/cell.rsum)*100)
sum(cc.resist$percent[cc.resist$Status == "Resistant_cell"])
sum(cc.resist$percent[cc.resist$Status == "Resistant_EV"])
#create dataframe for plotting using combined percentages
cc.resist_plotting <- aggregate(percent~Category+Status,cc.resist,sum)
cc.resist_plotting <- data.frame(cc.resist_plotting %>% 
  arrange(Status, desc(percent)) %>% 
  group_by(Status) %>% 
  mutate(position = 100 - (cumsum(percent) - percent/2)))
#arrange order and position on geom col
cc.resist_plotting$faucet <- 1:26
cc.resist_plotting$faucet <- factor(cc.resist_plotting$faucet, levels = unique(cc.resist_plotting$faucet))
#rename categories
name.list <- data.frame(Names = unique(cc.resist_plotting$Category))
name.list$cap <- str_to_title(name.list$Names)
name.list$cap[grep("Golgi", name.list$cap)] <- "Golgi Apparatus"
name.list$cap[grep("Extracellular", name.list$cap)] <- "Extracellular Region or Secreted"
name.list$cap[grep("Pmembrane", name.list$cap)] <- "Plasma Membrane"
name.list$cap[grep("Er", name.list$cap)] <- "ER"
cc.resist_plotting$Names <- name.list$cap[match(cc.resist_plotting$Category, name.list$Names)]
#calculate average score to select top 5 categories
cc.resist_plotting$avg.score <- NA
for (i in 1:length(unique(cc.resist_plotting$Status))) {
    df <- cc.resist[cc.resist$Status == unique(cc.resist_plotting$Status)[i],]
    for (ii in 1:length(unique(cc.resist_plotting$Category))) {
      cc.resist_plotting$avg.score[cc.resist_plotting$Status == unique(df$Status) & 
        cc.resist_plotting$Category == unique(cc.resist_plotting$Category)[ii]] <- mean(df$Score[df$Category == unique(cc.resist_plotting$Category)[ii]])
    }
}

#22RV1 ; sensitive
cc.normal <- cc_all.list[grep("Normal", cc_all.list$Status),]
exo.nsum = sum(cc.normal$input.val[cc.normal$Status == "Normal_EV"])
cell.nsum = sum(cc.normal$input.val[cc.normal$Status == "Normal_cell"])
cc.normal$percent <- ifelse(cc.normal$Status =="Normal_EV", (cc.normal$input.val/exo.nsum)*100, 
  (cc.normal$input.val/cell.nsum)*100)
sum(cc.normal$percent[cc.normal$Status == "Normal_cell"])
sum(cc.normal$percent[cc.normal$Status == "Normal_EV"])
#create dataframe for plotting using combined percentages
cc.normal_plotting <- aggregate(percent~Category+Status,cc.normal,sum)
cc.normal_plotting <- data.frame(cc.normal_plotting %>% 
  arrange(Status, desc(percent)) %>% 
  group_by(Status) %>% 
  mutate(position = 100 - (cumsum(percent) - percent/2)))
#arrange order and position on geom col
cc.normal_plotting$faucet <- 1:26
cc.normal_plotting$faucet <- factor(cc.normal_plotting$faucet, levels = unique(cc.normal_plotting$faucet))
#rename categories
name.list <- data.frame(Names = unique(cc.normal_plotting$Category))
name.list$cap <- str_to_title(name.list$Names)
name.list$cap[grep("Golgi", name.list$cap)] <- "Golgi Apparatus"
name.list$cap[grep("Extracellular", name.list$cap)] <- "Extracellular Region or Secreted"
name.list$cap[grep("Pmembrane", name.list$cap)] <- "Plasma Membrane"
name.list$cap[grep("Er", name.list$cap)] <- "ER"
cc.normal_plotting$Names <- name.list$cap[match(cc.normal_plotting$Category, name.list$Names)]
#calculate average score to select top 5 categories
cc.normal_plotting$avg.score <- NA
for (i in 1:length(unique(cc.normal_plotting$Status))) {
    df <- cc.normal[cc.normal$Status == unique(cc.normal_plotting$Status)[i],]
    for (ii in 1:length(unique(cc.normal_plotting$Category))) {
      cc.normal_plotting$avg.score[cc.normal_plotting$Status == unique(df$Status) & 
        cc.normal_plotting$Category == unique(cc.normal_plotting$Category)[ii]] <- mean(df$Score[df$Category == unique(cc.normal_plotting$Category)[ii]])
    }
}
cc.normal_plotting$avg.score <- -(cc.normal_plotting$avg.score)

#External dataset
#create empty list
cc.breast_plotting <- list(); cc.breast.processed <- list()
#loop by different different cell lines
for (i in 1:length(unique(phenotype_breast$cell_line))) {
    name <- unique(phenotype_breast$cell_line)[i]
    df <- cc.breast[grep(name, cc.breast$Status),]
    exo.rsum = sum(df$input.val[grep("EV", df$Status)])
    cell.rsum = sum(df$input.val[grep("Cell", df$Status)])
    df$percent <- ifelse(grepl("EV", df$Status), (df$input.val/exo.rsum)*100, 
        (df$input.val/cell.rsum)*100)
    #save categorized pathways
    cc.breast.processed[[name]] <- df
    if (sum(df$percent[grep("EV", df$Status)]) || sum(df$percent[grep("Cell", df$Status)])) {
        cc.breast_plotting[[name]] <- data.frame(aggregate(percent~Category+Status,df,sum) %>% 
          arrange(Status, desc(percent)) %>% group_by(Status) %>% 
          mutate(position = 100 - (cumsum(percent) - percent/2)))
        #save individual plotting dataframe for different cell lines
        cc.breast_plotting[[name]]$faucet <- factor(1:nrow(cc.breast_plotting[[name]]), levels = 1:nrow(cc.breast_plotting[[name]]))
        #rename categories
        name.list <- data.frame(Names = unique(cc.breast_plotting[[name]]$Category))
        name.list$cap <- str_to_title(name.list$Names)
        name.list$cap[grep("Golgi", name.list$cap)] <- "Golgi Apparatus"
        name.list$cap[grep("Extracellular", name.list$cap)] <- "Extracellular Region or Secreted"
        name.list$cap[grep("Pmembrane", name.list$cap)] <- "Plasma Membrane"
        name.list$cap[grep("Er", name.list$cap)] <- "ER"
        cc.breast_plotting[[name]]$Names <- name.list$cap[match(cc.breast_plotting[[name]]$Category, name.list$Names)]
        #calculate average score to select top 5 categories
        cc.breast_plotting[[name]]$avg.score <- NA
        for (ii in 1:length(unique(df$Status))) {
          dff <- df[df$Status == unique(df$Status)[ii],]
              for (iii in 1:length(unique(cc.breast_plotting[[name]]$Category))) {
                cc.breast_plotting[[name]]$avg.score[cc.breast_plotting[[name]]$Status == unique(dff$Status) & 
                cc.breast_plotting[[name]]$Category == unique(cc.breast_plotting[[name]]$Category)[iii]] <- mean(dff$Score[dff$Category == unique(cc.breast_plotting[[name]]$Category)[iii]])
              } 
        }              
    }
}
#merge all plotting dataframe
breast_df <- do.call(rbind,cc.breast_plotting)

#combine all categorized pathways
cc.22RV1 <- rbind(
        cc.normal[,c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes", "Status", "Score", "input.val", "GO.num", "Category", "percent")],
        cc.resist[,c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes", "Status", "Score", "input.val", "GO.num", "Category", "percent")]
      )
df_all <- do.call(rbind, cc.breast.processed)
#select only the pathways of cell lines with both EV and cell
pattern <- c("BT549", "MCF7", "HCC1954", "HCC1419")
df_all <- rbind(df_all[grep(paste(pattern, collapse = "|"), df_all$Status),], cc.22RV1)

#merge all aggregaed plotting dataframe
cc.all_plotting <- rbind(breast_df[grep(paste(pattern, collapse = "|"),breast_df$Status),], 
        cc.normal_plotting, cc.resist_plotting)
cc.all_plotting$label <- paste(cc.all_plotting$Status, cc.all_plotting$faucet, sep = "_")
#could use cc.all_plotting to plot all complete categories

#identify top 5 categories & new plotting position & percentage
cc.all_5_plotting <- data.frame(cc.all_plotting %>% group_by(Status) %>% arrange(desc(percent)) %>% slice(1:5))
other_plotting <- filter(cc.all_plotting, !cc.all_plotting$label %in% cc.all_5_plotting$label)
for (i in 1:length(unique(other_plotting$Status))) {
    df <- other_plotting[other_plotting$Status == unique(other_plotting$Status)[i],]
    new_df <- data.frame(Category = "Other", Status = unique(other_plotting$Status)[i],
      percent = sum(df$percent), position = df$position[1], faucet = df$faucet[1], 
      Names = "Other", avg.score = mean(df$avg.score), label = df$label[1], weighted.score = mean(df$weighted.score))
    if ("Other" %in% cc.all_5_plotting$Names[cc.all_5_plotting$Status == unique(other_plotting$Status)[i]]) {
         dff <- cc.all_5_plotting[cc.all_5_plotting$Status == unique(other_plotting$Status)[i],]
         ii <- dff[dff$Names == "Other",]
         cc.all_5_plotting$percent[cc.all_5_plotting$Status == unique(other_plotting$Status)[i] & 
                      cc.all_5_plotting$Names == "Other"] <- ii$percent + new_df$percent
         cc.all_5_plotting$position[cc.all_5_plotting$Status == 
              unique(other_plotting$Status)[i]] <- data.frame(aggregate(percent~Category+Status,
                     cc.all_5_plotting[cc.all_5_plotting$Status == unique(other_plotting$Status)[i],],sum) %>% 
               arrange(Status, desc(percent)) %>% group_by(Status) %>% 
               mutate(position = 100 - (cumsum(percent) - percent/2))) %>% pull(position)
    } else {
      cc.all_5_plotting <- rbind(cc.all_5_plotting, new_df)
    }
}
#check for sum for new aggregated plotting dataframe
for (i in 1:length(unique(cc.all_5_plotting$Status))) {
    df <- cc.all_5_plotting[cc.all_5_plotting$Status == unique(cc.all_5_plotting$Status)[i],]
    print(paste(sum(df$percent), unique(df$Status), sep = "_"))
}

#rename different cell lines/ conditions
cc.all_5_plotting$Status <- gsub("_Cell|_cell", " Cells", cc.all_5_plotting$Status)
cc.all_5_plotting$Status <- gsub("_EV", " Cells EVs", cc.all_5_plotting$Status)
cc.all_5_plotting$Status <- gsub("_MS", "", cc.all_5_plotting$Status)
cc.all_5_plotting$Status <- gsub("Normal", "22RV1", cc.all_5_plotting$Status)
cc.all_5_plotting$Status <- gsub("Resistant", "22RV1DR", cc.all_5_plotting$Status)
cc.all_5_plotting$Status <- gsub("hurex", "HuRex", cc.all_5_plotting$Status)
cc.all_5_plotting$Status <- gsub("Patient Cells EVs", "Patient Serum EVs", cc.all_5_plotting$Status)
cc.all_5_plotting$Status <- factor(cc.all_5_plotting$Status, 
          levels = c("BT549 Cells", "BT549 Cells EVs", "HCC1419 Cells", "HCC1419 Cells EVs", 
          "HCC1954 Cells", "HCC1954 Cells EVs", "MCF7 Cells", "MCF7 Cells EVs", 
          "22RV1 Cells", "22RV1 Cells EVs", "22RV1DR Cells", "22RV1DR Cells EVs", 
          "HuRex","Patient Serum EVs"))
#line break
cc.all_5_plotting$Names[grep("Extracellular", cc.all_5_plotting$Names)] <- "Extracellular\nRegion or Secreted"

### plot pie chart
#with labels
ggplot(cc.all_5_plotting, aes(x="",y = percent, fill = Names, group = faucet), color = "white") + 
    geom_bar(width = 1, stat = "identity", color = "black") + theme_void() + 
    theme(axis.text.x = element_blank()) + 
    geom_label_repel(aes(x = 1.45, y = position, label = sprintf("%0.2f%%", round(percent, digits = 2))), 
      show.legend = FALSE, max.overlaps = 40, size = 7) + coord_polar("y", start = 0) +
    labs(fill = "Cellular Component") +
    scale_fill_manual("legend", values = col[names(col) %in% unique(cc.all_5_plotting$Names)]) + 
    facet_wrap(~ Status, ncol = 2) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="right",
      legend.title=element_blank(), strip.text = element_text(size=rel(2.5)),
      legend.text=element_text(size=rel(1.5)))
#without labels
ggplot(cc.all_5_plotting, aes(x="",y = percent, fill = Names, group = faucet), color = "white") + 
    geom_bar(width = 1, stat = "identity", color = "black") + theme_void() + 
    theme(axis.text.x = element_blank()) + 
    #geom_label_repel(aes(x = 1.45, y = position, label = sprintf("%0.2f%%", round(percent, digits = 2))), 
    #  show.legend = FALSE, max.overlaps = 40, size = 7) + 
    coord_polar("y", start = 0) +
    labs(fill = "Cellular Component") +
    scale_fill_manual("legend", values = col[names(col) %in% unique(cc.all_5_plotting$Names)]) + 
    facet_wrap(~ Status, ncol = 2) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",
      legend.title=element_blank(), strip.text = element_text(size=rel(2.5)),
      legend.text=element_text(size=rel(1.5)))


### ++++++++++++++++++++++++++++++++++++
### FIGURE 3B : MASSPEC BOXPLOT
### ++++++++++++++++++++++++++++++++++++
#pathways that best represents each categories
key.terms <- c("focal adhesion (GO:0005925)",
              "nucleus (GO:0005634)",
              "intracellular organelle lumen (GO:0070013)",
              "early endosome (GO:0005769)",
              "cytoplasmic vesicle membrane (GO:0030659)",
              "mitochondrial inner membrane (GO:0005743)",
              "ribosome (GO:0005840)")
#separate key pathways from each cell lines
cc.all_boxplot <- df_all[df_all$Term %in% key.terms,]
cc.all_boxplot <- cc.all_boxplot[,c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes", "Status", "Score", "input.val", "GO.num", "Category", "percent")]
#rename categories
cc.all_boxplot$Status <- gsub("_Cell|_cell", " Cells", cc.all_boxplot$Status)
cc.all_boxplot$Status <- gsub("_EV", " Cells EVs", cc.all_boxplot$Status)
cc.all_boxplot$Status <- gsub("_MS", "", cc.all_boxplot$Status)
cc.all_boxplot$Status <- gsub("Normal", "22RV1", cc.all_boxplot$Status)
cc.all_boxplot$Status <- gsub("Resistant", "22RV1DR", cc.all_boxplot$Status)
cc.all_boxplot$Status <- factor(cc.all_boxplot$Status, levels = unique(cc.all_boxplot$Status))
cc.all_boxplot$Cell_EV <- ifelse(grepl("EVs|ex", cc.all_boxplot$Status), "EVs", "Cells")
cc.all_boxplot$Cell_lines <- sub(" .*", "", cc.all_boxplot$Status)
cc.all_boxplot$Score <- abs(cc.all_boxplot$Score)


### plot boxplot
ggplot(data=cc.all_boxplot, aes(x=Cell_EV , y=Score)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.7) + #ylim(0,20) +
    theme_classic() + #geom_jitter(width = 0.2, aes(color=Cell_lines)) + 
    facet_wrap(~Term, ncol = 2, scales = "free_y", dir = "v",labeller = label_wrap_gen(width = 18, multi_line = TRUE)) + 
    geom_point(aes(x=Cell_EV , y=Score, color=Cell_lines), position = position_dodge(0.5), size = 5) +
    #geom_line(aes(group=Cell_lines),color="slategray4",alpha=0.7,position = position_dodge(0.3)) + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                       paired=F, label="p.signif", hide.ns = FALSE, vjust = -0.2,
                       bracket.size = 1, size = 8) +
    scale_y_continuous(expand = c(.17,0)) + #geom_blank(aes(y=max.y)) +         
    labs(color = "Cell Lines", x = "", y = "-log10(FDR)") + 
    theme(axis.title = element_text(size=rel(2.5)), legend.title = element_text(size=rel(2.5)), 
          legend.text = element_text(size=rel(2.5)), 
          strip.text = element_text(size=rel(2.5)), axis.text = element_text(size=rel(2)))    



### ++++++++++++++++++++++++++++++++++++
### RNASEQ : PATIENT SERUM
### ++++++++++++++++++++++++++++++++++++
#load counts
allRNA.bowtie1.bestrata.m3n1.counts <- read.csv("Patient_RNAseq_raw_counts.csv", row.names = 1)

#make phenotype
patient_rnaseq_phenotype <- data.frame(sample = colnames(allRNA.bowtie1.bestrata.m3n1.counts),
									   patient = as.character(substr( colnames(allRNA.bowtie1.bestrata.m3n1.counts),2,4)))
patient_rnaseq_phenotype$technology[grep("DLD", patient_rnaseq_phenotype$sample)] <- "DLD"
patient_rnaseq_phenotype$technology[grep("UC", patient_rnaseq_phenotype$sample)] <- "UC"
patient_rnaseq_phenotype$condition[grep("Serum",patient_rnaseq_phenotype$sample)] <- "Serum"
patient_rnaseq_phenotype$variable <- paste(patient_rnaseq_phenotype$condition,patient_rnaseq_phenotype$technology,sep="_")
#double check colnames of raw counts match with sample names in phenotype
identical( colnames(allRNA.bowtie1.bestrata.m3n1.counts),as.character(patient_rnaseq_phenotype$sample))

#filtering
x <- allRNA.bowtie1.bestrata.m3n1.counts
cpm <- cpm(x)
keep.rows <- rowSums(cpm>10)>=1
rm(x,cpm)

#add symbol
allRNA.bowtie1.bestrata.m3n1.counts$geneSymbols <- mapIds(org.Hs.eg.db, keys = rownames(allRNA.bowtie1.bestrata.m3n1.counts), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
allRNA.bowtie1.bestrata.m3n1.counts$ensembl <- rownames(allRNA.bowtie1.bestrata.m3n1.counts)
#make separate count/gene dataframe
patient_rnaseq_counts <- allRNA.bowtie1.bestrata.m3n1.counts[,colnames(allRNA.bowtie1.bestrata.m3n1.counts) %in% patient_rnaseq_phenotype$sample]
patient_rnaseq_genes <- allRNA.bowtie1.bestrata.m3n1.counts[,c("geneSymbols", "ensembl")]

#make DGElist
dcpmx2 <- DGEList(counts=patient_rnaseq_counts[keep.rows,],group=patient_rnaseq_phenotype$condition, genes = patient_rnaseq_genes[keep.rows,])
dcpmx2 <- calcNormFactors(dcpmx2, method = "RLE") #method TMM
identical(rownames(dcpmx2$samples),as.character(patient_rnaseq_phenotype$sample))
patient_rnaseq_phenotype$norm.factors <- dcpmx2$samples$norm.factors
#make design matrix
design = model.matrix( ~ 0 + technology, data=patient_rnaseq_phenotype)
#voom normalization
v_m3n1_2 <- voom(dcpmx2, design, plot=TRUE, save.plot = TRUE)


### ++++++++++++++++++++++++++++++++++++
### MASSPEC : PATIENT SERUM
### ++++++++++++++++++++++++++++++++++++
#load counts
patient_masspec <- read.csv("Patient_Mass.Spec_raw_counts.csv", row.names = 1)
patient_masspec$uniprot <- rownames(patient_masspec)
patient_masspec$ensembl <- as.character(mapIds(org.Hs.eg.db, keys = as.character(patient_masspec$uniprot), column="ENSEMBL", keytype="UNIPROT", multiVals="first")) 
patient_masspec$symbol <- as.character(mapIds(org.Hs.eg.db, keys = as.character(patient_masspec$uniprot), column="SYMBOL", keytype="UNIPROT", multiVals="first")) 
#remove genes that have 0 counts across all samples
table(rowSums(patient_masspec[,grep("BioSample", colnames(patient_masspec))]) != 0)
patient_masspec <- patient_masspec[rowSums(patient_masspec[,grep("BioSample", colnames(patient_masspec))]) != 0,]
#make count files and gene files
patient_data <- patient_masspec[,grep("BioSample", colnames(patient_masspec))]
patient_genes <- patient_masspec[,!grepl("BioSample", colnames(patient_masspec))]

#filter
keep.masspec <- rowSums(patient_data>1)>=1
#create DGElist
dge_mspec <- DGEList(counts=patient_data[keep.masspec,], genes = patient_genes[keep.masspec,]) 
dge_mspec <- calcNormFactors(dge_mspec, method = "RLE") 
dge_mspec <- estimateCommonDisp(dge_mspec, verbose=TRUE)
dge_mspec <- estimateTagwiseDisp(dge_mspec)
#voom normalization
v_masspec_patient <- voom(dge_mspec, plot = TRUE,save.plot = TRUE) 

#make phenotype
phenotype_patient <- data.frame(samples = colnames(patient_data))
phenotype_patient$technology <- ifelse(grepl("UC", phenotype_patient$samples), "UC", "DLD")


### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### FIGURE 5D : PATIENT TRANSCRIPTOMIC HEATMAP
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#double check the column names of normalized counts and phenotype samples match
identical(colnames(v_m3n1_2$E), patient_rnaseq_phenotype$sample)

#use normalized data for heatmap matrix
generic_matrix <- v_m3n1_2$E
rownames(generic_matrix) <- v_m3n1_2$genes$gene_name
rownames(generic_matrix) <- paste(v_m3n1_2$genes$geneSymbols,1:nrow((generic_matrix)),sep="_")
#make annotation column
annotation_col <- data.frame (Condition = patient_rnaseq_phenotype$condition,
                              Patient = patient_rnaseq_phenotype$patient,
                              Technology = patient_rnaseq_phenotype$technology)
rownames(annotation_col) <- patient_rnaseq_phenotype$sample
#double check column names of heatmap matrix and rownames of annotation column match
table(colnames(generic_matrix) %in% rownames(annotation_col))
identical(colnames(generic_matrix),rownames(annotation_col))
#assign annotation colors
annotation_colors <- list ( Condition = c(Serum="red",Urine="blue"),
                            Technology = c (DLD="purple",UC="green") )
#make annotation row
ann_row <- data.frame( nanoDLD = rowMeans(generic_matrix[,which(annotation_col$Technology %in% 'UC')] ) ,
                       UC = rowMeans(generic_matrix[,which(annotation_col$Technology %in% 'DLD')] ) )
rownames(ann_row) <- rownames(generic_matrix)
#select top 30
iy <- order(ann_row$nanoDLD,decreasing = TRUE)[1:100][which(order(ann_row$nanoDLD,decreasing = TRUE)[1:100] %in% order(ann_row$UC,decreasing = TRUE)[1:100])]



### plot heatmap
pheatmap::pheatmap(generic_matrix[iy,], color = colorRampPalette(c("steelblue", "white","firebrick"))(255),  
                   cluster_cols = T, 
                   cluster_rows=T, 
                   main ='RNAseq Patients',
                   name='Ave.Expr.',
                   annotation_col = annotation_col,
                   annotation_row = ann_row,
                   annotation_colors = annotation_colors,
                   border_color="white", 
                   show_rownames = TRUE, 
                   show_colnames = FALSE, 
                   clustering_distance_rows ="euclidean", scale="row")

### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### FIGURE 5E : PATIENT PROTEOMIC HEATMAP
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#use normalized data for heatmap matrix
generic_matrix <- v_masspec_patient$E
rownames(generic_matrix) <- v_masspec_patient$genes$symbol 

#make annotation column
annotation_col <- data.frame (Condition = phenotype_patient$technology)
rownames(annotation_col) <- phenotype_patient$samples
#double check column names of heatmap matrix and rownames of annotation column match
identical(colnames(generic_matrix),rownames(annotation_col))
#relabel condition
annotation_col$Condition[ annotation_col$Condition %in% 'DLD'] <- 'Serum nanoDLD'
annotation_col$Condition[ annotation_col$Condition %in% 'UC'] <- 'Serum UC'

#make annotation row
ann_row <- data.frame( nanoDLD = rowMeans(generic_matrix[,which(annotation_col$Condition %in% 'Serum nanoDLD')]) ,
                       UC = rowMeans(generic_matrix[,which(annotation_col$Condition %in% 'Serum UC')]) )
rownames(ann_row) <- rownames(generic_matrix)
#select top 30
iy <- order(ann_row$nanoDLD,decreasing = TRUE)[1:30]

### plot heatmap
pheatmap::pheatmap(generic_matrix[iy], color = colorRampPalette(c("steelblue", "white","firebrick"))(255),    
                   cluster_cols = T, 
                   cluster_rows=T, 
                   main ='Masspec EVs',
                   name='Ave.Expr',
                   annotation_col = annotation_col,
                   annotation_row = ann_row,
                   #annotation_colors = annotation_colors,
                   border_color="white", 
                   show_rownames = TRUE, 
                   show_colnames = FALSE, 
                   clustering_distance_rows ="euclidean", scale="none")

