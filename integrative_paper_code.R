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
### FIGURE 4E - 1G: SPEARMAN CORRELATIONS
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


### ++++++++++++++++++++++++++++++++++++
### EV MASSPEC : 22RV1/22RV1DR
### ++++++++++++++++++++++++++++++++++++
mass_spec_results <- read.csv("EV_Mass.Spec_raw_counts.csv", row.names = T)
#Separate values from main table
mass_spec_data <- mass_spec_results
mass_spec_genes <- data.frame(uniprot_id = rownames(mass_spec_results))
#add gene symbol
mass_spec_genes$gene_symbol <- as.character(mapIds(org.Hs.eg.db, keys = as.character(mass_spec_genes$uniprot_id), column="SYMBOL", keytype="UNIPROT", multiVals="first")) 

### filter protein counts
table(rowSums(mass_spec_data>3)>=1)
keep.masspec <- rowSums(mass_spec_data>3)>=1

### Create a DGEList object
dge_mspec <- DGEList(counts=mass_spec_data[keep.masspec,], genes = mass_spec_genes[keep.masspec,]) 
dge_mspec <- calcNormFactors(dge_mspec, method = "RLE") 
dge_mspec <- estimateCommonDisp(dge_mspec, verbose=TRUE)
dge_mspec <- estimateTagwiseDisp(dge_mspec)

### voom normalization
v_masspec <- voom(dge_mspec, plot = TRUE,save.plot = TRUE) 
#design
masspec_meta <- data.frame(type=c("RESISTANT","RESISTANT","SENSITIVE","SENSITIVE"),sample=c("22RVIDR_1","22RVIDR_2","22RVI_1","22RVI_2"))
masspec_design = model.matrix( ~ 0 + type, data=masspec_meta)
colnames(masspec_design) <- levels(masspec_meta$type)
contr.matrix <- makeContrasts( "RESISTANT_vs_SENSITIVE" = RESISTANT - SENSITIVE,levels = colnames(masspec_design))
vfit <- lmFit(v_masspec, masspec_design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
summary(decideTests(efit))
#output results
results_masspec_DR_vs_SENS <- topTable(efit,n=Inf,coef=1)







