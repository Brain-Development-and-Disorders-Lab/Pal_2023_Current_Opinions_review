library(SingleCellExperiment)
library(scater)
library(scran)

#Load in ALM exon counts table
sparse.mat <- readSparseCounts("mouse_ALM_gene_expression_matrices_2018-06-14/mouse_ALM_2018-06-14_exon-matrix.csv",
                               sep=",", quote='"', row.names=TRUE, col.names=TRUE) #45768 genes and 10068 cells
#Rows contain EntrezID; columns contain cells

Gene.meta <- read.csv("mouse_ALM_gene_expression_matrices_2018-06-14/mouse_ALM_2018-06-14_genes-rows.csv",
                      header = TRUE)
rownames(Gene.meta) <- Gene.meta$gene_entrez_id
Sample.meta.ALM <- read.csv("mouse_ALM_gene_expression_matrices_2018-06-14/mouse_ALM_2018-06-14_samples-columns.csv",
                        header = TRUE)
rownames(Sample.meta.ALM) <- Sample.meta.ALM$sample_name

ALM.sce <- SingleCellExperiment(assays=list(counts=sparse.mat),
                                    colData=Sample.meta.ALM, rowData=Gene.meta)
ALM.sce <- ALM.sce[,ALM.sce$core_intermediate_call!='Not Calculated']
#9573 samples remaining as reported in Tasic 2018
rm(Gene.meta, Sample.meta.ALM, sparse.mat)

#Load in VISp exon counts table
sparse.mat <- readSparseCounts("mouse_VISp_gene_expression_matrices_2018-06-14/mouse_VISp_2018-06-14_exon-matrix.csv",
                               sep=",", quote='"', row.names=TRUE, col.names=TRUE)#45768 genes and 15413 cells

Gene.meta <- read.csv("mouse_VISp_gene_expression_matrices_2018-06-14/mouse_VISp_2018-06-14_genes-rows.csv",
                      header = TRUE)
rownames(Gene.meta) <- Gene.meta$gene_entrez_id
Sample.meta.VISp <- read.csv("mouse_VISp_gene_expression_matrices_2018-06-14/mouse_VISp_2018-06-14_samples-columns.csv",
                        header = TRUE)
rownames(Sample.meta.VISp) <- Sample.meta.VISp$sample_name

VISp.sce <- SingleCellExperiment(assays=list(counts=sparse.mat),
                                colData=Sample.meta.VISp, rowData=Gene.meta)
VISp.sce <- VISp.sce[,VISp.sce$core_intermediate_call!='Not Calculated'] #14249 samples remaining as reported in Tasic 2018
rm(Gene.meta, Sample.meta.VISp, sparse.mat)

#Combined ALM and VISp SCE objects
Full.SCE <- cbind(ALM.sce, VISp.sce)
rm(ALM.sce, VISp.sce)

#Subset glutamatergic cells
Glut.SCE <- Full.SCE[,Full.SCE$class=='Glutamatergic'] #11905 glutamatergic cells

#Run normalization, variance modelling, feature selection, PCA, TSNE
set.seed(1)
clusters1 <- quickCluster(Glut.SCE)
Glut.SCE <- computeSumFactors(Glut.SCE, cluster=clusters1)
Glut.SCE <- logNormCounts(Glut.SCE)
dec.all <- modelGeneVarByPoisson(Glut.SCE)
topHVG.all <- getTopHVGs(dec.all, n=2000) #2000 genes subsetted
Glut.SCE <- denoisePCA(Glut.SCE, technical=dec.all, subset.row=topHVG.all)
Glut.SCE <- runTSNE(Glut.SCE, dimred = "PCA")
saveRDS(Glut.SCE, file="Glutamatergic_cells.rds")

#Subset IT cells
IT.subset <- Glut.SCE[,Glut.SCE$subclass=='L2/3 IT'|Glut.SCE$subclass=='L4'|Glut.SCE$subclass=='L5 IT'|Glut.SCE$subclass=='L6 IT']
IT.subset$Subclass <- droplevels(IT.subset$Subclass)

#FigA: Color based on brain region dissected for scRNA-seq
FigA_meta <- read.csv("metadata/FigA_meta.csv", header=TRUE)
rownames(FigA_meta) <- FigA_meta$sample_name
FigA_meta <- FigA_meta[match(colnames(IT.subset), rownames(FigA_meta)),]
IT.subset$Brain_region <- as.factor(FigA_meta$FigA_label)
plotTSNE(IT.subset, colour_by='Brain_region') + scale_colour_manual(
  values = c("grey70","black","lightpink","red"))

#FigB: Color based on IT cell subclass
FigB_meta <- read.csv("metadata/FigB_meta.csv", header=TRUE)
rownames(FigB_meta) <- FigB_meta$sample_name
FigB_meta <- FigB_meta[match(colnames(IT.subset), rownames(FigB_meta)),]
IT.subset$Subclass <- as.factor(FigB_meta$FigB_label)
plotTSNE(IT.subset, colour_by='Subclass') + scale_colour_manual(
  values = c("black","cadetblue1","darkseagreen1",
             "lightpink","thistle2","red"))

#FigC: IT cells based on clusters
FigC_meta <- read.csv("metadata/FigC_meta.csv", header=TRUE)
rownames(FigC_meta) <- FigC_meta$sample_name
FigC_meta <- FigC_meta[match(colnames(IT.subset), rownames(FigC_meta)),]
IT.subset$Cluster <- as.factor(FigC_meta$FigC_label)
plotTSNE(IT.subset, colour_by='Cluster') + scale_colour_manual(
  values = c("black","gold","khaki1","yellow3","lightcyan3","lightcyan2",
             "skyblue","darkseagreen1","palevioletred1","plum1","mediumorchid3","purple",
             "mediumpurple","darksalmon","violet","maroon1","plum3","olivedrab1",
             "olivedrab3","springgreen","green","green4","cyan","deepskyblue",
             "royalblue1","skyblue1","turquoise","cyan4","blue","red"))

#FigD: Besides retrograde-labelled cells from ALM-c and VISp-C, also include cells labelled from ipsilateral cortex and contralateral striatum
FigD_meta <- read.csv("metadata/FigD_meta.csv", header=TRUE)
rownames(FigD_meta) <- FigD_meta$sample_name
FigD_meta <- FigD_meta[match(colnames(IT.subset), rownames(FigD_meta)),]
IT.subset$Retrograde <- as.factor(FigD_meta$FigD_label)
plotTSNE(IT.subset, colour_by='Retrograde') + scale_colour_manual(
  values = c("black","blue","cyan","gray90","gray90","gray90","gray90","gray90",
             "gray90","gray90","gray90","gray90","gray90","gray90",
             "gray90","gray90","gray90","gray90","gray90","gray90",
             "gray90","gray90","gray90","gray90","gray90","gray90",
             "gray90","gray90","gray90","gray90","gray90","red"))
