Args <- commandArgs()
print(Args)

library('Seurat')
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(clustree)

varGenes <- as.integer(Args[6])
CCAdim <- as.integer(Args[7])
PCAdim <- as.integer(Args[8])
#Resolution:
ResN <- Args[9]



PATH <- '/home/nifang/mmyy/project/WY_scRNA/seuratinteg'


outDir <- paste0(PATH,'/varGenes',Args[6],'_CCA',Args[7],'_PCA',Args[8],'_Res',Args[9])
if (!dir.exists(outDir)){
  dir.create(outDir)}
FilePrefix=paste0(outDir,'/NKcombined')

ctrl.data <- Read10X(data.dir='/home/nifang/mmyy/project/WY_scRNA/cellranger/con/outs/filtered_feature_bc_matrix/')
drug.data <- Read10X(data.dir='/home/nifang/mmyy/project/WY_scRNA/cellranger/ven/outs/filtered_feature_bc_matrix/')

#control
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "NK_ctrl", min.cells = 5)
ctrl$stim <- "CTRL"
ctrl[['percent.mt']] <- PercentageFeatureSet(object=ctrl, pattern='^MT-')
pdf(file=paste0(outDir,'/QC_ctrl.pdf'))
VlnPlot(ctrl,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
dev.off()
ctrl <- subset(ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 10)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = varGenes)
top10 <- head(VariableFeatures(ctrl), 10)
pdf(file=paste0(outDir,'/VariableFeatures_ctrl.pdf'))
plot1 <- VariableFeaturePlot(ctrl)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.off()

#drug
drug <- CreateSeuratObject(counts = drug.data, project = "NK_drug", min.cells = 5)
drug$stim <- "DRUG"
drug[['percent.mt']] <- PercentageFeatureSet(object=drug, pattern='^MT-')
pdf(file=paste0(outDir,'/QC_drug.pdf'))
VlnPlot(drug,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
dev.off()
drug <- subset(drug, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 10)
drug <- NormalizeData(drug, verbose = FALSE)
drug <- FindVariableFeatures(drug, selection.method = "vst", nfeatures = varGenes)
top10 <- head(VariableFeatures(drug), 10)
pdf(file=paste0(outDir,'/VariableFeatures_drug.pdf'))
plot1 <- VariableFeaturePlot(drug)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.off()


#Integration
NK.anchors <- FindIntegrationAnchors(object.list = list(ctrl, drug), dims = 1:CCAdim)
NK.combined <- IntegrateData(anchorset = NK.anchors, dims = 1:CCAdim)

DefaultAssay(NK.combined) <- "integrated"

DefaultAssay(NK) <- "RNA"

#Scale&PCA&Cluster
NK.combined <- ScaleData(NK.combined, verbose = FALSE)
NK.combined <- RunPCA(NK.combined, npcs = 40, verbose = FALSE)
print("Run PCA Done!")
pdf(file=paste0(outDir,'/ElbowPlot.PCA.pdf'))
ElbowPlot(object = NK.combined, ndims = 50)
dev.off()
NK.combined <- RunUMAP(NK.combined, reduction = "pca", dims = 1:PCAdim)
NK.combined <- FindNeighbors(NK.combined, reduction = "pca", dims = 1:PCAdim)

#resolution:
#saveRDS(NK.combined,paste0(outDir,'/NKcombined_no_Resolution.rds'))
#obj <- FindClusters(NK.combined, resolution = seq(0.1,0.7,by=0.1))
#pdf(file=paste0(outDir,'/Clustree_resolution.pdf'))
#clustree(obj)
#dev.off()

NK.combined <- FindClusters(NK.combined, resolution = as.numeric(ResN))

UMAP=as.data.frame(Embeddings(object = NK.combined, reduction = "umap"))
write.table(UMAP,file=paste0(FilePrefix,'.umap.txt'),sep='\t',quote=F)
print("Save UMAP Done!")

#Plot Dimplot
pdf(file=paste0(outDir,'/Dimplot_','condition.pdf'))
DimPlot(NK.combined, reduction = "umap", group.by = "stim")
dev.off()
pdf(file=paste0(outDir,'/Dimplot_','condition2.pdf'))
DimPlot(NK.combined, reduction = "umap", split.by = "stim")
dev.off()
pdf(file=paste0(outDir,'/Dimplot_','cluster.pdf'))
DimPlot(NK.combined, reduction = "umap", label = TRUE)
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_','genes.pdf'))
FeaturePlot(NK.combined, features = c("CD3D", "IL32", "IFI6", "KLRC2", "GNLY", "HLA-DQA1", "FCGR3A", "CD52", "FCER1G","SPON2","XCL1","MKI67"), min.cutoff = "q9")
dev.off()

#Save MetaData
MetaData=as.data.frame(NK.combined@meta.data)
write.table(MetaData,file=paste0(FilePrefix,'.MetaData.txt'),sep='\t',quote=F)

#SaveExpData
#Data=as.data.frame(as.matrix(GetAssayData(object = Aging.integrated)))
#write.table(Data,file=paste0(FilePrefix,'.dataNorm.txt'),sep='\t',quote=F)

#VarGenes=VariableFeatures(object = Aging.integrated)
#VarGeneData=Data[VarGenes,]
#write.table(VarGeneData,file=paste0(FilePrefix,'.VarGeneData.Integrate.txt'),sep='\t',quote=F)


###save file
NK.combined$Group <- paste(NK.combined$stim, NK.combined$seurat_clusters, sep = "_")
saveRDS(NK.combined,paste0(outDir,'/NKcombined.rds'))


#Marker Gene
DefaultAssay(NK.combined) <- "RNA"
nk.markers <- FindConservedMarkers(NK.combined, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
write.table(nk.markers,file=paste0(outDir,'/ClusterMarkers_0.txt'),quote =F,sep='\t')
nk.markers <- FindConservedMarkers(NK.combined, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
write.table(nk.markers,file=paste0(outDir,'/ClusterMarkers_1.txt'),quote =F,sep='\t')
nk.markers <- FindConservedMarkers(NK.combined, ident.1 = 2, grouping.var = "stim", verbose = FALSE)
write.table(nk.markers,file=paste0(outDir,'/ClusterMarkers_2.txt'),quote =F,sep='\t')
nk.markers <- FindConservedMarkers(NK.combined, ident.1 = 3, grouping.var = "stim", verbose = FALSE)
write.table(nk.markers,file=paste0(outDir,'/ClusterMarkers_3.txt'),quote =F,sep='\t')
nk.markers <- FindConservedMarkers(NK.combined, ident.1 = 4, grouping.var = "stim", verbose = FALSE)
write.table(nk.markers,file=paste0(outDir,'/ClusterMarkers_4.txt'),quote =F,sep='\t')
nk.markers <- FindConservedMarkers(NK.combined, ident.1 = 5, grouping.var = "stim", verbose = FALSE)
write.table(nk.markers,file=paste0(outDir,'/ClusterMarkers_5.txt'),quote =F,sep='\t')
nk.markers <- FindConservedMarkers(NK.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
write.table(nk.markers,file=paste0(outDir,'/ClusterMarkers_6.txt'),quote =F,sep='\t')

#DEGs
Idents(NK.combined) <- "Group"
Group <- FindMarkers(NK.combined, ident.1 = "DRUG_0", ident.2 = "CTRL_0", verbose = FALSE)
write.table(Group, file=paste0(outDir,'/DEG_0.txt'),quote =F,sep='\t')
Group <- FindMarkers(NK.combined, ident.1 = "DRUG_1", ident.2 = "CTRL_1", verbose = FALSE)
write.table(Group, file=paste0(outDir,'/DEG_1.txt'),quote =F,sep='\t')
Group <- FindMarkers(NK.combined, ident.1 = "DRUG_2", ident.2 = "CTRL_2", verbose = FALSE)
write.table(Group, file=paste0(outDir,'/DEG_2.txt'),quote =F,sep='\t')
Group <- FindMarkers(NK.combined, ident.1 = "DRUG_3", ident.2 = "CTRL_3", verbose = FALSE)
write.table(Group, file=paste0(outDir,'/DEG_3.txt'),quote =F,sep='\t')
Group <- FindMarkers(NK.combined, ident.1 = "DRUG_4", ident.2 = "CTRL_4", verbose = FALSE)
write.table(Group, file=paste0(outDir,'/DEG_4.txt'),quote =F,sep='\t')
Group <- FindMarkers(NK.combined, ident.1 = "DRUG_5", ident.2 = "CTRL_5", verbose = FALSE)
write.table(Group, file=paste0(outDir,'/DEG_5.txt'),quote =F,sep='\t')
Group <- FindMarkers(NK.combined, ident.1 = "DRUG_6", ident.2 = "CTRL_6", verbose = FALSE)
write.table(Group, file=paste0(outDir,'/DEG_6.txt'),quote =F,sep='\t')
Group <- FindMarkers(NK.combined, ident.1 = "DRUG_7", ident.2 = "CTRL_7", verbose = FALSE)
write.table(Group, file=paste0(outDir,'/DEG_7.txt'),quote =F,sep='\t')

#nk.markers <- FindConservedMarkers(NK.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
#write.table(nk.markers,file=paste0(outDir,'/ClusterMarkers_7.txt'),quote =F,sep='\t')


#ConserveMarker heatmap
#ConserveMarker <- read.table('/Users/ceci/Project/汪焱_NK_scRNAseq/Seurat/varGenes2000_CCA20_PCA20_Res0.18/ClusterMarkers_merge.txt', head=TRUE)
#Idents(NK.combined) <- "seurat_clusters"
#DefaultAssay(NK.combined) <- "RNA"
#pdf(file=paste0(outDir,'/ClusterMarkers_Top10_heatmap.pdf'))
#DoHeatmap(NK.combined, features = ConserveMarker$gene) + NoLegend()
#dev.off()

#VinPlot: gene exp in cluster
#pdf(file=paste0(outDir,'/GeneExp/VlnPlot.FCGR3A.pdf'))
#VlnPlot(NK.combined, features = c("FCGR3A"), slot = "data")
#dev.off()


gene <- 'ITGAL'
pdf(file=paste0(Dir,gene,'.pdf'))
FeaturePlot(object = NK, features = gene)+
scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdYlBu")))
dev.off()


pdf(file=paste0('/Users/ceci/Project/汪焱_NK_scRNAseq/Seurat/varGenes2000_CCA20_PCA20_Res0.18/Files_WY_NKscRNA/VlnPlot/VlnPlot.','FCGR3A.pdf'))
VlnPlot(NK, features = c("FCGR3A"), slot = "data")
dev.off()


markers.to.plot <- c("CD3D", "FCGR3A", "VMO1", "CCL2", "S100A9")

markers.to.plot <- c('ATP5F1E','ATP5MC2','ATP5PF','COX5B','COX6A1','COX6C','COX7B','COX7C','ND6','NDUFA4','NDUFAB1','NDUFB6','UQCRB','UQCRH','UQCR11','UQCR10','MPC2','GSTP1','PRDX1','PSMA4','PSMB1','PSMB2','PSME1','TXN','LAMTOR5','ELOB')



pdf(file=paste0(Dir,'DotPlot_ATP.pdf'))
DotPlot(NK, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8,
    split.by = "stim") + RotatedAxis()
dev.off()


Idents(NK_3) <- "stim"

pdf(file=paste0(Dir,'FCGR3A.pdf'),width=2.5, height=2.5)
VlnPlot(NK_3, features = c("FCGR3A"), slot = "data",  pt.size = 0, cols=c('red', 'orange'))+NoLegend()+theme(text = element_text(size = 12))
dev.off()



pdf(file=paste0(Dir,'VlnPlot_subset_gene1_16.pdf'),width=10, height=8)
VlnPlot(NK, features = g1,pt.size=0)
dev.off()

NK <- readRDS('/Users/ceci/Project/汪焱_NK_scRNAseq/Seurat/varGenes2000_CCA20_PCA20_Res0.18/NKcombined.rds')
NK_3 <- subset(NK, subset=seurat_clusters==3)
Idents(NK_3) <- "stim"

i<-"AP1G1"
pdf(file=paste0(Dir,i,".pdf"),width=2.5, height=2.5)
VlnPlot(NK_3, features = c(i), slot = "data",  pt.size = 0, cols=c("#6DA9DC", "#E96653"))+NoLegend()+theme(text = element_text(size = 12))
dev.off()


NK_all <- subset(NK, subset=seurat_clusters!=5)

Idents(NK_all) <- "seurat_clusters"
i<-"AP1G1"
pdf(file=paste0(Dir,i,".cluster.pdf"),width=5, height=2.5)
VlnPlot(NK_all, features = c(i), slot = "data",  pt.size = 0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 12))
dev.off()


NK <- readRDS('/Users/xixi/Project/WY_scRNA/varGenes2000_CCA20_PCA20_Res0.18/NKcombined.rds')
NKK <- subset(NK, subset=seurat_clusters!=5)



gene <- 'PNN'
pdf(file=paste0(Dir,gene,'.pdf'))
FeaturePlot(object = NKK, features = gene)+
scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdYlBu")))
dev.off()
gene <- 'BHLHE40'
pdf(file=paste0(Dir,gene,'.pdf'))
FeaturePlot(object = NKK, features = gene)+
scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdYlBu")))
dev.off()
gene <- 'GOLGA8A'
pdf(file=paste0(Dir,gene,'.pdf'))
FeaturePlot(object = NKK, features = gene)+
scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdYlBu")))
dev.off()
gene <- 'PRRC2C'
pdf(file=paste0(Dir,gene,'.pdf'))
FeaturePlot(object = NKK, features = gene)+
scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdYlBu")))
dev.off()
gene <- 'ANKRD36C'
pdf(file=paste0(Dir,gene,'.pdf'))
FeaturePlot(object = NKK, features = gene)+
scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdYlBu")))
dev.off()



i <- 'XCL1'
pdf(file=paste0(Dir,i,".cluster.pdf"),width=5, height=2.5)
VlnPlot(NKK, features = i, slot = "data",  pt.size = 0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 12))
dev.off()


part1 <- read.table('/Users/xixi/Project/WY_scRNA/varGenes2000_CCA20_PCA20_Res0.18/Files_WY_NKscRNA/未命名文件夹/subset_genes_part5.txt',header=TRUE)
gene_part1 <- part1$gene

pdf(file=paste0(Dir,'VlnPlot_subset_part5.pdf'),width=10, height=8)
VlnPlot(NK, features = gene_part1, pt.size=0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 8))
dev.off()

part1 <- read.table('/Users/xixi/Project/WY_scRNA/varGenes2000_CCA20_PCA20_Res0.18/Files_WY_NKscRNA/未命名文件夹/subset_genes_part6.txt',header=TRUE)
gene_part1 <- part1$gene

pdf(file=paste0(Dir,'VlnPlot_subset_part6.pdf'),width=6, height=5)
VlnPlot(NK, features = gene_part1, pt.size=0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 8))
dev.off()



part1 <- read.table('/Users/xixi/Project/WY_scRNA/varGenes2000_CCA20_PCA20_Res0.18/Files_WY_NKscRNA/未命名文件夹/drug_genes_part4.txt',header=TRUE)
gene_part1 <- part1$gene

pdf(file=paste0(Dir,'VlnPlot_drug_genes_part4.pdf'),width=10, height=9)
VlnPlot(NK_3, features = gene_part1, pt.size=0, cols=c("#6DA9DC", "#E96653"))+NoLegend()+theme(text = element_text(size = 8))
dev.off()


part1 <- read.table('/Users/xixi/Project/WY_scRNA/varGenes2000_CCA20_PCA20_Res0.18/Files_WY_NKscRNA/未命名文件夹/drug_genes_part5.txt',header=TRUE)
gene_part1 <- part1$gene

pdf(file=paste0(Dir,'VlnPlot_drug_ITGAL.pdf'),width=4, height=2)
VlnPlot(NK_3, features = 'ITGAL', pt.size=0, cols=c("#6DA9DC", "#E96653"))+NoLegend()+theme(text = element_text(size = 8))
dev.off()



pdf(file=paste0(Dir,'VlnPlot_KLRB1.pdf'))
VlnPlot(NK, features = 'KLRB1', pt.size=0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 8))
dev.off()

pdf(file=paste0(Dir,'VlnPlot_ITGAL.pdf'))
VlnPlot(NK, features = 'ITGAL', pt.size=0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 8))
dev.off()


Idents(NK_3) <- "stim"

NK_stim <- subset(NK, stim=='DRUG')
NK_ctrl<- subset(NK, stim=='CTRL')
pdf(file=paste0(Dir,'VlnPlot_ITGAL_NK_stim.pdf'))
VlnPlot(NK_stim, features = 'ITGAL', pt.size=0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 8))
dev.off()
pdf(file=paste0(Dir,'VlnPlot_ITGAL_NK_ctrl.pdf'))
VlnPlot(NK_ctrl, features = 'ITGAL', pt.size=0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 8))
dev.off()

Idents(NK) <- "seurat_clusters"
i<-"GNLY"
pdf(file=paste0(Dir,i,".seurat_clusters.pdf"),width=4, height=2)
VlnPlot(NK, features = c(i), slot = "data",  pt.size = 0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 10))
dev.off()

i<-"KIR2DL3"
pdf(file=paste0(Dir,i,".seurat_clusters.pdf"),width=4, height=2)
VlnPlot(NK, features = c(i), slot = "data",  pt.size = 0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 10))
dev.off()

i<-"KIR2DL4"
pdf(file=paste0(Dir,i,".seurat_clusters.pdf"),width=4, height=2)
VlnPlot(NK, features = c(i), slot = "data",  pt.size = 0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 10))
dev.off()
i<-"B3GAT1"
pdf(file=paste0(Dir,i,".seurat_clusters.pdf"),width=4, height=2)
VlnPlot(NK, features = c(i), slot = "data",  pt.size = 0, cols=c("#9DD1AC", "#60B6A9", "#455653", "#EB85A7", "#EFEFA3"))+NoLegend()+theme(text = element_text(size = 10))
dev.off()


markers.to.plot <- c('ATP5F1E','ATP5MC2','ATP5PF','COX5B','COX6A1','COX6C','COX7B','COX7C','ND6','NDUFA4','NDUFAB1','NDUFB6','UQCRB','UQCRH','UQCR11','UQCR10','MPC2','GSTP1','PRDX1','PSMA4','PSMB1','PSMB2','PSME1','TXN','LAMTOR5','ELOB')

markers.to.plot <- c('CD33','TIGIT','SIGLEC9','SIGLEC7','LAG3','KLRB1','KLRC1','KIR2DL1','LILRB1','KIR2DL3','KLRG1','KIR3DL1','KIR3DL2','CD300A')
pdf(file=paste0(Dir,'DotPlot_NKcellinhibitoryreceptor.pdf'),width=8, height=3.2)
DotPlot(NK, features = rev(markers.to.plot), cols = 'RdYlBu', dot.scale = 8) + RotatedAxis()
dev.off()

markers.to.plot <- c('IL17R'
,'IL12RA'
,'TOX'
,'IFNG'
,'CD58'
,'CD38'
,'SELL'
,'KIT'
,'NCAM1'
,'KLRD1'
,'B3GAT1'
,'EGR2'
,'CD7'
,'CD247'
,'ITGAL'
,'ITGAM'
,'EOMES'
,'ITGB2'
,'TYROBP'
,'CD8A'
,'CD69'
,'FOSL2'
,'HCST'
,'NKG7'
,'ELF4'
,'FCGR1G'
,'TBX21')
pdf(file=paste0(Dir,'DotPlot_NKcellmaturation.pdf'),width=8, height=3.2)
DotPlot(NK, features = rev(markers.to.plot), cols = 'RdYlBu', dot.scale = 8) + RotatedAxis()
dev.off()




pdf(file=paste0(Dir,'MarkerGene_Top10_DotPlot.pdf'),width=15, height=2.8)
DotPlot(NK, features = DF$gene, cols = 'RdYlBu', dot.scale = 6) + RotatedAxis()
dev.off()
