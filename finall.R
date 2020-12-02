#coldatafile
datTraits<- read.table(file = "clipboard", sep = "\t", header=TRUE)#for the WGCNA
m<-read.table(file = "clipboard", sep = "\t", header=TRUE)#for characteristics of patients
i<-read.table(file = "clipboard", sep = "\t", header=TRUE)#rnaseq id for each patient

library(dplyr)
m %>%
  filter(m$id %in% i$donor_id)
#bamfiles
filenames <- file.path("D:/OGT", paste0(m$Run, ".bam"))
file.exists(filenames)
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
library("GenomicFeatures")
TxDb <- makeTxDbFromGFF(file="E:/OGT/genes.gtf",format="gtf")
ebg <- exonsBy(TxDb, by = "gene")
library("GenomicAlignments")
library("BiocParallel")
multicoreWorkers()
se <- summarizeOverlaps(features = ebg, reads = bamfiles,
                        mode = "Union",
                        singleEnd = FALSE,
                        ignore.strand = TRUE,
                        fragments = TRUE)
#imputemissingdatausing KNN
library(VIM)
d<-kNN(m,k=3,imp_var = F)
#addingcoldata
dim(se)
colSums(assay(se))
colData(se) <- DataFrame(d)
#referencelevel
library("magrittr")
se$Dx %<>% relevel("Control")
se$Dx
library("DESeq2")
dds <- DESeqDataSet(se, design = ~Dx)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
#chosesamplevariance due to low number
vsd <- vst(dds, blind = T)

#heatmap
sampleDists <- dist(t(assay(vsd)))
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group, vsd$brain, sep = " - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
#checking covariates
#apoe4 and it was repeated for each one
#plotPCA(vsd, intgroup = c("Dx","apoe4"))
plotPCA(rld, intgroup = c("Dx","apoe4"))
library(ggplot2)
pcaData <- plotPCA(rld, intgroup=c("Dx", "apoe4"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Dx, shape=apoe4)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#batcheffect
library("sva")
dds <- DESeq(dds)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Dx, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$Dx, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + Dx
dds <- DESeq(ddssva)
resultsNames(dds)
res<-results(dds, contrast=c("Dx","Alzheimer","Control"))
summary(res)
plotMA(res, ylim=c(-2,2))

library("AnnotationDbi")
library("org.Hs.eg.db")
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = row.names(res),
                     column = "GENENAME",
                     keytype = "ENTREZID",
                     multiVals = "first")
res$symbol2 <- mapIds(org.Hs.eg.db,
                     keys = row.names(res),
                     column = "SYMBOL",
                     keytype = "ENTREZID",
                     multiVals = "first")
resOrdered <- res[order(res$padj),]
write.csv(resOrdered,"res.csv")

#after this step I used the online website for prediction of OGT then made a list of the genes of these proteins
#extract the value of logchange of O-GLcNAcyalted proteins from the results
#extractinglogexpression
keep<-c("100131755", "10040", "10137", "10194", "10410", "10645", "10716", 
        "10725", "10810", "10844", "10902", "10914", "11043", "11309", 
        "11325", "115548", "123591", "124402", "126248", "1415", "144108", 
        "151613", "152485", "154796", "154810", "157285", "1602", "1657", 
        "169792", "170302", "1874", "1875", "195827", "1981", "201191", 
        "2130", "221037", "22846", "23015", "23051", "23063", "23082", 
        "23122", "23139", "23144", "23164", "23177", "23210", "23268", 
        "23316", "23389", "23587", "23613", "250", "2521", "255967", 
        "25890", "26040", "26130", "26165", "26207", "26762", "284069", 
        "29123", "29127", "338773", "340596", "3566", "3654", "375449", 
        "387914", "4001", "4037", "404734", "4109", "4208", "4299", "440193", 
        "441452", "4774", "4798", "4850", "4948", "49854", "5080", "51230", 
        "51366", "51554", "5256", "5360", "53838", "54207", "5430", "54620", 
        "54682", "5469", "54815", "54819", "54823", "54882", "54949", 
        "552", "55578", "55619", "55689", "55691", "55728", "56301", 
        "56776", "56999", "57142", "57459", "57478", "57497", "575", 
        "57532", "57648", "57659", "57794", "5783", "5886", "5937", "5980", 
        "5990", "610", "6314", "63894", "6397", "64207", "6421", "6427", 
        "6431", "64324", "64375", "645369", "645961", "646851", "6664", 
        "6693", "6722", "6787", "6875", "6881", "6885", "6924", "7030", 
        "7090", "7133", "727830", "727905", "7456", "7536", "7545", "79576", 
        "79595", "79680", "79805", "79832", "79872", "79890", "79932", 
        "80129", "80209", "80223", "80315", "80829", "81545", "8208", 
        "8224", "83878", "83992", "84067", "84071", "8408", "84132", 
        "84243", "8464", "84706", "84717", "84722", "84930", "84951", 
        "84952", "8518", "85439", "85440", "85456", "8704", "8848", "90113", 
        "9113", "9149", "92211", "9278", "9369", "9423", "9451", "9609", 
        "9632", "9651", "9730", "9820", "9860", "9883", "9972", "20", 
        "118", "274", "1760", "2034", "3996", "4041", "4302", "4303", 
        "4337", "4661", "5116", "5331", "5394", "5451", "6657", "6909", 
        "7461", "8189", "8242", "8555", "9639", "9644", "9736", "9739", 
        "9794", "9846", "10298", "10454", "10743", "11016", "23037", 
        "23090", "23091", "23152", "23187", "23524", "23654", "25865", 
        "26205", "27146", "28996", "29994", "51564", "54467", "54754", 
        "55023", "55904", "57513", "57597", "57661", "58509", "64760", 
        "65125", "79447", "79778", "80115", "91369", "92154", "152273", 
        "155435", "283987", "441457", "10129", "10178", "10313", "10962", 
        "134957", "146223", "221895", "22866", "2290", "22941", "23035", 
        "23047", "23349", "23426", "2533", "25758", "25849", "26003", 
        "26037", "26115", "283106", "283209", "285513", "2904", "3781", 
        "3797", "388753", "390616", "4130", "4131", "427", "4584", "4747", 
        "5048", "54332", "54602", "55074", "55122", "57062", "57451", 
        "57512", "57583", "57620", "57692", "57700", "644150", "6854", 
        "6857", "7357", "79022", "79589", "79605", "79668", "79719", 
        "79772", "79884", "83593", "83851", "861", "8715", "8825", "8867", 
        "9856", "9892","28996", "23187", "20", "6909", "2034", "57513", "80115", 
        "23037", "58509", "4303", "4041", "274", "441457", "9639", "79778", 
        "9846", "57597", "57661", "23524", "4661", "23152", "51564", 
        "23090", "25865", "23654", "5116", "8555", "283987", "65125", 
        "6657", "26205", "5331", "8189", "9739", "64760", "4337", "10454", 
        "79447", "29994", "155435", "3996", "4302", "54467", "23091", 
        "55904", "5451", "5394", "27146", "55023", "8242", "118", "9736", 
        "54754", "7461", "1760", "9644", "92154", "10298", "91369", "11016", 
        "9794", "152273", "10743","100131755", "10040", "10137", "10194", "10410", 
        "10645", "10716", "10725", "10810", "10844", "10902", "10914", 
        "11043", "11309", "11325", "115548", "123591", "124402", "126248", 
        "1415", "144108", "151613", "152485", "154796", "154810", "157285", 
        "1602", "1657", "169792", "170302", "1874", "1875", "195827", 
        "1981", "201191", "2130", "221037", "22846", "23015", "23051", 
        "23063", "23082", "23122", "23139", "23144", "23164", "23177", 
        "23210", "23268", "23316", "23389", "23587", "23613", "250", 
        "2521", "255967", "25890", "26040", "26130", "26165", "26207", 
        "26762", "284069", "29123", "29127", "338773", "340596", "3566", 
        "3654", "375449", "387914", "4001", "4037", "404734", "4109", 
        "4208", "4299", "440193", "441452", "4774", "4798", "4850", "4948", 
        "49854", "5080", "51230", "51366", "51554", "5256", "5360", "53838", 
        "54207", "5430", "54620", "54682", "5469", "54815", "54819", 
        "54823", "54882", "54949", "552", "55578", "55619", "55689", 
        "55691", "55728", "56301", "56776", "56999", "57142", "57459", 
        "57478", "57497", "575", "57532", "57648", "57659", "57794", 
        "5783", "5886", "5937", "5980", "5990", "610", "6314", "63894", 
        "6397", "64207", "6421", "6427", "6431", "64324", "64375", "645369", 
        "645961", "646851", "6664", "6693", "6722", "6787", "6875", "6881", 
        "6885", "6924", "7030", "7090", "7133", "727830", "727905", "7456", 
        "7536", "7545", "79576", "79595", "79680", "79805", "79832", 
        "79872", "79890", "79932", "80129", "80209", "80223", "80315", 
        "80829", "81545", "8208", "8224", "83878", "83992", "84067", 
        "84071", "8408", "84132", "84243", "8464", "84706", "84717", 
        "84722", "84930", "84951", "84952", "8518", "85439", "85440", 
        "85456", "8704", "8848", "90113", "9113", "9149", "92211", "9278", 
        "9369", "9423", "9451", "9609", "9632", "9651", "9730", "9820", 
        "9860", "9883", "9972")
#9892
rr<-res[res$X1 %in% keep, ]
write.csv(rr,"ogtfnal.csv")
#gene_enrichment_analysis
library(clusterProfiler) #Does not include background genes so I need to know another one

#Molecular enrichment
sample_test <- enrichGO(rr$X1, pvalueCutoff=1, qvalueCutoff=1,'org.Hs.eg.db', ont="MF")
dotplot(sample_test, showCategory=30)
#Celular component
sample_test <- enrichGO(rr$X1, pvalueCutoff=1, qvalueCutoff=1,'org.Hs.eg.db', ont="CC")
dotplot(sample_test, showCategory=30)
#Biologic component
sample_test <- enrichGO(rr$X1, pvalueCutoff=1, qvalueCutoff=1,'org.Hs.eg.db', ont="BP")
dotplot(sample_test, showCategory=30)
#PPathwayenrichment
library(ReactomePA)
x <- enrichPathway(gene=rr$X1,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))
dotplot(x,showCategory=15)
#Diseaseenrichm
library(DOSE)
edo <- enrichDGN(down$E)
p1 <- emapplot(edo)
p2 <- emapplot(edo, pie_scale=1.5)
p3 <- emapplot(edo,layout="kk")
p4 <- emapplot(edo, pie_scale=1.5,layout="kk") 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
#####################################
#WCGNA
#extractOGTdatafromse
#openning the unormalized data(I could not create matrix from the se variable so I used FPKM)
fpkm_table_unnormalized <- read_csv("E://OGT/fpkm_table_unnormalized.csv")
#export gene names for annotation of well ID
rows_genes <- read_csv("E://OGT/rows-genes.csv")
#adding gene symbol to the data
dat<-cbind(rows_genes$gene_entrez_id,fpkm_table_unnormalized)
colnames(dff)<-dff[1,]
dff <- dff[-1, ]
colnames(dat)[2] = "wellid"
colnames(dat)[1] = "gene"
#remove the well id
dat$wellid<-NULL
#remove unneeded groups
tt<-t(dta_m)
tt<-as.data.frame(tt)
dput(names(dat))
#keep alzeheimer groups and control
t<-c("gene", "488395315", "496100277", "496100278", "496100279", 
     "496100283", "496100287", "496100288", "496100290", "496100294", 
     "496100295", "496100296", "496100297", "496100301", "496100303", 
     "496100304", "496100305", "496100306", "496100307", "496100310", 
     "496100313", "496100314", "496100315", "496100316", "496100317", 
     "496100318", "496100319", "496100320", "496100323", "496100324", 
     "496100327", "496100328", "496100329", "496100330", "496100332", 
     "496100333", "496100335", "496100336", "496100337", "496100338", 
     "496100339", "496100340", "496100341", "496100342", "496100343", 
     "496100344", "496100345", "496100346", "496100347", "496100349", 
     "496100351", "496100352", "496100353", "496100354", "496100355", 
     "496100356", "496100360", "496100361", "496100364", "496100365", 
     "496100366", "496100367", "496100368", "496100369", "496100370", 
     "496100371", "496100372", "496100373", "496100374", "496100375", 
     "496100376", "496100377", "496100378", "496100379", "496100381", 
     "496100382", "496100383", "496100384", "496100392", "496100393", 
     "496100395", "496100397", "496100400", "496100401", "496100402", 
     "496100416", "496100422")
df<-dat[t]
dff<-t(df)
#rownames
dff<-as.data.frame(dff)
dff <- data.frame(dff, row.names = 1)
#remove X in columnnames :(
colnames(dff)<-gsub("X","",colnames(dff))
gggg<-t(as.dataframe(dff))
gggg<-as.data.frame(gggg)
#keep genes that has O-glycNac signal
keep<-c("100131755", "10040", "10137", "10194", "10410", "10645", "10716", 
        "10725", "10810", "10844", "10902", "10914", "11043", "11309", 
        "11325", "115548", "123591", "124402", "126248", "1415", "144108", 
        "151613", "152485", "154796", "154810", "157285", "1602", "1657", 
        "169792", "170302", "1874", "1875", "195827", "1981", "201191", 
        "2130", "221037", "22846", "23015", "23051", "23063", "23082", 
        "23122", "23139", "23144", "23164", "23177", "23210", "23268", 
        "23316", "23389", "23587", "23613", "250", "2521", "255967", 
        "25890", "26040", "26130", "26165", "26207", "26762", "284069", 
        "29123", "29127", "338773", "340596", "3566", "3654", "375449", 
        "387914", "4001", "4037", "404734", "4109", "4208", "4299", "440193", 
        "441452", "4774", "4798", "4850", "4948", "49854", "5080", "51230", 
        "51366", "51554", "5256", "5360", "53838", "54207", "5430", "54620", 
        "54682", "5469", "54815", "54819", "54823", "54882", "54949", 
        "552", "55578", "55619", "55689", "55691", "55728", "56301", 
        "56776", "56999", "57142", "57459", "57478", "57497", "575", 
        "57532", "57648", "57659", "57794", "5783", "5886", "5937", "5980", 
        "5990", "610", "6314", "63894", "6397", "64207", "6421", "6427", 
        "6431", "64324", "64375", "645369", "645961", "646851", "6664", 
        "6693", "6722", "6787", "6875", "6881", "6885", "6924", "7030", 
        "7090", "7133", "727830", "727905", "7456", "7536", "7545", "79576", 
        "79595", "79680", "79805", "79832", "79872", "79890", "79932", 
        "80129", "80209", "80223", "80315", "80829", "81545", "8208", 
        "8224", "83878", "83992", "84067", "84071", "8408", "84132", 
        "84243", "8464", "84706", "84717", "84722", "84930", "84951", 
        "84952", "8518", "85439", "85440", "85456", "8704", "8848", "90113", 
        "9113", "9149", "92211", "9278", "9369", "9423", "9451", "9609", 
        "9632", "9651", "9730", "9820", "9860", "9883", "9972", "20", 
        "118", "274", "1760", "2034", "3996", "4041", "4302", "4303", 
        "4337", "4661", "5116", "5331", "5394", "5451", "6657", "6909", 
        "7461", "8189", "8242", "8555", "9639", "9644", "9736", "9739", 
        "9794", "9846", "10298", "10454", "10743", "11016", "23037", 
        "23090", "23091", "23152", "23187", "23524", "23654", "25865", 
        "26205", "27146", "28996", "29994", "51564", "54467", "54754", 
        "55023", "55904", "57513", "57597", "57661", "58509", "64760", 
        "65125", "79447", "79778", "80115", "91369", "92154", "152273", 
        "155435", "283987", "441457", "10129", "10178", "10313", "10962", 
        "134957", "146223", "221895", "22866", "2290", "22941", "23035", 
        "23047", "23349", "23426", "2533", "25758", "25849", "26003", 
        "26037", "26115", "283106", "283209", "285513", "2904", "3781", 
        "3797", "388753", "390616", "4130", "4131", "427", "4584", "4747", 
        "5048", "54332", "54602", "55074", "55122", "57062", "57451", 
        "57512", "57583", "57620", "57692", "57700", "644150", "6854", 
        "6857", "7357", "79022", "79589", "79605", "79668", "79719", 
        "79772", "79884", "83593", "83851", "861", "8715", "8825", "8867", 
        "9856", "9892","28996", "23187", "20", "6909", "2034", "57513", "80115", 
        "23037", "58509", "4303", "4041", "274", "441457", "9639", "79778", 
        "9846", "57597", "57661", "23524", "4661", "23152", "51564", 
        "23090", "25865", "23654", "5116", "8555", "283987", "65125", 
        "6657", "26205", "5331", "8189", "9739", "64760", "4337", "10454", 
        "79447", "29994", "155435", "3996", "4302", "54467", "23091", 
        "55904", "5451", "5394", "27146", "55023", "8242", "118", "9736", 
        "54754", "7461", "1760", "9644", "92154", "10298", "91369", "11016", 
        "9794", "152273", "10743","100131755", "10040", "10137", "10194", "10410", 
        "10645", "10716", "10725", "10810", "10844", "10902", "10914", 
        "11043", "11309", "11325", "115548", "123591", "124402", "126248", 
        "1415", "144108", "151613", "152485", "154796", "154810", "157285", 
        "1602", "1657", "169792", "170302", "1874", "1875", "195827", 
        "1981", "201191", "2130", "221037", "22846", "23015", "23051", 
        "23063", "23082", "23122", "23139", "23144", "23164", "23177", 
        "23210", "23268", "23316", "23389", "23587", "23613", "250", 
        "2521", "255967", "25890", "26040", "26130", "26165", "26207", 
        "26762", "284069", "29123", "29127", "338773", "340596", "3566", 
        "3654", "375449", "387914", "4001", "4037", "404734", "4109", 
        "4208", "4299", "440193", "441452", "4774", "4798", "4850", "4948", 
        "49854", "5080", "51230", "51366", "51554", "5256", "5360", "53838", 
        "54207", "5430", "54620", "54682", "5469", "54815", "54819", 
        "54823", "54882", "54949", "552", "55578", "55619", "55689", 
        "55691", "55728", "56301", "56776", "56999", "57142", "57459", 
        "57478", "57497", "575", "57532", "57648", "57659", "57794", 
        "5783", "5886", "5937", "5980", "5990", "610", "6314", "63894", 
        "6397", "64207", "6421", "6427", "6431", "64324", "64375", "645369", 
        "645961", "646851", "6664", "6693", "6722", "6787", "6875", "6881", 
        "6885", "6924", "7030", "7090", "7133", "727830", "727905", "7456", 
        "7536", "7545", "79576", "79595", "79680", "79805", "79832", 
        "79872", "79890", "79932", "80129", "80209", "80223", "80315", 
        "80829", "81545", "8208", "8224", "83878", "83992", "84067", 
        "84071", "8408", "84132", "84243", "8464", "84706", "84717", 
        "84722", "84930", "84951", "84952", "8518", "85439", "85440", 
        "85456", "8704", "8848", "90113", "9113", "9149", "92211", "9278", 
        "9369", "9423", "9451", "9609", "9632", "9651", "9730", "9820", 
        "9860", "9883", "9972")
datExpr<-dff[keep]
tttt<-as.data.frame(dat)
setDT(tttt, keep.rownames = TRUE)[]
#WCGNA analysis
library(WGCNA)
library(Biobase)
library(GEOquery)
library(limma)
library(MASS) # standard, no need to install
library(class) # standard, no need to install
library(cluster)
library(impute)# install it for imputing missing value
library(scatterplot3d) 
library(WGCNA)
datExpr<-t(rr)
colnames(datExpr) <- datExpr[1,]
datExpr <- datExpr[-1, ]
enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType ="signed")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),ylim=c(0,1));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
net = blockwiseModules(datExpr, power = 12,
                      networkType = "signed", minModuleSize = 10,
                      mergeCutHeight = 0.1,
                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                      saveTOMs = FALSE
                      ,deepSplit = 2,
                      saveTOMFileBase = "OGT",
                      verbose = 3, maxBlockSize = 100)

table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 3), "\n(",
                   signif(moduleTraitPvalue, 3), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = T,
               colors = blueWhiteRed(50),
               
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"),textMatrix = textMatrix)

#Choosing one module and Tau PROTEIN
# Define variable weight containing the weight column of datTrait
tau = as.data.frame(datTraits$ptau);
names(tau) = "Tau"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, tau, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(AS`), sep="")
names(GSPvalue) = paste("p.GS.", names(AS), sep="");
###################
module = "midnightblue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for tau2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

########################################################
#dendrogram
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, tau))
# Plot the relationships among the eigengenes and the trait
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
#######################
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 12);
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select module
module = "midnightblue";
# Select module probes
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,probeToGene = data.frame(rr$symbol, rr$symbol2) )
