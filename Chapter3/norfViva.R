## package installation
# Install and Load the BiocManager package
install.packages("BiocManager")
library(BiocManager)

install.packages("protr")

install.packages('svglite')

# Install the kebabs package
BiocManager::install("kebabs")

# Install the stringi package
install.packages("stringi")

# Install the Annotation Hub package
BiocManager::install("AnnotationHub")

# Install the protr package
BiocManager::install("protr")

# Install the ballgown package
BiocManager::install("ballgown")

# Install the ggplot2 package
install.packages("ggplot2")

# Install the corrplot package
install.packages("corrplot")

# Install the ggExtra package
install.packages("ggExtra")

# Install the plyr package
install.packages("plyr")

# Install the ggpubr package
install.packages("ggpubr")

# Install the caret package
install.packages("caret")

# Install the gridExtra package
install.packages("gridExtra")

# Install the pROC package
install.packages("pROC")

# Install the xgboost package
install.packages("xgboost")

# Install the data.table package
install.packages("data.table")

install.packages("drat", repos="https://cran.rstudio.com")
drat:::addRepo("dmlc")
install.packages("xgboost")
install.packages("phylotools")

## Package initialisation 
library(kebabs)
library(stringi)
library("protr")
library(ballgown)
library(ggplot2)
library(corrplot)
library(ggplot2)
library(ggExtra)
library(plyr)
library(ggpubr)
library(caret)
library(gridExtra)
library(pROC)
require(xgboost)
library(data.table)
library("Biostrings")
library(phylotools)



## data loading
# norfsDB 1.1
norfsDB = gffRead("nORFsDB.1.1.gtf")

norfsAttributesAAseq = getAttributeField(norfsDB$attributes, "AA_seq", attrsep = "; ")
norfsAttributesGeneId = getAttributeField(norfsDB$attributes, "gene_id", attrsep = "; ")
norfsAttributesStart = getAttributeField(norfsDB$attributes, "start_codon", attrsep = "; ")

nORFAAseq = gsub("[^[:alnum:]=\\.]", "", norfsAttributesAAseq)
nORFGeneIds = gsub("[^[:alnum:]=\\.]", "", norfsAttributesGeneId)
nORFStart = gsub("[^[:alnum:]=\\.]", "", norfsAttributesStart)

#
nORFGeneIds_missing = nORFGeneIds[which(nORFGeneIds %in% not_overlapping)]
nORFAAseq_missing = nORFAAseq[which(nORFGeneIds %in% not_overlapping)]
write.csv(nORFGeneIds_missing, "failed_esm_sequences.csv")


write.table(paste0(">",nORFGeneIds_missing,"\n",nORFAAseq_missing), file = "nORFs_missing.fasta", sep = "", row.names = F, col.names = F, quote = F)


#create .fasta for molecular dynamics
write.table(paste0(">",nORFGeneIds,"\n",nORFAAseq), file = "nORFs.fasta", sep = "", row.names = F, col.names = F, quote = F)

# random sequences
randomSeq = read.csv("1000RandomSequences.csv", header = TRUE)
randomSeq = randomSeq$aa

#check if all pdb files are converted
files = list.files(path = './pdb', all.files = TRUE, recursive = TRUE)
files <- sub(".pdb", "", files)
not_overlapping <- setdiff(nORFGeneIds,files)
print(not_overlapping)


# routine to generate random sequences instead
# seqLength <- runif(1000, min=50, max=1000)
# aa = genRandBioSeqs(seqType = "AA", 1000, seqLength,
#                    biostring = TRUE, seed=11)

# canonical ORFs 
# https://www.uniprot.org/proteomes/UP000005640
fastaFile <- readDNAStringSet("UP000005640_9606_20k.fasta")

seq_name = names(fastaFile)
df <- data.frame(seq_name, sequence)
head(df$seq_name)
cORFs = df$sequence


# Load the ggplot2 package
library(ggplot2)

# Count the number of occurrences of each start codon
codon_counts <- table(norfsAttributesStart)

# Create a data frame with the start codons and their counts
codon_data <- data.frame(codon=names(codon_counts), count=codon_counts)

# Create a bar plot of the codon counts
ggplot(codon_data, aes(x=codon, y=count)) +
  geom_col() +
  labs(x="Start Codon", y="Count") +
  theme(axis.text.x=element_text(angle=90, hjust=1))
## data processing 
brotoFunc <- function (x) {
  return(tryCatch(extractMoreauBroto(x), error=function(e) NULL))
}

triads <- function (x) {
  return(tryCatch(extractCTriad(x), error=function(e) NULL))
}

fAAC <- function (x) {
  return(tryCatch(extractAAC(x), error=function(e) NULL))
}

fDC <- function (x) {
  return(tryCatch(extractDC(x), error=function(e) NULL))
}

fCTDC <- function (x) {
  return(tryCatch(extractCTDC(x), error=function(e) NULL))
}

fCTDT <- function (x) {
  return(tryCatch(extractCTDT(x), error=function(e) NULL))
}

fCTDD <- function (x) {
  return(tryCatch(extractCTDD(x), error=function(e) NULL))
}

fQSO <- function (x) {
  return(tryCatch(extractSOCN(x), error=function(e) NULL))
}

fQSO <- function (x) {
  return(tryCatch(extractQSO(x), error=function(e) NULL))
}

fPAAC <- function (x) {
  return(tryCatch(extractPAAC(x), error=function(e) NULL))
}

fAPAAC <- function (x) {
  return(tryCatch(extractAPAAC(x), error=function(e) NULL))
}

fSOCN <- function (x) {
  return(tryCatch(extractSOCN(x), error=function(e) NULL))
}

##conversion to data frames
conversionFrame = function(x){
  a = data.frame(stri_list2matrix(lapply(x, unlist), byrow=TRUE, fill=NA))
  colnames(a) = names(x[[1]])
  return(a)
}

## analysis
pchemAnalysis = function(seq, type) { 

  mBroto = lapply(seq, brotoFunc)
  
  #extractCTriad
  mTriads = lapply(seq, triads)
  
  #AA composition
  mAAC = lapply(seq, fAAC)
  mDC = lapply(seq, fDC)
  
  #CTD
  mCTDC = lapply(seq, fCTDC)
  mCTDT = lapply(seq, fCTDT)
  mCTDD = lapply(seq, fCTDD)
  
  
  #quasi sequence order
  mSOCN = lapply(seq, fSOCN)
  mQSO = lapply(seq, fQSO)
  
  #pseudo amino acid composition
  mPAAC = lapply(seq, fPAAC)
  mAPAAC = lapply(seq, fAPAAC)
  
  #convert to dataframes
  Fmbroto = conversionFrame(mBroto)
  Fmtriads = conversionFrame(mTriads)
  FmAAC = conversionFrame(mAAC)
  FmDC = conversionFrame(mDC)
  FmCTDC = conversionFrame(mCTDC)
  FmCTDT = conversionFrame(mCTDT)
  FmCTDD = conversionFrame(mCTDD)
  FmSOCN = conversionFrame(mSOCN)
  FmQSO = conversionFrame(mQSO)
  FmPAAC = conversionFrame(mPAAC)
  FmAPAAC = conversionFrame(mAPAAC)
  
  #aggregate
  jointFrame = data.frame(type, seq, Fmbroto, Fmtriads, FmAAC, FmDC, FmCTDC, FmCTDT, FmCTDD, FmSOCN, FmQSO, FmPAAC, FmAPAAC)
  return(jointFrame)
}

## calculate

#nORFs
nORF = pchemAnalysis(nORFAAseq, "nORF")
nORF$id = nORFGeneIds
write.csv(nORF, "nORFs_w_id.csv")

# canonical
ORF = pchemAnalysis(cORFs, "ORF")
ORF$id = seq_name
ORF = na.omit(ORF)
write.csv(ORF, "ORFs_w_id.csv")

#random
Random = pchemAnalysis(randomSeq, "Random")
write.csv(Random, "Random_w_id.csv")

#Restore calculations
nORF = read.csv("nORFs_w_id.csv", header = TRUE)
ORF = read.csv("ORFs_w_id.csv", header = TRUE)
Random = read.csv("Random_w_id.csv", header = TRUE)

Random$X.1 = NULL
ORF$X.1 = NULL
nORF$X.1 = NULL

ORF$X = NULL
nORF$X = NULL
Random$X = NULL

nORFpure = nORF
nORFpure$id = NULL
nORFpure$type = NULL
nORFpure$X = NULL
nORFpure$seq = NULL
nORFpure$X.1 = NULL
non_numeric <- sapply(nORFpure, is.numeric)
non_numeric_cols <- names(which(!non_numeric))
nORFpureScale = scale(nORFpure)
nORFpureCor = cor(nORFpure)

#append Id seqs 
Random$id = "random"

#store current session
save(r, file="nORFData.RData")

#joinframe creation
jointFrame = rbind(nORF, ORF, Random)
corrFrame = rbind(nORF, ORF)

#get length of sequences
lengths = sapply(corrFrame$seq, nchar)

#prepare for tensorboard and corrplot
jointFrame_noNA = na.omit(corrFrame);
jointFrame_noNA$id = NULL
jointFrame_noNA$type = NULL
jointFrame_noNA$X = NULL
jointFrame_noNA$seq = NULL
jointFrame_noNA$X.1 = NULL
jointFrame_noNA = scale(jointFrame_noNA)
jointFrame_noNA$lengths = lengths

#corrplot
corrTable = cor(nORFpureScaleDF, method = "pearson" )


corrTableDF = as.data.frame(corrTable)

#identify those with > 0.3 correlation
names(jointFrame_noNA[(which(abs(corrTableDF$lengths) > 0.29))])
corNames = names(nORFpureScaleDF[(which(abs(corrTableDF$lengths) > 0.2139))])
#remove length
corNames = head(corNames, -1)

#remove those with > 0.3 correlation
nORFpureScaleDF[,corNames] = NULL

#plot groups
corrplot(corrTable, method = "color", order = "hclust", tl.cex = 0.7, addrect = 2)


#phastcons
BiocManager::install("phastCons100way.UCSC.hg38")
phast <- phastCons100way.UCSC.hg19



# load jointframe
dataset = na.omit(jointFrame)
removeRandom = dataset$type == "Random" 
dataset = dataset[!removeRandom,]
dataset$type[dataset$type == "nORF"] = 1
dataset$type[dataset$type == "ORF"] = 0
dataset$id = NULL
dataset$seq = NULL
dataset$length = NULL
dataset$X = NULL
dataset$X.1 = NULL

#equalising the amount of each group
#dataset = rbind(head(dataset[dataset$group == 1,], 10000), head(dataset[dataset$group == 0,], 10000))

#dataset <- as.data.frame(scale(dataset))

## 75% of the sample size
smp_size <- floor(0.70 * nrow(dataset))

## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(dataset)), size = smp_size)

train <- dataset[train_ind, ]
test <- dataset[-train_ind, ]

## create 2 test set (hidden and tuning test set)

## 75% of the sample size
smp_test_size <- floor(0.50 * nrow(test))

## set the seed to make your partition reproducible
test_ind <- sample(seq_len(nrow(test)), size = smp_test_size)

test_tuning <- test[test_ind, ]
test_hidden <- test[-test_ind, ]


trainLabel = train$type
trainLabel = as.numeric(trainLabel)
testTuningLabel = test_tuning$type
testHiddenLabel = test_hidden$type

train$type = NULL
test_tuning$type = NULL
test_hidden$type = NULL

train = as.matrix(train)
train = as.numeric(train)
test = as.matrix(test_tuning)
test_hidden = as.matrix(test_hidden)

bstSparse <- xgboost(data = train, label = trainLabel, max.depth = 5, eta = 0.025, gamma = 0.5, alpha = 3.5, lambda = 2.5, nrounds = 500, verbose = FALSE,
                    objective = 'binary:logistic', eval_metric = 'rmse', prediction = T )

bstSparse <- xgboost(data = train, label = trainLabel, max.depth = 10, eta = 0.1, gamma = 7.5, alpha = 2, lambda = 0.4, nrounds = 250, verbose = FALSE,
                     objective = 'reg:linear', eval_metric = 'auc', prediction = T )

bstSparse2 <- xgboost(data = train, label = trainLabel, nrounds = 100,  eval_metric = 'auc', prediction = T )
bstSparse = bstSparse2

importance_matrix <- xgb.importance(model = bstSparse,  )
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix,  top_n = 50)


it = which.max(bstSparse$evaluation_log$test_auc_mean)
best.iter = bstSparse$evaluation_log$iter[it]

preds <- predict(bstSparse,test)

prediction <- as.numeric(preds > 0.5)

err <- mean(as.numeric(preds > 0.5) != testTuningLabel)
print(paste("test-error=", err))


preds_hidden <- predict(bstSparse,test_hidden)

prediction <- as.numeric(preds_hidden > 0.5)

err <- mean(as.numeric(preds_hidden > 0.5) != testHiddenLabel)
print(paste("test-error=", err))




#install.packages("pROC")
plot(pROC::roc(response = testTuningLabel,
               predictor = prediction,
               levels=c(0, 1), ),
     lwd=3) 

#install.packages("caret")

confusionMatrix(as.factor(prediction), as.factor(testTuningLabel))


tune_grid <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.025, 0.05, 0.15, 0.3, 0.6),
  max_depth = c(2, 5, 10),
  gamma = c(0, 1, 3, 6, 10),
  min_child_weight = 1,
  colsample_bytree = 1,
  subsample = 1
)

tune_grid <- expand.grid(
  nrounds = c(500, 1000),
  eta = c(0.025),
  max_depth = c(2),
  gamma = c(0),
  min_child_weight = 1,
  colsample_bytree = 1,
  subsample = 1
)

tune_control <- caret::trainControl(
  method = "cv", # cross-validation
  number = 3, # with n folds 
  #index = createFolds(tr_treated$Id_clean), # fix the folds
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

xgb_tune <- caret::train(
  x = train,
  y = trainLabel,
  trControl = tune_control,
  tuneGrid = tune_grid,
  method = "xgbTree",
  
  verbose = TRUE
)

# helper function for the plots
tuneplot <- function(x, probs = .90) {
  ggplot(x) +
    theme_bw()
}

tuneplot(xgb_tune)

## top feature seperation plot 
transparentTheme(trans = .9)
featurePlot(x = iris[, 1:4], 
            y = iris$Species,
            plot = "density", 
            ## Pass in options to xyplot() to 
            ## make it prettier
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")), 
            adjust = 1.5, 
            pch = "|", 
            layout = c(4, 1), 
            auto.key = list(columns = 3))

featurePlot(x = iris[, 1:4], 
            y = iris$Species, 
            plot = "box", 
            ## Pass in options to bwplot() 
            scales = list(y = list(relation="free"),
                          x = list(rot = 90)),  
            layout = c(4,1 ), 
            auto.key = list(columns = 2))

featurePlot(tuningResults)



tuningResults$colsample_bytree = NULL
tuningResults$min_child_weight = NULL
tuningResults$subsample = NULL

corrTable = cor(tuningResults)
corrplot(corrTable, method = "color", order = "hclust", tl.cex = 0.7, addrect = 2, outline = TRUE)

## histogram for seperation
histogram = ggplot(a, aes(x = V535)) +
  geom_histogram(aes(color = group, fill = group), 
                 position = "stack", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#D22B2B")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#D22B2B"))

densitogram = ggplot(a, aes(x = V1027)) +
  geom_density(aes(color = group, fill = group), 
                 position = "stack", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#D22B2B")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#D22B2B"))

densitogram = ggplot(a, aes(x = V1027)) +
  geom_density(aes(color = group, fill = group), 
               position = "stack", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#D22B2B")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#D22B2B"))

plot1 <- histogram
plot2 <- densitogram
grid.arrange(plot1, plot2, ncol=2)


p <- ggplot(a, aes(x=V1027, y=V1099, color=group, fill=group)) + geom_point(alpha = 0.01)

ggMarginal(p, type="histogram")

## corrplot of features
featureMatrix = head(importance_matrix,100)
b = a[featureMatrix$Feature]
b$group = a$group 
b$group = NULL
b$length = a$length

corrTable = cor(b)
corrplot(corrTable, method = "square", order = "hclust", tl.cex = 0.7, addrect = 2, outline = TRUE)


corrplot(corrTable, 
         order = "hclust", tl.cex = 0.7, method="color", diag = FALSE, type = "lower",
         addrect = 4, outline = TRUE, 
          tl.col = "black", pch.col = 'grey20')

# feature importance matrix
sortedMatrix = importance_matrix
sortedMatrix$number = gsub("V", "",sortedMatrix$Feature)
sortedMatrix$number = as.integer(sortedMatrix$number)
sortedMatrix = sortedMatrix[order(sortedMatrix$number)]

mid<-mean(log(sortedMatrix$Frequency))

p <- ggplot(sortedMatrix, aes(x=number, y=log(Gain), color=log(Frequency))) + 
  geom_point(alpha = 0.6) + scale_color_gradient2(midpoint=mid, low="black", mid="blue",
                                                  high="red", space ="Lab" )
p

## get all columns with correlation over 0.4 of the top 100 features
featureMatrix = head(importance_matrix,100)
b = a[featureMatrix$Feature]
corrTable = cor(b)
maxN <- function(x, N=2){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}

res3 = apply(filterCorr, 2, function(x) maxN(x))
length(res3[res3 > 0.4])

## box plot comparison

p <- ggboxplot(a, x = "group", y = "V1027",
               color = "group", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")

compare_means(V1027 ~ group,  data = a)

ggboxplot(a, x = "group", y = "V1027",
          color = "group", palette = "jco")+
  stat_compare_means(method = "anova", label.y = 3.5)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "random")  


## plotting


plotTripletPlot = function(dataset, dimension) {
  
  sigHight1_nORF_ORF = max(dataset[dimension][dataset["group"] == "nORF" | dataset["group"] == "ORF"]) + 0.5
  sigHight2_nORF_Random = max(dataset[dimension][dataset["group"] == "nORF" | dataset["group"] == "Random"]) + 1.5
  sigHight3_ORF_Random = max(dataset[dimension][dataset["group"] == "ORF" | dataset["group"] == "Random"]) + 0.1
  

  scaleFUN <- function(x) sprintf("%.2f", x)
  scaleMin = round_any(c(min(dataset[dimension])-0.3) - 0.5, 0.5)
  maxScale = round_any(max(c(sigHight1_nORF_ORF, sigHight2_nORF_Random, sigHight3_ORF_Random)) + 2, 0.5)
  
  my_comparisons <- list(c("nORF", "ORF"), c("nORF", "Random"), c("ORF", "Random"))
  plot2 = ggboxplot(dataset, x = "group", y = dimension,  xlab=" ", size = 0.01,
            color = "group")+ scale_color_manual(values = c("#00AFBB", "#E7B800", "#D22B2B")) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#D22B2B")) + theme(axis.text.y = element_text(size = 10)) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.y = c(sigHight1_nORF_ORF, sigHight2_nORF_Random, sigHight3_ORF_Random), method = "t.test", ) +
    scale_y_continuous(limits = c(scaleMin, maxScale), 
              breaks = seq(scaleMin, maxScale, b=1), expand = c(0,0))
  
  plot1 = densitogram = ggplot(dataset, aes_string(x = dimension)) +
    geom_density(aes(color = group), 
                 position = "identity", alpha = 0.8, size=1) + scale_y_continuous(labels=scaleFUN) +
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#D22B2B")) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#D22B2B")) + theme(legend.position = "none", text = element_text(size = 15), axis.text.x=element_blank(),
                                                                           axis.title.x=element_blank())
  
  plot3 = densitogram = ggplot(dataset, aes_string(x = dimension)) +
    geom_density(aes(color = group, fill = group), 
                 position = "fill", alpha = 0.5) +
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#D22B2B")) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#D22B2B")) + 
    theme(legend.position = "none", text = element_text(size = 15)) + ylab("relative density")
  

  figure = grid.arrange(arrangeGrob(plot1, plot3), plot2, ncol=2, widths = c(2,1))
  #ggsave(file="nORFPlot4.svg", plot=figure, width=1200, height=400, units = "px", device = "svg", scale = TRUE, dpi = 120, )
  figure
}


plotTripletPlot(dataset2, "V1093")
plotTripletPlot(dataset2, "V532")
plotTripletPlot(dataset2, "V240")
plotTripletPlot(dataset2, "V534")
plotTripletPlot(dataset2, "V1122")
plotTripletPlot(dataset2, "V1119")
plotTripletPlot(dataset2, "V1007")
plotTripletPlot(dataset2, "V1011")
plotTripletPlot(dataset2, "V944")
plotTripletPlot(dataset2, "V1027")
plotTripletPlot(dataset2, "V702")

ggplot(a, aes(x="length", y="v240") ) +
  geom_density_2d()


#python: python esm/scripts/extract.py esm2_t36_3B_UR50D nORFs.fasta esmdata/ --repr_layers 0 35 36 --include mean per_tok

duplicated_columns <- duplicated(sapply(cNORF, digest))
colnames(my_df[duplicated_columns])


#pdb analysis
library(bio3d)
pdb = read.pdb("pdb/00ahH1.pdb")
mode = nma(pdb)
IFlex <- function(weights) {
       N <- length(weights)
       numerator <- sum(weights^2)
       denominator <- N
       IFlex <- sqrt(numerator/denominator)
       return(IFlex)
}



IFlex(mode$fluctuations)

#generate trajectory file
mktrj(modes2, file="test.pdb")

allPDB = list.files("pdb")
allNMA = list() 
for(i in 1:length(allPDB)){
  pdb = read.pdb(paste0("pdb/",allPDB[i]))
  mode = nma(pdb)
  print(mode$frequencies)
  allNMA[[allPDB[i]]] = mode
}

for(i in 82802:length(allPDB)){
  pdb = read.pdb(paste0("pdb/",allPDB[i]))
  mode = nma(pdb)
  print(mode$frequencies)
  allNMA[[allPDB[i]]] = mode
}

for(i in 92216:length(allPDB)){
  pdb = read.pdb(paste0("pdb/",allPDB[i]))
  tryCatch({
    mode = nma(pdb)
    print(mode$frequencies)
    allNMA[[allPDB[i]]] = mode
  }, error = function(e) {
    print(paste0("Error processing ", allPDB[i], ": ", e))
  })
}

library(bio3d)
library(tidyverse)

# Set the path to the folder containing PDB files
path <- "dssp"

# Create an empty dataframe to store the information
df <- data.frame(pdb_name = character(), sol_area = numeric(), sol_sd_area = numeric(), phi = numeric(), psi = numeric(), turn = numeric(), helix = numeric(), beta = numeric(), G = numeric())

# Loop through all DSSP files in the folder
for (file in list.files(path, pattern = "*.dssp")) {
  tryCatch({
  # Read in the DSSP file
  dssp <- read.dssp(file.path(path, file))
  
  # Extract the solvable area
  sol_area <- mean(dssp$acc)
  sol_sd_area <- sd(dssp$acc)
  
  # Extract the Ramachandran angles
  phi <- mean(dssp$phi)
  psi <- mean(dssp$psi)
  sse <- dssp$sse
  
  t_count <- sum(sse == "T")
  s_count <- sum(sse == "S")
  h_count <- sum(sse == "H")
  g_count <- sum(sse == "G")
  non_count <- sum(sse == " ")
  
  # Calculate the total number of SSEs
  total_count <- length(sse)
  
  # Calculate the relative proportion of each SSE
  t <- t_count / total_count
  s <- s_count / total_count
  h <- h_count / total_count
  g <- g_count / total_count
  non <- non_count / total_count
  
  # Extract the PDB name without the .dssp extension
  pdb_name <- gsub(".dssp", "", file)
  
  # Append the information to the dataframe
  df <- rbind(df, data.frame(id = pdb_name, sol_area = sol_area, sol_sd_area = sol_sd_area,  phi = phi, psi = psi, turn = t, helix = h, beta = s, G = g  ))
  }, error = function(e) {
    print(paste("Error in processing file:", file))
  })
}

#create sampled ramachandran plot
dsspFiles = list.files(path, pattern = "*.dssp")
dssFilesReduced = dsspFiles
length(dssFilesReduced)


dfReduced <- data.frame(phi = numeric(), psi = numeric())


for (file in dssFilesReduced) {
  tryCatch({
    # Read in the DSSP file
    dssp <- read.dssp(file.path(path, file))
    
    # Extract the Ramachandran angles
    phi <- dssp$phi
    psi <- dssp$psi

    # Append the information to the dataframe
    dfReduced <- rbind(dfReduced, data.frame(phi = phi, psi = psi))
  }, error = function(e) {
    print(paste("Error in processing file:", file))
  })
}


#plotting secondary structure
par(mfrow=c(2,1))
#avg Ramachandran Plot
ramachandranDF = data.frame(df$phi, df$psi)
ramachandranMA = as.matrix(ramachandranDF)
ramachandran(ramachandranMA, xBins = 150,
               yBins = 150, printLegend = TRUE, barePlot = FALSE, main = "Averaged Ramachandran Plot")
#subsampled Ramachandran Plot
ramachandranDFSub = data.frame(dfReduced$phi, dfReduced$psi)
ramachandranMASub = as.matrix(ramachandranDFSub)
ramachandran(ramachandranMASub, xBins = 150,
             yBins = 150, printLegend = TRUE, barePlot = FALSE, main = "nORF Ramachandran Plot")

# comparing relative secondary structure
data$id = NULL
data$phi = NULL
data$psi = NULL
df_long = gather(data, variable, value)
ggbetweenstats(
  data  = df_long,
  x     = variable,
  y     = value,
  title = "nORF Secondary Structure Distribution",
)

# boxplot secondary structures
library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)

data("penguins", package = "palmerpenguins")
penguins <- drop_na(penguins)

plt <- ggbetweenstats(
  data = penguins,
  x = species,
  y = bill_length_mm
)

# nORFs Kernel
corrFrameReduced<- corrFrame[,c("type", "prop3.G1.residue100", "VS577", "DAYM780201.lag30", "BHAR880101.lag29", "CIDH920105.lag10", "BIGC670101.lag30")]


# load jointframe
dataset = na.omit(corrFrameReduced)
removeRandom = dataset$type == "Random" 
dataset = dataset[!removeRandom,]
dataset$type[dataset$type == "nORF"] = 1
dataset$type[dataset$type == "ORF"] = 0
dataset$id = NULL
dataset$seq = NULL
dataset$length = NULL
dataset$X = NULL
dataset$X.1 = NULL

#equalising the amount of each group
#dataset = rbind(head(dataset[dataset$group == 1,], 10000), head(dataset[dataset$group == 0,], 10000))

#dataset <- as.data.frame(scale(dataset))

## 75% of the sample size
smp_size <- floor(0.70 * nrow(dataset))

## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(dataset)), size = smp_size)

train <- dataset[train_ind, ]
test <- dataset[-train_ind, ]

## create 2 test set (hidden and tuning test set)

## 75% of the sample size
smp_test_size <- floor(0.50 * nrow(test))

## set the seed to make your partition reproducible
test_ind <- sample(seq_len(nrow(test)), size = smp_test_size)

test_tuning <- test[test_ind, ]
test_hidden <- test[-test_ind, ]


trainLabel = train$type
# trainLabel = as.numeric(trainLabel)
testTuningLabel = test_tuning$type
testHiddenLabel = test_hidden$type

train$type = NULL
test_tuning$type = NULL
test_hidden$type = NULL
row.names(train) <- NULL

trainTemp = train

train = as.matrix(train)
#train = as.numeric(train)
test = as.matrix(test_tuning)
test_hidden = as.matrix(test_hidden)


bstSparse <- xgboost(data = train, label = trainLabel, max.depth = 10, eta = 0.1, gamma = 7.5, alpha = 2, lambda = 0.4, nrounds = 250, verbose = FALSE,
                     objective = 'reg:linear', eval_metric = 'auc', prediction = T )

bstSparse2 <- xgboost(data = train, label = trainLabel, nrounds = 100,  eval_metric = 'auc', prediction = T )
bstSparse = bstSparse2

importance_matrix <- xgb.importance(model = bstSparse,  )
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix,  top_n = 50)


it = which.max(bstSparse$evaluation_log$test_auc_mean)
best.iter = bstSparse$evaluation_log$iter[it]

## test set
preds <- predict(bstSparse,test)

prediction <- as.numeric(preds > 0.5)

err <- mean(as.numeric(preds > 0.5) != testTuningLabel)
print(paste("test-error=", err))

preds <- predict(bstSparse,test)

##hidden dataset

preds <- predict(bstSparse,test_hidden)
prediction <- as.numeric(preds > 0.5)

err <- mean(as.numeric(preds > 0.5) != testHiddenLabel)
print(paste("test-error=", err))

#make latex table
rownames(importance_matrix) <- NULL

# Convert dataframe to LaTeX table
print(xtable(importance_matrix, caption = "Gini Importance Matrix", align = c("c","c","c","c","c","c")), 
      include.rownames = FALSE, 
      floating = FALSE, 
      hline.after = c(-1, 0))

confusionMatrix(as.factor(prediction), as.factor(testHiddenLabel))



#creating nORF kernel
corrFrameReduced2 = corrFrameReduced$prop3.G1.residue100 + corrFrameReduced$VS577 + corrFrameReduced$DAYM780201.lag30 + corrFrameReduced$BHAR880101.lag29 + corrFrameReduced$CIDH920105.lag10 + corrFrameReduced$BIGC670101.lag30
corrFrameReduced2 = corrFrameReduced$prop3.G1.residue100 + corrFrameReduced$VS577 + corrFrameReduced$DAYM780201.lag30 + corrFrameReduced$BHAR880101.lag29 + corrFrameReduced$CIDH920105.lag10 + corrFrameReduced$BIGC670101.lag30


corrFrameReduced3 = log(corrFrameReduced2)


#PCA kernel extraction
pcaFrame = corrFrameReduced
pcaFrame$type = NULL
pca_result <- prcomp(pcaFrame, scale = TRUE)

# Extract the first principal component
first_pc <- pca_result$x[,1]
second_pc <- pca_result$x[,2]

pcaDF = data.frame(PCA1 = first_pc, PCA2 = second_pc, type = arr$type)



library(ggplot2)

# Create an example array
arr <- data.frame(x = lengths, y = corrFrameReduced3, type = corrFrameReduced$type)

plot(arr$x, arr$y, col = colors)
# Plot the array using ggplot2
ggplot(arr, aes(x = x, y = y, color = type)) + 
  geom_point() + 
  ggtitle("Array Plot by Type")+
  xlab("AA Length")+
  ylab("nORF-Kernel")+
  scale_color_manual(values = c("nORF" = "red", "ORF" = "blue"))+
  theme(legend.title=element_blank()) +
  theme(legend.position = "topright")


ggplot(pcaDF, aes(x = PCA1, y = PCA2, color = type, group= type)) + 
  geom_point(alpha = 0.3) + 
  ggtitle("Reduced Kernel Dimension PCA Plot")+
  xlab("PCA1")+
  ylab("PCA2")+
  theme_minimal() +
  scale_color_manual(values = c("nORF" = "red", "ORF" = "blue"))



#Divide the screen in 1 line and 2 columns
par(
  mfrow=c(1,2), 
  oma = c(0, 0, 2, 0)
) 

par(new = TRUE)

#Make the margin around each graph a bit smaller
par(mar=c(4,2,2,2))

# Histogram and Scatterplot
hist(arr$x,  main="" , breaks=30 , col=colors , xlab="height" , ylab="dist", xlim = c(0,5000))
plot(arr$x , arr$y,  main="" , pch=20 , cex=0.4 , col=colors  , xlab="primadur" , ylab="Ixos", xlim = c(0,5000))

#And I add only ONE title :
mtext("Primadur : Distribution and correlation with Ixos", outer = TRUE, cex = 1.5, font=4, col=rgb(0.1,0.3,0.5,0.5) )



