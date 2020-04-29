#code was derived from the code developed in Atul Butte's lab
#author: Bin Chen (date: 2015)
#this code is to generate ews resistant signature
#samr and rankprod are used (or maybe bayesian)
#probe annotation can be downloaded from AILLUN or from GEO website
#expression data canbe taken from GEO matrix or derived from raw data
#cutoff 

#pre install the following package
source("http://www.bioconductor.org/biocLite.R")
#biocLite("impute")
#biocLite("siggenes")
#biocLite("RankProd")
#biocLite("preprocessCore")
set.seed(97)


#user input
setwd("/Users/User/Documents/stanford/repurpose/ews/data/array")
platform <- "GPL570"
gse <- "GSE12102"
matrix.file <- "GSE12102_series_matrix.txt" 

outcome.file <- "GSE12102_outcome.txt"
control.label <- "good"
treatment.label <- "poor"

method_id <- 2 #2: rankprod, 3: siggenes sam, 4: samr

q_thresh <- 0.05 #fdr 
sig_cutoff <- 0.05 # p value

dz_sig.file <- paste("GSE12102_dz_sig_",method_id,".txt",sep="")

#######################
#functions#####
# Very crude for now
is_logged <- function(sample_frame) {
  mean_sample_var <- mean(apply(sample_frame,2,var,na.rm=T),na.rm=T)
  ifelse(mean_sample_var > 10,FALSE,TRUE)
}

stripNumberPrefix <-function(v, prefix) {
  # Fix this to only work if there is a prefix!  Right now, will error!
  # Definitely need to fix this up!!
  
  #return(as.numeric(as.matrix(data.frame(strsplit(v, paste("^", prefix, sep="") ))[2,] )))
  #use GSUB instead?  
  return(as.matrix(gsub(paste("^", prefix, sep="" ), "", v)))
}

getEntrezMappings <- function(gpl, con, DEBUG=0) {
  query <- paste("select probe, GeneID, Symbol from annot_gpl.gpl_", stripNumberPrefix(gpl, "GPL"), sep="" )
  if (DEBUG) {print(query)}
  GeneMappings <- dbGetQuery(con, query)
  dupProbes <- unique(GeneMappings$probe[duplicated(GeneMappings$probe)])
  GeneMappings <- GeneMappings[! GeneMappings$probe %in% dupProbes,]
  rownames(GeneMappings) <- GeneMappings$probe
  return(GeneMappings)	
}

getEntrezMappingsFromGEO <- function(gpl){
  library(Biobase)
  library(GEOquery)
  gpl <- getGEO(gpl, destdir=".")
  geneMappings <- Table(gpl)[,c("ID","GENE","GENE_SYMBOL","NAME")]
  names(geneMappings) <- c("probe","GeneID","Symbol","NAME")
  geneMappings <- subset(geneMappings, !is.na(GeneID))
  dupProbes <- unique(geneMappings$probe[duplicated(geneMappings$probe)])
  geneMappings <- geneMappings[! geneMappings$probe %in% dupProbes,] #remove duplicated probes
  rownames(geneMappings) <- geneMappings$probe
  return(geneMappings)
}
######################



#create control and treatment group
outcome=read.table(outcome.file,header=T,stringsAsFactors=F,sep="\t")
control=outcome$sample[outcome$outcome== control.label]
treatment=outcome$sample[outcome$outcome== treatment.label]

#create matrix
fromGEO <- TRUE
#if (fromGEO){
#  data <- exprs(getGEO(gse,GSEMatrix=TRUE)[[1]])
#}else{
data<-read.table(matrix.file,header=TRUE,sep="\t",na.strings="NA")
#remove duplicated probes
data <- data[!duplicated(data$ID_REF),]
#}

#creat comparison
comparison_frame=subset(data,select=c(control,treatment))
row.names(comparison_frame)=data[,1]
sample_class=c(rep(0,length(control)),rep(1,length(treatment)))

#remove probes if half of values are missed
complete_probes <- apply(comparison_frame,1,function(row) {
  ctl_vals <- row[sample_class==0]
  dz_vals <- row[sample_class==1]
  missing_ctl <- length(ctl_vals[is.na(ctl_vals)]) / length(ctl_vals)
  missing_dz <- length(dz_vals[is.na(dz_vals)]) / length(dz_vals)
  ifelse(((missing_ctl > 0.5) || (missing_dz > 0.5)),FALSE,TRUE)
})
comparison_frame <- comparison_frame[complete_probes,]

# deal negative values
library(preprocessCore)
colSummarizeMedian(as.matrix(comparison_frame))
for (i in 1:ncol(comparison_frame)){
  comparison_frame[,i]=comparison_frame[,i]+abs(min(comparison_frame[,i],na.rm=TRUE))+1 #transform negatives to postivies, as quantile normailzation will be applied, it would affect 
}
                    
# deal log
values_logged <- is_logged(comparison_frame)
if (!values_logged) {
  comparison_frame <- log2(comparison_frame+0.001) # add 1 to avoid log2(0)
  values_logged <- TRUE
}

#normalization
boxplot(comparison_frame)
comparison_frame1 <- normalize.quantiles(as.matrix(comparison_frame ))
row.names(comparison_frame1) <- row.names(comparison_frame)
colnames(comparison_frame1) <- colnames(comparison_frame)
comparison_frame <- comparison_frame1


if (method_id==2) {
  method <- 'RANKPROD_SINGLE'
} else if (method_id==3){
  method <- "SIGGENES_SAMR"
} else if (method_id==4){
  method <- "SAMR"
}

#probe names
genenames<-rownames(comparison_frame)

if (method == 'RANKPROD_SINGLE') {
  library(RankProd)
  
  # Evaluate the impact of na.rm = TRUE. Alex M seems to think it's OK
  RP_result <- RP(comparison_frame, sample_class,gene.names=genenames, num.perm = 100, logged = T, na.rm = TRUE, plot = FALSE, rand = 123)
  # Leave logged=FALSE because topGene() converts the fold-change incorrectly!
  siggenes <- topGene(RP_result,cutoff=sig_cutoff,method="pfp",logged=FALSE,gene.names=genenames)
  # Normalize the results across methods
  siggenes.result = list()
  siggenes.result$UP <- siggenes$Table1[,3:5]
  colnames(siggenes.result$UP) <- c("fold.change","q.value","p.value")
  siggenes.result$DOWN <- siggenes$Table2[,3:5]
  colnames(siggenes.result$DOWN) <- c("fold.change","q.value","p.value")
  # Since RANKPROD does the goofy condition 1 / condition 2, inverse the fold-change and convert to log if not logged. 
  if (values_logged) {
    siggenes.result$UP[,"fold.change"] <- -siggenes.result$UP[,"fold.change"]
    siggenes.result$DOWN[,"fold.change"] <- -siggenes.result$DOWN[,"fold.change"]
  } else {
    siggenes.result$UP[,"fold.change"] <- log2(1/siggenes.result$UP[,"fold.change"])
    siggenes.result$DOWN[,"fold.change"] <- log2(1/siggenes.result$DOWN[,"fold.change"])
  }
  siggenes.result$UP=data.frame(siggenes.result$UP,probe=rownames(siggenes.result$UP))
  siggenes.result$DOWN=data.frame(siggenes.result$DOWN,probe=rownames(siggenes.result$DOWN))
  
    
} else if (method == 'SIGGENES_SAMR') {
  # Using only the SAM implementation in SIGGENES, not EBAM
  library(siggenes)

  SAM_result <- sam(comparison_frame,sample_class,rand=123,R.unlog=T,gene.names=genenames) #q.version=1
  delta.table <- rbind(c(0,0,0),findDelta(SAM_result,fdr=0.9)) #fdr=0.9, too loose
  siggenes.table <- summary(SAM_result, delta.table[  dim(delta.table)[1],  1] );
  siggenes <- siggenes.table@mat.sig
  
  if( nrow(siggenes) > 0 ) {
    siggenes.result = list()
    siggenes.result$UP <- siggenes[siggenes$d.value > 0,c("R.fold","q.value","rawp")]
    colnames(siggenes.result$UP) <- c("fold.change","q.value","p.value")
    siggenes.result$DOWN <- siggenes[siggenes$d.value < 0,c("R.fold","q.value","rawp")]
    colnames(siggenes.result$DOWN) <- c("fold.change","q.value","p.value")
 
    
    siggenes.result$UP=data.frame(siggenes.result$UP,probe=rownames(siggenes.result$UP))
    siggenes.result$DOWN=data.frame(siggenes.result$DOWN,probe=rownames(siggenes.result$DOWN))
    
  }
  
} else if (method == 'SAMR') {
  library(samr)
  
  # class labels for SAMR are 1 or 2
  input_data <- list(
    x=data.matrix(comparison_frame),
    y=sample_class+1,
    geneid=rownames(comparison_frame),
    genenames=rownames(comparison_frame),
    logged2=T
  )
  
  samr.obj <- samr(input_data,resp.type="Two class unpaired",testStatistic="standard",nperms=100)
  delta.table <- samr.compute.delta.table(samr.obj)
  delta.index <- which.max( delta.table[,5] < 0.4 )
  delta=delta.table[delta.index,1]
  
  #replace data with input_data
  siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, input_data, delta.table)

  siggenes.result = list()
  
  if (!is.null(siggenes.table$genes.up)>0){
    UP <- as.data.frame(subset(siggenes.table$genes.up,select=c("Gene ID","Fold Change","q-value(%)")))
    colnames(UP) <- c("probe","fold.change","q.value")
    UP$p.value  <- NA    
    UP$probe <- as.character(UP$probe)
    UP$fold.change <- as.double(as.character(UP$fold.change))
    UP$q.value <- as.double(as.character(UP$q.value))
    siggenes.result$UP <- UP
  }
  if (!is.null(siggenes.table$genes.lo)){
    DOWN <- as.data.frame(subset(siggenes.table$genes.lo,select=c("Gene ID","Fold Change","q-value(%)")))
    colnames(DOWN) <- c("probe","fold.change","q.value")
    DOWN$p.value  <- NA    
    DOWN$probe <- as.character(DOWN$probe)
    DOWN$fold.change <- as.double(as.character(DOWN$fold.change))
    DOWN$q.value <- as.double(as.character(DOWN$q.value))
    siggenes.result$DOWN <- DOWN
  }
  
  
} 

######the annotations from GEO are twice bigger than those from AILLUN (which might be out of date)
egm <- getEntrezMappings(platform,con)
#egm2 <- getEntrezMappingsFromGEO(platform)
#egm <- if (nrow(egm1)>nrow(egm2)) egm1 else egm2 
##################################

#annotate probe by merging with GPL
genes.up <- merge(siggenes.result$UP,egm, by.x="probe",by.y="probe",sort=F,all.x=T)
genes.down <- merge(siggenes.result$DOWN,egm, by.x="probe",by.y="probe",sort=F,all.x=T)

#missing match
print(paste("missing matched up probes #", sum(is.na(genes.up$GeneID)) + sum(genes.up$GeneID==""),sep=""))
print(paste("missing matched down probes #", sum(is.na(genes.down$GeneID)) + sum(genes.down$GeneID==""),sep=""))

#removed missed match probes
genes.up <- genes.up[!is.na(genes.up$GeneID) & genes.up$GeneID!="",]
genes.down <- genes.down[!is.na(genes.down$GeneID) & genes.down$GeneID!="",]

genes.up[,"up_down"] <- "up"
genes.down[,"up_down"] <- "down"

genes.sig.up <- subset(genes.up,q.value<q_thresh, select=c("probe","GeneID","Symbol","fold.change","q.value","p.value","up_down"))
genes.sig.down <- subset(genes.down, q.value<q_thresh, select=c("probe","GeneID","Symbol","fold.change","q.value","p.value","up_down"))

#write disease signature
genes.sig.up <- genes.sig.up[order(-genes.sig.up$fold.change),]
genes.sig.down <- genes.sig.down[order(genes.sig.down$fold.change),]

#genes.sig.up[regexpr("KIT",genes.sig.up$Symbol)>0,]
#genes.sig.up[regexpr("ETV1",genes.sig.up$Symbol)>0,]

genes.sig.up <- genes.sig.up[1:min(nrow(genes.sig.up),200),]
genes.sig.down <- genes.sig.down[1:min(nrow(genes.sig.down),200),]
dz_sig <- rbind(genes.sig.up,genes.sig.down)

write.table(dz_sig,dz_sig.file,sep="\t",quote=F,row.names=F,col.names=T)

