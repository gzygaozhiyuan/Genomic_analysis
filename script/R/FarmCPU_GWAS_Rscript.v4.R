#!/PUBLIC/software/public/System/R-3.2.4/bin/Rscript
library("getopt")
library("bigmemory")
library("biganalytics")
library("compiler") #this library is already installed in R
#source("http://zzlab.net/GAPIT/gapit_functions.txt")
#source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
source("./gapit_functions.txt")
source("./FarmCPU_functions.txt")

#------------------------------------------------
#getting parameters
#------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help'      , 'h'      , 0, "logical",
	'genofile'  , 'geno'   , 1, "character",
	'phenofile' , 'pheno'  , 1, "character",
	'snpposinfo', 'info'   , 1, "character",
	'covariate' , 'cov'    , 2, "character",
	'outpath'   , 'out'    , 2, "character"
), byrow=TRUE, ncol=4);
spec
opt = getopt(spec);
opt
# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("Usage example: 
      1 Rscript xxx.R --genofile genotype.txt --phenofile phenotype.txt -- snpposinfo snpposinfo.txt --outpath /share/nas1/map 
      2 Rscript GOClassificationMap.R --infile in.dat --rotation 0.5 --outpath /share/nas1 --method log
      
      Options: 
      --help		-h	NULL	get this help
      --genofile	-geno   character	[forced]
      --phenofile	-pheno  character	[forced]
      --snpposinfo 	-info	character	[forced]
      --covariate       -cov    character       [optional]
      --outpath		-out	character	[optional]
      \n")
  q(status=1)
}

#if help was asked for print a friendly message
# and exit with a non-zero error code
if(!is.null(opt$help)){print_usage(spec)}
if ( is.null(opt$genofile) )  { print_usage(spec) }else {opt$genofile<-gsub("\\\\",replacement = "/",x = opt$genofile)}
genofile<-as.vector(opt$genofile)
if ( is.null(opt$phenofile) )  { print_usage(spec) }else {opt$phenofile<-gsub("\\\\",replacement = "/",x = opt$phenofile)}
phenofile<-as.vector(opt$phenofile)
if ( is.null(opt$snpposinfo) )  { print_usage(spec) }else {opt$snpposinfo<-gsub("\\\\",replacement = "/",x = opt$snpposinfo)}
snpposinfo<-as.vector(opt$snpposinfo)
if (is.null(opt$covariate))  {print_usage(spec)}else{opt$covariate<-gsub("\\\\",replacement="/",x=opt$covariate)}
covariate<-as.vector(opt$covariate)
if ( is.null(opt$outpath) )  { opt$outpath =getwd() }
outpath<-gsub("\\\\",replacement = "/",x = opt$outpath)
if(!file.exists(outpath)){dir.create(outpath,recursive = T)}

cat("#--------------------------------------------------------\n")
cat("#reading data\n")
cat("#--------------------------------------------------------\n")
if(file.exists(genofile)){
  genofiledata <- read.big.matrix((genofile),type="char", sep="\t", head = TRUE)
}else{
  stop("open genofile error:infile does not exist!")
}
genofiledata
if(file.exists(phenofile)){
  phenofiledata <- read.table((phenofile),sep="\t",head = TRUE)
}else{
  stop("open phenotype error:infile does not exist!")
}
dim(phenofiledata);head(phenofiledata)
if(file.exists(snpposinfo)){
  snpposinfodata <- read.table((snpposinfo),head = TRUE)
}else{
  stop("open snpposinfo error:infile does not exist!")
}
dim(snpposinfodata);head(snpposinfodata)
if(file.exists(covariate)){
	covariatedata <- read.table((covariate),head = TRUE,sep="\t")
	dim(covariatedata);head(covariatedata)
}else{
	print("open covariate error:covariate does not exist!")
}
dim(covariatedata);head(covariatedata)
setwd((outpath))
getwd()
#--------------------------------------------------------
# GWAS Analysis
#--------------------------------------------------------
TraitNo <- ncol(phenofiledata)-1
cat("TraitNo_Total:",TraitNo,"\n")
if(length(objects(pattern="covariatedata"))==1){
	for(i in 2:ncol(phenofiledata)){
		cat ("GWAS for ",i-1,"th:",colnames(phenofiledata)[i],"\n",sep="")
        	myFarmCPU <- FarmCPU(
                	Y=phenofiledata[,c(1,i)],
                	GD=genofiledata,
                	GM=snpposinfodata,
                	CV=covariatedata,
			cutOff=0.05
        	)
        	write.table(myFarmCPU$GWAS[,2:4],paste(outpath,"/","GWAS_FarmCPU_P_Results_",colnames(phenofiledata)[i],".txt",sep=""),row.names=F,quote=F,sep="\t")
	}
}else{
	for(i in 2:ncol(phenofiledata)){
		cat ("GWAS for ",i-1,"th:",colnames(phenofiledata)[i],"\n",sep="")
        	myFarmCPU <- FarmCPU(
                	Y=phenofiledata[,c(1,i)],
                	GD=genofiledata,
                	GM=snpposinfodata,
			cutOff=0.05
        	)
		write.table(myFarmCPU$GWAS[,2:4],paste(outpath,"/","GWAS_FarmCPU_Results_",colnames(phenofiledata)[i],".txt",sep=""),quote=F,row.names=F,sep="\t")
	}
}


#for(i in 2:TraitNo){
#	myFarmCPU <- FarmCPU(
#		Y=phenofiledata[,c(1,i)],
#		GD=genofiledata,
#		GM=snpposinfodata,
#		CV=covariatedata
#	)
#	write.table(myFarmCPU$GWAS,paste(outpath,"/","GWAS_FarmCPU_P_Results_",colnames(phenofiledata)[i],".txt",sep=""),quote=F)
#}


