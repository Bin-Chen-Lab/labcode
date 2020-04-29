#Author: Bin Chen (2015)
#Input a list of genes and a given cancer, compute their association with survival.
#user can change gene symbol, cancer name or molecular type.

library(cgdsr)
library(survival)
library(survminer)

mycgds = CGDS("https://www.cbioportal.org/")

test(mycgds)

# Get list of cancer studies at server
a = getCancerStudies(mycgds)
symbols = c("IFI27", "JUP", "LAX1", "HK3", "TNIP1", "GPAA1", "CTSB" ) #"ACTB", "CLTA", "ITGB1",

# Get available case lists (collection of samples) for a given cancer study
mycancerstudies = getCancerStudies(mycgds)[,1]

# Get data slices for a specified list of genes, genetic profile and case list
symbol_set = symbols #split(symbols, ceiling(seq_along(symbols)/90))

mycancerstudy = "brca_tcga_pub2015" #brca_metabric brca_tcga_pub2015
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
myclinicaldata = getClinicalData(mycgds,mycaselist)
#myclinicaldata = myclinicaldata[myclinicaldata$HER2_STATUS == "+",]
mygeneticprofiles = getGeneticProfiles(mycgds,mycancerstudy)[,1]
mygeneticprofile = mygeneticprofiles[grep("mrna", mygeneticprofiles)][1]
mrna = getProfileData(mycgds,c(as.character(symbols)) ,mygeneticprofile,mycaselist)

expr_avg = apply(mrna, 1, function(x) mean(x, na.rm = T))
expr_clinic = merge(data.frame(barcode = rownames(myclinicaldata), myclinicaldata), 
                    data.frame(barcode = rownames(mrna), infect_expr = expr_avg  ), by = "barcode" )
expr_clinic = subset(expr_clinic, !is.na(infect_expr))

expr_clinic$infect = ifelse(expr_clinic$infect_expr > quantile(expr_clinic$infect_expr)[3], "yes", "no")

expr_clinic = expr_clinic[!is.na(expr_clinic$OS_STATUS),]
expr_clinic$OS_STATUS_BIN = 1
expr_clinic$OS_STATUS_BIN[expr_clinic$OS_STATUS == "LIVING"] = 0

survdiff(Surv(OS_MONTHS, OS_STATUS_BIN) ~ infect, data=expr_clinic)

my.fit1 <- survfit(Surv(OS_MONTHS, OS_STATUS_BIN) ~ infect, data =  expr_clinic) 
plot(my.fit1)# uvm_tcga very strong
fit <- survival::coxph(Surv(OS_MONTHS, OS_STATUS_BIN) ~ infect, data=expr_clinic)
print(paste(mycancerstudy, anova(fit, test="Chisq")))

ggsurvplot(my.fit1, xlab = "months", legend.title = "", pval = T)

