##################################################################################################################
### 1, read in and prepare data (750 x 600)
##################################################################################################################
#update.packages(ask = FALSE)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("ggsurvfit")
BiocManager::install("Biobase")
BiocManager::install("GEOquery")
install.packages("writexl")
library(Biobase)
library(GEOquery)
library(writexl)
library(ggsurvfit)
getwd()

# Expression Data (Assay Data)
x <- read.table("assayData_Biopsy_Cibersort_n704_v3.txt", sep="\t", header = T, row.names = 1)
class(x) # Data frame
x <- as.matrix(x) 
class(x) # Matrix array

# Feature Data
f <- read.table("featureData_Biopsy_Cibersort_n704_v3.txt", sep="\t", header =T, row.names = 1)
dim(f)
class(f)

# Phenotype data
p <- read.table("Phenotype_n704_v3.txt", sep="\t", header =T, row.names = 1)
dim(p)
class(p)

library(dplyr)


####################################################################################################################
# Create Expression Set
####################################################################################################################
#BiocManager::install("Biobase")
library(Biobase)

eset <- ExpressionSet(assayData = x,
                      phenoData = AnnotatedDataFrame(p),
                      featureData = AnnotatedDataFrame(f))

#View the number of features/Proteins (rows) and samples (columns)
dim(eset)
# Access data from an ExpressionSet object
x <- exprs(eset) # you can retreive the expression matrix with the function `exprs`
f <- fData(eset) # Retrieve the feature data with `fData'
p <- pData(eset)# Retrieve the phenotype data with`pData`

hist(x)

# transform the expression data to Z scores                        
#x <- t(scale(t(x)))                                               
#hist(x)
# extract information of interest from the phenotype data (pdata)
idx <- which(colnames(pData(eset)) %in%
               c('CensOS','CensEFS', 'YearsOS', 'YearsEFS'))

metadata <- data.frame(pData(eset)[,idx],
                       row.names = rownames(pData(eset)))

# check that sample names match exactly between pdata and Z-scores 
all((colnames(x) == rownames(metadata)) == TRUE)

# create a merged pdata and Z-scores object
coxdata <- data.frame(metadata, t(x))

library(tibble)
coxdata_df <- tibble::rownames_to_column(coxdata, "CHIPID")

write_xlsx(coxdata_df,"YC_coxdata_df.xlsx")
# prepare phenotypes
coxdata$CensOS <- as.numeric(coxdata$CensOS) # Censor for the data
coxdata$YearsOS<- as.numeric(gsub('^KJX|^KJ', '', coxdata$YearsOS)) # Similar to DaysOS
coxdata$CensEFS <- as.numeric(coxdata$CensEFS) # Censor for the data
coxdata$YearsEFS <- as.numeric(gsub('^KJX|^KJ', '', coxdata$YearsEFS)) # Similar to DaysOS

#########################################################################################################
## Cutpoint (OS)
#########################################################################################################
library("survival")
library("survminer")
# surv_cutpoint(): Determine the optimal cutpoint for each variable using 'maxstat'.
res.cut_OS <-surv_cutpoint(
  coxdata,
  time = "YearsOS",
  event = "CensOS",
  variables = colnames(coxdata)[5:ncol(coxdata)],
  minprop = 0.05,
  progressbar = TRUE
)
options(max.print = 5000)        # Change global options
summary(res.cut_OS)

library(data.table)
res.cut_OS_df <- setDT(summary(res.cut_OS), keep.rownames = "Uniprot")[]
write_xlsx(res.cut_OS_df,"YC_RES_CUT_OS.xlsx")

# RES.CAT (categorize by High and Low expression from cutpoints)
res.cat_OS <- surv_categorize(res.cut_OS) ### FINE
head(res.cat_OS)
summary(res.cat_OS)
str(res.cat_OS)
write_xlsx(res.cat_OS,"YC_RESCAT.xlsx")

# 4. Fit survival curves and visualize (750 x 600)
{fit_OS <- survfit(Surv(YearsOS, CensOS) ~Neutrophils, data = res.cat_OS)
  a <- surv_pvalue(fit_OS)$pval
  b <- unname(summary(fit_OS)$table[,'records'])
  print(a)
  print(b)
  print(fit_OS)
  ggsurvplot(fit_OS, data = res.cat_OS, risk.table = TRUE, conf.int = F, pval = TRUE)
}

ggsurv_OS <- ggsurvplot(fit_OS, data = res.cat_OS, 
                        break.time.by = 5,
                        conf.int=F, pval=T, pval.size = 9, risk.table = "nrisk_cumevents", 
                        surv.median.line = "hv", # Specify median survival
                        tables.theme = theme_survminer(
                          font.main = c(20, "bold", "black"),
                          font.submain = c(18, "bold.italic", "purple"),
                          font.caption = c(18, "plain", "oCBX3ge"),
                          font.x = c(20, "bold.italic", "black"),
                          font.y = c(20, "bold.italic", "black"),
                          font.tickslab = c(18, "plain", "black")),
                        legend.labs=c("High","Low"), legend.title="% Population",  
                        palette=c("red3","blue"), xlab = "Time (Years)", 
                        title="             Neutrophils (OS)",
                        caption = "created with survminer",
                        font.title = c(35, "bold", "black"),
                        font.subtitle = c(25, "bold.italic", "purple"),
                        font.legend = c(25),
                        font.x = c(20, "bold.italic", "black"),
                        xlim = c(-1, 26), # to the graph
                        font.y = c(20, "bold.italic", "black"),
                        font.tickslab = c(20, "plain", "black"),
                        risk.table.fontsize = 5.5,
                        ########## risk table #########,
                        risk.table.height = 0.25)

ggsurv_OS


#####################################################################################################################
# Kaplan Meier Curves! NOT optimal cutpoint
#####################################################################################################################
fit_OS_HR <- survfit(Surv(YearsOS, CensOS) ~Neutrophils, data = coxdata)
str(fit_OS_HR)  
plot(fit_OS_HR)

survfit2(Surv(YearsOS, CensOS) ~ Neutrophils, data = coxdata) %>% 
  ggsurvfit() +
  labs(
    x = "Years (OS)",
    y = "Overall survival probability"
  ) + 
add_confidence_interval() +
  add_risktable()

summary(survfit(Surv(YearsOS, CensOS) ~ 1, data = coxdata), times = 1)

##test Neutrophils OS
coxdata$Neutrophils <- as.numeric(coxdata$Neutrophils)
str(coxdata)
coxdata$Neutrophils_group <- ifelse(coxdata$Neutrophils <= 0.123342974, "low", "high")
table(coxdata$Neutrophils_group)
fit <- survfit(Surv(YearsOS, CensOS) ~ Neutrophils_group, data = coxdata)
summary(fit)
ggsurvplot(fit, data = coxdata, pval = TRUE, risk.table = TRUE, conf.int = F)

ggsurvplot(fit, data = coxdata, 
                        break.time.by = 5,
                        conf.int=F, pval=T, pval.size = 9, risk.table = "nrisk_cumevents", 
                        surv.median.line = "hv", # Specify median survival
                        tables.theme = theme_survminer(
                          font.main = c(20, "bold", "black"),
                          font.submain = c(18, "bold.italic", "purple"),
                          font.caption = c(18, "plain", "oCBX3ge"),
                          font.x = c(20, "bold.italic", "black"),
                          font.y = c(20, "bold.italic", "black"),
                          font.tickslab = c(18, "plain", "black")),
                        legend.labs=c("High","Low"), legend.title="% Population",  
                        palette=c("red3","blue"), xlab = "Time (Years)", 
                        title="         Neutrophils (OS)",
                        caption = "created with survminer",
                        font.title = c(35, "bold", "black"),
                        font.subtitle = c(25, "bold.italic", "purple"),
                        font.legend = c(25),
                        font.x = c(20, "bold.italic", "black"),
                        xlim = c(-1, 26), # to the graph
                        font.y = c(20, "bold.italic", "black"),
                        font.tickslab = c(20, "plain", "black"),
                        risk.table.fontsize = 5.5,
                        ########## risk table #########,
                        risk.table.height = 0.25)

##test Monocytes
coxdata$Monocytes <- as.numeric(coxdata$Monocytes)
str(coxdata)
coxdata$Monocytes_group <- ifelse(coxdata$Monocytes <= 0.059031854, "low", "high")
table(coxdata$Monocytes_group)
fit <- survfit(Surv(YearsOS, CensOS) ~ Monocytes_group, data = coxdata)
summary(fit)
#ggsurvplot(fit, data = coxdata, pval = TRUE, risk.table = TRUE, conf.int = F)

ggsurvplot(fit, data = coxdata, 
           break.time.by = 5,
           conf.int=F, pval=T, pval.size = 9, risk.table = "nrisk_cumevents", 
           surv.median.line = "hv", # Specify median survival
           tables.theme = theme_survminer(
             font.main = c(20, "bold", "black"),
             font.submain = c(18, "bold.italic", "purple"),
             font.caption = c(18, "plain", "oCBX3ge"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black")),
           legend.labs=c("High","Low"), legend.title="% Population",  
           palette=c("red3","blue"), xlab = "Time (Years)", 
           title="           Monocytes (OS)",
           caption = "created with survminer",
           font.title = c(35, "bold", "black"),
           font.subtitle = c(25, "bold.italic", "purple"),
           font.legend = c(25),
           font.x = c(20, "bold.italic", "black"),
           xlim = c(-1, 26), # to the graph
           font.y = c(20, "bold.italic", "black"),
           font.tickslab = c(20, "plain", "black"),
           risk.table.fontsize = 5.5,
           ########## risk table #########,
           risk.table.height = 0.25)

##test Mast cells resting
coxdata$Mast_cells_resting <- as.numeric(coxdata$Mast_cells_resting)
str(coxdata)
coxdata$Mast_cells_resting_group <- ifelse(coxdata$Mast_cells_resting <= 0.044010462, "low", "high")
table(coxdata$Mast_cells_resting_group)
fit <- survfit(Surv(YearsOS, CensOS) ~ Mast_cells_resting_group, data = coxdata)
summary(fit)

#ggsurvplot(fit, data = coxdata, pval = TRUE, risk.table = TRUE, conf.int = F)

ggsurvplot(fit, data = coxdata, 
           break.time.by = 5,
           conf.int=F, pval=T, pval.size = 9, risk.table = "nrisk_cumevents", 
           surv.median.line = "hv", # Specify median survival
           tables.theme = theme_survminer(
             font.main = c(20, "bold", "black"),
             font.submain = c(18, "bold.italic", "purple"),
             font.caption = c(18, "plain", "oCBX3ge"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black")),
           legend.labs=c("High","Low"), legend.title="% Population",  
           palette=c("red3","blue"), xlab = "Time (Years)", 
           title="           Mast cells (OS)",
           caption = "created with survminer",
           font.title = c(35, "bold", "black"),
           font.subtitle = c(25, "bold.italic", "purple"),
           font.legend = c(25),
           font.x = c(20, "bold.italic", "black"),
           xlim = c(-1, 26), # to the graph
           font.y = c(20, "bold.italic", "black"),
           font.tickslab = c(20, "plain", "black"),
           risk.table.fontsize = 5.5,
           ########## risk table #########,
           risk.table.height = 0.25)

##test Macrophages_M2
coxdata$Macrophages_M2 <- as.numeric(coxdata$Macrophages_M2)

coxdata$Macrophages_M2_group <- ifelse(coxdata$Macrophages_M2 <= 0.073344083, "low", "high")

fit <- survfit(Surv(YearsOS, CensOS) ~ Macrophages_M2_group, data = coxdata)
summary(fit)

ggsurvplot(fit, data = coxdata, 
           break.time.by = 5,
           conf.int=F, pval=T, pval.size = 9, risk.table = "nrisk_cumevents", 
           surv.median.line = "hv", # Specify median survival
           tables.theme = theme_survminer(
             font.main = c(20, "bold", "black"),
             font.submain = c(18, "bold.italic", "purple"),
             font.caption = c(18, "plain", "oCBX3ge"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black")),
           legend.labs=c("High","Low"), legend.title="% Population",  
           palette=c("red3","blue"), xlab = "Time (Years)", 
           title="      Macrophages_M2 (OS)",
           caption = "created with survminer",
           font.title = c(35, "bold", "black"),
           font.subtitle = c(25, "bold.italic", "purple"),
           font.legend = c(25),
           font.x = c(20, "bold.italic", "black"),
           xlim = c(-1, 26), # to the graph
           font.y = c(20, "bold.italic", "black"),
           font.tickslab = c(20, "plain", "black"),
           risk.table.fontsize = 5.5,
           ########## risk table #########,
           risk.table.height = 0.25)

##not optimal cutpoint_T_cells_CD8
coxdata$T_cells_CD8 <- as.numeric(coxdata$T_cells_CD8)

coxdata$T_cells_CD8_group <- ifelse(coxdata$T_cells_CD8 <= 0.087008523, "low", "high")

fit <- survfit(Surv(YearsOS, CensOS) ~ T_cells_CD8_group, data = coxdata)
summary(fit)

ggsurvplot(fit, data = coxdata, 
           break.time.by = 5,
           conf.int=F, pval=T, pval.size = 9, risk.table = "nrisk_cumevents", 
           surv.median.line = "hv", # Specify median survival
           tables.theme = theme_survminer(
             font.main = c(20, "bold", "black"),
             font.submain = c(18, "bold.italic", "purple"),
             font.caption = c(18, "plain", "oCBX3ge"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black")),
           legend.labs=c("High","Low"), legend.title="% Population",  
           palette=c("red3","blue"), xlab = "Time (Years)", 
           title="         T_cells_CD8 (OS)",
           caption = "created with survminer",
           font.title = c(35, "bold", "black"),
           font.subtitle = c(25, "bold.italic", "purple"),
           font.legend = c(25),
           font.x = c(20, "bold.italic", "black"),
           xlim = c(-1, 26), # to the graph
           font.y = c(20, "bold.italic", "black"),
           font.tickslab = c(20, "plain", "black"),
           risk.table.fontsize = 5.5,
           ########## risk table #########,
           risk.table.height = 0.25)

##not optimal cutpoint_T_cells_CD4_naive
coxdata$T_cells_CD4_naive <- as.numeric(coxdata$T_cells_CD4_naive)

coxdata$T_cells_CD4_naive_group <- ifelse(coxdata$T_cells_CD4_naive <= 0.0000924591, "low", "high")

fit <- survfit(Surv(YearsOS, CensOS) ~ T_cells_CD4_naive_group, data = coxdata)
summary(fit)

ggsurvplot(fit, data = coxdata, 
           break.time.by = 5,
           conf.int=F, pval=T, pval.size = 9, risk.table = "nrisk_cumevents", 
           surv.median.line = "hv", # Specify median survival
           tables.theme = theme_survminer(
             font.main = c(20, "bold", "black"),
             font.submain = c(18, "bold.italic", "purple"),
             font.caption = c(18, "plain", "oCBX3ge"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black")),
           legend.labs=c("High","Low"), legend.title="% Population",  
           palette=c("red3","blue"), xlab = "Time (Years)", 
           title="    T_cells_CD4_naive (OS)",
           caption = "created with survminer",
           font.title = c(35, "bold", "black"),
           font.subtitle = c(25, "bold.italic", "purple"),
           font.legend = c(25),
           font.x = c(20, "bold.italic", "black"),
           xlim = c(-1, 26), # to the graph
           font.y = c(20, "bold.italic", "black"),
           font.tickslab = c(20, "plain", "black"),
           risk.table.fontsize = 5.5,
           ########## risk table #########,
           risk.table.height = 0.25)

##not optimal cutpoint_T_cells_gamma_delta
coxdata$T_cells_gamma_delta <- as.numeric(coxdata$T_cells_gamma_delta)

coxdata$T_cells_gamma_delta_group <- ifelse(coxdata$T_cells_gamma_delta <= 0.045859531, "low", "high")

fit <- survfit(Surv(YearsOS, CensOS) ~ T_cells_gamma_delta_group, data = coxdata)
summary(fit)

ggsurvplot(fit, data = coxdata, 
           break.time.by = 5,
           conf.int=F, pval=T, pval.size = 9, risk.table = "nrisk_cumevents", 
           surv.median.line = "hv", # Specify median survival
           tables.theme = theme_survminer(
             font.main = c(20, "bold", "black"),
             font.submain = c(18, "bold.italic", "purple"),
             font.caption = c(18, "plain", "oCBX3ge"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black")),
           legend.labs=c("High","Low"), legend.title="% Population",  
           palette=c("red3","blue"), xlab = "Time (Years)", 
           title=" T_cells_gamma_delta (OS)",
           caption = "created with survminer",
           font.title = c(35, "bold", "black"),
           font.subtitle = c(25, "bold.italic", "purple"),
           font.legend = c(25),
           font.x = c(20, "bold.italic", "black"),
           xlim = c(-1, 26), # to the graph
           font.y = c(20, "bold.italic", "black"),
           font.tickslab = c(20, "plain", "black"),
           risk.table.fontsize = 5.5,
           ########## risk table #########,
           risk.table.height = 0.25)

##not optimal cutpoint_T_cells_CD4_memory_activated
coxdata$T_cells_CD4_memory_activated <- as.numeric(coxdata$T_cells_CD4_memory_activated)

coxdata$T_cells_CD4_memory_activated_group <- ifelse(coxdata$T_cells_CD4_memory_activated <= 0.072091666, "low", "high")

fit <- survfit(Surv(YearsOS, CensOS) ~ T_cells_CD4_memory_activated_group, data = coxdata)
summary(fit)

ggsurvplot(fit, data = coxdata, 
           break.time.by = 5,
           conf.int=F, pval=T, pval.size = 9, risk.table = "nrisk_cumevents", 
           surv.median.line = "hv", # Specify median survival
           tables.theme = theme_survminer(
             font.main = c(20, "bold", "black"),
             font.submain = c(18, "bold.italic", "purple"),
             font.caption = c(18, "plain", "oCBX3ge"),
             font.x = c(20, "bold.italic", "black"),
             font.y = c(20, "bold.italic", "black"),
             font.tickslab = c(18, "plain", "black")),
           legend.labs=c("High","Low"), legend.title="% Population",  
           palette=c("red3","blue"), xlab = "Time (Years)", 
           title="T_cells_CD4_memory_activated (OS)",
           caption = "created with survminer",
           font.title = c(25, "bold", "black"),
           font.subtitle = c(25, "bold.italic", "purple"),
           font.legend = c(25),
           font.x = c(20, "bold.italic", "black"),
           xlim = c(-1, 26), # to the graph
           font.y = c(20, "bold.italic", "black"),
           font.tickslab = c(20, "plain", "black"),
           risk.table.fontsize = 5.5,
           ########## risk table #########,
           risk.table.height = 0.25)


#################################################################################################
## Cutpoint (EFS)
#################################################################################################
library("survival")
library(survminer)
# surv_cutpoint(): Determine the optimal cutpoint for each variable using 'maxstat'.
res.cut_EFS <-surv_cutpoint(
  coxdata,
  time = "YearsEFS",
  event = "CensEFS",
  variables = colnames(coxdata)[5:ncol(coxdata)],
  minprop = 0.05,
  progressbar = TRUE
)
options(max.print = 5000)        # Change global options
summary(res.cut_EFS)

res.cut_EFS_df <- setDT(summary(res.cut_EFS), keep.rownames = "Uniprot")[]
write_xlsx(res.cut_EFS_df,"YC_RES_CUT_EFS.xlsx")

# RES.CAT (categorize by High and Low expression from cutpoints)
res.cat_EFS <- surv_categorize(res.cut_EFS) ### FINE
head(res.cat_EFS)
summary(res.cat_EFS)
str(res.cat_EFS)
write_xlsx(res.cat_EFS,"YC_RESCAT_EFS.xlsx")

# 4. Fit survival curves and visualize (800 x 762)

 {fit_EFS <- survfit(Surv(YearsEFS, CensEFS) ~T_cells_gamma_delta, data = res.cat_EFS)
   a <- surv_pvalue(fit_EFS)$pval
   b <- unname(summary(fit_EFS)$table[,'records'])
   print(a)
   print(b)
   ggsurvplot(fit_EFS, data = res.cat_EFS, risk.table = TRUE, conf.int = F, pval = TRUE)
   }


ggsurv_EFS <- ggsurvplot(fit_EFS, data = res.cat_EFS, 
                        break.time.by = 5,
                        conf.int=F, pval=T, pval.size = 9, risk.table = "nrisk_cumevents", 
                        surv.median.line = "hv", # Specify median survival
                        tables.theme = theme_survminer(
                          font.main = c(20, "bold", "black"),
                          font.submain = c(18, "bold.italic", "purple"),
                          font.caption = c(18, "plain", "oCBX3ge"),
                          font.x = c(20, "bold.italic", "black"),
                          font.y = c(20, "bold.italic", "black"),
                          font.tickslab = c(18, "plain", "black")),
                        legend.labs=c("High","Low"), legend.title="% Population",  
                        palette=c("red3","blue"), xlab = "Time (Years)", 
                        title=" T_cells_gamma_delta (EFS)",
                        caption = "created with survminer",
                        font.title = c(35, "bold", "black"),
                        font.subtitle = c(25, "bold.italic", "purple"),
                        font.legend = c(25),
                        font.x = c(20, "bold.italic", "black"),
                        xlim = c(-1, 26), # to the graph
                        font.y = c(20, "bold.italic", "black"),
                        font.tickslab = c(20, "plain", "black"),
                        risk.table.fontsize = 5.5,
                        ########## risk table #########,
                        risk.table.height = 0.25)


ggsurv_EFS
getwd()


################################################################################################
### OS Hazard Ratios (OS cutpoint)
################################################################################################

library(data.table)
library(writexl)

Cell_Types <- read.csv("Cell Types.csv") #need to read in .csv file, not .txt file
Cell_Types <- t(Cell_Types) #Need to transform the data to get a vector with .csv

## To apply the univariate coxph function to multiple covariates at once for OS, type this:
univ_formulas_OS <- sapply(Cell_Types,
                           function(x) as.formula(paste('Surv(YearsOS, CensOS)~',x)))

univ_models_OS <- lapply(univ_formulas_OS, function(x){coxph(x, data = res.cat_OS)})

#Extract data
univ_results_OS <- lapply(univ_models_OS,
                          function(x){
                            x <- summary(x)
                            p.value <- signif(x$wald["pvalue"], digits = 3)
                            wald.test <- signif(x$wald["test"],digits = 2)
                            beta <- signif(x$coef[1],digits = 3); #coefficient beta
                            HR <- signif(x$coef[2],digits = 3); #exp(beta)
                            HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                            HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                            HR <- paste0(HR, " (",
                                         HR.confint.lower, "-",HR.confint.upper, ")")
                            res <- c(beta, HR, wald.test, p.value)
                            names(res) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                            return(res)
                            #return(exp(cbind(coef(x), confint(x))))
                          })
res_OS <- t(as.data.frame(univ_results_OS, check.names = FALSE))
res_OS_DF <- as.data.frame(res_OS)

res_OS_DF <- setDT(res_OS_DF, keep.rownames = "Cell Types")[]
getwd()
write_xlsx(res_OS_DF,"Hazard Ratios OS.xlsx")

################################################################################################
### EFS Hazard Ratios (EFS Cutpoint)
################################################################################################
library(data.table)
library(writexl)

Cell_Types <- read.csv("Cell Types.csv") #need to read in .csv file, not .txt file
Cell_Types <- t(Cell_Types) #Need to transform the data to get a vector with .csv

## To apply the univariate coxph function to multiple covariates at once for OS, type this:
univ_formulas_EFS <- sapply(Cell_Types,
                           function(x) as.formula(paste('Surv(YearsEFS, CensEFS)~',x)))

univ_models_EFS <- lapply(univ_formulas_EFS, function(x){coxph(x, data = res.cat_EFS)})
#Extract data
univ_results_EFS <- lapply(univ_models_EFS,
                          function(x){
                            x <- summary(x)
                            p.value <- signif(x$wald["pvalue"], digits = 3)
                            wald.test <- signif(x$wald["test"],digits = 2)
                            beta <- signif(x$coef[1],digits = 3); #coefficient beta
                            HR <- signif(x$coef[2],digits = 3); #exp(beta)
                            HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                            HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                            HR <- paste0(HR, " (",
                                         HR.confint.lower, "-",HR.confint.upper, ")")
                            res <- c(beta, HR, wald.test, p.value)
                            names(res) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                            return(res)
                            #return(exp(cbind(coef(x), confint(x))))
                          })
res_EFS <- t(as.data.frame(univ_results_EFS, check.names = FALSE))
res_EFS_DF <- as.data.frame(res_EFS)

res_EFS_DF <- setDT(res_EFS_DF, keep.rownames = "Cell Types")[]
getwd()
write_xlsx(res_EFS_DF,"Hazard Ratios EFS.xlsx")

################################################################################################
### OS Hazard Ratios (EFS Cutpoint)
################################################################################################

#Download EFS Cutpoint Table
class(res.cat_EFS)
res.cat_EFS_df <- data.frame(res.cat_EFS)
class(res.cat_EFS_df)
res.cat_EFS_df <- setDT(res.cat_EFS_df, keep.rownames = "CHIPID")[]
write_xlsx(res.cat_EFS_df,"EFS_Cutpoints.xlsx")

# Change the EFS data (CensEFS, YearsEFS) with OS data (CensOS, YearsOS) in Excel 
# Read back in .txt file of EFS cutpoints with OS data
EFS_Cutpoints_OS_Data <- read.table("EFS_Cutpoints_OS_Data.txt", sep="\t", header = T, row.names = 1)

## To apply the univariate coxph function to multiple covariates at once for OS, type this:
univ_formulas <- sapply(Cell_Types,
                            function(x) as.formula(paste('Surv(YearsOS, CensOS)~',x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = EFS_Cutpoints_OS_Data)})
#Extract data
univ_results <- lapply(univ_models,
                           function(x){
                             x <- summary(x)
                             p.value <- signif(x$wald["pvalue"], digits = 3)
                             wald.test <- signif(x$wald["test"],digits = 2)
                             beta <- signif(x$coef[1],digits = 3); #coefficient beta
                             HR <- signif(x$coef[2],digits = 3); #exp(beta)
                             HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                             HR.confint.upper<- signif(x$conf.int[,"upper .95"],3)
                             HR <- paste0(HR, " (",
                                          HR.confint.lower, "-",HR.confint.upper, ")")
                             res <- c(beta, HR, wald.test, p.value)
                             names(res) <- c("beta","HR (95% CI for HR)", "wald.test", "p.value")
                             return(res)
                             #return(exp(cbind(coef(x), confint(x))))
                           })
EFS_Cutpoints_OS_Data <- t(as.data.frame(univ_results, check.names = FALSE))
EFS_Cutpoints_OS_Data_DF <- as.data.frame(EFS_Cutpoints_OS_Data)

EFS_Cutpoints_OS_Data_DF <- setDT(EFS_Cutpoints_OS_Data_DF, keep.rownames = "Cell Types")[]
getwd()
write_xlsx(EFS_Cutpoints_OS_Data_DF,"Hazard Ratios OS from EFS Cutpoints.xlsx")
