# library(raters)
library(irr)
# library(ONEST)
library(combinat)
library(openxlsx)
# data('sp142_bin')
# dataSP = sp142_bin
library(ggplot2)
library(reshape2)


#######################################################################################
#######################################################################################
#Backend functions for ONEST technique
#######################################################################################
#######################################################################################

#Get functions. Change directories to where the R files are located.
source('D:/Documents/Rimm/CPS-Gastric_discordance/ONEST_core.R')
source('D:/Documents/Rimm/CPS-Gastric_discordance/agree_utils.R')

#######################################################################################
#######################################################################################
#Making ONEST plots from pathologist data
#######################################################################################
#######################################################################################

#Load data
setwd('D:/Documents/Rimm/CPS-Gastric_discordance')
score_sheet = which(grepl('\\CPS Score <>1\\b', getSheetNames("Hub-ALL_Gastric_CPS.ScoreComparison-02-02-22.xlsx")))
ptiles_sheet = which(grepl('\\CPS Ventiles\\b', getSheetNames("Hub-ALL_Gastric_CPS.ScoreComparison-02-02-22.xlsx")))
data_score = read.xlsx(xlsxFile = "Hub-ALL_Gastric_CPS.ScoreComparison-02-02-22.xlsx", sheet=score_sheet, fillMergedCells = TRUE, colNames = TRUE)
data_ptiles =  read.xlsx(xlsxFile = "Hub-ALL_Gastric_CPS.ScoreComparison-02-02-22.xlsx", sheet=ptiles_sheet, fillMergedCells = TRUE, colNames = TRUE)

#Replace categories with NA
data_score[data_score=='NA NO H&E'] = NA
data_score[data_score=='indeterminate'] = NA

data_ptiles[data_ptiles=='NA NO H&E'] = NA
data_ptiles[data_ptiles=='indeterminate'] = NA
data_ptiles[data_ptiles=='<1'] = 0

#Remove all NA columns
data_score = data_score[,colSums(is.na(data_score))<nrow(data_score)]
data_ptiles = data_ptiles[,colSums(is.na(data_ptiles))<nrow(data_ptiles)]

#Remove rows that contain more than 50% NA columns
data_score = data_score[rowSums(is.na(data_score))<as.integer(ncol(data_score)*0.5),]
data_ptiles = data_ptiles[rowSums(is.na(data_ptiles))<as.integer(ncol(data_ptiles)*0.5),]

#Make rownames the case names
rownames(data_score) = data_score$`Case#`
data_score$`Case#` = NULL

rownames(data_ptiles) = data_ptiles$`Case#`
data_ptiles$`Case#` = NULL

#Clean up categories
#CPS <1 or >=1
data_score[data_score == 0] = '<1'
data_score[data_score >= 1] = '>=1'
data_score[data_score == '>1'] = '>=1'

#Percentiles
data_ptiles[data_ptiles=='10 or 5'] = 10
#Replace NAs with data_score values
data_ptiles[is.na(data_ptiles)] = data_score[is.na(data_ptiles)]
data_ptiles[data_ptiles=='<1'] = 0
data_ptiles[data_ptiles=='>=1'] = 1


#######################################################################################
#######################################################################################
#Generate ONEST data for CPS score groupings
#######################################################################################
#######################################################################################

dir.create(file.path(getwd(),'ONEST_concord'))

#ONEST plots for 2 category CPS score (<1, >=1)
results = ONEST_plot(data_score, plotI=T, metric = 'OPA')
write.table(results$concord, 'ONEST_concord/CPS-gastric_concordance-opa.txt', sep='\t', row.names = F, quote = F)
write.table(results$stats, 'ONEST_concord/CPS-gastric_plotStats-opa.txt', sep='\t', row.names = F, quote = F)

consist = data.frame(results$modelData$consistency)
consist$path_number = 2:(nrow(consist)+1)
diff = data.frame(results$modelData$difference)
diff$path_number = 2:(nrow(consist))
mData = merge(consist, diff, by = 'path_number', all = TRUE)
write.table(mData, paste0('ONEST_concord/','CPS-gastric_modelData-opa.txt'), sep='\t', row.names = F, quote = F)


#ICC ONEST plot for percentile CPS score
results = ONEST_plot(data_ptiles, plotI=T, metric = 'icc')
write.table(results$concord, 'ONEST_concord/perCPS-gastric_concordance-icc.txt', sep='\t', row.names = F, quote = F)
write.table(results$stats, 'ONEST_concord/perCPS-gastric_plotStats-icc.txt', sep='\t', row.names = F, quote = F)

#Histogram of percentile scores, to view where the relative cutpoints are
# new_data=na.omit(as.vector(data.matrix(data_ptiles)))
# hist(new_data)
# d = density(new_data)
# plot(d)

# #Group by x<1, 1=< x <5, 5=< x 10,... every 5, every 10, every 20
# cats = c('<1', '>=1', 'every20')
# data_ptiles = data.frame(lapply(data_ptiles, as.numeric))
# 
# #Need to make data_backup. Changing the data from numeric to characters while evaluating
# #inequalities messes up the indexing for dataframe
# 
# data_backup = data_ptiles
# for(i in 1:length(cats)){
#   if(grepl('every',cats[i])){
#     every = as.integer(gsub('every','',cats[i]))
#     iter = as.integer(100/every)
#     start = 1
#     for(j in 0:iter){
#       if(j==0){
#         ineq = paste0('data_backup>',start,' & data_backup<=',every*(j+1))
#         group_label = paste0(start,'< x <=',every*(j+1))
#       } else { 
#         ineq = paste0('data_backup>',every*j,' & data_backup<=',every*(j+1))
#         group_label = paste0(every*j,'< x <=',every*(j+1))
#       }
#       
#       data_ptiles[eval(parse(text=ineq))] = group_label
#     }
#   } else {
#     ineq = paste0('data_backup',cats[i])
#     data_ptiles[eval(parse(text=ineq))] = cats[i]
#   }
#   
# }
# 
# results = ONEST_plot(data_ptiles, plotI=T, metric = 'OPA')
# write.table(results$concord, 'ONEST_concord/20perCPS-gastric_concordance-opa.txt', sep='\t', row.names = F, quote = F)
# write.table(results$stats, 'ONEST_concord/20perCPS-gastric_plotStats-opa.txt', sep='\t', row.names = F, quote = F)
# data_ptiles = data_backup


#LGper
#<1, 1-10, 10-20, 20-50, and 50-100
cats = c('>=0-<1', '>=1-<=10', '>10-<=20', '>20-<=50', '>50-<=100')
#<10, >10
cats = c('<10', '>=10')
#<20, >20
cats = c('<20', '>=20')
#SMper
#<1,1-20, >20
cats = c('<1', '>=1-<=20', '>20')

data_backup = data_ptiles
for(i in 1:length(cats)){
  if(grepl('-',cats[i])){
    bounds = strsplit(cats[i],'-')
    b_low = bounds[[1]][1]
    b_high = bounds[[1]][2]
    ineq = paste0('data_backup',b_low,' & data_backup',b_high)
    group_label = cats[i]
    data_ptiles[eval(parse(text=ineq))] = group_label
  } else {
    ineq = paste0('data_backup',cats[i])
    data_ptiles[eval(parse(text=ineq))] = cats[i]
  }
}

results = ONEST_plot(data_ptiles, plotI=T, metric = 'OPA')
write.table(results$concord, 'ONEST_concord/lt20perCPS-gastric_concordance-opa.txt', sep='\t', row.names = F, quote = F)
write.table(results$stats, 'ONEST_concord/lt20perCPS-gastric_plotStats-opa.txt', sep='\t', row.names = F, quote = F)
data_ptiles = data_backup

consist = data.frame(results$modelData$consistency)
consist$path_number = 2:(nrow(consist)+1)
diff = data.frame(results$modelData$difference)
diff$path_number = 2:(nrow(consist))
mData = merge(consist, diff, by = 'path_number', all = TRUE)
write.table(mData, paste0('ONEST_concord/','lt20perCPS-gastric_modelData-opa.txt'), sep='\t', row.names = F, quote = F)

#ONEST plots of cases with score of X only (by at least one rater)
# onlyComp = c(0, "1,2", 3)
# onlyComp = c(0,1,2,3)
# for(i in 1:length(onlyComp)){
#   f1 = onlyComp[i]
#   cat('Working on',f1,'only...','\n')
#   sub = categ_split(her2, type='only', f1 = f1)
#   results = ONEST_plot(sub, plotI=T, metric = 'OPA')
#   write.table(results$concord, paste0('ONEST_concord/HER2/',f1,'only_HER2_concordance-opa.txt'), sep='\t', row.names = F, quote = F)
#   write.table(results$stats, paste0('ONEST_concord/HER2/',f1,'only_HER2_plotStats-opa.txt'), sep='\t', row.names = F, quote = F)
# }

#ONEST plots of scores grouped by X vs not X
# score = c(0,1,2,3)
# score = c(0,'1,2',3)
# score = c(1,2)
# for(i in 1:length(score)){
#   f1 = score[i]
#   cat('Working on',f1,'vs not',f1,'...','\n')
#   sub = categ_split(her2low, type='not', f1 = f1)
#   results = ONEST_plot(sub, plotI=T, metric = 'OPA')
#   # write.table(results$concord, paste0('ONEST_concord/HER2/',f1,'vnot',f1,'_HER2_concordance-opa.txt'), sep='\t', row.names = F, quote = F)
#   # write.table(results$stats, paste0('ONEST_concord/HER2/',f1,'vnot',f1,'_HER2_plotStats-opa.txt'), sep='\t', row.names = F, quote = F)
#   
#   consist = data.frame(results$modelData$consistency)
#   consist$path_number = 2:(nrow(consist)+1)
#   diff = data.frame(results$modelData$difference)
#   diff$path_number = 2:(nrow(consist))
#   mData = merge(consist, diff, by = 'path_number', all = TRUE)
#   write.table(mData, paste0('ONEST_concord/HER2/',f1,'vnot',f1,'_HER2_modelData-opa.txt'), sep='\t', row.names = F, quote = F)
# }

#######################################################################################
#######################################################################################
#Plotting CPS ONEST plots
#######################################################################################
#######################################################################################

desiredPlots = c('CPS-gastric',
                 'LGperCPS-gastric',
                 'SMperCPS-gastric',
                 'lt10perCPS-gastric',
                 'lt20perCPS-gastric')

labels = c('CPS score (<1, >=1)', 
           'CPS score (<1, 1-10, 10-20, 20-50, and 50-100)',
           'CPS score (<1, 1-20, >20)',
           'CPS score (<10, >=10)',
           'CPS score (<20, >=20)')

desiredPlots = c('perCPS-gastric')
labels = c('Percent CPS score')

# datadir = 'ONEST_concord/Fleiss_Kappa'
# datadir = 'ONEST_concord/OPA'
datadir = 'ONEST_concord/ICC'

flist = list.files(datadir)

for (i in 1:length(desiredPlots)){
  #Find files with similar matching names
  print(paste0('Working on ', desiredPlots[i],' plots...'))
  loadFiles = flist[grep(paste0("^",desiredPlots[i]), flist)]
  concord = read.delim(paste0(datadir,'/',loadFiles[grep('concord',loadFiles)]), sep='\t')
  plot_data = read.delim(paste0(datadir,'/',loadFiles[grep('plotStats',loadFiles)]), sep='\t')
  ONEST_plot_fromData(concord, plot_data, name=labels[i], file=desiredPlots[i], ylab="ICC", color='blue', percent = F)
  # ONEST_plot_fromData(concord, plot_data, name=labels[i], file=desiredPlots[i], ylab="Fleiss' Kappa", color='red', percent = F)
  # ONEST_plot_fromData(concord, plot_data, name=labels[i], file=desiredPlots[i], ylab="Overall Percent Agreement", color='black', percent = T)
  if(length(loadFiles)==3){
    model_data = read.delim(paste0(datadir,'/',loadFiles[grep('modelData',loadFiles)]), sep='\t')
    ONEST_plotModel_fromData(model_data, name=labels[i], file=desiredPlots[i], percent=TRUE)
  }
}

#######################################################################################
#######################################################################################
#Plotting CPS gastric stacked bar graphs
#######################################################################################
#######################################################################################

staRes = perCaseBar(data_score, name='',file='CPS-gastric', legendTitle = 'CPS score', C=100)
staResPath = perPathBar(data_score, name='Percent of gastric CPS score assigned by each pathologist',file='CPS-gastric', legendTitle = 'CPS score')


#######################################################################################
#######################################################################################
#Creating table of inter-rater reliability metrics for HER2 IHC
#######################################################################################
#######################################################################################

group = c('4 cat', '3 cat Low',
          '0 only', '1 only', '2 only', '3 only', 'Low only',
          '0 vs not 0', 'Low vs not Low', '3 vs not 3')

her2 = data.frame(lapply(her2, as.numeric))
her2low = categ_split(her2, type='combined_all', f1 = 0, f2 = c(1,2))

metTab = data.frame(group)
metTab$OPA = NA
metTab$Fkappa = NA
metTab$ICC = NA
rownames(metTab) = metTab$group
library(raters)

metTab = fillTable(her2, metTab, "4 cat", noICC=F)
mher2low=her2low
mher2low[mher2low=='1,2']=1
mher2low[mher2low==3]=2
mher2low = data.frame(lapply(mher2low, as.numeric))
metTab = fillTable(mher2low, metTab, "3 cat Low", noICC=F)

onlyComp = c(0,1,2,3,'1,2')
for(i in 1:length(onlyComp)){
  f1 = onlyComp[i]
  cat('Working on',f1,'only...','\n')
  if(f1=='1,2'){
    sub = categ_split(mher2low, type='only', f1 = 1)
    metTab = fillTable(sub, metTab, 'Low only', noICC=F)
  }else{
    sub = categ_split(her2, type='only', f1 = f1)
    metTab = fillTable(sub, metTab, paste0(f1,' only'), noICC=F)
  }
}

score = c(0,'1,2',3)
for(i in 1:length(score)){
  f1 = score[i]
  if(f1=='1,2'){
    groupN = 'Low vs not Low'
  } else{
    groupN = paste0(f1, ' vs not ', f1)
  }
  cat('Working on',groupN,'...','\n')
  sub = categ_split(her2low, type='not', f1 = f1)
  sub[sub=='not0'] = 1
  sub[sub=='0'] = 0
  sub[sub=='not3'] = 0
  sub[sub=='3'] = 1
  if(f1!='1,2'){
    sub = data.frame(lapply(sub, as.numeric))
  }
  metTab = fillTable(sub, metTab, groupN, noICC=F)
}

write.table(metTab, 'HER2-discordance-metrics.txt', sep='\t', row.names=F, quote=F)