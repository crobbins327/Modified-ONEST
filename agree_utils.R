library(irr)
# library(ONEST)
library(combinat)
# library(openxlsx)
# data('sp142_bin')
# dataSP = sp142_bin
library(ggplot2)
library(reshape2)


opa = function(df){
  #Determine number of samples in df
  n = nrow(df)
  temp = rep(0,n)
  for (j in 1:n){
    comp = as.vector(unlist(df[j,]))
    #Find which samples have 1 unique category == the samples with complete agreement
    #NA counts as a unique value, and are not removed from calculation
    temp[j] = length(unique(comp))
  }
  #Find the overall percent agreement as total number with complete agreement over samples
  opaVal = sum(temp==1, na.rm = TRUE)/n 
  SE = sqrt(opaVal*(1-opaVal)/n)
  #95% CI calculation
  l_ci = opaVal - 1.96*SE
  u_ci = opaVal + 1.96*SE
  return(list(opval=opaVal, lower95=l_ci, upper95=u_ci))
}

categ_split = function(data, type='versus', f1 = 0, f2 = 1){
  #If f1 or f2 is a vector, then combine the values into a single group variable
  #e.g. f1 = c(1,2) then combine the variables into group "1,2"
  f1Name = f1
  f2Name = f2
  if (length(f1)>1){
    f1Name = paste0(f1, collapse=',')
    for(i in 1:length(f1)){
      data[data==f1[i]] = f1Name
    }
  }
  if (length(f2)>1){
    f2Name = paste0(f2, collapse=',')
    for(i in 1:length(f2)){
      data[data==f2[i]] = f2Name
    }
  }
  
  
  if (type=='only'){
    #Subset all cases that have been classified as f1
    includeV = rowSums(data==f1Name, na.rm = T) > 0
    sub = data[includeV, ]
    return(sub)
  }
  else if (type=='versus'){
    #If you want to compare 0 vs 1 or 0 vs 2, etc
    #Subset all cases that have been classified as f1 or f2
    includeV = rowSums(data==f1Name, na.rm = T) > 0 | rowSums(data==f2Name, na.rm = T) > 0
    sub = data[includeV, ]
    return(sub)
  } else if (type=='not'){
    #If you want to compare 0 vs not 0
    sub = data
    sub[sub!=f1Name]  = paste0('not',f1Name)
    return(sub)
  } else if (type=='combined_all'){
    #Takes lists of f1 and f2 and combines the factors into two separate groups
    #e.g. 0 vs 1,2(low) or 1,2(low) vs 3
    return(data)
  } else {
    cat(type,'not allowed for input variable type. Use either "only", "versus", "not", or "combined_all".')
  }
}

perCaseBar = function(data, name='Stacked barplot', file='Test',catLabs = NA, legendTitle = NA, C=100){
  #Replace NAs with a another label
  data[is.na(data)] = '_'
  #Get all categories
  cats = unique(as.vector(as.matrix(data)))
  cats = cats[order(cats)]
  #Remove NA category
  naI = grep('_',cats)
  if(length(naI)!=0){
    cats=cats[-naI]
  }
  if (is.na(catLabs)[1]){
    catLabs = cats
  }
  #Get number of raters
  m = ncol(data)
  #Remove samples with all NAs
  naCount = rowSums(data=='_')
  naI = which(naCount>=m)
  if(length(naI)!=0){
    data = data[-naI,]
  }
  #Calculate the percent observers for each category for each case
  catcos = c()
  catpers = c()
  #No repeated measurements by raters
  total = m
  for(j in 1:length(cats)){
    catco = paste0(as.character(cats[j]),' counts')
    catper = paste0(as.character(cats[j]),' percent')
    catcos[j] = catco
    catpers[j] = catper 
    data[,catco] = NA
    data[,catper] = NA
    for(i in 1:nrow(data)){
      data[i, catco] = sum(data[i,1:m]==cats[j])
      #Remove NA ratings from counts in calculating percentage
      naCount = sum(data[i,1:m]=='_')
      data[i, catper] = 100*data[i, catco]/(total-naCount)
    }
  }
  perOnly = data[,catpers]
  orderStr = "order("
  for (j in 1:length(catpers)){
    orderStr = paste0(orderStr,'perOnly[,',j,'], ')
  }
  orderStr = paste0(orderStr,'decreasing=T)')
  perOnly = perOnly[eval(parse(text=orderStr)),]
  perOnly$ID = rownames(perOnly)
  data$ID = rownames(data)
  data = data[match(perOnly$ID, data$ID),]
  perMelt = melt(perOnly[,c(catpers,"ID")], id.vars="ID", variable.name='cat_per')
  perMelt$order = 1:nrow(perMelt)
  # catLabs = c('Tumor', 'Stroma', 'None')
  catLabsP = catLabs[order(catpers, decreasing = T)]
  if(is.na(legendTitle)){
    legendTitle = 'Categories'
  }
  p = ggplot(data=perMelt, aes(x=reorder(ID,order), y=value, fill=factor(cat_per,levels=catpers[order(catpers, decreasing=T)]))) +
    geom_bar(position="stack", stat='identity', width=1) +
    scale_y_continuous(expand=c(0,0), breaks=seq(0,100,by=20), limits=c(0,100)) +
    ggtitle(name) +
    labs(fill=legendTitle) + 
    xlab('Cases') + ylab('% Observers') +
    scale_fill_discrete(labels = catLabsP, guide = guide_legend(reverse = TRUE)) +
    theme(text = element_text(size=14),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_text(size=14, color='black'),
          axis.title = element_text(size=16, color='black'),
          plot.title = element_text(size=16, color='black'),
          legend.title = element_text(color='black', size=16, face='bold'),
          legend.text = element_text(color='black', size=14),
          legend.position = 'bottom',
          axis.ticks.length = unit(1.4, 'mm'),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  p
  ggsave(paste0(file,'_stackedBarbyCases.jpg'), plot=p, device="jpeg", width=10, height=5, units="in", dpi=300)
  
  #Calculate number of concordant cases and discordant cases at each cutpoint
  #Calculate x % concordant cases for each category and store in vector
  concordN = c()
  for (j in 1:length(cats)){
    concordN[j] = sum(perOnly[,j]>=C)
    names(concordN)[j] = catLabs[j]
  }
  
  discordN = c()
  stopInd = 0
  for (j in 1:length(cats)-1){
    #Data frame is already ordered, 
    #Need to find the discordant cases between cats 1-2, 2-3, etc. 
    #Find the first instance where cat 2 == 100%
    startInd = stopInd + 1
    k=j+1
    stopInd = which(perOnly[, k]>=C)[1] - 1
    while(is.na(stopInd)&k<length(cats)){
      stopInd = which(perOnly[, k+1]>=C)[1] - 1
      k=k+1
    }
    if(!is.na(stopInd)&!is.na(startInd)){
      discordN[j] = sum(perOnly[startInd:stopInd, j]<C)
    }else{
      discordN[j] = NA
    }
    names(discordN)[j] = paste0(catLabs[j],'-',catLabs[k])
  }
  
  return(list(percent=data, concordN=concordN, discordN=discordN))
}

perPathBar = function(ratings, name='Percent of HER2 status assigned by each pathologist', file='Test',catLabs = NA, legendTitle = NA){
  ratings[is.na(ratings)] = '_'
  cats = unique(unlist(ratings))
  cats = cats[order(cats)]
  #Remove NA category
  naI = grep('_',cats)
  if(length(naI)!=0){
    cats=cats[-naI]
  }
  if (is.na(catLabs)[1]){
    catLabs = cats
  }
  #Get number of raters/pathologists
  m=ncol(ratings)
  #Remove samples with all NAs
  naCount = rowSums(ratings=='_')
  naI = which(naCount>=m)
  if(length(naI)!=0){
    ratings = ratings[-naI,]
  }
  #Sum the pathologists readings for all cases in each category
  for(i in 1:length(cats)){
    if(i==1){
      pBPath = data.frame(colSums(ratings==cats[i]))
    } else{
      pBPath[,i] = colSums(ratings==cats[i])
    }
  }
  colnames(pBPath) = cats
  # catcount = paste0(cats,'.counts')
  catLabsP = catLabs[order(cats, decreasing = T)]
  orderStr = "order("
  for (j in 1:length(cat)){
    orderStr = paste0(orderStr,'pBPath[,',j,'], ')
  }
  orderStr = paste0(orderStr,'decreasing=T)')
  pBPath = pBPath[eval(parse(text=orderStr)),]
  rownames(pBPath) = NULL
  # pBPath$raters = rownames(pBPath)
  pBPath$raterID = as.numeric(rownames(pBPath))
  #Convert counts to percent by dividing each row by case number
  pBPath[,1:length(cats)] = pBPath[,1:length(cats)]/rowSums(pBPath[,1:length(cats)])*100
  pMelt = melt(pBPath, id="raterID", variable.name = 'categ')
  pMelt$order = 1:nrow(pMelt)
  if(is.na(legendTitle)){
    legendTitle = 'Categories'
  }
  p = ggplot(data=pMelt, aes(x=reorder(raterID,order), y=value, fill=factor(categ,levels=cats[order(cats, decreasing=T)]))) +
    geom_bar(position="stack", stat='identity', width=1) +
    scale_y_continuous(expand=c(0,0), breaks=seq(0,100,by=20), limits=c(0,101)) +
    ggtitle(name) +
    labs(fill=legendTitle) + 
    xlab('Pathologist') + ylab('% Cases') +
    scale_fill_discrete(labels = catLabsP, guide = guide_legend(reverse = TRUE)) +
    theme(text = element_text(size=14),
          axis.text = element_text(size=14, color='black'),
          axis.title = element_text(size=16, color='black'),
          plot.title = element_text(size=16, color='black'),
          legend.title = element_text(color='black', size=16, face='bold'),
          legend.text = element_text(color='black', size=14),
          legend.position = 'bottom',
          axis.ticks.length = unit(1.4, 'mm'),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  p
  ggsave(paste0(file,'_stackedBarbyPath.jpg'), plot=p, device="jpeg", width=10, height=5, units="in", dpi=300)
  return(pathCount=pBPath)
}


fillTable = function(df, metTab, groupN, noICC=F){
  opvect = unlist(opa(df))
  opav = opvect[1]
  l_ci = opvect[2]
  u_ci = opvect[3]
  metTab[groupN, 'OPA'] = paste0(round(opav,4)*100,"(", round(l_ci,4)*100,",",round(u_ci,4)*100,")")
  # kap2 = kappam.fleiss(her2, detail=T)
  conT = concTab(df)
  kap = concordance(conT, test='Normal')
  metTab[groupN, 'Fkappa'] = paste0(round(unlist(kap$Fleiss[1]),3),"(", round(unlist(kap$Fleiss[2]),3),",", round(unlist(kap$Fleiss[3]),3),")")
  if(noICC){
    metTab[groupN, 'ICC'] = '--'
  }else{
    iccR = icc(df, model = "twoway", type = "agreement", unit = "single")
    metTab[groupN, 'ICC'] = paste0(round(unlist(iccR[7]), 3),"(", round(unlist(iccR[14]), 3),",", round(unlist(iccR[15]), 3),")")
  }
  return(metTab)
}

spAgree = function(ratings, categories=c(), weighting='identity' ){
  #Determine unique categories if none are supplied...
  ratings[is.na(ratings)] = '_'
  x = unique(unlist(ratings))
  x = x[order(x)]
  naI = grep('_',x)
  if(length(naI)!=0){
    x=x[-naI]
  }
  if(length(categories)==0){
    categories = x
  }
  #Get number of raters/pathologists
  r=ncol(ratings)
  #Remove samples with all NAs
  naCount = rowSums(ratings=='_')
  naI = which(naCount>=r)
  if(length(naI)!=0){
    data = ratings[-naI,]
  }
  
  n = nrow(ratings)
  
  q = length(categories)
  #Get weights from options
  if(weighting=='identity'){
    weights = diag(q)
  }else if(weighting=='linear'){
    
  }else if(weighting=='quadratic'){
    
  }else{
    cat("ERROR:", weighting, "is not an appropriate weighting option.")
  }
  #Create n x q matrix (rater counts in item by category matrix)
  r_ik = matrix(nrow=n, ncol=q)
  for(i in 1:q){
    r_ik[,i] = rowSums(ratings==categories[i])
  }
  #Weight n x q matrix based on weighting scheme
  rstar_ik = t(weights %*% t(r_ik))
  #Calculate category-specific agreement
  r_i = r_ik %*% matrix(1, nrow=q, ncol=1)
  r_i = matrix(rep(r_i, each=q), ncol=q)
  numerator = colSums(r_ik * (rstar_ik - 1))
  denominator = colSums(r_ik * (r_i-1))
  SA = numerator / denominator
  #Average observed agreement, sum of rater-rater pairs that agree over all possible rater-rater pairs across samples
  A = sum(numerator)/(n*r*(r-1))
  return(list(categories=categories, specific=SA, mean_observed=A))
}
