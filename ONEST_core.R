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

plotline = function(path, indi=1, color='red', metric='OPA'){
  n = nrow(path)
  m = ncol(path)
  a = rep(0,m-1)
  
  # Concordance for original data
  concord  = matrix(a,ncol=1)
  M = m-1
  if(metric=='OPA'){
    for (i in 1:M){
      I = i+1
      concord[i] = unlist(opa(path[,1:I]))[1]
    }
  }else if(metric=='fkappa'){
    for (i in 1:M){
      I = i+1
      concord[i] = unlist(kappam.fleiss(path[,1:I])[5])
    }
  }else if(metric=='icc'){
    #only works for continuous or ordinal data
    #will not work on subcategories?
    path = data.frame(lapply(path, as.numeric))
    for (i in 1:M){
      I = i+1
      iccR = icc(path[,1:I], model = "twoway", type = "agreement", unit = "single")
      concord[i] = as.numeric(unlist(iccR)[7])
    }
  }else{
    cat('ERROR: Metric',metric,'is not acceptable.')
    break
  }
  
  # Plot the data based on the indi status
  if (indi == 1){
    x_axis = 2:m
    plot(x_axis,concord, type = "o", col = color,lwd=1.8,xlim=c(2,m),ylim=c(0,1),
         xlab = 'Number of raters',ylab='Proportion of identical reads',main="Figure(1)")
  }
  
  C = concord
  return(C)
}

ONEST_plotModel_fromData = function(model_data, name, file, percent=TRUE, xBy=2){
  m = nrow(model_data)+1
  scale=1
  cons_ylab = 'Proportion of identitcal reads'
  if (percent){
    #don't multiply path number
    model_data[-1] = model_data[-1]*100
    scale=100
    cons_ylab = 'Overall Percent Agreement'
  }
  #should have path number
  consist = model_data[,c('path_number','consist_p', 'consist_low')]
  
  line_types = c("Model"=1,"95% lower bound"=2)
  color_codes = c("Model"='red',"95% lower bound"='blue') 
  q = ggplot(data = consist) + 
    geom_line(aes(x=path_number, y=consist_low, linetype="95% lower bound", color="95% lower bound"), size=1) + 
    geom_line(aes(x=path_number, y=consist_p, linetype="Model", color="Model"), size=1) + 
    geom_point(aes(x=path_number, y=consist_p), color='red', size=3, pch=21) + 
    scale_x_continuous(breaks= seq(2,m,by=xBy), limits= c(2,NA)) +
    xlab('Number of Raters in Group') + ylab(cons_ylab) + ggtitle(name) +
    labs(color='') +
    scale_linetype_manual(values=line_types) +
    scale_color_manual(values=color_codes, guide = guide_legend(reverse = TRUE)) +
    guides(linetype="none")+
    theme_linedraw()
  if (consist$consist_p[nrow(consist)]>0.75*scale){
    #Add legend below
    q = q + 
      theme(text=element_text(size=14),
            axis.text = element_text(size=14, color='black'),
            axis.title = element_text(size=16, color='black'),
            plot.title = element_text(size=16, color='black'), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = c(.82, .2),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.spacing.y = unit(0.1,'cm'),
            legend.key.size = unit(0.8, 'cm'),
            axis.ticks.length = unit(1.4, "mm"))
  } else{
    #Add legend above
    q = q + 
      theme(text=element_text(size=14),
            axis.text = element_text(size=14, color='black'),
            axis.title = element_text(size=16, color='black'),
            plot.title = element_text(size=16, color='black'), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = c(.82, .9),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.spacing.y = unit(0.1,'cm'),
            legend.key.size = unit(0.8, 'cm'),
            axis.ticks.length = unit(1.4, "mm"))
  }
  
  if (percent){
    q = q + scale_y_continuous(expand=c(0,0), breaks=seq(0,100,by=20), limits=c(0,100))
  } else{
    q = q + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,by=0.2), limits=c(0,1))
  }
  ggsave(paste0(file,'_model.jpg'), plot=q, device="jpeg", width=6.43, height=5.18, units="in", dpi=300)
  
  #should have path number
  diff = model_data[,c('path_number','diff_consist', 'diff_high')]
  diff = na.omit(diff)
  
  line_types = c("Difference trajectory"=1,"95% upper bound"=2)
  color_codes = c("Difference trajectory"='red',"95% upper bound"='blue') 
  p = ggplot(data = diff) + 
    geom_line(aes(x=path_number+1, y=-diff_consist, linetype="Difference trajectory", color="Difference trajectory"), size=1) + 
    geom_point(aes(x=path_number+1, y=-diff_consist), color='red', size=3, pch=21) + 
    geom_line(aes(x=path_number+1, y=diff_high, linetype="95% upper bound", color="95% upper bound"), size=1) + 
    scale_x_continuous(breaks= seq(3,m,by=xBy), limits=c(3,NA)) +
    xlab('Number of Raters in Group') + ylab('OPA Difference between points') + ggtitle(name) +
    labs(color='') +
    scale_linetype_manual(values=line_types) +
    scale_color_manual(values=color_codes, guide = guide_legend(reverse = TRUE)) +
    guides(linetype="none")+
    theme_linedraw()
  # p
  if (diff$diff_high[nrow(diff)]>0.75*scale){
    #Add legend below
    p = p + 
      theme(text=element_text(size=14),
            axis.text = element_text(size=14, color='black'),
            axis.title = element_text(size=16, color='black'),
            plot.title = element_text(size=16, color='black'), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = c(.8, .2),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.spacing.y = unit(0.1,'cm'),
            legend.key.size = unit(0.8, 'cm'),
            axis.ticks.length = unit(1.4, "mm"))
  } else{
    #Add legend above
    p = p + 
      theme(text=element_text(size=14),
            axis.text = element_text(size=14, color='black'),
            axis.title = element_text(size=16, color='black'),
            plot.title = element_text(size=16, color='black'), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = c(.8, .9),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.spacing.y = unit(0.1,'cm'),
            legend.key.size = unit(0.8, 'cm'),
            axis.ticks.length = unit(1.4, "mm"))
  }
  
  if (percent){
    # mval = ceiling(max(diff$diff_high)*1.1)
    mval = max(diff$diff_high)*1.1
    p = p + scale_y_continuous(expand=c(0,0), breaks=seq(0,mval,by=1), limits=c(0,mval))
  } else{
    mval = max(diff$diff_high)*1.1
    p = p + scale_y_continuous(expand=c(0,0), breaks=seq(0,mval,by=0.01), limits=c(0,mval))
  }
  p
  ggsave(paste0(file,'_diff-trajectory.jpg'), plot=p, device="jpeg", width=6.43, height=5.18, units="in", dpi=300)
}


ONEST_plot_fromData = function(concord, plot_data, name, file, ylab = 'Overall Percent Agreement', color='black', percent=TRUE){
  m = nrow(concord)+1
  scale=1
  if (percent){
    concord = concord*100
    plot_data = plot_data*100
    scale=100
  }
  concord$path_number = 2:m
  df = melt(concord[,c(1:100,ncol(concord))], id.vars="path_number", variable.name='iter')
  p = ggplot(data = df, aes(x=path_number, y=value, group=iter)) + geom_line(color=color, alpha=0.4, size=0.5) + geom_point(color=color, size=3, pch=21) + 
    scale_x_continuous(breaks= seq(2,m,by=2), limits=c(2,NA)) +
    xlab('Number of Raters in Group') + ylab(ylab) + ggtitle(name) +
    theme_linedraw() +
    theme(text=element_text(size=14),
          axis.text = element_text(size=14, color='black'),
          axis.title = element_text(size=16, color='black'),
          plot.title = element_text(size=16, color='black'),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position = "none",
          axis.ticks.length = unit(1.4, "mm"))
  if (percent){
    p = p + scale_y_continuous(expand=c(0,0), breaks=seq(0,100,by=20), limits=c(0,100))
  } else{
    p = p + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,by=0.2), limits=c(0,1))
  }
  # p
  ggsave(paste0(file,'_100sim.jpg'), plot=p, device="jpeg", width=6.43, height=5.18, units="in", dpi=300)
  
  line_types = c("Mean"=1,"95% CI"=2)
  plot_data[plot_data < 0] = 0
  q = ggplot(data = plot_data) + 
    geom_line(aes(x=2:m, y=mean, linetype="Mean"), color=color, size=1) + 
    geom_point(aes(x=2:m, y=mean), color=color, size=3, pch=21) + 
    geom_line(aes(x=2:m, y=quant_2.5, linetype="95% CI"), color=color, size=1) + 
    geom_line(aes(x=2:m, y=quant_97.5, linetype="95% CI"), color=color, size=1) + 
    scale_x_continuous(breaks= seq(2,m,by=2), limits=c(2,NA)) +
    xlab('Number of Raters in Group') + ylab(ylab) + ggtitle(name) +
    scale_linetype_manual(name='', values=line_types) +
    theme_linedraw()
  if (plot_data$quant_97.5[nrow(plot_data)]>0.75*scale){
    q = q + 
      theme(text=element_text(size=14),
            axis.text = element_text(size=14, color='black'),
            axis.title = element_text(size=16, color='black'),
            plot.title = element_text(size=16, color='black'), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = c(.9, .2),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.spacing.y = unit(0.1,'cm'),
            legend.key.size = unit(0.8, 'cm'),
            axis.ticks.length = unit(1.4, "mm"))
  } else{
    q = q + 
      theme(text=element_text(size=14),
            axis.text = element_text(size=14, color='black'),
            axis.title = element_text(size=16, color='black'),
            plot.title = element_text(size=16, color='black'), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = c(.9, .9),
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.spacing.y = unit(0.1,'cm'),
            legend.key.size = unit(0.8, 'cm'),
            axis.ticks.length = unit(1.4, "mm"))
  }
  
  if (percent){
    q = q + scale_y_continuous(expand=c(0,0), breaks=seq(0,100,by=20), limits=c(0,100))
  } else{
    q = q + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,by=0.2), limits=c(0,1))
  }
  # q
  ggsave(paste0(file,'_meanCI.jpg'), plot=q, device="jpeg", width=6.43, height=5.18, units="in", dpi=300)
  # dev.off()
}

ONEST_plot = function(data, plotI=TRUE, color='black', metric = 'OPA', skipModeling=FALSE){
  # Size of the data
  n = nrow(data)
  m = ncol(data)
  size_case = n
  size_rater = m
  M = m-1
  
  ## Step 1, plot the data
  # 1) Plot the agreement percentage in the order of columns in the input
  if (plotI){
    C = plotline(data,1,color,metric)
    lines(0,0,main="Figure(1)")
  }
  
  
  # 2) Plot 100 randomly chosen permutations
  print('Simluating 100 combinations...')
  concord = matrix(rep(0,100*M),ncol=100)
  for (i in 1:100) {
    concord[,i] = plotline(data[,sample(m)],0,color,metric)
  }
  
  if (plotI){
    x_axis = 2:m
    plot(x_axis,concord[,100], type = "o", col = color,lwd=1.8,xlim=c(2,m),ylim=c(0,1),main="Figure(2)",
         xlab = 'Number of raters',ylab='Proportion of agreement')
    for (i in 1:99){
      lines(x_axis,concord[,i], type = "o", col = color,lwd=1.8,xlim=c(2,m),ylim=c(0,1))
    }
  }
  
  # 3) Plot the empirical confidence interval
  print('Simluating 1000 combinations...')
  simN = 1000
  D = matrix(rep(0,simN*(m-1)),ncol=simN)
  D[,1:100]=concord
  for (i in 101:simN){
    D[,i]  = plotline(data[,sample(m)],0,color,metric)
  }
  
  print('Calculating plotting descriptive statistics...')
  # Calculate mean, 5 th, and 95 th percentiles
  plot_data = matrix(rep(0,3*(m-1)),ncol=3)
  for (j in 1:M){
    plot_data[j,1] = mean(D[j,])
    plot_data[j,2] = quantile(D[j,],0.025,names=FALSE,type=1)
    plot_data[j,3] = quantile(D[j,],0.975,names=FALSE,type=1)
  }
  colnames(plot_data) = c('mean','quant_2.5','quant_97.5')
  # plot the mean and empirical 95# CI
  if (plotI){
    x_axis = 2:m
    plot(x_axis,plot_data[,1], type = "o", col = color,lwd=1.8,xlim=c(2,m),ylim=c(0,1),
         xlab = 'Number of raters',ylab='Proportion of agreement',main="Figure(3)")
    lines(x_axis,plot_data[,2],lty=2, col = color,lwd=1.8,xlim=c(2,m),ylim=c(0,1))
    lines(x_axis,plot_data[,3],lty=2, col = color,lwd=1.8,xlim=c(2,m),ylim=c(0,1))
    legend("topright", c("Mean","95% empirical CI"),lty=1:2,pch = c(1,NA),cex = 0.6)
  }
  
  ## Step 2, model the data;
  
  # Input data = sp142_bin;
  
  #Can only model the data if it has 2 categories....
  cats = unique(as.vector(as.matrix(data)))
  cats = na.omit(cats)
  
  if (skipModeling){
    print('Skipping modeling')
    return(list(concord=D, stats=plot_data, modelData=c()))
  } else if(length(cats) != 2){
    print(paste('Cannot model data with', length(cats),'categories...'))
    print('Returning ONEST plot data without model data.')
    return(list(concord=D, stats=plot_data, modelData=c()))
  } else if(length(cats) == 2 & metric=='OPA'){
    print('Preparing to model data of 2 categories...')
  } else{
    print('Exception encountered... cannot model data. Cannot model ICC or Fleiss Kappa data. Check number of categories of input dataframe.')
    return(list(concord=D, stats=plot_data, modelData=c()))
  }
  
  #Convert categories to 0 and 1
  #Prefer to keep ordinal scale where 0 is negative and 1 is positive
  #"not" is considered as a 0 and is negative
  
  #Get the factors from cats and reorder them so that "not" and lower numbers always start first
  #Assumes ordinal scale is based on ordering of categories
  cats = cats[order(cats)]
  #Find if there is a "not", then reorder cats vector so that "not" is first
  if(sum(grepl('not',cats))>0){
    cats=cats[order(grepl('not',cats), decreasing=TRUE)]
  }
  
  #For each column in dataframe, convert it into factors and then numeric - 1 to create 0 and 1 values
  data_temp = data
  data_temp = data.frame(lapply(data_temp, function(x) as.numeric(factor(x, levels=cats))-1))
  rownames(data_temp) = rownames(data)
  
  # Total number of observations;
  N = n*m
  # count the number of nan;
  nan_num = 0
  for (i in 1:n){
    for (j in 1:m){
      if (is.na(data_temp[i,j])==1){
        nan_num = nan_num + 1
      }
    }
  }
  
  # total number of positve reads;
  K = sum(data_temp==1,na.rm = TRUE)
  # calculate p_plus and p_minus, the proportions of all positive and negative;
  p_plus  = 0; p_minus = 0;
  for (i in 1:n){
    if (sum(data_temp[i,]==1,na.rm = TRUE)==0){
      p_minus = p_minus + 1/n
    }
    
    if (sum(data_temp[i,]==1,na.rm = TRUE)==sum(data_temp[i,],na.rm = TRUE)+sum(data_temp[i,]==0,na.rm = TRUE)){
      p_plus = p_plus + 1/n
    }
  }
  # [7/68 19/68]
  # [p_minus p_plus]
  # p_plus  = 19/68; p_minus = 7/68;
  
  # Calculate the probability 'p'
  prop = (K/(N-nan_num)-p_plus)/(1-p_plus-p_minus)
  
  # Calculate the consistancy probability trajectory;
  consist_p = rep(0,m-1)
  diff_consist = rep(0,m-2)
  o=m-2
  for (i in 1:M){
    consist_p[i] = (1-p_plus-p_minus)*(prop^(i+1)+(1-prop)^(i+1))+p_plus+p_minus
  }
  # Calculate the difference trajectory;
  for (i in 1:o){
    diff_consist[i] = consist_p[i+1]-consist_p[i]
  }
  
  ## Construct the upper bound for the difference using the lower bound for
  # p_plus+p_minus;
  p_c = p_plus+p_minus
  p_c_low = p_c - 1.645*sqrt(p_c*(1-p_c)/n)
  
  # Calculate the trajectory;
  consist_low = rep(0,m-1)
  diff_low = rep(0,m-2)
  for (i in 1:M){
    consist_low[i]=(1-p_c_low)*(prop^(i+1)+(1-prop)^(i+1))+p_c_low
  }
  
  for (i in 1:o){
    diff_low[i] = consist_low[i+1]-consist_low[i]
  }
  
  if (plotI){
    # figure(5)
    plot(x_axis,consist_p, type = "o", col = "red",lwd=1.8,ylim = c(0,1),pch=2,xlab = 'Number of raters',ylab='Proportion of identical reads',main="Figure(5)")
    lines(x_axis,consist_low,lty=2,lwd=1.8,col = "blue")
    # figure(6)
    plot(1:o,-diff_consist, type = "o", col = "red",lwd=1.8,pch=2,xlab = 'Number of raters',ylab='Difference between the consist_p',main="Figure(6)")
    lines(1:o,-diff_low,lty=2,lwd=1.8,col = "blue")
  }
  
  diff_high = -diff_low
  p = prop
  lower_bound = plot_data[,2]
  upper_bound = plot_data[,3]
  mean = plot_data[,1]
  empirical = rbind(cbind(lower_bound,mean,upper_bound))
  consist = rbind(cbind(consist_p,consist_low))
  diff = rbind(cbind(diff_consist,diff_high))
  estimate = rbind(cbind(size_case,size_rater,p,p_plus,p_minus))
  modelResults = list(consistency=consist,difference=diff,estimates=estimate,empirical=empirical)
  
  return(list(concord=D, stats=plot_data, modelData=modelResults))
}