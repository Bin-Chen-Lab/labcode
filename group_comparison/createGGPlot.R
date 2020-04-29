#Author: Bin Chen (2015)
#the code was modified from Purvesh Khatri's code
#give two groups of data, compute their difference and visualize using boxplot

createGGPlot <- function(currPheno, measureVar, groupVars, main="GGPLot", ylab="pass ylab value", xlab="pass xlab value", alternative="less") {
  library("ggplot2")
  currPheno$group = as.factor(currPheno$group)
  temp = summarySE(currPheno, measurevar=measureVar, groupvars=groupVars)[, c(groupVars, measureVar, "se","max_group", "max_all")]
  temp_line = createLines(temp)
  temp_text = createText(currPheno, temp, alternative)
  a = ggplot() +
    geom_violin(data=currPheno, aes(x=group,y=score), fill='grey',trim=F) +
    geom_jitter(data=currPheno, aes(x=group,y=score, color=color), size = 2, position=position_jitter(width=jitterWidth, height=jitterHight)) +
    guides(color=FALSE) +
    xlab(xlab) +
    ylab(ylab) + 
    ggtitle(main) +
    theme(text = element_text(size=24)) +
    #theme(legend.position=c(0.5,0.9), legend.title=element_blank(), legend.text=element_text(size=24)) +
    theme(axis.text.x = element_text(size=24, angle = 90, hjust = 1), axis.text.y = element_text(size=24)) +
    geom_segment(data=temp,
                 aes(x=match(group,levels(group))-segmentWidth,
                     xend=match(group,levels(group))+segmentWidth,
                     y=score-se,yend=score-se),
                 col='black') +
    geom_segment(data=temp,
                 aes(x=match(group,levels(group))-segmentWidth,
                     xend=match(group,levels(group))+segmentWidth,
                     y=score+se,yend=score+se),
                 col='black') +
    geom_segment(data=temp,
                 aes(x=match(group,levels(group)),
                     xend=match(group,levels(group)),
                     y=score+se,yend=score-se),
                 col='black') +
    geom_segment(data=temp_line,
                 aes(x=x,xend=xend, y=y, yend=yend),
                 col='black') +
    geom_text(data=temp_text, aes(x=x, y=y, label = p_value), color="black", size=4) + 
    geom_point(data=temp, aes(x=group, y=score), color="black", size=1) 
  return (a)
}

createLines <- function(temp){
  #create lines for adding significance value
  a= data.frame(x=match(temp$group, levels(temp$group))[-length(temp$group)],
                xend= match(temp$group, levels(temp$group))[-1],
                y=temp$max_all[1] + temp$max_all[1]/10 +  match(temp$group, levels(temp$group))[-length(temp$group)]/10,
                yend=temp$max_all[1] +  temp$max_all[1]/10 + match(temp$group, levels(temp$group))[-length(temp$group)]/10)
  
  b = data.frame(x=match(temp$group, levels(temp$group))[-length(temp$group)],
                 xend= match(temp$group, levels(temp$group))[-length(temp$group)],
                 y=temp$max_all[1] + temp$max_all[1]/10 +  match(temp$group, levels(temp$group))[-length(temp$group)]/10,
                 yend=temp$max_all[1] +  temp$max_all[1]/10 + match(temp$group, levels(temp$group))[-length(temp$group)]/10 - 0.1)
  
  c = data.frame(x=match(temp$group, levels(temp$group))[-1],
                 xend= match(temp$group, levels(temp$group))[-1],
                 y=temp$max_all[1] + temp$max_all[1]/10 +  match(temp$group, levels(temp$group))[-length(temp$group)]/10,
                 yend=temp$max_all[1] +  temp$max_all[1]/10 + match(temp$group, levels(temp$group))[-length(temp$group)]/10 - 0.1)
  
  
  return(rbind(a, b, c))
}

createText <- function(currPheno, temp, alternative){
  p_value = sapply(1:(length(levels(currPheno$group)) -1), function(i){
    p= t.test(currPheno$score[currPheno$group == levels(currPheno$group)[i] ], currPheno$score[currPheno$group == levels(currPheno$group)[i +1] ], alternative)$p.value
    if (p > 0.05){
      "P>0.05"
    }else{
      paste("P=",format(p, scientific = T, digit=2), sep="")
    }
  })
  
  a= data.frame(x=match(temp$group, levels(temp$group))[-length(temp$group)] + 0.5,
                y=temp$max_all[1] + 1.5*temp$max_all[1]/10 +  match(temp$group, levels(temp$group))[-length(temp$group)]/10 )
  
  return(cbind(a, p_value))
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     max_group   = max     (xx[[col]], na.rm=na.rm)                     
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  datac$max_all <- max(datac$max_group)
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

segmentWidth = 0.025
jitterWidth = 0.1
jitterHight = 0
color = ""



a = rbind(data.frame(score = c(1,2,3,4,5), group="A"),
          data.frame(score = c(4,5,6,7,8), group="B"))
createGGPlot(a, measureVar="score", groupVars = "group", main = "", ylab = "", xlab = "", "two.sided")
