#Author: Bin Chen (2020)
#scatter plot for correlation analysis between two variables

library(ggplot2)

d1 = c(1, 2, 3, 4, 5)
d2 = c(1, 3, 4, 6, 7)
d3 = as.factor(c(0, 1, 0, 1, 1))
d = data.frame(d1, d2, d3)

t = cor.test(d$d1, d$d2)
l = lm(d1 ~ d2, d)

ggplot(d, aes(x= (d1  ), y = (d2  ), colour = d3 )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  annotate("text", label = paste("r=", format(summary(l)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(l)$`Pr(>F)`[1], digit=2), sep=""), x = 5, y = 10, size = 6, colour = "black") +
  annotate("text", label = paste("rho=", format(t$estimate, digit=2), ", P=", format(t$p.value, digit=3, scientific=T), sep=""), x = 5, y = 9.5, size = 6, colour = "black") +
  xlab("d1") + guides(shape=FALSE, size=FALSE) +
  ylab("d2") + coord_cartesian(xlim = c(0, 10), ylim=c(0, 10))
