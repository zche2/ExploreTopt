#### load ####
library("nlme")
library(ggpubr)
library(latex2exp)
library(ggplot2) 
library(ggthemes)

#### test mixing effect ####
#### load data ####
df <- read.csv("/Users/hzc/data/forAnova.csv",
               header = TRUE, 
               colClasses = c("factor", "numeric", "numeric", "numeric","numeric","factor","numeric","numeric"))

#### anova ####
rst <- aov(Topt_mean ~ Tgs_mean*vegtype, data = df)
summary(rst)

#### full model comparison 1）random effect structure 2) fixed effect structure ####
B1 <- gls(Topt_mean ~ 1 + Tgs_mean, data = df,
          method = "REML")
B2 <- lme(Topt_mean ~ 1 + Tgs_mean * vegtype, data = df,
          method = "REML", random = ~1|vegtype)
B3 <- lme(Topt_mean ~ 1 + Tgs_mean * vegtype, data = df,
          method = "REML", random = ~1 + Tgs_mean|vegtype)
AIC(B1,B2,B3)  # suggests B1

## try only Tgs_mean for fixed effect and random intercp (and slope)
B4 <- lme(Topt_mean ~ 1 + Tgs_mean, data = df,
          method = "REML", random = ~1|vegtype)
ctrl <- lmeControl(opt='optim');
B5 <- lme(Topt_mean ~ 1 + Tgs_mean, data = df,
          method = "REML", random = ~1 + Tgs_mean|vegtype,
          control=ctrl)
AIC(B1,B4,B5)


#### plot scatter ####
lgd <- c('EBF', 'ENF', 'DBF', 'MF', 'SHR', 'SAV', 'GRA', 'CRO', 'WET')
color <- c("#228B22","#32CD32","#00FF7F","#FFA500","#F0E68C","#FFA07A","#DB7093","#00008B","#A52A2A")
df$predict_byGroup <- fitted(B4, level = 1)
df$predict_All <- fitted(B4, level = 0)
df$vegtype <- factor(df$vegtype,levels = lgd)  # 调整排序的优先级
figpop1 <- ggplot()+
  geom_point(df, mapping=aes(x=Tgs_mean,y=Topt_mean,
                             group=vegtype,color=vegtype),shape=1)+ 
  geom_line(df, mapping=aes(x=Tgs_mean, y=predict_byGroup,
                            group=vegtype, color=vegtype),linewidth=.6)+
  labs(x=TeX("Mean annual $T_{gs\\ max}^{air}$ (˚C)"),
       y=TeX("$T_{opt}^{eco}$ (˚C)"))+
  scale_color_manual(values = color)+ 
  theme_base(base_size = 14)+
  theme(plot.background = element_rect(color = NA),
        legend.title = element_blank(),
        legend.position = c(.2, .8),
        axis.text = element_text(size=13),
        axis.title = element_text(size=18)
        )+
  guides(col=
    guide_legend(ncol=2)
  )+
  coord_fixed()


figpop1  #在窗口右下角预览一下图
