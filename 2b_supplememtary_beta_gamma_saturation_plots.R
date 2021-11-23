# supplementary materials: plot of beta and gamma for different I levels

require(viridis)
require(tidyverse)
require(ggpubr)

minI <- 10
maxI <- 2000





# epidemic values
I50 <- c(20,100,200,500)
betaMax <- 0.5
gammaMax <- 4


I <- seq(minI,maxI,by = 10)



df <- data.frame(I=I,rep(I50,each=length(I)/length(I50)))

df = tidyr::expand(df,I,I50)

df$beta <- betaMax * df$I/ (df$I + df$I50)
df$gamma <- gammaMax * df$I / (df$I + df$I50)
#df$I = df$I/maxI

df$label=""
labels=df %>% filter(I50==I) %>% mutate(label=paste0("I[50]==",I50))

df[df$I50==df$I,]=labels
df$beta_label = df$beta
df$gamma_label = df$gamma
df$i_label = df$I
df[df$I50==df$I,"beta_label"]<-c(.46,.42,.413,.385)
df[df$I50==df$I,"i_label"]<-c(80,350,750,1400)

#---- beta plot------

beta.plot = ggplot(df,aes(x=I,y=beta,col=factor(I50),group=factor(I50)))+geom_line(size=1.2) + 
  geom_hline(yintercept = betaMax,linetype="dotted",size=1.2,col="#777777")+
  annotate(x=0,y=betaMax+.015,label="beta[max]",parse=T,geom="text",size=6)+
  geom_text(aes(x=i_label,label=label,y=beta_label),parse=T,size=4.5)+
  theme_classic()+
  scale_color_viridis(discrete=T,end=.7)+
  theme(text=element_text(size=14),legend.position = "none")+
  labs(x="I",y=expression(beta(I)),parse=T)
beta.plot


# ----- GAMMA PLOT-------
df[df$I50==df$I,"gamma_label"]<-c(3.65,3.4,3.3,3.08)
gamma.plot = ggplot(df,aes(x=I,y=gamma,col=factor(I50),group=factor(I50)))+geom_line(size=1.2) + 
  geom_hline(yintercept = gammaMax,linetype="dotted",size=1.2,col="#777777")+
  annotate(x=0,y=gammaMax+.15,label="gamma[max]",parse=T,geom="text",size=6)+
  geom_text(aes(x=i_label,label=label,y=gamma_label),parse=T,size=4.5)+
  theme_classic()+
  scale_color_viridis(discrete=T,end=.7)+
  theme(text=element_text(size=14),legend.position = "none")+
  labs(x="I",y=expression(gamma(I)),parse=T)
gamma.plot

fig.supp.beta.gamma = ggarrange(beta.plot,gamma.plot,labels = "AUTO")
ggsave(fig.supp.beta.gamma,width=14,height=6,file="plots/Sbetagamma.pdf")
ggsave(fig.supp.beta.gamma,width=14,height=6,file="plots/Sbetagamma.png")
