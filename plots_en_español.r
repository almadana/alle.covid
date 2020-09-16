#plots en español - para GUIAD 
require(ggpubr)
# -------- FIGURA EFECTO ALLEE (1)--------------
###
allee1D_esp <- ggplot(rdf, aes(x = Infected, y = R, color = model)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#b33018", "#14b74b"), name = "Modelo SIR",
                     labels = c("con efecto Allee", "sin efecto Allee")) +
  geom_hline(yintercept = 1, size = 1, linetype = "dashed") +
  geom_segment(aes(x = point1, xend = point1, y = 0, yend = 2), color = "red",
               linetype = "dashed") +
  xlab("Proporción de infectados") +
  ylab(bquote('Número reproductivo ' ~R[e] )) +
  geom_point(data = equilibriumDf, size = 4, color = "black") +
  geom_segment(data = dplyr::filter(arrowDf, direction == "right"),
               aes(x = Ii, xend = If, y = R, yend = R), color = "#fa3d1b",
               arrow = arrow(length = unit(0.3, "cm")), size = 1.1) +
  geom_segment(data = dplyr::filter(arrowDf, direction == "left"),
               aes(x = Ii, xend = If, y = R, yend = R), color =  "#243faf",
               arrow = arrow(length = unit(0.3, "cm"), ends = "last"),
               size = 1.1) +
  #geom_text(aes(x = point1+50, y = 0.7), label = thresholdText, color = "black",
  #hjust = "inward", family = "sans", fontface = "plain", ) +
  annotate("text",x=point1+50,y=0.7,label="Equilibrio  inestable \n Umbral de epidemia",size=3,hjust="inward") +
  annotate("text",x=point2,y=0.7,label="Equilibrio estable",size=3,hjust="inward") +
  annotate("text",x=point1-100,y=2.2,label="Efecto Allee",size=3,hjust="inward") +
  #  geom_text(aes(x = point2, y = 0.7), label = immunityText, color = "black",
  #           hjust = "inward", family = "sans", fontface = "plain" ) +
  theme_classic() +scale_x_continuous(labels = c(0,.25,.5,.75,1)) + 
  theme(legend.position = c(.7,.7),text=element_text(size=12)) 


allee1D_esp
ggsave("allee1D_esp.pdf", allee1D_esp, width = 6, height = 4)
ggsave("allee1D_esp.png", allee1D_esp, width = 6, height = 4)



#fig 3 allee guiad -----  simulaciones con y sin allee

nPaneles=32

fig3AlleeGuiad.a <- dplyr::filter(simulations_fit, allee == TRUE & rep %in% sampleRepetitionsSM[1:nPaneles]) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#b33018",size=.6) +
  #  facet_wrap(~rep, scales = "free", ncol = 4) +
  facet_wrap(~rep, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  scale_y_continuous(limits = c(10, NA), trans = "log10") +
  #xlab("time (days)",size=10) +
  #ylab("Cumulative infected") +
  xlab(element_text("tiempo (días)",size=12)) +
  ylab(element_text("Infectados totales",size=12)) +
  ggtitle("con efecto Allee") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5,size=12),
        text = element_text(size=6)) +
  geom_line(aes(x=t,y=cumI_fit),color="pink")

fig3AlleeGuiad.a



fig3AlleeGuiad.b <- dplyr::filter(simulations_fit, allee == FALSE & rep %in% sampleRepetitionsSM[1:nPaneles]) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#14b74b",size=.6) +
  #  facet_wrap(~rep, scales = "free", ncol = 4) +
  facet_wrap(~rep, scales = "free", ncol = 4) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  scale_y_continuous(limits = c(10, NA), trans = "log10") +
  #xlab("time (days)",size=10) +
  #ylab("Cumulative infected") +
  xlab(element_text("tiempo (días)",size=12)) +
  ylab(element_text("Infectados totales",size=12)) +
  ggtitle("sin efecto Allee") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5,size=12),
        text = element_text(size=6)) +
  geom_line(aes(x=t,y=cumI_fit),color="green")

fig3AlleeGuiad.b

fig3AlleeGuiad = ggarrange(fig3AlleeGuiad.a,fig3AlleeGuiad.b,ncol = 2)

ggsave(fig3AlleeGuiad,file="fig3AlleeGuiad.pdf",width=6,height=6)
ggsave(fig3AlleeGuiad,file="fig3AlleeGuiad.png",width=6,height=6)



#------------ FIG 2 - ESPACIO DE PARÁMETROS -------------

plotea_espacio_esp <- function(inputDF, xName, reverse = FALSE){
  spacePlot <- ggplot(inputDF, aes(x = measure, y = Infected, fill = Rt_exp)) +
  geom_raster() +
  scale_fill_manual(values = c("#243faf", "#fa3d1b"),
                    name = element_blank(), labels = c("Re < 1 (Control)",
                                                       "Re > 1 (Epidemia)")) +
  scale_x_continuous(name = xName, expand = c(0,0)) +
  scale_y_continuous(name = "Proporción de infectados", expand = c(0,0),limits = c(0,.8)) +
  theme_bw()
}
  
flechas = data.frame(
  x.det = c(100,100,400,600,600,800),
  y.det=c(.1,.5,.3,.1,.7,.3),
  direction=rep(c("up","down"),each=3)
)

fntSize=4

detectedPlot_esp <- dplyr::rename(detectedLong, measure = maxDetected) %>% mutate(Infected =Infected/ Npop) %>% 
  plotea_espacio_esp(., "Capacidad de detección") +
  ggtitle("Sistema de trazado de casos") +
  theme(plot.title = element_text(hjust = 0.5)) +
#  geom_segment(data=flechas,aes(x=x.det,xend=x.det,y=y.det,yend=y.det+.1),arrow=arrow(length = unit(2,"mm"))) +
  annotate("text",x=200,y=.35,label="epidemia",size=fntSize,col="white")+
  annotate("text",x=850,y=.35,label="control",size=fntSize,col="white")
  
callsPlot_esp <- dplyr::rename(callsLong, measure = maxCalls) %>% mutate(Infected =Infected/ Npop) %>% 
  plotea_espacio_esp(., "Número máximo de llamadas") +
  ggtitle("Velocidad de trazado de contactos") +
  theme(plot.title = element_text(hjust = 0.5)) 


# Space plot for maximum number of links
linksPlot_es <- dplyr::rename(linksLong, measure = maxLinks) %>% mutate(Infected =Infected/ Npop) %>% 
  plotea_espacio_esp(., "Número máximo de contactos") +
  ggtitle("Reducción de contactos") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text",x=25,y=.35,label="epidemia",size=fntSize,col="white")+
  annotate("text",x=7,y=.35,label="control",size=fntSize,col="white")


# Space plot for pInfection 
infectionPlot_es <- dplyr::rename(infectionLong, measure = pInfection) %>%mutate(Infected =Infected/ Npop) %>% 
  plotea_espacio_esp(., "Probabilidad de transmisión en un contacto") +
  ggtitle("Distanciamiento físico e higiene") +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text",x=.26,y=.35,label="epidemia",size=fntSize,col="white")+
  annotate("text",x=.13,y=.35,label="control",size=fntSize,col="white")

# dinámica de espacio-fase, relajación de NPIs

x1 = 300
x2 = 580
y1 = .16
y2 = .3
y1bis = .05
xOffset = 25
yOffset= 15
colSquare = "yellow"
colLine = "white"
linSquare="dashed"
linLine = "dotted"
fontSize = 3
xLim = c(0,800)
yLim = c(0,.5)


npi.upper.plot = detectedPlot_esp + scale_x_continuous(breaks = c(0, 400, 800), 
                                                   labels = c("","",""),
                                                   expand = c(0,0),
                                                   limits = xLim) +
  scale_y_continuous(expand=c(0,0),limits=yLim)+
  ggtitle("") +
  labs(x="Fuerza de las InF",y="Proporción de infectados",title = "Relajamiento de las medidas")+
  
  theme(#axis.title.y=element_blank(),
    legend.key.size = unit(.5,"cm"),
    legend.text = element_text(size = 9),
    legend.position = "right")+
  annotate(geom="point",x=x1,y=y1,col=colSquare)+
  annotate(geom="point",x=x2,y=y1,col=colSquare)+
  annotate(geom="point",x=x1,y=y2,col=colSquare)+
  annotate(geom="point",x=x2,y=y2,col=colSquare) + 
  annotate(geom = "text",x=x2+xOffset,y1,label="italic(i)",size=fontSize,parse=T,hjust="outward",col=colSquare)+
  annotate(geom = "text",x=x2+xOffset,y2,label="italic(iv)",size=fontSize,parse=T,hjust="outward",col=colSquare)+
  annotate(geom = "text",x=x1-xOffset,y1,label="italic(ii)",size=fontSize,parse=T,hjust="outward",col=colSquare)+
  annotate(geom = "text",x=x1-xOffset,y2,label="italic(iii)",size=fontSize,parse=T,hjust="outward",col=colSquare)+
  annotate(geom="segment",x=x2,xend = x1+25,y=y1,yend=y1,col=colSquare,
           arrow = arrow(length = unit(0.2, "cm"), ends = "last",type="closed"),linetype=linSquare) +
  annotate(geom="segment",x=x1,xend = x2-25,y=y2,yend=y2,col=colSquare,
           arrow = arrow(length = unit(0.2, "cm"), ends = "last",type = "closed"),linetype=linSquare)+
  annotate(geom="segment",x=x1,xend = x1,y=y1,yend=y2-.02,col=colSquare,
           arrow = arrow(length = unit(0.2, "cm"), ends = "last",type = "closed"),linetype=linSquare)+
  annotate(geom="point",x=x1,y=y1bis,col=colLine)+
  annotate(geom = "text",x=x1-xOffset,y1bis,label="italic(v)",size=fontSize,parse=T,hjust="outward",col=colLine)+
  annotate(geom="segment",x=x2,xend = x1+xOffset,y=y1,yend=y1bis,col=colLine,
           arrow = arrow(length = unit(0.2, "cm"), ends = "last",type="closed"),linetype=linLine) +
  annotate("text",x=650,y=.05,label="control",size=fntSize,col="white")


npi.upper.plot



# # arrange plots into grid
# legendPlot <- detectedPlot + theme(legend.position = "top")
# legend <- gtable_filter(ggplotGrob(legendPlot), "guide-box")
# 
# plotGrid <- grid.arrange(legend,
#                          arrangeGrob(detectedPlot + theme(legend.position="none",
#                                                           axis.title.y=element_blank()),
#                                      callsPlot + theme(legend.position="none",
#                                                        axis.title.y=element_blank()),
#                                      infectionPlot + theme(legend.position="none",
#                                                            axis.title.y=element_blank()),
#                                      linksPlot + theme(legend.position="none",
#                                                        axis.title.y=element_blank()),
#                                      ncol=2, left = "Proportion infected"),
#                          heights=c(1, 10))
# 

fig2_guiad = ggarrange(detectedPlot_esp,linksPlot_es,infectionPlot_es,npi.upper.plot
,nrow = 4, common.legend = T,labels = "AUTO")

ggsave("phase_space_guiad.pdf", fig2_guiad, width = 4, height = 10)
ggsave("phase_space_guiad.png", fig2_guiad, width = 4, height = 10)


#-------- fig condados y países -------------

nExamples=10 # cuantos de cada
county_examples = 
  c(
    county_data_fit %>% filter(cociente>1)%>% group_by(Id) %>% slice(1) %>% ungroup() %>% 
      slice_max(cociente,n=nExamples) %>% pull(Id)
    ,
    county_data_fit %>% filter(cociente>1) %>% group_by(Id) %>% slice(1) %>% ungroup() %>% 
      slice_min(cociente,n=nExamples) %>% pull(Id)
  )


plotsCondados_esp <- dplyr::filter(county_data_fit, Id %in% county_examples) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#000080") +
  facet_wrap(~county.full, scales = "free", ncol = 5) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  #scale_y_continuous(limits = c(10, NA), trans = "log10") +
  scale_y_continuous(trans = "log10") +
  xlab("tiempo (días)") +
  ylab("Infectados totales") +
  ylab(element_blank()) +
  ggtitle("Condados de EE.UU.") +
  theme_classic() +
  theme(strip.background = element_rect(linetype = 0),
        #strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=9),
        strip.text.x = element_text(size=7)) +
  geom_line(aes(x=t,y=cumI_fit),color="sky blue")

plotsCondados_esp





country_examples = 
  c(
    country_data_fit %>% filter(cociente>1) %>% group_by(Id) %>% slice(1) %>% ungroup() %>% 
      slice_max(cociente,n=nExamples) %>% pull(Id)
    ,
    country_data_fit %>% filter(cociente>1) %>% group_by(Id) %>% slice(1) %>% ungroup() %>% 
      slice_min(cociente,n=nExamples) %>% pull(Id)
  )


plotsPaises_esp <- dplyr::filter(country_data_fit, Id %in% country_examples) %>%
  #plotsPaises <- dplyr::filter(country_data_fit, Id %in% 108) %>%
  mutate(country=ifelse(country %in% "Korea, South","South Korea",country)) %>% 
  #plotsPaises <- dplyr::filter(country_data_fit, Id %in% c(166,39)) %>%
  ggplot(., aes(x = t, y = cumI)) +
  geom_point(color = "#4020ab") +
  facet_wrap(~country, scales = "free", ncol = 5) +
  scale_x_continuous(breaks = c(1, 5, 25), trans = "log10") +
  #scale_y_continuous(limits = c(10, NA), trans = "log10") +
  scale_y_continuous(trans = "log10") +
  xlab("tiempo (días)") +
  ylab("Cumulative infected") +
  ylab(element_blank()) +
  ggtitle("Países y regiones") +
  theme_classic() +
  theme(strip.background = element_rect(linetype = 0),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=9),
        strip.text.x = element_text(size=7)) +
  geom_line(aes(x=t,y=cumI_fit),color="sky blue")

plotsPaises_esp

plot_condados_paises_esp = ggarrange(plotsCondados_esp,plotsPaises_esp,nrow = 2)
ggsave(plot_condados_paises_esp,file="fig_2_trayectorias_esp.pdf",width=6,height=8)
ggsave(plot_condados_paises_esp,file="fig_2_trayectorias_esp.png",width=6,height=8)
