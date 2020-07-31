#plots en español
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



#fig 3 allee guiad

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
