fig2 = ggarrange(plotsCondados ,
          plot_slopes , 
          plotsPaises,
          plot_slopes_p ,
          widths = c(4,3),
          labels = "auto")

fig2

ggsave(fig2,file="fig2.pdf",width=14,height=9)
ggsave(fig2,file="fig2.png",width=14,height=9)
