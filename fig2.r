fig2 = ggarrange(plotsCondados + labs(tag="a"),
          plot_slopes + labs(tag="b"), 
          plotsPaises+labs(tag="c"),
          plot_slopes_p + labs(tag="d"))

fig2

ggsave(fig2,file="fig2.pdf",width=14,height=9)
ggsave(fig2,file="fig2.png",width=14,height=9)
