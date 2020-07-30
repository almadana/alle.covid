#fig weighted evidence in counties and countries... dang!

View(c(bb$Weight.evid,
bb_p$Weight.evid))

weight_evid_df = data.frame(weight.evid=c(bb$Weight.evid,bb_p$Weight.evid),
                            geo= c( rep("U.S. counties",nrow(bb)), rep("countries & regions" ,nrow(bb_p)) )
                            )

relevel(weight_evid_df$geo,ref="U.S. counties")

weight_evid_df$weight.evid[is.na(weight_evid_df$weight.evid)]=0

# fig_weight_evid_ccr = 
#   weight_evid_df %>% 
#   ggplot(aes(x=weight.evid,fill=geo))+
#   geom_histogram(aes(y=..density..)) +
#   facet_grid(~geo)+
#   labs(x="Weighted evidence for breaking point",y="# of cases",fill="")+
#   scale_fill_manual(values = c("#4020ab", "#000080"))+
#   theme_classic() +
#   theme(strip.background = element_blank(),
#         plot.title = element_text(hjust = 0.5,size=10),
#         text = element_text(size=12),
#         legend.position = "none")  
# 
# 
# fig_weight_evid_ccr

fig_weight_evid_counties = 
  bb %>% 
  ggplot(aes(x=Weight.evid))+
  geom_histogram(fill="#000080") +
  labs(x="Weighted evidence for breaking point",y="# of cases",title="U.S. counties")+
  #scale_fill_manual(values = c("#4020ab", "#000080"))+
  theme_classic() +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size=10),
        text = element_text(size=12),
        legend.position = "none")  


fig_weight_evid_countries = 
  bb_p %>% 
  ggplot(aes(x=Weight.evid))+
  geom_histogram(fill="#4020ab") +
  labs(x="Weighted evidence for breaking point",y="# of cases",title="countries & regions")+
  #scale_fill_manual(values = c("#4020ab", "#000080"))+
  theme_classic() +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size=10),
        text = element_text(size=12),
        legend.position = "none")  

fig_weight_evid_ccr=ggarrange(fig_weight_evid_counties,fig_weight_evid_countries,labels = "auto",label.x = c(1,0),nrow=1)

ggsave(fig_weight_evid_ccr,file="figs_weit_evid_ccr.pdf",width=8,height=4)
ggsave(fig_weight_evid_ccr,file="figs_weit_evid_ccr.png",width=8,height=4)


