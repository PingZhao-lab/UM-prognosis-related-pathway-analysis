data <- data.frame(V=c(0.8125,0.6875),M=c("UM prognostic pathways","Differential pathways"),
                   C=c(rep("Signatures",2)))
data$M<-factor(data$M,levels = c("Differential pathways","UM prognostic pathways"))

ggplot(data,aes(x=M, y=V, fill=C))+
  geom_bar(alpha=1,stat='identity',position="dodge",colour="black",size=1)+
  scale_fill_manual("legend", values = c("Signatures" = "#E27069"))+
  labs(x = NULL, y = "Maximum accuracy")+
  theme(legend.position="none",text = element_text(color ="black",size=35),
        axis.text.x = element_text(color="black",size=30),
        axis.text.y = element_text(color="black",size=30),
        axis.ticks.x=element_line(color="black",size=2,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=1.5),
        axis.ticks.y=element_line(color="black",size=2,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=1.5),
        plot.title = element_text(size = 45),panel.background = element_blank(),
        strip.text.x = element_text(size = 45)
  )
