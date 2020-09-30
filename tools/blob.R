library("ggplot2")
png("tgs.png",width = 900, height = 900)
data=read.csv("tgs.txt",sep=',',header=T)
data[,'name']=factor(data[,'name'],levels=c('chr19','E.coli','R12','43-1A','44A','LB8','M714','Y319','Y983') )
ggplot(data,aes(x=kmer.coverage,y=assembly.coverage))+
    geom_point(aes(size=kmer.depth,color=name)) + 
    scale_color_manual(values=c('red','darkorange','gold','green','cyan','turquoise','darkcyan','royalblue','navy') ) +
    scale_x_log10(breaks=c(0.1,1,10,100)) + 
    scale_y_log10(breaks=c(0.1,1,10,100)) + 
    labs(x="k-mer completeness in reads" , y="k-mer completeness in assembly") + 
    theme(text = element_text(size = 20)) +
    guides(size=guide_legend(title="k-mer depth",order=1) , color=guide_legend(title='species')) 

dev.off()
