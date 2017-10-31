source('sampleqc.R')

metadata = read_metadata('all')

empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

metadata$low_contam = metadata$freemix < 0.075
scatter <- ggplot(metadata, aes(x = xhet, y = ycov, color = sex, alpha = low_contam)) + 
  geom_point(size=1.0) + 
  theme_bw() +
  theme(legend.position = c(0.85, .63),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))+
  scale_y_continuous("Normalized Y coverage") +
  scale_x_continuous("X heterozygosity")

#marginal density of x - plot on top
plot_top <- ggplot(metadata, aes(x = xhet)) + 
  geom_density(alpha=.5,fill="blue") + 
  theme_bw() +
  theme(axis.title.x = element_blank())

#marginal density of y - plot on the right
plot_right <- ggplot(metadata, aes(ycov)) + geom_density(alpha=.5,fill="red") + 
  coord_flip() + 
  theme_bw()
  
  
library(gridExtra)
pdf("XhetYcov_new_combined_sex_2.pdf",width=8,height=8)
grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()
png("XhetYcov_new_combined_sex_2.png",width=1500,height=1200,res=300)
grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()


