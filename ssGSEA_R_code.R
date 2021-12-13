# ssGSEA by cluster
clustering_info <-readxl::read_xlsx("/Users/davydeng/Desktop/EAS/CGGA.xlsx")
clustering_info <- clustering_info[!is.na(clustering_info$Subtype),]
clustering_info <- clustering_info %>% select(Subtype, ...1)
clustering_info <- clustering_info %>% rename(sample = ...1)
ssGSEA <- fread("/Users/davydeng/Desktop/EAS/ssgseaOut.txt")
ssGSEA_t <- as.data.frame(t(ssGSEA))
colnames(ssGSEA_t) <- ssGSEA_t[1,]
ssGSEA_t <- ssGSEA_t[-1,]
ssGSEA_t$sample <- rownames(ssGSEA_t)
ssGSEA_t <- right_join(clustering_info, ssGSEA_t, by='sample')
ssGSEA_t_dat <- melt(ssGSEA_t,measure.vars=colnames(ssGSEA_t)[3:72])
ssGSEA_t_dat$value <- as.numeric(ssGSEA_t_dat$value)


sapply(as.character(unique(ssGSEA_t_dat$variable)), function(x){
  ggplot(ssGSEA_t_dat[ssGSEA_t_dat$variable==x,], 
         aes(x= Subtype, y= value)) + geom_boxplot(position=position_dodge(1)) +
    theme_classic() + labs(x="Subtype", y = "value")+
    stat_compare_means(aes(label = format.pval(..p.adj..)),
                       label.y = c(seq(1.2, by = 0, length.out = 4)),
                       tip.length = 0, vjust = 0.5,
                       method = 'wilcox.test',ref.group='.all.')
  ggsave(paste0('EAS_cluster_',x,'.pdf'),width = 16, height = 20)
}
)

ssGSEA_df_wilcox <- as.data.frame(t(data.frame(sapply(as.character(unique(ssGSEA_t_dat$variable)), function(x){
  sapply(as.character(unique(ssGSEA_t_dat$Subtype)),function(y){
    wilcox.test(ssGSEA_t_dat[ssGSEA_t_dat$variable==x & ssGSEA_t_dat$Subtype == y,]$value,
                ssGSEA_t_dat[ssGSEA_t_dat$variable==x & ssGSEA_t_dat$Subtype != y,]$value)$p.value
  }
  )
}
))))


ssGSEA_df_wilcox$PL <- qvalue(ssGSEA_df_wilcox$PL, lambda=0)$qvalues
ssGSEA_df_wilcox$NF <- qvalue(ssGSEA_df_wilcox$NF, lambda=0)$qvalues
ssGSEA_df_wilcox$IM <- qvalue(ssGSEA_df_wilcox$IM, lambda=0)$qvalues
ssGSEA_df_wilcox$MB <- qvalue(ssGSEA_df_wilcox$MB, lambda=0)$qvalues

ssGSEA_df_wilcox <- -log10(ssGSEA_df_wilcox)
colnames(ssGSEA_df_wilcox)[2] <- 'NS'

write.csv(ssGSEA_df_wilcox, 'ssGSEA_df_wilcox.csv',row.names = TRUE)

