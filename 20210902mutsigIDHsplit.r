library(data.table)
metadD <- fread('/Users/fdubois/Dropbox (Partners HealthCare)/EAS_GBM/EAS_GBM_NonShare/data/20210701ancestry/20210701metaDexomes.txt')
mafCGC <- fread('/Users/fdubois/Dropbox (Partners HealthCare)/EAS_GBM/EAS_GBM_Nonshare/data/20210811mutsig/20210811_EASaGBM.maf')

metadD_IDHwt <- metadD[IDHstat == 'IDHwt']
mafCGC_IDHwt <- mafCGC[Tumor_Sample_Barcode %in% metadD_IDHwt$sampleID]
write.table(mafCGC_IDHwt,  '/Users/fdubois/Dropbox (Partners HealthCare)/EAS_GBM/EAS_GBM_Nonshare/data/20210902mutsig/20210902_EASaGBM_IDHwt.maf', quote = FALSE, sep = "\t", row.names = F)


mafCGC_IDHmut <- mafCGC[!(Tumor_Sample_Barcode %in% metadD_IDHwt$sampleID)]
write.table(mafCGC_IDHmut,  '/Users/fdubois/Dropbox (Partners HealthCare)/EAS_GBM/EAS_GBM_Nonshare/data/20210902mutsig/20210902_EASaGBM_IDHmut.maf', quote = FALSE, sep = "\t", row.names = F)
