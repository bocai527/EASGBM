### devtools::install_github("FunGeST/Palimpsest")
library(Palimpsest)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38) ### TCGA uses hg38
library(maftools)
library(reshape2)
library(ggplot2)
library(ggpubr)

laml = read.maf(maf = "/Users/davydeng/Desktop/evolution/CGGA_gbm_maftools_new.maf")
lam2 = read.maf(maf = "/Users/davydeng/Desktop/TCGA_eur_maftools.maf")

lam1_maf <- laml@data
lam1_maf <- lam1_maf %>% select(c(Tumor_Sample_Barcode, chrom, pos, ref, alt,
                                  n_F1R2, t_F1R2, t_alt_count, Hugo_Symbol,
                                  Variant_Type))
lam1_maf <- lam1_maf %>% dplyr::rename(Sample = Tumor_Sample_Barcode,
                                       Type = Variant_Type,
                                       CHROM = chrom,
                                       POS = pos,
                                       REF = ref, ALT= alt, Tumor_Varcount= t_alt_count,
                                       Tumor_Depth = t_F1R2, Normal_Depth = n_F1R2,
                                       Gene_Name = Hugo_Symbol)
levels(lam1_maf$Type) <- c(levels(lam1_maf$Type), "SNV") 
lam1_maf$Type[lam1_maf$Type == 'SNP']  <- "SNV" 
lam1_maf <- na.omit(lam1_maf)

vcf <- annotate_VCF(vcf = lam1_maf, ref_genome = BSgenome.Hsapiens.UCSC.hg19)
SBS_input <- palimpsest_input(vcf = vcf, Type = "SBS")

### De novo extraction
SBS_denovo_sigs <- NMF_Extraction(input_matrices = SBS_input,
                                  range_of_sigs = 2:10, nrun = 10,
                                  method = "lee")
compare_tab <- compare_results(reference_sigs = SBS_cosmic, 
                               extraction_1 = SBS_denovo_sigs)
resdir <- file.path("/Users/davydeng/Desktop/evolution")
pdf(file.path(resdir, "Cosine_Similarity_Heatmap.pdf"), width = 11, height = 10)
SBS_cosine_similarities <- deconvolution_compare(SBS_denovo_sigs, SBS_cosmic)
dev.off()
save(SBS_cosine_similarities, file = file.path(resdir, "Cosine_Similarity_matrix.RData"))



# Calculate and plot the exposure of the signatures across the series
SBS_col <- signature_colour_generator(rownames(SBS_denovo_sigs))
SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input, input_signatures = SBS_denovo_sigs,
                                        resdir = resdir, signature_colours = SBS_col, input_vcf = vcf)

pdf(file.path(resdir, "signature_content_plot.pdf"), width=15, height=10)
deconvolution_exposure(signature_colours = SBS_col, signature_contribution = SBS_signatures_exp)
dev.off()


### use the identified signatures from maftools
SBS_EAS_names <- c("SBS1","SBS5","SBS6",'SBS3', 'SBS10b','SBS11') 
SBS_EAS_sigs <- SBS_cosmic[rownames(SBS_cosmic) %in% SBS_EAS_names,]
sig_cols <- signature_colour_generator(rownames(SBS_EAS_sigs))
SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input,
                                        input_signatures = SBS_EAS_sigs,
                                        signature_colours = sig_cols,
                                        resdir = resdir)

pdf(file.path(resdir, "signature_content_plot_2.pdf"), width=15, height=10)
deconvolution_exposure(signature_contribution = SBS_signatures_exp,
                       signature_colours = sig_cols)
dev.off()



### The same procedure for EUR cohort
lam2_maf <- lam2@data
lam2_maf <- lam2_maf %>% select(c(Tumor_Sample_Barcode, Chromosome, Start_Position, 
                                  Reference_Allele, Tumor_Seq_Allele2,
                                  n_depth, t_depth, t_alt_count, Hugo_Symbol,
                                  Variant_Type))
lam2_maf <- lam2_maf %>% dplyr::rename(Sample = Tumor_Sample_Barcode,
                                       Type = Variant_Type,
                                       CHROM = Chromosome,
                                       POS = Start_Position,
                                       REF = Reference_Allele, 
                                       ALT= Tumor_Seq_Allele2, Tumor_Varcount= t_alt_count,
                                       Tumor_Depth = t_depth, Normal_Depth = n_depth,
                                       Gene_Name = Hugo_Symbol)
levels(lam2_maf$Type) <- c(levels(lam2_maf$Type), "SNV") 
lam2_maf$Type[lam2_maf$Type == 'SNP']  <- "SNV" 
lam2_maf <- na.omit(lam2_maf)
vcf_2 <- annotate_VCF(vcf = lam2_maf, ref_genome = BSgenome.Hsapiens.UCSC.hg38)
SBS_input_2 <- palimpsest_input(vcf = vcf_2, Type = "SBS")

### Skipped de novo
SBS_EAS_names <- c("SBS1","SBS5","SBS6",'SBS3', 'SBS10b','SBS11') ### SBS3, 10b, 11
SBS_EAS_sigs <- SBS_cosmic[rownames(SBS_cosmic) %in% SBS_EAS_names,]
sig_cols <- signature_colour_generator(rownames(SBS_EAS_sigs))
SBS_signatures_exp_2 <- deconvolution_fit(input_matrices = SBS_input_2,
                                          input_signatures = SBS_EAS_sigs,
                                          signature_colours = sig_cols,
                                          resdir = resdir)

pdf(file.path(resdir, "signature_content_plot_3.pdf"), width=15, height=10)
deconvolution_exposure(signature_contribution = SBS_signatures_exp_2,
                       signature_colours = sig_cols)
dev.off()


### check if cohorts have different distributions of signature probabilities
EAS_SBS <- as.data.frame(SBS_signatures_exp$sig_props)
EAS_SBS <- melt(EAS_SBS,  measure.vars=colnames(EAS_SBS))
EAS_SBS$supp <- 'EAS'
EUR_SBS <- as.data.frame(SBS_signatures_exp_2$sig_props)
EUR_SBS <- melt(EUR_SBS,  measure.vars=colnames(EUR_SBS))
EUR_SBS$supp <- 'EUR'
EAS_EUR_SBS <- rbind(EAS_SBS, EUR_SBS)

pdf(file.path(resdir, "Pamplimsest.pdf"), width=15, height=10)
ggplot(EAS_EUR_SBS, aes(x=variable, y=value, fill=supp)) + 
  geom_boxplot(position=position_dodge(1))+ theme_classic()+
  labs(x="SBS signature", y = "Probability")+
  stat_compare_means(aes(label = format.pval(..p.adj..)),
                     label.y = c(seq(1.2, by = 0, length.out = 4)),
                     tip.length = 0, vjust = 0.5,
                     method = "t.test")
dev.off()


### Check if tp53 gene has different signature distribution between cohorts
vcf_mut <- signature_origins(input = vcf, Type = "SBS",
                             input_signatures = SBS_EAS_sigs,
                             signature_contribution = SBS_signatures_exp)
vcf_mut_2 <- signature_origins(input = vcf_2, Type = "SBS",
                               input_signatures = SBS_EAS_sigs,
                               signature_contribution = SBS_signatures_exp_2)

vcf_mut_tp53_full <- vcf_mut %>% filter(Gene_Name=='TP53')
vcf_mut_2_tp53_full <- vcf_mut_2 %>% filter(Gene_Name=='TP53')

vcf_mut_tp53 <- vcf_mut_tp53_full[, c(18:23)]
vcf_mut_2_tp53 <- vcf_mut_2_tp53_full[,c(18:23)]
tp53_vcf <- melt(vcf_mut_tp53,  measure.vars=colnames(vcf_mut_tp53))
tp53_vcf_2 <- melt(vcf_mut_2_tp53,  measure.vars=colnames(vcf_mut_2_tp53))
tp53_vcf$supp <- 'EAS'
tp53_vcf_2$supp <- 'EUR'
EAS_EUR_tp53_SBS <- rbind(tp53_vcf, tp53_vcf_2)

pdf(file.path(resdir, "TP53_Pamplimsest.pdf"), width=15, height=10)
ggplot(EAS_EUR_tp53_SBS, aes(x=variable, y=value, fill=supp)) + 
  geom_boxplot(position=position_dodge(1))+ theme_classic()+
  labs(x="SBS signature", y = "Probability")+
  stat_compare_means(aes(label = format.pval(..p.adj..)),
                     label.y = c(seq(1.2, by = 0, length.out = 4)),
                     tip.length = 0, vjust = 0.5,
                     method = "t.test")
dev.off()
