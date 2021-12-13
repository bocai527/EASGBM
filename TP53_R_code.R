### P53_eur
library(maftools)
library(stringr)
Transactivation_domain <- c(1:50)
SH3_domain <- c(63:97)
DNA_binding <- c(102:292)
Tetramerization_reg <- c(323:393)

TCGA_eur <- read.maf("/Users/davydeng/Desktop/TCGA_eur_maftools.maf")
CGGA_gbm <- read.maf("/Users/davydeng/Desktop/CGGA_gbm_maftools_new.maf")


TP53 <- TCGA_eur@data[TCGA_eur@data$Hugo_Symbol == 'TP53',]
TP53_silent <- TCGA_eur@maf.silent[TCGA_eur@maf.silent$Hugo_Symbol == 'TP53',]

tp53_protein <- data.frame(V1=TP53$Protein_position, V2= TP53$Start_Position)
tp53_protein_silent <- data.frame(V1=TP53_silent$Protein_position, V2= TP53_silent$Start_Position)
tp53_protein <- rbind(tp53_protein,tp53_protein_silent)

tp53_protein$V1 <- as.numeric(str_extract(tp53_protein$V1, '\\d+'))
names(tp53_protein)[2] <- 'POS'
names(tp53_protein)[1] <- 'Protein_POS'


sum(tp53_protein %in% Transactivation_domain)
sum(tp53_protein %in% SH3_domain)
sum(tp53_protein %in% DNA_binding)
sum(tp53_protein %in% Tetramerization_reg)

### p53_eas

tp53_EAS_protein <- data.frame(V1=TP53_EAS$aaChange, V2= TP53_EAS$Start_Position)
tp53_EAS_protein_silent <- data.frame(V1=TP53_EAS_silent$aaChange, V2= TP53_EAS_silent$Start_Position)
tp53_EAS_protein <- rbind(tp53_EAS_protein,tp53_EAS_protein_silent)

tp53_EAS_protein$V1 <- as.numeric(str_extract(tp53_EAS_protein$V1, '\\d+'))
names(tp53_EAS_protein)[2] <- 'POS'
names(tp53_EAS_protein)[1] <- 'Protein_POS'



sum(tp53_EAS_protein %in% Transactivation_domain)
sum(tp53_EAS_protein %in% SH3_domain)
sum(tp53_EAS_protein %in% DNA_binding)
sum(tp53_EAS_protein %in% Tetramerization_reg)

p53_trans <- matrix(c(2, 23, 117, 54), nrow=2, dimnames=list(race=c('EUR','EAS'),
                                                             domain= c('trans','non-trans')))

fisher.test(p53_trans)

p53_SH3 <- matrix(c(1, 20, 118, 57), nrow=2, dimnames=list(race=c('EUR','EAS'),
                                                           domain= c('SH3','non-SH3')))
fisher.test(p53_SH3)

p53_DNA <- matrix(c(111, 34, 8, 43), nrow=2, dimnames=list(race=c('EUR','EAS'),
                                                           domain= c('DNA','non-DNA')))
fisher.test(p53_DNA)

p53_tetra <- matrix(c(5, 0, 114, 77), nrow=2, dimnames=list(race=c('EUR','EAS'),
                                                            domain= c('Tetra','non-Tetra')))
fisher.test(p53_tetra)