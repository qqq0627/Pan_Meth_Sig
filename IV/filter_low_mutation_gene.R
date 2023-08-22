############## Before intrumental variable analysis,we select those genes with more than 10 mutations of missensen and frame shift mutations, respectively, then choose the significant 
############## correlation analysis between gene mutations and each methylation signature 
#########################################################################################################################################################
#only missense
library(stringr)
library(ivpack)
load("~/hyper_hypo_coef_H.RData")  ## methlation signatures score
gene_matrix <- read.table("~/reference/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv.gz", header = T, sep = "\t") #20531 11070 gene expression
gene <- gene_matrix
gene$gene_id <- str_split(gene$gene_id, '[|]', simplify = T)[,1]
gene <- gene[!duplicated(gene$gene_id),]
rownames(gene) <- gene$gene_id
gene <- gene[,-1]
gene <- as.data.frame(t(gene))
gene <- log2(gene + 1)
genelist <- read.csv("~/intersect_genelist.csv")  ##genes with both RNA expression and mutation 

load("~/TCGA_snv_missense.RData")    ##snv_mis   
snv_domin <- snv_mis[colnames(snv_mis) %in% genelist$x]
snv_domin$sample_ID <- str_sub(rownames(snv_domin), 1, 12)
cof_snv <- merge(hyper_hypo_coef, snv_domin)             ######10273 15611

snv_cor_ge <- mclapply(colnames(cof_snv)[13:length(cof_snv)], function(p){
	sig <- lapply(colnames(cof_snv)[2:11], function(q){
		cat(paste(p,q), "\n")
		df <- cof_snv[c(q, "sample_ID", "cancertype", p)]
        df[,4] <- as.factor(df[,4])
		if(dim(table(df[,4])) == 1|length(which(df[,4] == 1)) < 10){
        re_p <- data.frame(sig = q,
        	    gene = p, 
        	    mut_rate = NA,
                p.values= NA) 
		}else{
            tes <- wilcox.test(df[,1] ~ df[,4])
            re_p <- data.frame(sig = q,
            	gene = p, 
        	    mut_rate = length(which(df[,4]== 1))/dim(df)[1],
                p.values= tes$p.value)}
        return(re_p)
    })		
}, mc.cores = 5)

snv_re <- lapply(1:length(snv_cor_ge), function(i){
	df <- do.call(rbind, snv_cor_ge[[i]])
	return(df)
	})

snv_mis_cor <- do.call(rbind, snv_re)
snv_mis_cor <- na.omit(snv_mis_cor)
re <- lapply(unique(snv_mis_cor$sig), function(i){
    df <- subset(snv_mis_cor, sig == i)
    df$fdr <- p.adjust(df$p.value, method = "BH")
    return(df)
})
snv_mis_cor <- do.call(rbind, re)

snv_mis_cor_signi <- subset(snv_mis_cor, fdr < 0.1)    
write.csv(snv_mis_cor_signi, file = "~/snv_mis_correlation_signature_gene_result_significant0.1.csv", row.names = F)



#############################################################################################################################################################################
#frame shift mutation 
genelist <- read.csv("~/intersect_genelist.csv")
gene_matrix <- read.table("~/reference/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv.gz", header = T, sep = "\t") #20531 11070
gene <- gene_matrix
gene$gene_id <- str_split(gene$gene_id, '[|]', simplify = T)[,1]
gene <- gene[!duplicated(gene$gene_id),]
rownames(gene) <- gene$gene_id
gene <- gene[,-1]
gene <- as.data.frame(t(gene))
gene <- log2(gene + 1)
load("~/hyper_hypo_coef_H.RData")
load("~/TCGA_frame_shift.RData")
snv_domin <- snv_fram[colnames(snv_fram) %in% genelist$x]
snv_domin$sample_ID <- str_sub(rownames(snv_domin), 1, 12)
cof_snv <- merge(hyper_hypo_coef, snv_domin)

#correlation
snv_cor_ge <- mclapply(colnames(cof_snv)[13:length(cof_snv)], function(p){
    sig <- lapply(colnames(cof_snv)[2:11], function(q){
        cat(paste(p,q), "\n")
        df <- cof_snv[c(q, "sample_ID", "cancertype", p)]
        df[,4] <- as.factor(df[,4])
        if(dim(table(df[,4])) == 1|length(which(df[,4] == 1)) < 10){
        re_p <- data.frame(sig = q,
                gene = p, 
                mut_rate = NA,
                p.values= NA) 
        }else{
            tes <- wilcox.test(df[,1] ~ df[,4])
            re_p <- data.frame(sig = q,
                gene = p, 
                mut_rate = length(which(df[,4]== 1))/dim(df)[1],
                p.values= tes$p.value)}
        return(re_p)
    })      
}, mc.cores = 5)

snv_re <- lapply(1:length(snv_cor_ge), function(i){
    df <- do.call(rbind, snv_cor_ge[[i]])
    return(df)
    })

snv_fram_cor <- do.call(rbind, snv_re)
snv_fram_cor <- na.omit(snv_fram_cor)
re <- lapply(unique(snv_fram_cor$sig), function(i){
    df <- subset(snv_fram_cor, sig == i)
    df$fdr <- p.adjust(df$p.value, method = "BH")
    return(df)
})
snv_fram_cor <- do.call(rbind, re)

snv_fram_cor_signi <- subset(snv_fram_cor, fdr < 0.1)
write.csv(snv_fram_cor_signi, file = "~/snv_frame_correlation_signature_gene_result_significant_fdr0.1.csv", row.names = F)
