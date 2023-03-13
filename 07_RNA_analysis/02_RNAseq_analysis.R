# Notes here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\RNAseq\genome_paper_revisions\gene_exp.r 

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(scales)
setwd("/ssd3/Mar_genome_analysis/RNAseq/paper_revisions")

# sample data
    samples_data <- read.table("samples_to_use.txt", header = TRUE) %>%
        mutate(new_name = paste(clam, tissue, sep = "_"))
    # samples_data_healthy <- samples_data %>% filter(heme != "BTN", clam != "HEL_4", clam != "HEL_5")
    samples <- samples_data %>% pull(name) %>% as.vector()
    # samples_healthy <- samples_data_healthy %>% pull(name) %>% as.vector()

# Gene name conversion (these are in same order)
    gene_names <- read.table("/ssd3/Mar_genome_analysis/genomes/maker/annotation/mya_gene_table.tsv") %>%
        select(mya_gene, hit)
    head(gene_names)

# polymerase genes
    # list of genes with hits in cosmic and uniprot
    gene_names2 <- read.table("/ssd3/Mar_genome_analysis/genomes/maker/annotation/mya_gene_table.tsv") %>%
        separate(hit, into = c("uniprot_gene", NA), sep = "_", remove = FALSE) %>%
        filter(uniprot_gene != "uncharacterized") %>%
        separate(uniprot_gene, into = c("uniprot_gene", NA), sep = "-lik")
    
    # polymerase genes
    pol_genes <- read.table("polymerase_genes.txt")%>% pull(V1) %>% as.vector()

    #intersections
    filter(gene_names2, uniprot_gene %in% pol_genes) # 9
    filter(gene_names2, cosmic_hit %in% pol_genes) # 16
    gene_names3 <- gene_names2 %>%
        separate(cosmic_hit, into = c("cosmic_hit", NA), sep = "_")
    filter(gene_names3, cosmic_hit %in% pol_genes) # 29
    filter(gene_names3, cosmic_hit %in% pol_genes) %>% select(cosmic_hit, uniprot_gene)

    # list of least stringent - can come back and cut later
    pol_mya_genes <- filter(gene_names3, cosmic_hit %in% pol_genes) %>% pull(hit) %>% as.vector()

# other genes involved in stability/repair/telomeres
    other_mya_genes <- read.table("other_genes.txt") %>% pull(V1) %>% as.vector()

# LOAD RNASEQ DATA
    for(sample in samples){
        file <- paste0("/ssd3/RNAseq/STAR_genome_output/",sample, "_ReadsPerGene.out.tab")
        file_load <- read.table(file, sep="\t", skip = 4, col.names = c("mya_gene", sample, "fwd", "rev")) %>%
            full_join(gene_names) %>%
            column_to_rownames(var="hit") %>%
            select(-fwd,-rev, -mya_gene)
        if(sample == samples[1]){
            print("first")
            rnaseq <- file_load 
        } else {
            print(sample)
            rnaseq <- cbind(rnaseq,file_load)
        }
    }
# rename with better columns
    colnames(rnaseq) <- samples_data$new_name
    head(rnaseq)
    #rnaseq_healthy <- rnaseq[,as.vector(samples_data_healthy$new_name)]

# Normalize rnaseq data frame
    norm_value <- max(colSums(rnaseq))/colSums(rnaseq)
    rnaseq_norm <- as.data.frame(t(t(as.matrix(rnaseq))*norm_value))
    head(rnaseq_norm)

# plot expression of select genes
    samples_data2 <- samples_data
    levels(samples_data2$heme) <- c(levels(samples_data2$heme), "MarBTN", "hemocyte", "solid tissue")
    samples_data2[samples_data2[, "heme"] == "BTN", "heme"] <- "MarBTN"
    samples_data2[samples_data2[, "heme"] == "heme", "heme"] <- "hemocyte"
    samples_data2[samples_data2[, "heme"] == "not", "heme"] <- "solid tissue"
    samples_data2
    rnaseq[other_mya_genes ,]
    rnaseq[pol_mya_genes,]
    select_genes_expr <- rnaseq_norm[c(pol_mya_genes, other_mya_genes),] %>%
        t() %>%
        cbind(samples_data2)
    select_genes_expr[select_genes_expr[, "heme"] == "not", "heme"] <- factor("tissue")
    gene_means <- select_genes_expr %>% 
        select(-name, -code, -clam, -tissue, -new_name) %>%
        group_by(heme) %>%
        summarise(across(everything(), list(mean)))
    gene_sds <- select_genes_expr %>% 
        select(-name, -code, -clam, -tissue, -new_name) %>%
        group_by(heme) %>%
        summarise(across(everything(), list(sd)))
    gene_upper <- (gene_means + gene_sds) %>%
        select(-heme)
    gene_lower <- (gene_means - gene_sds) %>%
        select(-heme)

# pdf("relative_expr_initial.pdf")
# for(i in c(pol_mya_genes, other_mya_genes)){
#     cat(i)
#     grouping1 <- group_by(select_genes_expr, heme) %>%
#         summarise(mean = mean(get(i)), sd = sd(get(i)))
#     plot1 <- grouping1 %>%
#         ggplot() +
#             #geom_bar(data = gene_means, aes(x = heme, y = get(paste0(i,"_1")))) +
#             #geom_dotplot(aes(x = as.factor(heme), y = as.numeric(get(i)))) + 
#             geom_bar(data = grouping1, aes(x = heme, y = mean), stat='identity') + # , fill = heme
#             geom_point(data = select_genes_expr, aes(x = heme, y = get(i))) +
#             ggtitle(i) +
#             xlab("BTN, hemocyte or tissue") +
#             ylab("reads per gene, normalized") +
#             theme_classic()
#     print(plot1)
# }
# dev.off()

select_genes_expr2 <- select_genes_expr %>%
    select(-name, -code, -clam, -tissue, -new_name) %>%
    pivot_longer(!heme, names_to = "gene", values_to = "expression")
select_genes_expr_grouped <- select_genes_expr2 %>%
    group_by(heme, gene) %>%
    summarise(mean = mean(expression), sd = sd(expression))


# pdf("relative_expr_2.pdf")
#     ggplot() +
#         geom_bar(data = filter(select_genes_expr_grouped, gene %in% pol_mya_genes),
#             aes(x = heme, y = mean), stat='identity') +
#         geom_errorbar(data = filter(select_genes_expr_grouped, gene %in% pol_mya_genes),
#             aes(x = heme, ymin = mean-sd, ymax = mean+sd)) +
#         geom_point(data = filter(select_genes_expr2, gene %in% pol_mya_genes),
#             aes(x = heme, y = expression)) +
#         ggtitle("Polymerase genes, relative") +
#         xlab("BTN, hemocyte or tissue") +
#         ylab("reads per gene, normalized") +
#         theme_classic()+
#         facet_wrap(~gene, ncol = 5)
#     ggplot() +
#         geom_bar(data = filter(select_genes_expr_grouped, gene %in% pol_mya_genes),
#             aes(x = heme, y = mean), stat='identity') +
#         geom_errorbar(data = filter(select_genes_expr_grouped, gene %in% pol_mya_genes),
#             aes(x = heme, ymin = mean-sd, ymax = mean+sd)) +
#         geom_point(data = filter(select_genes_expr2, gene %in% pol_mya_genes),
#             aes(x = heme, y = expression)) +
#         ggtitle("Polymerase genes, scaled") +
#         xlab("BTN, hemocyte or tissue") +
#         ylab("reads per gene, normalized") +
#         theme_classic()+
#         facet_wrap(~gene, ncol = 5, scales = "free")
#     ggplot() +
#         geom_bar(data = filter(select_genes_expr_grouped, gene %in% other_mya_genes),
#             aes(x = heme, y = mean), stat='identity') +
#         geom_errorbar(data = filter(select_genes_expr_grouped, gene %in% other_mya_genes),
#             aes(x = heme, ymin = mean-sd, ymax = mean+sd)) +
#         geom_point(data = filter(select_genes_expr2, gene %in% other_mya_genes),
#             aes(x = heme, y = expression)) +
#         ggtitle("Genome stability genes, relative") +
#         xlab("BTN, hemocyte or tissue") +
#         ylab("reads per gene, normalized") +
#         theme_classic() +
#         facet_wrap(~gene, ncol = 3)
#     ggplot() +
#         geom_bar(data = filter(select_genes_expr_grouped, gene %in% other_mya_genes),
#             aes(x = heme, y = mean), stat='identity') +
#         geom_errorbar(data = filter(select_genes_expr_grouped, gene %in% other_mya_genes),
#             aes(x = heme, ymin = mean-sd, ymax = mean+sd)) +
#         geom_point(data = filter(select_genes_expr2, gene %in% other_mya_genes),
#             aes(x = heme, y = expression)) +
#         ggtitle("Genome stability genes, scaled") +
#         xlab("BTN, hemocyte or tissue") +
#         ylab("reads per gene, normalized") +
#         theme_classic() +
#         facet_wrap(~gene, ncol = 3, scales = "free")
# dev.off()


# Deseq on BTN vs hemocytes (only those from same run)
    deseq_data_BTN <- DESeqDataSetFromMatrix(countData = rnaseq,
                                colData = samples_data,
                                design = ~ heme)
    deseq_run_BTN <- DESeq(deseq_data_BTN)
    resultsNames(deseq_run_BTN) 
    BTNvHeme <- results(deseq_run_BTN, name="heme_heme_vs_BTN") %>% as.data.frame()
    BTNvHTiss <- results(deseq_run_BTN, name="heme_not_vs_BTN") %>% as.data.frame()
    
    BTNvHeme[c(pol_mya_genes, other_mya_genes),] %>%
        select(log2FoldChange, pvalue, padj, stat)

    nrow(BTNvHeme[pol_mya_genes,]) # 29
    BTNvHeme[pol_mya_genes,] %>%
        select(log2FoldChange, pvalue, padj, stat) %>%
        arrange(padj)

    BTNvHeme[pol_mya_genes,] %>%
        filter(log2FoldChange < 0, padj < 0.05) %>%
        arrange(padj) %>%
        select(log2FoldChange, pvalue, padj, stat)
    BTNvHeme[pol_mya_genes,] %>%
        filter(log2FoldChange > 0, padj < 0.05)

    BTNvHTiss[pol_mya_genes,] %>%
        filter(log2FoldChange < 0, padj < 0.05)
    BTNvHTiss[pol_mya_genes,] %>%
        filter(log2FoldChange > 0, padj < 0.05)




 BTNvHeme[other_mya_genes,] %>%
        select(log2FoldChange, pvalue, padj, stat) %>%
        filter(padj < 0.05) %>%
        arrange(padj)

 BTNvHTiss[other_mya_genes,] %>%
        select(log2FoldChange, pvalue, padj, stat) %>%
        filter(padj < 0.05) %>%
        arrange(padj)

######################################################
# PCA and heir clustering


# DESEQ TO DETERMINE TISSUE-SPECIFIC EXPRESSION
    samples_data <- read.table("/ssd3/RNAseq/STAR_genome_output/28samples/sample_data.txt", header = TRUE) %>%
        mutate(new_name = paste(clam, tissue, sep = "_"))
    samples_data_healthy <- samples_data %>% filter(heme != "BTN", clam != "HEL_4", clam != "HEL_5")
    samples <- samples_data %>% pull(name) %>% as.vector()
    samples_healthy <- samples_data_healthy %>% pull(name) %>% as.vector()
    rnaseq_healthy <- rnaseq[,as.vector(samples_data_healthy$new_name)]

    # Blank df
        tissue_specific_top100 <- data.frame(gene = c(), tissue = c())
    #tissue list
        tissues <- samples_data_healthy$tissue %>% unique() %>% as.vector()
    for(i in tissues){
    # new column specififying target
        samples_data_healthy2 <- samples_data_healthy %>%
            mutate(target = get(i))
    # run deseq and extract top 100 genes expressed in each tissue
        tissue_specific_top100 <- DESeqDataSetFromMatrix(
                                    countData = rnaseq_healthy,
                                    colData = samples_data_healthy2,
                                    design = ~ clam + target) %>%
            DESeq() %>%
            results(contrast=c("target","not",i)) %>%
            as.data.frame() %>%
            arrange(stat) %>%
            mutate(tissue = i) %>%
            rownames_to_column(var = "gene") %>%
            select(gene, tissue) %>%
            head(n=100) %>%
            rbind(tissue_specific_top100)
    }
    tissue_specific_list <- unique(tissue_specific_top100$gene)
    rnaseq_norm_tissuespecific <- rnaseq_norm[tissue_specific_list,]

# Run deseq by all tissues
    samples_data <- read.table("samples_to_use.txt", header = TRUE) %>%
        mutate(new_name = paste(clam, tissue, sep = "_"))
    deseq_data <- DESeqDataSetFromMatrix(countData = rnaseq,
                                colData = samples_data,
                                design = ~ tissue)
    deseq_data
    deseq_run1 <- DESeq(deseq_data)



    # The palette with grey:
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # The palette with black:
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    size1=8
########################################################
## after running sections above, make output figures for paper
# Main figure: PCA, polymerase volcano, and p53/mortalin/BRCA1
    #PCA
    deseq_vsd <- vst(deseq_run1)
    pca_plot <- plotPCA(deseq_vsd, intgroup = c("tissue")) +
        theme_classic() +
        scale_colour_manual(values=cbPalette) +
        labs(color=NULL) +
        #ggtitle("vst method, defaults") +
        theme(
            aspect.ratio = 1,
            #plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size = size1, face = "bold"),
            axis.title = element_text(size = size1, face = "bold"),
            text = element_text(size = size1, face = "bold"),
            #legend.title = element_text("none"),
            legend.text = element_text(size = size1, face = "bold")
        )
    # polymerase volcano
    poly_volcano <- BTNvHeme[pol_mya_genes,] %>%
        ggplot(mapping = aes(x = -log2FoldChange, y = -log10(padj))) +
            geom_point() +
            xlab("Log2 Fold Change (BTN vs hemocytes)") +
            ylab("-Log10 adj p-value") +
            geom_hline(yintercept =-log10(0.05)) +        
            #ggtitle("Differential expression of polymerase genes (29)") +
            xlim(-8,8) +
            theme_classic() +
            theme(
                aspect.ratio = 0.95,
                #plot.title = element_text(hjust = 0.5),
                axis.text = element_text(size = size1, face = "bold"),
                axis.title = element_text(size = size1, face = "bold"),
                text = element_text(size = size1, face = "bold"),
                legend.title = element_text("none")
            )
    poly_volcano_skinny <- BTNvHeme[pol_mya_genes,] %>%
        ggplot(mapping = aes(x = -log2FoldChange, y = -log10(padj))) +
            geom_point() +
            xlab("Log2 Fold Change (BTN vs hemocytes)") +
            ylab("-Log10 adj p-value") +
            geom_hline(yintercept =-log10(0.05)) +        
            #ggtitle("Differential expression of polymerase genes (29)") +
            xlim(-8,8) +
            theme_classic() +
            theme(
                aspect.ratio = 0.5,
                plot.title = element_text(hjust = 0.5),
                axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 16, face = "bold"),
                text = element_text(size = 16, face = "bold"),
                legend.title = element_text("none")
            )
    # p53/mortalin/BRCA1
    three_genes <- filter(select_genes_expr_grouped, gene %in% c("P63-like_1", "BRCA1-like_1", "HSP6-like_1"))
    three_genes_plot <- ggplot() +
        geom_bar(data = three_genes,
            aes(x = heme, y = mean), stat='identity') +
        geom_errorbar(data = three_genes,
            aes(x = heme, ymin = mean-sd, ymax = mean+sd)) +
        geom_point(data = filter(select_genes_expr2, gene %in% c("P63-like_1", "BRCA1-like_1", "HSP6-like_1")),
            aes(x = heme, y = expression)) +
        xlab("BTN, hemocyte or tissue") +
        ylab("Reads per gene") +
        theme_classic() +
            theme(axis.title.x=element_blank()) +
        facet_wrap(~gene, ncol = 1, scales = "free")
    p53_plot <- ggplot(filter(select_genes_expr_grouped, gene =="P63-like_1")) +
        geom_bar(aes(x = heme, y = mean), stat='identity') +
        geom_errorbar(aes(x = heme, ymin = mean-sd, ymax = mean+sd)) +
        geom_point(data = filter(select_genes_expr2, gene =="P63-like_1"), aes(x = heme, y = expression), size = 0.5) +
        ylab("") +
        ggtitle("TP53") + 
        scale_y_continuous(label=scientific_format())+
        theme_classic() +
        theme(
                aspect.ratio = 1,
                axis.text = element_text(size = size1, face = "bold"),
                axis.title = element_text(size = size1, face = "bold"),
                #axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(size = size1, hjust = 0.5),
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                text = element_text(size = size1, face = "bold"),
                legend.title = element_text("none")
            )
    mort_plot <- ggplot(filter(select_genes_expr_grouped, gene =="HSP6-like_1")) +
        geom_bar(aes(x = heme, y = mean), stat='identity') +
        geom_errorbar(aes(x = heme, ymin = mean-sd, ymax = mean+sd)) +
        geom_point(data = filter(select_genes_expr2, gene =="HSP6-like_1"), aes(x = heme, y = expression), size = 0.5) +
        ylab("Reads/gene") +
        ggtitle("HSPA9") + 
        scale_y_continuous(label=scientific_format())+
        theme_classic() +
        theme(
                aspect.ratio = 1,
                axis.text = element_text(size = size1, face = "bold"),
                axis.title = element_text(size = size1, face = "bold"),
                #axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(size = size1, hjust = 0.5),
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                text = element_text(size = size1, face = "bold"),
                legend.title = element_text("none")
            )
    brca_plot <- ggplot(filter(select_genes_expr_grouped, gene =="BRCA1-like_1")) +
        geom_bar(aes(x = heme, y = mean), stat='identity') +
        geom_errorbar(aes(x = heme, ymin = mean-sd, ymax = mean+sd)) +
        geom_point(data = filter(select_genes_expr2, gene =="BRCA1-like_1"), aes(x = heme, y = expression), size = 0.5) +
        ylab("") +
        ggtitle("BRCA1") + 
        scale_y_continuous(label=scientific_format())+
        theme_classic() +
        theme(
                aspect.ratio = 1,
                axis.text = element_text(size = size1, face = "bold"),
                axis.title = element_text(size = size1, face = "bold"),
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(size = size1, hjust = 0.5),
                axis.title.x=element_blank(),
                text = element_text(size = size1, face = "bold"),
                legend.title = element_text("none")
            )


# Volcano - BTN vs tissue for polymerases
    poly_volcano_tissue <- BTNvHTiss[pol_mya_genes,] %>%
        ggplot(mapping = aes(x = -log2FoldChange, y = -log10(padj))) +
            geom_point() +
            xlab("Log2 Fold Change (BTN vs non-hemocyte tissues)") +
            ylab("-Log10 adj p-value") +
            geom_hline(yintercept =-log10(0.05)) +        
            #ggtitle("Differential expression of polymerase genes (29)") +
            xlim(-8,8) +
            theme_classic() +
            theme(
                aspect.ratio = 0.5,
                plot.title = element_text(hjust = 0.5),
                axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 16, face = "bold"),
                text = element_text(size = 16, face = "bold"),
                legend.title = element_text("none")
            )

# Positively selected genes and related genes
    possel_genes <- filter(select_genes_expr_grouped, gene %in% other_mya_genes[c(3, 11:13)])
    possel_genes_plot <- ggplot() +
        geom_bar(data = possel_genes,
            aes(x = heme, y = mean), stat='identity') +
        geom_errorbar(data = possel_genes,
            aes(x = heme, ymin = mean-sd, ymax = mean+sd)) +
        geom_point(data = filter(select_genes_expr2, gene %in% other_mya_genes[c(3, 11:13)]),
            aes(x = heme, y = expression)) +
        xlab("BTN, hemocyte or tissue") +
        ylab("Reads per gene") +
        theme_classic() +
            theme(axis.title.x=element_blank()) +
        facet_wrap(~gene, ncol = 2, scales = "free")
    possel_genes2 <- filter(select_genes_expr_grouped, gene %in% other_mya_genes[c(3, 11:15, 1, 2, 4, 5, 6)])
    possel_genes_plot2 <- ggplot() +
        geom_bar(data = possel_genes2,
            aes(x = heme, y = mean), stat='identity') +
        geom_errorbar(data = possel_genes2,
            aes(x = heme, ymin = mean-sd, ymax = mean+sd)) +
        geom_point(data = filter(select_genes_expr2, gene %in% other_mya_genes[c(3, 11:15, 1, 2, 4, 5, 6)]),
            aes(x = heme, y = expression)) +
        xlab("BTN, hemocyte or tissue") +
        ylab("Reads per gene") +
        theme_classic() +
            theme(axis.title.x=element_blank()) +
        facet_wrap(~gene, ncol = 3, scales = "free")

pdf("expression_supp_figures.pdf")
poly_volcano_skinny
poly_volcano_tissue
possel_genes_plot
possel_genes_plot2 
# Supp - Heir clustering
    annotations <- samples_data %>%
        column_to_rownames(var = "new_name") %>%
        select(tissue)
    rnaseq_norm_tissuespecific %>% 
        pheatmap(
        cluster_rows = FALSE,
        show_rownames = FALSE,
        cluster_cols = TRUE,
        clustering_distance_cols = "canberra",
        scale = "row",
        annotation_col = annotations,
        annotation_colors = cbPalette[1:7]
        )
dev.off()

# multi-panel main figure
pdf("expression_main_fig.pdf", width = 8, height = 3)
#grid.arrange(pca_plot, poly_volcano, three_genes_plot, nrow = 1)
plots <- list(pca_plot, poly_volcano, p53_plot, mort_plot, brca_plot)
#plots <- c(pca_plot, poly_volcano, three_genes_plot)
grid.arrange(
  grobs = plots,
  widths = c(3.35, 2.5, 0.9),
  layout_matrix = rbind(c(1,2,3),
                        c(1,2,3),
                        c(1,2,3),
                        c(1,2,4),
                        c(1,2,4),
                        c(1,2,4),
                        c(1,2,5),
                        c(1,2,5),
                        c(1,2,5),
                        c(1,2,5),
                        c(1,2,5))
)
dev.off()


