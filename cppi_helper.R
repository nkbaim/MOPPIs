#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Du Yang
# @Date: 2025-01-06 21:50:24
# @LastEditTime: 2025-01-10 09:04:22
# @LastEditors: Du Yang
# @Description: TODO
#


# 01. extract ppi -----------------------------
extract_interactome <- function(
    queried_gene,
    ppi_dir = "/mnt/dellfs/home/duyang/nsfc2021/01.pub_data/pina") {
    interactome <- import(paste0(ppi_dir, "/Homo_sapiens.symbols.sif"), format = "txt") %>%
        as_tibble()

    # step 1, 导出所有与target_node 存在互作节点
    net1 <- interactome %>%
        filter(source == queried_gene | target == queried_gene) %>%
        filter(source != target)

    net_nodes <- unique(c(net1$source, net1$target))

    # step 2, 导出激活点之间的链接
    net2 <- interactome %>% filter(source %in% net_nodes, target %in% net_nodes)
    net <- rbind(net1, net2) %>%
        filter(source != target)

    # step 3, 去重
    net3 <- net %>%
        mutate(
            edges = apply(., 1, function(x) {
                paste0(sort(x), collapse = "_")
            })
        ) %>%
        filter(!duplicated(edges))

    network <- list(
        queried_gene = queried_gene,
        nodes = data.frame(name = net_nodes),
        # links = net3 %>% mutate(source = match(source, net_nodes), target = match(target, net_nodes))
        links = net3
    )
    return(network)
}

# 02. context specific co-dependency ---------------------------
context_dependency <- function(
    network = NULL,
    cancer_context = "LUAD",
    essential_score = 0.5,
    dep_dir = "/mnt/dellfs/home/duyang/nsfc2021/01.pub_data/depmap/") {
    # 1. load gene score
    models <- import(paste0(dep_dir, "Model.csv")) %>%
        as_tibble() %>%
        filter(OncotreeCode == cancer_context) %>%
        select(ModelID)

    # 2. get context specific Essentiality
    essential <- import(paste0(dep_dir, "CRISPRGeneDependency.csv")) %>%
        as_tibble() %>%
        filter(V1 %in% models$ModelID) %>%
        set_names(gsub(" .*", "", colnames(.))) %>%
        select(-V1) %>%
        apply(., 2, function(x) {
            median(x)
        })

    nodes_essential <- network$nodes %>%
        inner_join(data.frame(
            name = names(essential),
            essential = essential
        )) %>%
        filter(essential > essential_score)

    # 3. get context specific co-dependency
    links_by_essential <- network$links %>%
        filter(
            source %in% nodes_essential$name,
            target %in% nodes_essential$name
        )

    dep_matrix <- import(paste0(dep_dir, "CRISPRGeneDependency.csv")) %>%
        as_tibble() %>%
        filter(V1 %in% models$ModelID) %>%
        set_names(gsub(" .*", "", colnames(.))) %>%
        select(-V1) %>%
        as.matrix()

    links_essential <- links_by_essential %>%
        mutate(
            cor = map2_dbl(source, target, ~ cor(
                dep_matrix[, .x],
                dep_matrix[, .y],
                method = "spearman", use = "pairwise.complete.obs"
            )),
            p_value = map2_dbl(source, target, ~ cor.test(
                dep_matrix[, .x],
                dep_matrix[, .y],
                method = "spearman", use = "pairwise.complete.obs"
            )$p.value)
        ) %>%
        mutate(
            p.adj = p.adjust(p_value)
        ) %>%
        filter(p.adj < 0.05)

    network$nodes_essential <- nodes_essential
    network$links_essential <- links_essential

    return(network)
}


# functions
mut_interaction <- function(
    network = NULL,
    data_source = "TCGA-LUAD",
    mut_data = NULL) {
    input <- import(mut_data, format = "txt")

    # 在这里整体的剔除高/突变样本;
    mut_count <- input %>%
        column_to_rownames("attrib_name") %>%
        apply(., 2, sum)

    removedsamples <- boxplot(mut_count)$out
    input <- input %>% select(-names(removedsamples))


    # net_edges <- apply(network$links, 1, function(x) paste0(sort(x), collapse = "_"))

    somatic_interactions <- somaticInteractions(input,
        top = 0.01, cancer_drivers = network$nodes$name
    ) %>%
        as_tibble() %>%
        select(gene1, gene2, Event, pValue) %>%
        mutate(
            data_source = cancer_context,
            ppa_type = "somatic_interaction"
        ) %>%
        mutate(
            edges = apply(., 1, function(x) {
                paste0(sort(c(x[1], x[2])), collapse = "_")
            })
        ) %>%
        filter(
            edges %in% network$links$edges
        ) %>%
        filter(
            !duplicated(edges)
        ) %>%
        mutate(p.adj = p.adjust(pValue)) %>%
        filter(p.adj < 0.05)


    network$somatic_interactions <- somatic_interactions

    return(network)
}



somaticInteractions <- function(df_mat, top = 0.05, cancer_drivers = NULL, pvalue = 0.05) {
    library(tidyverse)
    if (is.null(cancer_drivers)) {
        cancer_drivers <- rio::import("/mnt/dellfs/pub/data/CancerDriverGenes/Pancancer_driver_genes.txt")$Hugo_Symbol
    }

    mutMat <- df_mat %>%
        as_tibble() %>%
        filter(attrib_name %in% cancer_drivers) %>%
        column_to_rownames("attrib_name") %>%
        filter(
            apply(., 1, function(x) {
                sum(x) / length(x)
            }) >= top
        ) %>%
        t()

    interactions <- sapply(1:ncol(mutMat), function(i) {
        sapply(
            1:ncol(mutMat),
            function(j) {
                f <- try(fisher.test(mutMat[, i], mutMat[, j]), silent = TRUE)
                if (class(f) == "try-error") {
                    NA
                } else {
                    ifelse(f$estimate > 1, -log10(f$p.val), log10(f$p.val))
                }
            }
        )
    })
    oddsRatio <- oddsGenes <- sapply(1:ncol(mutMat), function(i) {
        sapply(
            1:ncol(mutMat),
            function(j) {
                f <- try(fisher.test(mutMat[, i], mutMat[, j]), silent = TRUE)
                if (class(f) == "try-error") {
                    f <- NA
                } else {
                    f$estimate
                }
            }
        )
    })

    rownames(interactions) <- colnames(interactions) <- rownames(oddsRatio) <- colnames(oddsRatio) <- colnames(mutMat)
    sigPairs <- which(x = 10^-abs(interactions) < 1, arr.ind = TRUE)
    sigPairs2 <- which(x = 10^-abs(interactions) >= 1, arr.ind = TRUE)
    if (nrow(sigPairs) < 1) {
        stop("No meaningful interactions found.")
    }
    sigPairs <- rbind(sigPairs, sigPairs2)
    sigPairsTbl <- data.table::rbindlist(lapply(
        X = seq_along(1:nrow(sigPairs)),
        function(i) {
            x <- sigPairs[i, ]
            g1 <- rownames(interactions[x[1], x[2], drop = FALSE])
            g2 <- colnames(interactions[x[1], x[2], drop = FALSE])
            tbl <- as.data.frame(table(apply(X = mutMat[, c(
                g1,
                g2
            ), drop = FALSE], 1, paste, collapse = "")))
            combn <- data.frame(t(tbl$Freq))
            colnames(combn) <- tbl$Var1
            pval <- 10^-abs(interactions[x[1], x[2]])
            fest <- oddsRatio[x[1], x[2]]
            d <- data.table::data.table(
                gene1 = g1, gene2 = g2,
                pValue = pval, oddsRatio = fest
            )
            d <- cbind(d, combn)
            d
        }
    ), fill = TRUE)
    sigPairsTbl <- sigPairsTbl[!gene1 == gene2]
    sigPairsTbl[is.na(sigPairsTbl)] <- 0
    sigPairsTbl$Event <- ifelse(test = sigPairsTbl$oddsRatio >
        1, yes = "Co_Occurence", no = "Mutually_Exclusive")
    sigPairsTbl$pair <- apply(
        X = sigPairsTbl[, .(gene1, gene2)],
        MARGIN = 1, FUN = function(x) {
            paste(sort(unique(x)),
                collapse = ", "
            )
        }
    )
    sigPairsTbl[, `:=`(event_ratio, `01` + `10`)]
    sigPairsTbl[, `:=`(event_ratio, paste0(`11`, "/", event_ratio))]
    # print(head(sigPairsTbl))
    # sigPairsTblSig <- sigPairsTbl %>%
    #     group_by(gene2) %>%
    #     mutate(
    #         p.adj = p.adjust(pValue, method = "BH")
    #     ) %>%
    #     ungroup() %>%
    #     filter(!duplicated(pair), p.adj < pvalue)
    return(sigPairsTbl)
}


coExp <- function(
    network = NULL, exp_data = NULL, rho_threshold = 0.2,
    pvalue = 0.05, data_type = "mRNA", data_source = "TCGA_LUAD") {
    df_mat <- import(exp_data, format = "txt")
    # update 2023-12-27, 去除空值；
    expMat <- df_mat %>%
        as_tibble() %>%
        filter(attrib_name %in% network$nodes$name) %>%
        filter(!duplicated(attrib_name)) %>%
        column_to_rownames("attrib_name") %>%
        mutate_all(
            function(x) {
                ifelse(x == 0, NA, x)
            }
        ) %>%
        filter(apply(., 1, function(x) {
            sum(is.na(x)) / length(x) == 0
        }))

    corr <- WGCNA::corAndPvalue(t(expMat), method = "spearman")

    # net_edges <- apply(network$links, 1, function(x) paste0(sort(x), collapse = "_"))

    corr.cor <- corr$cor %>%
        as.data.frame() %>%
        rownames_to_column("gene1") %>%
        pivot_longer(
            cols = -gene1,
            names_to = "gene2",
            values_to = "rho"
        )

    # update 2023-12-26, filter with FDR, use p.adj to replace p
    corr.pval <- corr$p %>%
        as.data.frame() %>%
        rownames_to_column("gene1") %>%
        pivot_longer(
            cols = -gene1,
            names_to = "gene2",
            values_to = "pval"
        )

    corr.df <- corr.cor %>%
        inner_join(corr.pval) %>%
        mutate(
            edges = apply(., 1, function(x) {
                paste0(sort(c(x[1], x[2])), collapse = "_")
            })
        ) %>%
        filter(
            edges %in% network$links$edges
        ) %>%
        filter(
            !duplicated(edges)
        ) %>%
        mutate(
            p.adj = p.adjust(pval, method = "BH")
        ) %>%
        filter(
            p.adj < pvalue,
            abs(rho) > rho_threshold
        ) %>%
        mutate(
            data_source = data_source,
            data_type = data_type
        )

    if (is.null(network$coExp)) {
        network$coExp <- corr.df
    } else {
        network$coExp <- rbind(network$coExp, corr.df)
    }
    return(network)
}

# annotations

# prognosis
prognosis_anno <- function(network = NULL, cancer_context = "LUAD") {
    network$node_prognosis <- import("/mnt/dellfs/home/duyang/nsfc2021/02.result/07.prognosis/survival_mrna_2020-09-24.txt") %>%
        as_tibble() %>%
        filter(cancer_type == paste0("TCGA-", cancer_context)) %>%
        filter(gene %in% network$nodes$name, pvalue < 0.05, cutoff == "High50_Low50")

    return(network)
}

drug_anno <- function(network = NULL) {
    network$node_drugs <- import("/mnt/dellfs/home/duyang/nsfc2021/01.pub_data/drugs/Drugables.json") %>%
        as_tibble() %>%
        filter(tier %in% c("T1", "T2", "T3", "T4", "T5")) %>%
        filter(gene %in% network$nodes$name)

    return(network)
}

cancer_dirver_anno <- function(network = NULL) {
    network$node_cancer_dirvers <- import("/mnt/dellfs/pub/data/CancerDriverGenes/Pancancer_driver_genes.txt") %>%
        as_tibble() %>%
        mutate(
            role = ifelse(
                OG == "YES" & TSG == "NO", "Oncogene",
                ifelse(OG == "NO" & TSG == "YES", "Tumor Suppressor",
                    ifelse(OG == "YES" & TSG == "YES", "Oncogene|Tumor Suppressor", "Not clear")
                )
            )
        ) %>%
        select(Hugo_Symbol, role) %>%
        set_names(c("gene", "role")) %>%
        filter(gene %in% network$nodes$name)

    return(network)
}


# functions
mut_assoc <- function(
    network = NULL,
    mut_data = "Human__TCGA_LUAD__WUSM__Mutation__GAIIx__01_28_2016__BI__Gene__Firehose_MutSig2CV.cbt",
    data_type = "mRNA",
    data_source = "TCGA-LUAD",
    exp_data = "rna.exp_gene_revised.txt",
    top = 0.05, cancer_drivers = NULL) {
    # net_edges <- apply(network$links, 1, function(x) paste0(sort(x), collapse = "_"))

    input <- import(mut_data, format = "txt")

    # 在这里整体的剔除高/突变样本;
    mut_count <- input %>%
        column_to_rownames("attrib_name") %>%
        apply(., 2, sum)

    removedsamples <- boxplot(mut_count)$out
    input <- input %>% select(-names(removedsamples))


    if (is.null(cancer_drivers)) {
        cancer_drivers <- rio::import("/mnt/dellfs/pub/data/CancerDriverGenes/Pancancer_driver_genes.txt")$Hugo_Symbol
    }

    mutMat <- input %>%
        as_tibble() %>%
        filter(attrib_name %in% cancer_drivers, attrib_name %in% network$nodes$name) %>%
        column_to_rownames("attrib_name") %>%
        filter(
            apply(., 1, function(x) {
                sum(x) / length(x)
            }) >= top
        ) %>%
        rownames_to_column("gene") %>%
        pivot_longer(
            cols = -gene,
            names_to = "sample",
            values_to = "mut_status"
        )

    if (nrow(mutMat) == 0) {
        network$mut_assoc <- NULL
        return(network)
    }

    # 表达数据
    expMat <- import(exp_data, format = "txt") %>%
        as_tibble() %>%
        filter(attrib_name %in% network$nodes$name) %>%
        filter(!duplicated(attrib_name)) %>%
        column_to_rownames("attrib_name") %>%
        mutate_all(
            function(x) {
                ifelse(x == 0, NA, x)
            }
        ) %>%
        filter(apply(., 1, function(x) {
            sum(is.na(x)) / length(x) == 0
        })) %>%
        rownames_to_column("geneExp") %>%
        pivot_longer(
            cols = -geneExp,
            names_to = "sample",
            values_to = "expression"
        )

    mut_assoc_res <- expMat %>%
        inner_join(mutMat) %>%
        group_by(gene, geneExp) %>%
        summarise(
            # mean of logFC
            logFC = ifelse(
                mean(expression[mut_status == 1] - expression[mut_status == 0], na.rm = TRUE),
                mean(expression[mut_status == 1], na.rm = T) - mean(expression[mut_status == 0], na.rm = TRUE)
            ),
            P.Value = tryCatch(
                t.test(
                    expression[mut_status == 1],
                    expression[mut_status == 0],
                )$p.value,
                error = function(e) NA
            )
        ) %>%
        set_names(
            c("Gene_mutated", "Gene_associated", "logFC", "P.Value")
        ) %>%
        mutate(
            data_type = data_type,
            data_source = data_source
        ) %>%
        ungroup() %>%
        mutate(
            edges = apply(., 1, function(x) {
                paste0(sort(c(x[1], x[2])), collapse = "_")
            })
        ) %>%
        filter(
            edges %in% network$links$edges
        ) %>%
        filter(
            !duplicated(edges)
        ) %>%
        group_by(Gene_mutated) %>%
        mutate(
            p.adj = p.adjust(P.Value)
        ) %>%
        filter(p.adj < 0.05) %>%
        ungroup()


    if (nrow(mut_assoc_res) == 0) {
        network$mut_assoc <- NULL
    } else {
        network$mut_assoc <- mut_assoc_res
    }


    return(network)
}


# use random walk, to prouning the network


# a unique function to perform Context-specific PPI analysis
cppi <- function(
    queried_gene = "NDRG1",
    cancer_context = "LUAD",
    data_source = "TCGA-LUAD",
    mut_data = NULL,
    mrna_data = NULL,
    protein_data = NULL,
    phospho_data = NULL) {
    # read data
    # 1. load ppi network, all the data is based on the ppi network
    net <- extract_interactome(queried_gene)

    # 2. load mutational data
    # e.g. "/mnt/dellfs/pub/data/LinkedOmics/2023/TCGA-LUAD/Human__TCGA_LUAD__WUSM__Mutation__GAIIx__01_28_2016__BI__Gene__Firehose_MutSig2CV.cbt"
    if (!is.null(mut_data)) {
        message("mutation interaction")
        net <- mut_interaction(net, mut_data = mut_data, data_source = data_source)
    }

    # mrna_data = paste0(tcga_dir, "TCGA-LUAD/Human__TCGA_LUAD__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct")
    # 3. coexp data
    if (!is.null(mrna_data)) {
        message("mrna coexp ")
        net <- coExp(net, exp_data = mrna_data, data_source = data_source, data_type = "mRNA")
    }

    if (!is.null(protein_data)) {
        message("protein coexp ")
        net <- coExp(net, exp_data = protein_data, data_source = data_source, data_type = "protein")
    }

    if (!is.null(phospho_data)) {
        message("phospho coexp ")
        net <- coExp(net, exp_data = phospho_data, data_source = data_source, data_type = "phosphoProtein")
    }

    # 4. mut association
    if (!is.null(mrna_data) & !is.null(mut_data)) {
        message("mut mrna assoc ")
        net <- mut_assoc(net, mut_data = mut_data, exp_data = mrna_data, data_source = data_source, data_type = "mRNA")
    }

    if (!is.null(protein_data) & !is.null(mut_data)) {
        message("mut protein assoc ")
        net <- mut_assoc(net, mut_data = mut_data, exp_data = protein_data, data_source = data_source, data_type = "protein")
    }

    if (!is.null(phospho_data) & !is.null(mut_data)) {
        message("mut phospho assoc ")
        net <- mut_assoc(net, mut_data = mut_data, exp_data = phospho_data, data_source = data_source, data_type = "phosphoProtein")
    }

    # 5. co-enssential
    net <- context_dependency(net, cancer_context = cancer_context)

    # 6. drug targets
    net <- drug_anno(net)
    net <- cancer_dirver_anno(net)

    # 7. prognosis markers
    net <- prognosis_anno(net, cancer_context = cancer_context)

    return(net)
}



# graph random walk with restart
require(igraph)

## 对network, 创建权重graph，
create_weighted_network <- function(network = NULL) {
    # edge multi-omics count, 2,3,4,5,6,7
    # node weight, survival: 5, enssential: 5, cancer driver: 3, drug target: 2
    network$omics_links <- rbind(
            network$somatic_interactions %>% select(edges),
            network$coExp %>% select(edges),
            network$mut_assoc %>% select(edges),
            network$links_essential %>% select(edges)
        ) %>%
        group_by(edges) %>%
        summarise(weight = n()) %>%
        ungroup() %>%
        left_join(network$links)

    network$omics_nodes <- network$nodes %>% 
        filter((name %in% network$omics_links$source) | (name %in% network$omics_links$target)) %>%
        mutate(weight = 1) %>%
        mutate(weight = ifelse(name %in% network$nodes_essential$name, weight + 5, weight)) %>%
        mutate(weight = ifelse(name %in% network$node_prognosis$gene, weight + 5, weight)) %>%
        mutate(weight = ifelse(name %in% network$node_cancer_dirvers$gene, weight + 3, weight)) %>%
        mutate(weight = ifelse(name %in% network$node_drugs$gene, weight + 2, weight))

    return(network)

}

create_graph_from_edges <- function(network = NULL) {
    # 检查输入数据
    if (length(source) != length(target)) {
        stop("Source and target vectors must have the same length.")
    }

    # 创建边列表
    edges <- data.frame(source = network$omics_links$source, target = network$omics_links$target)

    # 创建图
    G <- graph_from_data_frame(edges, directed = FALSE)


    V(G)$name
    V(G)$weight <- node_weights

    # 添加边权重（如果提供）
    if (!is.null(edge_weights)) {
        if (length(edge_weights) != ecount(G)) {
            stop("Edge weights must have the same length as the number of edges.")
        }
        E(G)$weight <- edge_weights
    }

    return(G)
}


random_walk_with_weights <- function(G, start_node, node_weights, restart_prob = 0.5, max_steps = 1000) {
    n_nodes <- vcount(G)
    visit_prob <- rep(0, n_nodes) # 初始化访问概率
    current_node <- start_node # 当前节点
    visit_prob[current_node] <- 1 # 起始节点的访问概率为 1

    for (step in 1:max_steps) {
        # 以 restart_prob 的概率返回到起始节点
        if (runif(1) < restart_prob) {
            current_node <- start_node
        } else {
            # 获取当前节点的邻居
            neighbors <- neighbors(G, current_node)
            if (length(neighbors) > 0) {
                # 根据节点权重选择下一个节点
                neighbor_weights <- node_weights[neighbors]
                current_node <- sample(neighbors, 1, prob = neighbor_weights)
            }
        }
        visit_prob[current_node] <- visit_prob[current_node] + 1
    }

    # 归一化访问概率
    visit_prob <- visit_prob / sum(visit_prob)
    return(visit_prob)
}

# 3. 根据 RWR 分数修剪网络
prune_network <- function(G, rwr_scores, threshold) {
    # 根据 RWR 分数修剪网络
    # 参数:
    #   G: 输入的图
    #   rwr_scores: 节点的 RWR 分数
    #   threshold: 修剪的阈值
    # 返回值: 修剪后的图
    nodes_to_keep <- V(G)[rwr_scores > threshold]
    pruned_G <- induced_subgraph(G, nodes_to_keep)
    return(pruned_G)
}

# 4. 可视化网络
plot_network <- function(G, title) {
    # 可视化网络
    # 参数:
    #   G: 输入的图
    #   title: 图的标题
    plot(G,
        vertex.size = 10, vertex.label.cex = 0.8, edge.arrow.size = 0.5,
        main = title, vertex.color = "lightblue", edge.color = "gray"
    )
}
