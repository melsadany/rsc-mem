################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
## function for correcting file path to adapt to file system: Argon, Topaz, IDAS, or personal
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
## this is Muhammad's wrapper functions, aes for viz, and other customized scripts
## It also loads libraries of interest and makes sure they're installed
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(patchwork)
library(Seurat);library(Signac)
registerDoMC(30)
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/rsc-mem")
setwd(project.dir)
################################################################################
################################################################################
################################################################################
# cell.colors <- c("#e76cf2", "#04aef3", "#04bd7c", "#a3a700", "#fa766d", "#cab2d6", "#b4df8b", "#a6cfe2")
# names(cell.colors) <- c("ExN", "InN", "Astro", "Oligo", "Microglia", "Endo","VLMC","OPC")

### colors from Budhaditya
class.colors <- c("ExN"="#E76BF3", 
                  "InN"="#00B0F6", 
                  "Astro"="#00BF7D",
                  "Oligo"="#A3A500",
                  "Microglia"="#F8766D",
                  "Endo"="#CAB2D6",
                  "VLMC"="#B2DF8A",
                  "OPC"="#A6CEE3")
subclass.colors <- c("L2_3"="#1F77B4FF","L4"="#E377C2FF","L5"="#D62728FF","L6"="#6677d4",
                     "NP SUB"= "#720ba9","L2_3 RSP"="#FF7F0EFF","L4 RSP"="#2CA02CFF",
                     "Pvalb"="#17BECFFF","Sst"="#b7950b","Lamp5"="#98DF8AFF", "Sncg"="#FF9896FF","Vip"="#411eee",
                     "Oligo"="#da5e0e", "Astro"="#9EDAE5FF","Microglia"="#6495ED","OPC"="#FFBB78FF", 
                     "Endo"="#8C564BFF","VLMC"="#C49C94FF")

## NEUROeSTIMator genes
ne.genes <- c("Arc", "Btg2", "Coq10b", "Crem", "Dusp1", "Dusp5", "Egr1", "Egr3", 
              "Fbxo33", "Fos", "Fosl2", "Gadd45g", "Gmeb2", "Grasp", "Junb", 
              "Nr4a1", "Nr4a2", "Nr4a3", "Per1", "Rgs2", "Sertad1", "Tiparp")
################################################################################
################################################################################
################################################################################
## read data
# exp 1
sor.hc <- read_rds("data/raw/RSC_integrated_xenium_new.rds");gc()
sor.hc.clean <- GetAssayData(JoinLayers(sor.hc[["RNA"]]))
sor.hc.df <- sor.hc.clean %>% as.data.frame()
sor.hc.meta <- sor.hc@meta.data %>% 
  rownames_to_column("cell") %>%
  select(cell, orig.ident, condition, celltype, class)
# exp 2
tau.ctrl <- read_rds("data/raw/Tau_RSC_integrated_xenium_new.rds");gc()
tau.ctrl.clean <- GetAssayData(JoinLayers(tau.ctrl[["RNA"]]))
tau.ctrl.df <- tau.ctrl.clean %>% as.data.frame()
tau.ctrl.meta <- tau.ctrl@meta.data %>% 
  rownames_to_column("cell") %>%
  select(cell, orig.ident, condition, celltype, class)
################################################################################
################################################################################
library(neuroestimator)
# exp 1
sor.hc.act <- neuroestimator(sor.hc.df, species = "mmusculus") 
# exp 2
tau.ctrl.act <- neuroestimator(tau.ctrl.df, species = "mmusculus") 

## combine and save
ne.res <- inner_join(sor.hc.meta, sor.hc.act %>% rownames_to_column("cell")) %>%
  mutate(experiment = "SOR-HC") %>%
  rbind(inner_join(tau.ctrl.meta, tau.ctrl.act %>% rownames_to_column("cell")) %>%
          mutate(experiment = "Tau-Ctrl"))

write_rds(ne.res, "data/derivatives/ne-results.rds", compress = "gz")
ne.res <- read_rds("data/derivatives/ne-results.rds")
################################################################################
################################################################################
################################################################################
## spatial viz
# exp 1
sor.hc <- AddMetaData(object = sor.hc, metadata = sor.hc.act, col.name = "neuroestimator_activity")
DefaultAssay(sor.hc) <- "Xenium"

pdf("figs/ne-spatial-viz-sor-hc.pdf", width = 12, height = 13)
ImageFeaturePlot(sor.hc, features = "neuroestimator_activity", 
                 fov = c("fov",paste0("fov.",c(2:8))),cols = c("gray95", "firebrick3"),
                 size =2, scale = "all",
                 border.size = NA, dark.background = F) 
dev.off()
# exp 2
tau.ctrl <- AddMetaData(object = tau.ctrl, metadata = tau.ctrl.act, col.name = "neuroestimator_activity")
DefaultAssay(tau.ctrl) <- "Xenium"

pdf("figs/ne-spatial-viz-tau-hc.pdf", width = 12, height = 13)
ImageFeaturePlot(tau.ctrl, features = "neuroestimator_activity", 
                 fov = c("fov",paste0("fov.",c(2:7))),cols = c("gray95", "firebrick3"),
                 size =1.5, scale = "all",
                 border.size = NA, dark.background = F) 
dev.off()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
ne.res %>% group_by(condition, class, experiment) %>%
  summarise(c = n()) %>% View()

## violin plots
##    by class
p1 <- ne.res %>%
  filter(experiment == "SOR-HC") %>%
  mutate(class = factor(class, levels = names(class.colors))) %>%
  ggplot(aes(x = condition, y = predicted_activity)) +
  geom_violin(aes(fill = class), show.legend = F) + 
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.3) +
  ggpubr::stat_compare_means(color = "red", label.y.npc = 1, size = 3) +
  ggh4x::facet_grid2(cols = vars(class), scales = "free") +
  bw.theme + scale_fill_manual(values = class.colors) +
  labs(x="", y="NEUROeSTIMator predicted activity")
p2 <- ne.res %>%
  filter(experiment != "SOR-HC") %>%
  mutate(class = factor(class, levels = names(class.colors))) %>%
  ggplot(aes(x = condition, y = predicted_activity)) +
  geom_violin(aes(fill = class), show.legend = F) + 
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.3) +
  ggpubr::stat_compare_means(color = "red", label.y.npc = 1,size=3) +
  ggh4x::facet_grid2(cols = vars(class), scales = "free") +
  bw.theme + scale_fill_manual(values = class.colors) +
  labs(x="", y="NEUROeSTIMator predicted activity")
pdf("figs/ne-activity-violin-plots.pdf", width = 14, height = 9)
wrap_plots(p1,p2, ncol = 1)
dev.off()

##    by subclass
p12 <- ne.res %>%
  filter(experiment == "SOR-HC") %>%
  mutate(celltype = factor(celltype, levels = names(subclass.colors))) %>%
  ggplot(aes(x = condition, y = predicted_activity)) +
  geom_violin(aes(fill = celltype), show.legend = F) + 
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.3) +
  ggpubr::stat_compare_means(color = "red", label.y.npc = 1, size = 3) +
  ggh4x::facet_grid2(cols = vars(celltype), scales = "free") +
  bw.theme + scale_fill_manual(values = subclass.colors) +
  labs(x="", y="NEUROeSTIMator predicted activity")
p22 <- ne.res %>%
  filter(experiment != "SOR-HC") %>%
  mutate(celltype = factor(celltype, levels = names(subclass.colors))) %>%
  ggplot(aes(x = condition, y = predicted_activity)) +
  geom_violin(aes(fill = celltype), show.legend = F) + 
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.3) +
  ggpubr::stat_compare_means(color = "red", label.y.npc = 1,size=3) +
  ggh4x::facet_grid2(cols = vars(celltype), scales = "free") +
  bw.theme + scale_fill_manual(values = subclass.colors) +
  labs(x="", y="NEUROeSTIMator predicted activity")
pdf("figs/ne-activity-violin-plots-subclass.pdf", width = 28, height = 9)
wrap_plots(p12,p22, ncol = 1)
dev.off()
################################################################################
################################################################################
################################################################################
exp <- unique(ne.res$experiment)
classes <- unique(ne.res$class)
subclasses <- unique(ne.res$celltype)
################################################################################
################################################################################
################################################################################
## for each celltype/class per experiment, build a model to predict condition
##    and get effect size
registerDoMC(4)
ne.LR.res <- foreach(ee = 1:length(exp), .combine = rbind) %dopar% {
  exp.n <- exp[ee]
  class.res <- foreach(cc = 1:length(classes), .combine = rbind) %dopar% {
    class.n <- classes[cc]
    
    df <- ne.res %>% filter(experiment == exp.n, class == class.n) %>%
      mutate(condition_binary = ifelse(condition %in% c("SOR","Tau"), 1,0),
             predicted_activity = scale(predicted_activity, T,T)[,1])
    coefs_table(glm(condition_binary ~ predicted_activity, 
                    data = df,family = binomial)) %>%
      mutate(experiment = exp.n, source = "class", name = class.n)
  }
  subclass.res <- foreach(sc = 1:length(subclasses), .combine = rbind) %dopar% {
    subclass.n <- subclasses[sc]
    
    df <- ne.res %>% filter(experiment == exp.n, celltype == subclass.n) %>%
      mutate(condition_binary = ifelse(condition %in% c("SOR","Tau"), 1,0),
             predicted_activity = scale(predicted_activity, T,T)[,1])
    coefs_table(glm(condition_binary ~ predicted_activity, 
                    data = df,family = binomial)) %>%
      mutate(experiment = exp.n, source = "subclass", name = subclass.n)
  }
  rbind(class.res, subclass.res)
}
p50 <- ne.LR.res %>%
  filter(x=="predicted_activity", source =="class") %>%
  mutate(sig = case_when(pval < 0.001 ~ "***",
                         pval < 0.01 ~ "**",
                         pval < 0.05 ~ "*")) %>%
  ggplot(aes(x = Estimate, y = factor(name,levels = names(class.colors)[8:1]), color = name,
             alpha = !is.na(sig))) +
  geom_vline(xintercept = 0, color = "red",linetype=2) +
  geom_point(show.legend = F) + 
  geom_text(aes(label = sig), vjust = 0, show.legend = F) +
  geom_errorbarh(aes(xmin = confin_min, xmax = confin_max), height = 0.2, 
                 show.legend = F) +
  scale_alpha_manual(values = c(0.3,1))+
  scale_color_manual(values = class.colors) +
  facet_wrap(~experiment, scales = "free") +
  bw.theme + labs(y = "", caption=paste0("*     pval < 0.05\n",
                                         "**    pval < 0.01\n",
                                         "***   pval < 0.001\n"))
p51 <- ne.LR.res %>%
  filter(x=="predicted_activity", source =="subclass") %>%
  mutate(sig = case_when(pval < 0.001 ~ "***",
                         pval < 0.01 ~ "**",
                         pval < 0.05 ~ "*")) %>%
  ggplot(aes(x = Estimate, y = factor(name,levels = names(subclass.colors)[18:1]), color = name,
             alpha = !is.na(sig))) +
  geom_vline(xintercept = 0, color = "red",linetype=2) +
  geom_point(show.legend = F) + 
  geom_text(aes(label = sig), vjust = 0, show.legend = F) +
  geom_errorbarh(aes(xmin = confin_min, xmax = confin_max), height = 0.2, 
                 show.legend = F) +
  scale_alpha_manual(values = c(0.3,1))+
  scale_color_manual(values = subclass.colors) +
  facet_wrap(~experiment, scales = "free") +
  bw.theme + labs(y = "", caption=paste0("*     pval < 0.05\n",
                                         "**    pval < 0.01\n",
                                         "***   pval < 0.001\n"))
pdf("figs/ne-activity-forest-plots.pdf", width = 8, height =14)
wrap_plots(p50,p51, ncol = 1,heights = c(1,2))
dev.off()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## there's a slight imbalance between count of cells per group/cluster
registerDoMC(cores = 3)
ne.bs.res <- foreach(cc = 1:length(classes), .combine = rbind) %dopar% {
  class.n <- classes[cc]
  exp.res <- foreach(ee = 1:length(exp), .combine = rbind) %dopar% {
    exp.n <- exp[ee]
    class.exp.df <- ne.res %>%
      filter(experiment == exp.n, class == class.n) %>%
      drop_na()
    sam.n <- min(table(class.exp.df$condition) %>% as.numeric())
    ii.bs <- foreach(ii = 1:100, .combine = rbind) %dopar% {
      class.exp.df.i <- class.exp.df %>%
        group_by(condition) %>%
        sample_n(size = sam.n, replace = T) %>%
        ungroup()
      conditions <- unique(class.exp.df.i$condition)
      null.d <- class.exp.df.i %>%
        mutate(condition = sample(rep(conditions, each=nrow(class.exp.df.i)/2)))
      
      data.frame(experiment = exp.n, class = class.n,
                 iteration = ii,condition=conditions,
                 mean_activity = c(mean(class.exp.df.i$predicted_activity[class.exp.df.i$condition==conditions[1]]),
                                   mean(class.exp.df.i$predicted_activity[class.exp.df.i$condition==conditions[2]])),
                 null_activity = c(mean(null.d$predicted_activity[null.d$condition==conditions[1]]),
                                   mean(null.d$predicted_activity[null.d$condition==conditions[2]])))
    }
    ii.bs %>%
      group_by(experiment, class, condition) %>%
      dplyr::summarise(confin_min = quantile(mean_activity, 0.025),
                       confin_max = quantile(mean_activity, 0.975),
                       mean_activity = mean(mean_activity),
                       null_confin_min = quantile(null_activity, 0.025),
                       null_confin_max = quantile(null_activity, 0.975),
                       null_activity = mean(null_activity))
  }
  return(exp.res)
}
write_rds(ne.bs.res, "data/derivatives/ne-results-bootstrapped.rds",compress = "gz")
ne.bs.res <- read_rds("data/derivatives/ne-results-bootstrapped.rds")


## plot
p3 <- ne.bs.res %>%
  filter(experiment == "SOR-HC") %>%
  pivot_longer(cols = -c(1:3)) %>%
  mutate(group = paste(condition, name, sep = "__")) %>%
  pivot_wider(names_from = group, values_from = value, id_cols = c(experiment , class)) %>%
  ggplot(aes(x = HC__mean_activity, y = SOR__mean_activity, color = class)) +
  geom_point() + geom_abline(slope = 1, linetype = 1, color = "red") +
  geom_errorbarh(aes(xmin = HC__confin_min, xmax = HC__confin_max)) +
  geom_errorbar(aes(ymin = SOR__confin_min, ymax = SOR__confin_max)) +
  xlim(0.83,1) + ylim(0.83,1) + scale_color_manual(values = cell.colors) +
  bw.theme + labs(x = "HC predicted activity", y = "SOR predicted activity")

p4 <- ne.bs.res %>%
  filter(experiment != "SOR-HC") %>%
  pivot_longer(cols = -c(1:3)) %>%
  mutate(group = paste(condition, name, sep = "__")) %>%
  pivot_wider(names_from = group, values_from = value, id_cols = c(experiment , class)) %>%
  ggplot(aes(x = Tau__mean_activity, y = Control__mean_activity, color = class)) +
  geom_point() + geom_abline(slope = 1, linetype = 1, color = "red") +
  geom_errorbarh(aes(xmin = Tau__confin_min, xmax = Tau__confin_max)) +
  geom_errorbar(aes(ymin = Control__confin_min, ymax = Control__confin_max)) +
  xlim(0.87,1) + ylim(0.87,1) +scale_color_manual(values = cell.colors) +
  bw.theme + labs(x = "Tau predicted activity", y = "Control predicted activity")

p5 <- ne.bs.res %>%
  filter(condition %in% c("SOR", "Control")) %>%
  pivot_longer(cols = -c(1:3)) %>%
  mutate(group = paste(condition, name, sep = "__")) %>%
  pivot_wider(names_from = group, values_from = value, id_cols = c(class)) %>%
  ggplot(aes(x = SOR__mean_activity, y = Control__mean_activity, color = class)) +
  geom_point() + geom_abline(slope = 1, linetype = 1, color = "red") +
  geom_errorbarh(aes(xmin = SOR__confin_min, xmax = SOR__confin_max)) +
  geom_errorbar(aes(ymin = Control__confin_min, ymax = Control__confin_max)) +
  xlim(0.87,1) + ylim(0.87,1) +scale_color_manual(values = cell.colors) +
  bw.theme + labs(x = "SOR predicted activity", y = "Control predicted activity")

wrap_plots(p3,p4,p5, nrow = 1)
ggsave2("figs/ne-activity-point-error-plots.png", width = 12, height = 5)
################################################################################
################################################################################
################################################################################
################################################################################
## show how many genes from the data are the main list of NE
table(ne.genes %in% rownames(sor.hc.df))
pdf("figs/neuroestimator-yes-no-genes.pdf", width = 3, height = 8)
data.frame(NE_gene = ne.genes) %>%
  mutate(present = NE_gene %in% rownames(sor.hc.df)) %>%
  ggplot(aes(x = "", y = NE_gene, fill = present)) +
  geom_tile() + scale_fill_manual(values = abstract.colors[c(4,2)]) +
  bw.theme + labs(x="", y = "NEUROeSTIMator genes")
dev.off()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
