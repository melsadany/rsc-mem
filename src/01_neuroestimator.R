################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
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
cell.colors <- c("#e76cf2", "#04aef3", "#04bd7c", "#a3a700", "#fa766d", "#cab2d6", "#b4df8b", "#a6cfe2")
names(cell.colors) <- c("ExN", "InN", "Astro", "Oligo", "Microglia", "Endo","VLMC","OPC")
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
ne.res %>% group_by(condition, class, experiment) %>%
  summarise(c = n()) %>% View()

p1 <- ne.res %>%
  filter(experiment == "SOR-HC") %>%
  ggplot(aes(x = condition, y = predicted_activity)) +
  geom_violin(aes(fill = class), show.legend = F) + geom_boxplot(width = 0.2, fill = "white") +
  ggpubr::stat_compare_means(color = "red", label.y.npc = 0.99, size = 3) +
  ggh4x::facet_grid2(cols = vars(class), scales = "free") +
  bw.theme + scale_fill_manual(values = cell.colors) +
  labs(x="", y="NEUROeSTIMator predicted activity")
p2 <- ne.res %>%
  filter(experiment != "SOR-HC") %>%
  ggplot(aes(x = condition, y = predicted_activity)) +
  geom_violin(aes(fill = class), show.legend = F) + geom_boxplot(width = 0.2, fill = "white") +
  ggpubr::stat_compare_means(color = "red", label.y.npc = 0.99,size=3) +
  ggh4x::facet_grid2(cols = vars(class), scales = "free") +
  bw.theme + scale_fill_manual(values = cell.colors) +
  labs(x="", y="NEUROeSTIMator predicted activity")

wrap_plots(p1,p2, ncol = 1)
ggsave2("figs/ne-activity-violin-plots.png", width = 14, height = 9)
################################################################################
################################################################################
## there's a slight imbalance between count of cells per group/cluster
classes <- unique(ne.res$class)
exp <- unique(ne.res$experiment)
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
