#MSC DDP simulation
#source("/Users/dokada/Dropbox/analysis/2025.4/msc_normhypo_simu0126.R") #Run2026.1.26
out_path <- "/Users/dokada/Desktop/work/msc_normhypo_simu0126/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#set seed
seed <- 100
set.seed(seed)

#load function
library(vioplot)
library(ggplot2)
source("/Users/dokada/Dropbox/analysis/2025.4/pseudo_msc_function.R") #Simulation関数を読み込み

#Set initial cell ages from young and old donor
n_generations <- 100
n_cells <- 10000
n_bad_cells <- rbinom(1, size=n_cells, p=0.05)
young_donor  <- c(rbeta(n_cells - n_bad_cells, shape1 = 3, shape2 = 15), rbeta(n_bad_cells, shape1 =15, shape2 = 3))
n_bad_cells <- rbinom(1, size=n_cells, p=0.05)
old_donor  <- c(rbeta(n_cells - n_bad_cells, shape1 = 1, shape2 = 1), rbeta(n_bad_cells, shape1 =15, shape2 = 3))#Visualize initial distributions

#Parameter setting for hypoxia analysis
hypo_div_t0 <- 0.5
hypo_death_t0 <- 0.75
hypo_age_increase <- 0.02


#Youngでの通常培養継続 vs低酸素培養切り替え
all_norm <- Simulation(young_donor, n_generations = n_generations)
all_mean_age <- sapply(all_norm$cell_age_gen, mean)
all_cell_inc_rates <- all_norm$cell_inc_rates
change_passage <- 6 #passage0 - passage5
est_cells_age <- all_norm$cell_age_gen[[change_passage]]
post_est_old <- Simulation(est_cells_age, n_generations = n_generations - change_passage + 1, age_increase=hypo_age_increase, div_t0=hypo_div_t0, death_t0=hypo_death_t0)
post_age <- sapply(post_est_old$cell_age_gen, mean)
kirikae_mean_age <- c(all_mean_age[1:change_passage], post_age[-1])
kirikae_cell_inc_rates <- c(all_norm$cell_inc_rates[1:change_passage], post_est_old$cell_inc_rates[-1])
mat_age <- matrix(NA, n_generations + 1, 2)
mat_age[1:length(all_mean_age),1] <- all_mean_age
mat_age[1:length(kirikae_mean_age),2] <- kirikae_mean_age
idx <- which(apply(mat_age, 1, function(v) any(!is.na(v))))
stopifnot(all(diff(idx) == 1))
png(paste0(out_path, "young_allnorm_switichhypo_age.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
xlab <- "Passage"
ylab <- "Mean Cell Age"
main <- ""
cols <- c("black", "red")
x_passage <- 0:(length(idx)-1)
matplot(x_passage, mat_age[idx,], type = "l", lwd = 4, lty = 1, pch = 19, col = cols, xlab = "", ylab = "", cex.axis = 4, cex.lab = 4, axes = TRUE)
abline(v = 5,  lty = "dotted", lwd = 4)
mtext(xlab, side = 1, line = 6, cex = 4)
mtext(ylab, side = 2, line = 6, cex = 4)
mtext(main, side = 3, line = 3, cex = 4)
dev.off()
#Cell counts
cell_inc_rates <- all_cell_inc_rates
all_cell_counts <- c()
cur_cell_counts <- n_cells
for(i in 1:length(cell_inc_rates)){
    cur_cell_counts <- cur_cell_counts * cell_inc_rates[i]
    all_cell_counts <- c(all_cell_counts, cur_cell_counts)
}
cell_inc_rates <- kirikae_cell_inc_rates
kirikae_cell_counts <- c()
cur_cell_counts <- n_cells
for(i in 1:length(cell_inc_rates)){
    cur_cell_counts <- cur_cell_counts * cell_inc_rates[i]
    kirikae_cell_counts <- c(kirikae_cell_counts, cur_cell_counts)
}
mat_count <- matrix(NA, n_generations, 2)
mat_count[1:length(all_cell_counts),1] <- all_cell_counts
mat_count[1:length(kirikae_cell_counts),2] <- kirikae_cell_counts
idx <- which(apply(mat_count, 1, function(v) any(!is.na(v))))
stopifnot(all(diff(idx) == 1))
png(paste0(out_path, "old_nochange_count.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
xlab <- "Passage"
ylab <- "log10(Cell count)"
main <- ""
cols <- c("black", "red")
x_passage <- 0:(length(idx)-1)
matplot(x_passage, log10(mat_count[idx,]), type = "l", lwd = 4, lty = 1, pch = 19, col = cols, xlab = "", ylab = "", cex.axis = 4, cex.lab = 4, axes = TRUE)
abline(v = 5,  lty = "dotted", lwd = 4)
mtext(xlab, side = 1, line = 6, cex = 4)
mtext(ylab, side = 2, line = 6, cex = 4)
mtext(main, side = 3, line = 3, cex = 4)
dev.off()


#保存するファイル名を指定 
sink(paste0(out_path, "session_info.txt")) 
print(sessionInfo()) 
sink()