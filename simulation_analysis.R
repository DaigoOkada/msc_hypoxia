out_path <- "/Users/dokada/Desktop/work/pseudo_msc_simu_int0817/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#set seed
seed <- 1000
set.seed(seed)

#load function
library(vioplot)
library(ggplot2)
source("/Users/dokada/Dropbox/analysis/2025.4/pseudo_msc_function.R")

#Set initial cell ages from young and old donor
n_cells <- 10000
n_bad_cells <- rbinom(1, size=n_cells, p=0.05)
young_donor  <- c(rbeta(n_cells - n_bad_cells, shape1 = 3, shape2 = 15), rbeta(n_bad_cells, shape1 =15, shape2 = 3))
n_bad_cells <- rbinom(1, size=n_cells, p=0.05)
old_donor  <- c(rbeta(n_cells - n_bad_cells, shape1 = 1, shape2 = 1), rbeta(n_bad_cells, shape1 =15, shape2 = 3))#Visualize initial distributions
png(paste0(out_path, "young_init_hist.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4), mgp = c(2.5, 2.5, 0))
hist(young_donor, breaks=50, main="", xlab="", ylab="", cex.axis=4, cex.lab=4, cex.main=4, xlim=c(0,1))
xlab =  "Cell Age"
ylab = "Frequency"
main = "Young"
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4, adj=0)
dev.off()
png(paste0(out_path, "old_init_hist.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4), mgp = c(2.5, 2.5, 0))
hist(old_donor, breaks=50, main="", xlab="", ylab="", cex.axis=4, cex.lab=4, cex.main=4, xlim=c(0,1))
xlab =  "Cell Age"
ylab = "Frequency"
main = "Old"
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4, adj=0)
dev.off()

#Young example
set.seed(1000)
control <- Simulation(young_donor, n_generations = 100, seed=1000)
donor_age <- "young"
mean_cell_age <- sapply(control$cell_age_gen, mean)
feature <- "age"
png(paste0(out_path, donor_age, "_", feature, ".example.png"), width=960, height=960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0)) 
xlab = "Passage"
ylab = "Mean Cell Age"
main = ""
plot(as.integer(names(mean_cell_age)), mean_cell_age, type="l", lwd=4, cex = 4, xlab="", ylab = "", cex.axis=4, cex.lab=4)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4)
dev.off()
#Cell counts
cell_inc_rates <- control$cell_inc_rates
cell_counts <- c()
cur_cell_counts <- n_cells
for(i in 1:length(cell_inc_rates)){
    cur_cell_counts <- cur_cell_counts * cell_inc_rates[i]
    cell_counts <- c(cell_counts, cur_cell_counts)
}
feature <- "count"
png(paste0(out_path, donor_age, "_", feature, ".example.png"), width=960, height=960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0)) 
xlab = "Passage"
ylab = "log(cell count)"
main = ""
plot(as.integer(names(cell_inc_rates)), log10(cell_counts), type="l", lwd=4, cex = 4, xlab="", ylab = "", cex.axis=4, cex.lab=4)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4)
dev.off()
#Passage vs Cell Age
pt_list <- control$cell_age_gen[1:10]
passage_levels <- names(pt_list)
png(paste0(out_path, donor_age, "_cellage_passage.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
vioplot(pt_list, names = passage_levels, col = "lightgray", border = "black", xlab = "", ylab = "", cex.axis=4)
xlab = "Passage"
ylab = "Cell Age"
main = ""
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4)
for (i in seq_along(pt_list)) {
    stripchart(pt_list[[i]], at = i, method = "jitter", pch = 16, col = rgb(0,0,0,0.5), vertical = TRUE, add = TRUE, cex=1.9)
}
dev.off()


#div, death function
div_prob_cal <- function(cell_ages, div_t0, div_k=10) {
    div_prob <- 1 / (1 + exp(div_k * (cell_ages - div_t0)))
    return(div_prob)
}
death_prob_cal <- function(cell_ages, death_t0, death_k=10) {
    death_prob <- 1 / (1 + exp(-death_k * (cell_ages - death_t0)))
    return(death_prob)
}
div_t0_can <- seq(0.1, 0.9, by=0.1)
death_t0_can <- seq(0.1, 0.9, by=0.1)
x_grid <- seq(0, 1, length.out=1000)
func_can <- c("Div", "Death")
for(i in 1:length(func_can)){
    if(func_can[i] == "Div"){
        mat <- matrix(NA, length(x_grid), length(div_t0_can))
        for(j in 1:ncol(mat)){
            mat[,j] <- div_prob_cal(cell_ages=x_grid, div_t0=div_t0_can[j])
        }
        groups <- div_t0_can
        legend_place <- "topright"
    }else{
        mat <- matrix(NA, length(x_grid), length(death_t0_can))
        for(j in 1:ncol(mat)){
            mat[,j] <- death_prob_cal(cell_ages=x_grid, death_t0=death_t0_can[j])
        }
        groups <- death_t0_can
        legend_place <- "topleft"
    }
    png(paste0(out_path, func_can[i], "_func_lowoxygen.png"), width = 960, height = 960)
    par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
    xlab <- "Cell Age"
    ylab <- "Probability"
    main <- func_can[i]
    cols <- c("black", "red", "blue", "orange", "purple", "gold", "brown", "deepskyblue3", "gray40")
    matplot(x_grid, mat, type = "l", lwd = 4, lty = 1, pch = 19, col = cols, xlab = "", ylab = "", cex.axis = 4, cex.lab = 4, axes = TRUE, ylim=c(0,1))
    mtext(xlab, side = 1, line = 6, cex = 4)
    mtext(ylab, side = 2, line = 6, cex = 4)
    mtext(main, side = 3, line = 3, cex = 4)
    legend(legend_place, legend = groups, col = cols, lty = 1, lwd = 4, cex = 1.9)
    dev.off()
}


#set parameter candidate
div_t0_can <- seq(0.1, 0.9, by=0.1)
death_t0_can <- seq(0.1, 0.9, by=0.1)
age_increase_can <- seq(0.01, 0.09, by=0.01)/2


#div_t0, death_t0 simulation
n_generations <- 100
age_increase <- 0.025 #fixed
para_can <- expand.grid(div_t0_can, death_t0_can)
colnames(para_can) <- c("div_t0", "death_t0")
result <- matrix(NA, nrow(para_can), 4)
colnames(result) <- c("Lifespan_Young", "Age_Young", "Lifespan_Old", "Age_Old")
for(i in 1:nrow(para_can)){
    tmp_div_t0 <- para_can[i,1]
    tmp_death_t0 <- para_can[i,2]
    control <- Simulation(young_donor, n_generations = n_generations, div_t0=tmp_div_t0, death_t0= tmp_death_t0, seed=seed)
    tmp <- sapply(control$cell_age_gen, mean)
    lifespan_y <- length(tmp) - 1
    quality_y <-  ifelse(lifespan_y > 0, tmp[ceiling(lifespan_y*0.5) + 1], NA)
    control <- Simulation(old_donor, n_generations = n_generations, div_t0=tmp_div_t0, death_t0=tmp_death_t0, seed=seed)
    tmp <- sapply(control$cell_age_gen, mean)
    lifespan_o <- length(tmp) - 1
    quality_o <- ifelse(lifespan_o > 0, tmp[ceiling(lifespan_o*0.5) + 1], NA)
    result[i,] <- c(lifespan_y, quality_y, lifespan_o, quality_o)
    cat(i, "\n")
}

#heatmap
for(i in 1:ncol(result)){
    para_can2 <- cbind(para_can, result[,i])
    if(i %in% c(1, 3)){
        color_limits <- c(0, n_generations + 1)
    }else{
        color_limits <- c(0, 1.01)
    }
    colnames(para_can2)[3] <- "score"
    p <- ggplot(para_can2, aes(x = factor(div_t0), y = factor(death_t0), fill = score)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "plasma", limits=color_limits) +
    labs(x = "div_t0", y = "death_t0", fill = "Score", title = colnames(result)[i]) +
    theme_minimal(base_size = 16) + 
    theme(
        axis.text.x = element_text(size = 25, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 30),
        plot.margin = margin(30, 30, 30, 30),
        legend.text   = element_text(size = 18),
        legend.title  = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"),
        plot.title     = element_text(size = 28, face = "bold", hjust = 0.5) 
    ) +
    coord_fixed()
    ggsave(filename = paste0(out_path, colnames(result)[i], ".div_death.png"), plot = p,width = 8, height = 8, units = "in", dpi = 300)
}


#age_increase x death_t0
para_can <- expand.grid(age_increase_can, death_t0_can)
colnames(para_can) <- c("age_increase", "death_t0")
result <- matrix(NA, nrow(para_can), 4)
colnames(result) <- c("Lifespan_Young", "Age_Young", "Lifespan_Old", "Age_Old")
for(i in 1:nrow(para_can)){
    tmp_age_increase <- para_can[i,1]
    tmp_death_t0 <- para_can[i,2]
    control <- Simulation(young_donor, n_generations = n_generations, age_increase= tmp_age_increase, death_t0= tmp_death_t0, seed=seed)
    tmp <- sapply(control$cell_age_gen, mean)
    lifespan_y <- length(tmp) - 1
    quality_y <-  ifelse(lifespan_y > 0, tmp[ceiling(lifespan_y*0.5) + 1], NA)
    control <- Simulation(old_donor, n_generations = n_generations, age_increase=tmp_age_increase, death_t0=tmp_death_t0, seed=seed)
    tmp <- sapply(control$cell_age_gen, mean)
    lifespan_o <- length(tmp) - 1
    quality_o <- ifelse(lifespan_o > 0, tmp[ceiling(lifespan_o*0.5) + 1], NA)
    result[i,] <- c(lifespan_y, quality_y, lifespan_o, quality_o)
    cat(i, "\n")
}


#heatmap
for(i in 1:ncol(result)){
    para_can2 <- cbind(para_can, result[,i])
    if(i %in% c(1, 3)){
        color_limits <- c(0, n_generations + 1)
    }else{
        color_limits <- c(0, 1.01)
    }
    colnames(para_can2)[3] <- "score"
    p <- ggplot(para_can2, aes(x = factor(age_increase), y = factor(death_t0), fill = score)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "plasma", limits=color_limits) +
    labs(x = " age_increase", y = "death_t0", fill = "Score", title = colnames(result)[i]) +
    theme_minimal(base_size = 16) + 
    theme(
        axis.text.x = element_text(size = 25, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 30),
        plot.margin = margin(30, 30, 30, 30),
        legend.text   = element_text(size = 18),
        legend.title  = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"),
        plot.title     = element_text(size = 28, face = "bold", hjust = 0.5) 
    ) +
    coord_fixed()
    ggsave(filename = paste0(out_path, colnames(result)[i], ".age_death.png"), plot = p,width = 8, height = 8, units = "in", dpi = 300)
}


#Parameter setting for hypoxia analysis
hypo_div_t0 <- 0.5
hypo_death_t0 <- 0.75
hypo_age_increase <- 0.02

#
set.seed(1000)
n_iters <- 10
lifespan <- rep(NA, n_iters)
for(i in 1:n_iters){
    tmp <- Simulation(old_donor, n_generations = n_generations)
    lifespan[i] <- length(tmp$cell_inc_rates) - 1
}
stopifnot(all(lifespan==0))


all_hypo <- Simulation(old_donor, n_generations = n_generations, age_increase=hypo_age_increase, div_t0=hypo_div_t0, death_t0=hypo_death_t0)
all_mean_age <- sapply(all_hypo$cell_age_gen, mean)
all_cell_inc_rates <- all_hypo$cell_inc_rates
change_passage <- 6 #passage0 - passage5
est_cells_age <- all_hypo$cell_age_gen[[change_passage]]
post_est_old <- Simulation(est_cells_age, n_generations = n_generations - change_passage + 1)
post_age <- sapply(post_est_old$cell_age_gen, mean)
kirikae_mean_age <- c(all_mean_age[1:change_passage], post_age[-1])
kirikae_cell_inc_rates <- c(all_hypo$cell_inc_rates[1:change_passage], post_est_old$cell_inc_rates[-1])
mat_age <- matrix(NA, n_generations + 1, 2)
mat_age[1:length(all_mean_age),1] <- all_mean_age
mat_age[1:length(kirikae_mean_age),2] <- kirikae_mean_age
idx <- which(apply(mat_age, 1, function(v) any(!is.na(v))))
stopifnot(all(diff(idx) == 1))
png(paste0(out_path, "old_nochange_age.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
xlab <- "Passage"
ylab <- "Mean Cell Age"
main <- ""
cols <- c("black", "red")
x_passage <- 0:(length(idx)-1)
matplot(x_passage, mat_age[idx,], type = "l", lwd = 4, lty = 1, pch = 19, col = cols, xlab = "", ylab = "", cex.axis = 4, cex.lab = 4, axes = TRUE)
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
mtext(xlab, side = 1, line = 6, cex = 4)
mtext(ylab, side = 2, line = 6, cex = 4)
mtext(main, side = 3, line = 3, cex = 4)
dev.off()



young_res <- Simulation(young_donor, n_generations = n_generations)
young_mean_age <- sapply(young_res$cell_age_gen, mean)
young_cell_inc_rates <- young_res$cell_inc_rates
mat_age <- matrix(NA, n_generations + 1, 2)
mat_age[1:length(young_mean_age),1] <- young_mean_age
mat_age[1:length(kirikae_mean_age),2] <- kirikae_mean_age
idx <- which(apply(mat_age, 1, function(v) any(!is.na(v))))
stopifnot(all(diff(idx) == 1))
png(paste0(out_path, "youngnormal_oldhypo_age.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
xlab <- "Passage"
ylab <- "Mean Cell Age"
main <- ""
cols <- c("black", "red")
x_passage <- 0:(length(idx)-1)
matplot(x_passage, mat_age[idx,], type = "l", lwd = 4, lty = 1, pch = 19, col = cols, xlab = "", ylab = "", cex.axis = 4, cex.lab = 4, axes = TRUE)
mtext(xlab, side = 1, line = 6, cex = 4)
mtext(ylab, side = 2, line = 6, cex = 4)
mtext(main, side = 3, line = 3, cex = 4)
dev.off()
#Cell counts
cell_inc_rates <- young_cell_inc_rates
young_cell_counts <- c()
cur_cell_counts <- n_cells
for(i in 1:length(cell_inc_rates)){
    cur_cell_counts <- cur_cell_counts * cell_inc_rates[i]
    young_cell_counts <- c(young_cell_counts, cur_cell_counts)
}
mat_count <- matrix(NA, n_generations, 2)
mat_count[1:length(young_cell_counts),1] <- young_cell_counts
mat_count[1:length(kirikae_cell_counts),2] <- kirikae_cell_counts
idx <- which(apply(mat_count, 1, function(v) any(!is.na(v))))
stopifnot(all(diff(idx) == 1))
png(paste0(out_path, "youngnormal_oldhypo_count.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
xlab <- "Passage"
ylab <- "log10(Cell count)"
main <- ""
cols <- c("black", "red")
x_passage <- 0:(length(idx)-1)
matplot(x_passage, log10(mat_count[idx,]), type = "l", lwd = 4, lty = 1, pch = 19, col = cols, xlab = "", ylab = "", cex.axis = 4, cex.lab = 4, axes = TRUE)
mtext(xlab, side = 1, line = 6, cex = 4)
mtext(ylab, side = 2, line = 6, cex = 4)
mtext(main, side = 3, line = 3, cex = 4)
#legend("topright", legend = change_gen_can, col = cols, lty = 1, lwd = 4, cex = 2)
dev.off()

#session info
sink(paste0(out_path, "session_info.txt")) 
print(sessionInfo()) 
sink()