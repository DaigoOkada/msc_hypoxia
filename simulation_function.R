#Simulation function
Simulation <- function(initial_cell_ages, n_generations = 100, age_increase=0.025, gate_generation = NULL,
                       div_t0=0.3, div_k=10, death_t0=0.5, death_k=10, seed = 123) {

    #Define function
    set.seed(seed)
    div_prob_cal <- function(cell_ages) {
        div_prob <- 1 / (1 + exp(div_k * (cell_ages - div_t0)))
        return(div_prob)
    }
    death_prob_cal <- function(cell_ages) {
        death_prob <- 1 / (1 + exp(-death_k * (cell_ages - death_t0)))
        return(death_prob)
    }

    #Parameter setting
    n_cells <- length(initial_cell_ages)
    p0_cell_ages <- initial_cell_ages
    cell_ages <- initial_cell_ages
    div_probs <- div_prob_cal(cell_ages)
    death_prob <- death_prob_cal(cell_ages)
    cell_age_gen <- list()
    cell_age_gen[[1]] <- p0_cell_ages
    cell_inc_rates <- c(1)

    #Siulation
    for(i in 1:n_generations){
    

        cur_cell_num <- length(cell_ages)
        div_probs <- div_prob_cal(cell_ages)
        div_hanteis <- sapply(div_probs, function(prob) rbinom(1, size=1, prob=prob))
        musume_cell_ages <- cell_ages[div_hanteis == 1] + abs(rnorm(sum(div_hanteis==1), mean=age_increase, sd=age_increase/10))
        cell_ages[div_hanteis == 1] <- cell_ages[div_hanteis == 1] + abs(rnorm(sum(div_hanteis==1), mean=age_increase, sd=age_increase/10)) #分裂するごとに親細胞の年齢が増加
        cell_ages <- c(cell_ages, musume_cell_ages)

        #Cell death
        death_probs <- death_prob_cal(cell_ages)
        death_hanteis <- sapply(death_probs, function(prob) rbinom(1, size=1, prob=prob))
        cell_ages <- cell_ages[death_hanteis == 0]
        r <-  length(cell_ages)/cur_cell_num
        if(r < 0.9){
            print(paste("Generation", i, ":", length(cell_ages), "cells, not enough cells to sample."))
            break
        }

        if(is.null(gate_generation)) {
            cell_ages <- sample(cell_ages, min(length(cell_ages), n_cells), replace = FALSE)  
        }else if(i %in% gate_generation) {
            cell_ages <- sort(cell_ages)[1:min(length(cell_ages), n_cells)]
        }else{
            cell_ages <- sample(cell_ages, min(length(cell_ages), n_cells), replace = FALSE)
        }
        cell_inc_rates <- c(cell_inc_rates, r)
        cell_age_gen[[i+1]] <- cell_ages
    }
    names(cell_age_gen) <- 0:(length(cell_age_gen) - 1)
    names(cell_inc_rates) <- 0:(length(cell_inc_rates) - 1)
    res <- list(cell_age_gen = cell_age_gen, cell_inc_rates = cell_inc_rates)
    return(res)
}