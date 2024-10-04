J # vector
ps # vector
dt_val # vector
ks # vector

dx # scalar
alpha # scalar
mu # scalar

update_bact_pop <- function(bact_pop, pop_size, dt_val) {
  L <- colSums(J * (ps * bact_pop)) * dx / (1 + pop_size)^alpha
  (bact_pop + dt_val * L) / (1 + dt_val * (mu + ks))
}

bact_pop <- matrix(nrow = length(mr_bact_pop), ncol = length(dt) + 1)
bact_pop[, 1] <- mr_bact_pop
mr_pop_size <- pop_size <- sum(mr_bact_pop)
for (i in seq_along(dt)) {    
  mr_bact_pop <- update_bact_pop(mr_bact_pop, mr_pop_size, dt[i])
  bact_pop[, i + 1] <- mr_bact_pop
  mr_pop_size <- sum(mr_bact_pop)
  pop_size <- c(pop_size, mr_pop_size)
}


# output: bact_pop and pop_size

#######################################################################################


