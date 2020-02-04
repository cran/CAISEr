## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(smoof))
suppressPackageStartupMessages(library(MOEADr))
suppressPackageStartupMessages(library(CAISEr))

### Build function names (instances: UF1 - UF7, dimensions 10 - 40)
fname   <- paste0("UF_", 1:7)
dims    <- c(10:40)
allfuns <- expand.grid(fname, dims, stringsAsFactors = FALSE)

# Assemble instances list
instances <- vector(nrow(allfuns), mode = "list")
for (i in 1:length(instances)){
  instances[[i]]$FUN <- paste0(allfuns[i,1], "_", allfuns[i,2])
}

### Build the functions listed above (so that they can be properly used)
for (i in 1:nrow(allfuns)){
  assign(x = instances[[i]]$FUN,
     value = MOEADr::make_vectorized_smoof(prob.name  = "UF",
                    dimensions = allfuns[i, 2],
                    id = as.numeric(strsplit(allfuns[i, 1], "_")[[1]][2])))
}

## -----------------------------------------------------------------------------
# Prepare algorithm function to be used in run_experiment():
myalgo <- function(type, instance){
  # Input parameters:
  #     - type (variant to use: "original", "original2", "moead.de" or "moead.de2")
  #     - instance (instance to be solved, e.g., instance = instances[[i]])
  # All other parameters are set internally

  ## Extract instance information to build the MOEADr problem format
  fdef  <- unlist(strsplit(instance$FUN, split = "_"))
  uffun <- smoof::makeUFFunction(dimensions = as.numeric(fdef[3]),
                                 id         = as.numeric(fdef[2]))
  fattr    <- attr(uffun, "par.set")
  prob.dim <- fattr$pars$x$len
  
  ## Build MOEADr problem list
  problem <- list(name = instance$FUN,
                  xmin = fattr$pars$x$lower,
                  xmax = fattr$pars$x$upper,
                  m    = attr(uffun, "n.objectives"))

  ## Load presets for the algorithm provided in input 'type' and 
  ## modify whatever is needed for this particular experiment
  de2 <- FALSE
  if (type == "moead.de2"){
    de2  <- TRUE
    type <- "moead.de"
  }
  algo.preset <- MOEADr::preset_moead(type)
  algo.preset$decomp$H <- 99 # <-- set population size
  algo.preset$stopcrit[[1]]$name <- "maxeval" # <-- type of stop criterion
  algo.preset$stopcrit[[1]]$maxeval <- 2000 * prob.dim # stop crit.
  poly.ind <- which(sapply(algo.preset$variation,
                           function(x){x$name == "polymut"}))
  algo.preset$variation[[poly.ind]]$pm <- 1 / prob.dim # <--- pm = 1/d
  if (de2){
    algo.preset$aggfun$name <- "pbi"
    algo.preset$aggfun$theta <- 5
    algo.preset$neighbors$name = "x"
  }
  
  ## Run algorithm on "instance"
  out <- MOEADr::moead(preset = algo.preset, problem = problem,
                       showpars = list(show.iters = "none"))

  ## Read reference data to calculate the IGD
  Yref  <- as.matrix(read.table(paste0("./inst/extdata/pf_data/",
                                       fdef[1], fdef[2], ".dat")))
  IGD = MOEADr::calcIGD(Y = out$Y, Yref = Yref)

  ## Return IGD as field "value" in the output list
  return(list(value = IGD))
}

## -----------------------------------------------------------------------------
# Assemble Algorithm.list. Notice that we need to provide an alias for each 
# method, since both algorithms have the same '$FUN' argument.
algorithms <- list(list(FUN   = "myalgo", 
                        alias = "Original 1", 
                        type  = "original"),
                   list(FUN   = "myalgo", 
                        alias = "Original 2", 
                        type  = "original2"),
                   list(FUN   = "myalgo", 
                        alias = "MOEAD-DE", 
                        type  = "moead.de"),
                   list(FUN   = "myalgo", 
                        alias = "MOEAD-DE2", 
                        type  = "moead.de2"))

## ---- eval=FALSE--------------------------------------------------------------
#  my.results <- run_experiment(instances  = instances,
#                               algorithms = algorithms,
#                               power = 0.8,      # Desired power: 80%
#                               power.target = "mean", # on average,
#                               d = 0.5,          # to detect differences greater
#                                                 # than 0.5 standard deviations
#                               sig.level = 0.05, # at a 95% confidence level.
#                               se.max = 0.05,    # Measurement error: 5%
#                               dif = "perc",     # on the paired percent
#                                                 # differences of means,
#                               method = "param", # calculated using parametric
#                                                 # formula.
#                               comparisons = "all.vs.all", # Compare all algorithms
#                                                           # vs all others,
#                               nstart = 15,      # start with 15 runs/algo/inst
#                               nmax   = 200,     # and do no more than 200 runs/inst.
#                               seed   = 1234,    # PRNG seed (for reproducibility)
#                               #
#                               # NOTICE: Using all but 1 cores. Change if needed
#                               ncpus  = parallel::detectCores() - 1)

## ---- echo=FALSE--------------------------------------------------------------
load("../inst/extdata/vignette_results.RData")

## ---- fig.align="center", fig.width=8, fig.height=8---------------------------
plot(my.results)

## ---- fig.align="center", fig.width=6, fig.height=10--------------------------
suppressPackageStartupMessages(library(car))

algopairs <- paste(my.results$data.summary$Alg1, 
                   my.results$data.summary$Alg2,
                   sep = " - ")

par(mfrow = c(3, 2))
for (i in seq_along(unique(algopairs))){
  tmp <- my.results$data.summary[algopairs == unique(algopairs)[i], ]
  car::qqPlot(tmp$Phi, 
            pch = 16, las = 1, main = unique(algopairs)[i],
            ylab = "observed", xlab = "theoretical quantiles")
}
par(mfrow = c(1, 1))

## ---- fig.align="center", fig.width=6, fig.height=8---------------------------
par(mfrow = c(3, 2))
for (i in seq_along(unique(algopairs))){
  tmp <- my.results$data.summary[algopairs == unique(algopairs)[i], ]
  boot.means <- CAISEr::boot_sdm(tmp$Phi, boot.R = 999)
  hist(boot.means, breaks = 30, main = unique(algopairs)[i], las = 1)
}
par(mfrow = c(1, 1))

## ---- fig.align="center", fig.width=6, fig.height=4---------------------------
df <- cbind(Comparison = algopairs, my.results$data.summary)

suppressPackageStartupMessages(library(ggplot2))
mp <- ggplot(df, aes(x = Comparison, y = Phi, fill = Comparison))
mp + 
  geom_violin(alpha = 0.6, show.legend = FALSE, scale = "width") + 
  geom_boxplot(alpha = 0, show.legend = FALSE, 
               outlier.shape = NA, width = .15) + 
  geom_point(shape = 16, col = "black", fill = "black", alpha = 0.6,
             position = position_jitter(width = .15)) + 
  geom_abline(slope = 0, intercept = 0, col = "red", lty = 2) + 
  ylab("Percent difference in IGD") + xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ---- fig.align="center", fig.width=6, fig.height=8---------------------------
ggplot(df, aes(x = Instance, y = Phi, colour = Comparison,
                     ymin = Phi - SE, ymax = Phi + SE)) + 
  geom_pointrange(show.legend = FALSE) + 
  geom_abline(slope = 0, intercept = 0, col = 1, lty = 2) + 
  facet_grid(Comparison ~ .) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  xlab("")

## ---- fig.align="center", fig.width=6, fig.height=8---------------------------
summary(my.results, test = "wilcoxon")

