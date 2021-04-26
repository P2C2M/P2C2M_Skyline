#' Run P2C2M.Skyline
#' @export
#' @description Assessing model adequacy for Bayesian Skyline Plots using posterior predictive simulation.
#' @param tree.file A ultrametric phylogenetic tree (NEXUS format).
#' @param log.file The log file resulting from a Bayesian Skyline analyzed in Tracer.
#' @param dir Set the working directory.
#' @param nrep The size of the null distribution, Default: 100.
#' @param path.to.ms Path to ms software.
#' @return A p-value that specifies whether the empirical dataset is or not a poor fit to the Skyline model.
#' @author Emanuel M. Fonseca, Drew J. Duckett, Filipe G. Almeida, Megan L. Smith, Maria Tereza C. Thomé, Bryan C. Carstens


P2C2M.Skyline <- function(tree.file,
                          log.file,
                          dir,
                          nrep=NULL,
                          path.to.ms){

  list.of.packages <- c("ape",
                        "phytools",
                        "phyclust",
                        "pegas")

  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

  suppressMessages(library(ape))
  suppressMessages(library(phyclust))
  suppressMessages(library(phytools))
  suppressMessages(library(pegas))

  if(is.null(nrep)){
    nrep <- 100
  }

  emp_tree <- read.nexus(tree.file)

  if(!is.ultrametric(emp_tree)){
    emp_tree <- force.ultrametric(emp_tree, method="nnls")
  }

  Time <- coalescent.intervals(emp_tree)$total.depth

  ntips <- length(emp_tree$tip.label)
  theta <- theta.tree(emp_tree)$theta

  emp_tree$edge.length<- emp_tree$edge.length/max(nodeHeights(emp_tree)[,2])

  ci <- sort(coalescent.intervals(emp_tree)$interval.length)

  calc <- 0
  for (i in 1:length(ci)){
    calc <- c(calc, (calc[i]+ci[i]))

    if (i == length(ci)){
      calc <- 1-calc
    }
  }

  calc <- calc[-1]
  ci_emp<- sum(calc/(rev(1:length(calc))))

  raw_data <- read.table(log.file,h=T, skip = 1)
  raw_data <- raw_data[complete.cases(raw_data),]

  Time <- raw_data$Time[length(raw_data[,1])]

  ci_sim <- data.frame(matrix(NA,nrep,1))
  setwd(dir)

  for (x in 1:nrep){

    if (x == 1){
      print(noquote("Running P2C2M.skyline"))
      progress <- txtProgressBar(min = 0, max = nrep, style = 3)
    }

    Ne <- numeric()

    for (i in 1:100){
      anc <- runif(1,raw_data$Lower[length(raw_data$Lower)],raw_data$Upper[length(raw_data$Upper)])
      curr <- runif(1,raw_data$Lower[1],raw_data$Upper[1])
      Ne <- c(Ne, anc/curr)
    }

    Ne <- mean(Ne)

    system(paste(path.to.ms, ntips,"1 -t", theta, "-eN", Time, Ne, "-T -seeds", sample(1:100000,1), sample(1:100000,1), sample(1:100000,1), "| tail -n +4 | grep -v //>treefile_sim.tre", sep=" "))

    sim_tree <- readLines(file.path(dir,"treefile_sim.tre"))[1]
    sim_tree <- read.tree(text=sim_tree)

    if(!is.ultrametric(sim_tree)){
      sim_tree <- force.ultrametric(sim_tree, method="nnls")
    }

    sim_tree$edge.length<- sim_tree$edge.length/max(nodeHeights(sim_tree)[,2])

    ci <- coalescent.intervals(sim_tree)$interval.length

    calc <- 0
    for (i in 1:length(ci)){
      calc <- c(calc, (calc[i]+ci[i]))

      if (i == length(ci)){
        calc <- 1-calc
      }
    }

    calc <- calc[-1]
    ci_sim[x,] <- sum(calc/(rev(1:length(calc))))


    setTxtProgressBar(progress, x)

    if (x ==nrep){
      close(progress)
    }
  }

  lower <- length(which(ci_sim <= ci_emp))
  upper <- length(which(ci_sim >= ci_emp))
  p.value <- min(lower,upper)*2/nrep


  if (p.value < 0.05){
    print("Your dataset does violate Skyline model")
    print(noquote(paste("p-value:", p.value)))
  } else {
    print("Your dataset does not violate Skyline model")
    print(noquote(paste("p-value:", p.value)))
  }

}
