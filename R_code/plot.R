#! /bin/Rscript


## -----------------------
## Libraries
## -----------------------
library(ggplot2)
library(plyr)
library(stringr)
library(tidyr)

## -----------------------
## Read in width for plot
## -----------------------
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
    x_lim   <- 5e6
    toPrint <- "all"
}else if(length(args) == 1){
    x_lim   <- as.integer(args[1])
    toPrint <- "all"
}else{
    x_lim   <- as.integer(args[1])
    toPrint <- args[2]
}

## -----------------------
## Read in dataset
## -----------------------

## List all Directories:
files <- list.files("../tests")
files <- files[!str_detect(files, ".sh")]

## Select what to print
if(toPrint != "all"){
    files <- files[tolower(files) == tolower(toPrint)]
}

## -----------------------
## List all subfiels
## -----------------------
all_files <- llply(files,
                  function(t) t <- list.files(paste0("../tests/", t)))
all_files <- llply(all_files, function(t) t <- t[!str_detect(t, ".txt")])

## -----------------------
## Adhoc tunning
## -----------------------

### SLM VARYING CG
all_files[[2]] <- all_files[[2]][c(5,  6, 7, 8)]
all_files[[1]] <- all_files[[1]][1]


## -----------------------
## All data
## -----------------------
all_data <- c()
for(i in 1:length(files)){
    for(j in 1:length(all_files[[i]])){
        if(length(all_files[[i]]) > 0){
            file <- paste("../tests",
                         files[i],
                         all_files[[i]][j],
                         sep = "/")
            aux_data<- read.csv(file, header = FALSE)
            names(aux_data) <- c("data_points", "precision")
            aux_data$params <- rep(str_replace(all_files[[i]][j], ".csv", ""),
                                  nrow(aux_data))
            all_data <- rbind(all_data, aux_data)
        }
    }
}

## -----------------------
## Plot
## -----------------------
ggplot(data = all_data,
       aes(x = data_points,
           y = precision,
           col = as.factor(params))) +
    geom_point(size = 3, alpha = .5) +
    geom_line(data = all_data,
              aes(x = data_points,
                  y = precision,
                  col = as.factor(params)),
              alpha = .3) +
    xlim(0, x_lim) +
    theme(panel.background = element_blank(),
          legend.title = element_text(colour = "#424242", face = "bold"),
          axis.title = element_text(colour = "#424242", face = "bold", size = 10),
          axis.text  = element_text(colour = "#424242", face = "bold", size = 7)) +
#    scale_colour_discrete(name   = "L-BFGS",
#                          labels = c("T = 2",
#                                     "T = 5",
#                                     "T = 15",
#                                     "T = 20"))+

#    scale_colour_discrete(name   = "Newton Truncado con Gradiente Conjugado y \nHessiana Submuestreada",
#                          labels = c("LBFGS  | T = 2",
#                                     "CG = 20 | SGI = 0 | SGM = 0 | Muestra = 10% ",
#                                     "CG = 20 | SGI = 0 | SGM = 0 | Muestra = 100% ",
#                                     "CG = 20 | SGI = 0 | SGM = 0 | Muestra = 20%  ",
#                                     "CG = 20 | SGI = 0 | SGM = 0 | Muestra = 80% "))+
   scale_colour_discrete(name   = "L-BFGS Estocásticamente inicializado y LBFGS",
                          labels = c("LBFGS  | T = 2",
                                     "SLM | CG = 2 | T = 2 | Muestra = 100% ",
                                     "SLM | CG = 2 | T = 5 | Muestra = 10%  ",
                                     "SLM | CG = 2 | T = 5 | Muestra = 50%  ",
                                     "SLM | CG = 5 | T = 2 | Muestra = 100%  "))+
    xlab("Puntos Explorados") + ylab("Precisión")

toPrint <- "NewtonSGVaryingSampleLBFG"
toPrint <- "LBFGSVaryingMemory"
toPrint <- "SLMVaryingAllLBFG"
## -----------------------
## Save plot
## -----------------------
ggsave(paste0("../graphs/final_graphs/", toPrint, "_plot.png"), width = 10)
