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
    x_lim <- 5e6
}else{
    x_lim <- as.integer(args[1])
}

## -----------------------
## Read in dataset
## -----------------------

## List all Directories:
files <- list.files("../tests")
files <- files[!str_detect(files, ".sh")]

## List all subfiels
all_files <- llply(files,
                  function(t) t <- list.files(paste0("../tests/", t)))
all_files <- llply(all_files, function(t) t <- t[!str_detect(t, ".txt")])

## All data
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
    geom_point() +
    geom_line(data = all_data,
              aes(x = data_points,
                  y = precision,
                  col = as.factor(params))) +
    xlim(0, x_lim) +
    theme(panel.background = element_blank())

## -----------------------
## Save plot
## -----------------------
ggsave("../graphs/all_plot.png", width = 10)
