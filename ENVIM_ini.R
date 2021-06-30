library(optparse, quietly = TRUE)

option_list <- list(
  make_option(c("-i","--pathENVIM"),type="character",
              help="Path of R file of ENVIM.R. If there is no input, we assume you use the same directory.",default = F),
  make_option(c("-a","--microtrain"), type="character",help="Path of microbial training data"),
  make_option(c("-b","--microtest"), type="character",help="Path of microbial testing data"),
  make_option(c("-c","--metabtrain"), type="character",help="Path of metabolite training data"),
  make_option(c("-d","--metabtest"), type="character",help="Path of metabolite testing data"),
  make_option(c("-e","--seed"), type="integer",
              help="Seed number for ENVIM to keep reproducible",default = 1234),
  make_option(c("-f","--foldrf"), type="integer",
              help="Fold number for cross-validated random forest model to find variable importance",default = 10),
  make_option(c("-g","--foldENVIM"), type="integer",
              help="Fold number for cross-validated ENVIM (elastic net model based on variable importance features) to find variable importance",default = 10),
  make_option(c("-o","--output"),type="character",
              help="Output directory. If there is no input, we assume you use the same directory",default = F)
)
parser <- OptionParser(usage="%prog [options] file", option_list=option_list)

args <- parse_args(parser)

if (args$pathENVIM == F){
  source(paste0(getwd(),"/ENVIM.R"))
  }else {
    source(args$pathENVIM)
    }

if (args$o == F){
  out = getwd()
}else {
  out = args$o
}

microbio.train = read.csv(args$microtrain,row.names = 1)
microbio.test = read.csv(args$microtest,row.names = 1)
metab.train = read.csv(args$metabtrain,row.names = 1)
metab.test = read.csv(args$metabtest,row.names = 1)

ENVIM(microbio.train = microbio.train,
      microbio.test = microbio.test,
      metab.train = metab.train,
      metab.test = metab.test,
      seed = args$seed,
      outputdirectory = out,
      fold_rf = args$foldrf,
      fold_ENVIM = args$foldENVIM)




