rm(list=ls(all=TRUE))

#setwd("~/repos/ArborWorkflows/")
#sowsear("diversitree&Arbor.R", type="Rmd")
#knit("diversitree&Arbor.Rmd")

rm(list=ls(all=TRUE))
require(sowsear)
require(knitr)
build.knitr <- function(dir, filebase){
  #rm(list=ls(all=TRUE))
  require(sowsear)
  setwd(dir)
  sowsear(paste(filebase, ".R",sep=""), type="Rmd")
  knit(paste(filebase, ".Rmd",sep=""))
}

unlink("~/repos/bayou/figure/*")
dirw <- "~/repos/bayou/"
build.knitr(dirw, "tutorial")
y
