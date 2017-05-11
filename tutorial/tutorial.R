## # Tutorial for the R package *treeplyr*
## The purpose of *treeplyr* is to provide functions for matching a tree to data, and to manipulate that data using 
## *dplyr*--all while maintaining a perfect match between the tree and the data. 
##
## ## I. Getting started
## First we'll load some packages and generate some data to analyze. We're going to create a somewhat messy dataset
## and use *treeplyr* to manipulate the data to fit our needs. 
require(treeplyr)
tree <- geiger::sim.bdtree(stop="taxa", n=50, seed=1)
dat <- data.frame(X1 = rnorm(50), X2 = rnorm(50), X3 = rnorm(50),
                  taxa=c(tree$tip.label[1:40], paste("s", 51:60, sep="")), 
                  D1=sample(c("Hello", "World"), 50, replace=TRUE), D2=rbinom(50,1,0.5), 
                  D3=0.3+rbinom(50, 4, c(0.1,0.1,0.5,0.3)), XNA1 = c(rep(NA, 10), rnorm(40)))

## Let's check and see what the tree and data look like:
tree
head(dat)

## Notice that the dataset we have created has a mix of trait types, with both discrete and continuous characters, some data with missing data,
## and species names buried in the middle of the data matrix. The function *make.treedata* will search the data table for the column with 
## the most matches to the tree, and automatically use this column for matching. It will also search the rownames. Either way, the command is
## quite simple:
td <- make.treedata(tree, dat)

## We can use summary to display information about our treedata object. 
summary(td)

## We can also use indices directly on the treedata object, but note that these drop the tree:
td[[1]]
td[['X1']]
td[1:10,1:2]
## For single brackets ('[]'), we can specify that we want to keep the tip labels:
td[1:10, 1, tip.label=TRUE]

## ## II. The Basics: Reorder, Select, Filter, & Mutate
## ### *reorder*
## The treedata object itself is made up of a list of two elements *$phy* giving the tree and *$dat* providing in the data. One operation
## that is relatively common is changing the ordering of the phylogeny. It's important to maintain a match between the tree and the data. 
td <- reorder(td, "postorder")

## ### *select*
## We can use *dplyr* functions *select*, *filter* and *mutate* directly on the treedata object. The *select* function allows you to choose
## which columns you want in the dataset. Any number of columns can be specified:
select(td, X1, D1)
select(td, 1:3)
select(td, 1, 4, 6)

## ### *filter*
## We can also use the filter function allows us to select only those rows that meet a specific critierion. 
## Multiple criteria can be used to limit the dataset even more. note that as the dataset is filtered, the tree
## is automatically pruned to reflect the datasets represented.
filter(td, X1 > 0, D1=="Hello", is.na(XNA1)==FALSE)
filter(td, X1 + X2 > 0 & D1 == "Hello")

## ### *mutate*
## The function *mutate* adds a new variable to the dataset. It may be a transformation of one or more existing variables.
## For example, we may wish to log transform a variable, or average two or more variables.
mutate(td, Xall = (X1+X2+X3)/3, D1.binary = as.numeric(D1)-1)

## ## III. Programmatic use
## Note that if you want to use *treeplyr* programatically, you may not want to use these functions, as they will
## not result in the desired behavior. For example:
mytraits <- c('X1', 'D1')
try(select(td, mytraits))

## The *treeplyr* function expects a column called "mytraits" rather than evaluating the variable *mytraits*. To 
## use programatically (as with *dplyr*), use the same functions followed by an underscore:
select_(td, .dots=as.list(mytraits))

## A similar operation can be done for functions such as filter and mutate:
criteria <- list("X1 > -1", "D1 == 'Hello'", "is.na(XNA1)==TRUE")
filter_(td, .dots=criteria)

## *treeplyr* is mostly just a wrapper that passes on functions to *dplyr*, so most of the *dplyr* functionality is still there. 
## All *treeplyr* does is make sure the tree and data stay matched through the course of your data manipulations. In other words,
## you can still combine select, mutate, and filter with lots of other nice *dplyr* functions. Here are some examples:
select(td, starts_with("D"))
select(td, ends_with("1"))
select(td, matches("NA"))
select(td, contains("NA"))
select(td, which(sapply(td$dat, type_sum)=="int"))

## Or dropping columns:
select(td, -matches("NA"))
select(td, -starts_with("X"))

## ## IV. Applying functions to treedata objects: treeply, treedply and tdapply
## In many cases, the user may simply want to split apart the treedata object after matching and proceed in their analyses as normal. 
## For example, we could measure phylogenetic signal in our trait *X1*:
phytools::phylosig(td$phy, td[['X1']])

## ### *treedply*
## You can also run it directly on the treedata object using the function treedply:
treedply(td, phytools::phylosig(phy, td[['X1']], "K"))

        
## Or multiple functions at once: 
treedply(td, list("K" = phytools::phylosig(phy, td[['X1']], "K"),
                  "lambda" = phytools::phylosig(phy, td[['X1']], "lambda"))
         )

## ### *forceFactor* & *forceNumeric*
## We can also use the function *tdapply* which calls the *apply* function, but allows inclusion of the phylogeny. This can be useful for 
## applying the same function over every column in our dataset. First though, we can use the functions *forceNumeric* and *forceFactor* to 
## make sure that every column is of a type that can be analyzed by a function like *phenogram* or *phylosig*.  These functions to force 
## the traits in the treedata object to be a factor or a continuous data, dropping those that cannot be converted. 
tdDiscrete <- forceFactor(td)
tdNumeric <- forceNumeric(td)

## We can further filter out missing data in the trait *XNA1*.
tdNumeric <- filter(tdNumeric, !is.na(XNA1))

## ### *tdapply*
## Then we can apply a function like *phenogram* to plot all the data (A table of NA's is produced 
## because of the request for output from the phenogram function that returns no values):
par(mfrow=c(2,3))
tdapply(tdNumeric, 2, phytools::phenogram, tree=phy, spread.labels=FALSE, ftype="off")

## Or you could fit all traits to a BM model and then pull out the sigsq parameter:
fitsBM <- tdapply(tdNumeric, 2, geiger::fitContinuous, phy=phy, model="BM")
sapply(fitsBM, function(x) x$opt$sigsq)

## Perhaps more elegantly, you could use pipes to chain all of these operations together:
td %>% filter(., !is.na(XNA1)) %>% forceNumeric(.) %>% tdapply(., 2, phytools::phylosig, tree=phy)

### ### *treeply*
## You can manipulate the tree as well, using the function *treeply*, which is meant for simple operations on the tree alone that may 
## or may not change the number of tips. For example, let's use the *geiger* function *rescale.phylo* to rescale the branches according 
## to an OU model with an alpha value of 10. 
td.OU10 <- treeply(td, geiger::rescale, model="OU", 10)
par(mfrow=c(1,2))
plot(td$phy)
plot(td.OU10$phy)

## Or we could drop tips from the tree (here we drop tips from 1 to 35).
treeply(td, drop.tip, c(1:35))

## ## V. Grouping
## One of the cool features about *dplyr* is that it allows you to group variables and perform analyses independently 
## on different groups with a single command. Currently, *treeplyr* only supports a single grouping variable at a time, 
## but allows similar functionality. In this example, we will group taxa by the trait *D1*. You could easily imagine 
## using this for a taxonomic level on the tree.
td.D1 <- group_by(td, D1)

## ### *summarize*/*summarise*
## What good does this do? Well, we can use *summarize* to apply functions to specific groups.
summarize(td.D1, mean(X1), sd(X1), mean(X2), sd(X2))

## But what about if our functions require a phylogeny? Well, we can do that too:
summarise(td.D1, ntips = length(phy$tip.label), 
            psig.X1 = phytools::phylosig(setNames(X1, phy$tip.label), tree=phy),
              psig.X2 = phytools::phylosig(setNames(X2, phy$tip.label), tree=phy))

## Note both British and American spellings of *summarize/summarise* work.You might also want to do something like find 
## the total branch length found in different groups of taxa:
summarise(td.D1, ntips = length(phy$tip.label), 
              totalTL = sum(phy$edge.length), varianceBL = var(phy$edge.length))

## Or you could fit models of trait evolution to different groups in the tree (but note this is a dumb way
## to build a summary table, as the fitContinuous function is run independently for each parameter).
summarise(td.D1, sigsq = geiger::fitContinuous(phy, setNames(X1, phy$tip.label))$opt$sigsq, 
                 root = geiger::fitContinuous(phy, setNames(X1, phy$tip.label))$opt$z0)

## ### *paint_clades*
## We can also create a new variable that paints clades according to a particular set of nodes or branches (assuming a postordered tree). 
## by using the function *paint_clades*. 
td.painted <- paint_clades(td, interactive=FALSE, type="nodes", ids=c(75, 66, 54, 48), plot=TRUE)

## Or alternatively, you can specify which clades you want to group interactively. In the case below, the user selects 4 clades 
## by clicking on the desired branches (i.e. when you run this script,YOU must click on 4 branches of your choosing 
## to move forward!!): 
##+ eval=FALSE
td.painted <- paint_clades(td, 4, interactive=TRUE)

## Now we group the phylogeny by the *clades* variable we just defined and calculate summary statistics for each group.
td.painted <- group_by(td.painted, clades)
summarise(td.painted, ntips = length(phy$tip.label), psig1 = phytools::phylosig(setNames(X1, phy$tip.label), tree=phy), 
          meanX1 = mean(X1), sdX1 = sd(X1), ntips =length(phy$tip.label))

## Currently, the functions *treeply*, *treedply* and *tdapply* do not work with grouped data frames, and will analyze the 
## entire dataset rather than specified subgroups. This will be added in future versions of *treeplyr*. 
