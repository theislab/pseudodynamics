rm(list=ls())
library(methods)
library(ggplot2)
dirRoot <- "your_root_directory_of_this_project"

### Load scanpy data
dfDMKlein15 <- read.csv(paste0(dirRoot, "dimreds/DiffMap_Klein15.csv"), row.names = NULL, header = FALSE)
dfDPTKlein15 <- read.csv(paste0(dirRoot, "dimreds/dpt_Klein15.csv"), header=T)
if(FALSE) { # check groups
    dfgplotTest <- as.data.frame(dfDMKlein15)
    colnames(dfgplotTest) <- paste0("DC", seq(1,dim(dfgplotTest)[2]))
    dfgplotTest <- cbind(dfgplotTest, dfDPTKlein15)
    dfgplotTest$time <- as.character(dfgplotTest$time)
    ggplot() + geom_point(data = dfgplotTest, aes(
        x = DC1, y = DC3, colour = time))
    ggplot() + geom_point(data = dfgplotTest, aes(
        x = DC1, y = DC3, colour = dpt_pseudotime))
}

dfDynamicsInput_DPT_klein15 <- data.frame(
    pseudodynamics_branching_region="nonbranching",
    pseudodynamics_branch_class="main",
    dpt_adjusted_groups="none",
    time=c(0,2,4,7)[as.numeric(dfDPTKlein15$time)],
    batch=paste0("timebatch_", as.character(dfDPTKlein15$time)),
    pseudotime=dfDPTKlein15$dpt_pseudotime,
    cell=paste0("c", seq(1, dim(dfDPTKlein15)[1])),
    dc1=dfDMKlein15[,1],
    dc2=dfDMKlein15[,2],
    dc3=dfDMKlein15[,3],
    dc4=dfDMKlein15[,4],
    stringsAsFactors = FALSE
)

# sanity test
if(FALSE) {
    ggplot() + geom_point(data = dfDynamicsInput_DPT_klein15, aes(
        x = dc1, y=dc3, colour = pseudotime))
    ggplot() + geom_point(data = dfDynamicsInput_DPT_klein15, aes(
        x = dc1, y=dc3, colour = time))
}

write.table(dfDynamicsInput_DPT_klein15, file=paste0(
    dirRoot, "input/pseudodynamics_input_Klein15_dpt_DPTpaperProcessedData.csv"),
            sep=",", quote = FALSE, row.names = FALSE)
