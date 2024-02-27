# Example script for O'Rourke (2024) Monte Carlo simulation tutorial for latent variable models

rm(list=ls())
install.packages("MplusAutomation")
install.packages("dplyr")
library(MplusAutomation)
library(dplyr)

## Create master directory for simulation study ##
# Set up directories to switch easily between home, work, etc. - choose 1
#homepath <- "C:/users/home/mc_sim_ex/"
#filepath <- homepath
workpath <- "C:/Users/holly/Dropbox (ASU)/Conferences-Workshops-Talks/2024 CSP/mc_sim_ex_csp/"
filepath <- workpath

# Create master directory on your machine
if (!dir.exists(filepath)) {
  dir.create(filepath, recursive = TRUE)
}

# Conditions to iterate over
# factors = # of indicators in each latent variable
# n_values = sample sizes
factors <- c(4,5)
n_values <- c(100,500)
reps <- 2

# Nested loop over factors and values
for (f in factors) {
  for (n in n_values) {
    
    # Create the subdirectory path for each condition
    sub_path <- paste0(filepath, "i = ",f, "/N = ", n)
    
    #Create actual subdirectory if not yet created
    if (!dir.exists(sub_path)) {
      dir.create(sub_path, recursive = TRUE)
    }
    
    ###############################    
    # Write data generation scripts
    ###############################  
    
    # Using writeLines to write text lines to a .inp file
    # sprintf() returns character objects containing a formatted combination of input values
    # sprintf() allows us to refer to our objects and factors within the loops
    # Within the sprintf function, %s refers to string value (text), %d refers to digit (numeric value)
    # the final line in sprintf() gives the order of the objects that are being referred to
    
    # Create the file to be written
    dgscr <- paste0(sub_path, "/i", f, "_n", n, ".inp")
    writeLines(sprintf('title: MC latent variable sim example for i=%d n=%d;
  montecarlo:
  	names = x1-x%d m1-m%d y1-y%d;
    seed = 11287 ;
  	nobs = %d;
  	nreps = %d;
    REPSAVE = ALL;
  	save = i%d_%d_*.dat;

    model montecarlo:

    	x by x1@1 x2-x%d*1;
    	m by m1@1 m2-m%d*.8;
    	y by y1@1 y2-y%d*.7;

    	x1-x%d*1;
      m1-m%d*1;
      y1-y%d*1;
    	x*1;
    	m*.5;
      y*.5;

    	m on x*.5;
    	y on m*.6 x*.05;

    model:

    	x by x1@1 x2-x%d*1;
    	m by m1@1 m2-m%d*.8;
    	y by y1@1 y2-y%d*.7;

    	x1-x%d*1;
      m1-m%d*1;
      y1-y%d*1;
    	x*1;
    	m*.5;
      y*.5;

    	m on x*.5;
    	y on m*.6 x*.05;
  output:
  	tech9;', 
                       f, n, 
                       f, f, f,
                       n,
                       reps,
                       f, n,
                       f, f, f,
                       f, f, f,
                       f, f, f,
                       f, f, f),
               con = dgscr)
    ## Simulate data ##
    runModels(sub_path)
    
    ## Create file suffix naming conventions for moving the files around ##
    # file suffix for the data gen .inp/.out files
    file_suffix_n <- paste0("i",f, "_n", n)
    # file suffix for the "list.dat" files that are created during runModels()
    file_suffix_list <- paste0("i",f, "_",n, "_list")
    
    # Move data gen .inp/.out & list files to "datagen" folder before creating analysis scripts
    datagen <- paste0(filepath, "datagen")
    if (!dir.exists(datagen)) {
      dir.create(datagen, recursive = TRUE)
    }
    file.copy(paste0(sub_path, "/", file_suffix_n, ".inp"), datagen)
    file.remove(paste0(sub_path, "/", file_suffix_n, ".inp"))
    file.copy(paste0(sub_path, "/", file_suffix_n, ".out"), datagen)
    file.remove(paste0(sub_path, "/", file_suffix_n, ".out"))
    file.copy(paste0(sub_path, "/", file_suffix_list, ".dat"), datagen)
    file.remove(paste0(sub_path, "/", file_suffix_list, ".dat"))
    
    ########################################################################
    # Write & save scripts that create analysis scripts for each replication
    ########################################################################
    
    # Create the file to be written; ensure suffix matches to code below
    ascr <- paste0(sub_path, "/i", f, "_n",n, "_a.inp")
    
    # Write the Mplus script to the file
    writeLines(sprintf('[[init]]
iterators = sample;
sample = 1:%d;
filename = "i%d_n%d_[[sample]].inp";
outputDirectory = "%s";
[[/init]]

TITLE: MC latent variable sim example for i=%d n=%d rep [[sample]];
DATA:
    File is i%d_%d_[[sample]].dat;

  VARIABLE:
      NAMES = x1-x%d m1-m%d y1-y%d;
      USEVARIABLES = x1-x%d m1-m%d y1-y%d;
ANALYSIS:
MODEL:
    	x by x1@1 x2-x%d*1;
    	m by m1@1 m2-m%d*.8;
    	y by y1@1 y2-y%d*.7;

    	x1-x%d*1;
      m1-m%d*1;
      y1-y%d*1;
    	x*1;
    	m*.5;
      y*.5;

    	m on x*.5 (a);
    	y on m*.6 (b);
    	y on x*.05 (cp);

  MODEL CONSTRAINT:
  NEW(ab);
ab = a*b;
  OUTPUT: TECH1;', 
                       reps, f, n, sub_path, 
                       f, n, f, n,
                       f, f, f,
                       f, f, f,
                       f, f, f,
                       f, f, f),
               con = ascr)

    # file suffix for the creation scripts 
    file_suffix_a <- paste0("i", f, "_n",n, "_a")
    
    ## Create .inp files for all replications in a subdirectory ##
    createModels(paste0(sub_path, "/", file_suffix_a, ".inp"))
    
    #Move createmodels .inp files out of subdirectory before running model scripts
    createm <- paste0(filepath,"/modelcreation")
    if (!dir.exists(createm)) {
      dir.create(createm, recursive = TRUE)
    }
    file.copy(paste0(sub_path, "/", file_suffix_a, ".inp"), createm)
    file.remove(paste0(sub_path, "/", file_suffix_a, ".inp")) 
    
    ## Run model scripts for all replications in subdirectory ##
    runModels(sub_path)
  } #this closes the n_values loop
} # this closes the factors loop

####################################################### PICK UP HERE EDITING CODE
  
#############################################
# SAVING OUTPUT ESTIMATES FOR FUTURE ANALYSES
#############################################

# This code saves raw parameter estimates, standard errors, p values, and confidence intervals

# pick up here with correct # of matcols

# Re-initialize loops because we need the attached if statements
# matcols = # of estimates in output; NOT tech1!
for (f in factors) {
  if (f == 4) {
    matcols <- 43
  } else if (f == 5) {
    matcols <- 52
  }
  for (n in n_values) {
    # Create references to condition subdirectories & file replications within these loops
    sub_path <- paste0(filepath, "i = ",f, "/N = ", n)
    file_suffix_n <- paste0("i",f, "_n", n)
    #Use MplusAuto to read in results from all output files in the referenced directory
    # Outputs are read in as lists by condition, with list length = reps
    output <- paste0("i",f, "_", n, "_output")
    assign(output, readModels(
      target = sub_path,  
      recursive = TRUE,
      what = "parameters",
      quiet = TRUE))
    
    #Create empty matrix with for all replications by condition
    # Create one per desired set of output
    # "matcols" = # of estimates from Mplus tech1
    mat_ests <- matrix(0, nrow = reps, ncol = matcols)
    mat_se <- matrix(0, nrow = reps, ncol = matcols)
    mat_p <- matrix(0, nrow = reps, ncol = matcols)
    
    # f1-f3 convert filepath to list naming conventions used by readModels()
    f1 <- gsub(":",".",filepath)
    f2 <- gsub("/",".",f1)
    f22 <- gsub(" ",".",f2)
    f23 <- gsub("\\(",".",f22)
    f24 <- gsub("\\)",".",f23)
    f3 <- gsub("-",".",f24)
    
    # Use list data to create data frames containing output from each replication by condition
    # get() pulls values from the list
    # the following line converts the list values to a numeric object, 
    #then to a vector, then transposes the vector such that rows are replications
    # the third line referring to the matrix adds each replication to a new row in the matrix
    # This process is repeated for each set of desired output
    for (i in 1:reps) {
      filename <- paste0(f3, "i...",f, ".N...", n, ".i", f, "_n", n, "_", i, ".out")
      ## Parameter estimates ##
      ests_1 <- get(output)[[filename]][["parameters"]][["unstandardized"]][["est"]]
      ests_2 <- t(as.matrix(unlist(ests_1)))
      mat_ests[i, ] <- as.numeric(ests_2)
      ## Standard errors ##
      se_1 <- get(output)[[filename]][["parameters"]][["unstandardized"]][["se"]]
      se_2 <- t(as.matrix(unlist(se_1)))
      mat_se[i, ] <- as.numeric(se_2)
      ## p values ##
      pval_1 <- get(output)[[filename]][["parameters"]][["unstandardized"]][["pval"]]
      pval_2 <- t(as.matrix(unlist(pval_1)))
      mat_p[i, ] <- as.numeric(pval_2)
      ## Pull column names
      estnames_1 <- get(output)[[filename]][["parameters"]][["unstandardized"]][["paramHeader"]]
      te_1 <- t(estnames_1)
      estnames_2 <- get(output)[[filename]][["parameters"]][["unstandardized"]][["param"]]
      te_2 <- t(estnames_2)
    } # this closes the reps loop
    #Create data frames from matrices
      ests <- as.data.frame(mat_ests)
      se <- as.data.frame(mat_se)
      pval <- as.data.frame(mat_p)
      #Create variable names
      names(ests) <- paste0(tolower(estnames_1),".",tolower(estnames_2))
      names(se) <- paste0("se.",tolower(estnames_1),".",tolower(estnames_2))
      names(pval) <- paste0("p.",tolower(estnames_1),".",tolower(estnames_2))
    
    # Save each output as file in condition subdirectories
    write.table(ests, file=paste0(sub_path, "/", file_suffix_n, "_ests.dat"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    write.table(se, file=paste0(sub_path, "/", file_suffix_n, "_stderrs.dat"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    write.table(pval, file=paste0(sub_path, "/", file_suffix_n, "_pvals.dat"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  } # this closes the n_values loop
} # this closes the factors loop


## Combine condition-level data into one dataset per outcome with all replications, for all conditions ##

# Create empty data frame for each output
simests <- data.frame()
ses <- data.frame()
pvs <- data.frame()

# Row counter
row1 <- 1

# Re-initialize the factor loops, because we no longer need the "if" statements in the loops above
for (f in factors) {
  for (n in n_values) {
    sub_path <- paste0(filepath, "i = ",f, "/N = ", n)
    file_suffix_n <- paste0("i",f, "_n", n)
    ## Parameter estimates ##
    filename_est <- paste0(sub_path, "/i", f, "_n", n, "_ests.dat")
    esttable <- read.table(filename_est, header = TRUE)
    simests <- bind_rows(simests, esttable)
    ## Standard errors ##
    filename_se <- paste0(sub_path, "/i", f, "_n", n, "_stderrs.dat")
    setable <- read.table(filename_se, header = TRUE)
    ses <- bind_rows(ses, setable)
    ## p values ##
    filename_p <- paste0(sub_path, "/i", f, "_n", n, "_pvals.dat")
    ptable <- read.table(filename_p, header = TRUE)
    pvs <- bind_rows(pvs, ptable)
    
    # create factor & rep variables for each dataset
    for (rep in 1:reps) {
      simests[row1, "rep"] <- rep
      simests[row1, "factor"] <- f
      simests[row1, "n"] <- n
      ses[row1, "rep"] <- rep
      ses[row1, "factor"] <- f
      ses[row1, "n"] <- n
      pvs[row1, "rep"] <- rep
      pvs[row1, "factor"] <- f
      pvs[row1, "n"] <- n
      # Add a row to the dataframe
      row1 <- row1 + 1
    } # this closes the reps loop
  } # this closes the n_values loop
} # this closes the factors loop

#Save estimates data as one large .dat file in the filepath main directory
write.table(simests, file=paste0(filepath, "/all_ests.dat"), row.names=FALSE, col.names=TRUE,sep="\t", quote=FALSE)
write.table(ses, file=paste0(filepath, "/stderrs.dat"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(pvs, file=paste0(filepath, "/pvals.dat"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

