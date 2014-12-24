

########################################################################
# Code to run simulations.
#
# Here is what I am doing:
#   1. Choose a population from the NMFS SNP baseline.
#   2. Compute the allele frequencies and then alter them by genetic drift
#      according to a certain amount of Fst (using the Dirichlet dsn, or, in the 
#      case of SNPs, the beta dsn).  
#   3. Create new genotypes from those allele frequencies.
#   4. Assign those genotypes back to the baseline, after having removed the
#      original focal population from the baseline.
#
# The goal is to see if fish tend to go back to the most closely related
# populations or at least back to the reporting unit.  


# this flag controls whether the simulations are totally done over (which can take a
# rather long time)
REDO_SIMS = TRUE

#### Load some libraries ####
library(digest)
library(gpiper)
library(stringr)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)



#### Get the baseline data from Dryad and read it in  ####

if( file.exists("baseline.zip") && digest("baseline.zip", algo = "sha1") == "ff9e7854e6ec22359e2f151f35db467e4faabe10") {
  message("You already have the correct baseline file.  No need to download again")
} else {
  message("Downloading the baseline file.")
  DryadURL <- "http://datadryad.org/bitstream/handle/10255/dryad.61648/Clemento_etal_2014_Chinook_SNP_baseline.zip"
  download.file(url = DryadURL, destfile = "baseline.zip")
}

# define some connections into the zip file
bz = unz(description = "baseline.zip", filename = "2010_Baseline_SNPset.txt")
rz = unz(description = "baseline.zip", filename = "ReportingUnits_SNPset_2010.txt")
sz = unz(description = "baseline.zip", filename = "SNP_names.txt")

baseline <- readLines(con = bz)
reporting_units <- readLines(con = rz)
snp_names <- readLines(con = sz)

# close the connections
invisible(lapply(list(bz, rz, sz), close))



#### Now process those inputs into some nicer formats ####

# pull together the baseline
glines <- baseline[sapply(strsplit(baseline, "\t"), length) > 180]
base_df <- read.table(textConnection(glines), sep = "\t", stringsAsFactors = FALSE, na.strings = "0")  # baseline data frame
names(base_df) <- c("ID", paste(rep(snp_names, each = 2), c("", ".1"), sep = ""))
# add a pop column
pop <- str_match(base_df$ID, "(^.*):")[,2]
base_df <- cbind(pop = pop, base_df, stringsAsFactors = FALSE)


if(FALSE) {  #the reporting unit file on dryad is not correct.  I leave this in here just to have it later
# and the reporting units
rspl <- strsplit(reporting_units, "\t")
rnames <- sapply(rspl, function(x) x[1])
rpops <- strsplit(sapply(rspl, function(x) x[2]), "[, ]+")
names(rpops) <- rnames
repu_df <- stack(rpops, stringsAsFactors = FALSE)  # reporting unit data frame
names(repu_df) <- c("pop", "repu")
repu_df$repu <- factor(repu_df$repu, levels <- unique(repu_df$repu))
# and, unfortunately, there were two duplicated populations in the repunits file, so we need to chuck those
repu_df <- repu_df[!duplicated(repu_df),]
# not only that, but it appears there are populations in the repunits file that don't appear in the
# baseline.  We will remove them from the repunits here
unique(base_df$pop) %in% repu_df$pop 
}

# actually we are just going to get the reporting units from the names of the pops
repu_df <- str_match(unique(base_df$pop), "^(.*)--.*") %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
  setNames(c("pop", "repu"))
repu_df$repu <- factor(repu_df$repu, levels = unique(repu_df$repu))


#### Now, let's get the allele frequencies for every population  ####


# to count up alleles, we run it through gsi_sim and then
# process the output into afreqs_list
cat(baseline, sep = "\n", file = "full_baseline.txt")
gsi_Run_gsi_sim(" -b full_baseline.txt")
dumpola <- readLines("GSI_SIM_Dumpola.txt")

# get the allele counts 
ac <- dumpola[str_detect(dumpola, "^IND_COLL_SUMMARY_popfreqsum: Pop")] %>%
  str_replace("^.*Pop=", "") %>%
  str_split("[ \t]") %>% 
  sapply(function(x) x[c(1,2,4)]) %>%
  t() %>%
  as.data.frame(stringsAsFactors = FALSE)

ac[,2] <- as.numeric(ac[,2])
ac[,3] <- as.numeric(ac[,3])

# and the allele names
an <- dumpola[str_detect(dumpola, "^IND_COLL_SUMMARY_popfreqsum: LOCUS")] %>%
  str_replace("^IND_COLL_SUMMARY_popfreqsum: LOCUS=", "") %>%
  str_split("\t") %>%
  sapply(., function(x) x) %>% 
  t()
rownames(an) <- an[,1]
an <- an[,-1]


# then make a list, with one component for each locus with the allele freqs as
# a matrix.  We add one observed gene copy of each allele to each population
# to avoid zeroes.
snp_name_factor <- factor(snp_names, levels = snp_names)  # otherwise split sorts things.  Such 
      # stupid friggin behavior that!  No wondery hadley wrote plyr.
afreqs_list <- split(ac, f = rep(snp_name_factor, each = length(unique(ac[,1])))) %>%
  lapply(., function(x) {rownames(x) <- x[,1]; x[,-1]}) %>%
  lapply(., function(x) {tt <- rowSums(x); 
                         x[,1] <- (x[,1] + 1) / (tt + 2);
                         x[,2] <- (x[,2] + 1) / (tt + 2); 
                         x})

# add the allele names on as the colnames of each list component:
afreqs_list <- lapply(names(afreqs_list), function(x) {
  y <- afreqs_list[[x]]
  colnames(y) <- an[x,]
  y
}) %>%
  setNames(., names(afreqs_list))


# and for one last hurrah here we need to change the name of NOregonCoast--Nehalem_
# to NOregonCoast--Nehalem_R.  It is named incorrectly on the POP line in the
# baseline file.
afreqs_list <- lapply(afreqs_list, function(x) {rownames(x)[27] <- "NOregonCoast--Nehalem_R"; x})


#### Some functions ####

# this function requires that afreqs_list
# has been appropriately defined in the global environment
# when you pass in pop (a string with the population name), this function returns a named list of loci
# with allele frequencies simulated assuming Wrights beta distribution
# and an f value.  
drifted_freq <- function(pop, f) {
  # compute parameters of the beta distribution
  beta_pars <- lapply(afreqs_list, function(x) {
    r = x[pop,]
    r * (1 - f) / f
  })
  
  # simulate the beta's
  ret <- lapply(beta_pars, function(x) {
    y <- rbeta(n = 1, shape1 = x[1, 1], shape2 = x[1, 2])
    setNames(c(y, 1-y), names(x))
  })
  
  ret
}



## This is a higher level function.  Again, it assumes that 
# afreqs_list has been appropriately
# defined in the global environment.  You pass it a population name
# (pop), a sample size for the simulated sample (sam_size), the number
# of replicated simulated samples to make (reps) and a drift parameter (f),
# the baseline data frame (baseline), and the repunits data frame (repunits).
# This function then does this:
# 1. simulates reps simulated data sets using drifted versions of pop's alle freqs.
#    These are returned in $samples which is a list of data frames
# 2. removes pop from the baseline and repunits and returns what remains of each of those
#    in $baseline and $repunits
sim_the_reps <- function(pop, reps, f, sam_size = 76, baseline = base_df, repunits = repu_df, noStrip = FALSE) {
  
  # get the drifted allele freqs
  dfr <- drifted_freq(pop, f)
  
  # now, make a list of simulated sample data frames
  tmp <- lapply(dfr, function(x) sample(names(x), size = sam_size * reps * 2, replace = TRUE, prob = x))
  simdat <- unlist(tmp) %>% 
    matrix(nrow = sam_size * reps) %>%  
    as.data.frame(stringsAsFactors = FALSE)
  
  # now split all those into reps separate data frames and add appropriate names to them
  samples <- split(x = simdat, f = rep(1:reps, each = sam_size)) %>%
    lapply(., function(x) {
      y <- setNames(x, names(baseline[-(1:2)]))
      rownames(y) <- paste("Sim", f, pop, 1:sam_size, sep = "_")
      y
      })
  
  # now strip pop from baselines and repunits
  if(noStrip == FALSE) {
    baseline_new <- baseline[baseline$pop != pop, ]
    repunits_new <- repunits[repunits$pop != pop, ]
  } else {
    baseline_new <- baseline
    repunits_new <- repunits
  }
  
  # and return the result
  list(samples =  samples, baseline = baseline_new, repunits = repunits_new)
  
}


# this function slurps up the gsi_sim results
slurp_gsi_sim <- function() {
  x <- read.table("rep_unit_pofz_full_em_mle.txt", header=T)
  
  # get the max posterior reporting unit and what that max posterior is
  max_repu <- names(x)[-(1:2)][max.col(x[,-(1:2)])]
  max_post <- apply(x[,-(1:2)], 1, max)
  
  # grab the sim-logls too
  y <- read.table("em_mixture_logl_summary.txt", header=T)
  
  
  # return it all as a big data frame
  data.frame(IndivName = x$IndivName, 
             max_repu, 
             max_post, 
             ObsLogL = y$ObsLogL,
             zScore = y$zScore,
             FractionSimmedGreater = y$FractionSimmedGreater,
             stringsAsFactors = FALSE)
}




## This function takes what gets returned by sim_the_reps (R)
#  and runs gsi_sim on it and slurps the results out nicely.
chuck_it_through_gsi_sim <- function(R) {
  
  # make the baseline:
  popf <- factor(R$baseline[, 1], levels = unique(R$baseline[, 1]))
  basein <- R$baseline
  rownames(basein) <- as.character(basein$ID)
  basein <- basein[,-(1:2)]
  basein[is.na(basein)] <- 0
  gPdf2gsi.sim(basein, pop.ize.them = popf)
  
  # make the reporting units file
  gsi_WriteReportingUnitsFile(R$repunits$pop, R$repunits$repu)
  
  lapply(R$samples, function(x) {
    # here is how we make a mixture file:
    gPdf2gsi.sim(x, outfile = "simmixfile.txt")
    
    # and then here we run it through gsi_sim
    gsi_Run_gsi_sim(" -b gsi_sim_file.txt  -t simmixfile.txt -r repunits.txt --mix-logl-sims 1000 0")
    
    # slurp out the results, find maxes and return as a data frame
    slurp_gsi_sim()
  })
  
}


# finally, here is a function that takes a population name and then runs it
# at a series of different f's and sends it all back as a data frame
run_a_pop <- function(pop,  f = c(.000001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.5 ), reps = 10, ...) {
  tmp <- lapply(f, function(x) {
    R <- sim_the_reps(pop, reps, x, ...)
    boing <- chuck_it_through_gsi_sim(R)
    output <- do.call(rbind, lapply(names(boing), function(x) cbind(rep = x, boing[[x]])) )
    output
    })
  
  do.call(rbind, tmp)
}


#### Now run this over all the different populations ####


# What I want to do now is run it for every population in the baseline for something
# like 10 reps at each of 5 different values of f.  Then summarize the output
# graphically.  One way would be to show the fraction of each rep data set that
# ended up being assigned to the Willamette (or lower columbia, etc).  Those could be boxplots for each
# and every baseline population.  Faceted over f.


# Here I get a table of populations that tells us how many pops are in each of the
# reporting units they are in
pops_in_repu_counts <- repu_df %>%
  tbl_df() %>%
  group_by(repu) %>%
  tally() %>%
  inner_join(., repu_df) %>%
  rename(number_of_pops_in_repu = n)

# I will use the names of the pops in that to cycle over:
pops_cycle  <- pops_in_repu_counts$pop
pops_cycle <- pops_cycle[!(pops_cycle %in% c("CohoSp--California_Coho"))]  # remove coho


if(REDO_SIMS == TRUE) {
  list_output <- lapply(pops_cycle, run_a_pop, reps = 10)
  names(list_output) <- pops_cycle
  saveRDS(list_output, "all_pops_sim_output.rds", compress = "xz")
  stop("That's as far as we go with REDO_SIMS == TRUE")
} else {
  list_output <- readRDS("all_pops_sim_output.rds")
}


#### Now, put the list_output into long format and prepare it for plotting ####

simdf <- do.call(rbind, list_output)
rownames(simdf) <- NULL

# now strip some values out of the IndivNames and add the to the data frame
tmp <- strsplit(as.character(simdf$IndivName), "_")
simdf$f_value <- as.numeric(sapply(tmp, "[", 2))
simdf$source_pop <- sapply(tmp, function(x) paste(x[-c(1,2, length(x))], collapse = "_")) %>%
  factor(., levels = pops_cycle)

# now, summarize with dplyr
simdf <- tbl_df(simdf)

# here is the number of fish from each source assigned to each max_repu
mr_counts <- simdf %>%
  group_by(source_pop, f_value, rep, max_repu) %>%
  tally()


# here we filter it to include only those fish assigned back to "WillametteR"
ppt_to_willam <- mr_counts %>%
  filter(max_repu == "WillametteR") %>%
  mutate(ppn = n / 76)


# now, let's make a quick plot:
ggplot(data = ppt_to_willam, aes(y = source_pop, x = ppn)) + 
  geom_point(color = "red") + 
  facet_wrap(~ f_value)


# OK. that all looks good, but there are two problems that must be dealt with:
# 1. I accidentally left reps = 2 in there, instead of 10
# 2. The two reps weren't actually different!  Something is wrong there---the reps look like
#    they might be perfect replicas of one another (the gsi-sim results are all the same.  I need
#    track that down!)