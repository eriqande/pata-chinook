

########################################################################
# Code to run simulations.
#
# Here is what I am doing:
#   1. Choose a population from the NMFS SNP baseline.
#   2. Compute the allele frequencies and then alter them by genetic drift
#      according to a certain amount of Fst (using the Dirichlet dsn).  
#   3. Create new genotypes from those allele frequencies.
#   4. Assign those genotypes back to the baseline, after having removed the
#      original focal population from the baseline.
#
# The goal is to see if fish tend to go back to the most closely related
# populations or at least back to the reporting unit.  



#### Load some libraries ####
library(digest)
library(gpiper)
library(stringr)
library(reshape2)
library(plyr)
library(dplyr)




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


# and the reporting units
rspl <- strsplit(reporting_units, "\t")
rnames <- sapply(rspl, function(x) x[1])
rpops <- strsplit(sapply(rspl, function(x) x[2]), "[, ]+")
names(rpops) <- rnames
repu_df <- stack(rpops, stringsAsFactors = FALSE)  # reporting unit data frame
names(repu_df) <- c("pop", "repu")



#### Now, let's get the allele frequencies for every population  ####

# To do this we make a tidy data frame
boing <- melt(base_df, id.vars = c("pop", "ID"))
alle <- rep("a", nrow(boing))
alle[str_detect(boing$variable, "\\.1$")] <- "b"
boing$copy <- alle
names(boing)[3:4] <- c("locus", "allele") 
gtidy <- tbl_df(boing[c(1,2,3,5,4)])
gtidy$locus <-str_replace(gtidy$locus, "\\.1$", "")
# make locus a factor in the order we want them to be in
gtidy$locus <- factor(gtidy$locus, levels = snp_names)

# here we can get allele freqs, but it sort of falls flat 
# with loci that are monomorphic within pops because we will be adding
# a prior and need to know how much prior to add.  
alle_freqs <- gtidy %>% 
  filter(!is.na(allele)) %>%
  group_by(pop, locus, allele) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n)) %>%
  dlply(.variables = c("pop", "locus"), .fun = function(x) setNames(x$freq, x$allele))

# actually I think we can rescue it by creating a prior for the 
# Dirichlet dsn that has the allele names in it
alle_names <- lapply(seq(3, ncol(base_df), by = 2), function(x) sort(unique(c(base_df[,x], base_df[,x+1]))))
alle_names <- setNames(alle_names, snp_names)
alle_prior <- lapply(alle_names, function(x) setNames(c(0.01, 0.01), x))


#### Some functions ####

# this function requires that alle_freqs, alle_priors, and snp_names
# have been appropriately defined in the global environment
# when you pass in pop, this function returns a named list of loci
# with allele frequencies simulated assuming Wrights beta distribution
# and and f value.  
drifted_freq <- function(pop, f) {
  # compute parameters of the beta distribution
  beta_pars <- lapply(names(alle_prior), function(x) {
    r = alle_prior[[x]]
    c = alle_freqs[[paste(pop, x, sep = ".")]]
    r <- r + c[names(r)]
    r[is.na(r)] <- 0.01  # total hack to put the prior back in there.
    r * (1 - f) / f
  })
  names(beta_pars) <- names(alle_prior)
  
  # simulate the beta's
  ret <- lapply(beta_pars, function(x) {
    y <- rbeta(n = 1, shape1 = x[1], shape2 = x[2])
    setNames(c(y, 1-y), names(x))
  })
  
  ret
}



# now I can do like this, which I should wrap up into a function
dfr <- drifted_freq("RogueR--Applegate_Cr", .05)

sam_size = 76
reps = 5

tmp <- lapply(dfr, function(x) sample(names(x), size = sam_size * reps * 2, replace = TRUE, prob = x))
simdat <- unlist(tmp) %>% matrix(nrow = sam_size * reps) %>% as.data.frame(stringsAsFactors = FALSE)
# now split all those into reps separate data frames
almost <- split(x = simdat, f = rep(1:reps, each = sam_size))
