

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

glines <- baseline[sapply(strsplit(baseline, "\t"), length) > 180]
base_df <- read.table(textConnection(glines), sep = "\t")  # baseline data frame


