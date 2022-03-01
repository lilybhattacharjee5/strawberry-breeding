# generate masked files for BEAGLE

# library imports
library(vcfR)
library(data.table)
library(readr)
library(dplyr)
library(glue)
library(lattice)
library(ggplot2)
library(foreach)
library(doParallel)
library(ranger)

# file / path names
root_path = "~/Documents/beagle/"
provided_data_path = glue(root_path, "new_provided_data")
beagle_inputs_path = glue(root_path, "new_inputs_both_het_hom_phase_all")
beagle_outputs_path = glue(root_path, "new_outputs_both_het_hom_phase_all")

dir.create(beagle_inputs_path)
dir.create(beagle_outputs_path)

dir.create(glue("{beagle_inputs_path}/unmasked"))
dir.create(glue("{beagle_outputs_path}/unmasked"))

dir.create(glue("{beagle_inputs_path}/masked"))
dir.create(glue("{beagle_outputs_path}/masked"))

dir.create(glue("{beagle_inputs_path}/masked_both"))
dir.create(glue("{beagle_outputs_path}/masked_both"))

snp_array_name = "Factorial_50K_SNP_Array_Genotypes.csv"
name_mapper_name = "organism_safe_name_mapper_updated_data.csv"
cam_map_name = "camMap.csv"

snp_array = read.csv(glue("{provided_data_path}/{snp_array_name}"))
cam_map = read.csv(glue("{provided_data_path}/{cam_map_name}"))
View(snp_array)

##### GENERATE BEAGLE UNMASKED INPUTS

# generate vcf file metadata section
metaData = c(
  "##fileformat=VCFv4.1",
  "##source=StrawberryBreeding",
  "##phasing=partial",
  "##FILTER=<ID=PASS,Description=PASS>",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# create hashmap between unsafe names (possibly containing _) & safe names e.g. Organism1 etc.
org_names = unique(tail(colnames(snp_array), -5))
safe_names = c()
count = 0
for (i in org_names) {
  safe_names = c(glue("Organism{count}"), safe_names)
  count = count + 1
}
name_mapper = data.frame(org_names, safe_names)
write.csv(name_mapper, glue("{beagle_inputs_path}/{name_mapper_name}"))

# add the cM position to each row of the snp array table
# joined_snp_cam = merge(snp_array, cam_map, by.x = "Probe_Set_ID", by.y = "marker")
joined_snp_cam = snp_array

# replace the org names in the cleaned_mm table with the safe names
c <- colnames(joined_snp_cam)
for (i in seq(1, nrow(name_mapper))) {
  curr_org_name <- name_mapper[i, "org_names"]
  curr_safe_name <- name_mapper[i, "safe_names"]
  c[c == curr_org_name] <- curr_safe_name
}
colnames(joined_snp_cam) <- c

# remove rows where the position or cM is NA
# filtered_joined_snp_cam = joined_snp_cam %>% filter(!is.na(cM)) %>% filter(!is.na(Position))
filtered_joined_snp_cam = joined_snp_cam %>% filter(!is.na(Position))
View(filtered_joined_snp_cam)

necessary_cols = "(Organism*)|(Probe_Set_ID)|(Chromosome)|(Position)|(REF$)|(ALT)"
filtered_joined_snp_cam = as.data.frame(filtered_joined_snp_cam %>% select(matches(necessary_cols)))
grouped_joined_snp_cam = filtered_joined_snp_cam %>% group_by(Chromosome) %>% arrange(Position, .by_group = TRUE)
filtered_joined_snp_cam = grouped_joined_snp_cam
View(filtered_joined_snp_cam)

# how many NAs are in the input? on average, how many NAs per column?
organism_genotype_data <- as.data.frame(select(filtered_joined_snp_cam, contains("Organism")))
print(glue("Number of NAs in dataset: {sum(is.na(organism_genotype_data))}"))
print(glue("Number of NAs per row on average: {sum(is.na(organism_genotype_data)) / (ncol(organism_genotype_data) - 1)}"))

# where are all of the NAs in the original dataset?
# replace cells: 0 => 0/0, 1 => 0/1, 2 => 1/1
filtered_joined_snp_cam_copy <- cbind(filtered_joined_snp_cam)
filtered_joined_snp_cam_copy[] <- lapply(filtered_joined_snp_cam, gsub, pattern = " ", replacement = "")
filtered_joined_snp_cam_copy[filtered_joined_snp_cam_copy == "0"] <- "0/0"
filtered_joined_snp_cam_copy[filtered_joined_snp_cam_copy == "1"] <- "0/1"
filtered_joined_snp_cam_copy[filtered_joined_snp_cam_copy == "2"] <- "1/1"

filtered_joined_snp_cam_copy_colnames <- colnames(filtered_joined_snp_cam_copy)
filtered_joined_snp_cam_copy <- filtered_joined_snp_cam_copy[, filtered_joined_snp_cam_copy_colnames[filtered_joined_snp_cam_copy_colnames %in% name_mapper$safe_names]]
filtered_joined_snp_cam_copy <- filtered_joined_snp_cam_copy[order(names(filtered_joined_snp_cam_copy))]
View(filtered_joined_snp_cam_copy)
na_locs <- data.frame(which(is.na(filtered_joined_snp_cam_copy), arr.ind = TRUE))
write_csv(na_locs, glue("{beagle_inputs_path}/na_locs.csv"))

generateBeagleFiles <- function(beagle_inputs_path, beagle_outputs_path, 
                                beagle_df, runNum, masked = FALSE, numMasked = 0, 
                                maskBoth = FALSE, maskedOriginal = FALSE) {
  subfolder = "unmasked"
  if (masked) {
    subfolder = "masked"
  }
  if (maskBoth) {
    subfolder = "masked_both"
  }
  if (maskedOriginal) {
    subfolder = "masked_original"
  }
  dir.create(glue("{beagle_inputs_path}/{subfolder}/{runNum}_inputs/"))
  dir.create(glue("{beagle_outputs_path}/{subfolder}/{runNum}_outputs/"))
  
  beagle_df_copy <- cbind(beagle_df)
  
  # create genotypes VCF file
  chrom_linkage_groups <- beagle_df_copy$Chromosome
  pos_vals <- beagle_df_copy$Position
  id_vals <- rep(NA, length(pos_vals))
  ref_vals <- beagle_df_copy$REF
  alt_vals <- beagle_df_copy$ALT
  qual_vals <- rep(99, length(pos_vals))
  filter_vals <- rep(NA, length(pos_vals))
  info_vals <- rep(NA, length(pos_vals))
  format_vals <- rep("GT", length(pos_vals))
  
  fixedData <- data.table(
    CHROM = chrom_linkage_groups,
    POS = pos_vals,
    ID = id_vals,
    REF = ref_vals,
    ALT = alt_vals,
    QUAL = qual_vals,
    FILTER = filter_vals,
    INFO = info_vals,
    FORMAT = format_vals
  )
  
  # replace cells: 0 => 0/0, 1 => 0/1, 2 => 1/1
  beagle_df_copy[] <- lapply(beagle_df_copy, gsub, pattern = " ", replacement = "")
  beagle_df_copy[beagle_df_copy == "0"] <- "0/0"
  beagle_df_copy[beagle_df_copy == "1"] <- "0/1"
  beagle_df_copy[beagle_df_copy == "2"] <- "1/1"
  
  beagle_df_copy_colnames <- colnames(beagle_df_copy)
  filtered_genotype_data <- beagle_df_copy[, beagle_df_copy_colnames[beagle_df_copy_colnames %in% name_mapper$safe_names]]
  filtered_genotype_data <- filtered_genotype_data[order(names(filtered_genotype_data))]
  
  # masking
  if (masked) {
    if (numMasked > 0) {
      allHetLocs <- data.frame(
        row = c(),
        col = c()
      )
      
      if (maskedOriginal) {
        # where are all of the NAs in the dataset?
        allHetLocs <- cbind(na_locs)
      } else {
        # randomly sample n hetLocs (for each organism)
        hetLocs <- data.frame(which(filtered_genotype_data == "0/1", arr.ind = TRUE))
        allLocs <- data.frame(which(!is.na(filtered_genotype_data), arr.ind = TRUE))
        
        for (i in seq(1, ncol(filtered_genotype_data))) {
          if (maskBoth) {
            hetLocsCurr = allLocs[allLocs$col == i, ]
            hetLocsSelected = hetLocsCurr[sample(nrow(hetLocsCurr), numMasked), ]
          } else {
            hetLocsCurr = hetLocs[hetLocs$col == i, ]
            hetLocsSelected = hetLocsCurr[sample(nrow(hetLocsCurr), numMasked), ]
          }
          
          allHetLocs = rbind(allHetLocs, hetLocsSelected)
        }
      }
      
      filtered_genotype_data_copy = cbind(filtered_genotype_data)
      
      for (h in seq(1, nrow(allHetLocs))) {
        currHetLoc <- allHetLocs[h, ]
        currX = currHetLoc$row
        currY = currHetLoc$col
        filtered_genotype_data_copy[currX, currY] = NA
      }
      
      # write masked hetLocs to file
      allHetLocs <- allHetLocs[order(allHetLocs$col), ]
      write.csv(allHetLocs, glue("{beagle_inputs_path}/{subfolder}/{runNum}_inputs/masked_locations.csv"))
      
      filtered_genotype_data = filtered_genotype_data_copy
    }
  }
  
  fixedData <- as.matrix(fixedData, row.names = NULL)
  genotypeData <- as.matrix(filtered_genotype_data, row.names = NULL)
  vcfData <- new("vcfR", meta = metaData, fix = fixedData, gt = genotypeData)
  write.vcf(vcfData, file = glue("{beagle_inputs_path}/{subfolder}/{runNum}_inputs/genotype.vcf.gz"))
}

# generate source of truth
generateBeagleFiles(beagle_inputs_path, beagle_outputs_path, filtered_joined_snp_cam, 1)

# phase source of truth here by running beagle -- OUTSIDE SCRIPT
# terminal command: for run in {1..1} ; do d=(./new_inputs_full_phase/unmasked/${run}_inputs/) ; pattern="${d}/genotype.vcf" ; files=( $pattern ) ; tr -d " " < ${files[0]} > ${d}/cleaned_genotypes.vcf ; echo ${d}/cleaned_genotypes.vcf; mkdir /Users/lilybhattacharjee/Documents/beagle/new_outputs_full_phase/unmasked/${run}_outputs/; java -Xmx1g -jar beagle.28Jun21.220.jar gt=$d/cleaned_genotypes.vcf out=/Users/lilybhattacharjee/Documents/beagle/new_outputs_full_phase/unmasked/${run}_outputs/phased_genotypes ; done ;

# convert output to unphased data
# read output VCF file
phased_genotypes = read.vcfR(glue("{beagle_outputs_path}/unmasked/1_outputs/phased_genotypes.vcf"))
phased_genotype_data = phased_genotypes@gt
fixed_data = as.data.frame(phased_genotypes@fix)
phased_genotype_data[phased_genotype_data == "0|0"] <- "0/0"
phased_genotype_data[phased_genotype_data == "0|1" | phased_genotype_data == "1|0"] <- "0/1"
phased_genotype_data[phased_genotype_data == "1|1"] <- "1/1"
View(phased_genotype_data)

print(glue("There are no NAs in the phased 'source of truth' -- {sum(is.na(phased_genotype_data))}"))

phased_genotype_data_df <- as.data.frame(phased_genotype_data)
phased_genotype_data_df$Probe_Set_ID = filtered_joined_snp_cam$Probe_Set_ID
phased_genotype_data_df$Chromosome = fixed_data$CHROM
phased_genotype_data_df$Position = fixed_data$POS
phased_genotype_data_df$REF = fixed_data$REF
phased_genotype_data_df$ALT = fixed_data$ALT

# remove some data by adding NAs
numMaskedArgs <- c(0, 2, 3, 5, 10, 15, 20, 50, 100, 200, 300, 500, 600, 700, 1000) # ~500-600 NAs / col on average
# will add 15, 20, 50, 100, 200, 300, 500, 600, 700, 1000

print("Generating masked het locations ONLY")
n.cores <- parallel::detectCores() - 1

# create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)

foreach::getDoParRegistered()

start = 1
foreach(
  i = seq(start, start + length(numMaskedArgs) - 1)
) %dopar% {
  library(glue)
  library(data.table)
  library(vcfR)
  generateBeagleFiles(beagle_inputs_path, beagle_outputs_path, phased_genotype_data_df, i, masked = TRUE, numMasked = numMaskedArgs[i - start + 1])
  print(glue("finished generating {numMaskedArgs[i - start + 1]} inputs"))
}

print("Generating masked het AND hom locations")
start = 1
foreach(
  i = seq(start, start + length(numMaskedArgs) - 1)
) %dopar% {
  library(glue)
  library(data.table)
  library(vcfR)
  generateBeagleFiles(beagle_inputs_path, beagle_outputs_path, phased_genotype_data_df, i, masked = TRUE, numMasked = numMaskedArgs[i - start + 1], maskBoth = TRUE)
  print(glue("finished generating {numMaskedArgs[i - start + 1]} inputs"))
}

print("Generating masked original locations")
generateBeagleFiles(beagle_inputs_path, beagle_outputs_path, phased_genotype_data_df, "original", masked = TRUE, numMasked = 1, maskedOriginal = TRUE)
print(glue("finished generating masked original inputs"))

# phase newly generated files by running beagle -- OUTSIDE SCRIPT (do this on Farm)
# terminal command: for run in {1..11} ; do d=(./new_inputs_full_phase/masked/${run}_inputs/) ; pattern="${d}/genotype.vcf" ; files=( $pattern ) ; tr -d " " < ${files[0]} > ${d}/cleaned_genotypes.vcf ; echo ${d}/cleaned_genotypes.vcf; mkdir /Users/lilybhattacharjee/Documents/beagle/new_outputs_full_phase/masked/${run}_outputs/; java -Xmx1g -jar beagle.28Jun21.220.jar gt=$d/cleaned_genotypes.vcf out=/Users/lilybhattacharjee/Documents/beagle/new_outputs_full_phase/masked/${run}_outputs/phased_genotypes ; done ;

##### CHECK BEAGLE OUTPUTS
# for masked (i.e. only het genotypes masked), check 2 metrics
# - what percent of *unmasked* het genotypes are phased correctly?
# - what percent of *masked* het genotypes are phased correctly?

# for masked_both (i.e. het and hom genotypes masked, potentially), check 3 metrics
# - what percent of het / hom masked genotypes are imputed properly as het / hom?
# - what percent of *unmasked* het genotypes are phased correctly?
# - what percent of *masked* het genotypes are phased correctly?

calculate_metrics <- function(mode) {
  correct_masked_het_ratios <- c()
  correct_unmasked_het_ratios <- c()
  correct_imputed_ratios <- c()
  
  for (run_i in seq(1, 15)) { #length(numMaskedArgs)
    het_phased_count <- 0
    total_het_count <- 0
    
    input_dirs <- list.dirs(path = glue("{beagle_inputs_path}/masking/{run_i}_inputs/"), full.names = TRUE)
    
    output_path <- glue("{beagle_outputs_path}/{mode}/{run_i}_outputs/phased_genotypes.vcf")
    input_path <- glue("{beagle_inputs_path}/{mode}/{run_i}_inputs/cleaned_genotypes.vcf")
    het_locs_path <- glue("{beagle_inputs_path}/{mode}/{run_i}_inputs/masked_locations.csv")
    source_of_truth_path <- glue("{beagle_outputs_path}/unmasked/1_outputs/phased_genotypes.vcf")
    
    curr_output <- read.vcfR(output_path)
    curr_input <- read.vcfR(input_path)
    source_of_truth <- read.vcfR(source_of_truth_path)
    
    curr_input_genotype_data <- as.data.frame(curr_input@gt)
    curr_output_genotype_data <- as.data.frame(curr_output@gt)
    curr_sot_genotype_data <- as.data.frame(source_of_truth@gt)
    
    curr_input_genotype_data <- curr_input_genotype_data[, seq(2, ncol(curr_input_genotype_data))]
    curr_output_genotype_data <- curr_output_genotype_data[, seq(2, ncol(curr_output_genotype_data))][, names(curr_input_genotype_data)]
    curr_sot_genotype_data <- curr_sot_genotype_data[, seq(2, ncol(curr_sot_genotype_data))][, names(curr_input_genotype_data)]
    
    na_pos = data.frame(which(is.na(curr_input_genotype_data), arr.ind = TRUE))
    
    curr_het_locs_data <- as.data.frame(data.table(
      row = c(),
      col = c()
    ))
    
    tryCatch({
      curr_het_locs_data <- as.data.frame(read.csv(het_locs_path))
      curr_het_locs_data <- curr_het_locs_data[order(curr_het_locs_data$col), ]
    },
    error = function(cond) {
    }
    )
    
    tot_rows <- seq(1, nrow(curr_output_genotype_data))
    
    correct_masked_het = 0
    total_masked_het = 0
    
    correct_unmasked_het = 0
    total_unmasked_het = 0
    
    correct_imputed = 0
    total_masked = 0
    
    for (i in seq(1, ncol(curr_output_genotype_data))) {
      masked_locs <- curr_het_locs_data[curr_het_locs_data$col == i, ]
      
      curr_col_unmasked_rows <- setdiff(tot_rows, masked_locs$row)
      curr_col_masked_rows <- masked_locs$row
      
      curr_unmasked_vals_sot <- curr_sot_genotype_data[curr_col_unmasked_rows, i]
      curr_masked_vals_sot <- curr_sot_genotype_data[curr_col_masked_rows, i]
      
      curr_unmasked_vals_out <- curr_output_genotype_data[curr_col_unmasked_rows, i]
      curr_masked_vals_out <- curr_output_genotype_data[curr_col_masked_rows, i]
      
      # # what percent of *unmasked* het genotypes are phased correctly?
      curr_unmasked_01 <- (curr_unmasked_vals_sot == "0|1") & (curr_unmasked_vals_out == "0|1")
      curr_unmasked_10 <- (curr_unmasked_vals_sot == "1|0") & (curr_unmasked_vals_out == "1|0")
      curr_correct_unmasked_het = sum(curr_unmasked_01) + sum(curr_unmasked_10)
      curr_total_unmasked_het = sum((curr_unmasked_vals_sot == "0|1")) + sum((curr_unmasked_vals_sot == "1|0"))
      
      # what percent of *masked* het genotypes are phased correctly?
      curr_masked_01 <- (curr_masked_vals_sot == "0|1") & (curr_masked_vals_out == "0|1")
      curr_masked_10 <- (curr_masked_vals_sot == "1|0") & (curr_masked_vals_out == "1|0")
      curr_correct_masked_het = sum(curr_masked_01) + sum(curr_masked_10)
      curr_total_masked_het = sum((curr_masked_vals_sot == "0|1")) + sum((curr_masked_vals_sot == "1|0"))
      
      # what percent of het / hom *masked* genotypes are imputed properly as het / hom?
      if (mode == "masked_both") {
        # het that are imputed as het
        curr_masked_het <- sum(((curr_masked_vals_out == "0|1") | (curr_masked_vals_out == "1|0")) & ((curr_masked_vals_sot == "0|1") | (curr_masked_vals_sot == "1|0")))
        # hom that are imputed as hom
        curr_masked_hom <- sum(((curr_masked_vals_out == "0|0") | (curr_masked_vals_out == "1|1")) & ((curr_masked_vals_sot == "0|0") | (curr_masked_vals_sot == "1|1")))
        
        curr_correct_imputed <- curr_masked_het + curr_masked_hom
        curr_total_masked = nrow(masked_locs)
        
        correct_imputed = correct_imputed + curr_correct_imputed
        total_masked = total_masked + curr_total_masked
      }
      
      correct_unmasked_het = correct_unmasked_het + curr_correct_unmasked_het
      total_unmasked_het = total_unmasked_het + curr_total_unmasked_het
      
      correct_masked_het = correct_masked_het + curr_correct_masked_het
      total_masked_het = total_masked_het + curr_total_masked_het
    }
    correct_unmasked_het_ratios <- c(correct_unmasked_het_ratios, correct_unmasked_het / total_unmasked_het)
    correct_masked_het_ratios <- c(correct_masked_het_ratios, correct_masked_het / total_masked_het)
    
    if (mode == "masked_both") {
      correct_imputed_ratios <- c(correct_imputed_ratios, correct_imputed / total_masked)
      print(glue("imputed ratios {correct_imputed_ratios}"))
    }
    
    print(glue("unmasked het ratios {correct_unmasked_het_ratios}"))
    print(glue("masked het ratios {correct_masked_het_ratios}"))
  }
  return(list(correct_unmasked_het_ratios, correct_masked_het_ratios, correct_imputed_ratios))
}

metrics <- calculate_metrics("masked")
correct_unmasked_het_ratios
correct_unmasked_het_ratios = metrics[1]
correct_masked_het_ratios = metrics[2]
correct_imputed_ratios = metrics[3]

print(correct_unmasked_het_ratios)
print(correct_masked_het_ratios)
print(correct_imputed_ratios)

metrics_both <- calculate_metrics("masked_both")
correct_unmasked_het_ratios_both = metrics_both[1]
correct_masked_het_ratios_both = metrics_both[2]
correct_imputed_ratios_both = metrics_both[3]

print(correct_unmasked_het_ratios_both)
print(correct_masked_het_ratios_both)
print(correct_imputed_ratios_both)

metrics_both[3]

# save older results
unmasked_het <- data.frame(numMaskedArgs, metrics[1], metrics_both[1])
names(unmasked_het) <- c("num_masked", "correct_unmasked_het_ratios", "correct_unmasked_het_ratios_both")
ggplot(data = unmasked_het, aes(x = num_masked)) + geom_line(aes(y = correct_unmasked_het_ratios), color = "blue") + geom_line(aes(y = correct_unmasked_het_ratios_both), color = "red") + geom_point(aes(y = correct_unmasked_het_ratios)) + geom_point(aes(y = correct_unmasked_het_ratios_both)) + ggtitle("Correctly Phased Unmasked Heterozygous Ratio vs. Number of Masked Unphased Genotypes") + xlab("Number of Masked Unphased Genotypes") + ylab("Correctly Phased Unmasked Heterozygous Ratio")
ggsave(glue("{beagle_outputs_path}/correct_unmasked_het_ratios.png"))

masked_het <- data.frame(numMaskedArgs, metrics[2], metrics_both[2])
names(masked_het) <- c("num_masked", "correct_masked_het_ratios", "correct_masked_het_ratios_both")
ggplot(data = masked_het, aes(x = num_masked)) + geom_line(aes(y = correct_masked_het_ratios), color = "blue") + geom_line(aes(y = correct_masked_het_ratios_both), color = "red") + geom_point(aes(y = correct_masked_het_ratios)) + geom_point(aes(y = correct_masked_het_ratios_both)) + ggtitle("Correctly Phased Masked Heterozygous Ratio vs. Number of Masked Unphased Genotypes") + xlab("Number of Masked Unphased Genotypes") + ylab("Correctly Phased Masked Heterozygous Ratio")
ggsave(glue("{beagle_outputs_path}/correct_masked_het_ratios.png"))

imputed <- data.frame(numMaskedArgs, metrics_both[3])
names(imputed) <- c("num_masked", "correct_imputed_ratios")
ggplot(data = imputed, aes(x = num_masked)) + geom_line(aes(y = correct_imputed_ratios), color = "red") + geom_point(aes(y = correct_imputed_ratios)) + ggtitle("Correctly Imputed Ratio vs. Number of Masked Unphased Genotypes") + xlab("Number of Masked Unphased Genotypes") + ylab("Correctly Imputed Ratio")
ggsave(glue("{beagle_outputs_path}/correct_imputed_ratios.png"))
