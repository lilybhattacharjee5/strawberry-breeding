# generate masked files for BEAGLE

# library imports
library(vcfR)
library(data.table)
library(readr)
library(dplyr)
library(glue)
library(lattice)
library(ggplot2)

# file / path names
root_path = "~/Documents/beagle/"
provided_data_path = glue(root_path, "new_provided_data")
beagle_inputs_path = glue(root_path, "new_inputs_by_organism")
beagle_outputs_path = glue(root_path, "new_outputs_by_organism")

dir.create(beagle_inputs_path)
dir.create(beagle_outputs_path)

dir.create(glue("{beagle_inputs_path}/masking"))
dir.create(glue("{beagle_outputs_path}/masking"))

snp_array_name = "Factorial_50K_SNP_Array_Genotypes.csv"
name_mapper_name = "organism_safe_name_mapper_updated_data.csv"
cam_map_name = "camMap.csv"

snp_array = read.csv(glue("{provided_data_path}/{snp_array_name}"))
cam_map = read.csv(glue("{provided_data_path}/{cam_map_name}"))
View(snp_array)

##### GENERATE BEAGLE MASKED INPUTS

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

# necessary_cols = "(Organism*)|(Probe_Set_ID)|(Chromosome.x)|(Position)|(REF$)|(ALT)|(cM)"
necessary_cols = "(Organism*)|(Probe_Set_ID)|(Chromosome)|(Position)|(REF$)|(ALT)"
filtered_joined_snp_cam = as.data.frame(filtered_joined_snp_cam %>% select(matches(necessary_cols)))
grouped_joined_snp_cam = filtered_joined_snp_cam %>% group_by(Chromosome) %>% arrange(Position, .by_group = TRUE)
filtered_joined_snp_cam = grouped_joined_snp_cam
View(filtered_joined_snp_cam)

generateMaskedFiles <- function(runNum, numMasked) {
  dir.create(glue("{beagle_inputs_path}/masking/{runNum}_inputs/"))
  dir.create(glue("{beagle_outputs_path}/masking/{runNum}_outputs/"))
  filtered_joined_snp_cam_copy <- cbind(filtered_joined_snp_cam)
  
  # create genotypes VCF file
  chrom_linkage_groups <- filtered_joined_snp_cam_copy$Chromosome
  pos_vals <- filtered_joined_snp_cam_copy$Position
  id_vals <- rep(NA, length(pos_vals))
  ref_vals <- filtered_joined_snp_cam_copy$REF
  alt_vals <- filtered_joined_snp_cam_copy$ALT
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
    INFO = info_vals
  )
  
  # replace cells: 0 => 0/0, 1 => 0/1, 2 => 1/1
  filtered_joined_snp_cam_copy[] <- lapply(filtered_joined_snp_cam_copy, gsub, pattern = " ", replacement = "")
  filtered_joined_snp_cam_copy[filtered_joined_snp_cam_copy == "0"] <- "0/0"
  filtered_joined_snp_cam_copy[filtered_joined_snp_cam_copy == "1"] <- "0/1"
  filtered_joined_snp_cam_copy[filtered_joined_snp_cam_copy == "2"] <- "1/1"
  
  filtered_joined_snp_cam_copy_colnames <- colnames(filtered_joined_snp_cam_copy)
  filtered_genotype_data <- filtered_joined_snp_cam_copy[, filtered_joined_snp_cam_copy_colnames[filtered_joined_snp_cam_copy_colnames %in% name_mapper$safe_names]]
  filtered_genotype_data$FORMAT <- rep("GT", length(pos_vals))
  filtered_genotype_data <- filtered_genotype_data[, c(length(colnames(filtered_genotype_data)), 1:length(colnames(filtered_genotype_data)) - 1)]

  # masking
  if (numMasked > 0) {
    # randomly sample n hetLocs (for each organism)
    allHetLocs <- data.frame(
      row = c(),
      col = c()
    )
    hetLocs <- data.frame(which(filtered_genotype_data == "0/1", arr.ind = TRUE))
    
    for (i in seq(2, ncol(filtered_genotype_data))) {
      hetLocsCurr = hetLocs[hetLocs$col == i, ]
      hetLocsSelected = hetLocs[sample(nrow(hetLocs), numMasked), ]
      
      allHetLocs = rbind(allHetLocs, hetLocsSelected)
    }
    
    filtered_genotype_data_copy = cbind(filtered_genotype_data)
    
    for (h in seq(1, nrow(allHetLocs))) {
      currHetLoc <- allHetLocs[h, ]
      currX = currHetLoc$row
      currY = currHetLoc$col
      filtered_genotype_data_copy[currX, currY] = NA
    }
    
    # write masked hetLocs to file
    write.csv(allHetLocs, glue("{beagle_inputs_path}/masking/{runNum}_inputs/masked_locations.csv"))
  }

  fixedData <- as.matrix(fixedData, row.names = NULL)
  genotypeData <- as.matrix(filtered_genotype_data_copy, row.names = NULL)
  vcfData <- new("vcfR", meta = metaData, fix = fixedData, gt = genotypeData)
  write.vcf(vcfData, file = glue("{beagle_inputs_path}/masking/{runNum}_inputs/genotype.vcf.gz"))
}

numMaskedArgs <- c(0, 2, 3, 5, 10, 15, 20, 50, 100, 200, 500, 1000, 2000, 5000) # c(100, 200, 500, 1000, 2000)
for (i in seq(1, length(numMaskedArgs))) {
   generateMaskedFiles(i, numMaskedArgs[i])
}

##### CHECK BEAGLE OUTPUTS

ratio <- c()

# combine the output vcf files
for (run_i in seq(1, length(numMaskedArgs))) {
  het_phased_count <- 0
  total_het_count <- 0
  
  input_dirs <- list.dirs(path = glue("{beagle_inputs_path}/masking/{run_i}_inputs/"), full.names = TRUE)
  
  output_path <- glue("{beagle_outputs_path}/masking/{run_i}_outputs/1_phased_genotypes.vcf")
  input_path <- glue("{beagle_inputs_path}/masking/1_inputs/cleaned_genotypes.vcf")
  het_locs_path <- glue("{beagle_inputs_path}/masking/{run_i}_inputs/masked_locations.csv")
  
  if (numMaskedArgs[run_i] == 0) {
    ratio <- c(ratio, 1)
    next
  }
  
  curr_output <- read.vcfR(output_path)
  curr_input <- read.vcfR(input_path)
  
  curr_fixed_data <- as.data.table(curr_output@fix)
  curr_unphased_genotype_data <- as.data.table(curr_input@gt)
  curr_genotype_data <- as.data.table(curr_output@gt)
  curr_het_locs_data <- as.data.frame(read.csv(het_locs_path))
  
  for (j in seq(1, nrow(curr_het_locs_data))) {
    idx <- curr_het_locs_data[j, ]
    curr_out <- as.vector(unname(curr_genotype_data[idx$row]))
    curr_in <- as.vector(unname(curr_unphased_genotype_data[idx$row]))
    
    curr_in_unphased <- curr_out[[idx$col]]
    curr_out_phased <- curr_out[[idx$col]]
    
    if (!is.na(curr_out_phased) & !is.na(curr_in_unphased)) {
      print("pair")
      print(curr_out_phased)
      print(curr_in_unphased)
      if (curr_out_phased == curr_in_unphased) {
        het_phased_count = het_phased_count + 1;
      }
    }
  }
  
  total_het_count = total_het_count + nrow(curr_het_locs_data)
  
  ratio <- c(ratio, het_phased_count / total_het_count)
  print(glue("correctly phased heterozygous ratio {het_phased_count / total_het_count}"))
}

print(ratio)

# graph correctly phased heterozygous ratio
heterozygous_phased_ratio <- data.frame(num_masked = numMaskedArgs * num_linkage_groups, correctly_phased_ratio = ratio)
head(heterozygous_phased_ratio)
ggplot(data = heterozygous_phased_ratio, aes(x = num_masked, y = correctly_phased_ratio)) + geom_line() + geom_point() + ggtitle("Correctly Phased Heterozygous Ratio vs. Number of Masked Unphased Genotypes") + xlab("Number of Masked Unphased Genotypes") + ylab("Correctly Phased Heterozygous Ratio")
ggsave(glue("{beagle_outputs_path}/masking/correctly_masked_ratio.png"))

sum(is.na(genotypeData)) / (nrow(genotypeData) * (ncol(genotypeData) - 1))
