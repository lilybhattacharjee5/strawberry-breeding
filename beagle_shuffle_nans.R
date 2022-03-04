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

load_data <- function(root_path, beagle_inputs_path, beagle_outputs_path) {
  # create inputs folder + subfolders
  dir.create(beagle_inputs_path)
  dir.create(glue("{beagle_inputs_path}/original"))
  dir.create(glue("{beagle_inputs_path}/masked"))
  
  # created outputs folder + subfolders
  dir.create(beagle_outputs_path)
  dir.create(glue("{beagle_outputs_path}/original"))
  dir.create(glue("{beagle_outputs_path}/masked"))
  
  snp_array_name = "Factorial_50K_SNP_Array_Genotypes.csv"
  cam_map_name = "camMap.csv"
  
  snp_array = read.csv(glue("{provided_data_path}/{snp_array_name}"))
  cam_map = read.csv(glue("{provided_data_path}/{cam_map_name}"))
  
  return(list(snp = snp_array, cam = cam_map))
}

map_safe_organism_names <- function(snp_array, beagle_inputs_path, name_mapper_name) {
  org_names = unique(tail(colnames(snp_array), -5))
  safe_names = c()
  count = 0
  for (i in org_names) {
    safe_names = c(glue("Organism{count}"), safe_names)
    count = count + 1
  }
  name_mapper = data.frame(org_names, safe_names)
  write.csv(name_mapper, glue("{beagle_inputs_path}/{name_mapper_name}"))
}

clean_snp_array <- function(snp_array, cam_map) {
  # add the cM position to each row of the snp array table
  # joined_snp_cam = merge(snp_array, cam_map, by.x = "Probe_Set_ID", by.y = "marker", all.x = TRUE)
  joined_snp_cam = cbind(snp_array)
  
  # replace the org names in the cleaned_mm table with the safe names
  c <- colnames(joined_snp_cam)
  for (i in seq(1, nrow(name_mapper))) {
    curr_org_name <- name_mapper[i, "org_names"]
    curr_safe_name <- name_mapper[i, "safe_names"]
    c[c == curr_org_name] <- curr_safe_name
  }
  colnames(joined_snp_cam) <- c
  
  # remove rows where the position or cM is NA
  filtered_joined_snp_cam = joined_snp_cam %>% filter(!is.na(Position))
  
  necessary_cols = "(Organism*)|(Probe_Set_ID)|(Chromosome)|(Position)|(REF$)|(ALT)|(cM)"
  filtered_joined_snp_cam = as.data.frame(filtered_joined_snp_cam %>% select(matches(necessary_cols)))
  grouped_joined_snp_cam = filtered_joined_snp_cam %>% group_by(Chromosome) %>% arrange(Position, .by_group = TRUE)
  filtered_joined_snp_cam = grouped_joined_snp_cam
  
  filtered_joined_snp_cam = filtered_joined_snp_cam[(!(duplicated(filtered_joined_snp_cam$Chromosome) & duplicated(filtered_joined_snp_cam$Position))),]
  
  return(filtered_joined_snp_cam)
}

find_nas_snp_array <- function(filtered_joined_snp_cam, beagle_inputs_path) {
  # where are all of the NAs in the original dataset?
  # replace cells: 0 => 0/0, 1 => 0/1, 2 => 1/1
  filtered_joined_snp_cam_copy <- cbind(filtered_joined_snp_cam)
  filtered_joined_snp_cam_copy[] <- lapply(filtered_joined_snp_cam_copy, gsub, pattern = " ", replacement = "")
  filtered_joined_snp_cam_copy[filtered_joined_snp_cam_copy == "0"] <- "0/0"
  filtered_joined_snp_cam_copy[filtered_joined_snp_cam_copy == "1"] <- "0/1"
  filtered_joined_snp_cam_copy[filtered_joined_snp_cam_copy == "2"] <- "1/1"
  
  filtered_joined_snp_cam_copy_colnames <- colnames(filtered_joined_snp_cam_copy)
  filtered_joined_snp_cam_copy <- filtered_joined_snp_cam_copy[, filtered_joined_snp_cam_copy_colnames[filtered_joined_snp_cam_copy_colnames %in% name_mapper$safe_names]]
  filtered_joined_snp_cam_copy <- filtered_joined_snp_cam_copy[order(names(filtered_joined_snp_cam_copy))]
  na_locs <- data.frame(which(is.na(filtered_joined_snp_cam_copy), arr.ind = TRUE))
  write_csv(na_locs, glue("{beagle_inputs_path}/na_locs.csv"))
  return(filtered_joined_snp_cam_copy)
}

swap_na_patterns <- function(beagle_inputs_path, na_locs) {
  num_orgs = length(unique(na_locs[["col"]]))
  mod_na_locs <- na_locs
  for (i in seq(0, num_orgs - 1)) {
    org_1_idx = i + 1
    org_2_idx = ((i + 1) %% num_orgs) + 1

    mod_na_locs$col[na_locs$col == org_1_idx] <- org_2_idx
  }
  return(mod_na_locs)
}

generate_beagle_files <- function(beagle_inputs_path, beagle_outputs_path, genotype_data, filtered_joined_snp_cam, name_mapper, na_locs, run_num, num_masked = 0, original = FALSE, all_masked = FALSE) {
  subfolder = "masked"
  if (original) {
    subfolder = "original"
  }
  if (all_masked) {
    subfolder = "all_masked"
  }
  
  dir.create(glue("{beagle_inputs_path}/{subfolder}/{run_num}_inputs/"))
  dir.create(glue("{beagle_outputs_path}/{subfolder}/{run_num}_outputs/"))
  
  beagle_df_copy <- cbind(filtered_joined_snp_cam)
  
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
  
  fixed_data <- data.table(
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
  
  meta_data = c(
    "##fileformat=VCFv4.1",
    "##source=StrawberryBreeding",
    "##phasing=partial",
    "##FILTER=<ID=PASS,Description=PASS>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
  )
  
  if (original || all_masked) {
    for (row in 1:nrow(na_locs)) {
      na_row <- na_locs[row, "row"]
      na_col <- na_locs[row, "col"]

      filtered_genotype_data[na_row, na_col] <- NA
    }
  } else {
    for (i in seq(1, length(unique(na_locs$col)))) {
      curr_na_locs = na_locs[na_locs$col == i, ]
      
      if (num_masked > nrow(curr_na_locs)) {
        random_na_locs = curr_na_locs[nrow(curr_na_locs), ]
      } else {
        random_na_locs = curr_na_locs[sample(nrow(curr_na_locs), num_masked), ]
      }
      
      for (j in seq(1, nrow(random_na_locs))) {
        curr_na_row = random_na_locs[j, "row"]
        curr_na_col = random_na_locs[j, "col"]
        
        filtered_genotype_data[curr_na_row, curr_na_col] <- NA
      }
    }
  }
  
  fixed_data <- as.matrix(fixed_data, row.names = NULL)
  genotype_data <- as.matrix(filtered_genotype_data, row.names = NULL)
  vcf_data <- new("vcfR", meta = meta_data, fix = fixed_data, gt = genotype_data)
  write.vcf(vcf_data, file = glue("{beagle_inputs_path}/{subfolder}/{run_num}_inputs/genotype.vcf.gz"))
}

phased_to_unphased <- function() {
  # read output VCF file
  phased_genotypes = read.vcfR(glue("{beagle_outputs_path}/original/1_outputs/phased_genotypes.vcf"))
  phased_genotype_data = phased_genotypes@gt
  fixed_data = as.data.frame(phased_genotypes@fix)
  phased_genotype_data[phased_genotype_data == "0|0"] <- "0/0"
  phased_genotype_data[phased_genotype_data == "0|1" | phased_genotype_data == "1|0"] <- "0/1"
  phased_genotype_data[phased_genotype_data == "1|1"] <- "1/1"
  
  print(glue("There are no NAs in the phased 'source of truth' -- {sum(is.na(phased_genotype_data))}"))
  
  phased_genotype_data_df <- as.data.frame(phased_genotype_data)
  phased_genotype_data_df$Chromosome = fixed_data$CHROM
  phased_genotype_data_df$Position = fixed_data$POS
  phased_genotype_data_df$REF = fixed_data$REF
  phased_genotype_data_df$ALT = fixed_data$ALT
  
  return(list(gt = phased_genotype_data, all = phased_genotype_data_df))
}

# RUN CODE
root_path = "~/Documents/beagle/"
provided_data_folder = "new_provided_data"
inputs_folder = "shuffle_nans_inputs"
outputs_folder = "shuffle_nans_outputs"

root_path = "~/Documents/beagle/"
provided_data_path = glue(root_path, provided_data_folder)
beagle_inputs_path = glue(root_path, inputs_folder)
beagle_outputs_path = glue(root_path, outputs_folder)

loaded = load_data(root_path, beagle_inputs_path, beagle_outputs_path)

snp_array = loaded$snp
cam_map = loaded$cam

name_mapper_name = "organism_safe_name_mapper.csv"
map_safe_organism_names(snp_array, beagle_inputs_path, name_mapper_name)

filtered_joined_snp_cam = clean_snp_array(snp_array, cam_map)
View(filtered_joined_snp_cam)

genotype_data = find_nas_snp_array(filtered_joined_snp_cam, beagle_inputs_path)

na_locs <- read.csv(glue("{beagle_inputs_path}/na_locs.csv"))
mod_na_locs <- swap_na_patterns(beagle_inputs_path,na_locs)
write.csv(mod_na_locs, glue("{beagle_inputs_path}/mod_na_locs.csv"))

# generate source of truth
name_mapper = read.csv(glue("{beagle_inputs_path}/{name_mapper_name}"))
generate_beagle_files(beagle_inputs_path, beagle_outputs_path, genotype_data, filtered_joined_snp_cam, name_mapper, na_locs, 1, original = TRUE)

# run beagle on source of truth

phased_results <- phased_to_unphased()
original_genotype <- phased_results$gt
original_unphased <- phased_results$all
View(original_unphased)

# all masked in modified locs
generate_beagle_files(beagle_inputs_path, beagle_outputs_path, original_genotype, original_unphased, name_mapper, mod_na_locs, 1, all_masked = TRUE)

num_masked_args = c(500, 300, 200, 100, 80, 50, 30, 20, 10, 5, 2)

print("Generating masked location beagle input files")

for (i in seq(1, length(num_masked_args))) {
  generate_beagle_files(beagle_inputs_path, beagle_outputs_path, original_genotype, original_unphased, name_mapper, mod_na_locs, i, num_masked = num_masked_args[i])
  print(glue("Generated files for {num_masked_args[i]}"))
}

# run beagle on farm on 