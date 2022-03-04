# simulate fake pedigree

# library imports
library(vcfR)
library(data.table)
library(readr)
library(dplyr)
library(glue)
library(lattice)
library(ggplot2)
library(stringr)
library(MoBPS)

set.seed(5)

# file / path names
root_path = "~/Documents/beagle/"
provided_data_path = glue(root_path, "new_provided_data")
beagle_inputs_path = glue(root_path, "simulated_inputs_rebred")
beagle_outputs_path = glue(root_path, "simulated_outputs_rebred")

dir.create(beagle_inputs_path)
dir.create(beagle_outputs_path)

# get the markers 
snp_array_name = "Factorial_50K_SNP_Array_Genotypes.csv"
snp_array = read.csv(glue("{provided_data_path}/{snp_array_name}"))

# grab (Probe_Set_ID, Chromosome, Position)
marker_data = snp_array[c("Probe_Set_ID", "Chromosome", "Position")]

# # how many founders have no NA SNPs? 0
# # which ones have the minimum number of NA SNPs?
# # grab 5 founders from the ordered set according to the number of NAs in the corresponding column
# genotype_data = snp_array[, seq(6, ncol(snp_array))]
# chosen_founders = names(sort(colSums(is.na(genotype_data)))[1:5])
# print(chosen_founders)
# chosen_founder_data = snp_array[, chosen_founders]
# View(chosen_founder_data)
# 
# # replace the 2's => 1|1, 1's => 0|1 or 1|0, 0's => 0|0
# chosen_founder_data[chosen_founder_data == 2] <- "1|1"
# chosen_founder_data[chosen_founder_data == 0] <- "0|0"
# chosen_founder_data[chosen_founder_data == 1] <- sample(c("0|1", "1|0"), 1)
# 
# # replace the Organism names with safe names
# counter = 0
# for (c in colnames(chosen_founder_data)) {
#   chosen_founder_data = rename(chosen_founder_data, !!glue("Organism_{counter}") := all_of(c))
#   counter = counter + 1
# }

# View(chosen_founder_data)
# 
# # use mobps to breed the founders over several generations


# # create an initial set of randomized founder organisms
# new_marker_data <- marker_data
# num_founders = 5
# counter = 0
# phased_genotypes = c("0|0", "0|1", "1|1", "1|0")
# for (i in seq(1, num_founders)) {
#   new_founder_data = sample(phased_genotypes, nrow(marker_data), replace = TRUE)
#   new_organism_name = glue("Organism_{counter}")
#   new_marker_data$new_organism = new_founder_data
#   new_marker_data = rename(new_marker_data, !!new_organism_name := "new_organism")
#   counter = counter + 1
# }

strsplit_pipe <- function(x) {
  return(strsplit(x, "\\|"))
}

# create_child <- function(x) {
#   return(sample(x, 1, replace = TRUE)[1])
# }
#  
# reproduce <- function(p1, p2, num_children, child_idx, marker_data, pedigree) {
#   p1_name = glue("Organism_{p1}")
#   p2_name = glue("Organism_{p2}")
#   
#   p1_data = marker_data[[p1_name]]
#   p2_data = marker_data[[p2_name]]
#   
#   p1_vals = sapply(p1_data, strsplit_pipe)
#   p2_vals = sapply(p2_data, strsplit_pipe)
#   
#   p1_child_vals = sapply(p1_vals, create_child)
#   p2_child_vals = sapply(p2_vals, create_child)
#   
#   for (i in seq(1, num_children)) {
#     p1_p2_child = paste0(p1_child_vals, "|", p2_child_vals)
#     new_organism_name = glue("Organism_{child_idx}")
#     marker_data$new_organism = p1_p2_child
#     marker_data = rename(marker_data, !!new_organism_name := "new_organism")
#     child_idx = child_idx + 1
#     pedigree[nrow(pedigree) + 1, ] = c(new_organism_name, p1_name, p2_name)
#   }
#   
#   return(list(marker_data, pedigree))
# }

# iterate through several generations, creating offspring from every selected pairing
num_generations = 2
num_parents = num_founders
num_total = num_founders
max_children = 4
pairings_per_generation = 0.8
curr_parents = seq(1, num_founders)

pedigree <- data.frame(
  organism = character(0),
  p1 = character(0),
  p2 = character(0)
)

for (i in seq(1, num_generations)) {
  num_pairings = floor(pairings_per_generation * num_parents * (num_parents - 1) / 2)
  possible_pairings = combn(1:num_parents, 2, simplify = FALSE)
  selected_pairings = sample(possible_pairings, num_pairings, replace = FALSE)
  num_children = 0
  generation_children = c()
  for (pairing in selected_pairings) {
    curr_p1 = pairing[1]
    curr_p2 = pairing[2]
    curr_num_children = sample(1:max_children, 1, replace = FALSE)[1]
    reproduction_results = reproduce(curr_p1, curr_p2, curr_num_children, num_total, new_marker_data, pedigree)
    new_marker_data = reproduction_results[[1]]
    pedigree = reproduction_results[[2]]
    generation_children = c(generation_children, seq(num_children, num_children + curr_num_children))
    num_children = num_children + curr_num_children
    num_total = num_total + curr_num_children
  }
  num_parents = num_children
  curr_parents = generation_children
}

View(new_marker_data)
View(pedigree)

write.csv(pedigree, glue("{beagle_inputs_path}/pedigree.csv"))

# randomly mask (i.e. set to NA) X genotypes per organism
num_masked = c(0, 2, 5, 10, 20, 50, 100) #, 500, 1000, 2000, 5000)

# create genotype VCF data structure
new_marker_data = new_marker_data %>% group_by(Chromosome) %>% arrange(Position, .by_group = TRUE)
metadata <- c(
  "##fileformat=VCFv4.1",
  "##source=StrawberryBreeding",
  "##phasing=partial",
  "##FILTER=<ID=PASS,Description=PASS>",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)
chrom_linkage_groups <- new_marker_data$Chromosome
pos_vals <- new_marker_data$Position
id_vals <- rep(NA, length(pos_vals))
ref_vals <- rep("A", length(pos_vals))
alt_vals <- rep("G", length(pos_vals))
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
fixed_data <- as.matrix(fixed_data, row.names = NULL)
genotype_data <- as.matrix(new_marker_data[, seq(4, ncol(new_marker_data))], row.names = NULL)
vcf_data <- new("vcfR", meta = metadata, fix = fixed_data, gt = genotype_data)
write.vcf(vcf_data, file = glue("{beagle_inputs_path}/phased_genotype.vcf.gz"))
View(new_marker_data)

new_marker_unphased <- data.frame(lapply(new_marker_data, function(x) {
  gsub("\\|", "/", x)
}))
new_marker_unphased <- data.frame(lapply(new_marker_unphased, function(x) {
  gsub("1/0", "0/1", x)
}))

run_num = 0
for (n in num_masked) {
  # masking
  filtered_genotype_data <- new_marker_unphased[, seq(4, ncol(new_marker_unphased))]
  
  dir.create(glue("{beagle_inputs_path}/{run_num}_inputs"))
  
  if (n > 0) {
    # randomly sample n hetLocs (for each organism)
    all_het_locs <- data.frame(
      row = c(),
      col = c()
    )
    het_locs <- data.frame(which(filtered_genotype_data == "0/1", arr.ind = TRUE))
    
    for (i in seq(1, ncol(filtered_genotype_data))) {
      het_locs_curr = het_locs[het_locs$col == i, ]
      het_locs_selected = het_locs[sample(nrow(het_locs), n), ]
      
      all_het_locs = rbind(all_het_locs, het_locs_selected)
    }
    
    filtered_genotype_data_copy <- filtered_genotype_data
    
    for (h in seq(1, nrow(all_het_locs))) {
      curr_het_loc <- all_het_locs[h, ]
      curr_x = curr_het_loc$row
      curr_y = curr_het_loc$col
      filtered_genotype_data_copy[curr_x, curr_y] = NA
    }
    
    filtered_genotype_data = filtered_genotype_data_copy
    
    # write masked hetLocs to file
    write.csv(all_het_locs, glue("{beagle_inputs_path}/{run_num}_inputs/masked_locations.csv"))
  }
  genotype_data <- as.matrix(filtered_genotype_data, row.names = NULL)
  vcf_data <- new("vcfR", meta = metadata, fix = fixed_data, gt = genotype_data)
  
  # save genotype + pedigree data
  write.vcf(vcf_data, file = glue("{beagle_inputs_path}/{run_num}_inputs/genotype.vcf.gz"))
  run_num = run_num + 1
  print(glue("{n} masked / organism done"))
}

# check the robustness of beagle outputs
ratio <- c()
# combine the output vcf files
for (run_i in seq(0, length(num_masked))) { 
  het_phased_count <- 0
  total_het_count <- 0
  
  input_dirs <- list.dirs(path = glue("{beagle_inputs_path}/{run_i}_inputs/"), full.names = TRUE)
  
  output_path <- glue("{beagle_outputs_path}/{run_i}_outputs_noimpute/1_phased_genotypes.vcf")
  input_path <- glue("{beagle_inputs_path}/phased_genotype.vcf")
  het_locs_path <- glue("{beagle_inputs_path}/{run_i}_inputs/masked_locations.csv")
  
  curr_output <- read.vcfR(output_path)
  curr_input <- read.vcfR(input_path)
  
  curr_fixed_data <- as.data.table(curr_output@fix)
  correct_genotype_data <- as.data.table(curr_input@gt)
  curr_genotype_data <- as.data.table(curr_output@gt)
  
  #View(correct_genotype_data)
  #View(curr_genotype_data)

  tryCatch(
    expr = {
      curr_het_locs_data <- as.data.frame(read.csv(het_locs_path))
    },
    error = function(e){ 
      curr_het_locs_data = as.data.frame(data.table(
        row = c(),
        col = c()
      ))
    }
  )
  
  for (j in seq(1, nrow(curr_het_locs_data))) {
    idx <- curr_het_locs_data[j, ]
    curr_in <- as.vector(unname(correct_genotype_data[idx$row]))
    curr_out <- as.vector(unname(curr_genotype_data[idx$row]))
    
    curr_in_unphased <- curr_in[[idx$col]]
    curr_out_phased <- curr_out[[idx$col]]
    
    if (!is.na(curr_out_phased) & !is.na(curr_in_unphased)) {
      if (curr_out_phased == curr_in_unphased) {
        het_phased_count = het_phased_count + 1;
      }
    }
  }
  
  total_het_count = total_het_count + nrow(curr_het_locs_data)
  
  ratio <- c(ratio, het_phased_count / total_het_count)
  print(glue("input run {run_i} correctly phased heterozygous ratio {het_phased_count / total_het_count} total phased correctly {het_phased_count} total het {total_het_count}"))
}

print(ratio)

# graph correctly phased heterozygous ratio
heterozygous_phased_ratio <- data.frame(num_masked = num_masked, correctly_phased_ratio = ratio)
head(heterozygous_phased_ratio)
ggplot(data = heterozygous_phased_ratio, aes(x = num_masked, y = correctly_phased_ratio)) + geom_line() + geom_point() + ggtitle("Correctly Phased Heterozygous Ratio vs. Number of Masked Unphased Genotypes") + xlab("Number of Masked Unphased Genotypes") + ylab("Correctly Phased Heterozygous Ratio")
ggsave(glue("{beagle_outputs_path}/correctly_masked_ratio.png"))

(500 * ncol(correct_genotype_data)) / (nrow(correct_genotype_data) * ncol(correct_genotype_data))
