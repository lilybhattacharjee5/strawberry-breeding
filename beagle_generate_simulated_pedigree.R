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

# file / path names
root_path = "~/Documents/beagle/"
provided_data_path = glue(root_path, "new_provided_data")
beagle_inputs_path = glue(root_path, "simulated_inputs")
beagle_outputs_path = glue(root_path, "simulated_outputs")

dir.create(beagle_inputs_path)
dir.create(beagle_outputs_path)

# get the markers 
snp_array_name = "Factorial_50K_SNP_Array_Genotypes.csv"
name_mapper_name = "organism_safe_name_mapper_updated_data.csv"
cam_map_name = "camMap.csv"

snp_array = read.csv(glue("{provided_data_path}/{snp_array_name}"))
# cam_map = read.csv(glue("{provided_data_path}/{cam_map_name}"))

View(snp_array)

# grab (Probe_Set_ID, Chromosome, Position)
marker_data = snp_array[c("Probe_Set_ID", "Chromosome", "Position")]
View(marker_data)

# create an initial set of randomized founder organisms
new_marker_data <- marker_data
num_founders = 5
counter = 0
unphased_genotypes = c("0/0", "0/1", "1/1")
for (i in seq(1, num_founders)) {
  new_founder_data = sample(unphased_genotypes, nrow(marker_data), replace = TRUE)
  new_organism_name = glue("Organism_{counter}")
  new_marker_data$new_organism = new_founder_data
  new_marker_data = rename(new_marker_data, !!new_organism_name := "new_organism")
  counter = counter + 1
}

strsplit_slash <- function(x) {
  return(strsplit(x, "/"))
}

create_child <- function(x) {
  return(sample(x, 1, replace = TRUE)[1])
}
 
reproduce <- function(p1, p2, num_children, child_idx, marker_data, pedigree) {
  p1_name = glue("Organism_{p1}")
  p2_name = glue("Organism_{p2}")
  
  p1_data = marker_data[[p1_name]]
  p2_data = marker_data[[p2_name]]
  
  p1_vals = sapply(p1_data, strsplit_slash)
  p2_vals = sapply(p2_data, strsplit_slash)
  
  p1_child_vals = sapply(p1_vals, create_child)
  p2_child_vals = sapply(p2_vals, create_child)
  
  for (i in seq(1, num_children)) {
    p1_p2_child = paste0(p1_child_vals, "/", p2_child_vals)
    new_organism_name = glue("Organism_{child_idx}")
    marker_data$new_organism = p1_p2_child
    marker_data = rename(marker_data, !!new_organism_name := "new_organism")
    child_idx = child_idx + 1
    pedigree[nrow(pedigree) + 1, ] = c(new_organism_name, p1_name, p2_name)
  }
  
  return(list(marker_data, pedigree))
}

# iterate through several generations, creating offspring from every selected pairing
num_generations = 1
num_parents = num_founders
num_total = num_founders
max_children = 4
pairings_per_generation = 0.8

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
  for (pairing in selected_pairings) {
    curr_p1 = pairing[1]
    curr_p2 = pairing[2]
    curr_num_children = sample(1:max_children, 1, replace = FALSE)[1]
    reproduction_results = reproduce(curr_p1, curr_p2, curr_num_children, num_total, new_marker_data, pedigree)
    new_marker_data = reproduction_results[[1]]
    pedigree = reproduction_results[[2]]
    num_children = num_children + curr_num_children
    num_total = num_total + curr_num_children
  }
  num_parents = num_children
}

View(new_marker_data)
View(pedigree)
