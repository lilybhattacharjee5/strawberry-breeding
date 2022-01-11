# generate masked files for BEAGLE

# library imports
library(vcfR)
library(data.table)
library(readr)
library(dplyr)
library(glue)
library(lattice)
library(ggplot2)

num_linkage_groups = 28

# file / path names
root_path = "~/Documents/beagle/"
provided_data_path = glue(root_path, "new_provided_data")
beagle_inputs_path = glue(root_path, "new_inputs")
beagle_outputs_path = glue(root_path, "new_outputs")

dir.create(beagle_inputs_path)
dir.create(beagle_outputs_path)

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
write.csv(name_mapper, glue("{root_path}/{name_mapper_name}"))

# add the cM position to each row of the snp array table
joined_snp_cam = merge(snp_array, cam_map, by.x = "Probe_Set_ID", by.y = "marker")

# replace the org names in the cleaned_mm table with the safe names
c <- colnames(joined_snp_cam)
for (i in seq(1, nrow(name_mapper))) {
  curr_org_name <- name_mapper[i, "org_names"]
  curr_safe_name <- name_mapper[i, "safe_names"]
  c[c == curr_org_name] <- curr_safe_name
}
colnames(joined_snp_cam) <- c

# remove rows where the position or cM is NA
nrow(joined_snp_cam)
nrow(snp_array)
filtered_joined_snp_cam = joined_snp_cam %>% filter(!is.na(cM)) %>% filter(!is.na(Position))

View(filtered_joined_snp_cam)

necessary_cols = "(Organism*)|(Probe_Set_ID)|(Chromosome.x)|(Position)|(REF$)|(ALT)|(cM)"
filtered_joined_snp_cam = filtered_joined_snp_cam %>% select(matches(necessary_cols))
View(filtered_joined_snp_cam)

# create genotypes VCF file
chrom_linkage_groups <- filtered_joined_snp_cam$Chromosome.x
pos_vals <- filtered_joined_snp_cam$Position
id_vals <- rep(NA, length(pos_vals))
ref_vals <- filtered_joined_snp_cam$REF
alt_vals <- filtered_joined_snp_cam$ALT
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

genotypeData <- as.matrix(fixedData, row.names = NULL)

# chromosome | rs# / snp identifier | genetic distance (morgans) | base-pair position
# Chromosome | SNP | cM / 100 | RR_Pos
plink_map <- data.frame(
  Chromosome = filtered_joined_snp_cam$Chromosome.x,
  SNP_ID = filtered_joined_snp_cam$Probe_Set_ID,
  Morgans = filtered_joined_snp_cam$cM / 100,
  Pos = filtered_joined_snp_cam$Position
)

# sort by ascending cM order
plink_map <- plink_map[order(plink_map$Morgans), ]

# modify the pos values to be consistent with morgans ordering
new_pos_vals <- plink_map$Pos
old_pos_vals <- pos_vals

for (j in seq(2, length(new_pos_vals))) {
  curr_val <- new_pos_vals[j]
  prev_val <- new_pos_vals[j - 1]
  if (prev_val >= curr_val) {
    new_pos_vals[j] <- prev_val + 1
  }
}
new_pos_matcher <- data.frame(old_pos_vals, new_pos_vals)

plink_map$Pos <- new_pos_vals

# create map file
write.table(plink_map, glue("{beagle_inputs_path}/strawberry.map"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

fixed_pos_vals <- fixedData$POS
for (j in seq(1, length(new_pos_vals))) {
  curr_fp_val <- fixed_pos_vals[j]

  # find curr_fp_val in new_pos_matcher
  match <- new_pos_matcher[new_pos_matcher$old_pos_vals == curr_fp_val, ][1, "new_pos_vals"]
  fixed_pos_vals[j] <- match
}
fixedData$POS <- fixed_pos_vals

# replace cells: 0 => 0/0, 1 => 0/1, 2 => 1/1
filtered_joined_snp_cam[] <- lapply(filtered_joined_snp_cam, gsub, pattern = " ", replacement = "")
filtered_joined_snp_cam[filtered_joined_snp_cam == "0"] <- "0/0"
filtered_joined_snp_cam[filtered_joined_snp_cam == "1"] <- "0/1"
filtered_joined_snp_cam[filtered_joined_snp_cam == "2"] <- "1/1"
View(filtered_joined_snp_cam)

filtered_joined_snp_cam_colnames <- colnames(filtered_joined_snp_cam)
filtered_genotype_data <- filtered_joined_snp_cam[, filtered_joined_snp_cam_colnames[filtered_joined_snp_cam_colnames %in% name_mapper$safe_names]]
filtered_genotype_data$FORMAT <- rep("GT", length(pos_vals))
filtered_genotype_data <- filtered_genotype_data[, c(length(colnames(filtered_genotype_data)), 1:length(colnames(filtered_genotype_data)) - 1)]

genotypeData <- as.matrix(filtered_genotype_data)
fixedData <- as.matrix(fixedData)

View(genotypeData)
View(fixedData)

vcfData <- new("vcfR", meta = metaData, fix = fixedData, gt = genotypeData)
write.vcf(vcfData, file = glue("{beagle_inputs_path}/unphased_genotypes.vcf.gz"))

# 
# generateMaskedFiles <- function(runNum, numMasked) {
#   # split genotypes by linkage group
#   chrGroupedAlleleFreqsData <- alleleFreqsData %>% group_by(Chr)
#   dir.create(glue("{beagle_inputs_path}/masking/{runNum}_inputs"), showWarnings = TRUE)
#   for (i in group_split(chrGroupedAlleleFreqsData)) {
#     
#     # write i to a separate file
#     currLinkageGroup = as.character(i[1, "Chr"])
#     dir.create(glue("{beagle_inputs_path}/masking/{runNum}_inputs/{currLinkageGroup}_inputs/"), showWarnings = TRUE)
#   }
#   
#   # create separate vcf files per "Group"
#   chrGroupedJoinedVcfData <- joined_vcf_table %>% group_by(Group)
#   for (i in group_split(chrGroupedJoinedVcfData)) {
#     currLinkageGroup = as.character(i[1, "Group"])
#     
#     # regenerate the vcf file format
#     chrom_linkage_groups <- i$Group
#     pos_vals <- i$RR_Pos
#     id_vals <- rep(NA, length(pos_vals))
#     ref_vals <- rep("A", length(pos_vals))
#     alt_vals <- rep("G", length(pos_vals))
#     qual_vals <- rep(99, length(pos_vals))
#     filter_vals <- rep(NA, length(pos_vals))
#     info_vals <- rep(NA, length(pos_vals))
#     
#     fixedData <- data.table(
#       CHROM = chrom_linkage_groups,
#       POS = pos_vals,
#       ID = id_vals,
#       REF = ref_vals,
#       ALT = alt_vals,
#       QUAL = qual_vals,
#       FILTER = filter_vals,
#       INFO = info_vals
#     )
#     
#     cleaned_mm_colnames <- colnames(joined_vcf_table)
#     filtered_genotype_data <- i[, cleaned_mm_colnames[cleaned_mm_colnames %in% name_mapper$safe_names]]
#     
#     filtered_genotype_data$FORMAT <- rep("GT", length(pos_vals))
#     filtered_genotype_data <- filtered_genotype_data[, c(length(colnames(filtered_genotype_data)), 1:length(colnames(filtered_genotype_data)) - 1)]
#     
#     # masking
#     hetLocs <- which(filtered_genotype_data == "0/1", arr.ind = TRUE)
#     # randomly sample n hetLocs
#     hetLocsSelected <- hetLocs[sample(nrow(hetLocs), numMasked), ]
#     filtered_genotype_data_copy <- cbind(filtered_genotype_data)
#     
#     for (h in seq(1, nrow(hetLocsSelected))) {
#       currHetLoc <- hetLocsSelected[h, ]
#       filtered_genotype_data_copy[currHetLoc["row"], currHetLoc["col"]] = NA
#     }
#     # write masked hetLocs to file
#     write.csv(hetLocsSelected, glue("{beagle_inputs_path}/masking/{runNum}_inputs/{currLinkageGroup}_inputs/masked_locations.csv"))
#     
#     genotypeData <- as.matrix(filtered_genotype_data_copy)
#     
#     # chromosome | rs# / snp identifier | genetic distance (morgans) | base-pair position
#     # Chromosome | Group | cM / 100 | RR_Pos
#     plink_map <- data.frame(
#       Chromosome = i$Group,
#       Group = i$Group, 
#       Morgans = i$cM / 100, 
#       Pos = i$RR_Pos
#     )
#     
#     # sort to put map positions in ascending order
#     plink_map <- plink_map[order(plink_map$Morgans), ]
#     
#     # modify the pos vals to be consistent with morgans ordering
#     new_pos_vals <- plink_map$Pos
#     old_pos_vals <- pos_vals
#     
#     for (j in seq(2, length(new_pos_vals))) {
#       curr_val <- new_pos_vals[j]
#       prev_val <- new_pos_vals[j - 1]
#       if (prev_val >= curr_val) {
#         new_pos_vals[j] <- prev_val + 1
#       }
#     }
#     new_pos_matcher <- data.frame(old_pos_vals, new_pos_vals)
#     
#     plink_map$Pos <- new_pos_vals
#     
#     write.table(plink_map, glue("{beagle_inputs_path}/masking/{runNum}_inputs/{currLinkageGroup}_inputs/strawberry.map"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
#     
#     fixed_pos_vals <- fixedData$POS
#     for (j in seq(1, length(new_pos_vals))) {
#       curr_fp_val <- fixed_pos_vals[j]
#       
#       # find curr_fp_val in new_pos_matcher
#       match <- new_pos_matcher[new_pos_matcher$old_pos_vals == curr_fp_val, ][1, "new_pos_vals"]
#       fixed_pos_vals[j] <- match
#     }
#     fixedData$POS <- fixed_pos_vals
#     
#     fixedData <- as.matrix(fixedData, row.names = NULL)
#     vcfData <- new("vcfR", meta = metaData, fix = fixedData, gt = genotypeData)
#     write.vcf(vcfData, file = glue("{beagle_inputs_path}/masking/{runNum}_inputs/{currLinkageGroup}_inputs/genotype_{currLinkageGroup}.vcf.gz"))
#   }
# }
# 
# numMaskedArgs <- c(2, 3, 5, 10, 15, 20, 50, 100, 200) # c(100, 200, 500, 1000, 2000)
# for (i in seq(1, length(numMaskedArgs))) {
#   generateMaskedFiles(i, numMaskedArgs[i])
# }

##### CHECK BEAGLE OUTPUTS

ratio <- c()

# combine the output vcf files
for (run_i in seq(1, length(numMaskedArgs))) {
  full_fixed_data <- data.table(
    CHROM = c(),
    POS = c(),
    ID = c(),
    REF = c(),
    ALT = c(),
    QUAL = c(),
    FILTER = c(),
    INFO = c()
  )
  
  full_masked <- data.table(
    row = c(),
    col = c()
  )
  
  het_phased_count <- 0
  total_het_count <- 0
  
  input_dirs <- list.dirs(path = glue("{beagle_inputs_path}/masking/{run_i}_inputs/"), full.names = TRUE)
  print(length(input_dirs))
  
  for (i in seq(1, num_linkage_groups)) {
    output_path <- glue("{beagle_outputs_path}/masking/{run_i}_outputs/{i}_phased_genotypes.vcf")
    input_path <- glue("{input_dirs[i + 1]}/cleaned_genotypes.vcf")
    het_locs_path <- glue("{input_dirs[i + 1]}/masked_locations.csv")
    
    curr_output <- read.vcfR(output_path)
    curr_input <- read.vcfR(input_path)
    
    curr_fixed_data <- as.data.table(curr_output@fix)
    curr_unphased_genotype_data <- as.data.table(curr_input@gt)
    curr_genotype_data <- as.data.table(curr_output@gt)
    curr_het_locs_data <- as.data.frame(read.csv(het_locs_path))
    
    for (j in nrow(curr_het_locs_data)) {
      idx <- curr_het_locs_data[j, ]
      curr_out <- as.vector(unname(curr_genotype_data[idx$row]))
      curr_in <- as.vector(unname(curr_unphased_genotype_data[idx$row]))
      
      curr_in_unphased <- curr_out[[idx$col]]
      curr_out_phased <- curr_out[[idx$col]]
      
      print(curr_out_phased)
      print(curr_in_unphased)
      if (!is.na(curr_out_phased)) {
        if (curr_out_phased == "0|1" || curr_out_phased == "1|0") {
          het_phased_count = het_phased_count + 1;
        }
      }
    }
    
    total_het_count = total_het_count + nrow(curr_het_locs_data)
    
    full_fixed_data <- rbindlist(list(full_fixed_data, curr_fixed_data))
    full_masked <- rbindlist(list(full_masked, curr_het_locs_data))
    if (i == 1) {
      full_genotype_data <- curr_genotype_data
    } else {
      full_genotype_data <- rbindlist(list(full_genotype_data, curr_genotype_data))
    }
  }
  
  ratio <- c(ratio, het_phased_count / total_het_count)
  print(glue("correctly phased heterozygous ratio {het_phased_count / total_het_count}"))
}

print(ratio)

# graph correctly phased heterozygous ratio
heterozygous_phased_ratio <- data.frame(num_masked = numMaskedArgs * num_linkage_groups, correctly_phased_ratio = ratio)
head(heterozygous_phased_ratio)
ggplot(data = heterozygous_phased_ratio, aes(x = num_masked, y = correctly_phased_ratio)) + geom_line() + geom_point() + ggtitle("Correctly Phased Heterozygous Ratio vs. Number of Masked Unphased Genotypes") + xlab("Number of Masked Unphased Genotypes") + ylab("Correctly Phased Heterozygous Ratio")


