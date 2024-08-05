

## Sample names with more than one strain
sample_name_nms <- c("MSC27", "MSC35", "MSC42", "MSC44", "MSC47", "MSC36", "MSC16", "MSC15") %>% sort()
# sample_name_nms <- c("MSC35",  "MSC47", "MSC36", "MSC16", "MSC15") %>% sort()
# sample_name_nms <- c("MSC47","MSC16") %>% sort()


# Get the optimum QCd K
opt_narrow_clusters_nms = list("MSC27" = 2, "MSC35" = c(2:4), "MSC42" = 2, "MSC44" = c(2:4), "MSC47" = c(5:7), "MSC36" = c(2:3), "MSC16" = c(2:3), "MSC15" = c(2:4)) %>% .[sample_name_nms]
# optimum_k = opt_narrow_clusters_nms

sample_name_k <- map2(sample_name_nms,
                      opt_narrow_clusters_nms,
                      ~ rep(.x, each = length(.y))) %>% unlist()

## sorted donor_opt_cluster names
sorted_names <- imap(opt_narrow_clusters_nms, ~ paste0(.y, "_", .x)) %>% unlist(., recursive =
                                                                                  F, use.names = F)

optimum_k <- opt_narrow_clusters_nms %>% unlist() %>% set_names(sorted_names)

# Get the pseudobulk groups/categories
grps <- c('strain')

##Alignment incorporation
algn_nms = c('minmap')
algn_nm = 'minmap'

gtyp_step = "pbulk_gtypes_preQC"
# gtyp_step = "pbulk_gtypes_GE_postQC"

sample_set <- "Mali2"
species <- "Pm"

mt_chrom_nm <- "PmUG01_MIT_v1"

## NB-DONT DELETE BELOW- RUN ONCE Get genome size like this cut -f1,2 fasta.fai > genome.size
# system("cut -f1,2 /lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/PlasmoDB-66_PmalariaeUG01_Genome.fasta.fai > /lustre/scratch126/tol/teams/lawniczak/users/jr35/Pm_genomes_gtfs_indexs/Pdb66_Genome.size")
chrm_sz <- "../../Pm_genomes_gtfs_indexs/Pdb66_Genome.size"
chrm_prfx <- "PmUG01_0|PmUG01_|_v1"