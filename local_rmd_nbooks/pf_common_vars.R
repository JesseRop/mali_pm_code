library(MetBrewer)
library(RColorBrewer)

## Common tables and variables 

##Strain cluster colours
strain_cols_n <-c("#910000", "#ff1a8d", "#007c71", "#0ea300", "#00deca", "#9a96c7", "#5d57a4", "#a8875b","#708db3", "#c39a00", "#ffd94d", "#9770b3", "#C9E0F1", "#a40000" , "#00b7a7", "#ffcd12","thistle","#a40000" , "#00b7a7", "#e7f1f9", 'black', '#050301','grey')
#, 'lavender', 'lightgray', 'grey'
names(strain_cols_n) <- c('SC2', 'SC1','SC5',  'SC4', 'SC3',  'SC6', 'SC7', 'SC8', 'SC9', 'SC10', 'SC11', 'SC12', 'S3_b', 'S1', 'S2', 'Singlet', 'Negative', 'SC1_3i','Doublet', 'Lab', 'Less5', 'PoorQC', NA)

dblt_cols <- c(brewer.pal(name = 'Set1', n = 8)[c(1:5,7:8)], 'lightblue')
names(dblt_cols) <- c('ScDfSt', 'ScDf', 'ScSt', 'DfSt', 'Sc', 'Df','St', 'Singlet')

cell_qc_filt_nm_lvls <- c("sc_neg", "lo_gn_cnt",  "hi_gn_cnt",  "hi_mtp", "ss_db",  "strn_db", "stg_db", "pass")
cell_qc_cols <- c(brewer.pal(name = 'Set1', n = 8)[c(1:5,7:8)], 'lightblue') %>% set_names(cell_qc_filt_nm_lvls)

dblt_ss_cols <- c(brewer.pal(name = 'Set1', n = 8)[c(1:3)], 'lightgrey')
names(dblt_ss_cols) <- c('Strain & stage', 'Strain', 'Stage', 'Singlet')

cluster_cols <- brewer.pal(8, name = 'Set2')[c(1:4,7,8)] %>% set_names(as.character(c(0:5)))

stg_refnd_lvls = c("Early ring", "early ring", "Late ring", "late ring", "Early trophozoite", "early trophozoite", "Late trophozoite", "late trophozoite", "Early schizont", "early schizont",  "Late schizont","late schizont",  "Asexual", "Gametocyte (committed)", "Gametocyte (developing)", "gametocyte (developing)", "Gametocyte (branching)", "Female (HE)", "Female", "Female (LE)", "Female LE", "Gametocyte (female)", "gametocyte (female)", "Gametocyte (early female)", "Gametocyte (late female)", "Gametocyte (male)", "Male (HE)", "Male", "Male (LE)", "Male LE", "Gametocyte (early male)", "Gametocyte (late male)", "gametocyte (male)", "Sexual", "Unassigned", "Failed qc", "Not assigned", "Atlas", "stgQC_fail")

# stg_lvls = c("Early ring", "Late ring", "Early trophozoite","Late trophozoite" , "Early schizont", "Late schizont",  "Gametocyte (developing)", "Female (HE)", "Female", "Female (LE)",  "Gametocyte (female)", "Male (HE)", "Male", "Male (LE)", "Gametocyte (male)", "Asexual", "Unassigned", "Failed QC", "Not assigned")


# stg_refnd_lvls2 = c('Early ring','Late ring','Early trophozoite','Late trophozoite','Early schizont','Late schizont', "Gametocyte (committed)","Gametocyte (developing)","Gametocyte (branching)", "Female (HE)", "Female", "Female (LE)", "Gametocyte (early female)","Gametocyte (late female)", "Male (HE)", "Male", "Male (LE)", "Gametocyte (early male)","Gametocyte (late male)", "Asexual", "Unassigned", "Failed QC", "Not assigned") 

sp_colors <- c('Early ring'='#D1EC9F',
               'Late ring'='#78C679',
               'Early trophozoite'='#FEEEAA',
               'Late trophozoite'='#FEB24C',
               'Early schizont'='#C9E8F1',
               'Late schizont'='#85B1D3',
               'Gametocyte (developing)'='thistle',
               'Gametocyte (female)'='#551A8B',
               'Gametocyte (male)'='mediumpurple',
               'Unassigned'='red', 
               'Not assigned'='#D6CFC7',
               'Female (LE)'='#208cf7',
               'Female (HE)'='#551A8B',
               'Male (LE)'='#a894d1',
               'Male (HE)'='mediumpurple',
               'Failed QC' = 'grey',
               'Gametocyte (committed)' = '#C399A2', 
               'Gametocyte (developing)' = '#66444C', 
               'Gametocyte (branching)' = '#66444C', 
               'Gametocyte (early female)' = '#749E89', 
               'Gametocyte (early male)' = '#7D87B2', 
               'Gametocyte (late female)' = '#4E6D58', 
               'Gametocyte (late male)' = '#41507B'
)


sp_fld_colors <- c('Late ring'= '#78C679',#met.brewer('Java')[2],
                   'Early trophozoite'='#FEEEAA',#met.brewer('Java')[3],
                   'Female LE'= met.brewer('Ingres')[5],
                   'Female'= '#551A8B', #met.brewer('Ingres')[7],
                   # 'Female (HE)'= '#551A8B', #met.brewer('Ingres')[7],
                   # 'Male (LE)'=met.brewer('Nizami')[6],
                   # 'Male (HE)'=met.brewer('Nizami')[8],
                   # 'Gametocyte (male)'=met.brewer('Austria')[5],
                   'Male LE'=met.brewer('Ingres')[7],#'#a2e19e',
                   'Male'= 'mediumpurple',#'#039911',
                   ## 'Male (HE)'= 'mediumpurple',#'#039911',
                   'Gametocyte (male)'='#3DED97',
                   'Gametocyte (female)'=met.brewer('Ingres')[6],
                   'Lab' = '#D3D3D3',
                   'Unassigned'='#DCDCDC',
                   'Failed QC' = 'grey',
                   'Not assigned' = 'grey',
                   'Atlas' = 'grey',
                   'Early ring'='#D1EC9F',
                   'Late trophozoite'='#FEB24C',
                   'Early schizont'='#C9E8F1',
                   'Late schizont'='#85B1D3',
                   'Gametocyte (developing)'='thistle',
                   'early female' = '#749E89',
                   'late female'= '#4E6D58',
                   "Asexual" = "#FBC81D",
                   NULL = 'lightgrey',
                   'Field unlabelled'='#D2C1DC'
)

cluster_cols <- brewer.pal(7, name = 'Set2') %>% set_names(as.character(c(0:6)))


## donor colors 
donor_cols_m = c(met.brewer("Thomas")[c(1)],c('#F39EC8','#99ceff', '#b9e589')) %>% set_names(paste0(c("MSC1", "MSC3", "MSC13", "MSC14" ), '.mrna'))
# donor_cols = met.brewer("Egypt")[c(2,4,1,3)] %>% set_names(paste0(c("MSC1", "MSC3", "MSC13", "MSC14" )))
donor_cols = c(met.brewer("Thomas")[c(1)],c('#F39EC8','#99ceff', '#b9e589')) %>% set_names(paste0(c("MSC1", "MSC3", "MSC13", "MSC14" )))

c_stg_cols= c('Asexual' = "#D7B7C7", 'Male' = "#C4CAA4", 'Female'="#6289A9", 'Developing' = "#C9F6E1", 'Unassigned'="lightgrey")


##Metadata from all p.falciparum genes downloaded from plasmoDB
Pf3D7_genes_plasmoDB <- read_csv("/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs/Pf3d7_Genes_Summary_270922.csv") %>%
  dplyr::rename("gene_id" = `Gene ID` ) %>%
  group_by(gene_id) %>%
  add_count(name='gn_iso') %>%
  ungroup() %>%
  mutate(across(gene_id, ~case_when(gn_iso>1 ~ source_id,
                                    TRUE ~ as.character(.)))) %>%
  mutate(across(gene_id, ~str_replace(.,"_", "-")))%>%
  dplyr::rename('Gene_name' = `Gene Name or Symbol`, 'Ka_Ks' = `NonSyn/Syn SNP Ratio All Strains`, 'TM_domains'=`# TM Domains`) %>%
  mutate(across(where(is.character), ~na_if(., "N/A")),
         Gene = coalesce(`Gene_name`, gene_id) )  %>%
  group_by(Gene) %>%
  add_count() %>%
  ungroup() %>%
  dplyr::mutate(across(Gene, ~case_when(n>1 ~ paste0(.,'_', gene_id),
                                        TRUE ~ as.character(.))))


## Irods decode file plus metadata
irods_id_sk21_mdata_cln <- read_csv("/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/raw/irods_id_sk21_mdata_cln.csv") %>%
  mutate(sample_nm = str_remove_all(sample_name, "_JC")) %>%
  arrange(sample_nm)

smpl_nms <-irods_id_sk21_mdata_cln$sample_nm


