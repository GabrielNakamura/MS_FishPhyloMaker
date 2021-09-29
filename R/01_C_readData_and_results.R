# read data ---------------------------------------------------------------

drainage_basins <- read.csv(here::here("data", "Drainage_Basins_Table.csv"), sep = ";")
occ_drainage <- read.csv(here::here("data", "Occurrence_Table.csv"), sep = ";")
occ_drainage$X2.Species.Name.in.Source <- gsub("[.]", "_", occ_drainage$X2.Species.Name.in.Source)
occ_drainage$X6.Fishbase.Valid.Species.Name <- gsub("[.]", "_", occ_drainage$X6.Fishbase.Valid.Species.Name)

# read results ------------------------------------------------------------

res_afrotropics <- readRDS(here::here("output", "phylo_afrotropics.rds"))
res_indomalay <- readRDS(here::here("output", "phylo_indomalay.rds"))
res_neartic <- readRDS(here::here("output", "phylo_neartic.rds"))
res_neotropic <- readRDS(here::here("output", "phylo_neotropic.rds"))

