devtools::install_github("GabrielNakamura/FishPhyloMaker", ref = "main", force = TRUE)
library(FishPhyloMaker)
data("fish_SAmerica")
taxn_data <- FishPhyloMaker::FishTaxaMaker(data = fish_SAmerica)
Loricariidae
Siluriformes
Anostomidae 
Characiformes
Anostomidae 
Characiformes
Anostomidae 
Characiformes
Serrasalmidae
Characiformes
Auchenipteridae
Siluriformes
Auchenipteridae
Siluriformes
Auchenipteridae
Siluriformes
insertions <- whichFishAdd(data = taxn_data)
phylo_fish_SAmerica <- FishPhyloMaker(data = taxn_data, return.insertions = T)
Serrasalmidae
Acanthicus
Schizodon
Liosomadoras
Hemiancistrus
Ctenolucius
Roeboexodon Exodon
Leporinus
Gymnotus
Loricariichthys
Anostomus
Furcodontichthys
Lasiancistrus
Leporinus
Schizodon Hypomasticus
Batrochoglanis
Pseudolithoxus
Leporinus
Anodus
Paralonchurus
Pachypops
Vandellia Paravandellia
Pachyurus
Anostomus
Sternarchorhamphus
Vandellia
Odontognathus
Exallodontus
Myleus
Gnathodolus
Trachelyopterus
Hemiodontichthys
Rhinelepis
Pseudohemiodon
Schizodon
Gnathodolus
Panaque Panaqolus
Pseudacanthicus Leporacanthicus
Hypopomidae
Sternarchogiton
Platyurosternarchus
Gymnocorymbus
Schizodon
Synbranchidae
Ageneiosus
Asterophysus
Trachelyopterus
Ageneiosus

dim(taxon_data)
