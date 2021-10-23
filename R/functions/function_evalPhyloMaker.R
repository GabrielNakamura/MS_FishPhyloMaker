eval_PhyloMaker <- function(k, 
                            tree_complete,
                            probs = 0.1,
                            insert.base.node = TRUE, 
                            return.insertions = TRUE, 
                            progress.bar = TRUE) 
{
  source(here::here("R","functions", "internal_filter_rank.R"))
  source(here::here("R", "functions", "internal_treedata_modif.R"))
  samp <- ceiling(length(tree_complete$tip.label) * probs)
  spp_samp <- tree_complete$tip.label[sample(1:length(tree_complete$tip.label),
                                             size = samp, replace = FALSE)]
  data_samp <- FishPhyloMaker::FishTaxaMaker(data = spp_samp, allow.manual.insert = FALSE)
  data_samp <- data_samp$All_info_fishbase[, c("user_spp", "Family", "Order")]
  colnames(data_samp) <- c("s", "f", "o")
  data_samp$o <- gsub("/.*", 
                      replacement = "", 
                      data_samp$o)
  spp_samp_return <- data_samp$s
  spp_samp_modif <- gsub("_.*", "_spp", data_samp$s)
  data_samp[, "s"] <- spp_samp_modif # species to be inserted
  spp_data <- 1:length(data_samp$s)
  names(spp_data) <- data_samp$s
  insert_spp <- treedata_modif(phy = tree_complete, data = spp_data, 
                               warnings = F)$nc$data_not_tree
  if (length(insert_spp) >= 1) {
    genus_in_tree <- sub("_.*", "", tree_complete$tip.label)[match(sub("_.*", 
                                                                     "", insert_spp), sub("_.*", "", tree_complete$tip.label))][!is.na(sub("_.*", 
                                                                                                                                         "", tree_complete$tip.label)[match(sub("_.*", "", 
                                                                                                                                                                              insert_spp), sub("_.*", "", tree_complete$tip.label))])]
    species_to_genus1 <- insert_spp[which(is.na(insert_spp[match(sub("_.*", 
                                                                     "", insert_spp), genus_in_tree)]) == FALSE)]
    if (length(species_to_genus1) >= 1) {
      if (progress.bar == TRUE) {
        pb_congeneric <- progress::progress_bar$new(format = "Adding congeneric species [:bar] :percent", 
                                                    total = length(species_to_genus1), clear = FALSE, 
                                                    width = 60, current = "<", incomplete = ">", 
                                                    complete = ">")
      }
      for (i in 1:length(species_to_genus1)) {
        tree_complete <- phytools::add.species.to.genus(tree = tree_complete, 
                                                      species = species_to_genus1[i])
        if (progress.bar == TRUE) {
          pb_congeneric$tick()
        }
      }
    }
    insert_spp2 <- treedata_modif(phy = tree_complete, 
                                  data = spp_data, warnings = F)$nc$data_not_tree
    if (length(insert_spp2) == 0) {
      data_final <- 1:length(as.character(data_samp$s))
      names(data_final) <- as.character(data_samp$s)
      tree_res <- tree_complete
      if (return.insertions == TRUE) {
        insertions <- rep("NA", nrow(data_samp))
        data_insertions <- cbind(data_samp, insertions)
        data_insertions[match(species_to_genus1, data_samp$s), 
                        "insertions"] <- rep("Congeneric_insertion", 
                                             length(species_to_genus1))
        spp_on_tree <- data_samp[-match(species_to_genus1, 
                                   data_samp$s), "s"]
        data_insertions[match(spp_on_tree, data_samp$s), 
                        "insertions"] <- rep("Present_in_Tree", length(spp_on_tree))
        list_res <- vector(mode = "list", length = 3)
        list_res[[1]] <- tree_res
        list_res[[2]] <- data_insertions
        list_res[[3]] <- spp_samp_return
        names(list_res) <- c("Phylogeny", "Insertions_data", "samp_spp")
        return(list_res)
      }
      else {
        return(tree_res)
      }
    }
    if (length(insert_spp2) >= 1) {
      data_exRound2 <- data_samp[match(insert_spp2, as.character(data_samp$s)), 
      ]
      rank_family2 <- unique(as.character(data_samp[match(insert_spp2, 
                                                     as.character(data_samp$s)), 2]))
      list_spp_step2 <- vector(mode = "list", length = length(rank_family2))
      for (i in 1:length(rank_family2)) {
        list_spp_step2[[i]] <- tryCatch(paste(ape::extract.clade(phy = tree_complete, 
                                                                 node = as.character(rank_family2[i]))$tip.label), 
                                        error = function(e) paste("noFamily", as.character(data_samp[which(rank_family2[i] == 
                                                                                                        data_samp$f), 1]), sep = "_"))
      }
      names(list_spp_step2) <- rank_family2
      data_exRound3 <- data_exRound2[!is.na(match(data_exRound2$f, 
                                                  names(which(unlist(lapply(lapply(list_spp_step2, 
                                                                                   function(x) which(sub("_.*", "", x) == "noFamily")), 
                                                                            function(y) length(y))) > 0)))), ]
      spp_family <- 1:nrow(data_exRound2)
      names(spp_family) <- data_exRound2$s
      spp_with_family <- names(which(unlist(lapply(lapply(list_spp_step2, 
                                                          function(x) which(sub("_.*", "", x) != "noFamily")), 
                                                   function(y) which(length(y) != 0))) > 0))
      spp_family_inTree <- list_spp_step2[match(spp_with_family, 
                                                names(list_spp_step2))]
      spp_to_add_round2 <- setdiff(data_exRound2$s, 
                                   data_exRound3$s)
      if (insert.base.node == TRUE) {
        pb_length <- unique(sub("_.*", "", as.character(spp_to_add_round2)))
        if (progress.bar == TRUE) {
          pb_insert_family_node <- progress::progress_bar$new(format = "Adding species to family node [:bar] :percent", 
                                                              total = length(pb_length), clear = FALSE, 
                                                              width = 60, current = "<", incomplete = ">", 
                                                              complete = ">")
        }
        count <- 0
        species_to_genus2 <- vector(mode = "list")
        while (length(spp_to_add_round2) >= 1) {
          count <- count + 1
          family_name <- data_samp[match(spp_to_add_round2[1], 
                                    data_samp$s), "f"]
          node_family <- which(c(tree_complete$tip.label, 
                                 tree_complete$node.label) == family_name)
          tree_complete <- phytools::bind.tip(tree = tree_complete, 
                                            tip.label = spp_to_add_round2[1], where = node_family, 
                                            position = 0)
          spp_data <- 1:length(spp_to_add_round2)
          names(spp_data) <- spp_to_add_round2
          insert_spp <- treedata_modif(phy = tree_complete, 
                                       data = spp_data, warnings = F)$nc$data_not_tree
          genus_in_tree <- sub("_.*", "", tree_complete$tip.label)[match(sub("_.*", 
                                                                           "", insert_spp), sub("_.*", "", tree_complete$tip.label))][!is.na(sub("_.*", 
                                                                                                                                               "", tree_complete$tip.label)[match(sub("_.*", 
                                                                                                                                                                                    "", insert_spp), sub("_.*", "", tree_complete$tip.label))])]
          if (length(genus_in_tree) >= 1) {
            species_to_genus <- insert_spp[which(is.na(insert_spp[match(sub("_.*", 
                                                                            "", insert_spp), genus_in_tree)]) == 
                                                   FALSE)]
            species_to_genus2[[count]] <- species_to_genus
            for (i in 1:length(species_to_genus)) {
              tree_complete <- phytools::add.species.to.genus(tree = tree_complete, 
                                                            species = species_to_genus[i])
            }
            insert_spp <- treedata_modif(phy = tree_complete, 
                                         data = spp_data, warnings = F)$nc$data_not_tree
          }
          spp_to_add_round2 <- setdiff(insert_spp, 
                                       data_exRound3$s)
          if (progress.bar == TRUE) {
            pb_insert_family_node$tick()
          }
        }
      }
      else {
        count <- 0
        species_to_genus2 <- vector(mode = "list")
        while (length(spp_to_add_round2) >= 1) {
          count <- count + 1
          family_name <- data_samp[match(spp_to_add_round2[1], 
                                    data_samp$s), "f"]
          spp_Byfamily_inTree <- as.character(unlist(spp_family_inTree[match(family_name, 
                                                                             names(spp_family_inTree))]))
          user_option_spp <- unique(sub("_.*", "", 
                                        as.character(unlist(spp_Byfamily_inTree))))
          local_to_add_spp <- readline(prompt = print_cat(print_cat = user_option_spp, 
                                                          spp = spp_to_add_round2[1], family_name))
          spp_user_opt <- unlist(strsplit(local_to_add_spp, 
                                          split = " "))
          if (length(spp_user_opt) > 2) {
            stop("/n Only two genus can be chosen to insert species")
          }
          if (length(spp_user_opt) < 1) {
            stop("/n At least one genus/family can be chosen to insert species")
          }
          if (any(is.na(match(spp_user_opt, user_option_spp) && 
                        is.na(match(spp_user_opt, family_name))))) {
            stop("/n Choose a validy genus/family to insert species")
          }
          if (length(spp_user_opt) == 2) {
            spp_to_add_tmp <- spp_Byfamily_inTree[match(spp_user_opt, 
                                                        sub("_.*", "", as.character(unlist(spp_Byfamily_inTree))))]
            node_btw_genus <- phytools::fastMRCA(tree = tree_complete, 
                                                 sp1 = spp_to_add_tmp[1], sp2 = spp_to_add_tmp[2])
            tree_complete <- phytools::bind.tip(tree = tree_complete, 
                                              tip.label = spp_to_add_round2[1], where = node_btw_genus, 
                                              position = 0)
          }
          if (length(spp_user_opt) == 1) {
            if (length(match(spp_user_opt, c(user_option_spp, 
                                             family_name))) != 1) {
              stop(paste("\n Check the spelling of Genus or Family in:", 
                         spp_user_opt, sep = " "))
            }
            if (any(spp_user_opt == family_name)) {
              node_family <- which(c(tree_complete$tip.label, 
                                     tree_complete$node.label) == family_name)
              tree_complete <- phytools::bind.tip(tree = tree_complete, 
                                                tip.label = spp_to_add_round2[1], where = node_family, 
                                                position = 0)
            }
            else {
              genus_nspp <- match(spp_user_opt, names(which(table(sub("_.*", 
                                                                      "", as.character(unlist(spp_Byfamily_inTree)))) > 
                                                              1)))
              if (!is.na(genus_nspp) == TRUE) {
                list_to_be_named <- list(spp_user_opt = spp_user_opt)
                name_list <- paste("MRCA", count, sep = "")
                names(list_to_be_named) <- name_list
                phylo_order2 <- ape::makeNodeLabel(phy = tree_complete, 
                                                   method = "u", nodeList = list_to_be_named)
                position_MRCA2 <- which(c(phylo_order2$tip.label, 
                                          phylo_order2$node.label) == name_list)
                check_name_node <- tree_complete$node.label[position_MRCA2 - 
                                                            length(tree_complete$tip.label)]
                if (check_name_node == family_name) {
                  list_to_be_named <- list(spp_user_opt = spp_user_opt)
                  name_list <- paste("MRCA", count, 
                                     sep = "")
                  names(list_to_be_named) <- name_list
                  tree_complete <- ape::makeNodeLabel(phy = tree_complete, 
                                                    method = "u", nodeList = list_to_be_named)
                  position_MRCA2 <- which(c(tree_complete$tip.label, 
                                            tree_complete$node.label) == name_list)
                  size_branch <- tree_complete$edge.length[sapply(position_MRCA2, 
                                                                function(x, y) which(y == x), y = tree_complete$edge[, 
                                                                                                                   2])]
                  tree_complete <- phytools::bind.tip(tree_complete, 
                                                    spp_to_add_round2[1], where = position_MRCA2, 
                                                    position = size_branch/2)
                  node_pos_new <- phytools::fastMRCA(tree = tree_complete, 
                                                     sp1 = spp_Byfamily_inTree[1], sp2 = spp_to_add_round2[1])
                  tree_complete$node.label[node_pos_new - 
                                           length(tree_complete$tip.label)] <- family_name
                }
                else {
                  tree_complete <- ape::makeNodeLabel(phy = tree_complete, 
                                                    method = "u", nodeList = list_to_be_named)
                  position_MRCA <- which(c(tree_complete$tip.label, 
                                           tree_complete$node.label) == name_list)
                  size_branch <- tree_complete$edge.length[sapply(position_MRCA, 
                                                                function(x, y) which(y == x), y = tree_complete$edge[, 
                                                                                                                   2])]
                  tree_complete <- phytools::bind.tip(tree_complete, 
                                                    spp_to_add_round2[1], where = position_MRCA, 
                                                    position = size_branch/2)
                }
              }
              if (is.na(genus_nspp) == TRUE) {
                tree_complete <- phytools::add.species.to.genus(tree = tree_complete, 
                                                              species = paste(sub("_.*", "", as.character(spp_user_opt)), 
                                                                              "toadd", sep = "_"))
                position_problem2 <- which(tree_complete$tip.label == 
                                             paste(sub("_.*", "", as.character(spp_user_opt)), 
                                                   "toadd", sep = "_"))
                tree_complete$tip.label[position_problem2] <- spp_to_add_round2[1]
              }
            }
          }
          spp_data <- 1:nrow(data_exRound2)
          names(spp_data) <- data_exRound2$s
          insert_spp <- treedata_modif(phy = tree_complete, 
                                       data = spp_data, warnings = F)$nc$data_not_tree
          genus_in_tree <- sub("_.*", "", tree_complete$tip.label)[match(sub("_.*", 
                                                                           "", insert_spp), sub("_.*", "", tree_complete$tip.label))][!is.na(sub("_.*", 
                                                                                                                                               "", tree_complete$tip.label)[match(sub("_.*", 
                                                                                                                                                                                    "", insert_spp), sub("_.*", "", tree_complete$tip.label))])]
          if (length(genus_in_tree) >= 1) {
            species_to_genus <- insert_spp[which(is.na(insert_spp[match(sub("_.*", 
                                                                            "", insert_spp), genus_in_tree)]) == 
                                                   FALSE)]
            species_to_genus2[[count]] <- species_to_genus
            for (i in 1:length(species_to_genus)) {
              tree_complete <- phytools::add.species.to.genus(tree = tree_complete, 
                                                            species = species_to_genus[i])
            }
            insert_spp <- treedata_modif(phy = tree_complete, 
                                         data = spp_data, warnings = F)$nc$data_not_tree
          }
          rank_family2 <- unique(as.character(data_samp[match(insert_spp, 
                                                         as.character(data_samp$s)), 2]))
          if (length(rank_family2) == 0) {
            data_final <- 1:length(as.character(data_samp$s))
            names(data_final) <- as.character(data_samp$s)
            tree_res <- tree_complete
            data_exRound3 <- NULL
            break
          }
          else {
            list_spp_step2 <- vector(mode = "list", 
                                     length = length(rank_family2))
            for (i in 1:length(rank_family2)) {
              list_spp_step2[[i]] <- tryCatch(paste(ape::extract.clade(phy = tree_complete, 
                                                                       node = as.character(rank_family2[i]))$tip.label), 
                                              error = function(e) paste("noFamily", 
                                                                        as.character(data_samp[which(rank_family2[i] == 
                                                                                                  data_samp$f), 1]), sep = "_"))
            }
            names(list_spp_step2) <- rank_family2
            spp_family <- 1:nrow(data_exRound2)
            names(spp_family) <- data_exRound2$s
            spp_with_family <- names(which(unlist(lapply(lapply(list_spp_step2, 
                                                                function(x) which(sub("_.*", "", x) != 
                                                                                    "noFamily")), function(y) which(length(y) != 
                                                                                                                      0))) > 0))
            spp_family_inTree <- list_spp_step2[match(spp_with_family, 
                                                      names(list_spp_step2))]
            spp_to_add_round2 <- setdiff(insert_spp, 
                                         data_exRound3$s)
          }
        }
      }
      if (is.null(data_exRound3) | dim(data_exRound3)[1] == 0) {
        if (return.insertions == TRUE) {
          data_final <- 1:length(as.character(data_samp$s))
          names(data_final) <- as.character(data_samp$s)
          tree_res <- tree_complete
          insertions <- rep("NA", nrow(data_samp))
          data_insertions <- cbind(data_samp, insertions)
          Present_in_tree <- data_samp$s[!is.na(match(data_samp$s, 
                                                 tree_complete$tip.label))]
          Congeneric_insertion <- species_to_genus1
          not_inserted <- data_samp$s[is.na(match(data_samp$s, 
                                             tree_res$tip.label))]
          Congeneric_round_family <- unlist(species_to_genus2)
          if(is.null(Congeneric_round_family) & length(not_inserted) == 0){
            data_insertions[match(insert_spp2, data_insertions$s), "insertions"] <- "Family_insertion"
          } else {
            Family_insertion <- insert_spp2[-match(c(Congeneric_round_family, 
                                                     not_inserted), insert_spp2)]
            data_insertions[match(Family_insertion, data_insertions$s), 
                            "insertions"] <- "Family_insertion"
          }
          
          data_insertions[match(Present_in_tree, data_insertions$s), 
                          "insertions"] <- "Present_in_Tree"
          data_insertions[match(Congeneric_insertion, 
                                data_insertions$s), "insertions"] <- "Congeneric_insertion"
          
          if (length(not_inserted) >= 1) {
            data_insertions[match(not_inserted, data_insertions$s), 
                            "insertions"] <- "Not_inserted"
          }
          if (length(Congeneric_round_family) >= 1) {
            data_insertions[match(Congeneric_round_family, 
                                  data_insertions$s), "insertions"] <- "Congeneric_Family_level"
          }
          list_res <- vector(mode = "list", length = 3)
          list_res[[1]] <- tree_res
          list_res[[2]] <- data_insertions
          list_res[[3]] <- spp_samp_return
          names(list_res) <- c("Phylogeny", "Insertions_data", "samp_spp")
          return(list_res)
        }
        else {
          return(tree_res)
        }
      }
      else {
        data_exRound3 <- data_exRound3[is.na(match(data_exRound3$o, 
                                                   order_rm_list)), ]
        rank_order_Round3 <- rank_order[match(data_exRound3$o, 
                                              rank_order)]
        families_round3 <- lapply(lapply(rank_order_Round3, 
                                         function(x) {
                                           fishbasedata[which(x == fishbasedata$Order), 
                                                        5]
                                         }), function(y) unique(y))
        names(families_round3) <- rank_order_Round3
        orders_round_3 <- unique(rank_order_Round3)
        families_round3_check <- unique(as.character(unlist(families_round3)))
        check_round3_insertion <- match(families_round3_check, 
                                        families_in_tree)
        check_allNA_round3 <- unique(check_round3_insertion)
        if (all(is.na(check_allNA_round3)) == TRUE) {
          no_represent_tree <- data_exRound3$s
          data_final <- 1:length(as.character(data_samp$s))
          names(data_final) <- as.character(data_samp$s)
          tree_res <- tree_complete
          if (return.insertions == TRUE) {
            insertions <- rep("NA", nrow(data_samp))
            species_to_genus2 <- unlist(species_to_genus2)
            data_insertions <- cbind(data_samp, insertions)
            Present_in_tree <- data_samp$s[!is.na(match(data_samp$s, 
                                                   tree_complete$tip.label))]
            Congeneric_insertion <- species_to_genus1
            not_inserted <- data_samp$s[is.na(match(data_samp$s, 
                                               tree_res$tip.label))]
            Congeneric_round_family <- species_to_genus2
            if(is.null(Congeneric_round_family) & length(not_inserted) == 0){
              data_insertions[match(insert_spp2, data_insertions$s), "insertions"] <- "Family_insertion"
            } else {
              Family_insertion <- insert_spp2[-match(c(Congeneric_round_family, 
                                                       not_inserted), insert_spp2)]
              data_insertions[match(Family_insertion, data_insertions$s), 
                              "insertions"] <- "Family_insertion"
            }
            data_insertions[match(Present_in_tree, 
                                  data_insertions$s), "insertions"] <- "Present_in_Tree"
            data_insertions[match(Congeneric_insertion, 
                                  data_insertions$s), "insertions"] <- "Congeneric_insertion"
            if (length(not_inserted) >= 1) {
              data_insertions[match(not_inserted, data_insertions$s), 
                              "insertions"] <- "Not_inserted"
            }
            if (length(Congeneric_round_family) >= 
                1) {
              data_insertions[match(Congeneric_round_family, 
                                    data_insertions$s), "insertions"] <- "Congeneric_Family_level"
            }
            list_res <- vector(mode = "list", length = 3)
            list_res[[1]] <- tree_res
            list_res[[2]] <- data_insertions
            list_res[[3]] <- spp_samp_return
            names(list_res) <- c("Phylogeny", "Insertions_data", "samp_spp")
            return(list_res)
          }
          return(tree_res)
        }
        order_add_round3 <- unique(data_exRound3$o)
        if (progress.bar == TRUE) {
          pb_insert_order_node_remaining <- progress::progress_bar$new(format = "Adding species to order node [:bar] :percent", 
                                                                       total = nrow(data_exRound3), clear = FALSE, 
                                                                       width = 60, current = "<", incomplete = ">", 
                                                                       complete = ">")
        }
        for (i in 1:length(order_add_round3)) {
          list_names_round3 <- list_order[[which(names(list_order) == 
                                                   order_add_round3[i])]]
          tree_complete <- ape::makeNodeLabel(tree_complete, 
                                            "u", nodeList = list(Ord_name = list_names_round3))
          tree_complete$node.label[which(tree_complete$node.label == 
                                         "Ord_name")] <- order_add_round3[i]
        }
        if (insert.base.node == TRUE) {
          for (i in 1:nrow(data_exRound3)) {
            node_order_pos <- which(c(tree_complete$tip.label, 
                                      tree_complete$node.label) == data_exRound3$o[i])
            tree_complete <- phytools::bind.tip(tree = tree_complete, 
                                              tip.label = data_exRound3$s[i], where = node_order_pos, 
                                              position = 0)
            pb_insert_order_node_remaining$tick()
          }
        }
        else {
          for (i in 1:length(families_round3)) {
            user_option_family <- tree_complete$node.label[match(unlist(families_round3[[i]]), 
                                                               tree_complete$node.label)[-which(match(match(unlist(families_round3[[i]]), 
                                                                                                          tree_complete$node.label), NA) == 1)]]
            local_to_add_spp_family <- readline(prompt = print_cat_family(print_cat = unlist(user_option_family), 
                                                                          spp = data_exRound3$s[i], data_exRound3$o[i]))
            family_user_opt <- unlist(strsplit(local_to_add_spp_family, 
                                               split = " "))
            if (length(family_user_opt) == 1) {
              if (family_user_opt == data_exRound3$o[i]) {
                node_order_pos <- which(c(tree_complete$tip.label, 
                                          tree_complete$node.label) == data_exRound3$o[i])
                tree_complete <- phytools::bind.tip(tree = tree_complete, 
                                                  tip.label = data_exRound3$s[i], where = node_order_pos, 
                                                  position = 0)
              }
              family_nspp <- length(list_family[[match(family_user_opt, 
                                                       names(list_family))]])
              if (family_nspp > 1) {
                position_family <- which(family_user_opt == 
                                           c(tree_complete$tip.label, tree_complete$node.label))
                size_branch_family <- tree_complete$edge.length[sapply(position_family, 
                                                                     function(x, y) which(y == x), y = tree_complete$edge[, 
                                                                                                                        2])]
                tree_complete <- phytools::bind.tip(tree_complete, 
                                                  data_exRound3$s[i], where = position_family, 
                                                  position = size_branch_family/2)
              }
              if (family_nspp == 1) {
                tree_complete <- phytools::add.species.to.genus(tree = tree_complete, 
                                                              species = paste(sub("_.*", "", as.character(list_family[[which(names(list_family) == 
                                                                                                                               family_user_opt)]][1])), "toadd", 
                                                                              sep = "_"))
                position_problem_family <- which(tree_complete$tip.label == 
                                                   paste(sub("_.*", "", as.character(list_family[[which(names(list_family) == 
                                                                                                          family_user_opt)]][1])), "toadd", 
                                                         sep = "_"))
                tree_complete$tip.label[position_problem_family] <- data_exRound3$s[i]
              }
            }
            if (length(family_user_opt) == 2) {
              spp_to_add_tmp_family1 <- ape::extract.clade(tree_complete, 
                                                           node = family_user_opt[1])$tip.label[1]
              spp_to_add_tmp_family2 <- ape::extract.clade(tree_complete, 
                                                           node = family_user_opt[2])$tip.label[1]
              node_btw_genus_family <- phytools::fastMRCA(tree = tree_complete, 
                                                          sp1 = spp_to_add_tmp_family1, sp2 = spp_to_add_tmp_family2)
              tree_complete <- phytools::bind.tip(tree = tree_complete, 
                                                tip.label = data_exRound3$s[i], where = node_btw_genus_family, 
                                                position = 0)
            }
          }
        }
        data_final <- 1:length(as.character(data_samp$s))
        names(data_final) <- as.character(data_samp$s)
        tree_res <- tree_complete
        if (return.insertions == TRUE) {
          insertions <- rep("NA", nrow(data_samp))
          data_insertions <- cbind(data_samp, insertions)
          species_to_genus2 <- unlist(species_to_genus2)
          Present_in_tree <- data_samp$s[!is.na(match(data_samp$s, 
                                                 tree_complete$tip.label))]
          Congeneric_insertion <- species_to_genus1
          not_inserted <- data_samp$s[is.na(match(data_samp$s, 
                                             tree_res$tip.label))]
          Congeneric_round_family <- species_to_genus2
          if(is.null(Congeneric_round_family) & length(not_inserted) == 0){
            data_insertions[match(insert_spp2, data_insertions$s), "insertions"] <- "Family_insertion"
          } else {
            Family_insertion <- insert_spp2[-match(c(Congeneric_round_family, 
                                                     not_inserted), insert_spp2, nomatch = 0)]
            data_insertions[match(Family_insertion, data_insertions$s), 
                            "insertions"] <- "Family_insertion"
          }
          data_insertions[match(Present_in_tree, data_insertions$s), 
                          "insertions"] <- "Present_in_Tree"
          data_insertions[match(Congeneric_insertion, 
                                data_insertions$s), "insertions"] <- "Congeneric_insertion"
          if (length(not_inserted) >= 1) {
            data_insertions[match(not_inserted, data_insertions$s), 
                            "insertions"] <- "Not_inserted"
          }
          if (length(Congeneric_round_family) >= 1) {
            data_insertions[match(Congeneric_round_family, 
                                  data_insertions$s), "insertions"] <- "Congeneric_Family_level"
          }
          if (dim(data_exRound3)[1] >= 1) {
            data_insertions[match(data_exRound3$s, 
                                  data_insertions$s), "insertions"] <- "Order_insertion"
          }
          list_res <- vector(mode = "list", length = 3)
          list_res[[1]] <- tree_res
          list_res[[2]] <- data_insertions
          list_res[[3]] <- spp_samp_return
          names(list_res) <- c("Phylogeny", "Insertions_data", "samp_spp")
          return(list_res)
        }
        else {
          return(tree_res)
        }
      }
    }
  }
  
}
