#####
# find_all_parents
#   find the parent links for all active nodes
#
#    Inputs:
#      mat      = subset of the link matrix indicating links for main effects
#      row_inds = current active links
#
#    Outputs:
#      parents  = list of indices of the parent links
#
#    Function Calls:
#      None
#
#    Called From:
#      add_links
#
#    Most Recent Author:
#      Mike Grosskopf
#
#    Most Recent Edit:
#      2015-12-15
#
#####
find_all_parents <- function(mat,row_inds){
  parents <- unique(unlist(apply(mat[row_inds,,drop=F],1,function(x){which(x == 1)})))
  parents
}

#####
# find_all_children
#   find the child links for all active nodes
#
#    Inputs:
#      mat      = subset of the link matrix indicating links for main effects
#      row_inds = current active links
#
#    Outputs:
#      children  = list of indices of the parent links
#
#    Function Calls:
#      None
#
#    Called From:
#      add_links
#
#    Most Recent Author:
#      Mike Grosskopf
#
#    Most Recent Edit:
#      2015-12-15
#
#####
find_all_children <- function(mat,row_inds){
  has_children <- unique(unlist(apply(mat[row_inds,,drop=F],1,function(x){which(x == -1)})))
  has_children
}

#####
# check_for_parent
#   ensure each potential parent link if valid by checking to see 
#   if the parent node has a link going into it
#
#    Inputs:
#      mat      = subset of the link matrix indicating links for main effects
#      row_inds = current active links
#
#    Outputs:
#      boolean indicated whether the link is valid (in the sense of above)
#
#    Function Calls:
#      None
#
#    Called From:
#      add_links
#
#    Most Recent Author:
#      Mike Grosskopf
#
#    Most Recent Edit:
#      2015-12-15
#
#####
check_for_parent <- function(mat,col_ind,active_parents){
  parent_ind       <- which( mat[,col_ind] == -1 )
  its_parents_inds <- which( mat[parent_ind,] == 1 )
  # Consider its parents active if it is a top node
  if (length(its_parents_inds) == 0) return(TRUE)
  
  any(its_parents_inds %in% active_parents)
}

#####
# add_interactions
#   Adds 1's to bit string if the interaction should be activated. This
#   is done when both main effects for the interaction are active
#
#    Inputs:
#      mat              = subset of the link matrix indicating links for main effects
#      active_links     = all active links in the model
#      main_effect_inds = list of indices of the main effects in bit string
#
#    Outputs:
#      indices of any interactions to be added to the bit string
#
#    Function Calls:
#      None
#
#    Called From:
#      add_link
#
#    Most Recent Author:
#      Mike Grosskopf
#
#    Most Recent Edit:
#      2015-12-15
#
#####
add_interactions <- function(mat,active_links,main_effect_inds){
  interaction_columns <- (1:ncol(mat))[-main_effect_inds]
  
  add_indicator  <- rep(FALSE, ncol(mat))
  
  for (i in interaction_columns) {
    parents <- which( mat[,i] == -1 )
    child   <- which( mat[,i] ==  1 )
    
    check_one <- any(apply(mat[-c(parents[1],child),active_links],2,function(x){all(x==0)}))
    check_two <- any(apply(mat[-c(parents[2],child),active_links],2,function(x){all(x==0)}))
    
    add_indicator[i] <- check_one && check_two
  }
  
  which(add_indicator)
}

#####
#  add_link:
#     Sequential link addition from initial small subgraph to full network
#
#    Inputs:
#      init_bit_string = bit string representing the currently active/inactive links in the subgraph
#      links_mat       = model$interMat from relevant fuzzy Logic model. Indicates all the links and which
#                            nodes they point between
#
#    Outputs:
#      new_b_string    = new bit string with new link activated. If the graph is complete
#                          the function prints a message and returns the previous bit string
#
#    Function Calls:
#      find_all_parents
#      find_all_children
#
#    Called From:
#      smc_test.R
#
#    Most Recent Author:
#      Mike Grosskopf
#
#    Most Recent Edit:
#      2015-10-30
#
#####
add_link <- function(init_bit_string,links_mat){
  main_effect_inds    <- 1:length(init_bit_string)
  new_bit_string      <- init_bit_string
  subgraph_bit_string <- as.logical(init_bit_string[main_effect_inds])
  mainEffects_mat     <- links_mat[,main_effect_inds]
  
  # Identify top nodes
  top_nodes      <- which(apply(mainEffects_mat,1,function(x){all(x != 1)}))
  
  # Find all nodes and links that are active in the main effects
  active_nodes   <- unique(c(apply(mainEffects_mat[,subgraph_bit_string,drop=F],2,function(x){which(x != 0)})))
  active_links   <- which(subgraph_bit_string)

  # Find the parents of all nodes in an active link
  active_parents <- find_all_parents(mainEffects_mat,active_nodes)
  # Find inactive links from a parent of an active node
  #     then check to ensure those potential parents have
  #     at least one active parent or are a top node
  potential_parents <- setdiff(active_parents, active_links)
  if (length(potential_parents) > 0){
    check_parents     <- sapply(potential_parents, check_for_parent, 
                                mat = mainEffects_mat,
                                active_parents = active_parents
    )
    potential_parents <- potential_parents[check_parents]
  }
  
  if( length(potential_parents) == 0 ) {
    
    active_children <- find_all_children(mainEffects_mat,unique(c(active_nodes,top_nodes)))

    if (all(active_children %in% active_links)){
      
      print( "No new nodes to add. Graph complete.")
      return(init_bit_string)
      
    } else {
      
      potential_children <- setdiff(union(active_children,active_parents), active_links)

      if (length(potential_children) == 1 ){
        new_link <- potential_children
      } else{
        new_link <- sample(potential_children, 1)
      }
      
    }
    
  } else {
    
    if (length(potential_parents) == 1 ){
      new_link <- potential_parents
    } else{
      new_link <- sample(potential_parents, 1)
    }
    
  }
  
  new_bit_string[new_link] <- 1
  
  active_links     <- c(active_links, new_link)
  # new_interactions <- add_interactions(links_mat,active_links,main_effect_inds)
  
  # new_bit_string[new_interactions] <- 1
  
  new_bit_string
}
