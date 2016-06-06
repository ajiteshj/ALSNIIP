iacal <- function(ip_file = " ", file_sep = " ", op_path = " "){
  
#Load packages
require(ProNet)
require(igraph)
require(abind)
require(sna)
require(ergm.count)
require(network)
require(tcltk)
require(tcltk2)

#Take user input and convert to graph object
input <- read.csv(ip_file , sep = file_sep, header = FALSE)
graph <- simplify(graph.data.frame(input, directed = FALSE))

#Initialise variables
percentage <- 0.1
major_cluster_count <- vector("numeric", 9)
temporary <- vector("numeric", 9)
difference <- vector("numeric", 9)
result <- 0

#sample until stopping criterion is satisfied
stop_index <- 1
while(percentage < 1){
  #Sample percentage amount of vertices and cluster 
  length <- vcount(graph) * percentage
  set.seed(123)
  subgraph <- extraction(graph, mode = "sample", sample.number = length)
  fastgreedy_results <- cluster_fast_greedy(subgraph)
  V(subgraph)$color <- fastgreedy_results$membership
  plot(subgraph, main = paste(percentage * 100, "% sample"), 
       vertex.color = V(subgraph)$color, layout = layout.fruchterman.reingold)
  visualise(object = subgraph, vcolor = V(subgraph)$color, 
            save_file = paste(op_path, "/sample/sample_", percentage * 100, ".eps", sep = ""))
  
  #Store only those clusters that have greater than 1% of total number of vertices
  #that were sampled in major_size_fastgreedy_results
  size_fastgreedy_results <- sizes(fastgreedy_results)
  cut_off <- 0.01*length
  major_size_fastgreedy_results <- size_fastgreedy_results[size_fastgreedy_results>cut_off]
  major_cluster_count[stop_index] <- length(major_size_fastgreedy_results)
  
  #Stopping criterion
  if(stop_index == 1){
    temporary[1] = major_cluster_count[stop_index]
  }else if(stop_index == 2){
    temporary[2] = major_cluster_count[stop_index]
    difference[1] = abs(temporary[2] - temporary[1])
  }else if(stop_index == 3){
    temporary[3] = major_cluster_count[stop_index]
    difference[2] = abs(temporary[3] - temporary[2])
    result = abs(difference[2] - difference[1])
  }else{
    temporary[stop_index] = major_cluster_count[stop_index]
    difference[stop_index-1] = abs(temporary[stop_index] - temporary[stop_index-1])
    if(abs(difference[stop_index-1] - difference[stop_index-2]) <= result){
      result = abs(difference[stop_index-1] - difference[stop_index-2])
      break
    }else{
      result = abs(difference[stop_index-1] - difference[stop_index-2])
    }
  }
  stop_index <- stop_index + 1
  percentage = percentage + 0.1
}

#Stabilised graph and its clusters
print("Stabilised sampled graph:")
print(subgraph)
print("Comunities in the graph:")
print(fastgreedy_results)

#Store those clusters in max_fastgreedy_results that have number of vertices greater
#than 1% of max_subgraph in final 
final_cut_off <- 0.01 * vcount(subgraph)
final <- fastgreedy_results[sizes(fastgreedy_results) > final_cut_off]
print("Major clusters:")
print(final)

#final_matrix and final_edge has all the clusters as matrices
final_matrix <- matrix(nrow = nrow(get.edgelist(graph)), ncol = 2)
for(final_matrix_cluster in 1:length(names(final))){
  if(final_matrix_cluster != 1){
    final_matrix <- abind(final_matrix, final_matrix, along = 3)
  }
}

final_edge <- matrix(NA, ncol=2)
final_index <- 1
for(final_names in names(final)){
  noneed <- V(subgraph)[membership(fastgreedy_results) != as.integer(final_names)]
  initial_cluster <- delete_vertices(subgraph, noneed)
  initial_cluster_edge <- get.edgelist(initial_cluster)
  final_matrix[1:nrow(initial_cluster_edge),,final_index] <- initial_cluster_edge
  final_index <- final_index + 1
  if(final_names == names(final)[1]){
    final_edge <- initial_cluster_edge
  }else{
    final_edge <- rbind(final_edge, initial_cluster_edge)
  }
}

#not_yet_edge has all the edges that have not yet been clustered
not_yet_edge <- matrix(NA, ncol = 2)
not_yet_edge <- rbind(get.edgelist(graph))
for(final_edge_row in 1:nrow(final_edge)){
  not_yet_edge_row <- 1
  while(not_yet_edge_row <= nrow(not_yet_edge)){
    if((final_edge[final_edge_row,1] == not_yet_edge[not_yet_edge_row,1]) && (final_edge[final_edge_row,2] == not_yet_edge[not_yet_edge_row,2])){
      not_yet_edge <- not_yet_edge[-not_yet_edge_row,]
      break
    }else{
      not_yet_edge_row <- not_yet_edge_row + 1
    }
  }
}

#for every edge in not_yet_edge
while(nrow(not_yet_edge) > 0){

  #not_yet_vertices has vertices of not_yet_graph
  not_yet_vertices <- V(graph.data.frame(not_yet_edge))

  #all_pageranks is a vector that stores the pagerank of not_yet in each cluster 
  all_pageranks <- vector("numeric", length(names(final)))

  #not_yet_random contains randomly selected node from not_yet and removed from not_yet
  not_yet_random <- sample(not_yet_vertices, 1)
  not_yet_random <- names(not_yet_random)
  print("Randomly selected vertex to be clustered:")
  print(not_yet_random)

  #edge_index contains the indices of all the edges corresponding to not_yet_random
  edge_index <- which(not_yet_edge == not_yet_random)
  correct_index <- which(edge_index > nrow(not_yet_edge))
  edge_index[correct_index] <- edge_index[correct_index] - nrow(not_yet_edge)
  edge_index <- unique(edge_index)

  #not_yet_random_edge contains edges corresponding to not_yet_random
  not_yet_random_edge <- matrix(NA, ncol = 2)
  if(NA %in% not_yet_random_edge){
    not_yet_random_edge <- rbind(not_yet_edge[edge_index,])
  }else{
    not_yet_random_edge <- rbind(not_yet_random_edge, not_yet_edge[edge_index,])
  }
  print("Edges corresponding to the randomly selected vertex:")
  print(not_yet_random_edge)

  #remove not_yet_random_edges from not_yet_edge
  not_yet_edge <- not_yet_edge[-edge_index,]
  
  if(class(not_yet_edge) == "character"){
    not_yet_edge <- t(as.matrix(not_yet_edge))  
  }
  
  #for all clusters in final_matrix
  for(final_matrix_index in 1:length(names(final))){
    for(final_matrix_row in 1:nrow(final_matrix[,,final_matrix_index])){
      if((NA %in% final_matrix[final_matrix_row,1,final_matrix_index]) == TRUE){
        stop_value <- final_matrix_row
        break
      }
    }
    
    #find the indexed cluster
    nonNA <- stop_value - 1
    cluster_index_edge <- matrix(NA, nrow = nonNA, ncol = 2)
    cluster_index_edge[1:nonNA, 1] <- final_matrix[1:nonNA, 1, final_matrix_index]
    cluster_index_edge[1:nonNA, 2] <- final_matrix[1:nonNA, 2, final_matrix_index]
    
    #add those edges in not_yet_edge to the indexed cluster if it has a link to a vertex in that cluster
    for(not_yet_index in 1:nrow(not_yet_random_edge)){
      if(((not_yet_random == not_yet_random_edge[not_yet_index,1]) && (not_yet_random_edge[not_yet_index,2] %in% cluster_index_edge) == TRUE) || ((not_yet_random == not_yet_random_edge[not_yet_index,2]) && (not_yet_random_edge[not_yet_index,1] %in% cluster_index_edge) == TRUE)){
        cluster_index_edge <- rbind(cluster_index_edge, not_yet_random_edge[not_yet_index,])
      }
    }

    #find pagerank of not_yet_random in the indexed cluster iff it is present 
    #else assign it as zero 
    if((not_yet_random %in% cluster_index_edge) == TRUE){
      cluster_index_graph <- graph.data.frame(cluster_index_edge, directed = TRUE)
      cluster_index_pagerank <- page_rank(cluster_index_graph)$vector
      all_pageranks[final_matrix_index] <- cluster_index_pagerank[[not_yet_random]] 
    }else{
      all_pageranks[final_matrix_index] <- 0
    }
  }
  print("PageRank of the randomly selected vertex in all the clusters:")
  print(all_pageranks)

  #if not_yet_random cannot be placed in any of the existing clusters then reject all the corresponding edges 
  #else add it to the identified cluster
  if(sum(all_pageranks) > 0){
    #find the cluster in which the not_yet element has highest pagerank
    index <- which(max(all_pageranks) == all_pageranks)
    print("Clustered to which the randomly selected vertex has been determined to belong:")
    print(index)
    
    #find the index of the last edge
    for(final_matrix_row in 1:nrow(final_matrix[,,index])){
      if((NA %in% final_matrix[final_matrix_row,1,index]) == TRUE){
        stop_value <- final_matrix_row
        break
      }
    }
    
    #find the indexed cluster
    nonNA <- stop_value - 1
    
    #add_ith_cluster_edge has edges in final_matrix[,,index]
    add_ith_cluster_edge <- matrix(NA, nrow = nonNA, ncol = 2)
    add_ith_cluster_edge[1:nonNA, 1] <- final_matrix[1:nonNA, 1, index]
    add_ith_cluster_edge[1:nonNA, 2] <- final_matrix[1:nonNA, 2, index]
    
    #add_final_edge contains those edges that are to be added
    add_final_edge <- matrix(NA, ncol = 2)
    for(add_not_yet_index in 1:nrow(not_yet_random_edge)){
      if((not_yet_random_edge[add_not_yet_index,1] %in% add_ith_cluster_edge) == TRUE || (not_yet_random_edge[not_yet_index,2] %in% add_ith_cluster_edge) == TRUE){
        if(NA %in% add_final_edge){
          add_final_edge <- rbind(not_yet_random_edge[add_not_yet_index,])
        }else{
          add_final_edge <- rbind(add_final_edge, not_yet_random_edge[add_not_yet_index,])
        }
      }
    }

    #find start and end index to add edges to final_matrix[,,index] and put add_final_edge there
    start_index <- stop_value
    end_index <- stop_value + nrow(add_final_edge) - 1
    final_matrix[start_index:end_index, 1, index] <- add_final_edge[1:nrow(add_final_edge), 1]
    final_matrix[start_index:end_index, 2, index] <- add_final_edge[1:nrow(add_final_edge), 2]
  }
}

#Convert all clusters to graph objects and find the influentials in each of them
for(each_cluster in 1:length(names(final))){
  for(final_matrix_row in 1:nrow(final_matrix[,,each_cluster])){
    if((NA %in% final_matrix[final_matrix_row,1,each_cluster]) == TRUE){
      stop_value <- final_matrix_row
      break
    }
  }
  
  #find the indexed cluster
  nonNA <- stop_value - 1
  real_cluster <- matrix(NA, ncol = 2)
  real_cluster <- rbind(final_matrix[1:nonNA, , each_cluster])
  cluster_graph <- graph.data.frame(real_cluster, directed = TRUE)
  plot(main = paste("Cluster", each_cluster, sep = " "), cluster_graph, 
       vertex.color = "gold", layout = layout.fruchterman.reingold)
  visualise(object = cluster_graph, 
            save_file = paste(op_path, "/cluster/cluster_", each_cluster, ".eps", sep = ""))
  
  #get adjacency matrix of the cluster
  cluster_graph_adj <- get.adjacency(cluster_graph, sparse = FALSE)
  
  #generate screeplot of the dataframe and plot
  in_ties <- sna::degree(cluster_graph_adj, gmode = "digraph", cmode = "indegree")
  names(in_ties) <- names(V(cluster_graph))
  desc_in_ties <- sort(in_ties, decreasing = TRUE)
  plot(desc_in_ties, main = paste("Scree plot of cluster", each_cluster, 
                                  sep = " "), xlab = "Actor", ylab = "In-Ties")
  legend('topright', c("Absolute Cut", "Fixed %", "Standard Deviation") , 
         lty=1, col=c('red', 'blue', 'green'), bty='n', cex=1)
  
  #absolute cut score
  cutoff_abs <- 0.8 * max(desc_in_ties)
  abline(h = cutoff_abs, col = "red") 
  abscut <- names(V(cluster_graph)[in_ties > cutoff_abs])
  
  #Fixed percentage of population
  cutoff_fix <- 0.2 * nrow(cluster_graph_adj)
  abline(v = cutoff_fix, col = "blue")
  fixcut <- names(head(desc_in_ties, cutoff_fix))
  
  #Standard deviation
  cutoff_sd <- mean(in_ties) + 2 * sd(in_ties)
  abline(h = cutoff_sd, col = "green")
  sdcut <- names(V(cluster_graph)[in_ties > cutoff_sd])
  
  #save plot as Screeplot<no>.png
  png(file = paste(op_path, "/screeplot/screeplot_", each_cluster, ".png", sep =""),
      height = 1600, width = 1000)
  plot(desc_in_ties, main = paste("Scree plot of cluster", each_cluster, 
                                  sep = " "), xlab = "Actor", ylab = "In-Ties")
  legend('topright', c("Absolute Cut", "Fixed %", "Standard Deviation") , 
         lty=1, col=c('red', 'blue', 'green'), bty='n', cex=1)
  abline(h = cutoff_abs, col = "red") 
  abline(v = cutoff_fix, col = "blue")
  abline(h = cutoff_sd, col = "green")
  dev.off()
  
  #color nodes of abscut and plot these iff there are influentials
  if(length(abscut) > 0){
    color <- vector("numeric", length = length(V(cluster_graph)))
    index <- 1
    while(index <= length(abscut)){
      c <- which(names(V(cluster_graph)) == abscut[index])
      color[c] <- 1
      index <- index + 1
    }
    plot(cluster_graph, main = paste("Influentials under Absolute Cut Score Method in cluster", each_cluster, sep = " "),
         layout = layout.fruchterman.reingold, vertex.color = ifelse(color==1, "red", "gold")) 
    legend('topright', c("Influentials", "Non-influentials"),
           col=c('red', 'gold'), bty='n', cex=1, pch=19)
    visualise(object = cluster_graph, vcolor = ifelse(color==1, "red", "gold"), 
              save_file = paste(op_path, "/abscut/abscut_", each_cluster, ".eps", sep = ""))
  }else{
    print(paste("No influentials have been generated by absolute cut method for cluster", each_cluster, sep = " "))
  }
  
  #color nodes of fixed% and plot these iff there are influentials
  if(length(fixcut) > 0){
    color <- vector("numeric", length = length(V(cluster_graph)))
    index <- 1
    while(index <= length(fixcut)){
      c <- which(names(V(cluster_graph)) == fixcut[index])
      color[c] <- 1
      index <- index + 1
    }
    plot(cluster_graph, main = paste("Influentials under Fixed Percentage of Population Method in cluster", each_cluster, sep = " "),
         layout = layout.fruchterman.reingold, vertex.color = ifelse(color==1, "dodgerblue", "gold")) 
    legend('topright', c("Influentials", "Non-influentials"),
           col=c('dodgerblue', 'gold'), bty='n', cex=1, pch=19)
    visualise(object = cluster_graph, vcolor = ifelse(color==1, "dodgerblue", "gold"), 
              save_file = paste(op_path, "/fixcut/fixcut_", each_cluster, ".eps", sep = ""))
  }else{
    print(paste("No influentials have been generated by fixed percentage of population method for cluster", each_cluster, sep = " "))
  }
  
  #color nodes of sd and plot these iff there are influentials  
  if(length(sdcut) > 0){
    color <- vector("numeric", length = length(V(cluster_graph)))
    index <- 1
    while(index <= length(fixcut)){
      c <- which(names(V(cluster_graph)) == fixcut[index])
      color[c] <- 1
      index <- index + 1
    }
    plot(cluster_graph, main = paste("Influentials under Standard Deviation Method in cluster", each_cluster, sep = " "),
         layout = layout.fruchterman.reingold, vertex.color = ifelse(color==1, "green", "gold")) 
    legend('topright', c("Influentials", "Non-influentials"),
           col=c('green', 'gold'), bty='n', cex=1, pch=19)
    visualise(object = cluster_graph, vcolor = ifelse(color==1, "green", "gold"), 
              save_file = paste(op_path, "/sdcut/sdcut_", each_cluster, ".eps", sep = ""))
  }else{
    print(paste("No influentials have been generated by standard deviation method for cluster", each_cluster, sep = " "))
  }
  
  #Random Permutation
  #Check the number of vertices in the cluster
  if(length(V(cluster_graph)) >= 5){  
    #Fit an ergm model for g1
    g_model <- ergm(cluster_graph_adj ~ edges)
    
    if(g_model$coef[[1]] != Inf){
      #Conditional simulation on outdegrees for the ergm model
      no_of_sim <- 1000
      g_model.sim <- simulate(g_model, constraints = ~odegrees, nsim = no_of_sim)
      
      #Store the indegrees of all the simulations in a vector
      in.deg.list <- vector(mode = "numeric", length = 0) 
      count <- 0
      
      for(sim_index in 1:no_of_sim){
        sim_mat <- as.matrix(g_model.sim[[sim_index]])
        sim_in_deg <- sna::degree(sim_mat, gmode = "digraph", cmode = "indegree")
        len <- length(sim_in_deg)
        
        count <- length(in.deg.list)
        min <- count + 1
        max <- count + len
        in.deg.list[min:max] <- sim_in_deg[1:len]
      }
      
      #histogram of the indegree alongwith the cut-off 
      q <- quantile(in.deg.list, prob = 0.95)
      hist(in.deg.list, 
           main = paste("Histogram of Random Permutation of cluster", each_cluster, 
                        sep = " "), xlab = "In-degree")
      abline(v = q, col = "navy")
      legend('topright', "alpha = 0.05", lty = 1, col = "navy",
             bty = 'n', cex = 1)
      
      #save plot as RandomHist<no>.png
      png(file = paste(op_path, "/randhist/randhist", each_cluster, ".png", sep = ""), 
          width = 1600, height = 1000)
      hist(in.deg.list, 
           main = paste("Histogram of Random Permutation of cluster", each_cluster, 
                        sep = " "), xlab = "In-degree")
      abline(v = q, col = "navy")
      legend('topright', "alpha = 0.05", lty = 1, col = "navy",
             bty = 'n', cex = 1)
      dev.off()
      
      cutoff <- q[[1]]
      data.in.deg <- sna::degree(cluster_graph_adj, gmode = "digraph", cmode = "indegree")
      
      names <- names(V(cluster_graph))
      names.list <- vector(mode = "character", length = 0)
      index.list <- vector(mode = "numeric", length = 0)
      x <- 1
      for(i in 1:length(data.in.deg)){
        if(data.in.deg[i] >= cutoff){
          names.list[x] <- names(V(cluster_graph)[i])
          index.list[x] <- i
          x = x + 1
        }
      }
      
      #color nodes of rand and plot these 
      color <- vector("numeric", length = length(V(cluster_graph)))
      index <- 1
      while(index <= length(names.list)){
        c <- which(names(V(cluster_graph)) == names.list[index])
        color[c] <- 1
        index <- index + 1
      }
      plot(cluster_graph, main = paste("Influentials under Random Permutation Method in cluster", each_cluster, sep = " "),
           layout = layout.fruchterman.reingold, 
           vertex.color = ifelse(color==1, "chocolate", "gold")) 
      legend('topright', c("Influentials", "Non-influentials"),
             col=c('chocolate', 'gold'), bty='n', cex=1, pch=19)
      visualise(object = cluster_graph, vcolor = ifelse(color==1, "chocolate", "gold"), 
                save_file = paste(op_path, "/randplot/randplot_", each_cluster, ".eps", sep = ""))
    }else{
      print(paste("The MLE coefficient is Inf. Cluster", each_cluster, 
                  "is saturated resulting in no influentials.", sep = " "))
    }
  }else{
    print(paste("The cluster", each_cluster, "has", length(names(V(cluster_graph))), "nodes.",
                sep = " "))
  }
}
}

#tkplot function
visualise <- function(object = " ", vcolor = "gold", save_file = " "){
  id <- tkplot(object, canvas.width = 1300, canvas.height = 680,
               vertex.color = vcolor, layout = layout.fruchterman.reingold)
  canvas <- tk_canvas(id)
  tkpostscript(canvas, file = save_file)
  tk_close(id)
}
