lpcal <- function(ip_file = " ", file_sep = ",", thresh = 0.05, op_file = " "){

#Load packages
library(igraph)
library(tcltk)
library(tcltk2)

#Read input file, convert to graph
input <- read.csv(ip_file ,header = FALSE, sep = file_sep)
g <- simplify(graph.data.frame(input, directed = FALSE))

#finding clusters using fast greedy approach
fc  <- cluster_fast_greedy(g, weights = E(g)$weight)
plot(fc, g, main = "Clusters", vertex.color = "gold", 
     layout = layout.fruchterman.reingold)

#save plot as clusters.png
png(file = "/home/mohan/Desktop/Execution/output/lp/clusters.png",
    height = 1600, width = 1000)
plot(fc, g, main = "Clusters", vertex.color = "gold", 
     layout = layout.fruchterman.reingold)
dev.off()

#finding link existence and non-existence
u <- (vcount(g) * (vcount(g) - 1)) / 2
le <- ecount(g) / u
lne <- (u - ecount(g)) / u
omega <- le / lne

#threshold for final similarity score
threshold <- thresh

#get adjacency matrix
adj <- get.adjacency(g, sparse = FALSE)

#declarations for computing score
tcn <- 0
wcn <- 0
ocn <- 0
temp <- 0
c1 <- 0
c2 <- 0
temp1<- 0
temp2 <- 0
temp3 <- 0

sim <- matrix(NA, ncol = 3)
prob <- matrix(NA, ncol = 3)

#computing similarity scores
for(i in 1:nrow(adj)){
  for(j in 1:ncol(adj)){
    if(i > j){
      if(adj[i,j]==0){
        tcn <- 0
        c1 <- 0
        c2 <- 0
        cn1 <- 0
        cn2 <- 0
        wcn <- 0
        ocn <- 0
        s <- 0
      
        #Computing total no of common neighbors - cocitation() gives the total number of common neighbors
        #In place of cocitation() other functions can be used which represent different similarity measures
        tcn <- cocitation(g,i)[j]
      
        #compute cluster id
        c1 <- membership(fc)[V(g) == i][[1]]
        c2 <- membership(fc)[V(g) == j][[1]]
      
        #if they belong to the same cluster
        if(c1==c2){
          for(l in 1:vcount(g)){
            cn1 <- cocitation(g,i)[l]
            cn2 <- cocitation(g,j)[l]
            if(cn1 != 0 && cn2 != 0){
              if(membership(fc)[V(g)==l][[1]] == c1){
                wcn <- wcn + 1
              }
            }
          }
        }else{
          wcn <- 0
        }
        ocn <- -(tcn - wcn)
        if(wcn !=0 && ocn !=0){
          s <- (wcn / ocn) * omega
        }
        if(s>0){ 
          if(NA %in% sim){
            sim <- rbind(c(i, j, s))
          }else{
            sim <- rbind(sim, c(i, j, s))
          }
        }
      }
    }
  }
}

#Sorting in decreasing order based on simialrity score 
sim <- sim[order(-sim[,3]),]

#Filtering and writing the filtered & sorted edges into a file
for(k in 1:nrow(sim)){
  if(sim[k,3] > threshold){
    if(NA %in% prob){
      prob <- rbind(sim[k,])
    }else{
      prob <- rbind(prob, sim[k,])
    }
  }
}

sink(op_file)
for(i in 1:nrow(prob)){
 cat(prob[i, 1], ",", prob[i, 2], ",", prob[i,3], "\n")
}
sink()
 
#Plotting the edges of original network in BLACK and predicted edges in RED
od <- read.csv(op_file, fill=T, header=F, sep=",")

#Plotting origianal edges in BLACK
od[,1]=as.character(od[,1])
od[,2]=as.character(od[,2])
E(g)$color <- "black"
 
#Plotting new predicted edges in RED
for( i in 1:length(readLines(op_file))){ 
  g <- g+edge(od[,1][i-1],od[,2][i-1],color="red") 
}
 
#Plotting the entire network including the predicted edges
id <- tkplot(g, edge.label=E(g)$Prob, vertex.color = "gold", canvas.width = 1300, 
       canvas.height = 680, layout = layout.fruchterman.reingold)
canvas <- tk_canvas(id)
tkpostscript(canvas, file = paste(strsplit(op_file, "[.]")[[1]][1], "_prob.eps", 
                                  sep = ""))
tk_close(id)

return(prob)
}
