wtcal<-function(ip_file = " ", file_sep = " ", isdirected = T, op_file = " "){

#Load package  
library(igraph)

#Read input file, convert to graph and plot
data<- read.csv(file = ip_file, header = FALSE, sep = file_sep)
graph<-simplify(graph.data.frame(data[,c(1,2)], directed = isdirected))

# plot(graph, main = "Network", vertex.color = "gold", 
#      layout = layout.fruchterman.reingold)
# 
# #save plot as network.png
# png(file = paste(strsplit(op_file, "[.]")[[1]][1], "_network.png", sep = ""),
#     height = 1600, width = 1000)
# plot(graph, main = "Network", vertex.color = "gold", 
#      layout = layout.fruchterman.reingold)
# dev.off()

#Get the edgelist, degree and vcount
el<-get.edgelist(graph)
deg<- degree(graph)
nk<-vcount(graph)

#wts contains weights assigned to each of the edges
wts<-vector("numeric", length = nrow(el))
for(i in 1:nrow(el)){
  nc1<- deg[[el[i,1]]]
  nc2<- deg[[el[i,2]]]
  wts[i]<- ((nc1+nc2)-1)/(nk-1)
}
wts<-round(wts,digits=4)
E(graph)$weight<-wts

#plot weighted network
# plot(graph, main = "Weighted Network", vertex.color = "gold", 
#      layout = layout.fruchterman.reingold, edge.label = E(graph)$weight)
# 
# #save plot as weighted.png
# png(file = paste(strsplit(op_file, "[.]")[[1]][1], "_weighted.png", sep = ""),
#     height = 1600, width = 1000)
# plot(graph, main = "Weighted Network", vertex.color = "gold", 
#      layout = layout.fruchterman.reingold, edge.label = E(graph)$weight)
# dev.off()

#op contains the edges along with the weights
op <- cbind(el, wts)

sink(op_file)
for(i in 1:nrow(op)){
  cat(op[i,1], " ", op[i,2], " ", op[i,3], "\n")
}
sink() 

return(wts)
}

