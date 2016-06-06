tsacal <- function(ip_file = " ", op_path = " "){

#load relevant packages
library(igraph)

#convert to igraph object with relevant headers
edges <- read.csv(ip_file, header = FALSE, sep = " ")
edges <- edges[order(edges[,3]),]
g <- graph.data.frame(edges[,c(1,2)], directed = FALSE)
E(g)$time <- edges[,3]
time <- edges[,3]

#generate color palette
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
YlOrBr.Lab <- colorRampPalette(YlOrBr, space = "Lab")
vcolor <- rev(YlOrBr.Lab(vcount(g)))

#where to begin
minimum <- min(time)
for(i in 1:length(time)){
  if(time[i] > minimum){
    start <- time[i]
    break
  }
}

#set layout and time
E(g)$weight <- ifelse(E(g)$time < start, 1, 0)
layout.old <- layout.fruchterman.reingold(g,params=list(weights=E(g)$weight))

#where to end
end <- max(time)

#compute increment
time_diff <- vector("numeric", length(time))
for(i in 2:length(time)){
  time_diff[i-1] <- time[i] - time[i-1]
}
occur <- table(time_diff)
occur <- sort(occur, decreasing = TRUE)
if(as.integer(names(occur[1])) == 0){
  incr <- as.integer(names(occur[2]))
}

#plot and save in png
png(file = paste(op_path, "/example%03d.png", sep = ""), width=1600,height=900)
for(ti in seq(start,end,incr)){
  E(g)$weight <- ifelse(E(g)$time < ti,1,0) 
  E(g)$color <- ifelse(E(g)$time < ti,"gray",rgb(0,0,0,0))
  V(g)$color <- ifelse(graph.strength(g)==0,rgb(0,0,0,0),vcolor)
  layout.new <- layout.fruchterman.reingold(g,params=list(niter=10,start=layout.old,weights=E(g)$weight,maxdelta=1))
  plot(g,layout=layout.new,vertex.size=1+2*log(graph.strength(g)),vertex.frame.color=V(g)$color,edge.width=1.5,asp=9/16,margin=-0.15)
  layout.old <- layout.new 
}
dev.off()

#run the movie
system("xdg-open /home/mohan/Desktop/Execution/output/tsa/dynamic/timer1.mov")
}

