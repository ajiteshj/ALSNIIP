data <- read.csv("/home/mohan/Desktop/ALSNIIP/datasets/ca-cit-HepPh/out.ca-cit-HepPh", 
                 header = FALSE, sep = " ")
data <- data[-c(1:3),]
wt_data <- data[,c(1,2)]
sink("/home/mohan/Desktop/ALSNIIP/results/wt/input/hep-ph_from_to.txt")
for(i in 1:nrow(wt_data)){
  wt_put <- paste(wt_data[i,1], wt_data[i,2], sep = " ")
  cat(wt_put, "\n")
}
sink()

tsa_data <- data[,c(1,2,4)]
sink("/home/mohan/Desktop/ALSNIIP/results/tsa/input/hep-ph_from_to_time.txt")
for(i in 1:nrow(tsa_data)){
  tsa_put <- paste(tsa_data[i,1], tsa_data[i,2], tsa_data[i,3], sep = " ")
  cat(tsa_put, "\n")
}
sink()
