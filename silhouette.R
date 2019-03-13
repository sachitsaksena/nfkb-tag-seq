silhouette <- function(pca) {
  #################
  # calculate within averages
  ##################
  # we just care about pc1 and pc2
  pca_cull <- as.data.frame(pca$x[,1:2])
  dist_mat <- as.data.frame(as.matrix(dist(pca_cull)))
  uniques <- unique(substr(colnames(dist_mat), 1, 1))
  
  # assign within average dataframe
  w = data.frame(name=colnames(dist_mat))
  within_averages <- list()
  for(i in seq(1, dim(dist_mat)[1], 4)){
    within_average=0
    tot=0
    # calculate the average of the 1st point in the set
    average1 = (dist_mat[i, i+1] +
                  dist_mat[i, i+2] +
                  dist_mat[i, i+3]) / 3
    # calculate the average of the 2nd point in the set
    average2 = (dist_mat[i, i] +
                  dist_mat[i, i+2] +
                  dist_mat[i, i+3]) / 3
    # calculate the average of the 3rd point in the set
    average3 = (dist_mat[i, i]+
                  dist_mat[i,i+1] +
                  dist_mat[i,i+3]) / 3
    # calculate the average of the 4th point in the set
    average4 = (dist_mat[i,i] + 
                  dist_mat[i, i+1] +
                  dist_mat[i, i+2]) / 3
    within_averages <- c(within_averages, average1, average2, average3, average4)
    }
  # append within averages
  w %>% mutate(within=within_averages) -> w
  w <- column_to_rownames(w, var="name")
  w$within <- as.numeric(w$within)
  # write within averages to file
  
  #################
  # calculate between averages
  ################
  
  b = data.frame(row.names = uniques)
  for (i in 1:dim(dist_mat)[1]){
    between_averages = list()
    current_name <- colnames(dist_mat)[i]
    # iterating through every other cluster
    for (j in seq(1, dim(dist_mat)[1], 4)){
      tot=0
      average=0
      # iterating through each point in a selected cluster
      for (k in 0:3){
        tot = tot + dist_mat[i, j+k]
      }
      # calculate average
      average=tot/4
      # append average to list of all the average distances to clusters
      between_averages <- c(between_averages, average)
    }
    # attaching the average distance of each point to each cluster in the dataset
    b[[current_name]] = with(b, between_averages)
    b[[current_name]] <- as.numeric(b[[current_name]])
  }
  # b <- as.data.frame(t(b))
  
  #####################
  # calculate silhouette values 
  ######################
  
  
  b_is = data.frame(row.names=colnames(dist_mat))
  # let's get our b(i)'s
  for (i in 1:dim(b)[2]){
    clust = substr(colnames(b[i]), 1,1)
    cs = array()
    for (j in 1:length(rownames(b[i]))){
      if (!as.character(clust)==rownames(b)[j]){
        cs <- c(cs, as.numeric(b[j,i]))
      }
    }
    b_i <- min(cs, na.rm=TRUE)
    b_is[[colnames(b[i])]] = b_i
  }
  # clean up the data frame 
  b_is <- unique(b_is)
  rownames(b_is) <- "b_i"
  # transpose the matrix
  b_is <- as.data.frame(t(b_is))
  
  # merge the a_i matrix and the b_i
  ab <- merge(w, b_is, by=0)
  ab <- column_to_rownames(ab, var="Row.names")
  
  ###############
  # calculate silhouette
  ##############
  s = data.frame(row.names=colnames(dist_mat))
  for (i in 1:dim(ab)[1]){
    s[[rownames(ab)[i]]] = (ab[i,2]-ab[i,1]) / max(ab[i,1], ab[i,2])
  }
  s <- as.data.frame(t(unique(s)))
  colnames(s) <- "silhouette"
  
  # calculate averages 
  s_avg = data.frame(row.names=uniques)
  cluster_averages = list()
  average = 0
  for (i in seq(1, dim(s)[1], 4)){
    average = (s[i, 1] + 
              s[i+1, 1] + 
              s[i+2, 1] +
              s[i+3, 1])/4
    cluster_averages <- c(cluster_averages, as.numeric(average))
  }
  
  # clean up data frame
  s_avg %>% rownames_to_column() %>% 
    mutate(cluster_averages) -> s_avg
  s_avg <- column_to_rownames(s_avg, var="rowname")
  
  return(s_avg)
}

