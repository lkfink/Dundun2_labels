# Custom functions for dundun w & w/out labels experiments
# LKF 2022
# lauren.fink@ae.mpg.de

# Function to test whether participants placed more stimuli into the left vs. right side
# input should be data frame with cols 1:30 (stimuli) and rows participants
# 1s represent left and 2s right
checkSideBias = function(d) 
{
  # get total number of stimuli placed on left vs right
  sum(d==1)
  sum(d==2)
  
  # Now use 1s for right
  d = d-1
  d$right = rowSums(d)
  d$left = 30-d$right # NOTE hard-coding 30 (num of stims in experiment)
  
  # t test for side bias
  tres = t.test(d$left, d$right, paired=TRUE)
  tres
}


# Function to plot heatmaps of raw data
plotHeatmap = function(d, exp)
{
  # Choose color representation based on experiment number
  # if(exp == 1){ourcolors = c("red", "blue")}
  # if(exp == 2){ourcolors = c(c("black", "grey"))}
  ourcolors = c(c("black", "grey"))
  
  # Create heatmap
  heatmap.2(data.matrix(d), 
            scale="none", 
            xlab="", ylab="", 
            col=ourcolors,
            dendrogram="none", #"col",
            Rowv=FALSE,
            trace="none",
            key=FALSE,
            density.info = "none",
            cexCol = 1,
            labRow = FALSE,
            Colv=FALSE
            # TODO figure out layout 
  )
}




# Function to create distance matrix based on the number of times each stimulus paired with each other stimulus
# Input = n x 30 data frame with participant sides for each stimulus (1 = left; 2 = right)
# Output = dissimilarity matrix

# first create empty data frame
get_groupings = function(d)
{
  # create empty matrix that will hold our results
  # keep same row and column names as original matrix
  df = matrix(nrow=30, ncol=30)
  rownames(df) = colnames(d)
  colnames(df) <- colnames(d)
  
  # find the number of times each row/col equal to each other
  for (ik in 1:30){
    for (ij in 1:30){
      if (ik == ij){
        df[ik,ij] = nrow(d) # max number of pairings
      }
      else {
        numpairs = length(d[d[,ik] == d[,ij], 1])
        df[ik,ij] = numpairs
      }
    }
  }
  
  # normalize matrix by dividing by max number of times stimuli could have been together
  df2 = df/nrow(d)
  
  # return matrix to user
  return(df2)
}



# Function to create dendograms
# Input: dissimilarity matric
# Output: hierarchical clustering results
get_dendro = function(dm){
  
  # Compute clusters
  hc.res <- eclust(dm, "hclust", hc_metric = "euclidean", 
                   hc_method = "ward.D2", graph = TRUE, nboot=500)
  
  # Create silhouette plot and calculate gap statistic
  fviz_silhouette(hc.res)
  fviz_gap_stat(hc.res$gap_stat)
  
  # Output number of clusters
  hc.res$nbclust
  
  return(hc.res)
  
}




# # Function for modeling acoustic features in relation to pca dimensions
# # Input: results of pca and table of acoustic features
# # Output: linear model results 
# getModelRes = function(pca.res, af){
#  
#   # Get all data into one table
#   pcadata = data.frame(pca.res$data$name, pca.res$data$x, pca.res$data$y, pca.res$data$coord) 
#   forlm = inner_join(pcadata, af, by = c("pca.res.data.name" = "stim"))
#   
#   # dim 1
#   # regresx = lm(pca.res.data.x ~ intensity_mean + intensity_meanDiff_betweenNotes + timbre_mean + IOI_mean + ratio_mean + pulseClarity + ampModSpectrum_peak, 
#   #              data = forlm)
#   
#   regresx = lm(pca.res.data.x ~ intensity_mean + intensity_meanDiff_betweenNotes + timbre_mean + IOI_mean + IOI_meanDiff_betweenNotes + ratio_mean + pulseClarity + ampModSpectrum_peak, 
#                data = forlm)
#   
#   print(summary(regresx))
#   print(car::vif(regresx))
#   print(tab_model(regresx, show.fstat=TRUE, show.obs=TRUE, pred.labels =
#               c("(Intercept)",  "intensity (mean)", "intensity (difference)", "timbre (mean)", "IOI (mean)", "IOI (difference)", "ratio (mean)", "pulse clarity", "Amp. mod. spectrum peak"),
#             dv.labels = "Dim. 1 position"))
# 
#   
#   # #dim 2
#   # regresy = lm(pca.res.data.y ~ intensity_mean + intensity_meanDiff_betweenNotes + timbre_mean + IOI_mean + IOI_meanDiff_betweenNotes + ratio_mean + pulseClarity + ampModSpectrum_peak,
#   #              data = forlm)
#   # print(summary(regresy))
#   # print(car::vif(regresy))
#   # print(tab_model(regresy, show.fstat=TRUE, show.obs=TRUE, pred.labels =
#   #             c("(Intercept)", "intensity (mean)", "intensity (difference)", "timbre (mean)", "IOI (mean)", "IOI (difference)", "ratio (mean)", "pulse clarity", "Amp. mod. spectrum peak"),
#   #           dv.labels = "Dim. 2 position"))
# 
# 
#   # return res x and y TODO
#   
# }


# Function for modeling acoustic features in relation to behavioural and acoustic pca dimensions
# Input: results of both pcas
# Output: linear model results 
getModelRes = function(pca.res, newaf){
  set.seed(123) # for reproducability
  
  # Get all data into one table
  pcadata = data.frame(pca.res$data$name, pca.res$data$x, pca.res$data$y, pca.res$data$coord) 
  forlm = inner_join(pcadata, newaf, by = c("pca.res.data.name" = "stim"))
  
  
  # Model to predict X (Dim 1) position
  regresx = lm(pca.res.data.x ~ PC1 + PC2 + PC3 + PC4, 
               data = forlm)
  
  print(summary(regresx))
  print(car::vif(regresx))
  print(tab_model(regresx, show.fstat=TRUE, show.obs=TRUE, pred.labels =
                    c("(Intercept)",  "PC1", "PC2", "PC3", "PC4"),
                  dv.labels = "Dim. 1 position"))
  
  # Check for overfit using cross-validation
  # Define training control
  train.control <- trainControl(method = "repeatedcv", 
                                number = 5, repeats = 10)
  # Train the model
  model_val_x <- train(pca.res.data.x~PC1 + PC2 + PC3 + PC4, data = forlm, method = "lm",
                 trControl = train.control)
  # Summarize the results
  print(model_val_x)
  
  
  # Model to predict Y (Dim 2) position
  regresy = lm(pca.res.data.y ~ PC1 + PC2 + PC3 + PC4, 
               data = forlm)
  
  print(summary(regresy))
  print(car::vif(regresy))
  print(tab_model(regresy, show.fstat=TRUE, show.obs=TRUE, pred.labels =
                    c("(Intercept)",  "PC1", "PC2", "PC3", "PC4"),
                  dv.labels = "Dim. 2 position"))
  
  # Check for overfit using cross-validation
  # Train the model
  model_val_y <- train(pca.res.data.y~PC1 + PC2 + PC3 + PC4, data = forlm, method = "lm",
                       trControl = train.control)
  print(model_val_y)
  
  # return model outputs in case we want to do anything else with them
  return(list(regresx, regresy, model_val_x, model_val_y))
}




