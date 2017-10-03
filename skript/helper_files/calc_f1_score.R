############################################
# F01 scores and hungarian helper function
############################################

# started https://stackoverflow.com/questions/8499361/easy-way-of-counting-precision-recall-and-f1-score-in-r
# initial organization from Chris Busch (http://cbuschnotes.tumblr.com/)
# substantially rewritten,  added Hungarian alg. from 'clue' package by Mark Robinson
# renamed calcF1Scores -> calc_f1_scores,
# changed code to transpose matrix; cluster are rows, ground truth is column


# # arguments:
# - cluster: cluster labels from algorithm
# - labels: true cluster labels
# (for both arguments: length = number of cells; names = cluster labels (integers))

calc_f1_score <- function(labels,cluster){
  {
    require(clue)
    #act and prd must be integers
    stopifnot(is.integer(labels))
    stopifnot(is.integer(cluster))
    
    # creat table with labels in rows and cluster in column
    tbl <- table(cluster=cluster, label=labels)
    # add additional column with max if number of column <= number of rows
    ifelse(ncol(tbl) >= nrow(tbl), tbl <- tbl, tbl <- cbind(tbl,"x"= array(max(tbl), dim=nrow(tbl))) )
    # compute hungarian assignement
    ha <- solve_LSAP(tbl, maximum=TRUE) # hungarian alg.
    df <- data.frame(cluster=as.integer(rownames(tbl)),
                     labels=as.integer(colnames(tbl)[ha[1:nrow(tbl)]]))
    df$tp <- df$fp <- df$fn <- df$f1 <- NA
    
    # compute tp, fp,fn and f1 score
    for ( i in seq_len(nrow(df)) ) {
      df$tp[i] <- sum( labels==df$labels[i] & cluster==df$cluster[i] )
      df$fp[i] <- sum( labels==df$labels[i] & cluster!=df$cluster[i] )
      df$fn[i] <- sum( labels!=df$labels[i] & cluster==df$cluster[i] )
      df$f1[i] <- with(df, (2*tp[i])/(2*tp[i]+fp[i]+fn[i]))  }
  }
  return(df)
}
