############################################
# F01 scores and hungarian helper function
############################################

# started https://stackoverflow.com/questions/8499361/easy-way-of-counting-precision-recall-and-f1-score-in-r
# initial organization from Chris Busch (http://cbuschnotes.tumblr.com/)
# substantially rewritten, renamed calcF1Scores -> calc_f1_scores, 
# added Hungarian alg. from 'clue' package
# input: prd are the predicted clusters, act is the "ground" truth
# changed code to transpose matrix if number cols is < nrow, have to rewrite the rest of the code
calc_f1_scores=function(prd,act){
  require(clue)
  #act and prd must be integers
  stopifnot(is.integer(act))
  stopifnot(is.integer(prd))
  
  
  tb <- table(act,prd)
  ifelse(ncol(tb) >= nrow(tb), tb <- tb, tb <- t(tb))
  ha <- solve_LSAP(tb, maximum=TRUE) # hungarian alg.
  
  df <- data.frame(act=as.integer(rownames(tb)),
                   prd=as.integer(colnames(tb)[ha[1:nrow(tb)]]))
  
  df$tp <- df$fp <- df$fn <- df$f1 <- NA
  for(i in seq_len(nrow(df))) {
    df$tp[i] <- sum( prd==df$prd[i] & act==df$act[i] )
    df$fp[i] <- sum( prd==df$prd[i] & act!=df$act[i] )
    df$fn[i] <- sum( prd!=df$prd[i] & act==df$act[i] )
    df$f1[i] <- with(df, (2*tp[i])/(2*tp[i]+fp[i]+fn[i]))
  }
  df
}



