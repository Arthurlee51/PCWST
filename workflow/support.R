#loss function for evaluations
loss.func = function(Pihat, Pi_star, n){
  out = sum((Pihat-Pi_star)^2)/(n^2-n)
  return(out)
}

#Get Pihat using the Bradley-Terry model.
BT_est.func = function(Y){
  n <- nrow(Y)
  btdata<- btdata(Y)
  #Fit bt data, set a=1 to ensure that mle is used.
  bt_result <- btfit(btdata,1)
  # Extract the coefficients
  coefficients <- coef(bt_result)
  
  #Reorder the coefficients to compute Pihat_bt
  # Get the names of the players (e.g., "Player 327", "Player 297", ...)
  if(is.list(coefficients)){
    if(length(coefficients)==2){
      coefficients <- c(coefficients[[1]],coefficients[[2]])
    }
    if(length(coefficients)==1){
      coefficients <- coefficients[[1]]
    }
  }
  player_names <- names(coefficients)
  
  # Extract the player numbers using regular expressions
  player_numbers <- as.numeric(sub("Player ", "", player_names))
  
  # Create a data frame with player numbers and their coefficients
  coeff_df <- data.frame(
    player = player_numbers,
    coefficient = coefficients
  )
  
  missing_ind <- setdiff(1:n, coeff_df$player)
  if(length(missing_ind)>0){
    #Add coefficient to this and set to 0.
    missing_data <- data.frame(cbind(missing_ind , rep(0, length(missing_ind))))
    names(missing_data) <- c("player", "coefficient")
    coeff_df <- rbind(coeff_df, missing_data )
  }
  
  
  # Order the data frame by player numbers
  coeff_df_ordered <- coeff_df[order(coeff_df$player), ]
  #Check for missing numbers
  
  coeff_ordered <- coeff_df_ordered$coefficient
  
  #Compute Pihat_bt 
  Pihat_bt <- g.func(matrix( coeff_ordered ,nrow=n, ncol=n ) - matrix(  coeff_ordered ,nrow=n, ncol=n, byrow=TRUE))
  
  return(list(Pihat_bt =  Pihat_bt ,coeff_bt = coeff_ordered)  )
}

#Function to generate simulated data. 
#logit function
g.func = function(x){
  out=1/(1 + exp(-x))
  #Make sure it does not explode
  out[x> 700]=1
  out[x< -800]=0
  return(out)
}


#===================================================================================================================
#Functions relevant to data analysis only.
#Function to remove players who only win/loses.
rm_extreme_data.func= function(data){
  # Identify players who win all games or lose all games
  # Players who win all games: only appear in the 'Winner' column, not in 'Loser'
  winners_only <- setdiff(data$Winner, data$Loser)
  # Players who lose all games: only appear in the 'Loser' column, not in 'Winner'
  losers_only <- setdiff(data$Loser, data$Winner)
  # Combine the list of players to exclude
  players_to_exclude <- unique(c(winners_only, losers_only))
  # Step 2: Filter the dataset to exclude these players
  if(length(players_to_exclude)>0){
    filtered_data <- data[-unique(c(which(data$Loser %in% losers_only),which(data$Winner %in% winners_only))), ]
    #Step 2.5: Perform filtering again for those who now lose only after removing the previous losers
    while(( length(setdiff(filtered_data$Winner, filtered_data$Loser))>0 ) || ( length(setdiff(filtered_data$Loser,filtered_data$Winner ))>0 ) ){
      winners_only <- setdiff(filtered_data$Winner, filtered_data$Loser)
      losers_only <- setdiff(filtered_data$Loser, filtered_data$Winner)
      players_to_exclude <- unique(c(winners_only, losers_only))
      filtered_data <- filtered_data[-unique(c(which(filtered_data$Loser %in% losers_only),which(filtered_data$Winner %in% winners_only))), ]
    }
  }else{
    filtered_data <- data
  }
  
  return(filtered_data)
}

GetW.func = function(data, n_players,players){
  W <- matrix(0, nrow = n_players, ncol = n_players)
  
  # Step 6: Populate the win-loss matrix 
  for (i in 1:nrow(data)) {
    winner <- data$Winner[i]
    loser <- data$Loser[i]
    
    # Get the corresponding player IDs
    winner_id <- which(players== winner)
    loser_id <- which(players== loser)
    
    # Update the win-loss matrix
    W[winner_id, loser_id] <- W[winner_id, loser_id] + 1
  }
  return(W)
}
