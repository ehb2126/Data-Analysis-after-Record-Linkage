set.seed(005)
setwd("/home/b/benda002/rFiles/bagging/Results_twitter_data/mse")
source("/home/b/benda002/rFiles/bagging/source_functions/testing_fn.R")

library(ranger)
library(rpart)
library(dplyr)


exact2 <-read.csv("exact2.csv")

exact2$gun_pred <- factor(exact2$gun_pred)
 exact2$rep_pred <- factor(exact2$rep_pred)
 exact2$male_pred <- factor(exact2$male_pred)
 exact2$older_pred <- factor(exact2$older_pred)
 exact2$ppp10012 <- as.numeric(exact2$ppp10012)
 
 simDat <- exact2%>%select( ppp10012, tweet_count , older_pred, male_pred, gun_pred , rep_pred)
 
 formula_simDat <- ppp10012 ~ tweet_count + older_pred + male_pred + gun_pred + rep_pred
 




n <- nrow(simDat)



M_naive_mse_error <- numeric(8)
M_naive_rf_mse_error <- numeric(8)
M1_mse_error <- numeric(8)
M_rf1_mse_error <- numeric(8)
M2_mse_error <- numeric(8)

Op_M_naive_mse_error <- numeric(8)
Op_M_naive_rf_mse_error <- numeric(8)
Op_M1_mse_error <- numeric(8)
Op_M_rf1_mse_error <- numeric(8)
Op_M2_mse_error <- numeric(8)

Op_mean_M_naive_mse_error <- numeric(8)
Op_mean_M_naive_rf_mse_error <- numeric(8)
Op_mean_M1_mse_error <- numeric(8)
Op_mean_M_rf1_mse_error <- numeric(8)
Op_mean_M2_mse_error <- numeric(8)





mse_error_sum <- as.data.frame(list(M_naive_mse_error = M_naive_mse_error,
                                    M_naive_rf_mse_error = M_naive_rf_mse_error,
                                    M1_mse_error = M1_mse_error,
                                    M_rf1_mse_error = M_rf1_mse_error,
                                    M2_mse_error = M2_mse_error,
                                    
                                    Op_M_naive_mse_error=Op_M_naive_mse_error, 
                                    Op_M_naive_rf_mse_error = Op_M_naive_rf_mse_error,
                                    Op_M1_mse_error = Op_M1_mse_error,
                                    Op_M_rf1_mse_error = Op_M_rf1_mse_error,
                                    Op_M2_mse_error = Op_M2_mse_error,
                                    
                                    Op_mean_M_naive_mse_error = Op_mean_M_naive_mse_error,
                                    Op_mean_M_naive_rf_mse_error = Op_mean_M_naive_rf_mse_error,
                                    Op_mean_M1_mse_error = Op_mean_M1_mse_error,
                                    Op_mean_M_rf1_mse_error = Op_mean_M_rf1_mse_error,
                                    Op_mean_M2_mse_error = Op_mean_M2_mse_error
))



M_naive_sum <- M_naive_rf_sum <- M1_sum <- M_rf1_sum <- M2_sum <- M_op_naive_sum <- M_op_naive_rf_sum <-
  M1_op_sum <- M_op_rf1_sum <- M2_op_sum <- M_op_mean_naive_sum <- M_op_mean_naive_rf_sum <- 
  M1_op_mean_sum <- M_op_mean_rf1_sum <- M2_op_mean_sum <- as.data.frame(
    cbind( numeric(n), numeric(n), numeric(n), numeric(n),
         numeric(n),numeric(n), numeric(n),numeric(n)))



max_iter = 100

for(iter in 1:max_iter){
  
  
  
  simDat <- model.frame(formula_simDat, simDat)
  
  
  
  
  Kappa <- c(0,.10,.15,.2,.25,.3,.35,.4)
  Alpha <- rep(.5, n)
  wghts <- rep(1, n)
  
  
  
  for(i in 1:8){
    
    kappa <- Kappa[i]
  
    simDat_pdat <- simDat
    simDat_pdat$ppp10012 <- y.perm(simDat$ppp10012, kappa)
    
    
    running_time_Alpha_op <- system.time(Alpha_op <- optimal_Alpha(formula_simDat, data = simDat_pdat))
    running_time_Alpha_mean <- system.time( Alpha_mean <- rep(optimal_mean_alpha(formula_simDat, data = simDat_pdat), n))
    
    running_time_bagging <- system.time(  m_naive <- bag_trees(formula_simDat, data = simDat_pdat,  wghts = wghts))
    M_naive_mse_error[i] <- mse(m_naive, simDat$ppp10012)
    M_naive_sum[,i] <- M_naive_sum[,i] + m_naive
    
    running_time_ranger <- system.time( m_naive_rf <- ranger(formula_simDat, simDat_pdat)$predictions)
    M_naive_rf_mse_error[i] <- mse(m_naive_rf, simDat$ppp10012)
    M_naive_rf_sum[, i] = M_naive_rf_sum[, i] + m_naive_rf
    
    running_time_adj1 <- system.time(m1 <- adj_bag_trees1(formula_simDat, data  = simDat_pdat,  ntrees=100, wghts,
                                                          Alpha))
    M1_mse_error[i] <- mse(m1, simDat$ppp10012)
    M1_sum[, i] = M1_sum[, i] + m1
    
    running_time_adj_ranger <- system.time(  m_rf1 <- adj_rf1(formula_simDat, data =  simDat_pdat,  Alpha))
    M_rf1_mse_error[i] <- mse(m_rf1, simDat$ppp10012)
    M_rf1_sum[, i] = M_rf1_sum[, i] + m_rf1
    
    running_time_adj2 <- system.time( m2 <- adj_bag_trees2(formula_simDat, data =  simDat_pdat,  ntrees=100, wghts,
                                                           Alpha))
    M2_mse_error[i] <- mse(m2, simDat$ppp10012)
    M2_sum[, i] = M2_sum[, i] + m2
    
    
    # Alpha adjustment using Alpha = Optimal Alpha
    
    
    running_time_alpha_adj <- system.time(Op_m_naive <- adj_alpha_method(formula_simDat, data =  simDat_pdat, Alpha = Alpha_op, 
                                                                         mu_ytildex = m_naive))
    Op_M_naive_mse_error[i] <- mse( Op_m_naive, simDat$ppp10012)
    M_op_naive_sum[, i] = M_op_naive_sum[, i] + Op_m_naive
    
    
    Op_m_naive_rf <- adj_alpha_method(formula_simDat, data =  simDat_pdat, 
                                      Alpha = Alpha_op, mu_ytildex = m_naive_rf)
    Op_M_naive_rf_mse_error[i] <- mse( Op_m_naive_rf, simDat$ppp10012)
    M_op_naive_rf_sum[, i] =  M_op_naive_rf_sum[, i] + Op_m_naive_rf
    
    
    Op_m1 <- adj_alpha_method(formula_simDat,data =  simDat_pdat, Alpha = Alpha_op, mu_ytildex = m1)
    Op_M1_mse_error[i] <- mse( Op_m1, simDat$ppp10012)
    M1_op_sum[, i] = M1_op_sum[, i] + Op_m1
    
    Op_m_rf1 <- adj_alpha_method(formula_simDat, data =  simDat_pdat, 
                                 Alpha = Alpha_op, mu_ytildex = m_rf1)
    Op_M_rf1_mse_error[i] <- mse( Op_m_rf1, simDat$ppp10012)
    M_op_rf1_sum[, i] = M_op_rf1_sum[, i] + Op_m_rf1
    
    Op_m2 <- adj_alpha_method(formula_simDat, data =  simDat_pdat,
                              Alpha = Alpha_op, mu_ytildex = m2)
    Op_M2_mse_error[i] <- mse( Op_m2, simDat$ppp10012)
    M2_op_mean_sum[, i] = M2_op_mean_sum[, i] + Op_m2
    
    
    # Alpha adjustment method using alpha_mean_optimal
    
    
    
    Op_mean_m_naive <- adj_alpha_method(formula_simDat, data =  simDat_pdat, Alpha = Alpha_mean, 
                                        mu_ytildex = m_naive)
    Op_mean_M_naive_mse_error[i] <- mse( Op_mean_m_naive, simDat$ppp10012)
    M_op_mean_naive_sum[, i] = M_op_mean_naive_sum[, i] + Op_mean_m_naive
    
    Op_mean_m_naive_rf <- adj_alpha_method(formula_simDat, data =  simDat_pdat, 
                                           Alpha = Alpha_mean, mu_ytildex = m_naive_rf)
    Op_mean_M_naive_rf_mse_error[i] <- mse( Op_mean_m_naive_rf, simDat$ppp10012)
    M_op_mean_rf1_sum[, i] = M_op_mean_rf1_sum[, i] + Op_mean_m_naive_rf
    
    
    Op_mean_m1 <- adj_alpha_method(formula_simDat,data =  simDat_pdat, 
                                   Alpha = Alpha_mean, mu_ytildex = m1)
    Op_mean_M1_mse_error[i] <- mse( Op_mean_m1, simDat$ppp10012)
    M1_op_mean_sum[, i] = M1_op_mean_sum[, i] + Op_mean_m1 
    
    Op_mean_m_rf1 <- adj_alpha_method(formula_simDat, data =  simDat_pdat, 
                                      Alpha = Alpha_mean, mu_ytildex = m_rf1)
    Op_mean_M_rf1_mse_error[i] <- mse( Op_mean_m_rf1, simDat$ppp10012)
    M_op_mean_rf1_sum[, i] = M_op_mean_rf1_sum[, i] + Op_mean_m_rf1
    
    Op_mean_m2 <- adj_alpha_method(formula_simDat, data =  simDat_pdat,
                                   Alpha = Alpha_mean, mu_ytildex = m2)
    Op_mean_M2_mse_error[i] <- mse(Op_mean_m2, simDat$ppp10012)
    M2_op_mean_sum[, i] = M2_op_mean_sum[, i] + Op_mean_m2 
    
  }
  
  
  
  table_simDat_mse_error <-  as.data.frame(list(M_naive_mse_error = M_naive_mse_error,
                                                   M_naive_rf_mse_error = M_naive_rf_mse_error,
                                                   M1_mse_error = M1_mse_error,
                                                   M_rf1_mse_error = M_rf1_mse_error,
                                                   M2_mse_error = M2_mse_error,
                                                   
                                                   Op_M_naive_mse_error=Op_M_naive_mse_error, 
                                                   Op_M_naive_rf_mse_error = Op_M_naive_rf_mse_error,
                                                   Op_M1_mse_error = Op_M1_mse_error,
                                                   Op_M_rf1_mse_error = Op_M_rf1_mse_error,
                                                   Op_M2_mse_error = Op_M2_mse_error,
                                                   
                                                   Op_mean_M_naive_mse_error = Op_mean_M_naive_mse_error,
                                                   Op_mean_M_naive_rf_mse_error = Op_mean_M_naive_rf_mse_error,
                                                   Op_mean_M1_mse_error = Op_mean_M1_mse_error,
                                                   Op_mean_M_rf1_mse_error = Op_mean_M_rf1_mse_error,
                                                   Op_mean_M2_mse_error = Op_mean_M2_mse_error))
  
  
  
  
  
  
  
  mse_error_sum <- table_simDat_mse_error + mse_error_sum
  
  
}



table_simDat_mse_error <- mse_error_sum/max_iter

table_M_naive_predictions <- M_naive_sum/max_iter
table_M_naive_rf_predictions <- M_naive_rf_sum/max_iter
table_M1_predictions <- M1_sum/max_iter
table_M2_predictions <- M2_sum/max_iter
table_M_rf1_predictions <- M_rf1_sum/max_iter
table_M_op_naive_predictions <- M_op_naive_sum/max_iter
table_M_op_naive_rf_predictions <- M_op_naive_rf_sum/max_iter
table_M1_op_predictions <- M1_op_sum/max_iter
table_M2_op_predictions <- M2_op_sum/max_iter
table_M_op_rf1_predictions <- M_op_rf1_sum/max_iter
table_M_op_mean_naive_predictions <- M_op_mean_naive_sum/max_iter
table_M_op_mean_naive_rf_predictions <- M_op_mean_naive_rf_sum/max_iter
table_M1_op_mean_predictions <- M1_op_mean_sum/max_iter
table_M2_op_mean_predictions <- M2_op_mean_sum/max_iter
table_M_op_mean_rf1_predictions <- M_op_mean_rf1_sum/max_iter

write.csv(table_M_naive_predictions, "table_M_naive_predictions_50_replicates.csv")
write.csv(table_M_naive_rf_predictions,"table_M_naive_rf_predictions_50_replicates.csv")
write.csv(table_M1_predictions,"table_M1_predictions_50_replicates.csv")
write.csv(table_M2_predictions,"table_M2_predictions_50_replicates.csv")
write.csv(table_M_rf1_predictions,"table_M_rf1_predictions_50_replicates.csv")

write.csv(table_M_op_naive_predictions, "table_M_op_naive_predictions_50_replicates.csv")
write.csv(table_M_op_naive_rf_predictions,"table_M_op_naive_rf_predictions_50_replicates.csv")
write.csv(table_M1_op_predictions,"table_M1_op_predictions_50_replicates.csv")
write.csv(table_M2_op_predictions,"table_M2_op_predictions_50_replicates.csv")
write.csv(table_M_op_rf1_predictions,"table_M_op_rf1_predictions_50_replicates.csv")

write.csv(table_M_op_mean_naive_predictions, "table_M_op_mean_naive_predictions_50_replicates.csv")
write.csv(table_M_op_mean_naive_rf_predictions,"table_M_op_mean_naive_rf_predictions_50_replicates.csv")
write.csv(table_M1_op_mean_predictions,"table_M1_op_mean_predictions_50_replicates.csv")
write.csv(table_M2_op_mean_predictions,"table_M2_op_mean_predictions_50_replicates.csv")
write.csv(table_M_op_mean_rf1_predictions,"table_M_op_mean_rf1_predictions_50_replicates.csv")



table_simDat_mse_error <- apply(table_simDat_mse_error,1, round, 2)


running_time <- as.data.frame(list( running_time_Alpha_op = running_time_Alpha_op[3], 
                                    running_time_Alpha_mean = running_time_Alpha_mean[3],
                                    running_time_bagging = running_time_bagging[3],
                                    running_time_ranger = running_time_ranger[3],
                                    running_time_adj1 = running_time_adj1[3],
                                    running_time_adj_ranger = running_time_adj_ranger[3],
                                    running_time_adj2 = running_time_adj2[3],
                                    running_time_alpha_adj = running_time_alpha_adj[3]
))




write.csv(table_simDat_mse_error,"table_simDat_mse_error_50_new.csv")



write.csv(running_time, "running_time_simData_50_new.csv")
