# addis controls the online FDR rate. It calculates the alpha value that is used
# to compare with the current p-value.
#
# alpha - The overall FDR rate. This is the FDR input from the main MRPC 
# funciton.
#
# tau - A number between 0 and 1. This value is used to determine if a p-value
# will be considered for testing. For example, if a p-value is greater than tau
# then it is discarded and no test will be performed.
#
# lambda - A number between 0 and 1. This valued is used to determine if a
# p-value is a candidate for rejection. For example, if a p-value is smaller
# than lambda then it can be rejected when testing the hypothesis (if the 
# p-value is smaller than alphai).
#
# iter - An integer. This is the current iteration of the MRPC algorithm.
#
# w0 - The initial "wealth" the addis algorithm starts out with. This is
# determined by the equation tau * lambda * alpha/2.
#
# pval - A vector of p-values up to the current iteration.
#
# alphai - A number between 0 and 1. This is the value that will be compared to
# the current p-value to determine whether or not it will be rejected.
#
# gammai - A vector containing the values of the infinite sequence sum 1 to inf
# 0.4374901658/(iteration^(1.6)).
#
# kappai - A vector containing the iterations at which each rejection occurs.
#
# kappai_star - A vector. Each element of this vector is the sum of the Si
# vector up to the index at which each rejection occurs.
#
# K - A scalar that denotes the number of rejections so far. This is the sum of
# the Ri (R in MRPC) vector or the length of the kappai vector.
#
# Ci - A vector indicating whether or not a p-value is a candidate for being
# rejected.
#
# Si - A vector indicating whether or not a p-value is selected for testing. In
# other words, whether or not the p-value was discarded.
#
# Ri - A vector indicating whether or not a hypothesis was rejected.
#
# Ci_plus - A vector of the number of times each kappai value was counted when
# calculating the alphai value.
#
# Ci_sum - A scalar. This is the sum of the Ci vector at each iteration.
#
# Si_sum - A scalar. This is the sum of the Si vector at each iteration.
#
# gammai_sum - A scalar. The sum of the gammai vector is calculated using the
# Si_sum and therefore has to be recalculated each time the Si_sum changes. This
# value is used in calculating the alphai_hat value at each iteration.
#
# normalizer - The value that ensures the gammai vector sums to one.
#
# exponent - The exponent of the p-series used to calculate each value of the
# gammai vector.
#
# Returns the alphai value calculated for the current iteration and the updated
# values for gammai, K, kappai, Ri, Si, Ci, Ci_plus, Ci_sum, and Si_sum.
#
addis <- function (alpha,
                   tau,
                   lambda,
                   iter,
                   w0,
                   pval,
                   alphai,
                   gammai,
                   kappai,
                   kappai_star,
                   K,
                   Ci,
                   Si,
                   Ri,
                   Ci_plus,
                   Ci_sum,
                   Si_sum,
                   gammai_sum,
                   normalizer,
                   exponent) {
  
  # Update all vectors and sums ----------
  
  # Calculate gamma for the current iteration.
  gammai[iter] <- normalizer / iter^exponent
  
  # If the current hypothesis was rejected then add the iteration # to the 
  # kappai vector.
  if (Ri[iter - 1] == 1) {
    
    # Increment K.
    K <- K + 1
    
    # Add the iteration # to the kappai vector.
    kappai[K] <- iter - 1
    
  }
  
  # Determine if the p-value will be tested.
  Si[iter - 1] <- pval[iter - 1] <= tau
  
  # Calculate the new sum of the Si vector.
  Si_sum <- Si_sum + Si[iter - 1]
  
  # Determine if the p-value is a candidate for rejection.
  Ci[iter - 1] <- pval[iter - 1] <= tau * lambda
  
  # Check if the current p-value is a candidate for being tested.
  if (Ci[iter - 1] == 1) {
    
    # Increase Ci_sum by one.
    Ci_sum <- Ci_sum + 1
    
  }
  
  # Calculate alpha depending on the number of discoveries ----------
  
  # Check how large K is and calculate alphai accordingly.
  if (K > 1) {
    
    # Check if the Kth element of kappai_star is a zero. If it is then calculate
    # the sum of the Si vector up to the index/iteration of the Kth rejection.
    if (kappai_star[K] == 0) {
      
      # Sum the number of selected p-values.
      kappai_star[K] <- sum(Si[seq_len(kappai[K])])
      
    }
    
    # Check if the value from the last iteration of the Si vector is a one. In
    # other words, was the p-value not discarded (pval[iter - 1] <= tau)?
    if (Si[iter - 1] == 1) {
      
      # Create a vector that is the length of the number of rejected hypotheses
      # minus one.
      Kseq <- seq_len(K - 1)
      
      # Check if the value from the last iteration of the Ci vector is a one. In
      # other words, is the p-value a candidate for rejection (pval[iter - 1]
      # <= tau * lambda)?
      if (Ci[iter - 1] == 1) {
        
        # Update the number of candidate p-values.
        Ci_plus[Kseq] <- Ci_plus[Kseq] + Ci[iter - 1]
        
      }
      
      # Sum the number of candidate for rejection p-values for the Kth element.
      Ci_plus[K] <- sum(Ci[seq(from = kappai[K] + 1,
                               to = max(iter - 1, kappai[K] + 1))])
      
      # Check the value of alphai from the last iteration and the p-value from
      # the current iteration. If alphai[iter - 1] = tau * lambda and
      # pval[iter] <= tau * lambda then there is no need to calculate the
      # gammai_sum because the alphai_hat value (which is calculated with the
      # gammai_sum value) will be larger than tau * lambda.
      if (!(alphai[iter-1] == tau * lambda && pval[iter-1] <= tau * lambda)) {
        
        # Sum the gammai sequence.
        gammai_sum <- sum(gammai[Si_sum
                                 - kappai_star[Kseq]
                                 - Ci_plus[Kseq] + 1])
        
        # Update the sum of the gammai sequence after the Kth element of the
        # Ci_plus vector has been updated. The Ci_plus vector is updated before
        # the if statement instead of after the gammai_sum value is initialized
        # to save time when the alphai value is larger than tau * lambda.
        gammai_sum <- (gammai_sum + gammai[Si_sum
                                           - kappai_star[K]
                                           - Ci_plus[K] + 1]
                       - gammai[Si_sum
                                - kappai_star[1]
                                - Ci_plus[1] + 1])
        
      }
      
    }
    
    # Check the value of alphai from the last iteration and the p-value from the
    # current iteration. If alphai[iter - 1] = tau * lambda and
    # pval[iter] <= tau * lambda then there is no need to calculate alphai_hat
    # because this value will be larger than tau * lambda.
    if (!(alphai[iter-1] == tau * lambda && pval[iter-1] <= tau * lambda)) {
      
      # Calculate the potential alphai value.
      alphai_hat <- (w0 * gammai[Si_sum
                                 - Ci_sum + 1]
                     + (tau * (1 - lambda) * alpha - w0)
                     * gammai[Si_sum
                              - kappai_star[1]
                              - Ci_plus[1] + 1]
                     + tau * (1 - lambda) * alpha * gammai_sum)
      
      # The following code will be run when alphai[iter - 1] = tau * lambda and
      # pval[iter] <= tau * lambda.
    } else {
      
      # Set alphai_hat to a value that will be larger than tau * lambda for any
      # values of tau and lambda. When alphai[iter - 1] = tau * lambda and
      # pval[iter] <= tau the gammai_sum will not be calculated to save time and
      # alphai_hat will need to be set to a large value.
      alphai_hat <- 1
      
    }
    
  } else if (K == 1) {
    
    # Check if the Kth element of kappai_star is a zero. If it is then calculate
    # the sum of the Si vector up to the index/iteration of the Kth rejection.
    if (kappai_star[K] == 0) {
      
      # Sum the number of selected p-values.
      kappai_star[K] <- sum(Si[seq_len(kappai[K])])
      
    }
    
    # Sum the number of candidate p-values.
    Ci_plus[1] <- sum(Ci[seq(from = kappai[K] + 1,
                             to = max(iter - 1, kappai[K] + 1))])
    
    # Calculate the potential alphai value.
    alphai_hat <- (w0 * gammai[Si_sum - Ci_sum + 1]
                   + (tau * (1 - lambda) * alpha - w0)
                   * gammai[Si_sum - kappai_star[K] - Ci_plus[1] + 1])
    
  } else {
    
    # Calculate the potential alphai value
    alphai_hat <- w0 * gammai[Si_sum - Ci_sum + 1]
    
  }
  
  # Choose the alpha value for the current iteration. This will be compared to
  # the p-value from this iteration at the beginning of the next iteration.
  alphai[iter] <- min(tau * lambda, alphai_hat)
  
  # Determine whether or not to reject the current hypothesis.
  Ri[iter] <- pval[iter] <= alphai[iter]
  
  # K will be zero until the first rejected hypothesis. We cannot subset the
  # kappai vector with 0 (i.e., kappai[[0]]). Until the first rejection we will
  # manually return the first element of the kappai vector.
  if (K == 0) {
    
    # Return the updated values for alphai, gammai, kappai, K, Ri, Si, Ci,
    # Ci_plus, Ci_sum, Si_sum, and kappai_star.
    return (list(alphai[iter],
                 gammai[iter],
                 K,
                 kappai[1],
                 Ri[iter],
                 Si[iter - 1],
                 Ci[iter - 1],
                 Ci_plus[1],
                 Ci_sum,
                 Si_sum,
                 kappai_star[1],
                 gammai_sum))
    
  } else {
    
    # Return the updated values for alphai, gammai, kappai, K, Ri, Si, Ci,
    # Ci_plus, Ci_sum, Si_sum, and kappai_star.
    return (list(alphai[iter],
                 gammai[iter],
                 K,
                 kappai[K],
                 Ri[iter],
                 Si[iter - 1],
                 Ci[iter - 1],
                 Ci_plus[1:K],
                 Ci_sum,
                 Si_sum,
                 kappai_star[K],
                 gammai_sum))
    
  }
  
}
