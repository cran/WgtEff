#' @title Calculate DEFF
#' @description Calculates design effect (DEFF)
#' @param x = weights vector (name of weights column)
#' @return Design effect (DEFF)
#' @examples DEFF(testweights$weights_column)
#' @references Design effect (DEFF) due to weighting => n * (sum(x^2) / sum(x)^2)
#' @export DEFF
DEFF=function(x){
  stopifnot(length(x)>2) #stops for bad inputs
  thedeff=((length(x)*sum(x^2))/(sum(x)^2)) #DEFF
  cat("Design effect (DEFF) =",thedeff)
}

#' @title Calculate weighting loss
#' @description Calculates weighting loss
#' @param x = weights vector (name of weights column)
#' @return Weighting loss
#' @examples WTGLOSS(testweights$weights_column)
#' @references Weighting loss => DEFF-1
#' @export WTGLOSS
WTGLOSS=function(x){
  stopifnot(length(x)>2) #stops for bad inputs
  thedeff=((length(x)*sum(x^2))/(sum(x)^2)) #WTGLOSS
  cat("Weighting loss =",(thedeff-1))
}

#' @title Calculate DEFT
#' @description Calculates root design effect (DEFT) 
#' @param x = weights vector (name of weights column)
#' @return Root design effect (DEFT) 
#' @examples DEFT(testweights$weights_column)
#' @references Root design effect (DEFT) => square root of DEFF
#' @export DEFT
DEFT=function(x){
  stopifnot(length(x)>2) #stops for bad inputs
  thedeft=sqrt(((length(x)*sum(x^2))/(sum(x)^2))) #DEFT
  cat("Root design effect (DEFT) =",thedeft)
}

#' @title Calculate ESS
#' @description Calculates effective sample size (ESS)
#' @param x = weights vector (name of weights column)
#' @return Effective sample size (ESS)
#' @examples ESS(testweights$weights_column)
#' @references Effective sample size (ESS) => sum(x)^2 / sum(x^2)
#' @export ESS
ESS=function(x){
  stopifnot(length(x)>2) #stops for bad inputs
  theess=((sum(x)^2)/(sum(x^2))) #ESS
  cat("Effective sample size (ESS) =",theess)
}

#' @title Calculate MOE
#' @description Calculates weighted margin of error (MOE)
#' @import stats
#' @param p = percentage for which MOE is calculated (optional, default is p = 50)
#' @param conf = level of confidence (optional, default is conf = 95)
#' @param N = population size (optional, used for finite population correction)
#' @param wtcol = Weights vector (name of weights column)
#' @return Weighted margin of error (MOE)
#' @examples MOE(N=3000, wtcol=testweights$weights_column)
#' @references Weighted margin of error (MOE) => unweighted MOE * DEFT
#' @export MOE
MOE=function(p=50, conf=95, N, wtcol){
  stopifnot(p<=100, p>=0, conf>=0, conf<=100, 
    N>=length(wtcol), length(wtcol)>2) #stops for bad inputs 
  if(missing(N)){fpc=1
  }
  if(!missing(N)){fpc=sqrt((N-length(wtcol))/(N-1))
  } #finite population correction
  tt=((100-((100-conf)/2))/100) #converts from two-tailed
  themoe=qnorm(tt)*sqrt(((p/100)*(1-(p/100)))/length(wtcol))*
    fpc*sqrt((((length(wtcol)*
    sum(wtcol^2))/(sum(wtcol)^2)))) #MoE formula
  cat("Margin of error (MoE) = +/-",themoe,"\n")
  cat("\n")
  cat("  **MoE calculated for percentages of",p,"percent, 
      at a",conf,"percent level of confidence")
}

#' @title Calculate Full Statistics
#' @description Calculates DEFF, weighting loss, DEFT, ESS, and MOE
#' @import stats
#' @param p = percentage for which MOE is calculated (optional, default is p = 50)
#' @param conf = level of confidence (optional, default is conf = 95)
#' @param N = population size (optional, used for finite population correction)
#' @param wtcol = Weights vector (name of weights column)
#' @return DEFF, weighting loss, DEFT, ESS, and MOE
#' @examples FULL(N=3000, wtcol=testweights$weights_column)
#' @export FULL
FULL=function(p=50, conf=95, N, wtcol){
  stopifnot(p<=100, p>=0, conf>=0, conf<=100, 
    N>=length(wtcol), length(wtcol)>2) #stops for bad inputs
  if(missing(N)){fpc=1
  }
  if(!missing(N)){fpc=sqrt((N-length(wtcol))/(N-1))
  } #finite population correction
  thedeff=(length(wtcol)*sum(wtcol^2))/(sum(wtcol)^2) #DEFF
  thedeft=sqrt(((length(wtcol)*sum(wtcol^2))/(sum(wtcol)^2))) #DEFT
  theess=((sum(wtcol)^2)/(sum(wtcol^2))) #ESS
  tt=((100-((100-conf)/2))/100) #converts from two-tailed
  themoe=qnorm(tt)*sqrt(((p/100)*(1-(p/100)))/length(wtcol))*
    fpc*thedeft #MoE formula
  #Print results
  cat("\n")
  cat("+++++++++++++++","\n")
  cat("\n")
  cat("For unweighted data:","\n")
  cat("\n")
  cat("- Margin of error (MoE) = +/-",themoe/thedeft,"\n")
  cat("- Sample size (n) =",length(wtcol),"\n")
  cat("\n")
  cat("For weighted data:","\n")
  cat("\n")
  cat("- Margin of error (MoE) = +/-",themoe,"\n")
  cat("- Effective sample size (ESS) =",theess,"\n")
  cat("- Design effect (DEFF) =",thedeff,"\n")
  cat("- Root design effect (DEFT) =",thedeft,"\n")
  cat("- Weighting loss =",thedeff-1,"\n")
  cat("\n")
  cat("  **MoE calculated for percentages of",p,"percent, 
      at a",conf,"percent level of confidence","\n")
  cat("\n")
  cat("+++++++++++++++","\n")
}
