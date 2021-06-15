#' German electricity prices and demand 
#'
#'This data set consist of electricity spot prices traded at the European Energy Exchange (EEX) in Leipzig and electricity demand reported by European Network of Transmission System Operators for Electricity from January 1st 2006 to  September 30th 2008,
#'excluding weekends and holidays. It was processed from its original version provided in the on-line supplementary materials of Liebl (2013).  
#'
#
#'
#' @format A list with 3 fields:
#' \describe{
#'   \item{price}{hourly evaluated  electricity spot prices in  638 days,  represented as a matrix of dimension 638 by 24.}
#'   \item{demand}{daily average of electricity demand in 638 days, represented as a vector.}
#'   \item{date}{dates represented as a Date object with each day displayed in the form of %d%b%y.}
#'   ...
#' }
#' @source Liebl, Dominik. "Modeling and forecasting electricity spot prices: A functional data perspective." The Annals of Applied Statistics 7.3 (2013): 1562-1592.
"GED"