
SeqFDR <- function(m,FDR,a=2,R)
{
  #Beta[m]=(6*FDR/pi^2)/m^2     #if a=2
  #pre-assigned level (FDR) that ensures FDR and mFDR remains below
  #Beta[m]=(FDR/1.202)/m^3    #if a=3
  #Beta[m]=(90*FDR/pi^4)/m^4   #if a=4
  if (m==1){
    Alpha <- (6*FDR/pi^2)/m^2   #for a=2
  }
  if(m>1){
    Alpha <- (6*FDR/pi^2)/m^2 *(sum(R[1:(m-1)])+1)
  } 
  return(Alpha)
}
