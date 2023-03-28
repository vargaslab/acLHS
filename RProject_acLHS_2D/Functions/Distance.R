Distance<-function (CoorX, CoorY) 
{
  geodata <- as.matrix(cbind(CoorX, CoorY))
  u <- as.vector(dist(geodata))
  DistMin <- min(u)
  DistMax <- max(u)
  Dist_Max_Min <- cbind(DistMax, DistMin)
  colnames(Dist_Max_Min) <- c("DMaxima", "DMinima")
  return(Dist_Max_Min)
}