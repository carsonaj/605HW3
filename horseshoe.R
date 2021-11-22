library(horseshoe)

HorseshoeMCMC = function(X,Y){
    startTime = Sys.time()
    InferenceResultList = horseshoe(Y,X)
    totalTime = Sys.time() - startTime
    InferenceResultList$CPUTime = totalTime
    return(InferenceResultList)
}

PostHorseshoeMCMC = function(InferenceResultList){
    beta_hats = abs(InferenceResultList$BetaHat)
    kmeans = kmeans(beta_hats, centers = c(min(beta_hats),max(beta_hats)))
    signal_id = c(1,2)[(kmeans$centers == max(kmeans$centers))]
    return(1 * (kmeans$cluster == signal_id))
}
