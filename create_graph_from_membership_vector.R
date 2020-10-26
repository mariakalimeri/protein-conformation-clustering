create_graph_from_membership_vector <- function(membershipVector) {
  # library(igraph)  # Needed library

  # First count the frequencies of each cluster in order to weight the nodes
  myHist <- hist(
    membershipVector, 
    seq(0.5, (max(membershipVector) + 0.5), 1), 
    plot = FALSE
  )
  nodeWeights <- myHist$counts

  # Create graph - NEW VERSION. It seems that the nodes are now 1-based, in the 
  # update??
  edgeList <- cbind(
    membershipVector[1:(length(membershipVector) - 1)], 
    membershipVector[2:(length(membershipVector))]
  )
  myG <- graph.edgelist(edgeList)

  # Simplify the graph by removing the loops
  waitingTimes <- is.loop(myG)
  myG <- simplify(myG, remove.multiple = FALSE, remove.loops = TRUE)

  # Count multiple edges and set them as an attribute to the newly created graph
  myG <- set.edge.attribute(myG, "weight", value = count.multiple(myG))
  myG <- simplify(myG, remove.multiple = TRUE)

  # This code is very sensitive to R and igraph updates. I added the following 
  # part after realizing that the code has suddenly stopped working correctly 
  # and after the simplification of the myG above it assings to the simplyfied 
  # edges the original weigth n times n = n^2 (new weight)
  myG <- set.edge.attribute(
    myG, 
    "weight", 
    value = sqrt(get.edge.attribute(myG, "weight"))
  )

  # Add nodeWeights as an attribute to the newly created graph
  V(myG)$weight <- nodeWeights

  # Create a version of the graph in the form of a weighted edgelist as MCL 
  # input
  weightedEdgeList <- get.edgelist(myG)
  weightedEdgeList <- cbind(
    weightedEdgeList, 
    get.edge.attribute(myG, "weight")
  )

  myG <- set.vertex.attribute(
    myG, 
    "label", 
    value = seq(1:max(membershipVector))
  )

  # write.graph(myG, file=paste("fnc", ".graphml", sep=""), format=c("graphml"))
  return(list(myG = myG, weightedEdgeList = weightedEdgeList))
}
