convert_to_kinetic_membership_vector <- function(initialMem){
    # This script reads the original membershipVector and the file created by mcl
    # algorithm and converts the initial trajectory to the CG trajectory.
    # The mcl output has been save from before with a standard title.

    # Read MCL output file with standard title
    mclOutput <- strsplit(readLines("mclOut.txt"), "\n")
    
    # Convert clusterfile to numeric vector
    mclOutputNew <- vector('list', length(mclOutput))
    for (i in 1:length(mclOutput)){
        mclOutputNew[[i]] <- as.numeric((strsplit(mclOutput[[i]], "\t"))[[1]])
    }

    newMem <- rep(0,length(initialMem))
    for (i in 1:length(initialMem)){
        for (j in 1:length(mclOutputNew)){
	        if (sum(mclOutputNew[[j]]==initialMem[i]) == 1){
          	    newMem[i] <- j
            }
        }
    }
    
    # Finally re-enumerate starting from 1
    numberOfClusters <- length(unique(newMem))
    superNewMem <- rep(0, length(initialMem))
    for (iko in 1:numberOfClusters){
        superNewMem[which(newMem==unique(newMem)[iko])] <- iko
    }
    return(superNewMem)
}
