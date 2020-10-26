cluster_function_gui <- function(trjfile,
                               cutoff,
                               fittingInds,
                               outputProtInds,
                               verbose = TRUE,
                               cell = FALSE) {

  # ============== PRE-EXISTING CODE by B. Grant==============================
  dcd.frame <- function(trj, head, cell) {
    if (head$charmm && head$extrablock) {
      a <- readBin(trj, "integer", 1, endian = head$end)
      u <- readBin(trj, "numeric",
        size = 8, n = (a / 8),
        endian = head$end
      )
      a <- readBin(trj, "integer", 1, endian = head$end)
    }
    if (head$nfixed == 0) {
      a <- readBin(trj, "integer", 1, endian = head$end)
      x <- readBin(trj, "numeric",
        size = 4, n = (a / 4),
        endian = head$end
      )
      a <- readBin(trj, "integer", 1, endian = head$end)
      a <- readBin(trj, "integer", 1, endian = head$end)
      y <- readBin(trj, "numeric",
        size = 4, n = (a / 4),
        endian = head$end
      )
      a <- readBin(trj, "integer", 1, endian = head$end)
      a <- readBin(trj, "integer", 1, endian = head$end)
      z <- readBin(trj, "numeric",
        size = 4, n = (a / 4),
        endian = head$end
      )
      a <- readBin(trj, "integer", 1, endian = head$end)
    }
    else {
    }
    if (head$charmm && head$four.dims) {
      a <- readBin(trj, "integer", 1, endian = head$end)
      seek(trj, where = a, origin = "current")
      a <- readBin(trj, "integer", 1, endian = head$end)
    }
    if (cell) {
      coords <- c(u[c(1, 3, 6)], (180 / pi) * acos(u[c(5, 4, 2)]))
    } else {
      coords <- as.vector(rbind(x, y, z))
    }
    class(coords) <- "xyz"
    return(coords)
  }

  if (!file.exists(trjfile)) {
    stop(paste("No input DCD file found with name:", trjfile))
  }
  trj <- file(trjfile, "rb")
  head <- dcd.header(trj, verbose)
  nframes <- head$nframe
  natoms <- head$natom

  if (verbose) {
    cat("Reading (x100)")
  }
  store <- NULL
  # =============== END OF PRE-EXISTING FUNCTIONS ============================

  # Read the first frame (outside the for loop) due to the linearity of leader
  curr.pos <- seek(trj, where = 0, origin = "current")
  # Store only current and previous frame for fitting
  #    coords <- matrix(NA, nrow = 2, ncol = natoms *  3)
  coords <- dcd.frame(trj, head, cell)
  store <- cbind(store, curr.pos)
  # Leaders start with the first frame
  leaders <- coords[outputProtInds]
  coordsFirstFrame <- leaders
  clusterIndex <- 1
  membershipVector <- clusterIndex
  leadersFramePositions <- 1
  if (length(leaders) == 0) {
    stop("No atoms are selected. Check that you chose the correct chain identifier. If your PDB file has no chain ID record for the protein you want to cluster, you must create one (e.g. with VMD) \n If all the above are correct, you may be helpless. Report this as a bug.")
  }
  xyzLeaders <- coordsFirstFrame

  for (i in 2:nframes) {
    curr.pos <- seek(trj, where = 0, origin = "current")
    if (curr.pos <= head$end.pos) {
      coords <- dcd.frame(trj, head, cell)
      coords <- coords[outputProtInds]

      # ========= MARIA's LEADER CLUSTERING ==============================
      ## FITTING ##
      fitted.coords <- fit.xyz(
        fixed = coordsFirstFrame, mobile = coords,
        fixed.inds = fittingInds,
        mobile.inds = fittingInds
      )

      ## RMSD Calculation of frame from all existing leaders
      rmsdValues <- rmsd(
        a = fitted.coords, b = leaders, a.inds = fittingInds,
        b.inds = fittingInds, fit = FALSE
      )
      if (min(rmsdValues) <= cutoff) {
        membershipVector[i] <- which(rmsdValues == min(rmsdValues))[1]
      } else {
        leaders <- rbind(leaders, fitted.coords)
        clusterIndex <- clusterIndex + 1
        membershipVector <- c(membershipVector, clusterIndex)
        leadersFramePositions <- c(leadersFramePositions, i)
        xyzLeaders <- rbind(xyzLeaders, fitted.coords)
      }
      # ========= END OF MARIA's LEADER CLUSTERING =======================

      if (verbose) {
        incProgress(1 / nframes, detail = paste("Reading frame", i))
      }
      store <- cbind(store, curr.pos)
    }
    else {
      print("Premature end of file")
      print(paste("  last frame:", i, "nframe:", head$nframe))
      break
    }
  }
  if (verbose) {
    cat("\n")
  }
  close(trj)

  # Get the upper limit of membership
  vectorCluster <- rep(1, nframes)
  for (fr in 2:nframes) {
    clus <- membershipVector[fr]
    if (clus > max(vectorCluster)) {
      vectorCluster[fr] <- clus
    } else {
      vectorCluster[fr] <- vectorCluster[fr - 1]
    }
  }
  return(list(membershipVector = membershipVector, leadersFramePositions = leadersFramePositions, vectorCluster = vectorCluster, nframes = nframes, xyzLeaders = xyzLeaders))
}

dcd.header <- function(trj, verbose = TRUE, ...) {
  end <- .Platform$endian
  check <- readBin(trj, "integer", 1, endian = end)
  if (check != 84) {
    if (end == "little") {
      end <- "big"
    }
    else {
      end <- "little"
    }
    check <- readBin(writeBin(check, raw()), "integer",
      1,
      endian = end
    )
    if (check != 84) {
      close(trj)
      stop("PROBLEM with endian detection")
    }
  }
  hdr <- readChar(trj, nchars = 4)
  cur.pos <- seek(trj, where = 1, origin = "end")
  end.pos <- seek(trj, where = cur.pos, origin = "start")
  icntrl <- readBin(trj, "integer", 20, endian = end)
  nframe <- icntrl[1]
  first <- icntrl[2]
  step <- icntrl[3]
  nstep <- icntrl[4]
  nfile <- nstep / step
  last <- first + (step * nframe)
  ndegf <- icntrl[8]
  nfixed <- icntrl[9]
  delta <- icntrl[10]
  cryst <- icntrl[11]
  block <- icntrl[12]
  vers <- icntrl[20]
  a <- readBin(trj, "integer", 1, endian = end)
  rm(icntrl)
  charmm <- FALSE
  extrablock <- FALSE
  four.dims <- FALSE
  if (vers != 0) {
    charmm <- TRUE
    if (cryst == 1) {
      extrablock <- TRUE
    }
    if (block == 1) {
      four.dims <- TRUE
    }
  }
  else {
    cur.pos <- seek(trj, where = 44, origin = "start")
    delta <- readBin(trj, "double", 1, endian = end)
    seek(trj, where = cur.pos, origin = "start")
  }
  a <- readBin(trj, "integer", 1, endian = end)
  ntitle <- readBin(trj, "integer", 1, endian = end)
  title <- NULL
  cur.pos <- seek(trj, where = NA)
  for (i in 1:ntitle) {
    ll <- try(title <- c(title, readChar(trj, 80)), silent = TRUE)
  }
  if (class(ll) == "try-error") {
    cur.pos <- seek(trj,
      where = (80 * ntitle + cur.pos),
      origin = "start"
    )
  }
  a <- readBin(trj, "integer", 1, endian = end)
  a <- readBin(trj, "integer", 1, endian = end)
  natom <- readBin(trj, "integer", 1, endian = end)
  a <- readBin(trj, "integer", 1, endian = end)
  if (nfixed != 0) {
    a <- readBin(trj, "integer", 1, endian = end)
    free.ind <- readBin(trj, "integer", (natom - nfixed),
      endian = end
    )
    a <- readBin(trj, "integer", 1, endian = end)
    print("FIXED ATOMS IN SIMULATION => CAN'T READ YET")
  }
  if (verbose) {
    cat(" NATOM =", natom, "\n")
    cat(" NFRAME=", nframe, "\n")
    # cat(" ISTART=", first, "\n")
    # cat(" last  =", last, "\n")
    # cat(" nstep =", nstep, "\n")
    # cat(" nfile =", nfile, "\n")
    # cat(" NSAVE =", step, "\n")
    # cat(" NDEGF =", ndegf, "\n")
    # cat(" version", vers, "\n")
  }
  header <- list(
    natom = natom, nframe = nframe, first = first,
    last = last, nstep = nstep, nfile = nfile, step = step,
    ndegf = ndegf, nfixed = nfixed, charmm = charmm,
    extrablock = extrablock, four.dims = four.dims, end.pos = end.pos,
    end = end
  )
}
