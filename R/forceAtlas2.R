forceAtlas2 <- 
function(mat, iterations = 100, linlog = FALSE, pos = NULL, nohubs = FALSE,
                              k = 400, gravity = 1, ks = 0.1, ksmax = 10, delta = 1, center = NULL,
                              tolerance = 0.1, dim = 2, ovs = NULL, ovd = NULL, plotstep = 0, seed = NULL, opts=NULL) {
  
  if (!is.matrix(mat)) stop("input should be a matrix")
  if (!isSymmetric(mat)) stop("input matrix is asymmetric")
  if (!is.null(ovs) & !is.null(ovd)) stop("only ovs or ovd parameter is allowed")
  if (!is.null(opts) & !is.list(opts)) stop("opts should be provided in a list format")

  A = unname(mat)

  #### center of gravity is by default set to the origin
  if(is.null(center)) center <- rep(0,dim)

  nnodes <- nrow(A)

  #### Binary will be an incidence matrix (0-not connected, 1-connected)
  Binary <- A
  Binary[Binary!=0] <- 1
  
  #### Deg will be a vector containing a degree of each node
  Deg <- rowSums(Binary)

  #### If there are no initial coordinates of points,
  #### they are chosen at random from 1000^dim space
  if (is.null(pos)) {
    difference <- 2000/((nnodes*dim)-1)
    if (!is.null(seed)) set.seed(seed)
    position <- matrix(sample(seq(-1000,1000,difference),nnodes*dim), nrow = nnodes, ncol = dim)
  } else {
    position <- pos
  }

  #### None of the nodes should be exactly at the center of gravity
  temp <- apply(position, 1, identical, center)
  if (any(temp)) position[temp,] <- center + 0.01
  rm(temp)

  #### Forces1 will be a matrix containing all the sums of forces acting on points at a step
  Forces1 <- matrix(0, nrow = dim, ncol = nnodes)

  #### displacement will be a matrix of points' movement at the current iteration
  # displacement <- matrix(0, nrow = dim, ncol = nnodes)

  # m <- nrow(position)
  # if (!is.null(ovd)) ov = seq(0, ovd, ovd/(iterations-1)) # linear growth
  # if (!is.null(ovd)) ov = exp(seq(0, log(ovd), log(ovd)/(iterations-1))) # exponential size growth
  if (!is.null(ovd)) ov = ovd * pbeta(seq(0,1, 1/(iterations-1)), 50,10) # beta distribution shape1 5 time larger than shape2
  if (!is.null(ovs)) ov = rep(ovs, iterations)

  for (iteration in 1:iterations) {

    if (!is.null(opts$progress)) sapply(iteration, opts$progress)

    # displacement <- displacement * 0
    #### Forces2 is the matrix of the forces from previous step
    #### Forces1 is the matrix of the forces from current step
    Forces2 <- Forces1
    # Forces1 <- matrix(, nrow = dim, ncol = 0)

    #### Calculate the Forces for each node
    ### Distance matrix between all nodes
    distances <- as.matrix(dist(position))
    distances[which(distances < 0.01)] <- 0.01 # we impose a minimum distance

    ### Each element of the list contains a matrix with the j = 1,2,..., dim dimension of the unitary vector 1
    mylist <- vector("list", dim)
    for (j in 1:dim){
      mylist[[j]] <- (tcrossprod(position[,j], rep(1,nnodes)) - tcrossprod(rep(1,nnodes), position[,j])) / distances
    }

    ### Prevent overlapping: if overlapping constant (ov) is provided; calculate d'(n1,n2); if d'(n1,n2)>0 use d' else use d to compute Fa and Fr
    # development note: the ov constant corresponds to the node size (arbitrary in R). We assume the node size is the same.
    if (!is.null(ovs) | !is.null(ovd)) {
      d <- distances - (2 * ov[iteration])
      diag(d) = 0
      # distances[which(d > 0)] <- d[which(d > 0)]
      pos = which(d > 0)
      neg = which(d < 0)

      Fa = d
      Fr = k * ((tcrossprod(rep(1,nnodes), Deg) + 1) * (tcrossprod(Deg, rep(1,nnodes)) + 1)) / d

      if (length(neg) > 0) {
        Fa[neg] = 0
        Fr[neg] = 100 * ((tcrossprod(rep(1,nnodes), Deg) + 1) * (tcrossprod(Deg, rep(1,nnodes)) + 1))[neg]
      }
    }
    else {
      ### Calculate the repulsion Force
      Fr <- k * ((tcrossprod(rep(1,nnodes), Deg) + 1) * (tcrossprod(Deg, rep(1,nnodes)) + 1)) / distances

      # The classical attraction force is just based on distance
      Fa <- distances
    }

    # The linlog mode calculates the attraction force as log(1+d(n1,n2))
    if (linlog) Fa <- log(1 + Fa)
    
    # Edge weights. The edges are weighted based on parameter delta. delta=0 implies no weight
    # development note: if delta = 0, then 0^0 = 1; all entries in the incidence matrix will == 1; multiply by Binary instead
    if (delta == 0) Fa <- Binary * Fa
    Fa <- (A^delta) * Fa 

    # Dissuade Hubs. This mode is meant to grant authorities (nodes with high indegree) a more central position than hubs (nodes with high outdegree)
    if (nohubs) Fa <- Fa/(tcrossprod(Deg, rep(1, nnodes)) + 1)

    ### Function to calculate the Attraction and Repulsion forces
    Farfunction <- function(x) rowSums(x*(Fr-Fa), na.rm=TRUE)
    ### And we aggregate it over all dimensions
    Far <- do.call(rbind, lapply(mylist, Farfunction))
    
    ### Unitary Vector 2, the directions between each point and the center
    uv2 <- apply(matrix(rep(center,nnodes),nrow=nnodes,byrow=T)-position, 1, function(x) x/sqrt(sum(x^2)))
    
    ### The gravity force
    # Fg <- uv2 * matrix(rep(gravity*(rowSums(A)+1), dim), nrow=dim, byrow=T)
    Fg <- uv2 * matrix(rep(gravity*(Deg+1), dim), nrow=dim, byrow=T)
    
    ### Forces 1 is the sum between all forces: Far (Fa + Fr) and Fg
    Forces1 <- Far + Fg
    Forces1 <- round(Forces1, 2) # use the first two decimals for the Forces.

    #### Swing is the vector of the swingings of all points
    swing <- abs(colSums((Forces1-Forces2)^2)^(1/2))
    Global_swing <- sum((Deg+1)*swing)

    # If the swing of all nodes is zero, then convergence is reached and we break.
    if (all(swing==0)) {
      if (!plotstep==0 & dim==2) {
        plot(position, main=paste0("iteration: ",iteration), xlab="", ylab="")
      }
      print(paste0("Convergence reached at step ", iteration))
      break
    }

    #### tra is the vector of the traction of all points
    tra <- abs(colSums((Forces1+Forces2)^2)^(1/2))/2
    Global_tra <- sum((Deg+1)*tra)

    #### Global speed calculation
    Global_speed <- tolerance * Global_tra/Global_swing
    #### speed is the vector of individual speeds of points
    speed <- ks * Global_speed / (1 + Global_speed * (swing)^(1/2))

    #### Imposing constrains on speed
    speed_constrain <- ksmax/abs(colSums((Forces1^2))^(1/2))
    speed <- ifelse(speed >= speed_constrain, speed_constrain, speed)

    #### calculating displacement and final position of points after iteration
    displacement <- Forces1 * t(matrix(rep(speed,dim), nrow=nnodes, ncol=dim))
    position <- position + t(displacement)

    #### Iteration plot. This is simply to see the evolution of the positions over iterations
    if (!plotstep==0 & dim==2) {
      if (iteration%%plotstep == 0) {
        plot(position, main=paste0("iteration: ",iteration), xlab="", ylab="")
      }
    }
  }

  return(position)
}
