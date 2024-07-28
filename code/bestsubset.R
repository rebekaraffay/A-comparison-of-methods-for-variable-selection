check.xy = function(x, y) {
  if (is.null(x) || !is.numeric(x)) stop("x must be a numeric matrix")
  if (is.null(y) || !is.numeric(y)) stop("y must be a numeric vector")
  if (length(y) == 0) stop("Must have length(y) > 0")
  if (nrow(x) != length(y)) stop("nrow(x) and length(y) must match")
  if (ncol(x) == 0) stop("Must have ncol(x) > 0")
  if (check.cols(x)) stop("x cannot have duplicate columns")
}

# Make sure that no two columms of A are the same (this works with
# probability one)

check.cols = function(A) {
  b = rnorm(nrow(A))
  a = sort(t(A)%*%b)
  return(any(diff(a)==0))
}

check.bool = function(b) {
  if (is.null(b) || length(b)!=1 || !is.logical(b))
    stop(paste(deparse(substitute(b)),"must be a Boolean"))
}

check.num = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a))
    stop(paste(deparse(substitute(a)),"must be a number"))
}

check.int = function(i) {
  if (is.null(i) || length(i)!= 1 || !is.numeric(i) || round(i) != i)
    stop(paste(deparse(substitute(i)),"must be an integer"))
}

check.pos.num = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a) || a<0)
    stop(paste(deparse(substitute(a)),"must be a positive number"))
}

check.pos.int = function(i) {
  if (is.null(i) || length(i)!= 1 || !is.numeric(i) || round(i) != i || i<1)
    stop(paste(deparse(substitute(i)),"must be a positive integer"))
}

check.num.01 = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a) || a<0 || a>1)
    stop(paste(deparse(substitute(a)),"must be a number between 0 and 1"))
}


# Special linear time order function, works only when x is
# a vector of integers

Order = function(x) {
  n = length(x)
  o = numeric(n)
  o[x] = Seq(1,n)
  return(o)
}

# Returns a sequence of integers from a to b if a <= b,
# otherwise nothing. You have no idea how important this
# function is...

Seq = function(a, b, ...) {
  if (a<=b) return(seq(a,b,...))
  else return(numeric(0))
}

# Returns the sign of x, with Sign(0) = 1

Sign = function(x) {
  return(-1+2*(x>=0))
}

# Truncate function
trunc = function(x, a, b) {
  return(ifelse(x>a,ifelse(x<b,2*x-a-b,b-a),a-b))
}

# Match duplicate row function

match.row = function(mat, a) {
  return(which(rowSums(abs(scale(mat,center=a,scale=F)))==0))
}

# Trace convenience function

Trace = function(mat) sum(diag(mat))

# Centering and scaling convenience function

standardize = function(x, y, intercept, normalize) {
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  
  if (intercept) {
    bx = colMeans(x)
    by = mean(y)
    x = scale(x,bx,FALSE)
    y = y-mean(y)
  } else {
    bx = rep(0,p)
    by = 0
  }
  if (normalize) {
    sx = sqrt(colSums(x^2))
    x = scale(x,FALSE,sx)
  } else {
    sx = rep(1,p)
  }

  return(list(x=x,y=y,bx=bx,by=by,sx=sx))
}

# Enlist function

enlist <- function (...) 
{
    result <- list(...)
    if ((nargs() == 1) & is.character(n <- result[[1]])) {
        result <- as.list(seq(n))
        names(result) <- n
        for (i in n) result[[i]] <- get(i)
    }
    else {
        n <- sys.call()
        n <- as.character(n)[-1]
        if (!is.null(n2 <- names(result))) {
            which <- n2 != ""
            n[which] <- n2[which]
        }
        names(result) <- n
    }
    result
}

#' Plot the results over several simulation settings.
#'
#' Plot results over several sets of simulations saved in files, where the same
#'   methods are run over different simulations settings.
#'
#' @param file.list Vector of strings that point to saved sim objects (each
#'   object produced by a call to \code{\link{sim.master}}).
#' @param row,col One of "beta", "rho", or "snr", indicating the variables for
#'   the rows and columns of plotting grid; note that row and col must be
#'   different.  Default is row="beta" and col="rho", so that the plotting grid
#'   displays the metric specified in what (default is relative test error, see
#'   below) versus the SNR, in a plotting grid with the coefficient types across
#'   the rows, and the correlation levels across the columns.
#' @param method.nums The indices of the methods that should be plotted. Default
#'   is NULL, in which case all methods are plotted.
#' @param method.names The names of the methods that should be plotted. Default
#'   is NULL, in which case the names are extracted from the sim objects.
#' @param what One of "error", "risk", "prop", or "nonzero", indicating whether
#'   to plot the relative test error, test proportion of variance explained, or
#'   number of nonzeros, for each method across the given SNR levels. When what
#'   is "prop", the x-axis is population proportion of variance explained,
#'   instead of SNR. Default is "error".
#' @param rel.to An index for a base method: when the what argument is "error"
#'   or "risk", the rel.to argument specifies the method with respect to which
#'   the relative error or relative risk is calculated. The default value of
#'   NULL means different things when what is equal to "error" and "risk". In
#'   the former case, the test error is calculated relative to the oracle; in
#'   the latter, the risk is calculated relative to the null model.
#' @param tuning Either "validation" or "oracle", indicating whether the tuning
#'   parameter for each method should be chosen according to minimizing
#'   validation error, or according to minimizing test error. Default is
#'   "validation".
#' @param type Either "ave" or "med", indicating whether the average or median
#'   of the test error (or number of nonzeros, if what is "nonzero") should be
#'   displayed. Default is "ave".
#' @param std Should standard errors be displayed (in parantheses)? When type
#'   is set to "med", the median absolute deviations are shown in place of the
#'   standard errors. Default is TRUE.
#' @param lwd,pch,main,ylim,legend.pos graphical parameters.
#' @param make.pdf Should a pdf be produced? Default is FALSE.
#' @param fig.dir,file.name The figure directory and file name to use, only
#'   when make.pdf is TRUE. Defaults are "." and "sim". (An extension of "pdf"
#'   is always appended to the given file name.)
#' @param w,h the width and height (in inches) for the plot, used only when
#'   make.pdf is TRUE. Defaults are 8 and 10, respectively.
#'
#' @export plot.from.file

plot.from.file = function(file.list,
                          row=c("beta","rho","snr"), col=c("rho","beta","snr"),
                          method.nums=NULL, method.names=NULL,
                          what=c("error","risk","prop","F","nonzero"), rel.to=NULL,
                          tuning=c("validation","oracle"), type=c("ave","med"),
                          std=TRUE, lwd=1, pch=19, main=NULL, ylim=NULL,
                          legend.pos=c("bottom","right","top","left","none"),
                          make.pdf=FALSE, fig.dir=".", file.name="sim",
                          w=8, h=10) {

  # Check for ggplot2 package
  if (!require("ggplot2",quietly=TRUE)) {
    stop("Package ggplot2 not installed (required here)!")
  }

  row = match.arg(row)
  col = match.arg(col)
  if (row==col) stop("row and col must be different")

  what = match.arg(what)
  tuning = match.arg(tuning)
  type = match.arg(type)
  legend.pos = match.arg(legend.pos)

  # Set the method numbers and names
  sim.obj = readRDS(file.list[1])
  if (is.null(method.nums)) method.nums = 1:length(sim.obj$err.test)
  if (is.null(method.names)) method.names =
                               names(sim.obj$err.test[method.nums])
  N = length(method.nums)

  # Set the base number and name
  if (is.null(rel.to)) {
    base.num = 0
    base.name = ifelse(what=="error","Bayes","null model")
  }
  else {
    base.num = which(method.nums==rel.to)
    base.name = tolower(method.names[base.num])
  }

  # Set the y-label
  ylab = switch(what,
                error=paste0("Relative test error (to ",base.name,")"),
                risk=paste0("Relative risk (to ",base.name,")"),
                prop="Proportion of variance explained",
                F="F classification of nonzeros",
                nonzero="Number of nonzeros")

  # Collect the y-variable from the file list
  yvec = ybar = beta.vec = rho.vec = snr.vec = c()
  for (i in 1:length(file.list)) {
    sim.obj = readRDS(file.list[i])
    beta.vec = c(beta.vec,rep(sim.obj$beta.type,N))
    rho.vec = c(rho.vec,rep(sim.obj$rho,N))
    snr.vec = c(snr.vec,rep(sim.obj$snr,N))

    z = sim.obj[[switch(what,
                        error="err.test",
                        risk="risk",
                        prop="prop",
                        F="F1",
                        nonzero="nzs")]]
    res = tune.and.aggregate(sim.obj, z)

    # For prop, F  and nonzero we ignore any request for a relative metric
    if (what=="prop" || what=="F" || what=="nonzero") {
      yvec = c(yvec,res[[paste0("z.",substr(tuning,1,3),".",type)]][method.nums])
      ybar = c(ybar,res[[paste0("z.",substr(tuning,1,3),".",
                                ifelse(type=="ave","std","mad"))]][method.nums])
    }

    # For err and risk we respect the request for a relative metric
    else {
      # First build the relative metric
      met = res[[paste0("z.",substr(tuning,1,3))]]#[method.nums]
      if (base.num == 0 && what=="error") denom = sim.obj$sigma^2
      else if (base.num == 0 && what=="risk") denom = sim.obj$risk.null
      else denom = met[[base.num]]
      z.rel = lapply(met, function(v) v / denom)
      # Now aggregate the relative metric
      res2 = tune.and.aggregate(sim.obj, z.rel, tune=FALSE)
      yvec = c(yvec,unlist(res2[[paste0("z.",type)]])[method.nums])
      ybar = c(ybar,unlist(res2[[paste0("z.",ifelse(type=="ave",
                                                        "std","mad"))]])[method.nums])
    }
  }
  # Set the x-variable and x-label
  xvec = snr.vec
  xlab = "Signal-to-noise ratio"

  # Set the y-limits
  if (is.null(ylim)) ylim = range(yvec-ybar, yvec+ybar)
  # Produce the plot
  beta.vec = factor(beta.vec)
  rho.vec = factor(rho.vec)
  snr.vec = factor(snr.vec)
  levels(beta.vec) = paste("Beta-type", levels(beta.vec))
  levels(rho.vec) = paste("Correlation", levels(rho.vec))

  dat = data.frame(x=xvec, y=yvec, se=ybar,
                   beta=beta.vec, rho=rho.vec, snr=snr.vec,
                   Method=factor(rep(method.names, length=length(xvec))))

  gp = ggplot(dat, aes(x=x,y=y,color=Method)) +
    xlab(xlab) + ylab(ylab) + coord_cartesian(ylim=ylim) +
    geom_line(lwd=lwd) + geom_point(pch=pch) +
    facet_grid(formula(paste(row,"~",col))) +
    theme_bw() + theme(legend.pos=legend.pos)
  if (!("snr" %in% c(row,col))) {
    # If SNR is being plotted on the x-axis in each plot, then define special
    # x-axis ticks and put the x-axis on a log scale
    snr.breaks = round(exp(seq(from=min(log(xvec)),
                               to=max(log(xvec)),length=4)),2)
    gp = gp + scale_x_continuous(trans="log", breaks=snr.breaks)
  }
  if (std) gp = gp + geom_errorbar(aes(ymin=y-se,ymax=y+se), width=0.02)
  if (what=="error") gp = gp + geom_line(aes(x=x, y=1+x), lwd=0.5,
                                         linetype=3, color="black")
  if (what=="prop") gp = gp + geom_line(aes(x=x, y=x/(1+x)), lwd=0.5,
                                        linetype=3, color="black")
  if (what =="nonzero") gp = gp + geom_line(aes(x=x, y=sim.obj$s), lwd=0.5,
                                            linetype=3, color="black")
  if (!is.null(main)) gp = gp + ggtitle(main)
  if (!is.null(ylim)) gp = gp + coord_cartesian(ylim=ylim)
  if (make.pdf) ggsave(sprintf("%s/%s.pdf",fig.dir,file.name),
                       height=h, width=w, device="pdf")
  else gp
}

get.stem = function(str.vec) {
  str.list = strsplit(str.vec, split="\\.")
  k = 0
  while (TRUE) {
    vec = c()
    for (i in 1:length(str.list)) {
      if (length(str.list[[i]]) < k+1) break
      vec = c(vec, str.list[[i]][k+1])
    }
    if (length(vec) < length(str.list)) break
    if (length(unique(vec)) > 1) break
    k = k+1
  }
  if (k == 0) stem = "foo"
  else stem = paste0(str.list[[1]][1:k], collapse=".")
  return(stem)
}


##############################


#' Tuning and aggregation function for sim object.
#'
#' Tune and aggregate a given metric across a set of simulations, stored in a
#'   sim object (produced by \code{\link{sim.master}}).
#'
#' @param obj The sim object.
#' @param z The metric to be tuned. Should be a list of lengh N, where N is the
#'   number of methods under consideration in the simulation, and in the ith
#'   element of the list should be a matrix of dimension nrep x m, where nrep is
#'   the number of repetitions in the simulation, and m is the number of tuning
#'   parameters for the ith method.
#' @param tune Either TRUE or FALSE, indicating whether the metric should be
#'   tuned; if FALSE, then it is only aggregated (for each value of the tuning
#'   parameter). Default is TRUE.
#'
#' @return A list with components z.ave, z.std, z.med, z.mad, z.val, z.ora,
#'   z.val.ave, z.val.std, z.val.med, z.val.mad, z.ora.ave, z.ora.std,
#'   z.ora.med, z.ora.mad.  The elements z.ave, z.std, z.med, z.mad are lists of
#'   length N, where N is the number of methods in the simulation; in the ith
#'   element of z.ave is a vector of length m, the number of tuning parameters
#'   for the ith method, giving the averages of the metric across the
#'   repetitions for each tuning parameter value; similarly z.std contains
#'   standard errors, z.med contains medians, and z.mad contains median absolute
#'   deviations.  The elements z.val, z.ora are matrices of dimension nrep x N,
#'   where nrep is the number of repetitions, whose rows contain the metric
#'   after it has been tuned on each repetition by selecting the tuning
#'   parameter for each method in one of two ways: validation tuning, where the
#'   tuning parameter is selected to minimize prediction error on the validation
#'   set, and oracle tuning, where the tuning parameter is selected to minimize
#'   average test error. The elements named z.val.xxx, where xxx is one of ave,
#'   std, med, or mad, are each vectors of length N, containing the average of
#'   the metric for each method under validation tuning; similarly for the
#'   postfixes std, med, and mad; and for the elements named z.ora.xxx.
#'
#' @export tune.and.aggregate

tune.and.aggregate = function(obj, z, tune=TRUE) {
  N = length(obj$err.test) # Number of methods
  nrep = nrow(obj$err.test[[1]]) # Number of repetitions

  # Aggregate
  z.ave = z.std = z.med = z.mad = vector(mode="list",length=N)
  for (j in 1:N) {
    z[[j]] = matrix(z[[j]], nrow=nrep) # Just in case it is not a matrix
    z.ave[[j]] = colMeans(z[[j]], na.rm=TRUE)
    z.std[[j]] = apply(z[[j]], 2, sd, na.rm=TRUE) /
      sqrt(colSums(!is.na(z[[j]])))
    z.med[[j]] = apply(z[[j]], 2, median, na.rm=TRUE)
    z.mad[[j]] = apply(z[[j]], 2, mad, na.rm=TRUE) /
      sqrt(colSums(!is.na(z[[j]])))
  }
  out = enlist(z.ave,z.std,z.med,z.mad)

  # Tune and aggregate
  if (tune) {
    z.val = z.ora = vector(mode="list",length=N)
    z.val.ave = z.val.std = z.val.med = z.val.mad = rep(NA,N)
    z.ora.ave = z.ora.std = z.ora.med = z.ora.mad = rep(NA,N)

    for (j in 1:N) {
      # Validation tuning
      z.val[[j]] = z[[j]][1:nrep + (obj$tun.val[,j]-1)*nrep]
      z.val.ave[j] = mean(z.val[[j]], na.rm=TRUE)
      z.val.std[j] = sd(z.val[[j]], na.rm=TRUE) / sqrt(sum(!is.na(z.val[[j]])))
      z.val.med[j] = median(z.val[[j]], na.rm=TRUE)
      z.val.mad[j] = mad(z.val[[j]], na.rm=TRUE) / sqrt(sum(!is.na(z.val[[j]])))

      # Oracle tuning
      z.ora[[j]] = z[[j]][,obj$tun.ora[j]]
      z.ora.ave[j] = mean(z.ora[[j]], na.rm=TRUE)
      z.ora.std[j] = sd(z.ora[[j]], na.rm=TRUE) / sqrt(sum(!is.na(z.ora[[j]])))
      z.ora.med[j] = median(z.ora[[j]], na.rm=TRUE)
      z.ora.mad[j] = mad(z.ora[[j]], na.rm=TRUE) / sqrt(sum(!is.na(z.ora[[j]])))
    }

    out = c(out, enlist(z.val,z.ora,z.val.ave,z.val.std,z.val.med,z.val.mad,
                        z.ora.ave,z.ora.std,z.ora.med,z.ora.mad))
  }

  return(out)
}

##############################

#' Print function for sim object.
#'
#' Summarize and print the results of a set of simulations, stored in an object
#'   of class sim (produced by \code{\link{sim.master}}).
#'
#' @param x The sim object.
#' @param what Either "error" or "risk", indicating whether the relative test
#'     error (to the Bayes error) or relative risk (to the null risk) should be
#'     displayed. Default is "error".
#' @param type Either "ave" or "med", indicating whether the average or median
#'   of the relative test error metric should be displayed. Default is "ave".
#' @param std Should standard errors be displayed (in parantheses)? When type
#'   is set to "med", the median absolute deviations are shown in place of the
#'   standard errors. Default is TRUE.
#' @param digits Number of digits to display. Default is 3.
#' @param ... Other arguments (currently not used).
#'
#' @export print.sim
#' @export

print.sim = function(x, what=c("error","risk"), type=c("ave","med"), std=TRUE,
                     digits=3, ...) {
  what = match.arg(what)
  type = match.arg(type)

  if (!is.null(x$call)) {
    cat("\nCall:\n")
    dput(x$call)
  }

  # Construct relative test error (to Bayes) or relative risk (to null)
  N = length(x$err.test) # Number of methods
  err.rel = vector(mode="list",length=N)
  for (j in 1:N) {
    if (what=="error") err.rel[[j]] = x$err.test[[j]] / x$sigma^2
    else err.rel[[j]] = x$risk[[j]] / x$risk.null
  }
  err.obj = tune.and.aggregate(x, err.rel)
  nzs.obj = tune.and.aggregate(x, x$nzs)

  cat("\nResults for tuning parameters chosen based on validation set:\n\n")
  if (type=="ave") {
    col1 = err.obj$z.val.ave
    col2 = nzs.obj$z.val.ave
    col1.std = err.obj$z.val.std
    col2.std = err.obj$z.val.std
  }
  else {
    col1 = err.obj$z.val.med
    col2 = nzs.obj$z.val.med
    col1.std = err.obj$z.val.mad
    col2.std = err.obj$z.val.mad
  }

  tab = round(cbind(col1,col2),digits)
  tab.std = round(cbind(col1.std,col2.std), digits)
  if (std) tab = matrix(paste0(tab," (",tab.std,")"),ncol=2)
  rownames(tab) = names(x$err.test)
  colnames(tab) = c(paste("Relative",ifelse(what=="error","test error","risk")),
                    "Number of nonzeros")
  print(tab,quote=F)

  cat("\nResults for tuning parameters chosen based on test set (oracle):\n\n")
    if (type=="ave") {
    col1 = err.obj$z.ora.ave
    col2 = nzs.obj$z.ora.ave
    col1.std = err.obj$z.ora.std
    col2.std = err.obj$z.ora.std
  }
  else {
    col1 = err.obj$z.ora.med
    col2 = nzs.obj$z.ora.med
    col1.std = err.obj$z.ora.mad
    col2.std = err.obj$z.ora.mad
  }

  tab = round(cbind(col1,col2),digits)
  tab.std = round(cbind(col1.std,col2.std), digits)
  if (std) tab = matrix(paste0(tab," (",tab.std,")"),ncol=2)
  rownames(tab) = names(x$err.test)
  colnames(tab) = c(paste("Relative",ifelse(what=="error","test error","risk")),
                    "Number of nonzeros")
  print(tab,quote=F)

  cat("\n")
  invisible()
}

#' Plot function for sim object.
#'
#' Plot the results of a set of simulations, stored in an object of class sim
#'   (produced by \code{\link{sim.master}}).
#'
#' @param x The sim object.
#' @param method.nums the indices of the methods that should be plotted. Default
#'   is to 1:length(x$err.test), which plots all methods.
#' @param method.names the names of the methods that should be plotted. Default
#'   is NULL, in which case the names are extracted from the sim object.
#' @param what Either "error" or "risk", indicating whether the relative test
#'     error (to the Bayes error) or relative risk (to the null risk) should be
#'     displayed. Default is "error".
#' @param type Either "ave" or "med", indicating whether the average or median
#'   of the relative test error metric should be displayed. Default is "ave".
#' @param std Should standard errors be displayed (in parantheses)? When type
#'   is set to "med", the median absolute deviations are shown in place of the
#'   standard errors. Default is TRUE.
#' @param lwd,pch,main,legend graphical parameters.
#' @param make.pdf Should a pdf be produced? Default is FALSE.
#' @param fig.dir,file.name The figure directory and file name to use, only
#'   when make.pdf is TRUE. Defaults are "." and "sim". (An extension of "pdf"
#'   is always appended to the given file name.)
#' @param w,h the width and height (in inches) for the plot, used only when
#'   make.pdf is TRUE. Defaults are 6 for both.
#'
#' @export plot.sim
#' @export

plot.sim = function(x, method.nums=1:length(x$err.test), method.names=NULL,
                    what=c("error","risk"), type=c("ave","med"), std=TRUE,
                    lwd=1, pch=19, main=NULL, legend=TRUE, make.pdf=FALSE,
                    fig.dir=".", file.name="sim", w=6, h=6) {

  # Check for ggplot2 package
  if (!require("ggplot2",quietly=TRUE)) {
    stop("Package ggplot2 not installed (required here)!")
  }

  what = match.arg(what)
  type = match.arg(type)
  if (is.null(method.names)) method.names = names(x$err.test[method.nums])
  ii = method.nums

  # Construct relative test error or relative risk
  N = length(x$err.test) # Number of methods
  err.rel = vector(mode="list",length=N)
  for (j in 1:N) {
    if (what=="error") err.rel[[j]] = x$err.test[[j]] / x$sigma^2
    else err.rel[[j]] = x$risk[[j]] / x$risk.null
  }
  err.obj = tune.and.aggregate(x, err.rel)
  nzs.obj = tune.and.aggregate(x, x$nzs)

  if (type=="ave") {
    xlist = nzs.obj$z.ave
    ylist = err.obj$z.ave
    ybars = err.obj$z.std
  }
  else {
    xlist = nzs.obj$z.med
    ylist = err.obj$z.med
    ybars = err.obj$z.mad
  }

  # Produce the plot
  dat = data.frame(x=unlist(xlist[ii]),
                   y=unlist(ylist[ii]),
                   se=unlist(ybars[ii]),
                   Method=factor(rep(method.names,
                                      lapply(xlist[ii],length))))

  gp = ggplot(dat, aes(x=x,y=y,color=Method)) +
    xlab("Number of nonzero coefficients") +
    ylab(ifelse(what=="error",
                "Relative test error (to Bayes)",
                "Relative risk (to null model)")) +
    geom_line(lwd=lwd) + geom_point(pch=pch) + theme_bw()
  if (std) gp = gp + geom_errorbar(aes(ymin=y-se,ymax=y+se), width=0.2)
  if (!is.null(main)) gp = gp + ggtitle(main)
  if (!legend) gp = gp + theme(legend.pos="none")
  if (make.pdf) ggsave(sprintf("%s/%s.pdf",fig.dir,file.name),
                       height=h, width=w, device="pdf")
  else gp
}

#' Print function for latex-style tables.
#'
#' Print a given table in format digestable by latex.
#'
#' @export print.tex

print.tex = function(tab, tab.se=NULL, digits=3, file=NULL, align="l") {
  tab = round(tab,digits)
  n = nrow(tab); m = ncol(tab)
  rownms = rownames(tab)
  colnms = colnames(tab)
  if (!is.null(tab.se)) {
    tab = matrix(paste0(tab, " (", round(tab.se,digits), ")"), ncol=m)
  }

  if (is.null(file)) file = ""
  cat("", file=file, append=F) # Wipe the file, or the console
  cat(paste0("\\begin{tabular}{|",paste0(rep(align,m+1),collapse="|"),"|}\n"),
      file=file, append=T)
  cat("\\hline\n", file=file, append=T)
  if (!is.null(colnms)) {
    for (j in 1:m) cat(paste0("& ", colnms[j]," "), file=file, append=T)
    cat("\\\\\n", file=file, append=T)
    cat("\\hline\n", file=file, append=T)
  }
  for (i in 1:n) {
    cat(paste0(rownms[i], " "), file=file, append=T)
    for (j in 1:m) cat(paste0("& ",tab[i,j]," "), file=file, append=T)
    cat("\\\\\n", file=file, append=T)
    cat("\\hline\n", file=file, append=T)
  }
  cat(paste0("\\end{tabular}\n"), file=file, append=T)
}