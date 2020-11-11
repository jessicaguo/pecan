##' @name fitA
##' @title fitA
##' @author Mike Dietze
##' @author Xiaohui Feng
##' @export
##' 
##' @param flux.data  data.frame of Licor data, concatenated by rows, and with a leading column 'fname' that is used to count the number of curves and match to covariates
##' @param cov.data   data.frame of covariate data. Column names used in formulas
##' @param model      list including at least 6 components: the fixed effects model for alpha (a.fixed) and Vcmax (V.fixed), the random effects for these (a.random, V.random), the variable used to match the gas-exchange and covariate data (match), and the number of MCMC interations (n.iter). Additional optional arguments: TPU = TRUE turns on TPU limitation; Temp == 'Bernacchi01' turns on the Bernacchi et al 2001 temperature correction. If this is turned on all parameters are estimated for 25C, otherwise no temperature correction is applied. Setting Temp = 'June2004' will turn on the June et al 2004 Funct Plant Biol temperature correction to Jmax. Note: these two corrections are not mutually exclusive, you can set Temp = c('June2004','Bernacchi2001')
##' @param pathway    character of either "C3" or "C4", which specifies which model script gets called to jags.model()
##' @param licor      character of either "6400" or "6800", which specifies whether column names need to be updated
##' 
##' Right now the fixed effects are specified as a string using the standard R lm formula syntax, but without the LHS variable (e.g. '~ SLA + chl + SLA:chl'). The tilde is optional. For random effects, the two options right now are just 'leaf' for leaf-level random effects and NULL. 'model' has a default that sets all effects to NULL (fit one curve to all data) and n.iter=1000.
##' 
fitA <- function(flux.data, cov.data = NULL, model = NULL, pathway, licor) {
  
  ##  TO-DO: 
  ##  Random effects using design matrix
  ##  Model selection
  ##  output variable selection: Pred Loss, WAIC?
  ##  function to do: multiple response curves
  ##  specify priors in model object
  ##  integrate with meta-analysis
  
  library(rjags)
  
  if (is.null(model)) {
    model <- list(V.fixed = NULL, V.random = NULL, 
                  a.fixed = NULL, a.random = NULL, 
                  k.fixed = NULL, k.random = NULL, 
                  r.random = NULL,
                  n.iter = 5000, match = "fname")
  }
  
  ## designate variables to monitor
  if(pathway == "C3"){
    out.variables <- c("r0", "vmax0", "alpha0", "Jmax0", "cp0", "tau", "pmean", "pA")} else if(pathway == "C4"){
        out.variables <- c("r0", "vmax0", "alpha0", "k0", "tau", "pmean", "pA")
                           #, "tau.r0")
      }

  
  ## Determine fixed and random effects
  V.fixed <- model$V.fixed
  V.random <- model$V.random
  a.fixed <- model$a.fixed
  a.random <- model$a.random
  k.fixed <- model$k.fixed
  k.random <- model$k.random
  r.random <- model$r.random
  
  if (is.null(model$match)) {
    model$match <- "fname"
  }
  
  ## convert flux.data to data.frame
  dat <- as.data.frame(flux.data)
  if(licor == "6800"){
    ind <- which(colnames(dat) %in% c("A", "Pci", "Qin"))
    colnames(dat)[ind] <- c("Photo", "Ci", "PARi")
  }
  
  
  ## obtain indexing variables
  id         <- dat[, model$match]
  n.curves   <- length(unique(id))
  curve.id   <- as.numeric(as.factor(id))
  curve.code <- tapply(as.character(id), curve.id, unique)
  
  ## if covariate data present, select only covariates associated with input flux.data
  if (!is.null(cov.data)) {
    ord <- match(curve.code, as.character(cov.data[, model$match]))
    cov.data <- cov.data[ord, ]
  }
  
  ## Adding design matrices for the Vcmax, alpha, and k fixed effects (if present)
  ## Vcmax design matrix
  if (is.null(V.fixed)) {
    XV <- NULL
  } else {
    if (is.null(cov.data)) {
      print("Vcmax formula provided but covariate data is absent:", V.fixed)
    }
    if (length(grep("~", V.fixed)) == 0) {
      V.fixed <- paste("~", V.fixed)
    }
    XV      <- with(cov.data, model.matrix(formula(V.fixed)))
    XV.cols <- colnames(XV)
    XV.cols <- XV.cols[XV.cols != "(Intercept)"]
    XV      <- as.matrix(XV[, XV.cols])
    colnames(XV) <- XV.cols
    Vcenter <- apply(XV, 2, mean, na.rm = TRUE)
    XV      <- t(t(XV) - Vcenter)
  }
  
  ## alpha design matrix
  if (is.null(a.fixed)) {
    Xa <- NULL
  } else {
    if (is.null(cov.data)) {
      print("alpha formula provided but covariate data is absent:", a.fixed)
    }
    a.fixed <- ifelse(length(grep("~", a.fixed)) == 0, paste("~", a.fixed), a.fixed)
    Xa      <- with(cov.data, model.matrix(formula(a.fixed)))
    Xa      <- as.matrix(Xa[, -which(colnames(Xa) == "(Intercept)")])
    acenter <- apply(Xa, 2, mean, na.rm = TRUE)
    Xa      <- t(t(Xa) - acenter)
  }
  
  ## k design matrix
  if (is.null(k.fixed)) {
    Xk <- NULL
  } else {
    if (is.null(cov.data)) {
      print("k formula provided but covariate data is absent:", k.fixed)
    }
    k.fixed <- ifelse(length(grep("~", k.fixed)) == 0, paste("~", k.fixed), k.fixed)
    Xk      <- with(cov.data, model.matrix(formula(k.fixed)))
    Xk      <- as.matrix(Xk[, -which(colnames(Xk) == "(Intercept)")])
    kcenter <- apply(Xk, 2, mean, na.rm = TRUE)
    Xk      <- t(t(Xk) - kcenter)
  }
  
## Specify between the C3 FvCB_RE model and the C4 Collatz_RE model
filename <- ifelse(pathway == "C3", "~/pecan/modules/photosynthesis/R/FvCB_RE.R", 
                   "~/pecan/modules/photosynthesis/R/Collatz_RE.R")
my.model <- readChar(filename, file.info(filename)$size)

## prepare list of data/constants/indices needed as model inputs
sel <- seq_len(nrow(dat))  #which(dat$spp == s)
if("Tleaf" %in% names(dat)){
  if(max(dat$Tleaf) < 100){ # if Tleaf in C, convert to K
    dat$Tleaf <- dat$Tleaf + 273.15
  } 
} else if (!"Tleaf" %in% names(dat)) {
    dat$Tleaf <- 25 + 273.15 ## if no Tleaf, assume 25C in Kelvin
    warning("No Leaf Temperature provided, setting to 25C\n",
            "To change add a column named Tleaf to flux.data data frame")
}

mydat <- list(an = dat$Photo[sel], 
              pi = dat$Ci[sel], 
              q = dat$PARi[sel],
              T = dat$Tleaf, 
              n = length(sel), 
              Kc = 46, ## Michaelis constant CO2 (Pa)
              Ko = 22000,  ## Michaelis constant O2  (Pa)
              po = 21000, ## partial pressure of O2  (Pa)
              rep = curve.id, 
              nrep = n.curves, 
              rd = abs(mean(dat$Photo[round(dat$PARi[sel], 2) == 0])))
                  

## TPU Limitation
if ("TPU" %in% names(model)) {
  if (model$TPU == TRUE) {
    my.model <- gsub(pattern = "#TPU", " ", my.model)
    out.variables <- c(out.variables, "tpu")
  }
}

## Temperature scaling
if ("Temp" %in% names(model)) {
  if ("Bernacchi01" %in% model$Temp) {
    my.model <- gsub(pattern = "##B01", " ", my.model)
  }
  if ("June2004" %in% model$Temp) {
    my.model <- gsub(pattern = "##J04", " ", my.model)
  }
}

## Adding fixed and random effect model formulations for Vcmax, alpha, and k
## Vcmax Formulas
Vformula <- NULL
# Designating random effects
if ("leaf" %in% V.random) {
  Vformula <- " + Vleaf[rep[i]]"
  my.model <- gsub(pattern = "#RLEAF.V", " ", my.model) # adds random effect priors
  out.variables <- c(out.variables, "tau.Vleaf", "Vleaf", "V") # adds "tau.Vleaf" and "Vleaf" to monitored variables
}
# Designating fixed effects based on covariate matrix XV
if (!is.null(XV)) {
  Vnames <- gsub(" ", "_", colnames(XV))
  Vformula <- paste(Vformula,
                    paste0("+ betaV", Vnames, "*XV[rep[i],", seq_along(XV), "]", collapse = " "))
  Vpriors <- paste0("     betaV", Vnames, "~dnorm(0,0.001)", collapse = "\n")
  my.model <- sub(pattern = "## Vcmax BETAs", Vpriors, my.model) # adds fixed effect priors
  mydat[["XV"]] <- XV # adds covariate matrix to the input data list
  out.variables <- c(out.variables, paste0("betaV", Vnames)) # adds fixed effects to monitored variables
}
# Updates model formula to account for random and/or fixed effects
if (!is.null(Vformula)) {
  my.model <- sub(pattern = "#VFORMULA", Vformula, my.model) 
} 

## alpha Formulas
Aformula <- NULL
# Designating random effects
if ("leaf" %in% a.random) {
  Aformula <- " + Aleaf[rep[i]]"
  my.model <- gsub(pattern = "#RLEAF.A", "", my.model) # adds random effect priors
  out.variables <- c(out.variables, "tau.Aleaf", "Aleaf", "A") # adds "tau.Aleaf" and "Aleaf" to monitored variables
}
# Designating fixed effects based on covariate matrix Xa
if (!is.null(Xa)) {
  Anames <- gsub(" ", "_", colnames(Xa))
  Aformula <- paste(Aformula, paste0("+ betaA", Anames, "*Xa[rep[i],", 1:ncol(Xa), 
                                     "]", collapse = " "))
  apriors <- paste0("betaA", Anames, "~dnorm(0,0.001)", collapse = "\n")
  my.model <- sub(pattern = "## alpha BETAs", apriors, my.model) # adds fixed effect priors
  mydat[["Xa"]] <- Xa # adds covariate matrix to the input data list
  out.variables <- c(out.variables, paste0("betaA", Anames)) # adds fixed effects to monitored variables
}
# Updates model formula to account for random and/or fixed effects
if (!is.null(Aformula)) {
  my.model <- sub(pattern = "#AFORMULA", Aformula, my.model)
}

## k Formulas
Kformula <- NULL
# Designating random effects
if ("leaf" %in% k.random) {
  Kformula <- " + Kleaf[rep[i]]"
  my.model <- gsub(pattern = "#RLEAF.K", "", my.model) # adds random effect priors
  out.variables <- c(out.variables, "tau.Kleaf", "Kleaf", "K") # adds "tau.Kleaf" and "Kleaf" to monitored variables
}
# Designating fixed effects based on covariate matrix Xk
if (!is.null(Xk)) {
  Knames <- gsub(" ", "_", colnames(Xk))
  Kformula <- paste(Kformula, paste0("+ betaK", Knames, "*Xk[rep[i],", 1:ncol(Xk), 
                                     "]", collapse = " "))
  kpriors <- paste0("betaK", Knames, "~dnorm(0,0.001)", collapse = "\n")
  my.model <- sub(pattern = "## k BETAs", kpriors, my.model) # adds fixed effect priors
  mydat[["Xk"]] <- Xk # adds covariate matrix to the input data list
  out.variables <- c(out.variables, paste0("betaK", Knames)) # adds fixed effects to monitored variables
}
# Updates model formula to account for random and/or fixed effects
if (!is.null(Kformula)) {
  my.model <- sub(pattern = "#KFORMULA", Kformula, my.model)
}

## r Formulas
Rformula <- NULL
# Designating random effects
if ("leaf" %in% r.random) {
  Rformula <- " * Rleaf[rep[i]]"
  my.model <- gsub(pattern = "#RLEAF.R", "", my.model) # adds random effect priors
  out.variables <- c(out.variables, "tau.Rleaf", "Rleaf", "R") # adds "tau.Rleaf" and "Rleaf" to monitored variables
}
# Updates model formula to account for random and/or fixed effects
if (!is.null(Rformula)) {
  my.model <- sub(pattern = "#RFORMULA", Rformula, my.model)
}

## Define initial conditions
if(pathway == "C3"){
  init <- list()
  init[[1]] <- list(r0 = 1.2, vmax0 = 39, alpha0 = 0.25, tau = 10, cp0 = 6, Jmax0 = 80)  ## tau.Vleaf=30,beta1=4, beta2=1,beta5=3,tau.Vmon=10,tpu=10,
  init[[2]] <- list(r0 = 1, vmax0 = 100, alpha0 = 0.2, tau = 20, cp0 = 4, Jmax0 = 150)  ##tau.Vleaf=20,beta1=1,beta2=1,beta5=-1,tau.Vmon=20,tpu=13,
  init[[3]] <- list(r0 = 2, vmax0 = 60, alpha0 = 0.28, tau = 20, cp0 = 5, Jmax0 = 60)  ##tau.Vleaf=100,beta1=1,beta2=2,beta5=2,tau.Vmon=3,tpu=20,
} else if(pathway == "C4"){
  init<-list()
  init[[1]]<-list(r0 = 0.5, vmax0 = 30, alpha0 = 0.03, tau = 10, k0 = 0.7*100000)   ## ,tau.Vleaf=300,tau.Kleaf=1e-10, beta1=4, beta2=1,beta5=3,tau.Vmon=10
  init[[2]]<-list(r0 = 0.75, vmax0 = 20, alpha0 = 0.07, tau = 20, k0 = 0.8*100000)    ## ,tau.Vleaf=200,tau.Kleaf=2e-9,beta1=1,beta2=1,beta5=-1,tau.Vmon=20
  init[[3]]<-list(r0 = 0.8, vmax0 = 15, alpha0 = 0.06, tau = 20, k0 = 0.2*1000000)    ## ,tau.Vleaf=100,tau.Kleaf=3e-8,beta1=1,beta2=2,beta5=2,tau.Vmon=3}
}
 
mc3 <- jags.model(file = textConnection(my.model), data = mydat, inits = init, n.chains = 3)

mc3.out <- coda.samples(model = mc3, variable.names = out.variables, n.iter = model$n.iter)

## split output
out         <- list(params = NULL, predict = NULL, model = my.model)
mfit        <- as.matrix(mc3.out, chains = TRUE)
pred.cols   <- union(grep("pA", colnames(mfit)), grep("pmean", colnames(mfit)))
chain.col   <- which(colnames(mfit) == "CHAIN")
out$predict <- mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
out$params  <- mat2mcmc.list(mfit[, -pred.cols])
return(out)
} # fitA


##' @name read_Licor
##' @title read_Licor
##' 
##' @author Mike Dietze
##' @export
##' 
##' @param filename  name of the file to read
##' @param sep       file delimiter. defaults to tab
##' @param ...       optional arguements forwarded to read.table
read_Licor <- function(filename, sep = "\t", ...) {
  fbase <- sub(".txt", "", tail(unlist(strsplit(filename, "/")), n = 1))
  print(fbase)
  full <- readLines(filename)
  ## remove meta=data
  start <- grep(pattern = "OPEN", full)
  skip <- grep(pattern = "STARTOFDATA", full)
  for (i in rev(seq_along(start))) {
    full <- full[-(start[i]:(skip[i] + 1 * (i > 1)))]  # +1 is to deal with second header
  }
  full <- full[grep("\t", full)]  ## skip timestamp lines
  dat <- read.table(textConnection(full), header = TRUE, blank.lines.skip = TRUE,
                    sep = sep, ...)
  fname <- rep(fbase, nrow(dat))
  dat <- as.data.frame(cbind(fname, dat))
  return(dat)
} # read_Licor


mat2mcmc.list <- function(w) {
  temp <- list()
  chain.col <- which(colnames(w) == "CHAIN")
  for (i in unique(w[, "CHAIN"])) {
    temp[[i]] <- as.mcmc(w[w[, "CHAIN"] == i, -chain.col])
  }
  return(as.mcmc.list(temp))
} # mat2mcmc.list
