# These are the R routines for the package mgcv (c) Simon Wood 2000

QT<-function(A) {

# This routine performs a QT factorization of matrix A. A QT factorization
# is of the form AQ=[0,T] where T is reverse lower triangular (top left zero)
# nrow(A)<ncol(A) and the first (ncol(A)-nrow(A)) columns of Q are the null
# space of A. Q is actually stored as Householder rotations in over A. 
# Specifically if u_i is row i of A then let H_i=(I-u_i u_i') then 
# Q = H_1 H_2 ... there are nrow(A) u_i's and H_i's. The last i-1 elements
# of u_i are always zero (i>=1). The main use of this routine is to find
# the null space of a constraint matrix C, suitable for input to mgcv().
# The return value is the over-written A. The factorization os performed
# using compiled code.

oo<-.C("RQT",as.double(A),as.integer(nrow(A)),as.integer(ncol(A)))
A<-matrix(oo[[1]],nrow(A),ncol(A))
A
}

GAMsetup<-function(G) {
# This function sets up the design matrix and penalty matrices for a GAM formulated with
# one dimensional cubic penalized regression splines. It takes a list G as argument which
# should have the following named elements:
#
# G$m the number of smooth terms in the model
# G$n the number of data to be modelled
# G$nsdf the number of user supplied columns of the design matrix for any parametric model parts
#        including the constant (if any)
# G$df an array of G$m integers specifying the maximum d.f. for each spline term
# G$x an array of G$n element arrays of data and (optionally) design matrix columns. The first
#     G$nsdf elements of G$x should contain the elements of the columns of the design matrix
#     corresponding to the parametric part of the model. The remaining G$m elements of G$x 
#     are the values of the covariates that are arguments of the spline terms.
#
# The function returns a list H containing all the above named elements plus the following:
#
# H$X the full design matrix.
# H$S an array of matrices containing the coefficients of the penalties. These are stored
#     in a compact form, so that H$S[i] is the smallest square matrix containing submatrix
#     conatining all the non-zero elements of S_i the ith penalty matrix. Element 0,0 of
#     H$S[i] is element off[i],off[i] of S_i, element 0,1 of H$S[i] is element off[i],off[i]+1
#     of S_i, and so on.
# H$off is an array of offsets, used to facilitate efficient storage of the penalty matrices
#       and to indicate where in the overall parameter vector the parameters of the ith 
#       spline reside (e.g. first parameter of ith spline is at p[off[i]]).
# H$C - the matrix of linear equality constraints on the parameters .  
# H$xp -matrix whose rows contain the covariate values corresponding to the parameters
#      of each spline - the splines are parameterized using their y- values at a series
#      of x values - these vectors contain those x values!

  q<-G$nsdf
  for (i in 1:G$m) { q<-q+G$df[i] }  # q stores total number of parameters
  X<-matrix(0,G$n,q)             # design matrix
  mdf<-max(G$df)
  xp<-matrix(0,G$m,mdf)          # covariate matrix
  S<-array(0,dim=c(G$m,mdf,mdf))   # array of wiggliness penalty matrices
  C<-matrix(0,G$m,q)             # fixed constraint matrix (one per spline term)
  off<-array(0,c(G$m))           # array for storing offsets for Wiggliness penalties
  o<-.C("RGAMsetup",as.double(X),as.double(C),as.double(S),as.double(xp),
        as.integer(off),as.double(G$x),as.integer(G$m),as.integer(G$n),as.integer(G$df),
        as.integer(G$nsdf)) # compiled code to set up matrices for model
  X<-matrix(o[[1]],G$n,q);
  C<-matrix(o[[2]],G$m,q);
  S<-array(o[[3]],dim=c(G$m,mdf,mdf));
  xp<-matrix(o[[4]],G$m,mdf);
  off<-array(o[[5]],c(G$m));
  G$X<-X;G$S<-S;G$off<-off;G$C<-C;G$xp<-xp
  G # ... and return
}


mgcv<-function(M) {

# Performs multiple smoothing parameter selection for Generalized ridge regression problems 
# M is a list which must containing the following named elements:
#
# M$y  - the response variable vector
# M$X  - the design matrix (declared as a matrix, with correct numbers of rows and columns,
#        as number of data and number of parameters are read from this)
# M$C  - matrix defining linear equality constraints on parameters (Cp=0). 
#        Number of rows is number of constraints.
# M$w  - weight vector (often proportional to inverse of variance)
# M$S  - A 3 dimensional array containing the penalty coefficient matrices. dim(M$S)[1], should
#        be the number of penalties. S[i] is the 2-d array containing the coeffs for the ith
#        penalty. S[i][j][k] contains element [j+M$off[i]][k+M$off[i]] of the ith penalty matrix
#        S_i. This slightly complicated way of storing things saves memory.
# M$off - array of offsets used as described for M$S and also indicates where in the parameter
#         vector (p, say) the parameters penalized by the ith penalty start.
# M$df  - array containing actual dimensions of non-zero part of S[i], i.e. S[i] contains only
#         an M$df square block of non-zero elements.
# M$fix - array containing TRUE if smooth i is to have fixed d.f. and FALSE otherwise
# M$sig2 - the error variance if it is known and to be used for UBRE - must be <=0 to use GCV.
#          Note the convention used here is that var(y_i) = sig2/w_i.
# M$sp   - the initial smoothing parameter estimates: if any of these are <=0 then automatic
#           initialisation is used
#
# The routine returns M with the following elements added (or reset):
#
# M$p    - the best fit parameter vector given the selected smoothing parameters
# M$sp   - the vector of smoothing parameters (\theta_i/\rho) estimated by the method
# M$sig2 - the estimate of the error variance (if GCV used)
# M$Vp   - the estimated covariance matrix of the parameters set Vp[1,1] <0 to not calculate
# M$edf  - the estimated degrees of freedom for the ith smooth term if Vp calculated
  
  if (is.null(M$sig2)) M$sig2<-0
  C.r<-nrow(M$C)          # number of equality constraints
  if (is.null(C.r)) C.r<-0
  q<-ncol(M$X)            # number of parameters
  
  n<-nrow(M$X)            # number of data
  # need to strip out penalties for fixed df smooths..... 
  k<-dim(M$S)[1]             # number of penalties (before stripping)
  if (is.null(k)) k<-0
  if (k!=0)  m<-k-sum(M$fix) # count penalty terms to include
  else m<-0
  S<-array(0,c(m,dim(M$S)[2],dim(M$S)[3])) # reduced smoothness array
  off<-array(0,m)  # reduced offset array
  df<-array(0,m)   # reduced penalty size array
  j<-1             # element index for off[], df and S[]    
  if (k!=0)
  for (i in 1:k)   # do actual stripping
  if (M$fix[i]==FALSE) 
  { S[j,,]<-M$S[i,,]
    off[j]<-M$off[i]
	df[j]<-M$df[i]
    j<-j+1
  }
  # deal with quantities that will be estimated
  p<-matrix(0,q,1)      # set up parameter vector
  Vp<-matrix(0.0,q,q)   # estimated covariance matrix
  edf<-array(0,m)       # estimated degrees of freedom
  Vp[1,1]<-1.0

  oo<-.C("mgcv",as.double(M$y),as.double(M$X),as.double(M$C),as.double(M$w),as.double(S),
         as.double(p),as.double(M$sp),as.integer(off),as.integer(df),as.integer(m),
         as.integer(n),as.integer(q),as.integer(C.r),as.double(M$sig2),as.double(Vp),
		 as.double(edf))
   
  p<-matrix(oo[[6]],q,1);
  sig2<-oo[[14]]
  Vp<-matrix(oo[[15]],q,q)
  sp<-matrix(oo[[7]])
  edf<-oo[[16]]
  # unpack results back to correct place in output (which includes fixed d.f. and free d.f. terms)
  M$sp<-array(0,k)
  M$edf<-array(0,k)
  j<-1
  if (k!=0)
  for (i in 1:k)
  { if (M$fix[i]==TRUE)
    { M$sp[i]<-0
	  M$edf[i]<- -1 # must be filled in outside this routine, with reference to constraints
    }
	else
	{ M$sp[i]<-sp[j]
	  M$edf[i]<-edf[j]
      j<-j+1
    }
  }
  M$Vp<-Vp
  M$p<-p;
  M$sig2<-sig2;
  M    # return list
}


SANtest<-function(n=100,sig2=-1) { 

# Example function that simulates data from a simple additive normal truth and then attempts
# to reconstruct it using a GAM constructed from penalized regression splines. 
# n is number of data. Magnitude of sig2 is error variance used for simulation.
# sig2<=0 => use GCV, otherwise use UBRE for smoothing parameter estimation.

  if (sig2<=0) noquote("GCV used") else noquote("UBRE used")
  if (n<60) 
  { n<-60
    warning("n reset to 60")
  }
# simulate some data to which GAM can be fitted......

  x<-runif(5*n,0,1)      # simulate covariates
  x<-array(x,dim=c(5,n)) # load into array
  pi<-asin(1)*2          
  y<-2*sin(pi*x[2,])  
  y<-y+exp(2*x[3,])-3.75887
  y<-y+0.2*x[4,]^11*(10*(1-x[4,]))^6+10*(10*x[4,])^3*(1-x[4,])^10-1.396 # simulate truth
  e<-rnorm(n,0,sqrt(abs(sig2)))  # noise to add to truth
  y<-y+e              # create data from truth + noise   
  w<-matrix(1,n,1)    # weight matrix
  

# display scatter plots of data against covariates.....

  par(mfrow=c(2,2))   # set display environment to 2 by 2
  plot(x[2,],y)      
  plot(x[3,],y)
  plot(x[4,],y)
  plot(x[5,],y)

# set up input list for GAMsetup and call GAMsetup.......

  x[1,]<-1 # set up dummy variable for intercept term
  G<-list(m=4,n=n,nsdf=1,df=c(15,15,15,15),x=x) # input list for GAMsetup
  H<-GAMsetup(G)  
  H$fix<-array(FALSE,H$m)
  H$y<-y;H$sig2<-sig2;H$w<-w # add data, variance and weights to mgcv input list
  H$sp<-array(-1,H$m)
  H<-mgcv(H)   # fit model
  xp<-H$xp
  p<-H$p;
  sig2<-H$sig2

  readline("Hit return to see fitted model and truth")
# plot reconstruction and truth....
  i<-c((H$off[1]+1):H$off[2]);yy<-p[i];yy<-yy-mean(yy);i<-i-H$off[1];xx<-xp[1,i];
  plot(xx,yy,type="l");t1<-2*sin(pi*xx);t1<-t1-mean(t1);lines(xx,t1); # t1 is truth
  i<-c((H$off[2]+1):H$off[3]);yy<-p[i];yy<-yy-mean(yy);i<-i-H$off[2];xx<-xp[2,i];
  plot(xx,yy,type="l");t2<-exp(2*xx)-3.75887;t2<-t2-mean(t2);lines(xx,t2);
  i<-c((H$off[3]+1):H$off[4]);yy<-p[i];yy<-yy-mean(yy);i<-i-H$off[3];xx<-xp[3,i];
  plot(xx,yy,type="l"); 
  t3<-0.2*xx^11*(10*(1-xx))^6+10*(10*xx)^3*(1-xx)^10-1.396;t3<-t3-mean(t3); lines(xx,t3);
  i<-c((H$off[4]+1):ncol(H$X));yy<-p[i];yy<-yy-mean(yy);i<-i-H$off[4];xx<-xp[4,i];
  plot(xx,yy,type="l");t4<-xx*0.0;lines(xx,t4);
  sig2 # the return value
}



########################################################################################
## The following provide a gam() function and associated gubbins for R......
########################################################################################

gam.parser<-function(formula)
# parses a gam formula of the generic form:
#   y~x0+x1+x3*x4 + s(x5)+ s(x6,x7) ....
# and returns 2 model formulae: 1 containing only the non-smooth stuff
# and the other relating to the smooths in a conventional model notation 
# without s(,,,)'s e.g:
#   y~x5+x6:x7
# On exit a list is returned (see ret below). The parametric part may also 
# include any model offsets.
# Smooth terms may include a max df as last argument e.g. s(x1,x2,x3,45).
#


{ ft<-as.character(formula)  # array of model formula components
  # unfortunately the above will truncate the rhs of the model formula if too long
  work.around<-as.character(deparse(formula[3]))
  ft[3]<-""    #  delete original, possibly defective, ft[3]
  for (i in 1:length(work.around)) ft[3]<-paste(ft[3],work.around[i],sep="")
  rm(work.around)           # ft[3] now contains rhs of formula + '()'
  sf<-paste(ft[2],"~",sep="") # smooth formula (as text)
  pf<-paste("~",sep="") # parametric formula (as text)
  i<-1       # character counter
  sm<-0      # 1 if currently outputting to smooth formula 0 for parametric formula
  suspend<-0 # 1 if output suspended because a "+" has occurred (could relate to para or smooth!) 
  smooth.started<-0 
  param.started<-0
  reading.df<-0
  n.smooths<-0
  lc<-" "
  c<-" "
  ndf<-c(-1)
  fix<-c(FALSE)
  options(warn=-1) # turn warnings off to prevent warning on is character a number test
  while(i <= nchar(ft[3])-2) # work through rhs of model (excluding final '()')
  { if (c!=" ") lc<-c # store last non-space
    c<-substring(ft[3],i,i) # next char from string
    if (c==" ") # then continue
    { i<-i+1
    } else
    { if (c=="+")    # stop o/p of parametric terms
      { suspend<-1
        i<-i+1
      } else # neither " " nor "+"
      { if (pmatch("s(",substring(ft,i,i+1),nomatch=0)) # then it's a smooth term
        { sm<-1     # so switch to smooth mode
          i<-i+2      
          if (smooth.started) sf<-paste(sf,"+",sep="")
          smooth.started<-1
          n.smooths<-n.smooths+1
		  fix[n.smooths]<-FALSE # default - d.f. not fixed
        } else # not the beginning of a smooth term
        { i<-i+1
          if (sm==0) # pass it to the parametric model
          { if (suspend&&param.started)
            { suspend<-0 
              pf<-paste(pf,"+",c,sep="")  # add term and "+" here
            } else    # add term
            { suspend<-0
              pf<-paste(pf,c,sep="")  # add term
            }
            param.started<-1
          } else # pass to smooth model unless it's ")"
          { if (c==")")   # end of smooth term
            { suspend<-0 
              sm<-0 
              if (reading.df)
              { # finish number and dump it to df[]
                if (fix[n.smooths]&&lc!="f"&&lc!="F") stop("Error in formula")
				ndf[n.smooths]<-as.numeric(dfs)     
                reading.df<-0
              } else
              { ndf[n.smooths]<-0 # no df given             
			  }
            } else
            { if (reading.df)
              { if (c=="|"||c=="f"||c=="F")
			    fix[n.smooths]<-TRUE
				else
				{ if (lc=="|"&&c!="f"&&c!="F") stop("Error in formula")
				  else
			      dfs<-paste(dfs,c,sep="")
                }
              } else
              { if (lc==","&& is.na(as.numeric(c))==FALSE ) # then this is start of df term
                { reading.df<-1
                  dfs<-c
                } else
                { if (lc==",") # replace with ":"
                  { sf<-paste(sf,":",sep="")
                    sf<-paste(sf,c,sep="")  
                  }
                  else 
                  { if (c!=",")
                    sf<-paste(sf,c,sep="")
                  } 
                }
              }
            }  
          }  
        }
      }
    }
  }
  options(warn=0) # turn warnings back on 
  if (param.started)
  { if (length(attr(terms(as.formula(pf)),"variables"))==1) # then formula contains only "1"'s or "-1"'s
    { if (length(grep("-",pf))>=1) # then no intercept
      param.started<-0  # there should be no parametric terms
      else pf <- "~ 1"    # make sure formula is unambiguous - there is an intercept only!
    } 
  } else # an intercept term is required
  { param.started<-1
    pf <- "~ 1"
  }
  # end of parse loop
  ret<-list(sfok=smooth.started,pfok=param.started,df=ndf,sftext=sf,pftext=pf,fix=fix)
  if (smooth.started) ret$sf<-as.formula(sf)
  if (param.started) ret$pf<-as.formula(pf)
  ret
}

gam.setup<-function(formula,data=list(),gam.call,predict=TRUE,parent.level=2,nsdf=-1)

# This gets the data referred to in the model formula, either from data frame "data"
# or from the level parent.level parent (e.g. when called from gam() this will be
# the parent that called gam(). 
# G$names[i] contains the names associated with each column of the design matrix,
# followed by the name of each smooth term - the constant is un-named! 
# Note that it is assumed that the default parent.level assumes that this function
# is called from the function that was called by the user! i.e. the user's calling 
# environment is the grandparent of the environment of gam.setup() 
# if nsdf==-1 then it assumed that the design matrix for the non-spline part of the 
# model is to be found here. Otherwise it is assumed that nsdf is the known 
# number of degrees of freedom for the non-spline part (including intercept), but
# the columns corresponding to the non-spline part of the model and the offset 
#  are set to zero.
# 
# NOTE: all.vars(formula) will return the names of all the variables refered to by 
# a formula...... use of this could make some of the code in this function somewhat
# simpler!

{ # now split the formula
  split<-gam.parser(formula) 
  dmiss<-missing(data)  
  if (dmiss) data<-sys.frame(sys.parent(n=parent.level))
  if (split$df[1]==-1) stop("You've got no smooth terms - use another model fitter")
  
  mt<-terms(split$sf)
  
  m<-length(split$df) # number of smooth terms
  G<-list(m=m,df=split$df)
  G$fix<-split$fix
  add.constant<-FALSE
  if (split$pfok)   # deal with the strictly parametric part of the model 
  { if (length(attr(terms(split$pf),"variables"))==1) # then the formula is ~ 1 and model.frame(), etc can't cope 
    { G$nsdf<-1
      add.constant<-TRUE  # must construct model matrix later when length of response is known
    } else
    { if (nsdf>=0)  # set G$nsdf to supplied value, but deal with X later
      { G$nsdf<-nsdf
      } else 
      { if (predict) # then user MUST have supplied all data in data frame and searching is easy
        { mf<-model.frame(split$pf,data)
        } else   # all sorts of stupid practices are allowed and must be catered for.... 
        { mf<-gam.call # assemble a call to model frame which will allow searching in the calling env. 
          mf[[1]]<-as.name("model.frame")      # replace "gam" with "model.frame" in call
          mf$family<-mf$weights<-mf$control<-mf$scale<-NULL
          mf$formula<-split$pf                 # set formula to parametric only part     
          # evaluate in the grandparent environment to get model frame......
          mf<-eval(mf,sys.frame(sys.parent(n = parent.level))) 
        }
        G$offset <- model.offset(mf)   # get the model offset (if any)
        if (!is.null(G$offset) && length(attr(terms(split$pf),"variables"))<=2) # then there is an offest but no terms other than "+1" or "-1"
        { if (length(grep("-",as.character(split$pf[2])))>=1) # then there is a "-1" term in formula
          X<-matrix(0,length(G$offset),0)
          else X <- model.matrix(split$pf,mf)  # then there is a constant in the model and model.matrix() can cope 
        }
        else X <- model.matrix(split$pf,mf)       # getting the model matrix
        G$nsdf <- dim(X)[2]                  # extract the total non- smooth df
      }        
    }	
  } else
  G$nsdf<-0
  
  xvars <- as.character(attr(mt, "variables"))[-1]
  if (!predict) # obtain the "y variable"  
  { # evaluates rhs of model using data and failing that the calling environment...
    G$y<-eval(terms(split$sf)[[2]],data,sys.frame(sys.parent(n = parent.level)))
    # checks that this has worked.....
    if (is.null(G$y)) stop(paste("Failed to find variable for",xvars[1]))
  }

  xfacs <- attr(mt, "factors")
  term.labels<-labels(mt)
 
  # now find number of data without looking at y (needed for predicting)
  
  label<-term.labels[1]

  vlist<-xvars[as.logical(xfacs[, label])]
 
  # need to know numer of rows in design matrix .....
  z<-eval(parse(text=vlist[1]),data,sys.frame(sys.parent(n = parent.level))) 
  #if (dmiss||is.null(data[[vlist[1]]])) z<-get(vlist[1],envir=sys.frame(sys.parent(n=parent.level)))
  #else z<-data[[vlist[1]]] code changed 22/2/01 

  G$n<-NROW(z)
  
  if (add.constant)          # then X must be created by hand for case where parametric model is ~1
  { X<-matrix(1,G$n,G$nsdf)
    colnames(X)<-"constant"
  }
 
  if (nsdf>=0) # then X still needs to be filled
  { X<-matrix(0,G$n,G$nsdf)
    colnames(X)<-rep("",G$nsdf)         # don't need them, but will fail below without!  
    G$offset<-array(0,G$n)
  }
  G$x<-array(0,dim=c(G$nsdf+m,G$n)) # explanatory data array
  G$names<-array("",G$nsdf+m)
  if (G$nsdf>0) for (k in 1:G$nsdf) 
  { G$x[k,]<-X[,k] # copy columns of parametric design matrix
    G$names[k]<-colnames(X)[k] # copy column names for later use as labels
  }
  if (G$nsdf>0) G$vnames<-c(colnames(X),term.labels) else
  G$vnames<-term.labels
  kk<-k<-G$nsdf+1
  for (label in term.labels)
  { vlist<-xvars[as.logical(xfacs[, label])]
	G$names[kk]<-"s("
    for (i in 1:length(vlist)) # produces a column for each variable in this smooth
    { z<-eval(parse(text=vlist[i]),data,sys.frame(sys.parent(n = parent.level)))
      #if (dmiss||is.null(data[[vlist[i]]]))  # then search the user's calling environment  
      #z<-get(vlist[i],envir=sys.frame(sys.parent(n=parent.level)))
      #else z<-data[[vlist[i]]]     # data should be in data frame  - old code changed 22/2/01
      if (is.null(z)) stop(paste("Failed to find variable",vlist[i]))  
      G$x[k,]<-z
      if (i<length(vlist)) G$names[kk]<-paste(G$names[kk],vlist[i],",",sep="")
      else G$names[kk]<-paste(G$names[kk],vlist[i],sep="")
      k<-k+1
    }
    if (length(vlist)>1) stop("Sorry - gam() can parse multidimensional smooths, but can't fit them yet")
    G$names[kk]<-paste(G$names[kk],")",sep="")
    kk<-kk+1
  }
  # now sort out max dfs per term (only if get.y - otherwise call from predict)
  if (!predict)
  { k<-0
    for (i in 1:G$m)
    { if (G$df[i]==2&&G$fix[i]==FALSE)
	  { G$fix[i]<-TRUE
	    warning("2 knot natural spline is a straight line, reset to fixed d.f.")
	  }
	  if (G$df[i]==1) 
          stop("You can't have a GAM component with one knot and zero degrees of freedom")
	  if (G$df[i]==0) G$df[i]<-10;
	  
	  if (length(unique(G$x[i+G$nsdf,]))<G$df[i]) # too many knots
	  { G$df[i]<-length(unique(G$x[i+G$nsdf,])) 
	    warning("Number of knots reset: must not exceed number of unique covariate values.")
	  }  
      k<-k+G$df[i] 
    }
    if (k>G$n) stop("You have specified more degrees of freedom than you have data")
  }
  G
}


gam<-function(formula,family=gaussian(),data=list(),weights=NULL,control=gam.control(),scale=0)

# Routine to fit a GAM to some data. The model is stated in the formula, which is then 
# parsed to figure out which bits relate to smooth terms and which to parametric terms.
# 

{ gam.call<-match.call()  # store the call to facilitate searching in gam.setup()
  # coerce family to correct type
  if (is.character(family)) 
  family <- get(family)
  if (is.function(family)) 
  family <- family()
  if (is.null(family$family)) 
  { print(family)
    stop("`family' not recognized")
  }
  if (missing(data)) G<-gam.setup(formula,gam.call=gam.call,predict=FALSE)
  else G<-gam.setup(formula,data,gam.call=gam.call,predict=FALSE)
 
  G<-GAMsetup(G) 
 
  if (is.null(G$offset)) G$offset<-rep(0,G$n)
  if (is.null(weights)) G$w<-rep(1,G$n) else G$w<-weights
  
  if (scale==0) 
  { if (family$family=="binomial"||family$family=="poisson") scale<-1 #ubre
    else scale <- -1 #gcv
  }

  G$sig2<-scale
  G$sp<-array(-1,G$m) # set up smoothing parameters for autoinitialization at first run

  object<-gam.fit(G,family=family,control=control)
  
  object$xp<-G$xp
  # return correct d.f. for fixed d.f. terms
  for (i in 1:G$m) if (G$fix[i]) object$edf[i]<-G$df[i]-1
  # now re-assign variable names to coefficients etc. 
  term.names<-array("",length(object$coefficients))
  #term.names[1]<-"(Intercept)"
  if (G$nsdf>0) for (i in 1:G$nsdf) term.names[i]<-G$names[i]
  i<-G$nsdf+1
  for (k in 1:G$m)
  for (j in 1:G$df[k])
  { term.names[i]<-paste(G$names[G$nsdf+k],".",as.character(j),sep="")
    i<-i+1
  }
  names(object$coefficients)<-term.names  # note - won't work on matrices!!
  object$formula<-formula
  object$x<-G$x
  object$call<-gam.call
  rownames(object$x)<-G$vnames
  names(object$x)<-NULL
  #object$X<-G$X
  class(object)<-"gam"

  object
}

print.gam<-function (x, ...) 
# default print function for gam objects
{ print(x$family)
  cat("Formula:\n")
  print(x$formula)
  cat("\nEstimated degrees of freedom:\n",x$edf,"  total = ",sum(x$edf)+x$nsdf,"\n")
  gcv<-x$df.null*x$sig2/(x$df.null-sum(x$edf)-x$nsdf)
  cat("\nGCV score: ",gcv,"\n")
}

gam.control<-function (epsilon = 1e-04, maxit = 30, trace = FALSE) 
# control structure for a gam
{   if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of epsilon must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace)
}

gam.fit<-function (G, start = NULL, etastart = NULL, 
    mustart = NULL, family = gaussian(), 
    control = gam.control()) 
# fitting function for a gam, modified from glm.fit
{
    #x <- as.matrix(x)
    #xnames <- dimnames(x)[[2]]
    #ynames <- names(y)
    conv <- FALSE
    nobs <- NROW(G$y)
    nvars <- NCOL(G$X) # check this needed
    y<-G$y # original data
    X<-G$X # original design matrix
    if (nvars == 0) stop("Model seems to contain no terms")
  
    weights<-G$w # original weights
   
    offset<-G$offset 

    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
    valideta <- family$valideta
    if (is.null(valideta)) 
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu)) 
        validmu <- function(mu) TRUE
    if (is.null(mustart)) 
        eval(family$initialize, sys.frame(sys.nframe()))
    
    if (NCOL(y) > 1) 
        stop("y must be univariate unless binomial")
    eta <- if (!is.null(etastart) && valideta(etastart)) 
        etastart
    else if (!is.null(start)) 
        if (length(start) != nvars) 
            stop(paste("Length of start should equal", nvars ))
        else as.vector(if (NCOL(G$X) == 1) 
            G$X * start
        else G$X %*% start)
    else family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta))) 
        stop("Can't find valid starting values: please specify some")
    devold <- sum(dev.resids(y, mu, weights))
    coefold <- start
    boundary <- FALSE
    scale<-G$sig2
    
    for (iter in 1:control$maxit) {
        good <- weights > 0
        varmu <- variance(mu)[good]
       
        if (any(is.na(varmu))) 
            stop("NAs in V(mu)")
        if (any(varmu == 0)) 
            stop("0s in V(mu)")
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val[good]))) 
            stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (mu.eta.val != 0) # note good modified here => must re-calc each iter
        if (all(!good)) {
            conv <- FALSE
            warning(paste("No observations informative at iteration", 
                iter))
            break
        }
   
        z<-G$y <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
        w<-G$w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
        G$w<-G$w^2 # this line is somewhat important
        G$X<-X[good,]  # truncated design matrix       
        ngoodobs <- as.integer(nobs - sum(!good))
        ncols <- as.integer(1)
        # must set G$sig2 to scale parameter or -1 here....
        G$sig2<-scale
        G<-mgcv(G) 
		
        start <- coef <- G$p
       
        eta[good] <- drop(X[good, , drop = FALSE] %*% start)
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace) 
            cat("Deviance =", dev, "Iterations -", iter, "\n")
        boundary <- FALSE
        if (is.na(dev) || any(is.na(coef))) {
            warning("Step size truncated due to divergence")
            ii <- 1
            while ((is.na(dev) || any(is.na(start)))) {
                if (ii > control$maxit) 
                  stop("inner loop 1; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2
                eta[good] <- drop(X[good, , drop = FALSE] %*% 
                  start)
                mu <- linkinv(eta <- eta + offset)
                dev <- sum(dev.resids(y, mu, weights))
            }
            boundary <- TRUE
            coef <- start
            if (control$trace) 
                cat("New Deviance =", dev, "\n")
        }
        if (!(valideta(eta) && validmu(mu))) {
            warning("Step size truncated: out of bounds.")
            ii <- 1
            while (!(valideta(eta) && validmu(mu))) {
                if (ii > control$maxit) 
                  stop("inner loop 2; can't correct step size")
                ii <- ii + 1
                start <- (start + coefold)/2
                eta[good] <- drop(X[good, , drop = FALSE] %*% 
                  start)
                mu <- linkinv(eta <- eta + offset)
            }
            boundary <- TRUE
            coef <- start
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("New Deviance =", dev, "\n")
        }
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
            conv <- TRUE
            break
        }
        else {
            devold <- dev
            coefold <- coef
        }
    }
    if (!conv) 
        warning("Algorithm did not converge")
    if (boundary) 
        warning("Algorithm stopped at boundary value")
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
        if (any(mu > 1 - eps) || any(mu < eps)) 
            warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (family$family == "poisson") {
        if (any(mu < eps)) 
            warning("fitted rates numerically 0 occurred")
    }
    
    residuals <- rep(NA, nobs)
    residuals[good] <- z - (eta - offset)[good]
    
    nr <- min(sum(good), nvars)

    wt <- rep(0, nobs)
    wt[good] <- w^2
   
    wtdmu <- sum(weights * y)/sum(weights)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok 

	list(coefficients = as.vector(coef), residuals = residuals, fitted.values = mu, 
        family = family, 
        linear.predictor = eta, deviance = dev,
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        #df.residual = resdf, 
        df.null = nulldf, y = y, converged = conv,sig2=G$sig2,edf=G$edf,
        boundary = boundary,sp = G$sp,df=G$df,nsdf=G$nsdf,Vp=G$Vp)
}


predict.gam<-function(object,newdata,type="link",se.fit=F,plot.call=FALSE, ...) {

# This function is used for predicting from a GAM. object is a gam object, newdata a dataframe to
# be used in prediction..... it's all done via a call to the compiled C routine RGAMpredict().
#
# Type == "link"     - for linear predictor
#      == "response" - for fitted values
#      == "terms"    - for individual terms on scale of linear predictor 
#
#  Steps are:
#  1. Construct x[i,j] - the explanatory data array for which predictions are required, this means
#     using the model formulae from the gam object, "object" to extract the model variables from 
#     the data-frame "newdata". Note np - the number of elements to predict for.
#  2. Convert the knot position data to a form usable by the C routine RGAMpredict(). 
#  3. Initialize storage space for eta and se.   
#  4. Call RGAMsetup()
#  5. Use eta and se to construct the returned vector, matrix or list.
#  6. Tidy up and return.  

  if (type!="link"&&type!="terms"&&type!="response")  
  { warning("Unknown type, reset to terms.")
    type<-"terms"
  }
  if (class(object)!="gam") stop("predict.gam can only be used to predict from gam objects")

  # get data from which to predict.....  
  
  if (missing(newdata)) 
  { m<-length(object$sp);n<-length(object$y)
    G<-list(x=object$x,nsdf=object$nsdf,m=m,n=n)
    no.data<-FALSE
  }
  else 
  { if (plot.call) # then it's a call from plot, and only spline parts are to be evaluated
    G<-gam.setup(object$formula,newdata,gam.call=object$call,predict=TRUE,nsdf=object$nsdf)
    else G<-gam.setup(object$formula,newdata,gam.call=object$call,predict=TRUE)    
    no.data<-FALSE
  }
  np<-G$n
  
  # now set up the other information for prediction.....
  control<-0
  if (type=="terms") 
  { control<-2
    eta<-array(0,dim=c(G$nsdf+G$m,np))
	se<-array(0,dim=c(G$nsdf+G$m,np))
  } else
  { eta<-array(0,dim=c(np))
    se<-array(0,dim=c(np))
  }
  if (se.fit&&object$Vp[1,1]<=0)
  { se.fit<-FALSE
    warning("No variance estimates available")
  }
  if (se.fit) control<-control+1
 
  # call compiled code for prediction ....... 
  o<-.C("RGAMpredict",as.double(object$xp),as.integer(G$nsdf),as.integer(object$df),
         as.integer(length(object$df)),as.double(G$x),as.integer(np),as.double(object$coefficients),
		 as.double(object$Vp),as.double(eta),as.double(se),as.integer(control))
  # now eta contains linear predictor (terms) and se may contain corresponding standard errors
   
  if (type=="terms")
  { eta<-array(o[[9]],c(G$nsdf+G$m,np));
    se<-array(o[[10]],c(G$nsdf+G$m,np));
  } else
  { eta<-array(o[[9]],c(np));
    se<-array(o[[10]],c(np));
    if (type=="response") # transform onto scale of data
    { fam<-object$family;linkinv<-fam$linkinv;dmu.deta<-fam$mu.eta
      se<-se*abs(dmu.deta(eta)) 
      eta<-linkinv(eta)
    }
  }
  if (se.fit)
  { H<-list(fit=eta,se.fit=se) }
  else
  { H<-eta}
 
  # tidy up? 
   
  H # ... and return
}

plot.gam<-function(x,rug=TRUE,se=TRUE,pages=0,scale=-1,n=100, ...)

# Create an appropriate plot for each smooth term of a GAM.....
# x is a gam object
# rug determines whether a rug plot should be added to each plot
# se determines whether twice standard error bars are to be added
# pages is the number of pages over which to split output - 0 implies that 
# graphic settings should not be changed for plotting
# scale -1 for same y scale for each plot
#        0 for different y scales for each plot
# n - number of x axis points to use for plotting each term

{ 
  fix.form<-function(formula)
  # internal function to replace "$" on rhs of formula with ".", to allow construction of 
  # data frame for prediction.....  
  { ft<-as.character(formula)
    work.around<-as.character(deparse(formula[3]))
    ft[3]<-""
    for (i in 1:length(work.around)) ft[3]<-paste(ft[3],work.around[i],sep="")
    form<-paste(ft[2],ft[1])
    for (i in 1:(nchar(ft[3])-2))
    { c<-substring(ft[3],i,i)
      if (c=="$") c<-"."
      form<-paste(form,c,sep="")
    }
    as.formula(form)
  }
  
  rename<-function (x)
  # internal function for renaming rows of x$x
  { nm<-rownames(x)
    for (j in 1:length(nm))
    { tmp<-nm[j]
      nm[j]<-""
      for (i in 1:nchar(tmp))
      { c<-substring(tmp,i,i)
        if (c=="$") c<-"."
        nm[j]<-paste(nm[j],c,sep="")
      }
    }
    nm
  }

  x$formula<-fix.form(x$formula)  # replace all "$" characters on rhs of formula
  rownames(x$x)<-rename(x$x)     # replace all "$" characters in row names of x$x
  
  m<-length(x$df) # number of smooth terms
  if (se && x$Vp[1,1]<=0) 
  { se<-FALSE
    warning("No variance estimates available")
  }
  # sort out number of pages and plots per page
  if (pages>m) pages<-m
  if (pages<0) pages<-0
  if (pages!=0)    # figure out how to display things
  { ppp<-m%/%pages
    if (m%%pages!=0) 
    { ppp<-ppp+1
      while (ppp*(pages-1)>=m) pages<-pages-1
      if (m%%pages) last.pages<-0 else last.ppp<-m-ppp*pages
    } 
    else last.ppp<-0
    # now figure out number of rows and columns
    c<-trunc(sqrt(ppp))
	if (c<1) c<-1
    r<-ppp%/%c
    if (r<1) r<-1
    while (r*c<ppp) r<-r+1
    while (r*c-ppp >c && r>1) r<-r-1
    while (r*c-ppp >r && c>1) c<-c-1 
    oldpar<-par(mfrow=c(r,c))
	ylim<-c(r,c)
  } else
  { ppp<-1}
  
  # now create dataframe to send to predict.gam() to return smooth terms
  xx<-data.frame(1:n)
  if (x$nsdf>0) 
  for (i in 1:x$nsdf) xx[i]<-seq(1:n)*0
  for (i in 1:m)
  { x0<-x$xp[i,1]       # min x
    x1<-x$xp[i,x$df[i]] # max x
    dx<-(x1-x0)/(n-1) 
	xx[i+x$nsdf]<-seq(x0,x1,dx)
  }
  names(xx)<-rownames(x$x)
  pl<-predict.gam(x,xx,type="terms",se,plot.call=TRUE)
 
  if (se)   # pl$fit and pl$se.fit
  { if (scale==-1)
    for (i in 1:m)
    { ul<-pl$fit[x$nsdf+i,]+2*pl$se.fit[x$nsdf+i,]
      ll<-pl$fit[x$nsdf+i,]-2*pl$se.fit[x$nsdf+i,]
	  if (i==1) ylim<-c(min(ll),max(ul))
	  else
	  { if (min(ll)<ylim[1]) ylim[1]<-min(ll)
	    if (max(ul)>ylim[2]) ylim[2]<-max(ul)
      }
    } 
    for (i in 1:m)
    { ul<-pl$fit[x$nsdf+i,]+2*pl$se.fit[x$nsdf+i,]
      ll<-pl$fit[x$nsdf+i,]-2*pl$se.fit[x$nsdf+i,]
	  if (scale==0) { ylim<-c(min(ll),max(ul))}
	  title<-paste("s(",rownames(x$x)[i+x$nsdf],",",as.character(round(x$edf[i],2)),")",sep="")
      if (interactive()&&i>1&&(i-1)%%ppp==0) readline("Press return for next page....")
      plot(xx[[x$nsdf+i]],pl$fit[x$nsdf+i,],type="l",xlab=rownames(x$x)[i+x$nsdf],ylim=ylim,ylab=title)
	  lines(xx[[x$nsdf+i]],ul,lty=2)
      lines(xx[[x$nsdf+i]],ll,lty=2)
	  rug(as.numeric(x$x[x$nsdf+i,]))
    }
  } else
  { if (scale==-1)
    for (i in 1:m)
	{ if (i==1) ylim<-range(pl[x$nsdf+i,])
	  else
      { if (min(pl[x$nsdf+i,])<ylim[1]) ylim[1]<-min(pl[x$nsdf+i,])
	    if (max(pl[x$nsdf+i,])>ylim[2]) ylim[2]<-max(pl[x$nsdf+i,])
	  }
	}
    for (i in 1:m)
    { if (interactive() && i>1&&(i-1)%%ppp==0) readline("Press return for next page....")
	  title<-paste("s(",rownames(x$x)[i+x$nsdf],",",as.character(round(x$edf[i],2)),")",sep="")
      if (scale==0) ylim<-range(pl[x$nsdf+i,])
	  plot(xx[[x$nsdf+i]],pl[x$nsdf+i,],type="l",,xlab=rownames(x$x)[i+x$nsdf],ylab=title,ylim=ylim)
      rug(as.numeric(x$x[x$nsdf+i,]))
    } 
  }
  if (pages>0) par(oldpar)
}


.First.lib <- function(lib, pkg) {
    library.dynam("mgcv", pkg, lib)
    cat("This is mgcv 0.5.3\n")
}










