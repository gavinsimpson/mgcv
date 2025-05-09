/* Source code for mgcv.dll/.so multiple smoothing parameter estimation code,
suitable for interfacing to R 

Copyright (C) 2000-2012 Simon N. Wood  simon.wood@r-project.org

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
(www.gnu.org/copyleft/gpl.html)

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
USA. */

#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tprs.h"
#include "mgcv.h"
#include "matrix.h"
#include "qp.h"
#include "general.h"

#define round(a) ((a)-floor(a) <0.5 ? (int)floor(a):(int) floor(a)+1)


//void RPackSarray(int m,matrix *S,double *RS)
/* Packs an array of matrices S[] into an R array RS in the manner described in RUnpackSarray
*/
//{ int start,i,j,k;
//  start=0;
//  for (k=0;k<m;k++)
//  { for (i=0;i<S[k].r;i++) for (j=0;j<S[k].c;j++) RS[start+i+S[k].r*j]=S[k].M[i][j];
//    start += S[k].r*S[k].c;
//  }
//}


matrix getD(matrix h,int nak)

/* the matrix mapping the value of the spline to the gradients at the knots.
   nak is true for 'not-a-knot' end conditions at the early end, otherwise
   'natural' end conditions are used. If there are only 2 knots then the spline
   is taken as a straight line if only 1 a constant. */

{ long i,j,n;
  matrix T,D,Res;
  n=h.r+1;
  T=initmat(n,n);D=initmat(n,n);Res=initmat(n,n);
  for (i=0;i<n;i++) for (j=0;j<n;j++)
  { T.M[i][j]=0.0;D.M[i][j]=0.0;}
  if (n==1L)
  { Res.M[0][0]=0.0;
  } else
  if (n==2L)
  { Res.M[0][0]=Res.M[1][0]=-1.0/h.V[0];
    Res.M[0][1]=Res.M[1][1]=1.0/h.V[0];
  } else
  { for (i=0;i<n;i++) T.M[i][i]=2.0;
    for (i=1;i<n-1;i++)
    { T.M[i][i-1]=h.V[i]/(h.V[i]+h.V[i-1]);
      T.M[i][i+1]=1.0-T.M[i][i-1];
      D.M[i][i-1]= -3.0*T.M[i][i-1]/h.V[i-1];
      D.M[i][i+1]=3.0*T.M[i][i+1]/h.V[i];
      D.M[i][i]= -(D.M[i][i+1]+D.M[i][i-1]);
    }
    if (!nak)
    { T.M[0][1]=1.0;D.M[0][0]= -3.0/h.V[0];D.M[0][1]= -D.M[0][0];}
    else
    { T.M[0][1]=2.0*(h.V[0]+h.V[1])/h.V[1];
      D.M[0][0]= -2.0*(3.0*h.V[0]+2.0*h.V[1])/
		(h.V[0]*(h.V[0]+h.V[1]));
      D.M[0][2]=2.0*h.V[0]*h.V[0]/
      (h.V[1]*h.V[1]*(h.V[0]+h.V[1]));
      D.M[0][1]= -D.M[0][0]-D.M[0][2];
    }
    T.M[n-1][n-2]=1.0;D.M[n-1][n-2]= -3.0/h.V[n-2];
    D.M[n-1][n-1]= -D.M[n-1][n-2];
    invert(&T);
    matmult(Res,T,D,0,0);
  }
  freemat(T);freemat(D);
  return(Res);
}

void MonoCon(matrix *A,matrix *b,matrix *x,int control,double lower,double upper ) 

/* gets matrices A and b for constraints of the form Ay>=b ensuring monotonic
   change  of the cubic spline interpolating (x_i,y_i) where h_i=x_{i+1}-x_i 
   
   control indicates type of constraints:
   up=control/4 - 0 for decrease, 1 for increase
   lo=(control-up*4)/2 - 1 for lower bound, 0 no lower bound
   hi=(control-up*4-lo*2) - 1 for upper bound, 0 no upper bound
   control = 4*up+2*lo+hi
*/

{ long i,j,n;
  int up,lo,hi;
  double m;
  matrix h,D;
  h=initmat(x->r-1,1L);
  n=h.r;
  for (i=0;i<n;i++) h.V[i]=x->V[i+1]-x->V[i];
  D=getD(h,0);
  up=control/4;control=control%4;
  lo=control/2;control=control%2;
  hi=control;
  if (up) m= -1.0; else m=1.0; 
  (*A)=initmat(4*n+hi+lo,n+1);
  for (i=0;i<n;i++)
  { for (j=0;j<n+1;j++)
    { if (j==i)
      { A->M[i][j]=(D.M[i][j]+3.0/h.V[i])*m;   /**not certain of d.M update**/
	    A->M[i+n][j]=(D.M[i+1][j]+3.0/h.V[i])*m;
	    A->M[i+2*n][j]=m;
	    A->M[i+3*n][j]= -D.M[i][j]*m;
      } else
      if (j==(i+1))
      { A->M[i][j]=(D.M[i][j]-3.0/h.V[i])*m;
	    A->M[i+n][j]=(D.M[i+1][j]-3.0/h.V[i])*m;
	    A->M[i+2*n][j]= -m;
	    A->M[i+3*n][j]= -D.M[i][j]*m;
      } else
      { A->M[i][j]=D.M[i][j]*m;
	    A->M[i+n][j]=D.M[i+1][j]*m;
	    A->M[i+2*n][j]=0.0;
	    A->M[i+3*n][j]= -D.M[i][j]*m;
      }
    }
  }
  *b = initmat(A->r,1L);
  if (lo)
  { for (j=0;j<n+1;j++) A->M[4*n][j]=0.0;
    if (up) A->M[4*n][0]=1.0; else A->M[4*n][n]=1.0;
    b->V[4*n]=lower;
  }
  if (hi)
  { for (j=0;j<n+1;j++) A->M[4*n][j]=0.0;
    if (up) A->M[4*n+lo][n]=-1.0; else A->M[4*n+lo][0]=-1.0;
    b->V[4*n+lo] = -upper;
  }
  freemat(D);
  freemat(h);
}


void getFS(double *x,int n,double *S,double *F) {
/* x contains ascending knot sequence for a cubic regression spline
   Routine finds wigglness penalty S and F such that F' maps function 
   values at knots to second derivatives. See Wood 2006 section 4.1.2.
   F and S are n by n. F is F' in 4.1.2 notation.
*/
  double *D,*ldB,*sdB,*h,*Di,*Di1,*Di2,*Fp,*Sp,a,b,c;
  int i,j,n1,n2;
  /* create knot spacing vector h */
  h = (double *)CALLOC((size_t)(n-1),sizeof(double));
  for (i=1;i<n;i++) h[i-1] = x[i]-x[i-1];

  /* create n-2 by n matrix D: D[i,i] = 1/h[i], D[i,i+1] = -1/h[i]-1/h[i+1]
     D[i,i+2] = 1/h[i+1], for i=0..(n-3). D is n-2 by n. */
  D = (double *)CALLOC((size_t)(n*(n-2)),sizeof(double));
  n1 = n-1;n2=n-2;
  for (Di=D,Di1=D+n2,Di2=Di1+n2,i=0;i<n2;i++,Di+=n1,Di1+=n1,Di2+=n1) {
    *Di = 1/h[i];*Di2 = 1/h[i+1];*Di1 = - *Di - *Di2;
  }
  /* create leading diagonal of B*/
  ldB = (double *)CALLOC((size_t)(n2),sizeof(double));
  for (i=0;i<n2;i++) ldB[i] = (h[i]+h[i+1])/3;
  sdB = (double *)CALLOC((size_t)(n2-1),sizeof(double));
  for (i=1;i<n2;i++) sdB[i-1] = h[i]/6;
  /* Now find B^{-1}D using LAPACK routine DPTSV (result in D) */
  F77_CALL(dptsv)(&n2,&n,ldB,sdB,D,&n2,&i);

  /* copy B^{-1}D into appropriate part of F */
  Di=D;
  for (i=0;i<n;i++) {
    Fp = F+i; /* point to row i of F */
    *Fp=0.0;Fp+=n;
    /* col i of D copied to row i of F */
    for (j=0;j<n2;j++,Fp+=n,Di++) *Fp = *Di; 
    *Fp=0.0;
  }

  /* now create D'B^{-1}D efficiently */
  a = 1/h[0];  /* row 0 */
  for (Sp=S,Di=D,i=0;i<n;i++,Sp+=n,Di+=n2) *Sp = *Di * a;
  if (n>3) {
    a = -1/h[0] - 1/h[1];b = 1/h[1]; /* row 1 */
    for (Sp=S+1,Di1=D+1,Di=D,i=0;i<n;i++,Sp+=n,Di+=n2,Di1+=n2) *Sp = *Di * a + *Di1 * b;
    for (j=2;j<n2;j++) { /* rows 2 to n-3 */
      a = 1/h[j-1];c = 1/h[j];b = -a -c; 
      for (Sp=S+j,Di=D+j-2,Di1 = D +j-1,Di2=D + j,i=0;i<n;i++,Sp+=n,Di+=n2,Di1+=n2,Di2+=n2) 
        *Sp = *Di * a + *Di1 * b + *Di2 * c;
    }
    j = n2; /* n-2 */
    a = 1/h[j-1]; b = -1/h[j-1] - 1/h[j]; /* row n-2 */
    for (Sp=S+n2,Di1=D+n2-1,Di=D+n2-2,i=0;i<n;i++,Sp+=n,Di+=n2,Di1+=n2) *Sp = *Di * a + *Di1 * b;
  } else { /* D' has only one column */
    a = -1/h[0] - 1/h[1]; 
    for (Sp=S+1,Di=D,i=0;i<n;i++,Sp+=n,Di+=n2) *Sp = *Di * a;
  }
  j = n2;
  a = 1/h[j]; /* row n-1 */
  for (Sp=S+n1,Di=D+n2-1,i=0;i<n;i++,Sp+=n,Di+=n2) *Sp = *Di * a;

  FREE(ldB);FREE(sdB);FREE(h);FREE(D);
} /* end of getFS*/


void crspl(double *x,int *n,double *xk, int *nk,double *X,double *S, double *F,int *Fsupplied) {
/* Routine to compute model matrix and optionally penalty matrix for cubic regression spline.
   * nk knots are supplied in an increasing sequence in xk. 
   * n data are in x (arbitrary order).
   * If Fsupplied!=0 then F' is matrix mapping function values at knots to second derivs,
     otherwise F and the penalty matrix S are computed and returned, along with X.         
*/
  int i,j=0,k,extrapolate,jup,jmid;
  double xlast=0.0,h=0.0,xi,kmax,kmin,ajm,ajp,cjm,cjp,*Fp,*Fp1,*Xp,xj,xj1,xik;
  if (! *Fsupplied) getFS(xk,*nk,S,F);
  kmax = xk[*nk-1];kmin = xk[0];
  for (i=0;i<*n;i++) { /* loop through x */
    xi = x[i];extrapolate=0;
    /* find interval containing x[i] */
    if (xi < kmin||xi>kmax) {
      extrapolate=1;
    } else if (i>0 && fabs(xlast-xi) < 2*h) { /* use simple direct search */
      while (xi <= xk[j] && j > 0) j--;
      while (xi > xk[j+1] && j < *nk-2) j++;
      /* next line should not be needed, except under dodgy use of 
         fpu registers during optimization... */
      if (j<0) j=0;
      if (j > *nk-2) j = *nk - 2; 
      /* now xk[j] <= x[i] <= xk[j+1] */ 
    } else { /* bisection search required */ 
      j=0;jup=*nk-1;
      while (jup-j>1) {
        jmid = (jup+j) >> 1; /* a midpoint */
        if (xi > xk[jmid]) j = jmid; else jup = jmid;
      }
      /* now xk[j] <= x[i] <= xk[j+1] */ 
    } /* end of bisection */
    
    /* knot interval containing x[i] now known. Compute spline basis */ 
    if (extrapolate) { /* x[i] is outside knot range */
      if (xi<kmin) {
        j = 0;
        h = xk[1] - kmin;
        xik = xi - kmin;
        cjm = -xik*h/3;
        cjp = -xik*h/6;
        Xp = X + i; /* ith row of X */
        for (Fp = F,Fp1 = F + *nk,k=0;k < *nk;k++,Xp += *n,Fp++,Fp1++) *Xp = cjm * *Fp + cjp * *Fp1 ;
        X[i] += 1 - xik/h;
        X[i + *n] += xik/h;
      } else { /* xi>kmax */
        j = *nk-1;
        h = kmax - xk[j-1];
        xik = xi - kmax;
        cjm= xik*h/6;
        cjp = xik*h/3;
        Xp = X + i; /* ith row of X */
        for (Fp1 = F+ j * *nk,Fp = Fp1 - *nk,k=0;k < *nk;k++,Xp += *n,Fp++) 
             *Xp = cjm * *Fp + cjp * *Fp1  ;
        X[i + *n * (*nk-2)] += - xik/h;
        X[i + *n * (*nk-1)] += 1+ xik/h;
      }
    } else { /* routine evaluation */
      xj = xk[j];xj1=xk[j+1];
      h = xj1-xj; /* interval width */
      ajm = (xj1 - xi);ajp = (xi-xj);
      cjm = (ajm*(ajm*ajm/h - h))/6;
      cjp = (ajp*(ajp*ajp/h - h))/6;
      ajm /= h;ajp /= h;
      
      Xp = X + i; /* ith row of X */

      for (Fp = F+ j * *nk, Fp1 = F+(j+1)* *nk,k=0;k < *nk;k++,Xp += *n,Fp++,Fp1++) 
        *Xp = cjm * *Fp + cjp * *Fp1;

      Xp = X + i + j * *n;
      *Xp += ajm; Xp += *n; *Xp += ajp;
    } 
    /* basis computation complete */
    xlast=xi;
  }

} /* end crspl */



void MinimumSeparation(double *x,int *n, int *d,double *t,int *m,double *dist) {
/* For each of n ppoints point x[i,] calculates the minimum Euclidian distance 
   to a point in m by d matrix t. These distances are stored in dist. 
*/
  int one=1,*ni;
  kdtree_type kd;
  kd_tree(t,m,d,&kd); /* build kd tree for target points */
  ni = (int *)CALLOC((size_t)*n,sizeof(int));
  k_newn_work(x,kd,t,dist,ni,n,m,d,&one);
  // for (i=0;i<*n;i++) {
  //  k = closest(&kd,t,x + i * *d,*m,&j,-1); /* index of nearest neighbour of x[i,] */
  //  dist[i] = xidist(x + i * *d,t,k,*d, *m); /* distance to this nearest neighbour */
  //}
  FREE(ni);
  free_kdtree(kd);
}

void MinimumSeparation_old(double *gx,double *gy,int *gn,double *dx,double *dy, int *dn,double *dist)
/* For each point gx[i],gy[i] calculates the minimum  Euclidian distance to a point in dx[], dy[].
   These distances are stored in dist. 
   Aweful routine:  O(gn * dn) cost!! 
*/

{ 
  double sep,xx,yy,*dum,*xdum,*ydum;
  int n,m;
  n = *gn;m = *dn;
  for (dum=dist;dum < dist + n; dum++,gx++,gy++)
  { xx= *gx - *dx;yy = *gy - *dy;*dum = xx*xx + yy*yy; /* first separation */
    for (xdum=dx+1,ydum=dy+1;xdum < dx + m;xdum++,ydum++)
    { xx= *gx - *xdum;yy = *gy - *ydum;sep = xx*xx + yy*yy; /* subsequent separations */
      if (sep < *dum) *dum = sep;
    }
    *dum = sqrt(*dum);
  }
} 



void RuniqueCombs(double *X,int *ind,int *r, int *c)

/* X is a matrix. This routine finds its unique rows and strips out the 
   duplicates. This is useful for finding out the number of unique covariate
   combinations present in a set of data. */

{ matrix B,Xd;
  int i,*ind1;
  B=Rmatrix(X,(long)(*r),(long)(*c));
  Xd=initmat(B.r,B.c+1);
  Xd.c--;mcopy(&B,&Xd);freemat(B);Xd.c++;
  for (i=0;i<Xd.r;i++) Xd.M[i][Xd.c-1]=(double)i;
  ind1=Xd_strip(&Xd);
  for (i=0;i<*r;i++) ind[i] = ind1[i]; /* copy index for return */
  Xd.c--; /* hide index array  */
  RArrayFromMatrix(X,Xd.r,&Xd);  /* NOTE: not sure about rows here!!!! */
  *r = (int)Xd.r; 
  freemat(Xd);FREE(ind1);

}

void RMonoCon(double *Ad,double *bd,double *xd,int *control,double *lower,double *upper,int *n)
/* obtains coefficient matrices for imposing monotonicity (and optionally bounds) on a 
   cubic regression spline with n knots located as specified in xd. 
   
   control indicates type of constraints:
   up=control/4 - 0 for decrease, 1 for increase
   lo=(control-up*4)/2 - 1 for lower bound, 0 no lower bound
   hi=(control-up*4-lo*2) - 1 for upper bound, 0 no upper bound
   control = 4*up+2*lo+hi
   
   lower and upper are the bounds to impose (ignored if control doesn't
   indicate that they should be used).

   Ad will have 4(n-1)+lo+hi rows and n columns
   bd will have 4(n-1)+lo+hi rows

*/ 

{ int i;
  matrix x,A,b;
  x=initmat((long)*n,1L);
  for (i=0;i<x.r;i++) x.V[i]=xd[i];
  MonoCon(&A,&b,&x,*control,*lower,*upper); 
  RArrayFromMatrix(Ad,A.r,&A);
  RArrayFromMatrix(bd,b.r,&b);
 
  freemat(x);freemat(A);freemat(b);  

}

