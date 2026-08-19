#include <R.h> 
#include <Rmath.h>
#include <Rinternals.h> // RK addition
#include <R_ext/RS.h>  // RK addition
#include <R_ext/Lapack.h> // RK addition
#include <R_ext/BLAS.h> // RK addition
#include <R_ext/Applic.h>  // ST addition
#include <R_ext/Rdynload.h>  // ST addition
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>




static int *checklist;
static double *tmplike;
static double *effects;
static double *residuals;
static int *pivot;
static double *work;
static double *qraux;
static double *coef;
static double *qr;
static double *Sumstats;




//To register C functions to be read within R!!!
//void R_init_EnvCpt(DllInfo *info){
//    R_registerRoutines(info, NULL, NULL, NULL, NULL);
//    R_useDynamicSymbols(info, TRUE);
//}








//Calculate summary statistics for the regression quadratic cost function.
void RegARQuadCost_SS(double *X, int *n, int *nc, double *SS, int *m){
  //X - data as double vector of length n*p   (passed from R, n*2(P+1))
  //n - number or rows of X
  //nc - number of columns of X (passed from R, 2(P+1))
  //SS - summary statistics vector of length (n+1)*m
  //m - number of rows of SS, m = p*(p+1)/2  (passed from R, (2P+2)*(2P+3)/2 = (P+1)(2p+3)
  
  
  int pos;       //position to fill in SS
  int i,j,l;     //loop indicies
  pos = 0;
  //Set first column of SS to zero
  for(i = 0; i < *m; i++){
    SS[pos] = 0;
    pos++;
  }
  
  //For each new row in X, calculate the unique cross prod of X (lower tri) and
  //   add to previous summary statistic.
  for(i = 0; i < *n; i++){
    for(j = 0; j < *nc; j++){ // nc is number of extended columns
      for(l = j; l < *nc; l++){
        SS[pos] = SS[pos - *m] + X[j * *n + i]*X[l * *n + i];
        // previous column of SS + X[i,j]*X[i,l]
        pos++;
      }
    }
  }
  
  
  return;
}




//Evaluate the regression quadratic cost function based on summary statistics
void RegARQuadCostFunc(double *SS, double *phi, int *m, int *n, int *P, int *start, int *end, 
                       double *cost, double *tol, int *error, double *scale, int *MBIC){
  //SS    - Summary statistics
  //phi   - fixed value of phi
  //m     - number of rows of SS
  //n     - number of columns of SS
  //P     - number of extended regressors, nb  m = (p+1)*(2p+3)    // previously m = (p+1)*(p+2)/2
  //start - index at start of segment
  //end   - index at end of segment
  //tol   - tolerence for lm
  //error - error index
  //scale - if>0, cost=-2loglike @ fixed scale/variance
  //        if=0, cost=-2logLik()
  //        if<0, cost=RSS
  //MBIC  - 1 if MBIC penalty used, 0 if not.
  
  int p;
  p= *P/2 -1; // number of regressors -1 as remove the data
  
  //-- memory allocation --
  //Summary statistics for segment
  double *Sumstats = calloc(*m, sizeof(double)); 
  if(Sumstats == NULL){
    *error = 101;
    goto err101;
  }
  //Lagged by 1 Summary statistics for segment
  double *Sumstatslag = calloc(*m, sizeof(double)); 
  if(Sumstatslag == NULL){
    *error = 201;
    goto err201;
  }
  //t(X) %*% X matrix (populated by Sumstats)
  double *XX = calloc(p * p, sizeof(double));  
  if(XX == NULL){
    *error = 102;
    goto err102;
  }
  double *XXlag = calloc(p * p, sizeof(double));  
  if(XXlag == NULL){
    *error = 202;
    goto err202;
  }
  double *XXcross = calloc(p * p, sizeof(double));  
  if(XXcross == NULL){
    *error = 302;
    goto err302;
  }
  //t(X) %*% y vector (populated by Sumstats)
  double *Xy = calloc(p, sizeof(double));
  if(Xy == NULL){
    *error = 103;
    goto err103;
  }
  double *Xylag = calloc(p, sizeof(double));
  if(Xylag == NULL){
    *error = 203;
    goto err203;
  }
  double *XycrossX = calloc(p, sizeof(double));
  if(XycrossX == NULL){
    *error = 303;
    goto err303;
  }
  double *Xycrossy = calloc(p, sizeof(double));
  if(Xycrossy == NULL){
    *error = 403;
    goto err403;
  }
  //coefficients reordered according to pivot (see below)
  double *beta = calloc(p, sizeof(double));
  if(beta == NULL){
    *error = 104;
    goto err104;
  }
  //copy of XX (converts to qr on return of dqrls)
  double *qr = calloc(p * p, sizeof(double));
  if(qr == NULL){
    *error = 105;
    goto err105;
  }
  //returns coefficient estimates of lm
  double *coef = calloc(p, sizeof(double));
  if(coef == NULL){
    *error = 106;
    goto err106;
  }
  //returns auxiliary information on qr decomposition
  double *qraux = calloc(p, sizeof(double));
  if(qraux == NULL){
    *error = 107;
    goto err107;
  }
  //working space
  double *work = calloc(2 * p, sizeof(double));
  if(work == NULL){
    *error = 108;
    goto err108;
  }
  //1:p (changes if cols of XX are permuted)
  int *pivot = calloc(p, sizeof(int));
  if(pivot == NULL){
    *error = 109;
    goto err109;
  }
  //copy of Xy (converts to y-Xb)
  double *residuals = calloc(p, sizeof(double));
  if(residuals == NULL){
    *error = 110;
    goto err110;
  }
  //copy of Xy (converts to orthog effects: t(q) %*% Xy)
  double *effects = calloc(p, sizeof(double));
  if(effects == NULL){
    *error = 111;
    goto err111;
  }
  //copy of Xy
  double *Xytmp = calloc(p, sizeof(double));
  if(Xytmp == NULL){
    *error = 112;
    goto err112;
  }
  
  
  int i, j, tri, skip;           //indices
  int rank = 0;            //returns estimated rank of system (XX)
  int ny = 1;              //number of columns of Xy (here, always 1)
  int nn = *end - *start;  //number of observations
  double RSS;              //residual sum of squares
  
  
  //Evaluate summary statistics for segment
  for(i = 0; i < *m; i++){
    //[yy, Xy, XX(lower.tri)]
    Sumstats[i] = SS[*end * *m + i] - SS[*start * *m + i]; // not end+1 as the column of zeroes is at index zero!!!
    Sumstatslag[i] = SS[(*end-1) * *m + i] - SS[(*start-1) * *m + i]; // just lagged by 1 on start and end
  } // includes the extended entries
  
  //Populate XX and Xy with values from Sumstats.
  //Also, make copies and initialise values for dqrls
  tri = 0;
  skip=p+1;               // skip the lagged at the start
  //  int endorg=(3*p*p +6*p+6)/2; // the end of the original set of p (note that the lagged ones are intersperced so this isn't (p+1)(p+2)/2)
  for(i = 0; i < p; i++){
    Xy[i] = Sumstats[i+1];
    Xytmp[i] = Xy[i];           //copy of Xy
    Xylag[i] = Sumstatslag[i+1];
    XycrossX[i] = Sumstats[p+1+ i+1];
    Xycrossy[i]=Sumstats[(3*p+2)+tri+skip-p-1]; // the i is in the tri
    coef[i] = 0.0;           //initialise to zero
    residuals[i] = Xy[i];    //copy of Xy
    effects[i] = Xy[i];      //copy of Xy
    pivot[i] = i+1;          //index from 1:p
    qraux[i] = 0.0;          //initialise to zero
    work[i] = 0.0;           //initialise first half to zero
    work[i + p] = 0.0;      //initialise second half to zero
    for(j = i; j < p; j++){
      // need to skip the extra entries for the extended
      XX[i * p + j] = Sumstats[p + 1 + tri + j +skip];
      XXlag[i * p + j] = Sumstatslag[p + 1 + tri + j +skip];
      qr[i * p + j] = XX[i * p + j];   //copy of XX
    }
    for(j=0;j<p;j++){
      XXcross[i * p + j] = Sumstats[p + 1 + tri + j +skip+p+1]; // not symmetric
    }
    skip+=p+1; // add another p+1 to skip the extended part
    tri += p - i - 1;
    for(j = 0; j < i; j++){
      XX[i * p + j] = XX[j * p + i];   //Copy to make XX symmetric
      XXlag[i * p + j] = XXlag[j * p + i];   //Copy to make XXlag symmetric
      qr[i * p + j] = XX[i * p + j];   //Copy of XX
    }
  }
  
  //Call to routine used by lm to solve system of linear equations
  F77_CALL(dqrls)(qr, &p, &p, Xytmp, &ny, tol, coef, residuals, effects, &rank ,
           pivot, qraux, work);
  
  //Extract coefficients in order
  for(i = 0; i < p; i++){
    beta[pivot[i]-1] = coef[i];
  }  
  
  
  //Calculate the quadratic cost
  RSS = Sumstats[0] + *phi* *phi*Sumstatslag[0] -2* *phi*Sumstats[p+1]; // standard + lagged + cross
  for(i = 0; i < p; i++){
    RSS += -2 * beta[i] * Xy[i] - 2* *phi * *phi*beta[i]*Xylag[i] + 2* *phi*beta[i]*XycrossX[i] + 2* *phi*beta[i]*Xycrossy[i]; // standard + lagged + cross
    if(*start==960){
      if(*end==1021){
        //Rprintf("%d, %d, %d,%d: %f \n", *start, *end, i,j, beta[i]);
      }
    }
    for(j = 0; j < p; j++){
      RSS += beta[i] * beta[j] * XX[i * p + j] + *phi* *phi*beta[i]*beta[j]*XXlag[i*p + j] - 2* *phi*beta[i]*beta[j]*XXcross[i*p + j]; // standard + lagged + cross
      if(*start==960){
        if(*end==1021){
          // Rprintf("%d, %d, %d,%d: %f, %f, %f \n", *start, *end, i,j, XX[i * p + j], XXlag[i*p + j], XXcross[i*p + j]);
        }
      }
    }
  }
  
  
  
  if(*scale == 0){
    // if(RSS<=0){ // RK to work out why RSS is negative!!!!! See testing file for example.
    //   RSS=0.00000000000000001; // 1e-16
    // }
    *cost = nn + nn*log(2 * M_PI * RSS) - nn*log(nn); //-2*logLik()
  }else if(*scale > 0){
    *cost = nn*log(2 * M_PI * *scale) + (RSS / *scale); //-2LL(sig2=scale)
  }else{
    *cost = RSS; //RSS
  }
  if(*MBIC==1){
    //*cost = RSS+log(nn); // extra log(length segment) for MBIC cost
  }
  // if(*start==960){
  //   Rprintf("end: %d, RSS: %f, cost: %f \n",*end,RSS,*cost);
  // }
  // //Free allocated memory
  free(Xytmp);
  err112: free(effects);
  err111: free(residuals);
  err110: free(pivot);
  err109: free(work);
  err108: free(qraux);
  err107: free(coef);
  err106: free(qr);
  err105: free(beta);
  err104: free(Xycrossy);
  err403: free(XycrossX);
  err303: free(Xylag);
  err203: free(Xy);
  err103: free(XXcross);
  err302: free(XXlag);
  err202: free(XX);
  err102: free(Sumstatslag);
  err201: free(Sumstats);
  err101: return;
  
  
}




//Find the minimum case
static void min_which(double *data, int *n, double *minval, int *minid){
  //data - values for which to find the minimum
  //n - number of items to search
  //minval - minimum value
  //minid - index from start where to find minimum
  int i;
  *minid = 0;
  *minval = data[*minid];
  for(i = 1; i < *n; i++){
    if(data[i] < *minval){
      *minid = i;
      *minval = data[i];
    }
  }
  return;
}








//Determine cpts for Normal regression under method PELT
void CptRegAR_Normal_PELT(double *data, double *phi, int *n, int *m, double *pen, int *cptsout,
                          int *error, double *shape, int *minseglen, double *tol, double *lastchangelike,
                          int *lastchangecpts, int *numchangecpts, int *MBIC, int *ncheck){
  
  
  //data           - vectorised matrix of size n x m
  //phi            - fixed value of phi
  //n              - number of records
  //m              - number of data points per record
  //pen            - penalty value
  //cptsout        - position of estimated changepoints
  //error          - Index stating if an error has occured (error=0 if ok)
  //shape          - known variance (if=0,use MLE, if<0, cost=RSS)
  //minseglen      - minimum segment length
  //lastchangelike - working space (cost value at last changepoint)
  //lastchangecpts - working space (index of last changepoint)
  //numchangecpts  - working space (number of cpts that have occured so far)
  //MBIC           - 1 if MBIC penalty, 0 if not
  
  
  int p = *m/2 - 1;   //number of regressors
  int np1 = *n + 1; //ncols of summary statistcs array 
  int size = (*m * (*m + 1)) * 0.5; //nrows of summary statistics array
  int nchecklist, nchecktmp;     //number of items in the checklist
  double minout;    //minimum cost value
  int tstar, i, j, start;  //indicies
  int whichout;     //index corresponding to the minimum cost value
  double segcost;   //cost over specified segment
  *error = 0;
  
  
  //working space: position of last changepoint for all 'active' branches
  int *checklist = (int *)calloc(np1, sizeof(int));
  if(checklist==NULL){
    *error = 1;
    goto err1;
  }
  //working space: updated cost wrt latest observation for all 'active' branches
  double *tmplike = (double *)calloc(np1, sizeof(double));
  if(tmplike == NULL){
    *error = 2;
    goto err2;
  }
  //Summary statistcs
  double *Sumstats = (double *)calloc(np1 * size, sizeof(double));
  if(Sumstats == NULL){
    *error = 3;
    goto err3;
  }
  
  //Evaluate the summary statistics
  RegARQuadCost_SS(data, n, m, Sumstats, &size);
  
  //Initialise
  lastchangelike[0] = -*pen;
  lastchangecpts[0] = 0;
  numchangecpts[0] = 0;
  for(j = 1; j <= *minseglen; j++){
    lastchangelike[j] = 0;
    lastchangecpts[j] = 0;
    numchangecpts[j] = 0;
  }
  //Evaluate cost for second minseglen
  start = 0;
  for(j = (*minseglen+1); j <= (2 * *minseglen); j++){
    RegARQuadCostFunc(Sumstats, phi, &size, &np1, m, &start, &j, lastchangelike + j,
                      tol, error, shape, MBIC);
    if(*error != 0){
      goto err4;
    }
    lastchangecpts[j] = 0;
    numchangecpts[j] = 1;
  }
  //setup check list;
  nchecklist = 2;
  checklist[0] = 0;
  checklist[1] = *minseglen+1;
  
  
  //progress through time series
  for(tstar = (2 * *minseglen)+1; tstar < np1; tstar++){
    R_CheckUserInterrupt(); //Has interrupred the R session? quits if true.
    
    
    for(i = 0; i < nchecklist; i++){
      //Evaluate cost&penalty based on 
      // total cost&pen at last cpt + cost over current segment + penalty
      start = checklist[i];  //last point of last segment
      RegARQuadCostFunc(Sumstats, phi, &size, &np1, m, &start, &tstar, &segcost,
                        tol, error, shape, MBIC);
      if(*error != 0){
        goto err4;
      }
      tmplike[i] = lastchangelike[start] + segcost + *pen;
      // if(checklist[i]==2080){
      //   Rprintf("tstar: %d, lclike: %f, segcost: %f, like: %f \n",tstar,lastchangelike[start],segcost,tmplike[i]);
      // }
    }
    // if(tstar>2141){
    //   if(tstar<2205){
    //     Rprintf("tstar: %d \n",tstar);
    //     for(i=0;i<nchecklist;i++){
    //       Rprintf(", %d",checklist[i]);
    //       // Rprintf("chk: %d, like: %f \n",checklist[i], tmplike[i]);
    //     }
    //     Rprintf("\n");
    //   }
    // }
    //Find and store branch with minimum cost&pen
    min_which(tmplike, &nchecklist, &minout, &whichout);
    lastchangelike[tstar] = minout;
    lastchangecpts[tstar] = checklist[whichout];
    numchangecpts[tstar] = numchangecpts[lastchangecpts[tstar]] + 1;
    
    //Prune out non-minimum branches
    nchecktmp = 0;
    for(i = 0; i < nchecklist; i++){
      // if((tstar-checklist[i])< (2* *minseglen)){
      //   checklist[nchecktmp] = checklist[i];
      //   nchecktmp++;
      //   continue;
      // }
      if(tmplike[i] <= (lastchangelike[tstar] +*pen)){
        checklist[nchecktmp] = checklist[i];
        nchecktmp++;
      }
    }
    nchecklist = nchecktmp;
    
    //Add new cpt to checklist
    checklist[nchecklist] = tstar - *minseglen;
    nchecklist++;
    ncheck[tstar]=nchecklist;
  }
  
  
  //Extract optimal changepoint set
  int ncpts = 0;
  int last = *n;
  while(last != 0){
    cptsout[ncpts] = last;
    last = lastchangecpts[last];
    ncpts++;
  }
  
  
  //Free allocated memory
  err4: free(Sumstats);
  err3: free(tmplike);
  err2: free(checklist);
  err1: return;
}








//Free allocated memory in case R session has been interrupred
void Free_CptRegAR_Normal_PELT(int *error){
  // Error code from CptRegAR_Normal_PELT C function, non-zero => error 
  if(*error==0){
    free((void *)checklist);
    free((void *)tmplike);
    free((void *)Sumstats);
  }
}








///////////////////////////////////////
///////// NEW VERSION /////////////////
///////////////////////////////////////




void RegARCostFunc(double *y, double *design, double *phi, int *n, int *p, int *start, 
                   int *end, double *like, double *tol, int *error, double *scale, int *MBIC,
                   double *coef, double *residuals, double *effects, int *rank, 
                   int *pivot, double *qraux, double *work){
  // Function to calculate the likelihood for the segment of data from start to end
  
  //y         - data
  //design    - vectorised design matrix
  //phi       - fixed value of phi
  //n         - number of rows of design
  //p         - number of columns of design
  //start     - index at start of segment (include in likelihood) - i.e. cpt +1
  //end       - index at end of segment (include in likelihood)
  //tol       - tolerence for lm
  //error     - error index
  //scale     - if>0, cost=-2loglike @ fixed scale/variance
  //            if=0, cost=-2logLik()
  //            if<0, cost=RSS
  //MBIC      - 1 if MBIC penalty used, 0 if not.
  // THE FOLLOWING VARIABLES ARE PASSED IN TO AVOID REPEATED MEMORY ALLOCATION
  //coef      - blank work vector for the coefficients 
  //residuals - blank work vector for the residuals
  //effects   - blank work vector for the effects 
  //rank      - blank work value for the rank
  //pivot     - blank work vector for the pviots
  //qraux     - blank work vector
  //work      - blank work matrix
  
  int nn=*end-*start;
  int i,j;
  
  int ny=1;
  
  //double tmplike[*n];
  double *qr;
  qr = (double *)calloc(*p*nn,sizeof(double));
  if (qr==NULL)   {
    *error = -1;
    goto errlike;
  }
  for(i=0;i<nn;i++){ // dqrls destroys the first argument so need a copy of design
    for(j=0;j<*p;j++){
      *(qr+i+(j*nn))=*(design+*start+i+(j * *n));
      // Rprintf("%d, %d, %f \n", i,j,*(qr+i+(j*nn)));
    }
  } // fills the nnxp matrix
  
  
  // First fit betas and get residuals from a standard OLS fit
  //Call to routine used by lm to solve system of linear equations
  //Rprintf("%d, %d, %d, %f, %d \n", nn,*p, ny,*tol, *rank);
  //Rprintf("%f, %f, %f, %f, %d, %f, %f",*(y+*start),*coef, *residuals, *effects,
  //        *pivot, *qraux, *work);
  F77_CALL(dqrls)(qr, &nn, p, y+*start, &ny, tol,
           coef, residuals, effects, rank, pivot, qraux, work);
  
  // Rprintf("Fit: %f, %f, %f, %f, %f, %f",*(y+*start),*coef, *(coef+1), *residuals, *effects,
  // *qraux);
  
  // use the residuals and filter by phi
  double sse=*residuals * *residuals; // starting with first value where resid[-1]=0
  for(i=1;i<nn;i++){ // calculate sum of squared residuals
    sse+=(*(residuals+i) - *phi * *(residuals+i-1)) * (*(residuals+i) - *phi * *(residuals+i-1));
  }
  // Rprintf("SSE: %f",sse);
  
  if(*scale == 0){
    *like = nn + nn*log(2 * M_PI * sse) - nn*log(nn); //-2*logLik()
  }else if(*scale > 0){
    *like = nn*log(2 * M_PI * *scale) + (sse / *scale); //-2LL(sig2=scale)
  }else{
    *like = sse; //RSS
  }
  if(*MBIC==1){
    *like = *like+log(nn); // extra log(length segment) for MBIC cost
  }
  
  
  free(qr);
  errlike: return;
}








void CptRegAR_PELT(double *y, double *design, double *phi, int *n, int *p, double *pen, int *cptsout,
                   int *error, double *shape, int *minseglen, double *tol, double *lastchangelike,
                   int *lastchangecpts, int *numchangecpts, int *MBIC, int *ncheck){
  
  //y              - vector length n
  //design         - vectorised matrix of size n x p
  //phi            - fixed value of phi
  //n              - number of records
  //p              - number of regressors per element of y
  //pen            - penalty value
  //cptsout        - position of estimated changepoints
  //error          - Index stating if an error has occured (error=0 if ok)
  //shape          - known variance (if=0,use MLE, if<0, cost=RSS)
  //minseglen      - minimum segment length
  //lastchangelike - working space (cost value at last changepoint)
  //lastchangecpts - working space (index of last changepoint)
  //numchangecpts  - working space (number of cpts that have occured so far)
  //MBIC           - 1 if MBIC penalty, 0 if not
  
  int np1 = *n + 1; //ncols of summary statistcs array 
  int nchecklist, nchecktmp;     //number of items in the checklist
  double minout;    //minimum cost value
  int tstar, i, j, start;  //indicies
  int whichout;     //index corresponding to the minimum cost value
  double segcost;   //cost over specified segment
  *error = 0;
  int rank =0;
  
  //working space: position of last changepoint for all 'active' branches
  int *checklist = (int *)calloc(np1, sizeof(int));
  if(checklist==NULL){
    *error = 1;
    goto err1;
  }
  //working space: updated cost wrt latest observation for all 'active' branches
  double *tmplike = (double *)calloc(np1, sizeof(double));
  if(tmplike == NULL){
    *error = 2;
    goto err2;
  }
  //returns coefficient estimates of lm
  double *coef = calloc(*p, sizeof(double));
  if(coef == NULL){
    *error = 3;
    goto err3;
  }
  //returns auxiliary information on qr decomposition
  double *qraux = calloc(*p, sizeof(double));
  if(qraux == NULL){
    *error = 4;
    goto err4;
  }
  //working space
  double *work = calloc(2 * *p, sizeof(double));
  if(work == NULL){
    *error = 5;
    goto err5;
  }
  //1:p (changes if cols of XX are permuted)
  int *pivot = calloc(*p, sizeof(int));
  if(pivot == NULL){
    *error = 6;
    goto err6;
  }
  //copy of Xy (converts to y-Xb)
  double *residuals = calloc(*n, sizeof(double));
  if(residuals == NULL){
    *error = 7;
    goto err7;
  }
  //copy of Xy (converts to orthog effects: t(q) %*% Xy)
  double *effects = calloc(*n, sizeof(double));
  if(effects == NULL){
    *error = 8;
    goto err8;
  }
  
  //Initialise
  lastchangelike[0] = -*pen;
  lastchangecpts[0] = 0;
  numchangecpts[0] = 0;
  for(j = 1; j <= *minseglen; j++){
    lastchangelike[j] = 0;
    lastchangecpts[j] = 0;
    numchangecpts[j] = 0;
  }
  //Evaluate cost for second minseglen
  start = 0;
  for(j = *minseglen; j < (2 * *minseglen); j++){
    RegARCostFunc(y, design, phi, n, p, &start, &j, lastchangelike + j,
                  tol, error, shape, MBIC,
                  coef, residuals, effects, &rank, pivot, qraux, work);
    if(*error != 0){
      goto err9;
    }
    lastchangecpts[j] = 0;
    numchangecpts[j] = 1;
  }
  //setup check list;
  nchecklist = 2;
  checklist[0] = 0;
  checklist[1] = *minseglen;
  
  
  //progress through time series
  for(tstar = 2 * *minseglen; tstar < np1; tstar++){
    R_CheckUserInterrupt(); //Has interrupted the R session? quits if true.
    
    
    
    
    for(i = 0; i < nchecklist; i++){
      //Evaluate cost&penalty based on
      // total cost&pen at last cpt + cost over current segment + penalty
      start = checklist[i];  //last observation of the previous segment but correct
      // index for y and design (as index starts at 0 but tstar
      //starts at 1)
      RegARCostFunc(y, design, phi, n, p, &start, &tstar, &segcost,
                    tol, error, shape, MBIC,
                    coef, residuals, effects, &rank, pivot, qraux, work);
      if(*error != 0){
        goto err9;
      }
      tmplike[i] = lastchangelike[start] + segcost + *pen;
    }
    //Find and store candidate with minimum cost&pen
    min_which(tmplike, &nchecklist, &minout, &whichout);
    lastchangelike[tstar] = minout;
    lastchangecpts[tstar] = checklist[whichout];
    numchangecpts[tstar] = numchangecpts[lastchangecpts[tstar]] + 1;
    
    
    
    
    //Prune out non-minimum changes
    nchecktmp = 0;
    for(i = 0; i < nchecklist; i++){
      // if((tstar-checklist[i])< (2* *minseglen)){
      //   checklist[nchecktmp] = checklist[i];
      //   nchecktmp++;
      //   continue;
      // }
      if(tmplike[i] <= (lastchangelike[tstar] +*pen)){
        checklist[nchecktmp] = checklist[i];
        nchecktmp++;
      }
    }
    nchecklist = nchecktmp;
    
    
    
    
    //Add new cpt to checklist
    checklist[nchecklist] = tstar - *minseglen;
    nchecklist++;
    ncheck[tstar]=nchecklist;
  }
  
  // start = checklist[0];
  // tstar=5;
  // RegARCostFunc(y, design, phi, n, p, &start, &tstar, &segcost,
  //                                 tol, error, shape, MBIC,
  //                                 coef, residuals, effects, &rank, pivot, qraux, work);
  // Rprintf("Null like: %f",segcost);  
  
  //Extract optimal changepoint set
  int ncpts = 0;
  int last = *n;
  while(last != 0){
    cptsout[ncpts] = last;
    last = lastchangecpts[last];
    ncpts++;
  }
  
  //Free allocated memory
  err9: free(effects);
  err8: free(residuals);
  err7: free(pivot);
  err6: free(work);
  err5: free(qraux);
  err4: free(coef);
  err3: free(tmplike);
  err2: free(checklist);
  err1: return;
}




//Free allocated memory in case R session has been interrupted
void Free_CptRegAR_PELT_new(int *error){
  // Error code from CptRegAR_Normal_PELT C function, non-zero => error 
  if(*error==0){
    free((void *)qr);
    free((void *)effects);
    free((void *)residuals);
    free((void *)pivot);
    free((void *)work);
    free((void *)qraux);
    free((void *)coef);
    free((void *)checklist);
    free((void *)tmplike);
  }
}
































//Determine cpts for Normal regression under method AMOC
void CptRegAR_Normal_AMOC(double *data, double *phi, int *n, int *m, double *pen,
                          int *error, double *shape, int *minseglen, double *tol, int *tau,
                          double *nulllike, double *taulike, double *tmplike, int *MBIC){
  
  
  //data           - vectorised matrix of size n x m
  //phi            - fixed value of phi
  //n              - number of records
  //m              - number of data points per record
  //pen            - penalty value
  //error          - Index stating if an error has occured (error=0 if ok)
  //shape          - known variance (if=0,use MLE, if<0, cost=RSS)
  //minseglen      - minimum segment length
  //tau            - estimated single changepoint position
  //nulllike       - estimated cost of no changepoints
  //taulike        - estimated cost of single changepoint
  //tmplike        - working memory: store all single changepoint costs
  //MBIC           - 1 if MBIC penalty, 0 if not
  
  
  int p = *m/2 - 1;   //number of regressors
  int np1 = *n + 1; //ncols of summary statistcs array 
  int size = (*m * (*m + 1)) * 0.5; //nrows of summary statistics array
  int tstar;
  double seg1cost, seg2cost;
  int zero;
  int neval;
  *error = 0;
  
  
  //Summary statistcs
  double *Sumstats = (double *)calloc(np1 * size, sizeof(double));
  if(Sumstats == NULL){
    *error = 1;
    goto err1;
  }
  //Evaluate the summary statistics
  RegARQuadCost_SS(data, n, m, Sumstats, &size);
  
  
  //cost with no changepoints
  zero = 0;
  RegARQuadCostFunc(Sumstats, phi, &size, &np1, m, &zero, n, nulllike, tol, error, shape, MBIC);
  if(*error != 0){
    goto err2;
  }
  
  
  //cost with at most one changepoint
  neval = 0;
  for(tstar = *minseglen; tstar <= (*n - *minseglen); tstar++){  //??
    R_CheckUserInterrupt(); //Has interrupred the R session? quits if true.
    RegARQuadCostFunc(Sumstats, phi, &size, &np1, m, &zero, &tstar, &seg1cost,
                      tol, error, shape, MBIC);
    if(*error != 0){
      goto err2;
    }
    RegARQuadCostFunc(Sumstats, phi, &size, &np1, m, &tstar, n, &seg2cost,
                      tol, error, shape, MBIC);
    if(*error != 0){
      goto err2;
    }
    tmplike[tstar-1] = seg1cost + seg2cost;
    neval++;
  }
  //Which tstar returns the minimum cost
  min_which(tmplike + *minseglen - 1, &neval, taulike, tau);
  *tau += *minseglen;
  
  
  //Free allocated memory
  err2: free(Sumstats);
  err1: return;
}




//Free allocated memory in case R session has been interrupred
void Free_CptRegAR_Normal_AMOC(int *error){
  // Error code from CptRegAR_Normal_AMOC C function, non-zero => error 
  if(*error==0){
    free((void *)Sumstats);
  }
}
