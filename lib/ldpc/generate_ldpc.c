/* Copyright (c) 2000, 2001, 2006 by Radford M. Neal and Peter Junteng Liu
 *
 * Permission is granted for anyone to copy, use, modify, or distribute this
 * program and accompanying programs and documents for any purpose, provided
 * this copyright notice is retained and prominently displayed, along with
 * a note saying that the original programs are available from Radford Neal's
 * web page, and note is made of any changes made to the programs.  The
 * programs and documents are distributed without any warranty, express or
 * implied.  As the programs were written for research purposes only, they have
 * not been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own risk.
 *
 * Modified by Wai Lok Lai (amagicube@gmail.com), 2016.
 * Outputs the LDPC generated into MATLAB through the MEX functionality.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rand.h"
#include "alloc.h"
#include "intio.h"
#include "open.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"
#include "rcode.h"
#include "distrib.h"
#include "mex.h"

typedef enum { 
  Evencol, 	/* Uniform number of bits per column, with number specified */
  Evenboth 	/* Uniform (as possible) over both columns and rows */
} make_method;


void make_ldpc (int, int, int, make_method, distrib *, int, int);
int *column_partition (distrib *, int);


/* begin anguslai's contributions */

void mexFunction(
     int nlhs,       mxArray *plhs[],
     int nrhs, const mxArray *prhs[]
) {
  /* Declare variable */
  int M,N, NN;
  mwSize nzmax, oldnzmax;
  mwIndex *irs,*jcs,j,k;
  char *T;
  int Q, seed, evc, no4cycle, verbose;
  double *sr;
  double percent_sparse;
  distrib *d;
  make_method method;
  mod2entry *e;

  mwSize dopesize;
  double *dopepr;
  mwIndex cumdope;
  
  /* Check for proper number of input and output arguments */
  if (nrhs != 9) {
    mexErrMsgIdAndTxt( "MATLAB:ldpc_generate_n4c:invalidNumInputs",
            "Exactly 9 input arguments required.");
  }
  if (nlhs > 1) {
    mexErrMsgIdAndTxt( "MATLAB:ldpc_generate_n4c:maxlhs",
            "Exactly 1 output required");
  }
  
  /* Get the size and pointers to input data */
  dopesize = mxGetN(prhs[5]);
  dopepr = mxGetPr(prhs[5]);

  M  = (int) *(mxGetPr(prhs[0])); 
  NN  = (int) *(mxGetPr(prhs[1]));
  T  = mxArrayToString(prhs[2]); 
  seed = (int) *(mxGetPr(prhs[4]));

  evc = (int) *(mxGetPr(prhs[6]));
  no4cycle = (int) *(mxGetPr(prhs[7]));
  verbose = (int) *(mxGetPr(prhs[8]));

  N = NN - (int) dopesize;
  d = distrib_create(T);
  method = evc ? Evencol : Evenboth;
  /* Check for some problems. */
  
  if (distrib_max(d)>M) { 
    mexErrMsgIdAndTxt( "MATLAB:ldpc_generate_n4c:invalidDistribution",
            "At least one checks per bit (d) is greater than total checks (M)");
  }
  
  if (distrib_max(d)==M && N>1 && no4cycle) { 
    mexErrMsgIdAndTxt( "MATLAB:ldpc_generate_n4c:no4cycTooManyChecks",
            "Can't eliminate cycles of length four with this many checks per bit");
  }

  /* Allocate space for sparse matrix. Use ceil
   * to cause it to round up.
   */
  percent_sparse = atof(T) / M ; 
  nzmax=(mwSize)ceil((double)M*(double)N*percent_sparse);
  mexPrintf("LDPC %d...\n",nzmax);

  plhs[0] = mxCreateSparse(M,NN,nzmax,0);
  sr  = mxGetPr(plhs[0]);
  irs = mxGetIr(plhs[0]);
  jcs = mxGetJc(plhs[0]);

  /* Make the parity check matrix. */
  make_ldpc(M, N, seed,method,d,no4cycle,verbose);
  
  /* Copy nonzeros */
  cumdope = 0;
  k = 0;
  for (j = 0; j<mod2sparse_cols(H)+dopesize; j++) { 
    jcs[j] = k;
    if (dopesize && !(cumdope >= dopesize) && j == ((int) dopepr[cumdope])-1) { /*check dopesize>0 to avoid access array, matlab is 1 index, so -1*/
      cumdope++;
      jcs[j] = k; /*skip the column if it is in dopelist*/
    } else {
      e = mod2sparse_first_in_col(H,j-cumdope);
      while (!mod2sparse_at_end(e)) { 
        /* Check to see if non-zero element will fit in
         * allocated output array.  If not, increase percent_sparse
         * by 0.1/M, recalculate nzmax, and augment the sparse array
         */
        if (k>=nzmax){
          oldnzmax = nzmax;
          percent_sparse += .1/((double)M);
          nzmax = (mwSize)ceil((double)M*(double)N*percent_sparse);
          /* increase by at least 5 spots*/
          if (oldnzmax + 5 >= nzmax) {
            nzmax = oldnzmax + 5;
          }
          if (verbose)
            mexPrintf("... Matrix full, expanding to %d slots...\n", nzmax);

          mxSetNzmax(plhs[0], nzmax);
          mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
          mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(mwIndex)));
          
          sr  = mxGetPr(plhs[0]);
          irs = mxGetIr(plhs[0]);
        }

        sr[k] = 1;
        irs[k] = (int)mod2sparse_row(e);
        e = mod2sparse_next_in_col(e);
        k++;
      }
    }
  }
  jcs[j] = k;
  if (verbose)
    mexPrintf("... Done.\n");
  
  mod2sparse_free(H);

  return;
}

/* end anguslai's contributions */


/* CREATE A SPARSE PARITY-CHECK MATRIX.  Of size M by N, stored in H. */

void make_ldpc ( 
    int M,
    int N,
    int seed,		/* Random number seed */
    make_method method,	/* How to make it */
    distrib *d,		/* Distribution list specified */
    int no4cycle,		/* Eliminate cycles of length four? */
    int verbose
) {
  mod2entry *e, *f, *g, *h;
  int added, uneven, elim4, all_even, n_full, left;
  int i, j, k, t, z, cb_N;
  int *part, *u;
  
  rand_seed(10*seed+1);
  
  H = mod2sparse_allocate(M,N);
  part = column_partition(d,N);
  
  /* Create the initial version of the parity check matrix. */
  
  switch (method) {
    case Evencol: {
      z = 0;
      left = part[z];
      
      for (j = 0; j<N; j++) { 
        while (left==0) { 
          z += 1;
          if (z>distrib_size(d)){ 
            abort();
          }
          left = part[z];
        }
        for (k = 0; k<distrib_num(d,z); k++) { \
          do { 
            i = rand_int(M);
          } while (mod2sparse_find(H,i,j));
          mod2sparse_insert(H,i,j);
        }
        left -= 1;
      }
      
      break;
    }
    
    case Evenboth: {
      cb_N = 0;
      for (z = 0; z<distrib_size(d); z++) { 
        cb_N += distrib_num(d,z) * part[z];
      }
      
      u = chk_alloc (cb_N, sizeof *u);
      
      for (k = cb_N-1; k>=0; k--) { 
        u[k] = k%M;
      }
      
      uneven = 0;
      t = 0;
      z = 0;
      left = part[z];
      
      for (j = 0; j<N; j++) {
        while (left==0) { 
          z += 1;
          if (z>distrib_size(d)) { 
            abort();
          }
          left = part[z];
        }
        
        for (k = 0; k<distrib_num(d,z); k++) {
          for (i = t; i<cb_N && mod2sparse_find(H,u[i],j); i++) ;
          
          if (i==cb_N) { 
            uneven += 1;
            do { 
              i = rand_int(M);
            } while (mod2sparse_find(H,i,j));
            mod2sparse_insert(H,i,j);
          } else { 
            do { 
              i = t + rand_int(cb_N-t);
            } while (mod2sparse_find(H,u[i],j));
            mod2sparse_insert(H,u[i],j);
            u[i] = u[t];
            t += 1;
          }
        }
        
        left -= 1;
      }
      
      if (uneven>0) { 
        if (verbose)
          mexPrintf("... Had to place %d checks in rows unevenly\n",uneven);
      }
      break;
    }
    default: abort();
  }
  
  /* Add extra bits to avoid rows with less than two checks. */
  
  added = 0;
  
  for (i = 0; i<M; i++) { 
    e = mod2sparse_first_in_row(H,i);
    if (mod2sparse_at_end(e)) { 
      j = rand_int(N);
      e = mod2sparse_insert(H,i,j);
      added += 1;
    }
    e = mod2sparse_first_in_row(H,i);
    if (mod2sparse_at_end(mod2sparse_next_in_row(e)) && N>1) { 
      do { 
        j = rand_int(N);
      } while (j==mod2sparse_col(e));
      mod2sparse_insert(H,i,j);
      added += 1;
    }
  }
  
  if (added>0) { 
    if (verbose)
      mexPrintf(
          "... Added %d extra bit-checks to make row counts at least two\n",
          added);
  }
  
  /* Add extra bits to try to avoid problems with even column counts. */
  
  n_full = 0;
  all_even = 1;
  for (z = 0; z<distrib_size(d); z++) { 
    if (distrib_num(d,z)==M) { 
      n_full += part[z];
    }
    if (distrib_num(d,z)%2==1) { 
      all_even = 0;
    }
  }
  
  if (all_even && N-n_full>1 && added<2) { 
    int a;
    for (a = 0; added+a<2; a++) { 
      do { 
        i = rand_int(M);
        j = rand_int(N);
      } while (mod2sparse_find(H,i,j));
      mod2sparse_insert(H,i,j);
    }
    if (verbose)
      mexPrintf(
          "... Added %d extra bit-checks to try to avoid problems from even column counts\n",
          a);
  }
  
  /* Eliminate cycles of length four, if asked, and if possible. */
  
  if (no4cycle) {
    elim4 = 0;
    
    for (t = 0; t<10; t++) { 
      k = 0;
      for (j = 0; j<N; j++) { 
        for (e = mod2sparse_first_in_col(H,j);
            !mod2sparse_at_end(e);
            e = mod2sparse_next_in_col(e)) { 
          for (f = mod2sparse_first_in_row(H,mod2sparse_row(e));
              !mod2sparse_at_end(f);
              f = mod2sparse_next_in_row(f)) { 
            if (f==e) continue;
            for (g = mod2sparse_first_in_col(H,mod2sparse_col(f));
                !mod2sparse_at_end(g);
                g = mod2sparse_next_in_col(g)) { 
              if (g==f) 
                continue;
              for (h = mod2sparse_first_in_row(H,mod2sparse_row(g));
                  !mod2sparse_at_end(h);
                  h = mod2sparse_next_in_row(h)) { 
                if (mod2sparse_col(h)==j) { 
                  do { 
                    i = rand_int(M);
                  } while (mod2sparse_find(H,i,j));
                  mod2sparse_delete(H,e);
                  mod2sparse_insert(H,i,j);
                  elim4 += 1;
                  k += 1;
                  goto nextj;
                }
              }
            }
          }
        }
        nextj: ;
      }
      if (k==0) break;
    }
    
    if (elim4>0) { 
      if (verbose)
        mexPrintf(
            "... Eliminated %d cycles of length four by moving checks within column\n",
            elim4);
    }
    
    if (t==10) { 
      if (verbose)
        mexPrintf(
            "... Couldn't eliminate all cycles of length four in 10 passes\n");
    }
  }
}


/* PARTITION THE COLUMNS ACCORDING TO THE SPECIFIED PROPORTIONS.  It
 * may not be possible to do this exactly.  Returns a pointer to an
 * array of integers containing the numbers of columns corresponding
 * to the entries in the distribution passed. */

int *column_partition( 
    distrib *d,		/* List of proportions and number of check-bits */
    int n			/* Total number of columns to partition */
) {
  double *trunc;
  int *part;
  int cur, used;
  int i, j;
  
  trunc = chk_alloc (distrib_size(d), sizeof(double));
  part = chk_alloc (distrib_size(d), sizeof(int));
  
  used = 0;
  for (i = 0; i<distrib_size(d); i++) { 
    cur = floor(distrib_prop(d,i)*n);
    part[i] = cur;
    trunc[i] = distrib_prop(d,i)*n - cur;
    used += cur;
  }
  
  if (used>n) { 
    abort();
  }
  
  while (used<n) { 
    cur = 0;
    for (j = 1; j<distrib_size(d); j++) { 
      if (trunc[j]>trunc[cur]) { 
        cur = j;
      }
    }
    part[cur] += 1;
    used += 1;
    trunc[cur] = -1;
  }
  
  free(trunc);
  return part;
}
