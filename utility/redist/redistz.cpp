#include "mex.h"
#include "redist.hpp"
#include "defs.h"
#include "array2d.hpp"
#include "toolbox.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Matlab interface:
       [v,cpx,cpy] = redist(u,width,flag,dx,dy)
       where u is m x n 2D level set function to be redistanced,
             width is half-width of tubular region for redistancing,
  */

  double *u, *v, *cpx, *cpy, *cpz;
  double dx, dy, dz;
  int width, flag, m, n, N, ii;
  int *bnd;
  /* u holds the input array */
  u = mxGetPr(prhs[0]);
  width = static_cast<int>(mxGetScalar(prhs[1]));
  flag = static_cast<int>(mxGetScalar(prhs[2]));
  dx = mxGetScalar(prhs[3]);
  dy = mxGetScalar(prhs[4]);
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  N = m*n;
  
  plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m,n,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(m,n,mxREAL);
  v = mxGetPr(plhs[0]);
  cpx = mxGetPr(plhs[1]);
  cpy = mxGetPr(plhs[2]);
  // copy data into Array2D<double> type
  Array2D<double> uu(u,m,n,dx,dy);
  Redist r(uu,dx,dy,width,flag);
  r.redistance();
  // handle outputs
  r.dump_u(v);
  r.dump_cp(cpx,cpy);

}

