#include "mex.h"

void usage() {
  printf("/*\n");
  printf("* =============================================================\n");
  printf("* function g = statistical_inefficiency_mex(A_t, B_t);\n");
  printf("*\n");
  printf("* Compute the statistical inefficiency.\n");
  printf("* \n");
  printf("* arguments;\n");
  printf("*  A(t) is snapshot t of observable A.\n");
  printf("*  B(t) is snapshot t of observable B.\n");
  printf("* \n");
  printf("* This is a MEX-file for MATLAB.  \n");
  printf("* =============================================================\n");
  printf("*/\n\n");
}

/* If you are using a compiler that equates NaN to zero, you must
 * compile this example using the flag -DNAN_EQUALS_ZERO. For 
 * example:
 *
 *     mex -DNAN_EQUALS_ZERO findnz.c  
 *
 * This will correctly define the IsNonZero macro for your
   compiler. */
#if NAN_EQUALS_ZERO
#define IsNonZero(d) ((d) != 0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d) != 0.0)
#endif

#define max(a,b) ( ((a) < (b)) ? (b) : (a) )
#define min(a,b) ( ((a) < (b)) ? (a) : (b) )
#define SQR(x) ((x)*(x))

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  int T;
  int ndim;
  const int * dim_array;
  int ntau;
  int tau_max;
  double * A_t;
  double * B_t;
  double C;
  double * g;
  int tau;
  int tau_step;
  int n;
  int t0;
  double EA; /* Expectation of A over trajectory. */
  double EB; /* Expectation of B over trajectory. */
  double EAB; /* Expectation of AB over trajectory. */

  /* Check for proper number of input and output arguments. */
  if (nrhs != 2 || nlhs != 1) {
    usage();
    mexErrMsgTxt("Wrong number of input/output arguments.");
  } 

  /* Get the number of snapshots provided. */
  ndim = mxGetNumberOfDimensions(prhs[0]);
  dim_array = mxGetDimensions(prhs[0]);
  T = max(dim_array[0], dim_array[1]);
  /* printf("T = %d\n", T); */

  /* Get observable. */
  A_t = (double *)mxGetPr(prhs[0]);
  B_t = (double *)mxGetPr(prhs[1]);
  
  /* Allocate space for output. */
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  g = (double *)mxGetPr(plhs[0]);  
  *g = 0.5;

  /* Compute expectation of A and A^2 over trajectory. */
  EA = 0.0;
  EB = 0.0;
  EAB = 0.0;
  for(t0 = 0; t0 < T; t0++) {
    EA += A_t[t0];
    EB += B_t[t0];
    EAB += A_t[t0] * B_t[t0];
  }
  EA /= (double)T;
  EB /= (double)T;
  EAB /= (double)T;

  /* printf("EA = %e, EB = %e, EAB = %e, EAB-EA*EB = %e\n", EA, EB, EAB, EAB - EA*EB); */

  /* Return if there is a problem. */
  if(EAB - EA*EB == 0)
    return;

  /* Compute unnormalized time-correlation function. */
  tau_step = 1;
  for(tau = 1; tau < T; tau += tau_step)
    {
      /*  Compute TCF. */
      C = 0.0;

      for(t0 = 0; t0 < T - tau; t0++) {
	C += (A_t[t0] - EA) * (B_t[t0+tau] - EB) + (B_t[t0] - EB) * (A_t[t0+tau] - EA);
      }

      C /= 2.0 * (double)(T - tau);
      C /= EAB - EA*EB;

      if(C > 0) {
	if(tau == 1)
	  *g += C * (1. - (double)tau/T);
	else if(tau == 2)
	  *g += C * (1. - (double)tau/T) * 1.5;
	else
	  *g += C * (1. - (double)tau/T) * (1 + 0.5*(tau_step-1) + 0.5 * (tau_step));
      }
      else
	break;

      if(tau > 1)
	tau_step++;
    }
  
  *g *= 2.0;

  /* printf("Computed out to %d/%d\n", tau, T); */
}
