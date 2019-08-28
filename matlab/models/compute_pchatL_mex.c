#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "float.h"

/*
 * compute_pchatL_mex.c
 *
 *COMPUTE_PCHATL Compute probability of responding Left.
 *
 * ================ INPUT VARIABLES ====================
 * NC: number of signed contrast levels. [scalar] (integer)
 * NX: number of grid points. [scalar] (integer)
 * NP: number of probability points. [scalar] (integer)
 * LOGLIKEDIFF_T: log likelihood of Left minus Right. [Nx,Nc] (double)
 * LOGPRIOR_ODDS: log prior odds. [1,Np] (double)
 * SOFTMAX_ETA: softmax inverse temperature. [scalar] (double)
 * SOFTMAX_BIAS: softmax bias. [scalar] (double)
 * W: grid weight vector. [1,Nx] (double)
 * 
 * ================ OUTPUT VARIABLES ==================
 * PCHATL: probability of responding left. [Nc,Np] (double)
 *
 * This is a MEX-file for MATLAB.
 * Template C code generated on 26-Aug-2019 with MEXXER v0.2 
 * (https://github.com/lacerbi/mexxer).
 */

/* Set ARGSCHECK to 0 to skip argument checking (for minor speedup) */
#define ARGSCHECK 1

void compute_pchatL( double *PChatL, int Nc, int Nx, int Np, double *loglikediff_t, double *logprior_odds, double softmax_eta, double softmax_bias, double *W )
{

    
    int ip,ic,ix;
    double dhat,Inf,lodds,PCx_soft;
    
    Inf = mxGetInf();   // Value of Inf
	
    
    if (softmax_eta == Inf) {        

        for (ip = 0; ip < Np; ip++) {
            
            lodds = *(logprior_odds++);

            for (ic = 0; ic < Nc; ic++, PChatL++) {
                
                *PChatL = 0.;

                dhat = *(loglikediff_t++) + lodds;                    
                PCx_soft = (dhat > -softmax_bias) ? 1.0 : 0.0;
                if (dhat == -softmax_bias) PCx_soft = 0.5;                    
                *PChatL += 0.5 * PCx_soft * *(W++);
                
                for (ix = 1; ix < Nx-1; ix++) {

                    dhat = *(loglikediff_t++) + lodds;                    
                    PCx_soft = (dhat > -softmax_bias) ? 1.0 : 0.0;
                    if (dhat == -softmax_bias) PCx_soft = 0.5;                    
                    *PChatL += PCx_soft * *(W++);
                    
                }

                dhat = *(loglikediff_t++) + lodds;                    
                PCx_soft = (dhat > -softmax_bias) ? 1.0 : 0.0;
                if (dhat == -softmax_bias) PCx_soft = 0.5;                    
                *PChatL += 0.5 * PCx_soft * *(W++);
                
                W -= Nx;
                
            }

            loglikediff_t -= Nx*Nc;
        }
    }	
    else {
        
        for (ip = 0; ip < Np; ip++) {
            
            lodds = *(logprior_odds++);

            for (ic = 0; ic < Nc; ic++, PChatL++) {
                
                *PChatL = 0.;

                dhat = *(loglikediff_t++) + lodds;
                PCx_soft = 1./(1. + exp(-softmax_eta*(dhat + softmax_bias)));
                if (dhat == -softmax_bias) PCx_soft = 0.5;
                *PChatL += 0.5 * PCx_soft * *(W++);
                
                for (ix = 1; ix < Nx-1; ix++) {

                    dhat = *(loglikediff_t++) + lodds;
                    PCx_soft = 1./(1. + exp(-softmax_eta*(dhat + softmax_bias)));
                    if (dhat == -softmax_bias) PCx_soft = 0.5;                    
                    *PChatL += PCx_soft * *(W++);
                    
                }

                dhat = *(loglikediff_t++) + lodds;
                PCx_soft = 1./(1. + exp(-softmax_eta*(dhat + softmax_bias)));
                if (dhat == -softmax_bias) PCx_soft = 0.5;                    
                *PChatL += 0.5 * PCx_soft * *(W++);
                
                W -= Nx;
                
            }

            loglikediff_t -= Nx*Nc;
        }    
        
        
        
        
    }
    
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *PChatL, *loglikediff_t, *logprior_odds, softmax_eta, softmax_bias, *W;
	int Nc, Nx, Np;
#if ( ARGSCHECK==0 )
#else /* ( ARGSCHECK!=0 ) */ 
	mwSize *dims_loglikediff_t, *dims_logprior_odds, *dims_W;
#endif /* ( ARGSCHECK!=0 ) */ 

	/*  check for proper number of arguments */
	/* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
	   within an if statement, because it will never get to the else
	   statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks
	   you out of the MEX-file) */
	if ( nrhs<8 || nrhs>8 )
		mexErrMsgIdAndTxt( "MATLAB:compute_pchatL:invalidNumInputs",
			"Eight inputs required.");
	if ( nlhs<1 || nlhs>1 )
		mexErrMsgIdAndTxt( "MATLAB:compute_pchatL:invalidNumOutputs",
			"One outputs required.");

	/* Get first input (NC, scalar int) */
	Nc = (int) mxGetScalar(prhs[0]);

	/* Get second input (NX, scalar int) */
	Nx = (int) mxGetScalar(prhs[1]);

	/* Get third input (NP, scalar int) */
	Np = (int) mxGetScalar(prhs[2]);

	/* Get fourth input (LOGLIKEDIFF_T, Nx-by-Nc double) */
	loglikediff_t = (double*) mxGetPr(prhs[3]);

	/* Get fifth input (LOGPRIOR_ODDS, 1-by-Np double) */
	logprior_odds = (double*) mxGetPr(prhs[4]);

	/* Get sixth input (SOFTMAX_ETA, scalar double) */
	softmax_eta = (double) mxGetScalar(prhs[5]);

	/* Get seventh input (SOFTMAX_BIAS, scalar double) */
	softmax_bias = (double) mxGetScalar(prhs[6]);

	/* Get eighth input (W, 1-by-Nx double) */
	W = (double*) mxGetPr(prhs[7]);

	/* Check sizes of input arguments (define ARGSCHECK to 0 above to skip) */
#if ( ARGSCHECK==1 )
		if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetN(prhs[0])*mxGetM(prhs[0])!=1) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:NcNotScalar", "Input NC must be a scalar.");

		if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetN(prhs[1])*mxGetM(prhs[1])!=1) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:NxNotScalar", "Input NX must be a scalar.");

		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || (mxGetN(prhs[2])*mxGetM(prhs[2])!=1) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:NpNotScalar", "Input NP must be a scalar.");

		dims_loglikediff_t = (mwSize*) mxGetDimensions(prhs[3]);
		if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
				mexErrMsgIdAndTxt("MATLAB:compute_pchatL:loglikediff_tNotReal", "Input LOGLIKEDIFF_T must be real.");
		if ( dims_loglikediff_t[0] != ((mwSize) (Nx)) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:loglikediff_tWrongSize", "The first dimension of input LOGLIKEDIFF_T has the wrong size (should be Nx).");
		if ( dims_loglikediff_t[1] != ((mwSize) (Nc)) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:loglikediff_tWrongSize", "The second dimension of input LOGLIKEDIFF_T has the wrong size (should be Nc).");

		dims_logprior_odds = (mwSize*) mxGetDimensions(prhs[4]);
		if ( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) )
				mexErrMsgIdAndTxt("MATLAB:compute_pchatL:logprior_oddsNotReal", "Input LOGPRIOR_ODDS must be real.");
		if ( dims_logprior_odds[0] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:logprior_oddsWrongSize", "The first dimension of input LOGPRIOR_ODDS has the wrong size (should be 1).");
		if ( dims_logprior_odds[1] != ((mwSize) (Np)) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:logprior_oddsWrongSize", "The second dimension of input LOGPRIOR_ODDS has the wrong size (should be Np).");

		if ( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || (mxGetN(prhs[5])*mxGetM(prhs[5])!=1) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:softmax_etaNotScalar", "Input SOFTMAX_ETA must be a scalar.");

		if ( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || (mxGetN(prhs[6])*mxGetM(prhs[6])!=1) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:softmax_biasNotScalar", "Input SOFTMAX_BIAS must be a scalar.");

		dims_W = (mwSize*) mxGetDimensions(prhs[7]);
		if ( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) )
				mexErrMsgIdAndTxt("MATLAB:compute_pchatL:WNotReal", "Input W must be real.");
		if ( dims_W[0] != ((mwSize) (1)) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:WWrongSize", "The first dimension of input W has the wrong size (should be 1).");
		if ( dims_W[1] != ((mwSize) (Nx)) )
			mexErrMsgIdAndTxt("MATLAB:compute_pchatL:WWrongSize", "The second dimension of input W has the wrong size (should be Nx).");
#endif /* ( ARGSCHECK==1 ) */ 

	/* Pointer to first output (PCHATL, Nc-by-Np double) */
	plhs[0] = mxCreateDoubleMatrix((mwSize) (Nc), (mwSize) (Np), mxREAL);
	PChatL = mxGetPr(plhs[0]);

	/* Call the C subroutine */
	compute_pchatL(PChatL, Nc, Nx, Np, loglikediff_t, logprior_odds, softmax_eta, softmax_bias, W);

}
