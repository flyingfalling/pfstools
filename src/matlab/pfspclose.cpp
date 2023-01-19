#include "compatibility.h"

#include <stdio.h>
#include <stdlib.h>

#include <config.h>

#if defined(_WIN32) || defined(_WIN64) 
#include <io.h>
#include <process.h>
#else
#include <unistd.h>
#endif

#include "mex.h"
#include "mex_utils.h"

#define SCRIPT_NAME "pfspclose"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{     
  /* Check for proper number of arguments. */
  if (nrhs != 1 && !is_mex_scalar( prhs[0] ) ) 
    mexErrMsgTxt( SCRIPT_NAME ": Expecting one parameter - file descriptor." );

  if( mxGetN(prhs[0]) != 2 )
    mexErrMsgTxt( SCRIPT_NAME ": File descriptor not created with pfspopen." );
  
   double *p_fid = mxGetPr( prhs[0] );

   // cast 64-bit number into a pointer   
   FILE *fh  = (FILE*)*((unsigned long long*)(p_fid+1));
  
    if( pclose( fh ) == -1 )
		mexErrMsgTxt( SCRIPT_NAME ": pclose has failed.");
  
  return;
}
