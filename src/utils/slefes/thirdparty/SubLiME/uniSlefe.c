/* ======================================================================
 * File     : uniSlefe.c
 * Purpose  : To compute the upper/lower bounds for any given 
 *             univariate function
 *
 * Author(s): Xiaobin Wu (xwu@cise.ufl.edu)
 * Date     : 11/04/2002
 * ======================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "type.h"
#include "SubLiME.h"

// The bound for the basis functions
REAL* upper_bound[MAXDEG+1][MAXPCS+1];  // pre-tabulated upper bounds 
REAL* lower_bound[MAXDEG+1][MAXPCS+1];  // pre-tabulated lower bounds 

// flags show that if bounds for a certain degree/piece 
// has already been loaded  into memory.
int   loaded[MAXDEG+1][MAXPCS+1];   

void InitBspBounds();   // in bspslefe.c

char *sublime_path;

/* ---------------- codes start here ---------------------- */

/* --------------------------------------------------------
 * Initialize the flags for the uni-variate bounds
 *
 * For compilers that initiate variables to zeros by default,
 * this function has no effect.
 */
void InitUniBounds()
{
	int i,j;
	for(i=0; i<=MAXDEG; i++)
		for(j=0; j<=MAXPCS; j++)
	{
        loaded[i][j]=0;
		upper_bound[i][j] = lower_bound[i][j] = NULL;
	}
}
 

/* --------------------------------------------------------
 * Initialize all the bounds
 *
 */
void InitBounds() {

	sublime_path = getenv("SUBLIMEPATH");
	if(sublime_path == NULL) {
		printf("Environment variable SUBLIMEPATH not found.\n");
		exit(0);
	}

	if(sublime_path == NULL) {
		printf("Environment varaible SUBLIMEPATH not found.\n");
		exit(-1);
	}
	InitUniBounds();
	InitBspBounds();
}

/* -------------------------------------------
 * Load the bounds from the input file.
 *  Input:
 *         deg:  degree of the Bezier function
 *         seg:  number of linear pieces in the slefe
 */
int loadUniRange(int deg, int seg)
{
    char filename[255];  // name of the file contains the bounds
    FILE* rangefile;     // file handler for the bound data

    REAL *P, *M;         // pointer to the memory saving the bounds
   	                     //   (P for plus and M for Minus)
	int bas, pts;        // bas: number of base functions
	                     // pts: number of breaking points in the bound

    int i,j;

    bas = deg-1;   // number of base functions 
	pts = seg+1;   // number of breaking points

	// allocate the memory
	P = (REAL*) malloc( (bas*pts) *sizeof(REAL));
	M = (REAL*) malloc( (bas*pts) *sizeof(REAL));

    sprintf(filename, "%s/range/unirange-%d_%d.asc", sublime_path, deg, seg);
    //sprintf(filename, "range/unirange-%d_%d.best", deg+1, seg+1);

    if ((rangefile = fopen(filename, "r")) == NULL)
    {
		printf("The bound file %s is not found.\n", filename);
        return 1;  // return error if file not exist 
	}

	// read in the bounds
    for (i=0;i<bas; i++) {
      for (j=0;j<pts; j++)
         fscanf(rangefile, "%lf", &P[i*pts+j]);

      for (j=0;j<pts; j++)
         fscanf(rangefile, "%lf", &M[i*pts+j]);
	}

    fclose(rangefile);

	// store them into the bound array
	upper_bound[deg][seg] = P;
	lower_bound[deg][seg] = M;

	loaded[deg][seg] = 1;  // set the flag so that bounds for
                           //  this deg/seg won't be loaded again

    return 0;  // successful return
}


/* --------------------------------------------------------
 * compute the slefe for a univariate function
 *
 * Input: 
 *         coeff:  bezier coefficient array
 *        stride:  stride between two coefficients in the input array
 *           deg:  degree of the bezier function
 *           seg:  number of linear pieces in the slefe
 *  stride_slefe:  stride between two results in the output arrays
 *
 * Output: 
 *         upper:  the upper slefe for the function
 *         lower:  the lower slefe for the function
 */
int uniSlefe ( REAL *coeff, int stride, int deg, int seg, 
                       REAL* upper, REAL* lower, int stride_slefe)
{
    int bas = deg-1;   // num of basis (also num of 2nd differences)
	int pts = seg+1;   // num of breaking points in slefe

	REAL *P, *M;       // pointer to the memory saving the bounds
    REAL D2b[MAXDEG];  // 2nd difference array
	REAL left_end, right_end;

    int i,j;

	// load the bounds if necessary
	if(!loaded[deg][seg]) {    
		if(loadUniRange(deg, seg)) return 1;  // return error if not succeed 
	}

	// get the pointer from the bound array
	P = upper_bound[deg][seg];   
	M = lower_bound[deg][seg];

	// compute the 2nd difference of the input function
    for (i=0;i<bas;i++)
        D2b[i] = coeff[i*stride] -2*coeff[(i+1)*stride] + coeff[(i+2)*stride];

	// compute the function bounds according to P, M, and D2b
	
	left_end  = coeff[0*stride];   // save two end points
	right_end = coeff[deg*stride];

	// compute each break point
    for(i=0;i<=seg;i++) {
        double u = (double)i/seg;
		int    loc = i*stride_slefe;  // location at the result array

	    // 1. initialize to linear average of two end points
		if(upper!=NULL)
			upper[loc] = (1-u)*left_end + u*right_end;
		if(lower!=NULL)
			lower[loc] = (1-u)*left_end + u*right_end;

	    // 2. add contribution from every coefficients
        for(j=0;j<bas;j++) {
			if(D2b[j]>0) {
		        if(upper!=NULL)
                    upper[loc] += P[j*pts+i]* D2b[j];
		        if(lower!=NULL)
                    lower[loc] += M[j*pts+i]* D2b[j];
		    } else {
		        if(upper!=NULL)
                    upper[loc] += M[j*pts+i]* D2b[j];
		        if(lower!=NULL)
                    lower[loc] += P[j*pts+i]* D2b[j];
			}
		}
    }

	return 0;  // successful return
}

