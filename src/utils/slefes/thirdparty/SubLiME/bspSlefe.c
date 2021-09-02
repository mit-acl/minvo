#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include "type.h"
#include "SubLiME.h"

REAL* BspUpperBound[MAXDEG+1][MAXPCS+1];  // pre-tabulated upper bounds 
REAL* BspLowerBound[MAXDEG+1][MAXPCS+1];  // pre-tabulated lower bounds 

// Boolean flags show that if bounds for a certain degree/piece 
// has already been loaded  into memory.
int BspLoaded[MAXDEG+1][MAXPCS+1];   

extern char *sublime_path;

/* ---------------- codes start here ---------------------- */

/* --------------------------------------------------------
 * Initialize the flags for the bounds
 */
void InitBspBounds()
{
	int i,j;
	for(i=0; i<=MAXDEG; i++)
		for(j=0; j<=MAXPCS; j++)
	{
        BspLoaded[i][j]=0;
		BspUpperBound[i][j] = BspLowerBound[i][j] = NULL;
	}
}


/* -------------------------------------------
 * Load the bounds from the input file.
 *  Input:
 *         seg:  number of linear pieces on each side of the slefe
 */
int loadBspRange(int deg, int seg)
{
    char bndfilename[255];
    FILE* bndfile;
    REAL* P = (REAL*) malloc( sizeof(REAL) * (seg+1));
    REAL* M = (REAL*) malloc( sizeof(REAL) * (seg+1));
	int i;

    // Load bounds from the data file (ascii format)
    sprintf(bndfilename, "%s/range/bsprange-%d_%d.asc", sublime_path, deg, seg);
    if( (bndfile = fopen(bndfilename, "r")) == NULL)
    {
		printf("The bound file %s is not found.\n", bndfilename);
		return 1;  // return error if file not exist 
	}

    for(i=0;i<=seg;i++)
        fscanf(bndfile, "%lf", &P[i]);
    for(i=0;i<=seg;i++)
        fscanf(bndfile, "%lf", &M[i]);
    fclose(bndfile);


	// store them into the bound array
	BspUpperBound[deg][seg] = P;
	BspLowerBound[deg][seg] = M;

	BspLoaded[deg][seg] = 1;  // set the flag so that bounds for
                           //  this seg won't be loaded again

    return 0;  // successful return
}



/* --------------------------------------------------------
 * compute the slefe for a univariate bspline function (uniform knots)
 *
 * Input: 
 *     coeff: the coefficients of the spline 
 *     cpnum:  number of control points
 *       deg:  degree of the bspline
 *    stride:  stride between two coefficients in the input array
 *       seg:  number of segments on each interval [t^*_k , t^*_{k+1}]
 *                    grid.
 *  stride_slefe:  stride between two results in the output arrays
 *
 * Output: 
 *         upper:  the upper slefe for the function
 *         lower:  the lower slefe for the function
 *    1-d array with length: (cpnum-3)*seg+1
 *    (there are cpnum-2 second differences, so there are cpnum-3 segments 
 *      where each of the segments has seg pieces. 
 *	    Therefore there will be (cpnum-3)*seg+1 points.)
 */
int bspSlefe(REAL* coeff, int cpnum, int deg, int stride, int seg, 
			  REAL* upper, REAL* lower, int slefe_stride)
{
    int i,j;
	REAL* P, *M;

	if(!BspLoaded[deg][seg]) {    
		if(loadBspRange(deg, seg)) return 1;  // return error if not succeed 
	}

	// get the pointer from the bound array
	P = BspUpperBound[deg][seg];   
	M = BspLowerBound[deg][seg];

    // for each control point 
    for(i=1;i<cpnum-2;i++) 
    {
        REAL d2k  = coeff[(i-1)*stride] -2*coeff[i*stride]   + coeff[(i+1)*stride];
        REAL d2kp = coeff[i*stride]   -2*coeff[(i+1)*stride] + coeff[(i+2)*stride];

        for(j=0;j<=seg;j++) // for each gridpoint
        {
            REAL u=(REAL)j/seg; 
            REAL ctl_poly = (1-u)*coeff[i*stride] + u*coeff[(i+1)*stride];
            REAL h;

            int  pos = (i-1)*seg+j; // position at the result array

            h = max(d2k,0)*P[j] + min(d2k,0)*M[j] +
                max(d2kp,0)*P[seg-j] + min(d2kp,0)*M[seg-j];

            //printf("h= %f\n", h);
            // take the maximum at the joint vertices 
            if(j==0 && i!=1)
               upper[pos*(slefe_stride)] = max(upper[pos*(slefe_stride)],ctl_poly+h);
            else
               upper[pos*(slefe_stride)] = ctl_poly+h;
            // for deg=2 we can simply use 
            // upper[pos] = ctl_poly+h
            // is this true in general?

            h = min(d2k,0)*P[j] + max(d2k,0)*M[j] +
                min(d2kp,0)*P[seg-j] + max(d2kp,0)*M[seg-j];

            //printf("h= %f\n", h);
            // take the minimum at the joint vertices 
            if(j==0 && i!=1)
                lower[pos*(slefe_stride)] = min(lower[pos*(slefe_stride)], ctl_poly+h);
            else
                lower[pos*(slefe_stride)] = ctl_poly+h;
        }
    }

	return 0;
}


