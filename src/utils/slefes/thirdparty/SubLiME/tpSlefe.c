/* ======================================================================
 * File     : tpSlefe.c
 * Purpose  : To compute the upper/lower bounds for any given 
 *             tensor-product bi-variate function
 *
 * Author(s): Xiaobin Wu (xwu@cise.ufl.edu)
 * Date     : 03/04/2003
 * ======================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "type.h"
#include "SubLiME.h"

REAL ubuffer[(MAXDEG+1)*(MAXPCS+1)];
REAL lbuffer[(MAXDEG+1)*(MAXPCS+1)];

/* --------------------------------------------------------
 * compute the slefe for a tensor-product bi-variate function
 *
 * Input: 
 *          coeff:  Bezier coefficient array
 *        strideu: 
 *        stridev:  strides between two coefficients in the input arrays
 *     degu, degv:  degrees of the Bezier function
 *     pcsu, pcsv:  number of linear pieces in the slefe
 *  stride_slefeu:  
 *  stride_slefev:  stride: between two results in the output arrays
 *
 * Output: 
 *          upper:  the upper slefe for the function
 *          lower:  the lower slefe for the function
 */
int tpSlefe ( REAL *coeff, int strideu, int stridev, 
			int degu, int degv, int pcsu, int pcsv, 
			REAL* upper, REAL* lower, 
			int stride_slefeu, int stride_slefev) 
{
	int i;

	// The enclosure is computed by tensoring two uni-variate bounds.
	// We divide the process into two steps,
	// first, v direction. then,  u direction
	
	// step 1. compute bounds for each row of the control points
	//         bounds needed for degree = degv
	//
    for(i=0;i<=degu;i++) {
	    uniSlefe(&coeff[i*strideu], stridev, degv, pcsv, 
			&ubuffer[i*(MAXPCS+1)], &lbuffer[i*(MAXPCS+1)], 1);
	}
	// step 2. compute the bounds for the result from step 1, 
	//        but for each column, use bound for degree = dg1
	for(i=0;i<=pcsv;i++) {
		uniSlefe(&ubuffer[i], (MAXPCS+1), degu, pcsu, 
			&upper[i*stride_slefev], NULL, stride_slefeu);
		uniSlefe(&lbuffer[i], (MAXPCS+1), degu, pcsu, 
			NULL, &lower[i*stride_slefev], stride_slefeu);
    }

	return 0;   // success
}
