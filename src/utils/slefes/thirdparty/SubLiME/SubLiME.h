#ifndef SUBLIME_H_2002_11_4
#define SUBLIME_H_2002_11_4

/* =================================================================
 * This package provides functions involving slefe construction
 *
 * For example:
 *     uniSlefe(): compute the slefe for a univariate function
 *      tpSlefe(): compute the slefe for a tensor-product bi-variate function
 *
 * Arthor(s):
 *      Xiaobin Wu (xwu@cise.ufl.edu)
 *
 * Copyright (c) SurfLab, University of Florida, 2002
 * =================================================================
 */ 

#ifdef __cplusplus
extern "C" {
#endif

#ifndef REAL
#define REAL double
#endif

/* --------------------------------------------------------
 * Initialize the bounds
 *
 * --- Important!! -----
 *  Call this functions once before you call any of the other 
 *  functions in the library.
 */
void InitBounds();

/* --------------------------------------------------------
 * compute the slefe for a univariate function
 *
 * Input: 
 *         coeff:  bezier coefficient array
 *        stride:  stride between two coefficients in the input array
 *           deg:  degree of the bezier function
 *           seg:  number of linear segments in the slefe
 *  stride_slefe:  stride between two results in the output arrays
 *
 * Output: 
 *         upper:  the upper slefe for the function
 *         lower:  the lower slefe for the function
 */
int uniSlefe ( REAL *coeff, int stride, int deg, int seg, 
                       REAL* upper, REAL* lower, int stride_slefe);

/* --------------------------------------------------------
 * compute the slefe for a tensor-product bi-variate function
 *
 * Input: 
 *          coeff:  bezier coefficient array
 *        strideu: 
 *        stridev:  strides between two coefficients in the input arrays
 *     degu, degv:  degrees of the bezier function
 *     segu, segv:  number of linear segments in the slefe
 *  stride_slefeu:  
 *  stride_slefev:  stride: between two results in the output arrays
 *
 * Output: 
 *          upper:  the upper slefe for the function
 *          lower:  the lower slefe for the function
 */
int tpSlefe ( REAL *coeff, int strideu, int stridev, 
				int degu, int degv, int segu, int segv, 
				REAL* upper, REAL* lower, 
				int stride_slefeu, int stride_slefev); 

/* --------------------------------------------------------
 * compute the slefe for a univariate bspline function (uniform knots)
 *
 * Input: 
 *     coeff: the coefficients of the spline 
 *     cpnum:  number of control points
 *       deg:  degree of the bspline
 *    stride:  stride between two coefficients in the input array
 *       seg:  number of segments on each interval [t^*_k , t^*_{k+1}]
 *                    gridding
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
			  REAL* upper, REAL* lower, int slefe_stride);

#ifdef __cplusplus
}
#endif

#endif /* SUBLIME_H_2002_11_4 */
