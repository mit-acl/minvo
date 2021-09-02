=================================================================
                 SubLiME package
                   Version 1.1
                 June. 20, 2003

 Copyright (C)  2002 SurfLab, Univ. of Florida ,
 See "COPYRIGHT.txt" for the copyright information.

 For the latest version and up-to-date information, please visit
 http://www.cise.ufl.edu/research/SurfLab/

 Authors  : Xiaobin Wu (xwu@cise.ufl.edu), 
            Jorg Peters(jorg@cise.ufl.edu)
 Contact  : Xiaobin Wu (xwu@cise.ufl.edu)

==================================================================





-------------------------------------------------------------------------------

(I) User's Guide:

------------------------------------------------------------------------------



 Overview
 ---------

   This package includes: 

   Routines:
      univariate     Bezier     slefe computation (uniSlefe())
      tensor product Bezier     slefe computation (tpSlefe ())
      univariate     bspline    slefe computation (bspSlefe())

   Data:
      bounds for univariate Bezier basis 
          (range/unirange%d_%d.asc, 
          where first integer is degree and second integer is the number of segments)
      bounds for univariate bspline basis 
          (range/bsprange%d_%d.asc, 
          where first integer is degree and second integer is the number of segments)


 Pre-Installation (For both unix and Windows)
 ----------------
     (i)   unzip the downloaded file

           Now you should have a SubLiME/ directory created.

     (ii)  Set up an environment variable named SUBLIMEPATH 
           to the SubLiME/ directory 

           for example, in csh:

           % setenv SUBLIMEPATH /usr/local/lib/SubLiME

           or in Windows

           My Computer -> Properties -> Advanced -> New System Variables
           Variable Name:  SUBLIMEPATH
           Variable Value: C:\download\Library\SubLiME


 Installation and Compilation (For unix)
 -----------------------------

     (i) Compile and use the library
        % make
        then include the SubLiME head file (SubLiME.h) in your source code
        and link the lib file libSubLiME.a to your executable.
        
     (ii) Check the example
         % make example
         % ./uniexample  
	  or % ./bspexample 

         * Hint -- using mouse to drag the control points *

 Installation and Compilation (For Windows)
 -----------------------------

     (i)  Open SubLiME.dsw using Visual Studio

     (ii) Compile project SubLiME to get a library file named: SubLiME.lib

     (iii) Compile and run projects 'uniexample' and 'bspexample' to 
           view the slefe of a bezier curve and a uniform bspline curve, respectively.

           * Hint -- using mouse to drag the control points *

     (iv)  To use the library in your own projects, include "SubLiME.h" and link "SubLiME.lib".
           (You either have to copy these two files into your project directory or
            set up the proper input directory in your project settings)

 Files:
 -------
      README        : this file
      SubLiME.h, 
      libSubLiME.a, : the library files (For Unix)
      SubLiME.lib,  : the library files (For Windows)
      range/*       : data file
      *.c           : source code


---------------------------------------------------------------------------------

(II) Programmer's Guide -- references of the functions supported by the library.

---------------------------------------------------------------------------------


  (1)
  /* --------------------------------------------------------
     void InitBounds();
    
     Initialize the bounds
    
     Important!! -----
      Call this functions once before you call any of the other 
      functions in the library.
     -------------------------------------------------------- */

  (2)
  /* --------------------------------------------------------
     int uniSlefe ( REAL *coeff, int stride, int deg, int seg, 
                         REAL* upper, REAL* lower, int stride_slefe);
    
     compute the slefe for a univariate function
    
     Input: 
             coeff:  Bezier coefficient array
            stride:  stride between two coefficients in the input array
               deg:  degree of the Bezier function
               seg:  number of linear segments in the slefe
      stride_slefe:  stride between two results in the output arrays
    
     Output: 
             upper:  the upper slefe for the function
             lower:  the lower slefe for the function
    ---------------------------------------------------------- */

  (3)
  /* --------------------------------------------------------
     int tpSlefe ( REAL *coeff, int strideu, int stridev, 
    			int degu, int degv, int segu, int segv, 
     			REAL* upper, REAL* lower, 
    			int stride_slefeu, int stride_slefev); 
    
     compute the slefe for a tensor-product bi-variate function
    
     Input: 
              coeff:  Bezier coefficient array
            strideu: 
            stridev:  strides between two coefficients in the input arrays
         degu, degv:  degrees of the Bezier function
         segu, segv:  number of linear segments in the slefe
      stride_slefeu:  
      stride_slefev:  stride: between two results in the output arrays
    
     Output: 
              upper:  the upper slefe for the function
              lower:  the lower slefe for the function
     -----------------------------------------------------------  */

  (4) 
  /* --------------------------------------------------------
     int bspSlefe(REAL* coeff, int cpnum, int deg, int stride, int seg, 
	              REAL* upper, REAL* lower, int slefe_stride)

     compute the slefe for a univariate bspline function (uniform knots)

     Input: 
              coeff: the coefficients of the spline 
              cpnum:  number of control points
                deg:  degree of the bspline
             stride:  stride between two coefficients in the input array
                seg:  number of segments on each interval [t^*_k , t^*_{k+1}]
                      grid.
       stride_slefe:  stride between two results in the output arrays

  Output: 
          upper:  the upper slefe for the function
          lower:  the lower slefe for the function
     1-d array with length: (cpnum-3)*seg+1
     (there are cpnum-2 second differences, so there are cpnum-3 segments 
       where each of the segments has seg pieces. 
 	    Therefore there will be (cpnum-3)*seg+1 points.)
    ---------------------------------------------------------- */

Last Modified:
June. 20th 2003
