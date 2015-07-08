#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "dma.h"
#include "SteffenInterp.h"

using namespace std;


double signum(double n)
{
	double sign = 0.0;
	
	if (n >= 0.0) {sign = 1.0;}
	else {sign = -1.0;}
	
	return sign;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This is stupid. I couldn't compile using min and max algorithm...so i had to create one!
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double findMin(double a, double b)
{
	double min = 0.0;
	if (a<b) {min = a;}
	if (a>b) {min = b;}
	return min;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void calculateSteffenInterpolant(int N,double *x, double *y,double **c)
{
/*
!**** Calculate a piecewise cubic interpolant of a tabulated function (x_i,y_i) that is monotonic on
!**** each interval (x_i, x_i+1) for i=1,2,...,N-1. On each interval, the interpolant is of the form
!**** 
!****    f_i(x) = c_(1,i)*(x-x_i)^3 + c_(2,i)*(x-x_i)^2 
!****           + c_(3,i)*(x-x_i) + c_(4,i)
!**** 
!**** The interpolant and its first 3 derivatives can be evaluated by calling the SUBROUTINE 
!**** evaluateSteffenInterpolant. Note however that the interpolation scheme only guarantees continuity
!**** of the function and its first derivative, not of the higher order derivatives. 
!**** The method is described in detail in 
!****    Steffen, M., 1990, A&A 239, 443
!**** 
!**** Input parameters
!**** ----------------
!****    N: The number of data points
!****    x: 1-D array of size N containing the independent variable
!****    y: 1-D array of size N containing the dependent variable
!**** 
!**** Output parameters
!**** -----------------
!****    c: 2-D array of size 4x(N-1) containing the coefficients 
!****       of the interpolant on each interval (x_i,x_i+1) for 
!****       i=1,2,...,N-1
!**** 
!**** Coded by Bart Willems on 08/28/06
!**** Rewritten in c++ by Francesca Valsecchi 10/30/08
!****
*/
//	cout << "....calculateSteffenInterpolant.... " << endl;

    double 	p;
    double	*dy, *h, *s;
    	dy = dfun(N);
    	h = dfun(N-1);
    	s = dfun(N-1);
// Determine the sizes h_i of the interval (x_i,x_i+1) and the slope of the secants through (x_i,y_i) 
// and (x_i+1,y_i+1) for i=1,2,...,N-1 [Eqs. (6) and (7) in Steffen (1990)]
// HOW DOES THE LAST ELEMENTS LOOK LIKE?
	for (int i=0; i<N-1; i++)
     {
         h[i] = x[i+1]-x[i];
         s[i] = (y[i+1]-y[i])/h[i];
     }
// Determine the derivative at (x_i,y_i) for i=1,2,...,N [Eqs. (11), (26), (27) in Steffen (1990)]
    p = s[0]*(1.0+(h[0]/(h[0]+h[1]))) - s[1]*h[0]/(h[0]+h[1]);
    dy[0] = (signum(p)+signum(s[0]))*
    		findMin(fabs(s[0]),0.5*fabs(p));
	
	for (int i=1; i<N-1; i++)
	{
		p = (s[i-1]*h[i]+s[i]*h[i-1])/(h[i-1]+h[i]);
        dy[i] = (signum(s[i-1]) + signum(s[i]))*
        		findMin(fabs(s[i-1]),findMin(fabs(s[i]),0.5*fabs(p)));
	}
    p = s[N-2]*(1.0+(h[N-2]/(h[N-2]+h[N-3]))) - s[N-3]*h[N-2]/(h[N-2]+h[N-3]);
	dy[N-1] = (signum(p) + signum(s[N-2]))*
			findMin(fabs(s[N-2]),0.5*fabs(p));

// Determine the coefficients of the cubic interpolants on 
// the intervals (x_i,x_i+1) for i=1,2,...,N-1
// [Eqs. (2)-(5) in Steffen (1990)]
	for(int i=0; i<N-1; i++)
    {
		c[0][i] = (dy[i]+dy[i+1]-2.0*s[i])/(h[i]*h[i]);
        c[1][i] = (3.0*s[i]-2.0*dy[i]-dy[i+1])/h[i];
        c[2][i] = dy[i];
        c[3][i] = y[i];
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void evaluateSteffenInterpolant(int N,double *x, double **c, double x0,double *f)
{
/*
!**** Evaluate the Steffen interpolant and its first 3 derivatives at x0. On each interval 
!**** [x_i,x_i+1] with i=1,2,...,N-1, the interpolant is of the form
!**** 
!****    f_i(x) = c_(1,i)*(x-x_i)^3 + c_(2,i)*(x-x_i)^2 
!****           + c_(3,i)*(x-x_i) + c_(4,i)
!**** 
!**** The coefficients c_(1,i), c_(2,i), c_(3,i), c_(4,i) are determined by calling the SUBROUTINE 
!**** calculateSteffenInterpolant. Note however that the interpolation scheme only guarantees continuity 
!**** of the function and its first derivative, not of the higher order derivatives. The method is 
!**** described in detail in
!**** 
!****    Steffen, M., 1990, A&A 239, 443
!**** 
!**** Input parameters
!**** ----------------
!****    N: The number of data points
!****    x: 1-D array of size N containing the independent variable
!****    c: 2-D array of size 4x(N-1) containing the coefficients 
!****       of the interpolant on each interval (x_i,x_i+1) with
!****       i=1,2,...,N-1
!****    x0: The point at which the interpolant and its derivatives
!****        is to be evaluated
!**** 
!**** Output parameters
!**** -----------------
!****    f: 1-D array of size 4 containing the value of the interpolant
!****       and its first 3 derivatives at x0 
!****
!****       -> f(0) contains the function value
!****       -> f(1) contains the 1st derivative
!****       -> f(2) contains the 2nd derivative
!****       -> f(3) contains the 3rd derivative
!**** 
!**** Coded by Bart Willems on 08/28/06
!**** Rewritten in c++ by Francesca Valsecchi 10/30/08
!****
*/
	int i = 0;
// Check if x0 is in the range of the (x_i) data set
	if (x0 < x[0] || x0 > x[N-1])
	{
		cout << "FATAL ERROR: Extrapolation not allowed in " 
			 << "SUBROUTINE evaluateSteffenInterpolant!'" << endl;	         	
	   	exit(EXIT_FAILURE);
    }
// Locate x0 in the (x_i) data set
	for (int j=0; j<N-1; j++)
    {
         i = j;
         if (x0 < x[j+1]) {break;}         
    }
// Evaluate the interpolant and its first 3 derivatives at x0
      f[0] = c[0][i]*pow((x0-x[i]),3.0) + c[1][i]*pow((x0-x[i]), 2.0) + 
      		 c[2][i]*(x0-x[i]) + c[3][i];
      f[1] = 3.0*c[0][i]*pow((x0-x[i]),2.0) + 2.0*c[1][i]*(x0-x[i]) + c[2][i];
      f[2] = 6.0*c[0][i]*(x0-x[i]) + 2.0*c[1][i];
      f[3] = 6.0*c[0][i];
}

