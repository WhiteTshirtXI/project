/* Random Number generators */
#include <math.h>

double rand_num(long *idum_ptr, double range1,double range2)
{
/* From Press et al. "Numerical Recipes": */
/* Long period (>2e18) random number generator of L'Ecuyer with Bays-Durham shuffle */
/* and added safeguards. Returns a uniform random deviate between range1 and range2 (exclusive of */
/* the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter */
/* idum between successive deviates in a sequence. RNMX should approximate the largest floating */
/* value that is less than 1. */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS) /*closest to 1 number of type double*/


  long j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;

  if (*idum_ptr <= 0)	
    {
      /*initialize */
      if (-(*idum_ptr) < 1)
	{
	  *idum_ptr = 1; /*avoid idum==0 */
	}
      else
	{
	  *idum_ptr = -(*idum_ptr);
	}
      idum2 = (*idum_ptr);

      for (j=NTAB+7; j>=0; j--)	 /*load the shuffle table after 8 warm-ups */
	{
	  k = (*idum_ptr)/IQ1;
	  *idum_ptr = IA1*(*idum_ptr-k*IQ1)-IR1*k;
	  if (*idum_ptr < 0)
	    {
	      *idum_ptr += IM1;
	    }

	  if (j < NTAB)
	    {
	      iv[j] = *idum_ptr;
	    }
	}
      iy = iv[0];
    }


  /*Start here when not initializing */
  k = (*idum_ptr)/IQ1;					    
  *idum_ptr = IA1*(*idum_ptr-k*IQ1)-k*IR1;  /*Compute idum= IA1*idum%IM1 without overflows */
  if (*idum_ptr < 0)
    {
      *idum_ptr += IM1;
    }
  k = idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2;          /*Compute idum2= IA2*idum%IM2 without overflows */
  if(idum2 < 0)
    {
      idum2 += IM2;
    }

  j = iy/NDIV;
  iy = iv[j]-idum2;
  iv[j] = *idum_ptr;
  if (iy<1)
    {
      iy += IMM1;
    }

  /*this routine NEVER returns range2, only numbers below it.*/
  if( (temp = AM*iy) > RNMX)
    {
      return range1+(range2-range1)*RNMX;  
    }
  else
    {
      return range1+(range2-range1)*temp;
    }
}




