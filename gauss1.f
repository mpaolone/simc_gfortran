	real*8 function gauss1(nsigmax)

**************************************************************************
* This subroutine generates a random number distributed about a Gaussian
* centered at zero with a standrd deviation of 1.  The gaussian is
* truncated at nsigmax.  Use as a function like: VAL' = VAL + WIDTH*X
*
* Algorithm is from PHYSICS LETTERS B V.204, P.83 ''Review of particle
* properties", Particle Data Group. (with nsigmax cutoff added).
**************************************************************************

	implicit none

	real*8 u1,u2,v1,v2,s,nsigmax
	real*8 grnd

1	u1 = grnd()
	u2 = grnd()
	v1 = 2.0*u1-1.0
	v2 = 2.0*u2-1.0
	s  = v1**2+v2**2

	if (s.gt.1. .or. s.eq.0) goto 1

	gauss1 = v1*sqrt(-2.*log(s)/s)  ! <--want a natural log here

	if (abs(gauss1).gt.nsigmax) goto 1	!truncate at nsigmax

	return
	end
	
**************************************************************************
	
	real*8 function cauchy1(nsigmax)
	
**************************************************************************
* This subroutine generates a random number distributed about a Cauchy
* distribution centered at zero with a standerd deviation of 1.  It is
* truncated at nsigmax.  Use as a function like: VAL' = VAL + WIDTH*X
**************************************************************************

	implicit none

	real*8 u1,u2,v1,v2,s,nsigmax
	real*8 grnd

1	u1 = grnd()
	u2 = grnd()
	v1 = 2.0*u1-1.0
	v2 = 2.0*u2-1.0
	s  = v1**2+v2**2

	if (s.gt.1. .or. s.eq.0) goto 1
	
	cauchy1 = u1/u2;
	
	if (abs(cauchy1).gt.nsigmax) goto 1	!truncate at nsigmax
	
	return
	end
