module gammaf90
  use mtmod
  implicit none

contains
 real(8) function gammadev(alpha,beta)  result(g)
!!$ Generates pseudo-random numbers from a Gamma distribution with parameters alpha and beta, 
!!$ such that the mean of the distribution, mu = alpha * beta. 
!!$ (Different notations are used for the second parameter, sometimes defined as lambda = 1/beta)
!!$ The algorithms in this routine are from "Random Number Generation and Monte Carlo Methods", by James E. Gentle
!!$ (2nd edition, Springer-Verlag, New York, 2003)
   implicit none
   real(8), intent(in) :: alpha, beta
   real(8) :: am1, b, u1, u2, v, t, x, y
   integer :: i, j

   am1 = alpha - 1.0d0
   g = -1

   if( alpha .gt. 1.0d0 ) then  
      do while( g .le. 0.0d0 ) 
         !u1 = grnd()
         !u2 = grnd()
         call random_number(u1)
         call random_number(u2)
         v = (alpha - 1.0d0/(6.0d0*alpha))*u1/(am1*u2)
         if( 2.0d0 * (u2 - 1.0d0) /am1 + v + 1.0d0/v .le. 2.0d0 ) then
            g = am1 * v
         else if( 2*log(u2)/am1 - log(v) + v .le. 1.0d0 ) then
            g = am1 * v
         end if
      end do
   else ! alpha <= 1
      t = 0.07 + 0.75*sqrt(1.0d0 - alpha)
      b = 1.0d0 + exp(-t)*alpha/t
      do while( g .le. 0.0d0 ) 
         !u1 = grnd()
         !u2 = grnd()
         call random_number(u1)
         call random_number(u2)
         v = b * u1
         if( v .le. 1 ) then
            x = t * v**(1.0d0/alpha)
            if( u2 .lt. (2.0d0-x)/(2.0d0+x) ) then
               g = x
            else if( u2 .lt. exp(-x) ) then
               g = x
            end if
         else
            x = - log(t*(b-v)/alpha)
            y = x/t
            if( u2*(alpha + y*(1.0d0-alpha)) .le. 1.0d0 ) then
               g = x
            else if( u2 .lt. y**am1 ) then
               g = x
            end if
         end if
      end do
   end if
   g = beta * g
   return
 end function gammadev
end module gammaf90
