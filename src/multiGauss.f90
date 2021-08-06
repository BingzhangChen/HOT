module mGf90
  use cholf90
  use gasdevf90
  implicit none

contains
 real function multigauss(R,mu,n) result(G)
!!$ Generates a one dimensional array of length n containing multivariate Gaussian pseudo-random numbers
!!$ where R is the lower-triangular Cholesky factor of the covariance matrix of the desired distribution 
!!$ and mu is a one dimensional array of length n, containing the mean value for each random variable. 
!!$ R must be a one dimensional array, containing the lower-triangular matrix stored by row, starting from the 
!!$ uppermost, leftmost entry (first row, first column). 
   implicit none
   integer, intent(in) :: n
   real, intent(in) :: R(n*(n+1)/2)
   real, intent(in) :: mu(n)
   real :: Nu(n)
   dimension :: G(n)
   integer :: i, j
   
   if( n*(n+1)/2 .ne. size(R) ) then
      write(6,*) ' n*(n+1)/2 != size(R) in multiGauss '
      write(6,*) ' Cholesky factor matrix has size ', size(R) 
      write(6,*) ' but you are requesting a vector of size ',n
      stop
   else
      
!!$ generate an array of independent Gaussian variates
      do j = 1, n
         Nu(j) = gasdev()
      end do

!!$ start with the desired mean, and add the product of 
!!$ [the lower-triangular Cholesky factor of the covariance matrix] x [ the vector of iid Gaussian variates]
      do i = 1, n
         G(i) = mu(i)
         do j =  1,i
            G(i) = G(i) + R( i*(i-1)/2 + j ) * Nu(j)
         end do
      end do

   end if
   return
 end function multigauss
end module mGf90
