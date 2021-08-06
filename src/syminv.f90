module syminvf90
  use cholf90
contains
  subroutine syminv ( a, n, c, w, nullty, ifault )
    !*****************************************************************************80
    !
    !! SYMINV computes the inverse of a symmetric matrix.
    !
    !  Modified:
    !
    !    11 February 2008
    !
    !  Author:
    !
    !    FORTRAN77 version by Michael Healy
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    Michael Healy,
    !    Algorithm AS 7:
    !    Inversion of a Positive Semi-Definite Symmetric Matrix,
    !    Applied Statistics,
    !    Volume 17, Number 2, 1968, pages 198-199.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A((N*(N+1))/2), a positive definite matrix stored 
    !    by rows in lower triangular form as a one dimensional array, in the sequence
    !    A(1,1),
    !    A(2,1), A(2,2),
    !    A(3,1), A(3,2), A(3,3), and so on.
    !
    !    Input, integer ( kind = 4 ) N, the order of A.
    !
    !    Output, real ( kind = 8 ) C((N*(N+1))/2), the inverse of A, or generalized
    !    inverse if A is singular, stored using the same storage scheme employed
    !    for A.  The program is written in such a way that A and U can share storage.
    !
    !    Workspace, real ( kind = 8 ) W(N).
    !
    !    Output, integer ( kind = 4 ) NULLTY, the rank deficiency of A.  If NULLTY is zero,
    !    the matrix is judged to have full rank.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error indicator.
    !    0, no error detected.
    !    1, N < 1.
    !    2, A is not positive semi-definite.
    !
    implicit none
    
    integer ( kind = 4 ), intent(in) :: n
    
    real    ( kind = 8 ), intent(in)  :: a((n*(n+1))/2)
    real    ( kind = 8 ), intent(out) :: c((n*(n+1))/2)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) icol
    integer ( kind = 4 ), intent(out) :: ifault
    integer ( kind = 4 ) irow
    integer ( kind = 4 ) j
    integer ( kind = 4 ) jcol
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) mdiag
    integer ( kind = 4 ) ndiag
    integer ( kind = 4 ) nn
    integer ( kind = 4 ) nrow
    integer ( kind = 4 ), intent(out) :: nullty
    real    ( kind = 8 ) w(n)
    real    ( kind = 8 ) x
    
    ifault = 0
    
    if ( n <= 0 ) then
       ifault = 1
       return
    end if
    
    nrow = n
    !
    !  Compute the Cholesky factorization of A.
    !  The result is stored in C.
    !
    nn = ( n * ( n + 1 ) ) / 2
    
    call cholesky ( a, n, nn, c, nullty, ifault )
    
    if ( ifault /= 0 ) then
       return
    end if
    !
    !  Invert C and form the product (Cinv)' * Cinv, where Cinv is the inverse
    !  of C, row by row starting with the last row.
    !  IROW = the row number, 
    !  NDIAG = location of last element in the row.
    !
    irow = nrow
    ndiag = nn
    
    do
       !
       !  Special case, zero diagonal element.
       !
       if ( c(ndiag) == 0.0D+00 ) then
          
          l = ndiag
          do j = irow, nrow
             c(l) = 0.0D+00
             l = l + j
          end do
          
       else
          
          l = ndiag
          do i = irow, nrow
             w(i) = c(l)
             l = l + i
          end do
          
          icol = nrow
          jcol = nn
          mdiag = nn
          
          do
             
             l = jcol
             
             if ( icol == irow ) then
                x = 1.0D+00 / w(irow)
             else
                x = 0.0D+00
             end if
             
             k = nrow
             
             do while ( irow < k )
                
                x = x - w(k) * c(l)
                k = k - 1
                l = l - 1
                
                if ( mdiag < l ) then
                   l = l - k + 1
                end if
                
             end do
             
             c(l) = x / w(irow)
             
             if ( icol <= irow ) then
                exit
             end if
             
             mdiag = mdiag - icol
             icol = icol - 1
             jcol = jcol - 1
             
          end do
          
       end if
       
       ndiag = ndiag - irow
       irow = irow - 1
       
       if ( irow <= 0 ) then
          exit
       end if
       
    end do
    
    return
  end subroutine syminv
end module syminvf90
