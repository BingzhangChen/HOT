subroutine matmuls(x,v,p)
USE PARAM_MOD, only: NPar
implicit none
real, intent(in) :: v(NPar)
real, intent(in) :: x(NPar,NPar)
real, intent(out):: p(NPar)

integer :: n

n = size(v,1)
!    SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
call DSYMV('u',n,1d0,x,n,v,1,0d0,p,1)
end subroutine matmuls

