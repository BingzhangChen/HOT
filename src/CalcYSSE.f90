!-----------------------------------------------------------------------
!Major computation is called here. 
! Calculates model output and Sum of Squares
SUBROUTINE CalcYSSE(pars,MODVAL,OBSVAL,ssqe)
USE MOD_1D
implicit none
! parameter sets for calculating sum of squares
real, intent(in)  :: pars(NPar)    
real, intent(OUT) :: modval(ANobs)
! observational values
real, intent(out) :: obsval(ANobs)

! Minimal and maximal values for transformation
real              :: max_y,min_y
real, allocatable :: OBS_data(:), MOD_data(:)
! Sum of squares for each data type
real, intent(OUT) :: ssqe(NDTYPE * Nstn)          
integer, parameter:: logtransform = NO
integer, parameter:: backlog      = YES
integer           :: don,i,k,j,kk(NDTYPE)

obsval(:) = 0d0
modval(:) = 0d0
! Make the model read in parameters generated in the main program:      
params(:) = pars(:)

call Timestep

! Match the final year result of the model output to the observational data 
! The same type of data at different stations should be different data types

ssqe(:) = 0d0 
  ! Normalize all the data (including both obs. and mod.) before calculating ssqe 
  !(to give equal weight to each data point)
kk(:) = 0
   k  = 0
DO j = 1, Nstn
   Do i = 1, NDTYPE
     allocate(OBS_data(NDPTS(i,j)))
     allocate(MOD_data(NDPTS(i,j)))
   
     if (i .eq. itNO3) then   ! Nitrate
       OBS_data = TINData((kk(i)+1):(kk(i)+nrow(1,j)),3)
       MOD_data = TINout( (kk(i)+1):(kk(i)+nrow(1,j)),1)
     elseif (i .eq. itCHL) then  ! CHL
       OBS_data = CHLData((kk(i)+1):(kk(i)+nrow(2,j)),3)
       MOD_data = CHLout((kk(i)+1):(kk(i)+nrow(2,j)),1)
     elseif (i .eq. itNPP) then  ! NPP
       OBS_data = NPPData((kk(i)+1):(kk(i)+nrow(3,j)),3)
       MOD_data = NPPout((kk(i)+1):(kk(i)+nrow(3,j)),1)
     elseif (i .eq. itPON) then  ! PON
       OBS_data = PONData((kk(i)+1):(kk(i)+nrow(4,j)),3)
       MOD_data =  PONout((kk(i)+1):(kk(i)+nrow(4,j)),1)
     elseif (i .eq. itDFE .and. do_IRON) then ! DFe
       OBS_data = DFeData((kk(i)+1):(kk(i)+nrow(itDFe,j)),3)
       MOD_data =  DFeout((kk(i)+1):(kk(i)+nrow(itDFe,j)),1)
     elseif (i .eq. itPO4 .and. N2fix) then   ! DIP
       OBS_data = PO4Data((kk(i)+1):(kk(i)+nrow(itPO4,j)),3)
       MOD_data =  PO4out((kk(i)+1):(kk(i)+nrow(itPO4,j)),1)
     elseif (i .eq. itPOP .and. N2fix) then   ! POP
       OBS_data = POPData((kk(i)+1):(kk(i)+nrow(itPOP,j)),3)
       MOD_data =  POPout((kk(i)+1):(kk(i)+nrow(itPOP,j)),1)
     elseif (i .eq. itDIA .and. N2fix) then   ! DIA
       OBS_data = DIAData((kk(i)+1):(kk(i)+nrow(itDIA,j)),3)
       MOD_data =  DIAout((kk(i)+1):(kk(i)+nrow(itDIA,j)),1)
     elseif (INCLUDESIZE) then  ! Size-fractionated Chl
       OBS_data = SizeData((kk(i)+1):(kk(i)+nrow(NDTYPE-3,j)), i-(NDTYPE-6))
       MOD_data =  Sizeout((kk(i)+1):(kk(i)+nrow(NDTYPE-3,j)), i-(NDTYPE-4))
     else
       print *, 'Errors in selecting data types! Quit!'
       stop
     endif
   
     if (INCLUDESIZE .and. i .ge. NDTYPE-3) then
        ! The percentages of size-fractions just between 0 and 1
        max_y = 1d0
        min_y = 0d0
     else
        max_y =     maxval(OBS_data,1)
        min_y = max(minval(OBS_data,1),0d0) ! Must be positive
     endif
   
   ! Transform (square root) both model and obs. data and normalize between 0 and 1
     call transform(NDPTS(i,j),OBS_data,min_y,max_y,logtransform,  &
          obsval((k+1):(k+NDPTS(i,j)))  )
   
     call transform(NDPTS(i,j),MOD_data,min_y,max_y,logtransform,  &
          modval((k+1):(k+NDPTS(i,j)))  )

     ! Calculate SSqE:
     do don = (k+1),(k+NDPTS(i,j))
        ssqe(i+(j-1)*NDTYPE) = ssqe(i+(j-1)*NDTYPE) &
                             + (modval(don)-obsval(don))**2 
     enddo
     
     ! Convert back to absolute values
     call transform( NDPTS(i,j),   modval((k+1):(k+NDPTS(i,j))),      &
          min_y,max_y,backlog,     modval((k+1):(k+NDPTS(i,j)))  )

     call transform( NDPTS(i,j),   obsval((k+1):(k+NDPTS(i,j))),      &
          min_y,max_y,backlog,     obsval((k+1):(k+NDPTS(i,j)))  )
  
     k = k + NDPTS(i,j)
     deallocate(OBS_data)
     deallocate(MOD_data)
     kk(i) = kk(i) + NDPTS(i,j)  !Update station index
   Enddo
ENDDO
end subroutine CalcYSSE
