MODULE MOD_1D
USE PARAM_MOD
USE bioCOM_MOD
implicit none

! Grid parameters
real, private, parameter :: hmax   = 250d0 ! Total water depth
real, private, parameter :: thetaS = 2d0   ! surface stretching parameter
real, private, parameter :: zero   = 0d0

!how many seconds of one year
integer, private, parameter :: y_per_s = INT(d_per_s*y_per_d), &

! Number of vertical points of NO3,Temp, and Aks profile
                      N_NO3   = 37,               &  
                      N_PO4   = 37,               &
                      N_Temp  = 57,               &
                      N_w     = 40,               &
                      N_par   = 1 ,               &
                      N_Dust  = 1,                &
                      N_wstr  = 1

integer, private, parameter :: N_PHY = 19  !Why this number? 


integer, private      :: N_Aks(Nstn)              ! Station dependent

! Forcing data time indices 
real, private, target :: obs_time_temp(NFobs(etemp), Nstn)
real, private, target :: obs_time_NO3( NFobs(eNO3) , Nstn)
real, private, target :: obs_time_PO4( NFobs(ePO4) , Nstn)
real, private, target :: obs_time_Aks( NFobs(eAks) , Nstn)
real, private, target :: obs_time_w(   NFobs(ew  ) , Nstn)
real, private, target :: obs_time_par( NFobs(ePAR) , Nstn)
real, private, target :: obs_time_wstr(NFobs(ewstr), Nstn)
real, private, target :: obs_time_Dust(NFobs(eDust), Nstn)
real, private, target :: obs_time_fer( NFobs(eFer) , Nstn)

! Forcing data
real, private :: obs_NO3( N_NO3 ,  1+NFobs(eNO3 ))
real, private :: obs_PO4( N_PO4 ,  1+NFobs(ePO4 ))
real, private :: obs_fer( N_fer ,  1+NFobs(efer ))
real, private :: obs_Temp(N_Temp,  1+NFobs(etemp))
real, private :: obs_PAR( N_par ,  1+NFobs(ePAR ))
real, private :: obs_w(   N_w   ,  1+NFobs(ew   ))
real, private :: obs_Dust(N_Dust,  1+NFobs(eDust))
real, private :: obs_wstr(N_wstr,  1+NFobs(ewstr))
real, private, allocatable :: obs_Aks(:,:) ! (N_Aks ,  1+NFobs(eAks ))

! Initial phytoplankton data:
real, private :: PHY_ini(N_PHY, 2) = 0.1

! Bottom data for NO3, PO4, and fer
real, private, target :: NO3_bot(1,NFobs(eNO3), Nstn)
real, private, target :: PO4_bot(1,NFobs(ePO4), Nstn)
real, private, target :: fer_bot(1,NFobs(eFer), Nstn)
real, private, allocatable :: VarsBom(:,:)

real, pointer         :: pb(:), pc(:,:)
integer, private      :: ncff  ! a scratch integer for the dimensions of matrix

! Vertically interpolated temperature and Aks at each obs timing
real, private, target :: VTemp(   nlev, NFobs(etemp), Nstn)
real, private, target :: VAks(0:nlev,NFobs(eAks ), Nstn)
real, private, target :: Vw(  0:nlev,NFobs(ew )  , Nstn)
real, private, target :: VPAR(1,     NFobs(ePAR) , Nstn)
real, private, target :: VDust(1,    NFobs(eDust), Nstn)
real, private, target :: Vwstr(1,    NFobs(ewstr), Nstn)

logical, public  :: savefile    = .FALSE.
logical, public  :: INCLUDESIZE = .FALSE.

! New calculated state variables
real, private, allocatable    :: ww(:,:)  
!!$$-----------------------------------------------------------------
! The declaration of the data part:
! length of string for labels 
integer, parameter   :: LabelLen = 7 

! NDTYPE is the number of types of data (or observations), 
integer              :: NDTYPE
character(LabelLen), allocatable :: DataLabel(:)

! Number of observations (data points) of each data type
integer, allocatable :: NDPTS(:,:)
integer              :: TNobs(Nstn), ANobs

! Indeces for data type:
integer, parameter   :: itNO3 = 1, itCHL = 2, itNPP = 3
integer, parameter   :: itPON = 4
integer  :: itPO4 = 5,  itPOP = 6, itDIA = 7
integer  :: itP10 = 8
integer  :: itP03 = 9
integer  :: itP01 = 10
integer  :: itP_1 = 11
integer  :: itDFe = 5

! The number of days for model run, needs to change at different stages during MCMC

! The data arrays are read from csv files. 
! Dimensions of the input data
integer, allocatable :: ncol(:,:)
integer, allocatable :: nrow(:,:)
integer, allocatable :: Tnrow(:)

! Declare the dimensions of the data:
real,    allocatable ::  TINData(:,:)
real,    allocatable ::  PO4Data(:,:)
real,    allocatable ::  CHLData(:,:)
real,    allocatable ::  NPPData(:,:)
real,    allocatable ::  PONData(:,:)
real,    allocatable ::  POPData(:,:)
real,    allocatable ::  DIAData(:,:)
real,    allocatable ::  DFeData(:,:)
real,    allocatable :: SizeData(:,:)
real,    allocatable ::  OBSData(:,:)

! The model output (final year) to match with observational data: 
real,    allocatable ::  TINout(:,:) 
real,    allocatable ::  CHLout(:,:) 
real,    allocatable ::  NPPout(:,:) 
real,    allocatable ::  PONout(:,:) 
real,    allocatable ::  PO4out(:,:) 
real,    allocatable ::  POPout(:,:) 
real,    allocatable ::  DFeout(:,:) 
real,    allocatable ::  DIAout(:,:) 
real,    allocatable :: Sizeout(:,:) 

! DOY and Depth for Observational data assembled as a single matrix 
real,    allocatable ::  OBS_DOY(:), OBS_Depth(:)
! Data label for OBS data
character(LabelLen), allocatable :: OBS_Label(:)

! Initial profile of NO3 and PO4:
real, private        :: NO3(nlev,1,Nstn) = 0d0
real, private        :: PO4(nlev,1,Nstn) = 0d0
real, private        :: PH0(nlev,1,Nstn) = 0.1d0

! Initial profile of fer:
real, private        :: fer(nlev,1,Nstn) = 0d0

! interpolation results for particles
real, private, ALLOCATABLE :: zp(:), zout(:,:), dat_in(:,:)

CONTAINS
!-----------------------------------------------------------------------
! Call other subroutines to read and set-up data, as needed 
! This subroutine is called before the real model run
SUBROUTINE SETUP_OBSDATA
implicit none
character(LEN=20)    :: TIN_OBS_file(Nstn)
character(LEN=20)    :: PO4_OBS_file(Nstn)
character(LEN=20)    :: CHL_OBS_file(Nstn)
character(LEN=20)    :: NPP_OBS_file(Nstn)
character(LEN=20)    :: PON_OBS_file(Nstn)
character(LEN=20)    :: POP_OBS_file(Nstn)
character(LEN=20)    :: DFe_OBS_file(Nstn)
character(LEN=20)    :: DIA_OBS_file(Nstn)
character(LEN=20)    :: SIZE_OBS_file(Nstn)
real,    allocatable :: DOY(:), Depth(:)
integer              :: k, i,oi,j
integer              :: NN  ! Number of data types in observational data
integer              :: NL  ! for counting of OBS_DOY et al.
integer, allocatable :: kk(:)

SELECT CASE(MODEL_ID)
  case(GMK98_Lag)
    NDTYPE      = 4  !TIN, CHL, PP, PON
    INCLUDESIZE = .FALSE.
    N2fix       = .FALSE.
  case(NPclosure)
     NDTYPE      = 3  !TIN, CHL, PP
     INCLUDESIZE = .FALSE.
  case(NPZDFix, NPZclosure, NPPZDD,EFTPPDD,Geidersimple,GeiderDroop, EFTsimple, NPZDFixIRON, GeidsimIRON,EFT2sp,NPZD2sp, EFTsimIRON)
     ! Data types must be the same for different stations.
     ! But can be different for different models
     if (trim(Stn(1)) .eq. 'SEATS' .and. Nstn .eq. 1) then
     	NDTYPE   = 3  !TIN, CHL, PP
     else
     	NDTYPE   = 4  !TIN, CHL, PP, PON
     endif

     INCLUDESIZE = .FALSE.
     N2fix       = .FALSE.
  case(NPZDN2)
     NDTYPE      = 7  !TIN, CHL, PP, PON, PO4, POP, DIA
     INCLUDESIZE = .FALSE.
     N2fix       = .TRUE.
 
  case(NPZDdisc, EFTdisc,EFTcont, Geiderdisc,NPZDCONT,CITRATE3)
     if (do_IRON) then
       NDTYPE = 9 !TIN, CHL, PP, PON, DFe, CHL>10, CHL3-10,CHL1-3,CHL<1
     else
       NDTYPE = 8 
     endif
     INCLUDESIZE = .TRUE.
     if (trim(Stn(1))  .eq. 'HOT' .and. Nstn .eq. 1) then
          NDTYPE = 5  !TIN, CHL, PP, PON, DFe
     INCLUDESIZE = .FALSE.
     endif
  case default
     stop 'Model option incorrect! Quit!'
ENDSELECT

ALLOCATE(DataLabel(NDTYPE),  STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating DataLabel ***"
DataLabel(itNO3) = 'TIN'
DataLabel(itCHL) = 'CHL'
DataLabel(itNPP) = 'NPP'
if (NDTYPE .ge. 4) DataLabel(itPON) = 'PON'
if(DO_IRON) then 
   DataLabel(itDFe) = 'DFe'
   itDFe            = itPON + 1       
endif
if (Model_ID .eq. NPZDN2) then
    if (DO_IRON) then
      itPO4          = itDFe + 1
    else
      itPO4          = itPON + 1
    endif
    itPOP            = itPO4 + 1
    itDIA            = itPOP + 1
    DataLabel(itPO4) = 'PO4'
    DataLabel(itPOP) = 'POP'
    DataLabel(itDIA) = 'DIA'
endif

If (INCLUDESIZE .and. (.not. N2fix)) then
  if (DO_IRON) then
    itP10          = itDFe + 1
  else
    itP10          = itPON + 1
  endif
  itP03            = itP10 + 1
  itP01            = itP03 + 1
  itP_1            = itP01 + 1
  DataLabel(itP10) = 'P10'
  DataLabel(itP03) = 'P03'
  DataLabel(itP01) = 'P01'
  DataLabel(itP_1) = 'P_1'
  allocate(nrow(NDTYPE-3,Nstn), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Error in allocating nrow ***"
  allocate(ncol(NDTYPE-3,Nstn), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Error in allocating ncol ***"
Elseif (INCLUDESIZE .and. N2fix) then
  if (DO_IRON) then
    itPO4          = itDFe + 1
  else
    itPO4          = itPON + 1
  endif
  itPOP            = itPO4 + 1
  itP10            = itPOP + 1
  itP03            = itP10 + 1
  itP01            = itP03 + 1
  itP_1            = itP01 + 1
  DataLabel(itPO4) = 'PO4'
  DataLabel(itPOP) = 'POP'
  DataLabel(itP10) = 'P10'
  DataLabel(itP03) = 'P03'
  DataLabel(itP01) = 'P01'
  DataLabel(itP_1) = 'P_1'

  allocate(nrow(NDTYPE-3,Nstn), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Error in allocating nrow ***"
  allocate(ncol(NDTYPE-3,Nstn), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Error in allocating ncol ***"

ELSE
  allocate(nrow(NDTYPE,Nstn), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Error in allocating nrow ***"
  allocate(ncol(NDTYPE,Nstn), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Error in allocating ncol ***"
ENDIF

ALLOCATE(NDPTS(NDTYPE,Nstn),  STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Error in allocating NDPTS ***"

DO j = 1, Nstn
  if(taskid==0) write(6,*) 'Station name: ', trim(Stn(j))
! Assign the data dimension, must be consistent with external file:
  if (trim(Stn(j)) .eq. 'S1') then
      N_Aks(j)       = 40
      NDPTS(itNO3,j) = 1127
      NDPTS(itCHL,j) = 426
      NDPTS(itNPP,j) = 128
      if (NDTYPE .ge. 4) NDPTS(itPON,j) = 32
      if (DO_IRON) NDPTS(itDFe,j) = 168
      if (INCLUDESIZE) then
         NDPTS(itP10,j) = 166
         NDPTS(itP03,j) = NDPTS(itP10,j)
         NDPTS(itP01,j) = NDPTS(itP10,j)
         NDPTS(itP_1,j) = NDPTS(itP10,j)
      endif
  else if (trim(Stn(j)) .eq. 'K2') then
      N_Aks(j)       = 25
      NDPTS(itNO3,j) = 1186
      NDPTS(itCHL,j) = 470
      NDPTS(itNPP,j) = 112
      if (NDTYPE .ge. 4) NDPTS(itPON,j) = 29
      if (DO_IRON) NDPTS(itDFe,j) = 168
      if (INCLUDESIZE) then
         NDPTS(itP10,j) = 143
         NDPTS(itP03,j) = NDPTS(itP10,j)
         NDPTS(itP01,j) = NDPTS(itP10,j)
         NDPTS(itP_1,j) = NDPTS(itP10,j)
      endif
  else if (trim(Stn(j)) .eq. 'HOT') then
      N_Aks(j)       = 41
      NDPTS(itNO3,j) = 6962
      if (Model_ID .eq. NPZDN2) then
         NDPTS(itPO4,j) = 8072
         NDPTS(itPOP,j) = 2850
         NDPTS(itDIA,j) = 23
      endif
      NDPTS(itCHL,j) = 9281
      NDPTS(itNPP,j) = 1904
      if (NDTYPE .ge. 4) NDPTS(itPON,j) = 2852
      if (DO_IRON) NDPTS(itDFe,j) = 168
  else if (trim(Stn(j)) .eq. 'SEATS') then
      N_Aks(j)       = 41
      NDPTS(itNO3,j) = 41
      if (Model_ID .eq. NPZDN2) then
         NDPTS(itPO4,j) = 0
         NDPTS(itPOP,j) = 2850
         NDPTS(itDIA,j) = 23
      endif
      NDPTS(itCHL,j) = 41
      NDPTS(itNPP,j) = 384
      if (NDTYPE .ge. 4) NDPTS(itPON,j) = 2806
      if (DO_IRON) NDPTS(itDFe,j) = 168
  else
      write(6,*) 'Station number incorrect! Stop!'
      stop
  endif
  TNobs(j)  = sum(NDPTS(:,j))
  ANobs     = sum(TNobs(:))   !All data points of all stations

  ! Assign obs. data matrix: 
  IF (INCLUDESIZE) THEN
     ncff = NDTYPE-3
     do i = 1,ncff
        nrow(i,j)=NDPTS(i,j)
        if (i < ncff) then
           ncol(i,j)=3
        else
           ncol(i,j)=6
        endif
     enddo
  ELSE
     do i = 1,NDTYPE
        nrow(i,j)=NDPTS(i,j)
        ncol(i,j)=4
     enddo
  ENDIF
  if (taskid==0) &
  write(6,'(A31,1x,A3,1x,A1,1x,I5)') 'Total number of observations of',Stn(j),'=',TNobs(j)
Enddo  ! ==> End of assigning data numbers at different stations
  
! Total number of observations of each data type for all stations
NN = size(nrow,1)
ALLOCATE(Tnrow(NN))
DO i = 1, NN
   Tnrow(i)=sum(nrow(i,:))
ENDDO

!Start to assign TIN,Chl, and NPP data:
allocate( TINData(Tnrow(itNO3),ncol(itNO3,1)), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating TINData ***"

allocate( CHLData(Tnrow(itCHL),ncol(itCHL,1)), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating CHLData ***"

allocate( NPPData(Tnrow(itNPP),ncol(itNPP,1)), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating NPPData ***"

if (NDTYPE .ge. 4) then
   allocate( PONData(Tnrow(itPON),ncol(itPON,1)), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating PONData ***"
endif

if (do_IRON) then
   allocate( DFeData(Tnrow(itDFe),ncol(itDFe,1)), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating DFeData ***"
endif

if (N2fix) then
   allocate( PO4Data(Tnrow(itPO4),ncol(itPO4,1)), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating PO4Data ***"
   allocate( POPData(Tnrow(itPOP),ncol(itPOP,1)), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating POPData ***"
   allocate( DIAData(Tnrow(itDIA),ncol(itDIA,1)), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating DIAData ***"
endif

allocate( OBSData(ANobs  ,3 ), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating OBSData ***"

! Initialize observational data:
TINData(:,:) = 0d0
CHLData(:,:) = 0d0
NPPData(:,:) = 0d0
if (NDTYPE .ge. 4) PONData(:,:) = 0d0

if (N2fix) then
    PO4Data(:,:) = 0.
    POPData(:,:) = 0.
    DIAData(:,:) = 0.
endif

if (do_IRON) DFeData(:,:) = 0d0
OBSData(:,:) = 0d0

! Initialize size data:
if (INCLUDESIZE) then
   ncff = NDTYPE - 3
   allocate(SizeData(Tnrow(ncff),ncol(ncff,1)), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating SizeData ***"
   SizeData(:,:) = 0d0
endif

! Initialize observed DOY and Depth
! This is for writing out output
allocate(OBS_DOY(  ANobs))
allocate(OBS_Depth(ANobs))
allocate(OBS_Label(ANobs))

OBS_DOY(:)   = 0d0
OBS_Depth(:) = 0d0
OBS_Label(:) = ' NA    '

allocate(kk(NN))
kk(:) = 0
NL    = 0
Do j = 1, Nstn
  ! Read observational data:
  TIN_OBS_file(j) = trim(Stn(j))//'_DIN.dat' 
  PO4_OBS_file(j) = trim(Stn(j))//'_DIP.dat' 
  POP_OBS_file(j) = trim(Stn(j))//'_POP.dat' 
  CHL_OBS_file(j) = trim(Stn(j))//'_CHL.dat'
  NPP_OBS_file(j) = trim(Stn(j))//'_PP.dat'
  PON_OBS_file(j) = trim(Stn(j))//'_PON.dat'
  DFe_OBS_file(j) = trim(Stn(j))//'_DFe.dat'
  DIA_OBS_file(j) = trim(Stn(j))//'_DIA.dat'
  
  !The first column of the data file must be DOY
  !The second column must be depth
  !The third column must be the real data
  call Readcsv(TIN_OBS_file(j),nrow(1,j),ncol(1,j), TINData((kk(1)+1):(kk(1)+nrow(1,j)),:))
  call Readcsv(CHL_OBS_file(j),nrow(2,j),ncol(2,j), CHLData((kk(2)+1):(kk(2)+nrow(2,j)),:))
  call Readcsv(NPP_OBS_file(j),nrow(3,j),ncol(3,j), NPPData((kk(3)+1):(kk(3)+nrow(3,j)),:))

  if (NDTYPE .ge. 4) & 
  call Readcsv(PON_OBS_file(j),nrow(4,j),ncol(4,j), &
               PONData((kk(4)+1):(kk(4)+nrow(4,j)),:))

  if (N2fix) then
     call Readcsv(PO4_OBS_file(j),nrow(itPO4,j),ncol(itPO4,j), PO4Data((kk(itPO4)+1):(kk(itPO4)+nrow(itPO4,j)),:))
     call Readcsv(POP_OBS_file(j),nrow(itPOP,j),ncol(itPOP,j), POPData((kk(itPOP)+1):(kk(itPOP)+nrow(itPOP,j)),:))
     call Readcsv(DIA_OBS_file(j),nrow(itDIA,j),ncol(itDIA,j), DIAData((kk(itDIA)+1):(kk(itDIA)+nrow(itDIA,j)),:))
  endif

  if (do_IRON) &
  call Readcsv(DFe_OBS_file(j),nrow(itDFe,j),ncol(itDFe,j), DFeData((kk(itDFe)+1):(kk(itDFe)+nrow(itDFe,j)),:))

  if (INCLUDESIZE) then
     ! Use percentage data
     SIZE_OBS_file(j) = trim(Stn(j))//'_size_Perc.dat'
     call Readcsv(SIZE_OBS_file(j),nrow(ncff,j),ncol(ncff,j),SizeData((kk(ncff)+1):(kk(ncff)+nrow(ncff,j)),:))
  endif

  ! Give DOY and Depth for all the obs. data
  k = 0
  do i = 1, NDTYPE
    allocate(  DOY(NDPTS(i,j)))
    allocate(Depth(NDPTS(i,j)))

    if (i .eq. itNO3) then   ! Nitrate
      DOY   = TINData((kk(1)+1):(kk(1)+nrow(1,j)),1)
      Depth = TINData((kk(1)+1):(kk(1)+nrow(1,j)),2)
    else if (i .eq. itCHL) then   ! CHL
      DOY   = CHLData((kk(2)+1):(kk(2)+nrow(2,j)),1)
      Depth = CHLData((kk(2)+1):(kk(2)+nrow(2,j)),2)
    else if (i .eq. itNPP) then   ! NPP
      DOY   = NPPData((kk(3)+1):(kk(3)+nrow(3,j)),1)
      Depth = NPPData((kk(3)+1):(kk(3)+nrow(3,j)),2)
    else if (i .eq. itPON) then  ! PON
      DOY   = PONData((kk(4)+1):(kk(4)+nrow(4,j)),1)
      Depth = PONData((kk(4)+1):(kk(4)+nrow(4,j)),2)
    else
      if (do_IRON) then
        if (i .eq. itDFe) then  !DFe
           DOY   = DFeData((kk(itDFe)+1):(kk(itDFe)+nrow(itDFe,j)),1)
           Depth = DFeData((kk(itDFe)+1):(kk(itDFe)+nrow(itDFe,j)),2)
        elseif (N2fix) then
          if (i .eq. itPO4) then   ! PO4
             DOY   = PO4Data((kk(itPO4)+1):(kk(itPO4)+nrow(itPO4,j)),1)
             Depth = PO4Data((kk(itPO4)+1):(kk(itPO4)+nrow(itPO4,j)),2)
          else if (i .eq. itPOP) then   ! POP
             DOY   = POPData((kk(itPOP)+1):(kk(itPOP)+nrow(itPOP,j)),1)
             Depth = POPData((kk(itPOP)+1):(kk(itPOP)+nrow(itPOP,j)),2)
          else if (i .eq. itDIA) then   ! DIA
             DOY   = DIAData((kk(itDIA)+1):(kk(itDIA)+nrow(itDIA,j)),1)
             Depth = DIAData((kk(itDIA)+1):(kk(itDIA)+nrow(itDIA,j)),2)
          else   ! Size-fractionated Chl
             DOY   = SizeData((kk(ncff)+1):(kk(ncff)+nrow(ncff,j)),1)
             Depth = SizeData((kk(ncff)+1):(kk(ncff)+nrow(ncff,j)),2)
          endif
        else  ! Size-fractionated Chl
            DOY   = SizeData((kk(ncff)+1):(kk(ncff)+nrow(ncff,j)),1)
            Depth = SizeData((kk(ncff)+1):(kk(ncff)+nrow(ncff,j)),2)
        endif ! <== if (i .eq. itDFe)
      else ! .not. do_IRON
        if (N2fix) then
          if (i .eq. itPO4) then   ! PO4
             DOY   = PO4Data((kk(itPO4)+1):(kk(itPO4)+nrow(itPO4,j)),1)
             Depth = PO4Data((kk(itPO4)+1):(kk(itPO4)+nrow(itPO4,j)),2)
          else if (i .eq. itPOP) then   ! POP
             DOY   = POPData((kk(itPOP)+1):(kk(itPOP)+nrow(itPOP,j)),1)
             Depth = POPData((kk(itPOP)+1):(kk(itPOP)+nrow(itPOP,j)),2)
          else if (i .eq. itDIA) then   ! DIA
             DOY   = DIAData((kk(itDIA)+1):(kk(itDIA)+nrow(itDIA,j)),1)
             Depth = DIAData((kk(itDIA)+1):(kk(itDIA)+nrow(itDIA,j)),2)
          else   ! Size-fractionated Chl
             DOY   = SizeData((kk(ncff)+1):(kk(ncff)+nrow(ncff,j)),1)
             Depth = SizeData((kk(ncff)+1):(kk(ncff)+nrow(ncff,j)),2)
          endif
        else  ! Size-fractionated Chl
           DOY   = SizeData((kk(ncff)+1):(kk(ncff)+nrow(ncff,j)),1)
           Depth = SizeData((kk(ncff)+1):(kk(ncff)+nrow(ncff,j)),2)
        endif
      endif
    endif

    do oi = 1, NDPTS(i,j)
         OBS_DOY(NL+k+oi) =   DOY(oi)
       OBS_Depth(NL+k+oi) = Depth(oi)
       OBS_Label(NL+k+oi) = DataLabel(i)
    enddo
    k = k + NDPTS(i,j)
    DEALLOCATE(DOY)
    DEALLOCATE(Depth)
  ENDDO

  !Update kk (the index in TINData etc.)
  do oi = 1, NN
     kk(oi) = kk(oi) + nrow(oi,j)
  enddo
  !Update NL (for OBS_DOY et al.)
  NL = NL + k
  if (taskid==0) write(6,'(A27,1x,A5)') 'Observation data set up for', Stn(j)
Enddo
END SUBROUTINE Setup_OBSdata

!========================================================
SUBROUTINE MODEL_SETUP
! Called at the beginning of the program, only once
implicit none
integer            :: i,j,k
real               :: a(1)  !a scratch vector with length 1
! Input files
character(LEN=20)  :: forcfile(TNFo), forcfiletime(TNFo)
character(LEN=20)  :: PHYfile

logical            :: there  !Check whether file exists

!the fraction of a time step in one day
dtdays = dtsec/d_per_s

if (taskid==0) &
write(6,'(A13,1x,F12.3,A12)') 'Timestepping: ',dtdays,'of one day.'

! Setup grid
call setup_grid

do j = 1, Nstn
   ! Forcing files:
   do i = 1, TNFo
      forcfile(i)     = trim(Stn(j))//'_'//trim(LabelForc(i))//'.dat'
      forcfiletime(i) = trim(Stn(j))//'_'//trim(LabelForc(i))//'_time.dat'
   enddo
   
   ! Use Taketo's Aks data only for K2:
   if (trim(Stn(j)) .eq. 'K2') then
      forcfile(eAks)     = trim(Stn(j))//'_'//trim(LabelForc(eAks))//'T.dat'
      forcfiletime(eAks) = trim(Stn(j))//'_'//trim(LabelForc(eAks))//'_timeT.dat'
   endif
   
   ! Read NO3 data:
   call Readcsv(forcfile(eNO3), N_NO3, size(obs_NO3,2), obs_NO3) 
   
   ! Interpolate initial NO3:
   call gridinterpol(N_NO3,1,obs_NO3(:,1),obs_NO3(:,2),                 &
                     nlev, Z_r, NO3(:,1,j)) 

   ! Get bottom data for NO3:
   ! Interpolate NO3 to the bottom depth:
   a(1) = Z_w(0)
   call gridinterpol(N_NO3,NFobs(eNO3),obs_NO3(:,1),               &
                     obs_NO3(:,2:(NFobs(eNO3)+1)),                 &
                     1, a, NO3_bot(1,:,j)) 

   ! Read initial PHYTO data:
   PHYfile = trim(Stn(j))//'_'//'PHY.dat'

   !Inquire whether PHYfile exists
   INQUIRE( FILE=PHYfile, EXIST=THERE )
   if (there) then
      call Readcsv(PHYfile, N_PHY, 2, PHY_ini) 

      ! Interpolate initial PHY:
      call gridinterpol(N_PHY,1,PHY_ini(:,1),PHY_ini(:,2),                 &
                        nlev, Z_r, PH0(:,1,j)) 

   endif

   if (N2fix) then
      ! Read PO4 data:
      call Readcsv(forcfile(ePO4), N_PO4, size(obs_PO4,2), obs_PO4) 
      
      ! Interpolate initial PO4:
      call gridinterpol(N_PO4,1,obs_PO4(:,1),obs_PO4(:,2),                 &
                        nlev, Z_r, PO4(:,1,j)) 

      ! Get bottom data for PO4:
      call gridinterpol(N_PO4,NFobs(ePO4),obs_PO4(:,1),               &
                     obs_PO4(:,2:(NFobs(ePO4)+1)),                    &
                     1, a, PO4_bot(1,:,j)) 
   endif

   if (do_IRON) then
     ! Read fer data:
     call Readcsv(forcfile(efer), N_fer, size(obs_fer,2), obs_fer) 
     
     ! Interpolate initial fer:
     call gridinterpol(N_fer,1,obs_fer(:,1),obs_fer(:,2),                 &
                       nlev, Z_r, fer(:,1,j)) 

     ! Get bottom data for fer:
     call gridinterpol(N_fer,NFobs(efer),obs_fer(:,1),               &
                    obs_fer(:,2:(NFobs(efer)+1)),                    &
                    1, a, fer_bot(1,:,j)) 

   endif
  
   ! Calculate obs. time indices in seconds:
   ! Read obs. time file:
   call Readcsv(forcfiletime(ePAR), 1,NFobs(ePAR), obs_time_par(:,j)) 

   ! Convert obs_time to the unit of seconds:
   obs_time_par(:,j)  = obs_time_par(:,j) *3d1*dble(d_per_s)

   call Readcsv(forcfiletime(eNO3), 1,NFobs(eNO3), obs_time_NO3(:,j)) 
   obs_time_NO3(:,j)  = obs_time_NO3(:,j) *3d1*dble(d_per_s)

   call Readcsv(forcfiletime(etemp),1,NFobs(etemp),obs_time_temp(:,j)) 
   obs_time_temp(:,j) = obs_time_temp(:,j)*3d1*dble(d_per_s)

   call Readcsv(forcfiletime(eAks), 1,NFobs(eAks), obs_time_Aks(:,j)) 
   obs_time_Aks(:,j)  = obs_time_Aks(:,j) *3d1*dble(d_per_s)

   !call Readcsv(forcfiletime(ew  ), 1,NFobs(ew)  , obs_time_w(:,j)) 
   !obs_time_w(:,j)    = obs_time_w(:,j)   *3d1*dble(d_per_s)
   
   if (N2fix) then
      call Readcsv(forcfiletime(ePO4), 1,NFobs(ePO4), obs_time_PO4(:,j)) 
      obs_time_PO4(:,j)  = obs_time_PO4(:,j) *3d1*dble(d_per_s)

      ! Read wind stress data:
      call Readcsv(forcfiletime(ewstr), 1,NFobs(ewstr), obs_time_wstr(:,j)) 
      obs_time_wstr(:,j)  = obs_time_wstr(:,j) *3d1*dble(d_per_s)
   endif

   if (do_IRON) then  ! Read dust deposition time data:
      call Readcsv(forcfiletime(eDust), 1,NFobs(eDust), obs_time_Dust(:,j)) 
      obs_time_Dust(:,j) = obs_time_Dust(:,j) *3d1*dble(d_per_s)  ! Dust deposition

      if (bot_bound .eq. Dirichlet) then  ! Read iron time
         call Readcsv(forcfiletime(eFer), 1,NFobs(eFer), obs_time_fer(:,j)) 
         obs_time_fer(:,j) = obs_time_fer(:,j) *3d1*dble(d_per_s)  ! Bottom fer
      endif
   endif

   ! Read external w data:
   !call Readcsv(forcfile(ew), size(obs_w,1), size(obs_w,2), obs_w) 
   
   ! Interpolate external w data:
   
   !subroutine gridinterpol(N,cols,obs_z,obs_prof,nlev_,model_z,model_prof)
   !call gridinterpol(N_w,NFobs(ew),obs_w(:,1),       &
   !      obs_w(:,2:(NFobs(ew)+1)),                   &
   !      nlev+1, Z_w, Vw(:,:)) 
   
   ! Read external PAR data:
   call Readcsv(forcfile(ePAR), size(obs_PAR,1),size(obs_PAR,2),obs_PAR) 

   ! No need to vertically interpolate PAR data:
   do k = 1, (size(obs_PAR,2)-1)
      VPAR(1,k,j) = obs_PAR(1,k+1)
   enddo

   If (do_IRON) then
      ! Read external Dust data:
      call Readcsv(forcfile(eDust),size(obs_Dust,1),size(obs_Dust,2),obs_Dust) 

      ! No need to vertically interpolate Dust data:
      do k = 1, (size(obs_Dust,2)-1)
         VDust(1,k,j) = obs_Dust(1,k+1)
      enddo
   Endif

   If (N2fix) then
      ! Read external wind stress data:
      call Readcsv(forcfile(ewstr),size(obs_wstr,1),size(obs_wstr,2),obs_wstr) 

      ! No need to vertically interpolate wstr data:
      do k = 1, (size(obs_wstr,2)-1)
         Vwstr(1,k,j) = obs_wstr(1,k+1)
      enddo
   Endif

   ! Read external Temp data:
   call Readcsv(forcfile(etemp),size(obs_Temp,1),size(obs_Temp,2),obs_Temp) 
   
   ! Interpolate external Temp data:
   call gridinterpol(N_Temp, NFobs(etemp),obs_Temp(:,1),                &
         obs_Temp(:,2:(NFobs(etemp)+1)),                                &
         nlev, Z_r, VTemp(:,:,j)) 

   ! Read external Aks data:
   allocate(obs_Aks(N_Aks(j), 1+NFobs(eAks)),  STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating obs_Aks ***"
   call Readcsv(forcfile(eAks), N_Aks(j), size(obs_Aks,2), obs_Aks) 
   
   ! Interpolate external Aks data:
   call gridinterpol(N_Aks(j),NFobs(eAks),obs_Aks(:,1),                 &
         obs_Aks(:,2:(NFobs(eAks)+1)),                                  &
         nlev+1, Z_w, VAks(:,:,j)) 
   
   savefile = .FALSE.
   if(taskid==0) write(6,*) 'Model setup finished for ', Stn(j)
   deallocate(obs_Aks)
Enddo
END SUBROUTINE Model_setup
!========================================================
! This subroutine does the main heavy lifting work in the 1D model
! and give the output to match with the observational data
! Must be called after initialization
SUBROUTINE Timestep
implicit none
real,    parameter  :: cnpar      = 0.6
real,    parameter  :: Taur(nlev) = 1D12  !Relaxation time
real,    parameter  :: Vec0(nlev) = 0d0   !Vectors of zero

! Local scratch variables
integer  :: i,k,nm,m,j,jj, Nstep, current_DOY,DOY
integer  :: current_year = 0
integer  :: current_hour = 0
integer  :: j_ = 0   !A scratch variable
integer  :: NP = 0   !The number of particles
logical  :: match = .FALSE.
INTEGER,         ALLOCATABLE :: add_death(:)
INTEGER,         ALLOCATABLE ::    ideath(:)
Type (Particle), ALLOCATABLE ::        Zd(:)

! Counting of the index of each data type
real     :: cff,current_sec
real     :: I_0(1), dust0(1), Aks(0:nlev),w(0:nlev)
real     :: depth(1)
real     :: Vars1(nlev),Vars2(nlev),ww_(0:nlev)
real, allocatable  :: a(:,:)  
integer, parameter :: mode0   = 0
integer, parameter :: mode1   = 1
real,    parameter :: mon2sec = 2592D3

character(LEN=20)  :: outfile, par_file
integer, PARAMETER :: parunit = 35

!End of declaration

! For showing whether some unconstrained outputs are realistic or not
if (.not. allocated(TINout)) allocate(TINout(size(TINData,1),1),STAT = AllocateStatus)
if (AllocateStatus /= 0) STOP "Problem in allocating TINout!"
TINout(:,:) = 0d0

if (.not. allocated(CHLout)) allocate(CHLout(size(CHLData,1),1),STAT = AllocateStatus)
if (AllocateStatus /= 0) STOP "Problem in allocating CHLout!"
CHLout(:,:) = 0d0

if (.not. allocated(NPPout)) allocate(NPPout(size(NPPData,1),1),STAT = AllocateStatus)
if (AllocateStatus /= 0) STOP "Problem in allocating NPPout!"
NPPout(:,:) = 0d0

If (NDTYPE .ge. 4) then
  if (.not. allocated(PONout)) allocate(PONout(size(PONData,1),1),STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Problem in allocating PONout!"
  PONout(:,:) = 0d0
Endif

if (do_IRON) then
  if (.not. allocated(DFeout)) allocate(DFeout(size(DFeData,1),1),STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "Problem in allocating DFeout!"
  DFeout(:,:) = 0.
endif

if (N2fix) then
   if (.not. allocated(PO4out)) allocate(PO4out(size(PO4Data,1),1),STAT = AllocateStatus)
   if (AllocateStatus /= 0) STOP "Problem in allocating PO4out!"
   PO4out(:,:) = 0.
   if (.not. allocated(POPout)) allocate(POPout(size(POPData,1),1),STAT = AllocateStatus)
   if (AllocateStatus /= 0) STOP "Problem in allocating POPout!"
   POPout(:,:) = 0.
   if (.not. allocated(DIAout)) allocate(DIAout(size(DIAData,1),1),STAT = AllocateStatus)
   if (AllocateStatus /= 0) STOP "Problem in allocating DIAout!"
   DIAout(:,:) = 0.
endif
if (INCLUDESIZE) then
   if (.not. allocated(Sizeout)) allocate(Sizeout(size(SizeData,1),4),STAT = AllocateStatus)
   if (AllocateStatus /= 0) STOP "Problem in allocating Sizeout!"
   Sizeout(:,:) = 0d0
endif

DO jj = 1, Nstn

  ! Initialize initial NO3:
  Vars(iNO3,:) = NO3(:,1,jj)

   !For lagrangian model, initialize individuals
  IF (Model_ID .eq. GMK98_Lag) THEN
     allocate(p_PHY(N_ini),  STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Problem in allocating p_PHY ***"

     allocate(AllIDs(N_ini), stat = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Problem in allocating AllIDs ***"

     do k = 1, N_ini
        p_PHY(k)%ID = k        

        !Compute rx (randomly distributed between Z_w(0) and Z_w(nlev)

        CALL RANDOM_NUMBER(cff)
        p_PHY(k)%rz = cff *(Z_w(nlev)-Z_w(0)) + Z_w(0)

        ! Find iz
        DO j_ = 1, nlev
           IF (Z_w(j_-1) .le. p_PHY(k)%rz .AND. Z_w(j_) .gt. p_PHY(k)%rz) THEN
              p_PHY(k)%iz = j_
              EXIT            
           ENDIF
        ENDDO

        AllIDs(k)   = k
     enddo
  ENDIF 

  ! Initialize other variables in the eulerian model:
  DO k = 1, nlev

     if (NVAR > iDET) then
        do i = (iDET+1), NVAR
           Vars(i,k) = 1D-2
        enddo
     endif

     if (Model_ID .ne. GMK98_Lag) then
        do i = 1,NPHY
           ! Initialize initial PHYTO:
           Vars(iPHY(i), k)     = PH0(k,1,jj)/float(NPHY)
           if (Model_ID==GeiderDroop) then
              Vars(iPHYC(i), k) = Vars(iPHY(i),k) * 106./16.
              Vars(iCHL(i),  k) = Vars(iPHYC(i),k)*  12./50.
           endif
        enddo
     endif

     Vars(iZOO,k) = 0.1

     !Following Verity et al. AME (1996)
     if (iDET  > 0) Vars(iDET,k) = Vars(iZOO,k) 
     if (iZOO2 > 0) Vars(iZOO2,k)= Vars(iZOO,k) * .5

  ENDDO

  if (do_IRON) then
    ! Initialize initial fer:
     Vars(ifer,:) = fer(:,1,jj)
     do k = 1, nlev
        Vars(iDETFe,k) = Vars(iDET,k)*Fe_N
     enddo
  endif

   ! Initialize initial PO4:
  if (iPO4 > 0) Vars(iPO4,:)=PO4(:,1,jj)

  ! Initialize diazotrophs:
  if (Model_ID .eq. NPZDN2) then
      do k=1,nlev
        do i = 1,NPHY
           Vars(iPHY(i), k) = 0.01/16d0/float(NPHY)
        enddo
        Vars(iZOO,k) =0.01/16d0
        Vars(iDIA, k)=1D-6/16d0
        Vars(iDETp,k)=1D-6
      enddo
  endif

  ! Initialize output data:
  Varout(:,:) = 0D0
  do i = 1,NVAR
     Varout(i,:)= Vars(i,:)
  enddo

  if (bot_bound .eq. Dirichlet) then
     ! Initialize Bottom values of Vars:
     allocate(VarsBom(1,NVAR),  STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Problem in allocating VarsBom ***"
     VarsBom(:,:)        = 1D-10
     VarsBom(1, iDET)    = 3D-2   !Based on PON at 500 m at HOT

     if (Model_ID .eq. NPZDN2 .and. trim(Stn(jj)) .eq. 'HOT') then 
        ! Initialize bottom values of PHY, ZOO, DET based on PON and POP at HOT:
        VarsBom(1, iZOO)    = 2D-3/3d0 ! Unit: uM P   
        VarsBom(1, iPHY(1)) = 2D-3/3d0 ! Unit: uM P
        VarsBom(1, iDETp)   = 2D-3/3d0 ! Unit: uM P
     endif

     if (iDETFe > 0) then
        VarsBom(1, iDETFe)  = VarsBom(1,iDET)*Fe_N  !Unit: nM
     endif
  endif

  ! Sinking rate
  allocate(ww(0:nlev,NVsinkterms),  STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Problem in allocating ww ***"

  ! Initialize sinking rate (UNIT: m/s !):
  ww(:,:) = 0d0
  ncff    = size(ww,2)
  cff     = exp(params(iwDET))
  do k = 0,nlev-1

   !Phytoplankton no sinking 
   !Detritus sinking rate (convert to UNIT: m/s)
    if (Model_ID==NPPZDD .OR. Model_ID==EFTPPDD) then
       ww(k,ncff-1)= -cff/dble(d_per_s) 
       ww(k,ncff)  = -exp(Params(iwDET2))/dble(d_per_s) 
    elseif (Model_ID==NPZDN2 .or. Model_ID==NPZDcont .or. &
            Model_ID==CITRATE3) then
       ww(k,ncff-1)=-cff/dble(d_per_s) 
       ww(k,ncff)  =ww(k,ncff-1)
    else
       ww(k,ncff) = -cff/dble(d_per_s) 
    endif
  enddo

  ! Create output files:
  if (savefile .and. taskid == 0) then

    ! Merge all data into one output file (Oct. 15 2016)
     outfile = trim(Stn(jj))//'.out'

    ! Create data out file
    open (unit=10, file = outfile, status = 'replace')
    write(10, 1) 'Type','Timestep','Day', 'Hour', Z_r
  endif
  
1 format(4(A10,1x), <nlev  >(2x,F12.5))
3 format(  A10,1x,  <nlev+1>(2x,F12.5))
  
  ! Number of time steps
  Nstep = NDays*INT(d_per_s)/INT(dtsec) 

  ! 'START TIME STEPPING'
  DO it = 1, Nstep+1
  
  ! Calculate current timing (zero is starting time):
     current_sec = float(it-1)*dtsec

  !For each time step, read in external environmental data%%%%%%%%%%%%%%%

    ! Interpolate Temp data throughout the water column:
     pb => obs_time_temp(:,jj)
     ncff = size(pb,1)
     pc => VTemp(:,:,jj)
     call time_interp(int(current_sec),ncff, nlev,pb,pc,Temp)

    ! Interpolate temporal PAR data: 
    ! pb => obs_time_par(:,jj)
    ! ncff = size(pb,1)
    ! pc => VPAR(:,:,jj)
    ! call time_interp(int(current_sec),ncff,1,pb,pc,I_0)
    
    ! Convert the unit of par to W m-2:
    ! cff = I_0(1)/0.4

    ! Calculate diurnal changes of light
    ! Calculate time of the day
     sec_of_day = mod(nint(current_sec), nint(d_per_s))

    ! Calculate DOY
     ! Calculate model time in days:
     current_day = int(current_sec/d_per_s)
   
     ! Calculate current year
     current_year= int(current_day/y_per_d)
 
     ! Calculate DATE OF the YEAR (DOY)
     current_DOY = mod(current_day, y_per_d)

     ! Calculate current hour
     current_hour = int(current_sec/s_per_h)
     current_hour = mod(current_hour, h_per_d)

     !Calculate surface light
     cff = diurnal_light(current_DOY, Stn_lat(jj), sec_of_day)

    ! Calculate light field:
     Vars1(:) = Varout(oCHLt,:)
     call Calculate_PAR(cff, nlev, Hz, Vars1, PAR, PAR_w)

     if (bot_bound .eq. Dirichlet) then
        !  Interpolate bottom NO3 data:
        pb => obs_time_NO3(:,jj)
        pc => NO3_bot(:,:,jj)
        call time_interp(int(current_sec),NFobs(eNO3),1,pb,pc,VarsBom(:,iNO3))
     endif

     IF (DO_IRON) THEN
        ! Interpolate temporal Dust data: 
        pb => obs_time_Dust(:,jj)
        ncff = size(pb,1)
        pc => VDust(:,:,jj)
        call time_interp(int(current_sec),ncff,1,pb,pc, dust0)
  
     ! Iron atmospheric deposition:
     ! Soluble iron deposition unit: kg/m2/s.
     ! so need to convert into nM at each time step
     ! Deposition = Dust*10^12/56*dtsec*surface_area/surface_grid_volume
       cff= dust0(1)*1D9/55.85/Hz(nlev)*dtsec

       ! added dissolved Fe (nM/d) on top grid: 
       Varout(odstdep,nlev)=cff/dtsec*d_per_s
       
       Vars(ifer,nlev) = Vars(ifer,nlev) + cff

       if (bot_bound .eq. Dirichlet) then
          !  Interpolate bottom fer data:
          pb => obs_time_fer(:,jj)
          pc => fer_bot(:,:,jj)
          call time_interp(int(current_sec),NFobs(eFer),1,pb,pc,VarsBom(:,ifer))
       endif
     ENDIF

     if (N2fix) then
        ! Interpolate temporal wstr data
        pb => obs_time_wstr(:,jj)
        ncff = size(pb,1)
        pc => Vwstr(:,:,jj)
        call time_interp(int(current_sec),ncff,1,pb,pc, wstr0)

        if (bot_bound .eq. Dirichlet) then
           ! Interpolate bottom PO4 data:
           pb => obs_time_PO4(:,jj)
           pc => PO4_bot(:,:,jj)
           call time_interp(int(current_sec),NFobs(ePO4),1,pb,pc,VarsBom(:,iPO4))
        endif
     endif

    ! Interpolate Aks data throughout the water column:
     pb => obs_time_Aks(:,jj)
     ncff = size(pb,1)
     pc => VAks(:,:,jj)
     call time_interp(int(current_sec), ncff, nlev+1,pb, pc, Aks)
     nullify(pb)
     nullify(pc) 

    ! Update PAR, Temp, NO3, and Aks associated with each particle for Lagrangian model
     if (Model_ID .eq. GMK98_Lag) then 
        allocate(zp(1), stat=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Problem in allocating zp ***"
        zp(1) = 0.
        allocate(zout(1,1), stat=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Problem in allocating zout ***"
        zout(1,1) = 0.
        
        allocate(dat_in(0:nlev, 1), stat=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Problem in allocating dat_in ***"
        dat_in(:,1) = PAR_w(:)  !Observed profile of PAR

        do k = 1, size(p_PHY) 
           zp(1) = p_PHY(k)%rz
           call gridinterpol(nlev+1, 1, Z_w, dat_in, 1, zp, zout)
           p_PHY(k)%PAR = zout(1,1)

        enddo

        DEALLOCATE(dat_in)

        allocate(dat_in(nlev, 1), stat=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Problem in allocating dat_in ***"
        dat_in(:,1) = Temp(:)  !Observed profile of temperature
        do k = 1, size(p_PHY) 
           zp(1) = p_PHY(k)%rz
           call gridinterpol(nlev, 1, Z_r, dat_in, 1, zp, zout)
           p_PHY(k)%temp = zout(1,1)
        enddo

        dat_in(:,1) = Vars(iNO3,:)  !Observed profile of NO3
        do k = 1, size(p_PHY) 
           zp(1) = p_PHY(k)%rz
           call gridinterpol(nlev, 1, Z_r, dat_in, 1, zp, zout)
           p_PHY(k)%NO3 = zout(1,1)
        enddo

        DEALLOCATE(dat_in)
        DEALLOCATE(zout)
        DEALLOCATE(zp)
     endif

    ! Biological rhs:
     call BIOLOGY

 
    ! Interpolate w data throughout the water column:
    ! call time_interp(int(current_sec), size(obs_time_w,1), nlev+1,&
    !      obs_time_w, Vw, w)
    w(0   ) = 0d0
    w(nlev) = 0d0
    w(:)    = 0d0

    IF (mod(it, nsave) .EQ. 1) THEN

      ! Check whether the values are valid:
      do j = 1,NVAR
         do k = 1, nlev
            if( Vars(j,k) < 0d0 .or. Vars(j,k) .ne. Vars(j,k)) then
                write(6,*) 'At Year ', current_year
                write(6,*) 'At day ',  current_day
                write(6,*) 'At hour ', current_hour
                write(6,*) 'WARNING! The variable ',trim(Labelout(j+ow)), &
                          ' is ', Vars(j,k), ' at depth ',Z_r(k)
                write(6,*) 'Stop this model run and go to next iteration!'
                return

            endif
         enddo 
      enddo
  
      ! Save data to output files:
      if (savefile .and. taskid == 0) then

      ! Save Temp data into the output file:
         write(10, 200) trim(Labelout(oTemp)), it, current_day, current_hour, Temp

      ! Save PAR data into the output file:
         write(10, 200) trim(Labelout(oPAR)),  it, current_day, current_hour, PAR_int
  
      ! Save Aks data into the output file:
         write(10, 200) trim(Labelout(oAks)),  it, current_day, current_hour, Aks(1:nlev)

      ! Save w data into the output file:
      !   write(9+ow  , 4) it, current_day, w
  
      ! Save state variables and diagnostics:
         do i=(ow+1),(Nout+ow)
            Vars1(:) = Varout(i-ow,:)
            write(10, 200) trim(Labelout(i)),  it, current_day, current_hour, Vars1
         enddo
200   format(A10,1x,I10,1x,I7,1x, I7, <nlev>(2x, 1pe12.4))

      ! Save particles for Lagrangian model (save every year)
         if (Model_ID .eq. GMK98_Lag) then
                        ! If it is the beginning of a new day, create a new file
            if (current_hour .eq. 0) then

               !Name the par_file
               write(par_file, 100) trim(Stn(jj)),'Ind_Y', current_year,&
                                    '_D', current_DOY

               open(unit=parunit, file=par_file, status='replace')
               write(parunit, 299) 'Timestep','Day','Hour','ID','Grid', &
                                   'Z', 'C', 'N', 'Chl', 'NO3', 'PAR',  &
                                   'Temp' 
            else
		         OPEN(parunit, file=par_file, status='old', action='write',&
                                                        position='append')
            endif

100 format(A3,A5,I1,A2,I0)
299 format(12(A8, 1x))

            do i = 1, size(p_PHY)
               write(parunit, 300) it,                                  &
                                 current_day, current_hour,p_PHY(i)%ID, &
                                 p_PHY(i)%iz, p_PHY(i)%rz, p_PHY(i)%C,  &
                                 p_PHY(i)%N,  p_PHY(i)%Chl,p_PHY(i)%NO3,&
                                 p_PHY(i)%PAR,p_PHY(i)%temp
            enddo
            close(parunit)

         endif
    ENDIF  !==> End of saving results
    
  !-------------------------------------------------------------------------
    ! Calculate model outputs (final year) to match with obs. data
   If ((NDays-current_day) .le. y_per_d .and. sec_of_day .eq. 0) Then
 
     ! Loop through the DOY of the observed data
     Do i=1, sum(nrow(:,jj))
        ! Calculate the index for the data
        if (jj .eq. 1) then
           nm = i
        else
           nm = 0
           do j = 1, (jj-1)
              nm = nm + sum(NDPTS(:,j))
           enddo
           nm = nm + i
        endif
        DOY      = INT(min(OBS_DOY(nm), dble(y_per_d)))
        depth(1) = -abs(OBS_Depth(nm))  !Consistent with ROMS convention

        IF (DOY .eq. current_DOY) then

          allocate(a(nlev,1))
          a(:,:) = zero
           
          if (i .le. nrow(1,jj)) then
             if (jj .eq. 1) then
                nm = i
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(1, j)
                enddo  !Obtain the total number of obs. for all previous stations
                nm = nm + i
             endif

           ! Calculate TIN output:
           a(:,1) = Varout(iNO3,:)

           call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1), TINout(nm,1))       
           TINout(nm,1) = max(TINout(nm,1), zero)
          elseif (i .le. (nrow(1,jj)+nrow(2,jj)) ) then
             if (jj .eq. 1) then
                nm = i- nrow(1,jj)
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(2, j)
                enddo
                nm = nm + i - nrow(1,jj)
             endif

          ! Calculate CHL output:
           a(:,1) = Varout(oCHLt,:)
           call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1), CHLout(nm,1))   
           CHLout(nm,1) = max(CHLout(nm,1), zero)
          elseif (i .le. (nrow(1,jj)+nrow(2,jj)+nrow(3,jj)) ) then
             if (jj .eq. 1) then
                nm = i- nrow(1,jj) - nrow(2,jj)
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(3, j)
                enddo
                nm = nm + i - nrow(1,jj) - nrow(2,jj)
             endif

          ! Calculate NPP output:
           a(:,1) = Varout(oNPP,:)   ! Carbon based PP
           call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                NPPout(nm,1))  

           NPPout(nm,1) = max(NPPout(nm,1), zero)

          elseif (i .le. (nrow(1,jj)+nrow(2,jj)+nrow(3,jj)+nrow(4,jj))) then
             if (jj .eq. 1) then
                nm = i- nrow(1,jj) - nrow(2,jj) - nrow(3,jj)
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(4, j)
                enddo
                nm = nm + i - nrow(1,jj) - nrow(2,jj) - nrow(3,jj)
             endif

          ! Calculate PON output:
           a(:,1) = Varout(oPON,:)   ! Total PON
           call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                PONout(nm,1))  

           PONout(nm,1) = max(PONout(nm,1), zero)
          elseif (do_IRON .and. (i .le. (nrow(1,jj)+nrow(2,jj)+nrow(3,jj)+nrow(4,jj)+nrow(5,jj)))) then

             if (jj .eq. 1) then
                nm = i- nrow(1,jj) - nrow(2,jj) - nrow(3,jj) - nrow(4,jj)
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(5, j)
                enddo
                nm = nm + i - nrow(1,jj) - nrow(2,jj) - nrow(3,jj)-nrow(4,jj)
             endif

          ! Calculate DFe output:
           a(:,1) = Varout(oFER,:)   ! Total DFe
           call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                DFeout(nm,1))  

          elseif (.not. N2fix) then  ! Separate NPZDN2 and NPZDCONT models

           goto 1000   ! <== for models without N2 fixation but with size

          else
           If (.not. do_IRON) then
             ! For models with N2fixation but without iron
             if (i.le.(nrow(1,jj)+nrow(2,jj)+nrow(3,jj)+nrow(4,jj)+nrow(5,jj))) then
              if (jj .eq. 1) then
                 nm = i- nrow(1,jj) - nrow(2,jj) - nrow(3,jj)-nrow(4,jj)
              else
                 nm = 0
                 do j = 1, (jj-1)
                    nm = nm + nrow(5, j)
                 enddo
                 nm = nm + i - nrow(1,jj) - nrow(2,jj) - nrow(3,jj) - nrow(4,jj)
              endif

          !   Calculate PO4 output:
              a(:,1) = Varout(oPO4,:)   ! PO4
              call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                   PO4out(nm,1))  

             elseif (&
              i.le.(nrow(1,jj)+nrow(2,jj)+nrow(3,jj)+nrow(4,jj)+nrow(5,jj)+nrow(6,jj))) then
               if (jj .eq. 1) then
                  nm = i- nrow(1,jj) - nrow(2,jj) - nrow(3,jj)-nrow(4,jj)-nrow(5,jj)
               else
                  nm = 0
                  do j = 1, (jj-1)
                     nm = nm + nrow(6, j)
                  enddo
                  nm = nm + i - nrow(1,jj) - nrow(2,jj) - nrow(3,jj) - nrow(4,jj)-nrow(5,jj)
               endif

          !   Calculate POP output:
              a(:,1) = Varout(oPOP,:)   ! POP
              call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                   POPout(nm,1))  

             else
               if (jj .eq. 1) then
                  nm = i- nrow(1,jj)-nrow(2,jj)-nrow(3,jj)-nrow(4,jj)-nrow(5,jj)-nrow(6,jj)
               else
                  nm = 0
                  do j = 1, (jj-1)
                     nm = nm + nrow(7, j)
                  enddo
                  nm=nm+i-nrow(1,jj)-nrow(2,jj)-nrow(3,jj)- nrow(4,jj)-nrow(5,jj)-nrow(6,jj)
               endif

              ! Calculate DIA output:
               a(:,1) = Varout(oDIA,:)   ! POP
               call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                    DIAout(nm,1))  
             endif ! <== End of selection of special data types in NPZDN2 model
           Else
           ! For models with N2fixation and iron
             if (i.le.(nrow(1,jj)+nrow(2,jj)+nrow(3,jj)+nrow(4,jj)+nrow(5,jj)+nrow(6,jj))) then
              if (jj .eq. 1) then
                 nm = i- nrow(1,jj) - nrow(2,jj) - nrow(3,jj)-nrow(4,jj) -nrow(5,jj)
              else
                 nm = 0
                 do j = 1, (jj-1)
                    nm = nm + nrow(6, j)
                 enddo
                 nm = nm + i - nrow(1,jj) - nrow(2,jj) - nrow(3,jj) - nrow(4,jj) - nrow(5,jj)
              endif

          !   Calculate PO4 output:
              a(:,1) = Varout(oPO4,:)   ! PO4
              call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                   PO4out(nm,1))  

             elseif (&
              i.le.(nrow(1,jj)+nrow(2,jj)+nrow(3,jj)+nrow(4,jj)+nrow(5,jj)+nrow(6,jj)+nrow(7,jj))) then
               if (jj .eq. 1) then
                  nm = i- nrow(1,jj) - nrow(2,jj) - nrow(3,jj)-nrow(4,jj)-nrow(5,jj)-nrow(6,jj)
               else
                  nm = 0
                  do j = 1, (jj-1)
                     nm = nm + nrow(7, j)
                  enddo
                  nm = nm + i - nrow(1,jj) - nrow(2,jj) - nrow(3,jj) - nrow(4,jj)-nrow(5,jj)-nrow(6,jj)
               endif

          !   Calculate POP output:
              a(:,1) = Varout(oPOP,:)   ! POP
              call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                   POPout(nm,1))  

             else
               if (jj .eq. 1) then
                  nm = i- nrow(1,jj)-nrow(2,jj)-nrow(3,jj)-nrow(4,jj)-nrow(5,jj)-nrow(6,jj)-nrow(7,jj)
               else
                  nm = 0
                  do j = 1, (jj-1)
                     nm = nm + nrow(8, j)
                  enddo
                  nm=nm+i-nrow(1,jj)-nrow(2,jj)-nrow(3,jj)- nrow(4,jj)-nrow(5,jj)-nrow(6,jj)-nrow(7,jj)
               endif

              ! Calculate DIA output:
               a(:,1) = Varout(oDIA,:)   ! POP
               call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                    DIAout(nm,1))  
             endif ! <== End of selection of special data types in NPZDN2 model

           Endif  ! <== endof both do_IRON and N2 fixation

1000      if(INCLUDESIZE) then

             ncff = NDTYPE-4  !Number of data types excluding size

             !Calculate the total number of observations before size
             j_   = 0
             do j = 1,ncff
                j_= j_+nrow(j,jj)
             enddo

             if (jj .eq. 1) then
                nm = i-j_
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(ncff+1, j)
                enddo
                nm = nm + i - j_
             endif

          !  Calculate Size-fractionated Chl output:
             do j = 1, 4
                a(:,1) = Varout(oCHLs(j),:) 
                call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1), &
                     Sizeout(nm,j))
             enddo
            endif
          endif  ! ==> End of selecting data types
          deallocate(a)
        Endif !==> End of matching DOY
       Enddo  !==> End of loop of selecting obs. data
   Endif !==> End of checking if final year
  !---------------------------------------------------------------------
 ENDIF  !==> END of daily work
4 format(I10,1x,I7, <nlev+1>(2x,1pe12.3))

  ! Pass the new state variables to Vars
  do j = 1,NVAR
     Vars(j,:)= Varout(j,:)
  enddo
  
  ! Diffusion:

  ! For Lagrangian model, run random walk for each particle
  IF (Model_ID .eq. GMK98_Lag) THEN
     
     NP = size(p_PHY)
     Do j = 1, NP
        CALL LAGRANGE(nlev, Z_w, Aks, w, p_PHY(j)%iz, p_PHY(j)%rz)

        !record the indexes of the lost particles or if the particle's C is below the critical value
	!Do not remove this individual now, but record this index
	!Update ideath

        if (p_PHY(j)%iz < 0 .or. p_PHY(j)%C < Cmin) then
	     if (allocated(ideath)) then
	        allocate(add_death(size(ideath) + 1), stat=allocatestatus)
	        IF (allocateStatus /= 0) STOP "*** Problem in allocating add_death***"
	        add_death(1:size(ideath))   = ideath
	        add_death(size(ideath) + 1) = j
	        call move_alloc(add_death, ideath)   !Update ideath and add_death has been deallocated
	    else
	        allocate(ideath(1))
	        ideath(1)  = j 
	    endif
        endif
     Enddo  !End of lagrangian motions

     !Remove the lost and dead particles from the community
	  IF (allocated(ideath) .and. size(ideath) > 0) THEN

	    !Construct a new Zd with length of size(Z)-length(ideath) 
			allocate(Zd(NP - size(ideath)), stat=allocateStatus)
			IF (allocateStatus /= 0) STOP "*** Problem in allocating Zd ***"

			j_ = 0   !Counting the index in Zd
			do j = 1, NP
				match = .false.
				do m = 1, size(ideath)
					if (j .eq. ideath(m)) then
						match = .TRUE.
                  exit
               endif
               if (j_ .eq. (NP-size(ideath))) exit
            enddo  !<-- End of finding whether this particle had been lost or not

				if (.not. match) then !Include this individual in the new community Zd
					j_     = j_ + 1
					Zd(j_) = p_PHY(j)
				endif
			enddo	!End of constructing Zd

			call move_alloc(Zd, p_PHY)  !Update p_PHY and Zd is deallocated
			deallocate(ideath)
      ENDIF  !End of removing lost and dead particles

  ENDIF  !<-- End of Lagragian model

  do j = 1,NVAR
     Vars1(:) = Vars(j,:)

     ! At surface, assume zero flux  (Neumann boundary condition)
     ! For Nutrient, use Dirichlet boundary condition at the bottom
     IF (j .eq. iNO3 .or. j .eq. ifer .or. j .eq. iPO4) THEN
        SELECTCASE (bot_bound) ! Select the type of boundary condition at bottom
        CASE (Dirichlet)

          !! At bottom, assume constant values obtained from observation (Dirichlet boundary condition)
          call diff_center(nlev,dtsec,cnpar,1,Hz, Neumann, Dirichlet, &
                            zero, VarsBom(1,j),Aks,Vec0,Vec0,Taur,Vars1,Vars1,Vars2)

        CASE (Neumann)
        ! Zero flux at bottom
          call diff_center(nlev,dtsec,cnpar,1,Hz, Neumann, Neumann, &
                        zero, zero, Aks,Vec0,Vec0,Taur,Vars1,Vars1,Vars2)
        CASE DEFAULT
          stop 'The type of bottom boundary condition incorrect!'
        END SELECT

     ELSE
       ! For other tracers, use closed boundary condition
       ! Zero flux at bottom
        call diff_center(nlev,dtsec,cnpar,1,Hz, Neumann, Neumann, &
                      zero, zero, Aks,Vec0,Vec0,Taur,Vars1,Vars1,Vars2)

     ENDIF

     ! Save diffusion fluxes (normalized to per day)
     Varout(oD_VARS(j),:) = (Vars2(:) - Vars1(:))/dtdays
  
     ! Update the state variables:
     Vars(j,:) = Vars2(:)
  enddo
  
  ! Sinking:
  do j = 1,NVsinkterms
     ww_(:)   = ww(:,j)
     Vars2(:) = Vars(Windex(j),:)

     SELECT CASE (bot_bound)
     case(Neumann)  ! closed at bottom (Conserve total N mass)
        call adv_center(nlev,dtsec,Hz,Hz,ww_(:),1,1,zero,zero,    6,mode1,Vars2(:))
     case(Dirichlet)
        ! Open bottom boundary
        call adv_center(nlev,dtsec,Hz,Hz,ww_(:),1,2,zero,Vars2(1),6,mode1,Vars2(:))
     case default
        stop "The boundary conditions incorrect! STOP!"
     ENDSELECT
     Vars(Windex(j),:) = Vars2(:)
  enddo
  
  ! Vertical advection (due to w, unit: m/s):
  ! Now do not consider vertical w
  !do j = 1,NVAR
  !   call adv_center(nlev,dtsec,Hz,Hz,w(:),  &
  !                   1,1,Yup,Ydw,6,mode0,Vars(j,:))
  !enddo
  
  ENDDO    ! ==> End of time stepping
  
  ! Close files:
  if (savefile .and. taskid == 0) then
     do i = 1, (Nout+ow)
        close (unit=9+i)
     enddo
  endif
  if(allocated(ww))      deallocate(ww)
  if(allocated(VarsBom)) deallocate(VarsBom)
  if(allocated(p_PHY))   deallocate(p_PHY)
ENDDO  ! ==> End of Stn

300   format(I10,1x, 2(I7, 1x), I0, 1x, I2, 1x, 7(1pe12.4, 1x))
END subroutine Timestep

!========================================================
subroutine setup_grid
implicit none

! Local scratch variable:
integer                :: i
real                   :: sc_r, C_sig

Z_w(0) = -hmax

!Following Song and Haidvogel (1994). sinh is the hyperbolic sin function
do i   = 1,nlev
    sc_r  = (float(i-nlev) - 0.5)/float(nlev)
   C_sig  = sinh(thetaS*sc_r)/sinh(thetaS)      ! -1 < C_sig < 0
   Z_r(i) = C_sig*hmax
    sc_r  = (float(i-nlev))/float(nlev)
   C_sig  = sinh(thetaS*sc_r)/sinh(thetaS)      ! -1 < C_sig < 0
   Z_w(i) = C_sig*hmax
    Hz(i) = Z_w(i) - Z_w(i-1)
enddo
end subroutine setup_grid
!===========
subroutine Calculate_PAR(I_0, nlev_, Hz, Chl, PAR, PARw)
  !top level is nlev_, bottom layer is 1, following ROMS convention
  implicit none
  real   , intent(in) :: I_0
  integer, intent(in) :: nlev_    ! Total number of vertical layers
  real   , intent(in) :: Hz(nlev_), Chl(nlev_)
  real   , intent(out):: PAR(nlev_),PARw(0:nlev_)
  integer             :: i
  real                :: par0, attn    ! Scratch variable
  real   , parameter  :: kw = 0.04, kc = 0.025

  par0         = I_0   !Light at the grid surface
  parw(nlev_)  = I_0
  do i = nlev_,1,-1

     attn   = exp(-0.5*Hz(i)*(kw + kc*Chl(i)))
     PAR(i) = par0*attn
     par0   = PAR(i)*attn
     PARw(i-1) = par0
  enddo
end subroutine Calculate_PAR
End MODULE MOD_1D
