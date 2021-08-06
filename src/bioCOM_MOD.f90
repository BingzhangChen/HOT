MODULE BIOCOM_MOD
! This module contains some common functions and variables for all biological models
IMPLICIT NONE

! MPI variables:
integer            :: numtasks, ierr, MPIRUN, taskid

! Options of biological models
integer, parameter :: EFTdisc     = 1
integer, parameter :: EFTcont     = 2
integer, parameter :: Geiderdisc  = 3
integer, parameter :: NPZDcont    = 4
integer, parameter :: EFTsimple   = 5
integer, parameter :: Geidersimple= 6
integer, parameter :: NPZDFix     = 7
integer, parameter :: NPZDdisc    = 8
integer, parameter :: EFTsimIRON  = 9
integer, parameter :: NPZDFixIRON = 10
integer, parameter :: GeidsimIRON = 11
integer, parameter :: NPZDdiscFe  = 12
integer, parameter :: EFTdiscFe   = 13
integer, parameter :: EFTDarwin   = 14
integer, parameter :: EFT2sp      = 15
integer, parameter :: NPZD2sp     = 16
integer, parameter :: NPPZDD      = 17
integer, parameter :: EFTPPDD     = 18
integer, parameter :: NPZDN2      = 19
integer, parameter :: CITRATE3    = 20
integer, parameter :: GeiderDroop = 21
integer, parameter :: NPclosure   = 22
integer, parameter :: NPZclosure  = 23
integer, parameter :: TOPT_MOD    = 24
integer, parameter :: GMK98_Lag   = 25  !Lagrangian model using GMK98 model

! Bottom boundary condition:
integer, parameter :: Dirichlet   = 0
integer, parameter :: Neumann     = 1

integer            :: NDays = 1080  !Total number of days for simulation
! Indices for output variables
integer, parameter :: oTemp  = 1,oPAR=2,oAks=3,oDust=4,ow=5

! Parameters for phytoplankton size fractional Chl
real, parameter :: pi      = 3.1415926535897932384633D0
real, parameter :: eps     = 1d-20
real, parameter :: PMU_min = log(1d1*pi/6d0*0.6**3)
real, parameter :: PMU_max = log(1d1*pi/6d0*4d1**3)
real, parameter :: PMU_1   = log(1d1*pi/6d0)
real, parameter :: PMU_3   = log(1d1*pi/6d0*3d0**3)
real, parameter :: PMU_10  = log(1d1*pi/6d0*1d1**3)
integer  :: AllocateStatus

! Number of vertical layers
integer, parameter :: nlev = 40  

real,    parameter ::d_per_s = 864d2 ! how many seconds in one day
real,    parameter ::s_per_h = 3600d0 ! how many seconds in one hour
integer, parameter ::h_per_d = 24      ! how many hours in one day

!how many days of one year
integer, parameter ::y_per_d = 365

!Time of the day
integer  :: sec_of_day   = 0
integer  :: current_day  = 0
integer  :: N_MLD     ! vertical grid index at the bottom of MLD
real     :: Temp(nlev), PAR(nlev), dtdays, Ntot, PARavg, wstr0(1) 

!Tracers for NPP incubation
real     :: NO3e(nlev) = 0.
real     :: PHYe(nlev) = 0.
real     :: PHYCe(nlev)= 0.
real     :: ZOOe(nlev) = 0.
real     :: DETe(nlev) = 0.
real     :: CHLe(nlev) = 0.

! PAR at w points
real     :: PAR_w(0:nlev) = 0.

! Daily integrated PAR
real     :: PAR_int(nlev) = 0.

! Count numbers for averaging
integer  :: count(nlev) = 0

integer  :: it  !The current step

real     :: DFe(nlev)                           ! Dissolved iron concentration
real     :: Z_r(1:nlev), Z_w(0:nlev), Hz(nlev)  ! Grid variables
real     :: I_zero
real     :: rhoChl_L = 0.5d0  !rhoChl calculated for the end of the preceding light period  

!  Indices for external forcing variables
integer, parameter :: etemp      = 1
integer, parameter :: eNO3       = 2
integer, parameter :: eAks       = 3
integer, parameter :: ew         = 4
integer, parameter :: ePAR       = 5
integer, parameter :: eDust      = 6
integer, parameter :: eFer       = 7
integer, parameter :: ePO4       = 8
integer, parameter :: ewstr      = 9
integer, parameter :: TNFo       = ewstr ! TOtal number of forcings

! Total of observation times in forcing data
character(LEN=5), parameter :: LabelForc(TNFo) &
     = (/'temp ','NO3  ','Aks  ','wROMS','par  ', 'solfe','fer  ', 'PO4  ', 'wstr '/)

integer,          parameter :: NFobs(TNFo)     & 
     = (/    12,   12, 365,     12,  12,      12,   12,    12, 12 /)

integer  :: iZOO2,iDET2,iDETp,iDETFe, iPMU, iVAR, iPO4,iDIA
integer  :: oZOO, oZOO2,oDET, oFER, oZ2N, oD2N, oPHYt
integer  :: oPO4, oPOP, oDIA, oDIAu,oDETp, oDET2, oDETFe, iwDET2
integer  :: itau,igb,oFescav,odstdep,ifer
!integer  :: oPMU, oVAR
!integer  :: odmudl,odgdl,od2mu,od2gdl,odmudT,od2mudT2,odmudI,od2mudI2  
!integer  :: oMESg,oMESgMIC,odgdl1,odgdl2,od2gdl1,od2gdl2,odVAR

! Indices for parameters used in DRAM
!integer  :: imu0,iaI0,igmax,iKN,iKP,iKPHY,iKPnif,iLnifp,iKFe,iRDN_N,iRDN_P
!integer  :: ialphamu,ibetamu,ialphaKN,irhom,iDp, iIopt, ibeta
!integer  :: imu0B, iaI0B, iA0N2, iRL2,iKN2,ibI0B
!integer  :: iVTR,iVTRL,iVTRT,iVTRI,,od3mu,od4mu
!integer  :: izetaN,izetaChl, iaI0_C 
!integer  :: ialphaI,iA0N,ialphaA,ialphaG,ialphaK, ialphaFe
!integer  :: iQ0N,ialphaQ,iPenfac,iLref,iwDET,irdN,imz

integer   :: oCHLs(4) = 0   ! Four size fractions of CHL
! Some common Model parameters:
real, parameter :: PMUmax =1.5D1, VARmax=50D0
real            :: Femin  =0.02,K0Fe  =0.8, alphaFe=0.14
real, parameter :: GGE    =0.3, unass =0.24
real, parameter :: Fe_N   =0.0265 ! Fe:Nitrogen molar ratio (nmol : umol)
real, parameter :: RMchl0 =0.1
real, parameter :: Ep     =0.4,  Ez    =0.65 
real            :: alphamu=0.2,  betamu=-0.01
real, parameter :: zetaChl=0.6,  zetaN =0.8
real            :: thetamax = 0.63, thetamin=0.02  !Unit: gChl/molC
integer, parameter :: nutrient_uptake=1
real :: KFe    =0.08     !Unit: nM. Gregg (2003) used ~0.1 nM for phyto

! Size and weights for discrete models
real, allocatable          :: PMU_(:)
real, allocatable          :: wtCHL(:,:)  ! weight for each size class
logical, parameter :: kill_the_winner=.TRUE.
logical, public    :: N2fix          =.FALSE.
logical            :: singlerun      =.FALSE.
integer, parameter :: namlst=8

! Declaration for phyto particles
TYPE Particle

    ! Particle ID
    integer :: ID = 1
    ! Grid indices for particles (range from nlev to 1)
    integer :: iz = 1
    
    ! Z coordinates for particles
    real    :: rz = 0.

       ! associated PAR 
    real    :: PAR = 0.01

    ! associated Temperature
    real    :: temp= 20.

    ! associated nutrient
    real    :: NO3 = 0.1

    ! cellular carbon content (pmol; assuming a 1 micron cell)
    real    :: C   = 0.02

    ! cellular nitrogen content (pmol; assuming a 1 micron cell)
    real    :: N   = 0.02/106.*16.

    ! Cellular Chl content (pg Chl)
    real    :: Chl = 0.02 * 12/50

END TYPE Particle

Type (Particle), ALLOCATABLE :: p_PHY(:)

INTEGER, ALLOCATABLE         :: AllIDs(:)  !Recording all ID names

!Initial number of individuals
INTEGER(kind=8), PARAMETER   :: N_ini = 1000

! cellular carbon content threshold for division (pmol)
real,    PARAMETER :: Cdiv= 0.04

! Subsistence cellular carbon content (pmol), below which the cell will die
real,    PARAMETER :: Cmin = 0.01

! Some common local variables:
real :: tf_p, tf_z

! Day length
real :: Day_len = 0d0

! PAR at noon
real :: I_noon_ = 0d0


! Common Indices (may not be used, but needed when compiling)
integer            :: iMTo, iVTo, iMIo
integer            :: iVIo

integer, parameter :: Yes = 1, No = 0

CONTAINS

!========================================================
pure real function TEMPBOL(Ea,tC)
implicit none
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
real, intent (in) :: Ea, tC
! boltzman constant constant [ eV /K ]
real, parameter   :: kb = 8.62d-5, Tr = 15D0

TEMPBOL = exp(-(Ea/kb)*(1D0/(273.15 + tC)-1D0/(273.15 + Tr)))
return 
end function TEMPBOL

!====================================================
real function ScaleTrait( logsize, star, alpha ) 
implicit none
real, intent(IN) :: logsize, star, alpha

! Calculate the size-scaled value of a trait
! for the given log (natural, base e) of cell volume as pi/6*ESD**3 (micrometers). 

ScaleTrait = star * exp( alpha * logsize )

return
end function ScaleTrait
!====================================================
real function PenAff( logsize, alpha, Pfac, lmin ) 
implicit none
real, intent(IN) :: logsize, alpha, Pfac, lmin 

!A 'penalty' function to reduce the value of affinity for nutrient at very small cell sizes
!in order to avoid modeling unrealistically small cell sizes.  This is needed because affnity
!increases with decreasing cell size, which means that under low-nutrient conditions, without
!such a penalty, unrealistically small cell sizes could be predicted.
!This penalty function becomes zero at logsize = lmin.   
   
  PenAff = 1.0 - exp(Pfac*alpha*(logsize - lmin))
end function PenAff
!====================================================
pure real function grazing(Hollingtype, Ksat, Prey)
implicit none
real,    intent(in) :: Ksat, Prey
integer, intent(in) :: Hollingtype

grazing = 0.
! kp relates to the capture coefficient
SELECT CASE(Hollingtype)
  ! Holling Type I
  case (1)
    grazing = min(Prey/2.0/Ksat,1.0)
  ! Holling Type II
  case (2)
    grazing = Prey/(Ksat + Prey)  
  ! Holling Type III
  case (3) 
    grazing = min(Prey*Prey/(Ksat*Ksat + Prey*Prey), 1D0)
 ! Ivlev
  case (4)
 !To be consistent with other functions  
    grazing = 1d0-exp(-log(2d0)*Prey/Ksat) 

END SELECT
return
end function grazing
!===================================================
!Functions deriving the first and second derivatives of X/Y ~ L 
pure real function dY_Xdl(Y, X, dYdl, dXdl) 
implicit none
real, intent(in)    :: Y, X, dYdl, dXdl
   dY_Xdl = dYdl/X - Y/X**2*dXdl
end function

pure real function d2Y_Xdl2(Y, X, dYdl, dXdl, d2Ydl2, d2Xdl2) 
implicit none
real, intent(in)    :: Y, X, dYdl, dXdl, d2Ydl2, d2Xdl2
 d2Y_Xdl2 = d2Ydl2/X - 2.*dYdl*dXdl/(X**2) - Y/(X**2)*d2Xdl2 +  &
                                          2.*Y/(X**3)*(dXdl**2)
end function

! Light ~ growth function from Edwards et al. L&O 2015
pure real function mu_Edwards2015(I_, lnIopt, mumax, alpha) 
implicit none
! lnIopt is log(Iopt)
real, intent(in)    :: I_, lnIopt, mumax, alpha
mu_Edwards2015 =  1./( I_/(alpha*exp(2.*lnIopt)) + 1./mumax -2./(alpha*exp(lnIopt))  &
            + 1./(I_*alpha)  ) 
end function mu_Edwards2015

!Function following Anderson et al. GMD 2015 to calculate PAR below the ocean surface
! Noon PAR, W m-2

pure real function noonparnow(jday, latradians)
implicit none
integer, intent(in) :: jday      !Day of year
real,    intent(in) :: latradians!Latitude in radianns

real            :: Rv        !Earth's radius vector (accounts for eccentricity of the earth's orbit)
real, parameter :: solarconst = 1368.0!Solar constant (W m-2)
real, parameter :: e0    = 12.0  !partial pressure of water vapour in the atmosphere (mb)
real, parameter :: parrac= 0.43  !Fraction of solar radiation
real, parameter :: albedo= 0.04  !Ocean albedo
real, parameter :: clouds= 6d0   !cloud fraction in oktas
real            :: cfac     !Effect of clouds on atmospheric transmission

real            :: Iclear   !Short-wave irradiance at the ocean surface on a clear day

real            :: coszen    !Cosine of solar zenith angle

real            :: declin    !Solar declination angle (in radians)

real            :: zen       !solar zenith angle

real            :: Rvector, Inoon

declin = 23.45*sin(2.*pi*(284. + dble(jday))*0.00274)*pi/180.     ! solar declination angle
coszen = sin(latradians)*sin(declin)+cos(latradians)*cos(declin)  ! cosine of zenith angle
zen    = acos(coszen)*180./pi                                     !zenith angle, degrees
Rvector= 1d0/sqrt(1. + 0.033*cos(2.*pi*dble(jday)*0.00274))       ! Earth's radius vector
Iclear = solarconst*coszen**2/(Rvector**2)/(1.2*coszen+e0*(1.0+coszen)*0.001+0.0455) ! irradiance at ocean surface, clear sky
cfac   = (1.-0.62*clouds*0.125 + 0.0019*(90.-zen))                     ! cloud factor (atmospheric transmission)
Inoon  = Iclear*cfac*(1.-albedo)                                    ! noon irradiance: total solar
noonparnow = parrac*Inoon   
return

end function noonparnow

! ----------------------------------------------------------------------- #
! Calculation of day length as function of day of year and latitude       #
! ----------------------------------------------------------------------- #
pure real function FNdaylcalc(jday,latradians) result(daylnow)
implicit none
integer, intent(in) :: jday      !Day of year
real,    intent(in) :: latradians!Latitude in radianns

real                :: declin    !Solar declination angle (in radians)

declin = 23.45*sin(2.*pi*(284. + jday)*0.00274)*pi/180.      ! solar declination angle
daylnow= 2.*acos(-1. * tan(latradians)*tan(declin))*12./pi   ! hours
end function FNdaylcalc

REAL FUNCTION surface_par(Inoon, timeofday, daylength)
IMPLICIT NONE
REAL,    INTENT(IN) :: Inoon         ! Surface par at noon Unit: W m-2
INTEGER, INTENT(IN) :: timeofday     ! Unit: second
REAL,    INTENT(IN) :: daylength     ! Unit: hours

real,    parameter  :: Tnoon = 43200d0 ! seconds at noon
real    :: Trise   !Unit: seconds
real    :: Tset    !Unit: seconds
real    :: Day_len_sec  ! Convert day length into seconds

Day_len_sec = daylength*s_per_h

Trise = Tnoon - Day_len_sec/2d0 
Tset  = Tnoon + Day_len_sec/2d0 

! Initialize surface PAR based on diel cycle
if (dble(timeofday) .le. Trise .or. dble(timeofday) .ge. Tset) then
    surface_par = 0d0
else
    surface_par = sin(pi*(dble(timeofday) - Trise)/Day_len_sec)*Inoon
endif

return
END FUNCTION surface_par

real function diurnal_light(jday,latradians, time_of_day) result(y)
IMPLICIT NONE
integer, intent(in) :: jday
real,    INTENT(in) :: latradians
integer, INTENT(in) :: time_of_day

if (time_of_day .eq. 0) then
   !Daylength
   Day_len = FNdaylcalc(jday, latradians)
   
   !I_noon
   I_noon_ = noonparnow(jday, latradians)
endif

!Diel light
y = surface_par(I_noon_, time_of_day, Day_len)

return
end function diurnal_light

real function rexp(u) result(y)
implicit none

real, intent(in) :: u   !Rate parameter
real             :: rnd !a random number between 0 and 1

y = 0.d0

if (u < 0.d0) then
        write(6, *) 'The rate parameter has to be positive!'
        return
endif

call random_number(rnd)

y = -log(rnd)/u
return
END function rexp

Integer function sample(nn, x) result(y)
!A function to randomly sample 1 elements from a vector of nn integers
implicit none
integer, intent(in) :: nn
integer, intent(in) :: x(nn)

real    :: r = 0.d0
integer :: i

call random_number(r)

if (r .le. 1./dble(nn)) then
	y = x(1)
else if (r .gt. (1. - 1./dble(nn))) then
	y = x(nn)
else
	do i = 2, (nn-1)
		if (r .gt. dble(i-1)/dble(nn) .and. r .le. dble(i)/dble(nn)) then
			y = x(i)
			return
		endif
	enddo
endif

return
END function sample

END Module BIOCOM_MOD
