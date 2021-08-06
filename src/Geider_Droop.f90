SUBROUTINE BIOLOGY
!Add phyto. Carbon into existing Geider model
USE PARAM_MOD
implicit none
integer :: k 
real :: NPPc_(nlev)  = 0.  !C-based phytoplankton production (mg C m-3 d-1)
real :: NPPn_(nlev)  = 0.  !N-based phytoplankton production (mmol N m-3 d-1)
real :: NPPchl_(nlev)= 0.  !Chl-based phytoplankton production (mgChl m-3 d-1)
real :: Graz_(nlev)  = 0.  !phytoplankton loss due to grazing (mmol N m-3 d-1)
real :: NPPe_(nlev)  = 0.  !Incubation-based NPP (in bottles)
real :: IntPAR_(nlev)= 0.  !Integrated PAR over 1 day
real :: PHYCavg(nlev)= 0.  !Averaged PHY carbon over 1 day
real :: PHYNavg(nlev)= 0.  !Averaged PHY N over 1 day
real ::  CHLavg(nlev)= 0.  !Averaged Chl over 1 day
real :: Lno3avg(nlev)= 0.  !Averaged Lno3  over 1 day
real :: SI_avg(nlev) = 0.  !Averaged SI  over 1 day
real :: NO3, PHY, PHYC, ZOO, DET, CHL, mu, QN, theta
real :: dCdt      = 0.
real :: dChldt    = 0.
real :: dCdte     = 0.
real :: Lno3(nlev)= 0.
real :: g_        = 0.  !Scrath variable for grazing rate
real :: SI_(nlev) = 0.
!End of declaration

DO k = 1, nlev
   NO3 = Vars(iNO3,    k)
   PHYC= Vars(iPHYC(1),k)
   PHY = Vars(iPHY(1), k)  !Phyto. N
   ZOO = Vars(iZOO,    k)
   CHL = Vars(iCHL(1), k)
   DET = Vars(iDET,    k)

   ! For calculating NPP, we need to integrate over 24 hours
   If ((NDays-current_day) .le. y_per_d) Then

      !A new time interval starts for output
      IF (mod(it, nsave) .EQ. 0) THEN
 
         PAR_int(k)         = IntPar_(k)
         Varout(oNPPr,k)    = NPPc_(k)  !NPP of the past day; this is real NPP (mg C d-1 L-1)
         Varout(omuC(1),k)  = Varout(oNPPr,k)/PHYCavg(k)/12.  !Averaged Carbon-based Specific growth rate
         Varout(oGraz(1),k) = Graz_(k)/PHYNavg(k)
         Varout(omuN(1),k)  = NPPn_(k)/PHYNavg(k)     !Averaged N-based Specific growth rate
         Varout(omuChl(1),k)= NPPchl_(k)/Chlavg(k)   !Averaged Chl-based Specific growth rate
         Varout(oQN(1), k)  = PHYNavg(k)/PHYCavg(k)
         Varout(othe(1),k)  =  Chlavg(k)/PHYCavg(k)
         Varout(oSI(1),k)   = SI_avg(k)
         Varout(oLno3(1),k) = Lno3avg(k)

         PHYNavg(k) = 0d0      !Reset PHYNavg
         Chlavg(k)  = 0d0      !Reset Chlavg
         PHYCavg(k) = 0d0      !Reset PHYCavg
         NPPchl_(k) = 0d0      !Reset NPPchl_
         NPPn_(k)   = 0d0      !Reset NPPn_
         NPPc_(k)   = 0d0      !Reset NPPc_
         Graz_(k)   = 0d0      !Reset average grazing rate
         SI_avg(k)  = 0d0      !Reset SI_avg
         Lno3avg(k) = 0d0      !Reset Lno3avg
         IntPar_(k) = 0d0      !Reset IntPar_

         count(k)   = 0        !Reset count

         !During simulation of the final year, Take water samples into a bottle
         !No detritus sinking; no mixing
         Varout(oNPP,k) = NPPe_(k) !NPP in the bottle of the past day
         NPPe_(k)       = 0d0      !Reset NPPe_

         !Reset bottles
         NO3e(k)        = NO3
         CHLe(k)        = CHL
         PHYe(k)        = PHY
         ZOOe(k)        = ZOO
         DETe(k)        = DET
         PHYCe(k)       = PHYC
      endif

      call rhs_NPZDGeiderDroop(Temp(k), PAR(k), NO3e(k), PHYe(k),       &
                               PHYCe(k), CHLe(k), ZOOe(k),DETe(k),      &
                               Lno3(k), SI_(k),                         &
                               dCdte, mu, dChldt, g_)

      NPPe_(k) = NPPe_(k) + PHYCe(k)*dCdte*12.*dtdays   !Unit: mgC m-3 d-1

   Endif

!RK4
!k1   <- k2 <- k3 <- k4 <- Y
!k1   <- ydot(Y,pars)*dtdays
!k2   <- ydot(Y+k1,pars)*dtdays
!k3   <- ydot(Y+k2,pars)*dtdays
!k4   <- ydot(Y+k3,pars)*dtdays
!Y    <- Y + (k1+2*(k2+k3)+k4)/6

   ! Main subroutine for biological reactions (rhs)
   call rhs_NPZDGeiderDroop(Temp(k), PAR(k),NO3, PHY, PHYC, CHL,ZOO,DET,&
                            Lno3(k), SI_(k), dCdt, mu, dChldt, g_)

   !Need the output rates to calculate the following quantities
   If ((NDays-current_day) .le. y_per_d)  then                              
      NPPc_(k)   = NPPc_(k)   + PHYC*dCdt*12.*dtdays     !Real NPP, Unit: mgC m-3 d-1
      NPPn_(k)   = NPPn_(k)   + PHY*mu       *dtdays     !N-based production, Unit: mmolN m-3 d-1
      NPPchl_(k) = NPPchl_(k) + CHL*dChldt   *dtdays     !Chl-based production,Unit: mgChl m-3 d-1
      Graz_(k)   = Graz_(k)   + PHY*g_       *dtdays     !Zooplankton grazing, Unit: mmolN m-3 d-1
      IntPar_(k) = IntPar_(k) + PAR(k)*dtdays

      !Update Lno3avg
      Lno3avg(k) = (Lno3avg(k)*dble(count(k)) + Lno3(k))/dble(count(k)+1)

      !Update SI_avg
      SI_avg(k)  = (SI_avg(k)*dble(count(k)) + SI_(k))/dble(count(k)+1)

      !Update PHYCavg
      PHYCavg(k) = (PHYCavg(k)*dble(count(k)) + PHYC)/dble(count(k)+1)

      !Update PHYNavg
      PHYNavg(k) = (PHYNavg(k)*dble(count(k)) + PHY)/dble(count(k)+1)

      !Update Chlavg
      CHLavg(k)  = ( CHLavg(k)*dble(count(k)) + CHL)/dble(count(k)+1)

      count(k)   = count(k) + 1
   Endif

   !Update 6 tracers
   Varout(iCHL(1),k) = CHL
   Varout(oCHLt,  k) = CHL
   Varout(iDET,   k) = DET
   Varout(iPHY(1),k) = PHY
   Varout(iPHYC(1),k)= PHYC
   Varout(iZOO,   k) = ZOO
   Varout(iNO3,   k) = NO3
   Varout(oPON,k)    = PHY+ZOO+DET
ENDDO
return
END SUBROUTINE BIOLOGY

SUBROUTINE rhs_NPZDGeiderDroop(Temp_, PAR_, NO3, PHY, PHYC, CHL,ZOO,DET,&
                               Lno3_, SI_, dCdt_, dNdt_, dChldt_,graz_)
USE PARAM_MOD, only: params, TEMPBOL, Ez, igmax, iKp, tf_z, iRdn, imz 
USE PARAM_MOD, only: grazing, GGE, unass, dtdays
implicit none

real, intent(in)     :: Temp_, PAR_
real, intent(inout)  :: NO3, PHY, PHYC, CHL, ZOO, DET
real, intent(out)    :: Lno3_, SI_, graz_
real, intent(out)    :: dNdt_      !Specific Nitrogen uptake rate (1/d) 
real, intent(out)    :: dCdt_      !Specific Carbon   uptake rate (1/d) 
real, intent(out)    :: dChldt_    !Specific Chl synthesis rate (1/d) 

! Regeneration constant from detritus to DIN (d-1)
real :: RDN   = 0.1  
real :: gmax  = 0d0 !Maximal zooplankton grazing rate (d-1)

real :: INGES,gbar,EGES,Zmort,RES, QN_, theta_ 
real :: pp_ZP, pp_NZ, pp_ND, pp_DZ 
!End of declaration

gmax = exp(params(igmax))
RDN  = exp(params(iRDN))

 QN_ = PHY/PHYC !Cellular N:C ratio (mol N: mol C)

! Define phytoplankton Chl-to-carbon ratio (gChl/molC)
theta_= CHL/PHYC

! Main phyto. subroutine
call PHY_GeiderDroop(Temp_, PAR_, NO3, QN_, theta_, &
                     Lno3_, SI_, dCdt_, dNdt_, dChldt_)

! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
tf_z   = TEMPBOL(Ez,Temp_)
gbar   = grazing(3, exp(params(iKp)), PHY)
INGES  = ZOO*gmax*tf_z*gbar
Zmort  = ZOO*ZOO*exp(params(imz))*tf_z  !Mortality term for ZOO

!Zooplankton excretion rate (-> DOM)
RES = INGES*(1d0-GGE-unass)

!ZOOPLANKTON EGESTION (-> POM)
EGES = INGES*unass

! For production/destruction matrix:
pp_ND = RDN*DET*tf_z   
pp_NZ = ZOO*RES        
pp_DZ = ZOO*EGES+Zmort 
pp_ZP = ZOO*INGES      

DET   = DET + dtdays*(pp_DZ - pp_ND)
NO3   = NO3 + dtdays*(pp_ND + pp_NZ - PHY * dNdt_)

graz_ = pp_ZP/PHY

! Equation for PHY. nitrogen (decouple nutrient uptake and photosynthesis)
PHY  = PHY  + dtdays*(PHY*dNdt_ - pp_ZP)

ZOO  = ZOO  + dtdays*(pp_ZP     - pp_DZ - pp_NZ)

! Equation for PHY. carbon
PHYC = PHYC + dtdays*(PHYC*dCdt_- pp_ZP/QN_)

CHL  = CHL  + dtdays*(CHL*dChldt_ - pp_ZP/QN_*theta_)
return
end subroutine rhs_NPZDGeiderDroop

SUBROUTINE PHY_GeiderDroop(Temp, PAR, NO3, QN, theta, Lno3, SI, dCdt, dNdt, dChldt)
!PARAM_MOD: a module containing the input parameters
!TEMPBOL:   Arrhenius temperature function following MTE
!Ep: activation energy of photosynthesis
!iaI0: index for the parameter of slope of P-I curve
!ithemax: index for the parameter of maximal N:Chl ratio
!imu00: index for the parameter of PCref
!iKN: index for the parameter of half-saturation constant of N uptake
!tf_p: a scratch variable for storing the temperature effect
!iRchl: index for the parameter of maintenance cost of Chl synthesis
!rchoChl_L: a scratch variable for storing the previous rhoChl value 
!QN: N:C ratio
!theta: Chl:C ratio
!Lno3: output -- nutrient limitation index
!SI: output -- light limitation index
!dCdt: output -- carbon specific growth/photosynthesis rate
!dNdt: output -- N specific uptake rate
!dChldt: output -- Chl specific synthesis rate

USE PARAM_MOD, only: params, TEMPBOL, Ep, iaI0, ithemax, imu00, iKN, tf_p, iRchl,rhoChl_L

implicit none

real, intent(in)  :: Temp, PAR, NO3, QN, theta

! Minimal and maximal N:C ratio
real, parameter   :: QNmin = 0.04, QNmax = 0.18, dQN = 0.14  

! Specific Net Nitrogen uptake rate
real, intent(out) :: dNdt 

! Specific Carbon based photosynthetical rate (d-1)
real, intent(out) :: dCdt

! Specific Chl based synthesis rate (d-1)
real, intent(out) :: dChldt

! Indices for light and nutrient limitation
real, intent(out) :: Lno3, SI

real :: PCmax     = 0d0
real :: thetaNmax = 0d0     !Maximal Chl:N ratio (gChl/molN)

real, PARAMETER   :: zeta = 0d0!3.0    !cost of biosynthesis (molC molN-1)
real, PARAMETER   :: Rc   = 0d0!0.025  !Basic respiration rate (d-1)
real, PARAMETER   :: RN   = 0d0!0.025  !Basic respiration rate (d-1)
real              :: RChl = 0.01       !Basic respiration rate (d-1)

! Gross DIN uptake rate by phytoplankton (molN/molC/d)
real :: VCN = 0.

real :: Pref= 0.  !Reference maximal photosynthetical rate  (per d)

real :: alphaChl = 0d0
real :: Vcref    = 0d0
real :: PC       = 0d0
real :: rhochl   = 0d0
!End of declaration

!Temperature coefficient
tf_p = TEMPBOL(Ep, Temp)  

alphaChl = exp(params(iaI0))
RChl = exp(params(iRchl))
Pref = exp(params(imu00))

thetaNmax = exp(params(ithemax))

! DIN uptake rate by phytoplankton (Ward 2017)
Vcref= Pref * QNmax
VCN  = Vcref* NO3/(NO3 + exp(params(iKN)))*(QNmax-QN)/dQN*tf_p

! Specific N uptake rate
dNdt = VCN/QN - RN*tf_p

! N limitation
Lno3 = (QN - QNmin)/dQN

!Maximal photosynthesis rate (regulated by QN)
PCmax = Pref *tf_P*Lno3

!The light limitation index (fpar)
SI = 1d0 - exp(-exp(params(iaI0)) * PAR * theta / PCmax)

PC = PCmax * SI

! Photosynthesis rate (d-1)
dCdt = PC - zeta*VCN - Rc*tf_p

!If dark, assume that rhochl equaled the value calculated for the end of the preceding light period 
if (PAR .le. 0d0) then
   rhochl   = rhoChl_L
else
   rhochl   = thetaNmax*PC/alphaChl/theta/PAR
   rhoChl_L = rhochl
endif

dChldt = rhochl*VCN/theta - RChl*tf_p

return
end subroutine PHY_GeiderDroop