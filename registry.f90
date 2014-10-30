module registry
  use wrf2gem_subs, only: read_values, read_values_const, mem_error,          &
                          freq_divisible, open_wrf_file, close_wrf_file
  use gempak,       only: write_gempak
  use diagnostics,  only: vint, hydro_interp, mix_ratio_to_spec_hum,          &
                          temp_from_theta_p, surface_pres, density, sea_pres, &
                          lifted_index, mean_rh, unstag_hght, cape,           &
                          sfcape => sbcape, mlcape, lcl_lfc, storm_motion,    &
                          storm_rel_hel, calcdbz, theta_phi, pvor
  use diagnostics,  only: Grav, Eps
  use maths,        only: mean, cen_diff_stag, cen_diff, uneven_deriv,        &
                          calc_mfpsi
  use maths,        only: Half, One
  use partition_wind, only: partition, balance_lhs, balance_lhs_exp,          &
                            balance_lhs_tlm, balance_lhs_adj, balance_rhs_simp,  &
                            balance_rhs_simp_exp, balance_rhs_simp_tlm, balance_rhs_simp_adj,  &
                            imbalance_simp_adj, pv_exp, pv_tlm, pv_adj, imbalance_exp,  &
                            solve_phi_psi
  use current_kind, only: cp, Quarter
  implicit none
  save
  private
  public :: init_reg, output_var

  ! The currently programmed WRF variables are:
  ! 1  - HGT: Terrain height (m)
  ! 2  - T2: 2-m temperature (K)
  ! 3  - Q2: 2-m water vapor mixing ratio (kg/kg)
  ! 4  - ZNU: Eta coordinate on mass (half) levels
  ! 5  - PB: Base state pressure (Pa)
  ! 6  - P: Perturbation pressure (Pa)
  ! 7  - T: Perterbation potential temperature (K)
  ! 8  - U: U component of wind (m/s)
  ! 9  - V: V component of wind (m/s)
  ! 10 - W: W component of wind (m/s)
  ! 11 - U10: U at 10 m (m/s)
  ! 12 - V10: V at 10 m (m/s)
  ! 13 - QVAPOR: Water vapor mixing ratio (kg/kg)
  ! 14 - QCLOUD: Cloud water mixing ratio (kg/kg)
  ! 15 - QICE: Ice mixing ratio (kg/kg)
  ! 16 - RAINC: Accumulated cumulus precip (mm)
  ! 17 - RAINNC: Accumulated grid-scale precip (mm)
  ! 18 - PH: Perturbation geopotential (m^2/s^2)
  ! 19 - PHB: Base-state geopotential (m^2/s^2)
  ! 20 - ZNW: Eta coordinate on w (full) levels
  ! 21 - LU_INDEX: Land use category
  ! 22 - MU: Perturbation dry air mass in column (Pa)
  ! 23 - MUB: Base-state dry air mass in column (Pa)
  ! 24 - TH2: 2-m potential temperature (K)
  ! 25 - QRAIN: Rain water mixing ratio (kg/kg)
  ! 26 - QSNOW: Snow mixing ratio (kg/kg)
  ! 27 - QGRAUP: Graupel mixing ratio (kg/kg)
  ! 28 - PBLH: PBL height (m)
  ! 29 - SST: Sea surface temperature (K)
  ! 30 - TSK: Surface skin temperature (K)
  ! 31 - SWDOWN: Downward short wave flux at ground surface (W/m^2)
  ! 32 - GLW: Downward long wave flux at ground surface (W/m^2)
  ! 33 - GRDFLX: Ground heat flux (W/m^2)
  ! 34 - HFX: Upward heat flux at the surface (W/m^2)
  ! 35 - QFX: Upward moisture flux at the surface (W/m^2)
  ! 36 - LH: Latent heat flux at the surface (W/m^2)
  ! 37 - MAPFAC_M: Map scale factor on mass grid
  ! 38 - MAPFAC_U: Map scale factor on u-grid
  ! 39 - MAPFAC_V: Map scale factor on v-grid
  ! 40 - F: Coriolis sine latitude term (s^-1)

  ! The currently programmed diagnostics are:
  ! 101 - Surface pressure (hPa)
  ! 102 - Pressure at 2 m (hPa)
  ! 103 - Sea level pressure (hPa)
  ! 104 - Convective precipitation accumulated over one hour (mm)
  ! 105 - Total precipitation accumulated over one hour (mm)
  ! 106 - Total precipitation accumulated over three hours (mm)
  ! 107 - Total precipitation accumulated over six hours (mm)
  ! 108 - Total accumulated precipitation (mm)
  ! 109 - Pressure (hPa)
  ! 110 - Potential temperature (K)
  ! 111 - Unstaggered U component of wind (m/s)
  ! 112 - Unstaggered V component of wind (m/s)
  ! 113 - Cloud water + Ice mixing ratio (kg/kg)
  ! 114 - Unstaggered W component of wind (m/s)
  ! 115 - 2-m temperature including at initial time (K)
  ! 116 - 2-m water vapor mixing ratio including at initial time (kg/kg)
  ! 117 - U at 10 m including at initial time (m/s)
  ! 118 - V at 10 m including at initial time (m/s)
  ! 119 - Lifted Index using lowest 4 levels (K)
  ! 120 - Mean relative humidity (%) in 850-500-hPa layer
  ! 121 - Geopotential height (m) on w (full) levels
  ! 122 - 2-m specific humidity (kg/kg)
  ! 123 - Specific humidity (kg/kg)
  ! 124 - Temperature (K)
  ! 125 - Geopotential height (m) on half levels
  ! 126 - Dry air mass in column (Pa)
  ! 127 - Most unstable CAPE w/ virt. temp. correction (J/kg)
  ! 128 - Lifted parcel level of the most unstable parcel (m)
  ! 129 - Surface-based CAPE w/ virt. temp. correction (J/kg)
  ! 130 - Surface-based CIN w/ virt. temp. correction (J/kg)
  ! 131 - Mixed-layer CAPE w/ virt. temp. correction (J/kg)
  ! 132 - Mixed-layer CIN w/ virt. temp. correction (J/kg)
  ! 133 - Mixed-layer LCL (m)
  ! 134 - Mixed-layer LFC (m)
  ! 135 - U component of Bunkers storm motion (m/s)
  ! 136 - V component of Bunkers storm motion (m/s)
  ! 137 - 0-1-km storm relative helicity (m2/s2)
  ! 138 - Column integrated precipitable water (cm)
  ! 139 - Reflectivity (dBZ)
  ! 140 - Composite reflectivity (dBZ)

  ! The currently programmed diagnostics interpolated to pressure are:
  ! 501 - Temperature (K)
  ! 502 - U component of wind (m/s)
  ! 503 - V component of wind (m/s)
  ! 504 - W component of wind (m/s)
  ! 505 - Geopotential height (m)
  ! 506 - Specific humidity (kg/kg)

  ! Encode GEMPAK vertical coordinate IDs
  integer, parameter :: HAGL = 1279738184, None = 0, Sgma = 4, Pres = 1

  ! Set up enumeration for the available variables
  integer, parameter :: LastWRF = 40, FirstDiag = 101, LastDiag = 260,  &
                        FirstPresDiag = 501, LastPresDiag = 506
  integer, parameter :: TerHght   = 1,   SfcTemp   = 2,   SfcMixR   = 3,  &
                        EtaMass   = 4,   PBase     = 5,   PPert     = 6,  &
                        ThetaPert = 7,   UStag     = 8,   VStag     = 9,  &
                        WStag     = 10,  SfcU      = 11,  SfcV      = 12,  &
                        MixR      = 13,  MixRCloud = 14,  MixRIce   = 15,  &
                        ConvPre   = 16,  StabPre   = 17,  GeopPert  = 18,  &
                        GeopBase  = 19,  EtaW      = 20,  LandUse   = 21,  &
                        MuPert    = 22,  MuBase    = 23,  SfcTheta  = 24,  &
                        MixRRain  = 25,  MixRSnow  = 26,  MixRGraup = 27,  &
                        PBLHgt    = 28,  SeaTemp   = 29,  SkinTemp  = 30,  &
                        SWGrd     = 31,  LWGrd     = 32,  GrdFx     = 33,  &
                        HeatFx    = 34,  MoisFx    = 35,  LaHeat    = 36,  &
                        MapFacM   = 37,  MapFacU   = 38,  MapFacV   = 39,  &
                        Coriolis  = 40
  integer, parameter :: SfcP      = 101, SfcP2m    = 102, SeaLevP   = 103,  &
                        ConvPre1h = 104, TotPre1h  = 105, TotPre3h  = 106,  &
                        TotPre6h  = 107, TotPre    = 108, Pressure  = 109,  &
                        Theta     = 110, UUnstag   = 111, VUnstag   = 112,  &
                        CWtr      = 113, WUnstag   = 114, SfcTempF0 = 115,  &
                        SfcMixRF0 = 116, SfcUF0    = 117, SfcVF0    = 118,  &
                        LiftIdx   = 119, MeanRH    = 120, GeopHghtS = 121,  &
                        SfcSpHum  = 122, SpHum     = 123, TmpK      = 124,  &
                        GeopHght  = 125, DryAir    = 126, MUCAPE    = 127,  &
                        MULPL     = 128, SBCAPE    = 129, SBCIN     = 130,  &
                        MixCAPE   = 131, MixCIN    = 132, MixLCL    = 133,  &
                        MixLFC    = 134, UStorm    = 135, VStorm    = 136,  &
                        StormRHel = 137, PreWater  = 138, Refl      = 139,  &
                        CompRefl  = 140
  integer, parameter :: MCap      = 201, MScript1  = 202, MScript2  = 203,  &
                        Diverg    = 204, DScript   = 205, DDDt      = 206,  &
                        EtaDot    = 207, TermA2    = 208, TermA3    = 209,  &
                        TermA4    = 210, TermA5    = 211, TermB1    = 212,  &
                        TermB2    = 213, TermB3    = 214, TermC     = 215,  &
                        TermD     = 216, MS1Unstag = 217, MS2Unstag = 218,  &
                        TermE     = 219, B1B2CD    = 220, ABCD      = 221,  &
                        Stream    = 222, VelocPot  = 223, DivPsi    = 224,  &
                        MS1PsiUnstag = 225, MS2PsiUnstag = 226,  &
                        MS1ChiUnstag = 227, MS2ChiUnstag = 228,  &
                        BalLHS    = 229, BalRHSSimp= 230, ImbalSimp = 231,  &
                        BalLHSExp = 232, BalLHSTLM = 233, BalLHSPert= 234,  &
                        BalLHSAdj     = 235, BalRHSSimpExp= 236, BalRHSSimpTLM= 237,  &
                        BalRHSSimpPert= 238, BalRHSSimpAdj= 239, ImbalSimpPert= 240,  &
                        ImbalSimpTLM  = 241, ImbalSimpAdj = 242, ThetaFromPhi = 243,  &
                        PotVort   = 244, PVExp     = 245, PVTLM     = 246,  &
                        PVPert    = 247, PVAdj     = 248, CostGradPhi = 249,  &
                        InvertPhi = 250, InvertPsi = 251

  integer, parameter :: TmpKPres     = 501, UUnstagPres  = 502,  &
                        VUnstagPres  = 503, WUnstagPres  = 504,  &
                        GeopHghtPres = 505, SpHumPres    = 506

  ! This structure contains information about the various WRF variables.
  type wrf_var
     integer :: dimType,  &       ! 0-z,t; 1-x,y,t; 2-x,y,z,t; 3-xs,y,z,t;
                                  ! 4-x,ys,z,t; 5-x,y,zs,t
                lev1, lev2,  &    ! GEMPAK level parameters. Relevant only for
                                  ! dimType 1
                vertCordID        ! GEMPAK vertical coordinate ID
     character(len=8) :: wrfName  ! WRF variable name
     character(len=4) :: gemName  ! GEMPAK variable name
     logical          :: const    ! T - Fixed with time; F - Varies with time
  end type wrf_var

  ! This is an analogous structure for diagnosed variables.
  type diag_var
     integer           :: dimType, lev1, lev2, vertCordID, diagID
     character(len=12) :: gemName
  end type diag_var

  ! And finally a structure for diagnosed variables interpolated to pressure.
  type pres_diag_var
     integer           :: interpType,  &  ! 1-linear interp; 2-log interp
                          diagID
     character(len=12) :: gemName
  end type pres_diag_var

  type(wrf_var),  dimension(LastWRF)                         :: wrfVar
  type(diag_var), dimension(FirstDiag:LastDiag)              :: diagVar
  type(pres_diag_var), dimension(FirstPresDiag:LastPresDiag) :: presDiagVar
  integer :: timesPer12hr, pBot, pTop, dp, np
  
contains  ! ===================================================================

  ! init_reg initializes the registry, which contains the vital parameters for
  ! each variable or diagnostic, with the number of outputs per 12 hours (freq)
  ! and pressure interpolation information (pb, pt, deltap) provided as input.
  subroutine init_reg(freq, pb, pt, deltap)
    integer, intent(in) :: freq, pb, pt, deltap

    wrfVar(TerHght)   = wrf_var(1, 0, -1, None, "HGT", "HGHT", .true.)
    wrfVar(SfcTemp)   = wrf_var(1, 2, -1, HAGL, "T2", "TMPK", .false.)
    wrfVar(SfcMixR)   = wrf_var(1, 2, -1, HAGL, "Q2", "MIXR", .false.)
    wrfVar(EtaMass)   = wrf_var(0, -1, -1, -1, "ZNU", "SGMA", .true.)
    wrfVar(PBase)     = wrf_var(2, -1, -1, Sgma, "PB", "PBAR", .true.)
    wrfVar(PPert)     = wrf_var(2, -1, -1, Sgma, "P", "PP", .false.)
    wrfVar(ThetaPert) = wrf_var(2, -1, -1, Sgma, "T", "TMPO", .false.)
    wrfVar(UStag)     = wrf_var(3, -1, -1, Sgma, "U", "USTG", .false.)
    wrfVar(VStag)     = wrf_var(4, -1, -1, Sgma, "V", "VSTG", .false.)
    wrfVar(WStag)     = wrf_var(5, -1, -1, Sgma, "W", "WSTG", .false.)
    wrfVar(SfcU)      = wrf_var(1, 10, -1, HAGL, "U10", "UREL", .false.)
    wrfVar(SfcV)      = wrf_var(1, 10, -1, HAGL, "V10", "VREL", .false.)
    wrfVar(MixR)      = wrf_var(2, -1, -1, Sgma, "QVAPOR", "MIXR", .false.)
    wrfVar(MixRCloud) = wrf_var(2, -1, -1, Sgma, "QCLOUD", "WC", .false.)
    wrfVar(MixRIce)   = wrf_var(2, -1, -1, Sgma, "QICE", "WI", .false.)
    wrfVar(ConvPre)   = wrf_var(1, 0, -1, None, "RAINC", "CTOT", .false.)
    wrfVar(StabPre)   = wrf_var(1, 0, -1, None, "RAINNC", "ETOT", .false.)
    wrfVar(GeopPert)  = wrf_var(5, -1, -1, Sgma, "PH", "HP", .false.)
    wrfVar(GeopBase)  = wrf_var(5, -1, -1, Sgma, "PHB", "HBAR", .true.)
    wrfVar(EtaW)      = wrf_var(0, -1, -1, -1, "ZNW", "SGMA", .true.)
    wrfVar(LandUse)   = wrf_var(1, 0, -1, None, "LU_INDEX", "LUC", .true.)
    wrfVar(MuPert)    = wrf_var(1, 0, -1, None, "MU", "MUP", .false.)
    wrfVar(MuBase)    = wrf_var(1, 0, -1, None, "MUB", "MUB", .true.)
    wrfVar(SfcTheta)  = wrf_var(1, 2, -1, HAGL, "TH2", "THTA", .false.)
    wrfVar(MixRRain)  = wrf_var(2, -1, -1, Sgma, "QRAIN", "WR", .false.)
    wrfVar(MixRSnow)  = wrf_var(2, -1, -1, Sgma, "QSNOW", "WS", .false.)
    wrfVar(MixRGraup) = wrf_var(2, -1, -1, Sgma, "QGRAUP", "WG", .false.)
    wrfVar(PBLHgt)    = wrf_var(1, 0, -1, None, "PBLH", "PBLH", .false.)
    wrfVar(SeaTemp)   = wrf_var(1, 0, -1, None, "SST", "SST", .false.)
    wrfVar(SkinTemp)  = wrf_var(1, 0, -1, None, "TSK", "TSK", .false.)
    wrfVar(SWGrd)     = wrf_var(1, 0, -1, None, "SWDOWN", "SWD", .false.)
    wrfVar(LWGrd)     = wrf_var(1, 0, -1, None, "GLW", "LWD", .false.)
    wrfVar(GrdFx)     = wrf_var(1, 0, -1, None, "GRDFLX", "GRDF", .false.)
    wrfVar(HeatFx)    = wrf_var(1, 0, -1, None, "HFX", "HFX", .false.)
    wrfVar(MoisFx)    = wrf_var(1, 0, -1, None, "QFX", "QFX", .false.)
    wrfVar(LaHeat)    = wrf_var(1, 0, -1, None, "LH", "LH", .false.)
    wrfVar(MapFacM)   = wrf_var(1, 0, -1, None, "MAPFAC_M", "MFM", .true.)
    wrfVar(MapFacU)   = wrf_var(3, 0, -1, None, "MAPFAC_U", "MFU", .true.)
    wrfVar(MapFacV)   = wrf_var(4, 0, -1, None, "MAPFAC_V", "MFV", .true.)
    wrfVar(Coriolis)  = wrf_var(1, 0, -1, None, "F", "F", .true.)    

    diagVar(SfcP)      = diag_var(1, 0, -1, None, SfcP, "PRES")
    diagVar(SfcP2m)    = diag_var(1, 2, -1, HAGL, SfcP2m, "PRES")
    diagVar(SeaLevP)   = diag_var(1, 0, -1, None, SeaLevP, "PMSL")
    diagVar(ConvPre1h) = diag_var(1, 0, -1, None, ConvPre1h, "C01M")
    diagVar(TotPre1h)  = diag_var(1, 0, -1, None, TotPre1h, "P01M")
    diagVar(TotPre3h)  = diag_var(1, 0, -1, None, TotPre3h, "P03M")
    diagVar(TotPre6h)  = diag_var(1, 0, -1, None, TotPre6h, "P06M")
    diagVar(TotPre)    = diag_var(1, 0, -1, None, TotPre, "PTOT")
    diagVar(Pressure)  = diag_var(2, -1, -1, Sgma, Pressure, "PRES")
    diagVar(Theta)     = diag_var(2, -1, -1, Sgma, Theta, "THTA")
    diagVar(UUnstag)   = diag_var(2, -1, -1, Sgma, UUnstag, "UREL")
    diagVar(VUnstag)   = diag_var(2, -1, -1, Sgma, VUnstag, "VREL")
    diagVar(CWtr)      = diag_var(2, -1, -1, Sgma, CWtr, "CWTR")
    diagVar(WUnstag)   = diag_var(2, -1, -1, Sgma, WUnstag, "W")
    diagVar(SfcTempF0) = diag_var(1, 2, -1, HAGL, SfcTempF0, "TMPK")
    diagVar(SfcMixRF0) = diag_var(1, 2, -1, HAGL, SfcMixRF0, "MIXR")
    diagVar(SfcUF0)    = diag_var(1, 10, -1, HAGL, SfcUF0, "UREL")
    diagVar(SfcVF0)    = diag_var(1, 10, -1, HAGL, SfcVF0, "VREL")
    diagVar(LiftIdx)   = diag_var(1, 0, -1, None, LiftIdx, "LFT4")
    diagVar(MeanRH)    = diag_var(1, 850, 500, Pres, MeanRH, "RELH")
    diagVar(GeopHghtS) = diag_var(5, -1, -1, Sgma, GeopHghtS, "HGHT")
    diagVar(SfcSpHum)  = diag_var(1, 2, -1, HAGL, SfcSpHum, "SPFH")
    diagVar(SpHum)     = diag_var(2, -1, -1, Sgma, SpHum, "SPFH")
    diagVar(TmpK)      = diag_var(2, -1, -1, Sgma, TmpK, "TMPK")
    diagVar(GeopHght)  = diag_var(2, -1, -1, Sgma, GeopHght, "HGHT")
    diagVar(DryAir)    = diag_var(1, 0, -1, None, DryAir, "MU")
    diagVar(MUCAPE)    = diag_var(1, 0, -1, None, MUCAPE, "MUCAPE")
    diagVar(MULPL)     = diag_var(1, 0, -1, None, MULPL, "LPL")
    diagVar(SBCAPE)    = diag_var(1, 0, -1, None, SBCAPE, "SBCAPE")
    diagVar(SBCIN)     = diag_var(1, 0, -1, None, SBCIN, "SBCIN")
    diagVar(MixCAPE)   = diag_var(1, 0, -1, None, MixCAPE, "MLCAPE")
    diagVar(MixCIN)    = diag_var(1, 0, -1, None, MixCIN, "MLCIN")
    diagVar(MixLCL)    = diag_var(1, 0, -1, None, MixLCL, "LCL")
    diagVar(MixLFC)    = diag_var(1, 0, -1, None, MixLFC, "LFC")
    diagVar(UStorm)    = diag_var(1, 0, -1, None, UStorm, "USTM")
    diagVar(VStorm)    = diag_var(1, 0, -1, None, VStorm, "VSTM")
    diagVar(StormRHel) = diag_var(1, 0, 1000, HAGL, StormRHel, "SRH")
    diagVar(PreWater)  = diag_var(1, 0, -1, None, PreWater, "PW")
    diagVar(Refl)      = diag_var(2, -1, -1, Sgma, Refl, "REFL")
    diagVar(CompRefl)  = diag_var(1, 0, -1, None, CompRefl, "CREF")

    diagVar(MCap)      = diag_var(1, 0, -1, None, MCap, "MCAP")
    diagVar(MScript1)  = diag_var(3, -1, -1, Sgma, MScript1, "MSCR1")
    diagVar(MScript2)  = diag_var(4, -1, -1, Sgma, MScript2, "MSCR2")
    diagVar(Diverg)    = diag_var(2, -1, -1, Sgma, Diverg, "DIVG")
    diagVar(DScript)   = diag_var(2, -1, -1, Sgma, DScript, "DSCR")
    diagVar(DDDt)      = diag_var(2, -1, -1, Sgma, DDDt, "DDDT")    
    diagVar(EtaDot)    = diag_var(2, -1, -1, Sgma, EtaDot, "EDOT")
    diagVar(TermA2)    = diag_var(2, -1, -1, Sgma, TermA2, "A2")    
    diagVar(TermA3)    = diag_var(2, -1, -1, Sgma, TermA3, "A3")
    diagVar(TermA4)    = diag_var(2, -1, -1, Sgma, TermA4, "A4")
    diagVar(TermA5)    = diag_var(2, -1, -1, Sgma, TermA5, "A5")    
    diagVar(TermB1)    = diag_var(2, -1, -1, Sgma, TermB1, "B1")
    diagVar(TermB2)    = diag_var(2, -1, -1, Sgma, TermB2, "B2")    
    diagVar(TermB3)    = diag_var(2, -1, -1, Sgma, TermB3, "B3")
    diagVar(TermC)     = diag_var(2, -1, -1, Sgma, TermC, "C")    
    diagVar(TermD)     = diag_var(2, -1, -1, Sgma, TermD, "D")
    diagVar(MS1Unstag) = diag_var(2, -1, -1, Sgma, MS1Unstag, "MSCR1")
    diagVar(MS2Unstag) = diag_var(2, -1, -1, Sgma, MS2Unstag, "MSCR2")
    diagVar(TermE)     = diag_var(2, -1, -1, Sgma, TermE, "E")    
    diagVar(B1B2CD)    = diag_var(2, -1, -1, Sgma, B1B2CD, "B1B2CD")    
    diagVar(ABCD)      = diag_var(2, -1, -1, Sgma, ABCD, "ABCD")    
    diagVar(Stream)    = diag_var(2, -1, -1, Sgma, Stream, "STRFUN")     
    diagVar(VelocPot)  = diag_var(2, -1, -1, Sgma, VelocPot, "VELPOT")
    diagVar(DivPsi)    = diag_var(2, -1, -1, Sgma, DivPsi, "DIVPSI")    
    diagVar(MS1PsiUnstag) = diag_var(2, -1, -1, Sgma, MS1PsiUnstag, "MS1PSI") 
    diagVar(MS2PsiUnstag) = diag_var(2, -1, -1, Sgma, MS2PsiUnstag, "MS2PSI")
    diagVar(MS1ChiUnstag) = diag_var(2, -1, -1, Sgma, MS1ChiUnstag, "MS1CHI") 
    diagVar(MS2ChiUnstag) = diag_var(2, -1, -1, Sgma, MS2ChiUnstag, "MS2CHI")
    diagVar(BalLHS)    = diag_var(2, -1, -1, Sgma, BalLHS, "LHS")
    diagVar(BalRHSSimp) = diag_var(2, -1, -1, Sgma, BalRHSSimp, "RHS")    
    diagVar(ImbalSimp) = diag_var(2, -1, -1, Sgma, ImbalSimp, "IMB")
    diagVar(BalLHSExp) = diag_var(2, -1, -1, Sgma, BalLHSExp, "LHSE")
    diagVar(BalLHSTLM) = diag_var(2, -1, -1, Sgma, BalLHSTLM, "LHSTLM")
    diagVar(BalLHSPert) = diag_var(2, -1, -1, Sgma, BalLHSPert, "LHSP")    
    diagVar(BalLHSAdj) = diag_var(2, -1, -1, Sgma, BalLHSAdj, "LHSADJ")    
    diagVar(BalRHSSimpExp) = diag_var(2, -1, -1, Sgma, BalRHSSimpExp, "RHSE")
    diagVar(BalRHSSimpTLM) = diag_var(2, -1, -1, Sgma, BalRHSSimpTLM, "RHSTLM")
    diagVar(BalRHSSimpPert) = diag_var(2, -1, -1, Sgma, BalRHSSimpPert, "RHSP")
    diagVar(BalRHSSimpAdj) = diag_var(2, -1, -1, Sgma, BalRHSSimpAdj, "RHSADJ")
    diagVar(ImbalSimpPert) = diag_var(2, -1, -1, Sgma, ImbalSimpPert, "IMBP")
    diagVar(ImbalSimpTLM) = diag_var(2, -1, -1, Sgma, ImbalSimpTLM, "IMBTLM")
    diagVar(ImbalSimpAdj) = diag_var(2, -1, -1, Sgma, ImbalSimpAdj, "IMBADJ")
    diagVar(ThetaFromPhi) = diag_var(2, -1, -1, Sgma, ThetaFromPhi, "THPHI")    
    diagVar(PotVort)   = diag_var(2, -1, -1, Sgma, PotVort, "PV")
    diagVar(PVExp)     = diag_var(2, -1, -1, Sgma, PVExp, "PVE")    
    diagVar(PVTLM)     = diag_var(2, -1, -1, Sgma, PVTLM, "PVTLM")    
    diagVar(PVPert)    = diag_var(2, -1, -1, Sgma, PVPert, "PVP")
    diagVar(PVAdj)     = diag_var(2, -1, -1, Sgma, PVAdj, "PVADJ")    
    diagVar(CostGradPhi) = diag_var(2, -1, -1, Sgma, CostGradPhi, "CGPHI")
    diagVar(InvertPhi) = diag_var(2, -1, -1, Sgma, InvertPhi, "PHII")
    diagVar(InvertPsi) = diag_var(2, -1, -1, Sgma, InvertPsi, "PSII")   

    presDiagVar(TmpKPres)     = pres_diag_var(1, TmpKPres, "TMPK")
    presDiagVar(UUnstagPres)  = pres_diag_var(1, UUnstagPres, "UREL")
    presDiagVar(VUnstagPres)  = pres_diag_var(1, VUnstagPres, "VREL")
    presDiagVar(WUnstagPres)  = pres_diag_var(1, WUnstagPres, "W")
    presDiagVar(GeopHghtPres) = pres_diag_var(2, GeopHghtPres, "HGHT")! 1 or 2?
    presDiagVar(SpHumPres)    = pres_diag_var(1, SpHumPres, "SPFH")

    timesPer12hr = freq

    pBot = pb
    pTop = pt
    dp = deltap
    np = (pBot - pTop) / dp + 1
  end subroutine init_reg

  ! ===========================================================================

  ! output_var writes a particular WRF variable or diagnostic to the GEMPAK
  ! file associated with gemid.  ncid contains the file handle for the WRF
  ! history file from which the WRF data is being read.  The WRF model
  ! dimensions are specified by nx, ny, nz.  var is the ID for the requested
  ! diagnostic.  ind is the current index of histFiles, so that histFiles(ind)
  ! is the pathname for the GEMPAK history file currently being read.  times 
  ! is an array of output times in GEMPAK form, which corresponds to all of the
  ! times in the history file.
  subroutine output_var(ncid, gemid, nx, ny, nz, var, ind, histFiles, times)
    integer, intent(in) :: ncid, gemid, nx, ny, nz, var, ind
    character(len=80), dimension(:), intent(in) :: histFiles
    character(len=15), dimension(:), intent(in) :: times

    character(len=*), parameter :: fmt1 = "(2x,a,i3,a)"

    write (*,fmt1) "Writing field ", var, " to GEMPAK..."

    ! Call the appropriate subroutine to produce the required output for GEMPAK
    select case (var)
    case (1:LastWRF)
       if (wrfVar(var)%dimType == 1 .or. wrfVar(var)%dimType == 2 .or.  &
            wrfVar(var)%dimType == 5) then
          call simple_out(ncid, gemid, nx, ny, nz, wrfVar(var), times)
       else
          print *, "It does not make sense to output this variable to GEMPAK. &
               &Skipping!"
       end if
    case (FirstDiag:LastDiag)
       if (diagVar(var)%dimType == 1 .or. diagVar(var)%dimType == 2 .or.  &
            diagVar(var)%dimType == 5) then
          call diag_out(ncid, gemid, nx, ny, nz, ind, histFiles,  &
               diagVar(var), times)
       else
          print *, "This diagnostic is not appropriate for output to GEMPAK."
          print *, "It's either staggered in x or y, or not a function of x &
               &and y."
          print *, "Skipping!"
       end if
    case (FirstPresDiag:LastPresDiag)
       call pres_diag_out(ncid, gemid, nx, ny, nz, presDiagVar(var), times)
    case default
       write (*,fmt1) "Variable ", var, " is unknown. Skipping!"
    end select
  end subroutine output_var

  ! ===========================================================================

  ! simple_out reads the variable var in from the WRF history file associated
  ! with ncid and writes it out to the GEMPAK file associated with gemid.
  ! Other inputs to the subroutine are the WRF model dimensions (nx,ny,nz) and
  ! the array of times in GEMPAK form to be output, which corresponds to all of
  ! the times in the history file.
  subroutine simple_out(ncid, gemid, nx, ny, nz, var, times)
    integer,                         intent(in) :: ncid, gemid, nx, ny, nz
    type(wrf_var),                   intent(in) :: var
    character(len=15), dimension(:), intent(in) :: times
     
    real, dimension(:,:,:,:), allocatable :: out4d
    real, dimension(:,:,:),   allocatable :: out3d
    real, dimension(:,:),     allocatable :: elvl
    integer                               :: kMax, eta, nt, status, t, k, lev
    character(len=20)                     :: gdat1
    character(len=12)                     :: parm

    ! Initialization
    select case (var%dimType)
    case (1)  ! x,y,t
       kMax = 0
    case (2)  ! x,y,z,t
       kMax = nz
       eta = EtaMass
    case (5)  ! x,y,zs,t
       kMax = nz + 1
       eta = EtaW
    case default
       ! GEMPAK only understands variables that depend on (unstaggered) x and y
       print *, "Invalid dimType in simple_out!"
       stop
    end select
    
    nt = size(times)
    parm = var%gemName
    
    ! Read variable, and, if necessary, the vertical levels it's at
    if (var%dimType == 1) then
       allocate (out3d(nx,ny,nt), stat=status)
       call mem_error(status, 1, "simple_out")
       call read_values(ncid, var%wrfName, out3d)
    else
       allocate (out4d(nx,ny,kMax,nt), elvl(kMax,nt), stat=status)
       call mem_error(status, 1, "simple_out")
       call read_values(ncid, var%wrfName, out4d)
       call read_values(ncid, wrfVar(eta)%wrfName, elvl)
    end if

    ! Write the variable to GEMPAK
    if (var%const) nt = 1  ! Output constant variables only at initial time
    do t = 1, nt
       gdat1 = times(t)
       if (var%dimType == 1) then
          call write_gempak(gemid, out3d(:,:,t), nx, ny, gdat1, var%lev1,  &
               var%lev2, var%vertCordID, parm)
       else
          do k = 1, kMax
             lev = nint(elvl(k,1) * 10000)
             call write_gempak(gemid, out4d(:,:,k,t), nx, ny, gdat1, lev,  &
                  var%lev2, var%vertCordID, parm)
          end do
       end if
    end do
        
    ! Clean up
    if (allocated(out3d)) then
       deallocate (out3d, stat=status)
    else
       deallocate (out4d, elvl, stat=status)
    end if
    call mem_error(status, 2, "simple_out")
  end subroutine simple_out
  
  ! ===========================================================================

  ! diag_out is a driver routine that produces the diagnostic var based on data
  ! from the WRF history file associated with ncid.  The diagnostic is written
  ! out to the GEMPAK file associated with gemid.  Other inputs to the
  ! subroutine are the WRF model dimensions (nx,ny,nz), the current index (ind)
  ! of histFiles, so that histFiles(ind) is the pathname for the GEMPAK history
  ! file currently being read, and the array of times in GEMPAK form to be
  ! output, which corresponds to all of the times in the history file.
  subroutine diag_out(ncid, gemid, nx, ny, nz, ind, histFiles, var, times)
    integer,                         intent(in) :: ncid, gemid, nx, ny, nz, ind
    character(len=80), dimension(:), intent(in) :: histFiles
    type(diag_var),                  intent(in) :: var
    character(len=15), dimension(:), intent(in) :: times
    
    real, dimension(:,:,:,:), allocatable :: out4d
    real, dimension(:,:,:),   allocatable :: out3d
    real, dimension(:,:),     allocatable :: elvl
    integer           :: kMax, eta, nt, status, prev, diff, t, k, lev
    character(len=20) :: gdat1
    logical           :: ok
    
    ! Initialization
    select case (var%dimType)
    case (1)  ! x,y,t
       kMax = 0
    case (2)  ! x,y,z,t
       kMax = nz
       eta = EtaMass
    case (5)  ! x,y,zs,t
       kMax = nz + 1
       eta = EtaW
    case default
       print *, "Invalid dimType in diag_out!"
       stop
    end select

    nt = size(times)
    ok = .true.
  
    if (var%dimType == 1) then
       allocate (out3d(nx,ny,nt), stat=status)
       call mem_error(status, 1, "diag_out")
    else
       allocate (out4d(nx,ny,kMax,nt), elvl(kMax,nt), stat=status)
       call mem_error(status, 1, "diag_out")
       call read_values(ncid, wrfVar(eta)%wrfName, elvl)
    end if

    ! Compute appropriate diagnostic
    select case (var%diagID)
    case (SfcP)
       out3d = diag_SfcP(ncid, nx, ny, nz, nt)
    case (SfcP2m)
       out3d = diag_SfcP2m(ncid, nx, ny, nz, nt)
    case (SeaLevP)
       out3d = diag_SeaLevP(ncid, nx, ny, nz, nt)
    case (Pressure)
       out4d = diag_Pressure(ncid, nx, ny, nz, nt)
    case (Theta)
       out4d = diag_Theta(ncid, nx, ny, nz, nt)
    case (UUnstag)
       out4d = diag_UUnstag(ncid, nx, ny, nz, nt)
    case (VUnstag)
       out4d = diag_VUnstag(ncid, nx, ny, nz, nt)
    case (CWtr)
       out4d = diag_CWtr(ncid, nx, ny, nz, nt)
    case (WUnstag)
       out4d = diag_WUnstag(ncid, nx, ny, nz, nt)
    case (SfcTempF0)
       out3d = diag_SfcTempF0(ncid, nx, ny, nz, nt)
    case (SfcMixRF0)
       out3d = diag_SfcMixRF0(ncid, nx, ny, nz, nt)
    case (SfcUF0)
       out3d = diag_SfcUF0(ncid, nx, ny, nz, nt)
    case (SfcVF0)
       out3d = diag_SfcVF0(ncid, nx, ny, nz, nt)
    case (LiftIdx)
       out3d = diag_LiftIdx(ncid, nx, ny, nz, nt)
    case (MeanRH)
       out3d = diag_MeanRH(ncid, nx, ny, nz, nt)
    case (GeopHghtS)
       out4d = diag_GeopHghtS(ncid, nx, ny, nz, nt)
    case (SfcSpHum)
       out3d = diag_SfcSpHum(ncid, nx, ny, nz, nt)
    case (SpHum)
       out4d = diag_SpHum(ncid, nx, ny, nz, nt)
    case (TmpK)
       out4d = diag_TmpK(ncid, nx, ny, nz, nt)
    case (GeopHght)
       out4d = diag_GeopHght(ncid, nx, ny, nz, nt)
    case (DryAir)
       out3d = diag_DryAir(ncid, nx, ny, nt)

    case (MCap)
       out3d = diag_MCap(ncid, nx, ny, nt)
    case (MScript1)
       out4d = diag_MScript1(ncid, nx, ny, nz, nt)
    case (MScript2)
       out4d = diag_MScript2(ncid, nx, ny, nz, nt)       
    case (Diverg)
       out4d = diag_Diverg(ncid, nx, ny, nz, nt)
    case (DScript)
       out4d = diag_DScript(ncid, nx, ny, nz, nt)
    case (DDDt)
       out4d = diag_DDDt(ncid, nx, ny, nz, nt)
    case (EtaDot)
       out4d = diag_EtaDot(ncid, nx, ny, nz, nt)       
    case (TermA2)
       out4d = diag_TermA2(ncid, nx, ny, nz, nt)              
    case (TermA3)
       out4d = diag_TermA3(ncid, nx, ny, nz, nt)       
    case (TermA4)
       out4d = diag_TermA4(ncid, nx, ny, nz, nt)
    case (TermA5)
       out4d = diag_TermA5(ncid, nx, ny, nz, nt)       
    case (TermB1)
       out4d = diag_TermB1(ncid, nx, ny, nz, nt)              
    case (TermB2)
       out4d = diag_TermB2(ncid, nx, ny, nz, nt)
    case (TermB3)
       out4d = diag_TermB3(ncid, nx, ny, nz, nt)       
    case (TermC)
       out4d = diag_TermC(ncid, nx, ny, nz, nt)       
    case (TermD)
       out4d = diag_TermD(ncid, nx, ny, nz, nt)       
    case (MS1Unstag)
       out4d = diag_MS1Unstag(ncid, nx, ny, nz, nt)       
    case (MS2Unstag)
       out4d = diag_MS2Unstag(ncid, nx, ny, nz, nt)       
    case (TermE)
       out4d = diag_TermE(ncid, nx, ny, nz, nt)        
    case (B1B2CD)
       out4d = diag_B1B2CD(ncid, nx, ny, nz, nt)
    case (ABCD)
       out4d = diag_ABCD(ncid, nx, ny, nz, nt)
    case (Stream)
       out4d = diag_Stream(ncid, nx, ny, nz, nt)
    case (VelocPot)
       out4d = diag_VelocPot(ncid, nx, ny, nz, nt)       
    case (DivPsi)
       out4d = diag_DivPsi(ncid, nx, ny, nz, nt)
    case (MS1PsiUnstag)
       out4d = diag_MS1PsiUnstag(ncid, nx, ny, nz, nt)
    case (MS2PsiUnstag)
       out4d = diag_MS2PsiUnstag(ncid, nx, ny, nz, nt)
    case (MS1ChiUnstag)
       out4d = diag_MS1ChiUnstag(ncid, nx, ny, nz, nt)
    case (MS2ChiUnstag)
       out4d = diag_MS2ChiUnstag(ncid, nx, ny, nz, nt)       
    case (BalLHS)
       out4d = diag_BalLHS(ncid, nx, ny, nz, nt)
    case (BalRHSSimp)
       out4d = diag_BalRHSSimp(ncid, nx, ny, nz, nt)       
    case (ImbalSimp)
       out4d = diag_ImbalSimp(ncid, nx, ny, nz, nt)     
    case (BalLHSExp)
       out4d = diag_BalLHSExp(ncid, nx, ny, nz, nt)
    case (BalLHSTLM)
       out4d = diag_BalLHSTLM(ncid, nx, ny, nz, nt)
    case (BalLHSPert)
       out4d = diag_BalLHSPert(ncid, nx, ny, nz, nt)
    case (BalLHSAdj)
       out4d = diag_BalLHSAdj(ncid, nx, ny, nz, nt)      
    case (BalRHSSimpExp)
       out4d = diag_BalRHSSimpExp(ncid, nx, ny, nz, nt)
    case (BalRHSSimpTLM)
       out4d = diag_BalRHSSimpTLM(ncid, nx, ny, nz, nt)
    case (BalRHSSimpPert)
       out4d = diag_BalRHSSimpPert(ncid, nx, ny, nz, nt)
    case (BalRHSSimpAdj)
       out4d = diag_BalRHSSimpAdj(ncid, nx, ny, nz, nt)          
    case (ImbalSimpPert)
       out4d = diag_ImbalSimpPert(ncid, nx, ny, nz, nt)    
    case (ImbalSimpTLM)
       out4d = diag_ImbalSimpTLM(ncid, nx, ny, nz, nt)
    case (ImbalSimpAdj)
       out4d = diag_ImbalSimpAdj(ncid, nx, ny, nz, nt)       
    case (ThetaFromPhi)
       out4d = diag_ThetaFromPhi(ncid, nx, ny, nz, nt)       
    case (PotVort)
       out4d = diag_PotVort(ncid, nx, ny, nz, nt)          
    case (PVExp)
       out4d = diag_PVExp(ncid, nx, ny, nz, nt)   
    case (PVTLM)
       out4d = diag_PVTLM(ncid, nx, ny, nz, nt)   
    case (PVPert)
       out4d = diag_PVPert(ncid, nx, ny, nz, nt)   
    case (PVAdj)
       out4d = diag_PVAdj(ncid, nx, ny, nz, nt)      
    case (CostGradPhi)
       out4d = diag_CostGradPhi(ncid, nx, ny, nz, nt)       
    case (InvertPhi)
       out4d = diag_InvertPhi(ncid, nx, ny, nz, nt)       
    case (InvertPsi)
       out4d = diag_InvertPsi(ncid, nx, ny, nz, nt)       

    case default
       write (*,"(a,i3,a)") "WARNING: Field ", var%diagID,  &
            " unknown in diag_out!"
    end select

    ! Write the diagnostic out to GEMPAK
    if (ok) then
       do t = 1, nt
          gdat1 = times(t)
          if (var%dimType == 1) then
             call write_gempak(gemid, out3d(:,:,t), nx, ny, gdat1, var%lev1,  &
                  var%lev2, var%vertCordID, var%gemName)
          else
             do k = 1, kMax
                lev = nint(elvl(k,1) * 10000)
                call write_gempak(gemid, out4d(:,:,k,t), nx, ny, gdat1, lev,  &
                     var%lev2, var%vertCordID, var%gemName)
             end do
          end if
       end do
    end if

    ! Clean up
    if (allocated(out3d)) then
       deallocate (out3d, stat=status)
    else
       deallocate (out4d, elvl, stat=status)
    end if
    call mem_error(status, 2, "diag_out")
  end subroutine diag_out

  ! ===========================================================================

  ! pres_diag_out is a driver routine that produces the diagnostic var
  ! interpolated to pressure surfaces based on data from the WRF history file
  ! associated with ncid.  The diagnostic is written out to the GEMPAK file
  ! associated with gemid.  Other inputs to the subroutine are the WRF model
  ! dimensions (nx,ny,nz) and the array of times in GEMPAK form to be output,
  ! which corresponds to all of the times in the history file.  np is obtained
  ! through use association.
  subroutine pres_diag_out(ncid, gemid, nx, ny, nz, var, times)
    integer,                         intent(in) :: ncid, gemid, nx, ny, nz
    type(pres_diag_var),             intent(in) :: var
    character(len=15), dimension(:), intent(in) :: times

    real, dimension(:,:,:,:), allocatable :: tk
    real, dimension(nx,ny,nz+1,size(times)) :: outNoInterp, p
    real, dimension(nx,ny,np,size(times))   :: outInterp
    integer                                 :: nt, t, j, i, status, k, lev
    character(len=20)                       :: gdat1
    
    nt = size(times)
    
    ! Compute appropriate diagnostic
    select case (var%diagID)
    case (TmpKPres)
       outNoInterp(:,:,1,:) = diag_SfcTempF0(ncid, nx, ny, nz, nt)
       outNoInterp(:,:,2:nz+1,:) = diag_TmpK(ncid, nx, ny, nz, nt)
    case (UUnstagPres)
       outNoInterp(:,:,1,:) = diag_SfcUF0(ncid, nx, ny, nz, nt)
       outNoInterp(:,:,2:nz+1,:) = diag_UUnstag(ncid, nx, ny, nz, nt)
    case (VUnstagPres)
       outNoInterp(:,:,1,:) = diag_SfcVF0(ncid, nx, ny, nz, nt)
       outNoInterp(:,:,2:nz+1,:) = diag_VUnstag(ncid, nx, ny, nz, nt)
    case (WUnstagPres)
       outNoInterp(:,:,1,:) = 0
       outNoInterp(:,:,2:nz+1,:) = diag_WUnstag(ncid, nx, ny, nz, nt)
    case (GeopHghtPres)
       call read_values(ncid, wrfVar(TerHght)%wrfName, outNoInterp(:,:,1,:))
       outNoInterp(:,:,2:nz+1,:) = diag_GeopHght(ncid, nx, ny, nz, nt)
    case (SpHumPres)
       outNoInterp(:,:,1,:) = diag_SfcSpHum(ncid, nx, ny, nz, nt)
       outNoInterp(:,:,2:nz+1,:) = diag_SpHum(ncid, nx, ny, nz, nt)
    case default
       write (*,"(a,i3,a)") "WARNING: Field ", var,  &
            " unknown in pres_diag_out!"
    end select
    
    ! Need pressure for interpolation
    p(:,:,1,:) = diag_SfcP(ncid, nx, ny, nz, nt)
    p(:,:,2:nz+1,:) = diag_Pressure(ncid, nx, ny, nz, nt)

    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             call vint(outNoInterp(i,j,:,t), p(i,j,:,t), pBot, -dp,  &
                  var%interpType, outInterp(i,j,:,t))
          end do
       end do
    end do

    ! Interpolate underground for heights
    if (var%diagID == GeopHghtPres) then
       allocate (tk(nx,ny,nz,nt), stat=status)
       call mem_error(status, 1, "pres_diag_out")
       tk = diag_TmpK(ncid, nx, ny, nz, nt)

       call hydro_interp(outNoInterp(:,:,1,:), p(:,:,1,:), tk(:,:,1,:), pBot, &
            -dp, outInterp)
       
       deallocate (tk, stat=status)
       call mem_error(status, 2, "pres_diag_out")
    end if
    
    ! Write the diagnostic out to GEMPAK
    do t = 1, nt
       gdat1 = times(t)
       do k = 1, np
          lev = pBot - (k-1) * dp
          call write_gempak(gemid, outInterp(:,:,k,t), nx, ny, gdat1, lev,   &
               -1, Pres, var%gemName)
       end do
    end do
  end subroutine pres_diag_out

  ! ===========================================================================

  ! The next group of functions produce the various diagnostics.  They all take
  ! as input the WRF history file handle (ncid), the WRF domain size 
  ! (nx, ny, and sometimes nz), and the number of times contained in the
  ! history file (nt).

  ! diag_SfcP diagnoses the surface pressure (in hPa).
  function diag_SfcP(ncid, nx, ny, nz, nt) result (ps)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: ps

    character(len=9), parameter :: funName = "diag_SfcP"

    real, dimension(:,:,:,:), allocatable :: p, w, t, ht
    real, dimension(nx,ny,2,nt)           :: ht12
    real, dimension(nx,ny,nt)             :: p1, q1, sfcq, t1
    integer                               :: status
    
    ! Read in pressure, saving only lowest level
    allocate (p(nx,ny,nz,nt), stat=status)
    call mem_error(status, 1, funName//" p")
    p = diag_Pressure(ncid, nx, ny, nz, nt, si=.true.)  ! Do calculations in Pa
    p1 = p(:,:,1,:)
    deallocate (p, stat=status)
    call mem_error(status, 2, funName//" p")

    ! Handle specific humidity
    allocate (w(nx,ny,nz,nt), stat=status)
    call mem_error(status, 1, funName//" w")
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    q1 = mix_ratio_to_spec_hum(w(:,:,1,:))
    deallocate (w, stat=status)
    call mem_error(status, 2, funName//" w")
    sfcq = mix_ratio_to_spec_hum(diag_SfcMixRF0(ncid, nx, ny, nz, nt))

    ! Handle temperature
    allocate (t(nx,ny,nz,nt), stat=status)
    call mem_error(status, 1, funName//" t")
    t = diag_Theta(ncid, nx, ny, nz, nt)
    t1 = temp_from_theta_p(t(:,:,1,:), p1)
    deallocate (t, stat=status)
    call mem_error(status, 2, funName//" t")

    ! Handle geopotential height
    allocate (ht(nx,ny,nz+1,nt), stat=status)
    call mem_error(status, 1, funName//" ht")
    ht = diag_GeopHghtS(ncid, nx, ny, nz, nt)
    ht12 = ht(:,:,1:2,:)
    deallocate (ht, stat=status)
    call mem_error(status, 2, funName//" ht")

    ps = .01 * surface_pres(ht12, p1, q1, t1, sfcq)  ! Convert to hPa
  end function diag_SfcP

  ! ===========================================================================

  ! diag_SfcP2m diagnoses the pressure (in hPa) 2 meters above the surface.
  ! The density is approximated by using surface pressure and 2-m temperature.
  function diag_SfcP2m(ncid, nx, ny, nz, nt) result (ps2m)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: ps, t2, q, rhog, ps2m
    
    ps = diag_SfcP(ncid, nx, ny, nz, nt) * 100  ! In Pa
    t2 = diag_SfcTempF0(ncid, nx, ny, nz, nt)
    q = diag_SfcMixRF0(ncid, nx, ny, nz, nt)
    rhog = Grav * density(ps, t2, q)

    ! Use hydrostatic equation
    ps2m = .01 * (ps - 2*rhog)  ! In hPa
  end function diag_SfcP2m

  ! ===========================================================================

  ! diag_SeaLevP diagnoses the sea level pressure (in hPa).
  function diag_SeaLevP(ncid, nx, ny, nz, nt) result (slp)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: slp
    
    real, dimension(nx,ny,nz+1,nt) :: ht
    real, dimension(nx,ny,nz,nt)   :: th, p, w
    
    ! Read in WRF data
    ht = diag_GeopHghtS(ncid, nx, ny, nz, nt)
    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt, si=.true.)  ! Do calculations in Pa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)  

    ! Compute sea level pressure
    slp = .01 * sea_pres(ht, th, p, mix_ratio_to_spec_hum(w))  ! Convert to hPa
  end function diag_SeaLevP

  ! ===========================================================================
  
  ! diag_Pressure diagnoses pressure. By default, the output is in hPa, but
  ! if si is true, the output is in Pa.
  function diag_Pressure(ncid, nx, ny, nz, nt, si) result (p)
    integer, intent(in)           :: ncid, nx, ny, nz, nt
    logical, intent(in), optional :: si
    real, dimension(nx,ny,nz,nt)  :: p

    real, dimension(nx,ny,nz,nt) :: pb
    logical                      :: hPa
    
    ! Set appropriate output unit flag
    if (present(si)) then
       hPa = .not. si
    else
       hPa = .true.
    end if

    ! Read in WRF data
    call read_values(ncid, wrfVar(PBase)%wrfName, pb)
    call read_values(ncid, wrfVar(PPert)%wrfName, p)

    ! Calculate pressure
    p = p + pb
    if (hPa) p = .01 * p
  end function diag_Pressure
  
  ! ===========================================================================

  ! diag_Theta diagnoses potential temperature (K).
  function diag_Theta(ncid, nx, ny, nz, nt) result (theta)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: theta!, q
    
    call read_values(ncid, wrfVar(ThetaPert)%wrfName, theta)
    theta = theta + 300  ! WRF stores theta as theta-300
  end function diag_Theta

  ! ===========================================================================

  ! diag_UUnstag diagnoses the wind component (m/s) in the x-coordinate
  ! direction at the mass points (grid box centers).
  function diag_UUnstag(ncid, nx, ny, nz, nt) result (u)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: u

    real, dimension(nx+1,ny,nz,nt) :: us

    ! Read staggered wind
    call read_values(ncid, wrfVar(UStag)%wrfName, us)

    ! Unstagger it
    u(1:nx,:,:,:) = .5 * (us(1:nx,:,:,:) + us(2:nx+1,:,:,:))
  end function diag_UUnstag
  
  ! ===========================================================================

  ! diag_VUnstag diagnoses the wind component (m/s) in the y-coordinate
  ! direction at the mass points (grid box centers).
  function diag_VUnstag(ncid, nx, ny, nz, nt) result (v)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: v

    real, dimension(nx,ny+1,nz,nt) :: vs

    ! Read staggered wind
    call read_values(ncid, wrfVar(VStag)%wrfName, vs)

    ! Unstagger it
    v(:,1:ny,:,:) = .5 * (vs(:,1:ny,:,:) + vs(:,2:ny+1,:,:))
  end function diag_VUnstag
  
  ! ===========================================================================

  ! diag_CWtr diagnoses the cloud water mixing ratio (kg/kg), including ice.
  function diag_CWtr(ncid, nx, ny, nz, nt) result (cw)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: cw

    real, dimension(nx,ny,nz,nt) :: cloud, ice
    
    call read_values(ncid, wrfVar(MixRCloud)%wrfName, cloud)
    call read_values(ncid, wrfVar(MixRIce)%wrfName, ice)
    cw = cloud + ice
  end function diag_CWtr

  ! ===========================================================================
  
  ! diag_WUnstag diagnoses the wind component (m/s) in the vertical direction
  ! at the mass points (grid box centers).
  function diag_WUnstag(ncid, nx, ny, nz, nt) result (w)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: w

    real, dimension(nx,ny,nz+1,nt) :: ws

    ! Read staggered wind
    call read_values(ncid, wrfVar(WStag)%wrfName, ws)

    ! Unstagger it
    w(:,:,1:nz,:) = .5 * (ws(:,:,1:nz,:) + ws(:,:,2:nz+1,:))
  end function diag_WUnstag

  ! ===========================================================================

  ! diag_GeopHghtS diagnoses geopotential height (m) on w (full) levels.
  function diag_GeopHghtS(ncid, nx, ny, nz, nt) result (height)
    integer, intent(in)            :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz+1,nt) :: height
    
    real, dimension(nx,ny,nz+1,nt) :: gpb, gpp
    
    call read_values(ncid, wrfVar(GeopBase)%wrfName, gpb)
    call read_values(ncid, wrfVar(GeopPert)%wrfName, gpp)
    height = (gpb + gpp) / Grav
  end function diag_GeopHghtS

  ! ===========================================================================

  ! diag_SfcTempF0 diagnoses 2-m temperature (K), including the initial time.
  function diag_SfcTempF0(ncid, nx, ny, nz, nt) result (temp)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: temp

    real, dimension(nx,ny,nz,nt) :: potTemp, p
    real, dimension(nx,ny)       :: tk1

    ! For initial time, calculate temperature in lowest grid box.
    potTemp = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt, si = .true.)  ! p in Pa
    tk1 = temp_from_theta_p(potTemp(:,:,1,1), p(:,:,1,1))

    ! Read in WRF 2-m temperature, but WRF doesn't include this variable at
    ! initial time, so assign the temperature in the lowest grid box to it.
    call read_values(ncid, wrfVar(SfcTemp)%wrfName, temp)
    temp(:,:,1) = tk1
  end function diag_SfcTempF0

  ! ===========================================================================

  ! diag_SfcMixRF0 diagnoses 2-m mixing ratio (kg/kg), including the initial
  ! time.
  function diag_SfcMixRF0(ncid, nx, ny, nz, nt) result (w2m)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: w2m

    real, dimension(nx,ny,nz,nt) :: w

    ! Read in WRF 2-m mixing ratio, but WRF doesn't include this variable at
    ! initial time, so assign the mixing ratio in the lowest grid box to it.
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    call read_values(ncid, wrfVar(SfcMixR)%wrfName, w2m)
    w2m(:,:,1) = w(:,:,1,1)
  end function diag_SfcMixRF0

  ! ===========================================================================

  ! diag_SfcUF0 diagnoses the 2-m wind component (m/s) in the x-coordinate
  ! direction, including the initial time.
  function diag_SfcUF0(ncid, nx, ny, nz, nt) result (u2m)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: u2m
    
    real, dimension(nx,ny,nz,nt) :: u
    
    ! Read in WRF 2-m u-component wind, but WRF doesn't include this variable
    ! at initial time, so assign the u-component wind in the lowest grid box to
    ! it.
    u = diag_UUnstag(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(SfcU)%wrfName, u2m)
    u2m(:,:,1) = u(:,:,1,1)
  end function diag_SfcUF0

  ! ===========================================================================

  ! diag_SfcVF0 diagnoses the 2-m wind component (m/s) in the y-coordinate
  ! direction, including the initial time.
  function diag_SfcVF0(ncid, nx, ny, nz, nt) result (v2m)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: v2m

    real, dimension(nx,ny,nz,nt) :: v
    
    ! Read in WRF 2-m v-component wind, but WRF doesn't include this variable
    ! at initial time, so assign the v-component wind in the lowest grid box to
    ! it.
    v = diag_VUnstag(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(SfcV)%wrfName, v2m)
    v2m(:,:,1) = v(:,:,1,1)
  end function diag_SfcVF0

  ! ===========================================================================

  ! diag_LiftIdx diagnoses the Lifted Index (K).
  function diag_LiftIdx(ncid, nx, ny, nz, nt) result (li)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: li

    real, dimension(nx,ny,nz,nt) :: th, p, w
    integer :: i, j, t

    ! Read in the necessary data from WRF output.
    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    
    ! Calculate lifted index, one profile at a time.
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             li(i,j,t) = lifted_index(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t))
          end do
       end do
    end do
  end function diag_LiftIdx

  ! ===========================================================================
  
  ! diag_MeanRH diagnoses the mean relative humidity (%) in the 850-500 mb
  ! layer.
  function diag_MeanRH(ncid, nx, ny, nz, nt) result (mrh)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: mrh

    real, dimension(nx,ny,nz,nt) :: th, p, w
    real, dimension(nx,ny,nt)    :: ps
    integer :: i, j, t
    
    ! Read in the necessary data from WRF output.
    th = diag_Theta(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt)  ! p in hPa
    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    ps = diag_SfcP(ncid, nx, ny, nz, nt)  ! ps in hPa
    
    ! Calculate mean RH, one profile at a time.
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             mrh(i,j,t) = mean_rh(th(i,j,:,t), p(i,j,:,t), w(i,j,:,t),  &
                  ps(i,j,t), 850., 500.)
          end do
       end do
    end do
  end function diag_MeanRH

  ! ===========================================================================

  ! diag_SfcSpHum diagnoses 2-m specific humidity (kg/kg).
  function diag_SfcSpHum(ncid, nx, ny, nz, nt) result (q2m)
    integer, intent(in)       :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nt) :: q2m

    q2m = mix_ratio_to_spec_hum(diag_SfcMixRF0(ncid, nx, ny, nz, nt))
  end function diag_SfcSpHum

  ! ===========================================================================
  
  ! diag_SpHum diagnoses specific humidity (kg/kg).
  function diag_SpHum(ncid, nx, ny, nz, nt) result (q)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: q

    real, dimension(nx,ny,nz,nt) :: w

    call read_values(ncid, wrfVar(MixR)%wrfName, w)
    q = mix_ratio_to_spec_hum(w)
  end function diag_SpHum

  ! ===========================================================================

  ! diag_TmpK diagnoses temperature (K).
  function diag_TmpK(ncid, nx, ny, nz, nt) result (tk)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: tk

    tk = temp_from_theta_p(diag_Theta(ncid, nx, ny, nz, nt),  &
         diag_Pressure(ncid, nx, ny, nz, nt, si=.true.))
  end function diag_TmpK

  ! ===========================================================================

  ! diag_GeopHght diagnoses geopotential height (m) on half (theta) levels.
  function diag_GeopHght(ncid, nx, ny, nz, nt) result (hght)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: hght

    real, dimension(:,:), allocatable :: etaTemp
    real, dimension(nx,ny,nz+1,nt)    :: hghtStag
    real, dimension(nx,ny,nt)         :: mu
    real, dimension(1,nt)             :: pt  ! really a constant
    real, dimension(nz+1)             :: znw
    real, dimension(nz)               :: znu
    integer :: status, i, j, t

    hghtStag = diag_GeopHghtS(ncid, nx, ny, nz, nt)
    mu = diag_DryAir(ncid, nx, ny, nt)

    allocate(etaTemp(nz+1,nt), stat=status)
    call mem_error(status, 1, "diag_GeopHght")
    call read_values(ncid, wrfVar(EtaW)%wrfName, etaTemp)
    znw = etaTemp(:,1)
    deallocate(etaTemp, stat=status)
    call mem_error(status, 2, "diag_GeopHght")

    allocate(etaTemp(nz,nt), stat=status)
    call mem_error(status, 1, "diag_GeopHght")
    call read_values(ncid, wrfVar(EtaMass)%wrfName, etaTemp)
    znu = etaTemp(:,1)
    deallocate(etaTemp, stat=status)
    call mem_error(status, 2, "diag_GeopHght")

    call read_values(ncid, "P_TOP", pt)

    ! Unstagger height by column.
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             hght(i,j,:,t) =  &
                  unstag_hght(hghtStag(i,j,:,t), znw, znu, mu(i,j,t), pt(1,1))
          end do
       end do
    end do
  end function diag_GeopHght

  ! ===========================================================================

  ! diag_DryAir diagnoses the dry air mass in the column (Pa).
  function diag_DryAir(ncid, nx, ny, nt) result (mu)
    integer, intent(in)       :: ncid, nx, ny, nt
    real, dimension(nx,ny,nt) :: mu
    
    real, dimension(nx,ny,nt) :: mub
    
    call read_values(ncid, wrfVar(MuBase)%wrfName, mub)
    call read_values(ncid, wrfVar(MuPert)%wrfName, mu)
    mu = mub + mu
  end function diag_dryAir

  ! ===========================================================================

  function diag_MCap(ncid, nx, ny, nt) result(mCap)
    integer, intent(in)       :: ncid, nx, ny, nt
    real, dimension(nx,ny,nt) :: mCap

    real, parameter :: p0 = 100000.

    real    :: pt

    call read_values_const(ncid, "P_TOP", pt)
    mCap = diag_dryAir(ncid, nx, ny, nt) / (p0 - pt)
!!$    print *, minval(mCap(:,:,13))
!!$    print *, mean(mCap(:,:,13))
!!$    print *, maxval(mCap(:,:,13))
  end function diag_MCap

  ! ===========================================================================

  function diag_MScript1(ncid, nx, ny, nz, nt) result(ms1)
    integer, intent(in)            :: ncid, nx, ny, nz, nt
    real, dimension(nx+1,ny,nz,nt) :: ms1

    real, dimension(nx+1,ny,nz,nt) :: us   
    real, dimension(nx,ny,nt)   :: mCap
    
    call read_values(ncid, wrfVar(UStag)%wrfName, us)
    mCap = diag_MCap(ncid, nx, ny, nt)
!    mCap = spread(diag_MCap(ncid, nx, ny, nt),3,nz)

    ms1(2:nx,:,:,:) = us(2:nx,:,:,:) * Half *  &
         spread(mCap(1:nx-1,:,:) + mCap(2:nx,:,:), 3, nz)
    ms1((/1,nx+1/),:,:,:) = us((/1,nx+1/),:,:,:) *  &
         spread(mCap((/1,nx/),:,:), 3, nz)
  end function diag_MScript1
  
  ! ===========================================================================

  function diag_MScript2(ncid, nx, ny, nz, nt) result(ms2)
    integer, intent(in)            :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny+1,nz,nt) :: ms2
    
    real, dimension(nx,ny+1,nz,nt) :: vs   
    real, dimension(nx,ny,nt)   :: mCap
    
    call read_values(ncid, wrfVar(VStag)%wrfName, vs)
    mCap = diag_MCap(ncid, nx, ny, nt)
!    mCap = spread(diag_MCap(ncid, nx, ny, nt),3,nz)

    ms2(:,2:ny,:,:) = vs(:,2:ny,:,:) * Half *  &
         spread(mCap(:,1:ny-1,:) + mCap(:,2:ny,:), 3, nz)
    ms2(:,(/1,ny+1/),:,:) = vs(:,(/1,ny+1/),:,:) *  &
         spread(mCap(:,(/1,ny/),:), 3, nz)
  end function diag_MScript2

  ! ===========================================================================

  function diag_MS1Unstag(ncid, nx, ny, nz, nt) result (ms1u)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: ms1u
   
    real, dimension(nx+1,ny,nz,nt) :: ms1

    ms1 = diag_MScript1(ncid, nx, ny, nz, nt)

    ms1u = Half * (ms1(:nx,:,:,:) + ms1(2:,:,:,:))
  end function diag_MS1Unstag

  ! ===========================================================================

  function diag_MS2Unstag(ncid, nx, ny, nz, nt) result (ms2u)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: ms2u
    
    real, dimension(nx,ny+1,nz,nt) :: ms2
    
    ms2 = diag_MScript2(ncid, nx, ny, nz, nt)

    ms2u = Half * (ms2(:,:ny,:,:) + ms2(:,2:,:,:))
  end function diag_MS2Unstag
  
  ! ===========================================================================

  function diag_Diverg(ncid, nx, ny, nz, nt) result(div)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: div

    real, dimension(nx,ny,nz,nt) :: u, v
    real, dimension(nx,ny)   :: dudx, dvdy, mfm
    real    :: rdx
    integer :: t, k

    u = diag_UUnstag(ncid, nx, ny, nz, nt)
    v = diag_VUnstag(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)
    call read_values_const(ncid, "RDX", rdx)
    
    do t = 1, nt
       do k = 1, nz
          dudx(2:nx-1,:) = rdx * Half * mfm(2:nx-1,:) *  &
               ((u(3:nx,:,k,t) - u(1:nx-2,:,k,t)) -  &
               u(2:nx-1,:,k,t) * (mfm(3:nx,:) - mfm(1:nx-2,:)))
          dudx((/1,nx/),:) = rdx * mfm((/1,nx/),:) *  &
               ((u((/2,nx/),:,k,t) - u((/1,nx-1/),:,k,t)) -  &
               u((/1,nx/),:,k,t) * (mfm((/2,nx/),:) - mfm((/1,nx-1/),:))) 
          dvdy(:,2:ny-1) = rdx * Half * mfm(:,2:ny-1) *  &
               ((v(:,3:ny,k,t) - v(:,1:ny-2,k,t)) -  &
               v(:,2:ny-1,k,t) * (mfm(:,3:ny) - mfm(:,1:ny-2)))
          dvdy(:,(/1,ny/)) = rdx * mfm(:,(/1,ny/)) *  &
               ((v(:,(/2,ny/),k,t) - v(:,(/1,ny-1/),k,t)) -  &
               v(:,(/1,ny/),k,t) * (mfm(:,(/2,ny/)) - mfm(:,(/1,ny-1/))))
          div(:,:,k,t) = dudx + dvdy
       end do
    end do


!!$    do t = 1, nt
!!$       do k = 1, nz
!!$          dudx(2:nx-1,:,k,t) = rdx * Half *  &
!!$               (mfm(3:nx,:) * u(3:nx,:,k,t) - mfm(1:nx-2,:) * u(1:nx-2,:,k,t))
!!$          dudx((/1,nx/),:,k,t) = rdx * (mfm((/2,nx/),:) * u((/2,nx/),:,k,t) - &
!!$               mfm((/1,nx-1/),:) * u((/1,nx-1/),:,k,t))
!!$          dvdy(:,2:ny-1,k,t) = rdx * Half *  &
!!$               (mfm(:,3:ny) * v(:,3:ny,k,t) - mfm(:,1:ny-2) * v(:,1:ny-2,k,t))
!!$          dvdy(:,(/1,ny/),k,t) = rdx * (mfm(:,(/2,ny/)) * v(:,(/2,ny/),k,t) - &
!!$               mfm(:,(/1,ny-1/)) * v(:,(/1,ny-1/),k,t))
!!$       end do
!!$    end do
    
!    div = spread(spread(mfm**2,3,nz),4,nt) * (dudx + dvdy)

!          div(2:nx-1,2:ny-1,k,t) = mfm**2 * rdx * Half * (  &
!               u(3:nx,2:ny-1,k,t) - u(1:nx-2,2:ny-1,k,t) +  &
!               v(2:nx-1,3:ny,k,t) - v(2:nx-1,1:ny-1,k,t))
!          div(1,
!          div(:,:,k,t) = mfm**2 * (cen_diff_stag(us(:,:,k,t)/mfu, rdx, 1) + cen_diff_stag(vs(:,:,k,t)/mfv, rdx, 2))
!          div(:,:,k,t) = mfm**2 * (cen_diff_stag(us(:,:,k,t), mfu, rdx, 1) + rdx * (mfv(:,2:ny+1)*vs(:,2:ny+1,k,t) - mfv(:,1:ny)*vs(:,1:ny,k,t)))
          
!          div(:,:,k,t) = mfm**2 * rdx * (mfu(2:nx+1,:)*us(2:nx+1,:,k,t)  &
!               + mfv(:,2:ny+1)*vs(:,2:ny+1,k,t) - mfu(1:nx,:)*us(1:nx,:,k,t)  &
!               - mfv(:,1:ny)*vs(:,1:ny,k,t))
  end function diag_Diverg

  ! ===========================================================================

  function diag_DScript(ncid, nx, ny, nz, nt) result(dScr)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: dScr
    
    real, dimension(nx+1,ny,nz,nt) :: mScr1
    real, dimension(nx,ny+1,nz,nt) :: mScr2 
    real, dimension(nx+1,ny) :: mfu
    real, dimension(nx,ny+1) :: mfv
    real, dimension(nx,ny)   :: mfm
    real    :: rdx
    integer :: t, k
 
    mScr1 = diag_MScript1(ncid, nx, ny, nz, nt)
    mScr2 = diag_MScript2(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)
    call read_values_const(ncid, "RDX", rdx)
    
    do t = 1, nt
       do k = 1, nz
          dScr(:,:,k,t) = mfm**2 * (  &
               cen_diff_stag(mScr1(:,:,k,t)/mfu, rdx, 1) +  &
               cen_diff_stag(mScr2(:,:,k,t)/mfv, rdx, 2))
       end do
    end do

!    print *, maxval(abs(dScr(:,:,:,13))), maxloc(abs(dScr(:,:,:,13)))
!    print *, mean(abs(dScr(:,:,9,13)))
  end function diag_DScript

!!$!    ms1((/1,nx+1/),:,:,:) = us((/1,nx+1/),:,:,:) * spread(mCap((/1,nx/),:,:),dim=3)
!!$!    ms1(2:nx,:,:,:) = .5 * us(2:nx,:,:,:) * spread(mCap(1:nx-1,:,:)+mCap(2:nx,:,:), dim=3)
!!$
!!$    do t = 1, nt
!!$       do k = 1, nz
!!$          do j = 1, ny
!!$             ms1(1,j,k,t) = us(1,j,k,t) * mCap(1,j,t)
!!$             do i = 2, nx
!!$                ms1(i,j,k,t) = us(i,j,k,t) * .5 * (mCap(i-1,j,t) + mCap(i,j,t))
!!$             end do
!!$             ms1(nx+1,j,k,t) = us(nx+1,j,k,t) * mCap(nx,j,t)
!!$          end do
!!$       end do
!!$    end do
    
  ! ===========================================================================

  function diag_DDDt(ncid, nx, ny, nz, nt) result(dd)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: dd
   
    real, parameter :: dtr = One / 3600.

    real, dimension(nx,ny,nz,nt) :: dScr

    dScr = diag_DScript(ncid, nx, ny, nz, nt)
    
    dd = cen_diff(dScr, dtr, 4)
    print *, maxval(abs(dd(:,:,29,5))), maxloc(abs(dd(:,:,29,5)))
    print *, mean(abs(dd(:,:,29,5)))
  end function diag_DDDt

  ! ===========================================================================

  function diag_EtaDot(ncid, nx, ny, nz, nt) result(ed)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: ed
    
    real, dimension(nx,ny,nz+1,nt) :: z
    real, dimension(nz+1)          :: eta

    z = diag_GeopHghtS(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, eta)

    ed = diag_WUnStag(ncid, nx, ny, nz, nt) *  &
         spread(spread(spread(eta(2:)-eta(:nz),1,nx),2,ny),4,nt) /  &
         (z(:,:,2:,:) - z(:,:,:nz,:))
!!$    print *, maxval(abs(ed(:,:,:,13))), maxloc(abs(ed(:,:,:,13)))
!!$    print *, mean(abs(ed(:,:,3,13)))
  end function diag_EtaDot

  ! ===========================================================================
  
  function diag_TermA2(ncid, nx, ny, nz, nt) result(a2)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: a2

    real, dimension(nz)          :: eta

    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, eta)
    
    a2 = diag_EtaDot(ncid, nx, ny, nz, nt) *  &
         uneven_deriv(diag_DScript(ncid, nx, ny, nz, nt), eta, 3)

    print *, maxval(abs(a2(:,:,29,5))), maxloc(abs(a2(:,:,29,5)))
    print *, mean(abs(a2(:,:,29,5)))
  end function diag_TermA2

  ! ===========================================================================

  function diag_TermA3(ncid, nx, ny, nz, nt) result(a3)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: a3

    real, dimension(nx+1,ny,nz,nt) :: ms1
    real, dimension(nx,ny+1,nz,nt) :: ms2
    real, dimension(nx,ny,nz,nt)   :: eDot
    real, dimension(nx,ny)         :: mfm
    real, dimension(nz)            :: eta
    real                           :: rdx

    ms1 = diag_MScript1(ncid, nx, ny, nz, nt)
    ms2 = diag_MScript2(ncid, nx, ny, nz, nt)
    eDot = diag_EtaDot(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, eta)
    call read_values_const(ncid, "RDX", rdx)
    
    a3 = spread(spread(mfm,3,nz),4,nt) * Half * (  &
         cen_diff(eDot, rdx, 1) *  &
         uneven_deriv(ms1(:nx,:,:,:) + ms1(2:,:,:,:), eta, 3) +  &
         cen_diff(eDot, rdx, 2) *  &
         uneven_deriv(ms2(:,:ny,:,:) + ms2(:,2:,:,:), eta, 3) )

    print *, maxval(abs(a3(:,:,29,5))), maxloc(abs(a3(:,:,29,5)))
    print *, mean(abs(a3(:,:,29,5)))
  end function diag_TermA3

  ! ===========================================================================
  
  function diag_TermA4(ncid, nx, ny, nz, nt) result(a4)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: a4

    real, parameter :: dtr = One / 3600.
    
    a4 = -diag_DScript(ncid, nx, ny, nz, nt) *  &
         cen_diff(spread(log(diag_MCap(ncid, nx, ny, nt)),3,nz), dtr, 4)
    
    print *, maxval(abs(a4(:,:,29,5))), maxloc(abs(a4(:,:,29,5)))
    print *, mean(abs(a4(:,:,29,5)))
  end function diag_TermA4

  ! ===========================================================================

  function diag_TermA5(ncid, nx, ny, nz, nt) result(a5)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: a5

    real, parameter :: dtr = One / 3600.
  
    real, dimension(nx+1,ny,nz,nt) :: ms1
    real, dimension(nx,ny+1,nz,nt) :: ms2
    real, dimension(nx,ny,nz,nt)   :: dlnmdt
    real, dimension(nx,ny)         :: mfm
    real                           :: rdx

    ms1 = diag_MScript1(ncid, nx, ny, nz, nt)
    ms2 = diag_MScript2(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, "RDX", rdx)

    dlnmdt = cen_diff(spread(log(diag_MCap(ncid, nx, ny, nt)),3,nz), dtr, 4)
    a5 = spread(spread(mfm,3,nz),4,nt) * Half *  &
         (-(ms1(:nx,:,:,:) + ms1(2:,:,:,:)) * cen_diff(dlnmdt, rdx, 1) -  &
         (ms2(:,:ny,:,:) + ms2(:,2:,:,:)) * cen_diff(dlnmdt, rdx, 2))

    print *, maxval(abs(a5(:,:,29,5))), maxloc(abs(a5(:,:,29,5)))
    print *, mean(abs(a5(:,:,29,5)))
  end function diag_TermA5

  ! ===========================================================================

  function diag_TermB1(ncid, nx, ny, nz, nt) result(b1)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: b1
    
    real, dimension(nx+1,ny,nz,nt) :: us, ms1
    real, dimension(nx,ny+1,nz,nt) :: vs, ms2
    real, dimension(nx,ny,nz,nt)   :: dvdx, dudy, dms2dx, dms1dy
    real, dimension(nx,ny)         :: mfm
    real    :: rdx
    integer :: k, t

    call read_values(ncid, wrfVar(UStag)%wrfName, us)
    ms1 = diag_MScript1(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(VStag)%wrfName, vs)
    ms2 = diag_MScript2(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, "RDX", rdx)

    dvdx = Half * cen_diff(vs(:,:ny,:,:) + vs(:,2:,:,:), rdx, 1)
    dudy = Half * cen_diff(us(:nx,:,:,:) + us(2:,:,:,:), rdx, 2)
    dms2dx = Half * cen_diff(ms2(:,:ny,:,:) + ms2(:,2:,:,:), rdx, 1)
    dms1dy = Half * cen_diff(ms1(:nx,:,:,:) + ms1(2:,:,:,:), rdx, 2)
    do t = 1, nt
       do k = 1, nz
          b1(:,:,k,t) = mfm**2 * (cen_diff_stag(us(:,:,k,t), rdx, 1) *  &
               cen_diff_stag(ms1(:,:,k,t), rdx, 1) + dvdx(:,:,k,t) *  &
               dms1dy(:,:,k,t) + dudy(:,:,k,t) * dms2dx(:,:,k,t) +  &
               cen_diff_stag(vs(:,:,k,t), rdx, 2) *  &
               cen_diff_stag(ms2(:,:,k,t), rdx, 2))
       end do
    end do

    print *, maxval(abs(b1(:,:,29,5))), maxloc(abs(b1(:,:,29,5)))
    print *, mean(abs(b1(:,:,29,5)))
  end function diag_TermB1

  ! ===========================================================================

  function diag_TermB2(ncid, nx, ny, nz, nt) result(b2)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: b2
  
    real, dimension(nx+1,ny,nz,nt) :: ms1
    real, dimension(nx,ny+1,nz,nt) :: ms2
    real, dimension(nx,ny,nz,nt)   :: u, v, lnm, vdglm, mfm4d
    real, dimension(nx,ny)         :: mfm
    real :: rdx

    ms1 = diag_MScript1(ncid, nx, ny, nz, nt)
    ms2 = diag_MScript2(ncid, nx, ny, nz, nt)
    u = diag_UUnstag(ncid, nx, ny, nz, nt)
    v = diag_VUnstag(ncid, nx, ny, nz, nt)
    lnm = spread(log(diag_MCap(ncid, nx, ny, nt)),3,nz)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, "RDX", rdx)
   
    mfm4d = spread(spread(mfm,3,nz),4,nt)
    vdglm = mfm4d * (u * cen_diff(lnm, rdx, 1) + v * cen_diff(lnm, rdx, 2))

    b2 = -Half * mfm4d * (  &
         (ms1(:nx,:,:,:) + ms1(2:,:,:,:)) * cen_diff(vdglm, rdx, 1) +  &
         (ms2(:,:ny,:,:) + ms2(:,2:,:,:)) * cen_diff(vdglm, rdx, 2) )

    print *, maxval(abs(b2(:,:,29,5))), maxloc(abs(b2(:,:,29,5)))
    print *, mean(abs(b2(:,:,29,5)))
  end function diag_TermB2

  ! ===========================================================================

  function diag_TermB3(ncid, nx, ny, nz, nt) result(b3)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: b3

    real, dimension(nx,ny,nz,nt) :: u, v, dScr, lnm
    real, dimension(nx,ny)       :: mfm
    real :: rdx
   
    u = diag_UUnstag(ncid, nx, ny, nz, nt)
    v = diag_VUnstag(ncid, nx, ny, nz, nt)
    dScr = diag_DScript(ncid, nx, ny, nz, nt)
    lnm = spread(log(diag_MCap(ncid, nx, ny, nt)),3,nz)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, "RDX", rdx)

    b3 = spread(spread(mfm,3,nz),4,nt) *  &
         (u * cen_diff(dScr, rdx, 1) + v * cen_diff(dScr, rdx, 2) -  &
         dScr * (u * cen_diff(lnm, rdx, 1) + v * cen_diff(lnm, rdx, 2)))
    print *, maxval(abs(b3(:,:,29,5))), maxloc(abs(b3(:,:,29,5)))
    print *, mean(abs(b3(:,:,29,5)))
  end function diag_TermB3

  ! ===========================================================================

  function diag_TermC(ncid, nx, ny, nz, nt) result(c)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: c

    real, dimension(nx+1,ny,nz,nt) :: ms1
    real, dimension(nx,ny+1,nz,nt) :: ms2    
    real, dimension(nx,ny,nz,nt)   :: mfm4d, ms1u, ms2u
    real, dimension(nx,ny)         :: mfm, f
    real :: rdx
    
    ms1 = diag_MScript1(ncid, nx, ny, nz, nt)
    ms2 = diag_MScript2(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, "RDX", rdx)    

    mfm4d = spread(spread(mfm,3,nz),4,nt)
    ms1u = Half * (ms1(:nx,:,:,:) + ms1(2:,:,:,:))
    ms2u = Half * (ms2(:,:ny,:,:) + ms2(:,2:,:,:))
    c = -spread(spread(f,3,nz),4,nt) *  &
         (mfm4d * cen_diff(ms2u, rdx, 1) - mfm4d * cen_diff(ms1u, rdx, 2) -  &
         ms2u * cen_diff(mfm4d, rdx, 1) + ms1u * cen_diff(mfm4d, rdx, 2)) 

    print *, maxval(abs(c(:,:,29,5))), maxloc(abs(c(:,:,29,5)))
    print *, mean(abs(c(:,:,29,5)))  
  end function diag_TermC

  ! ===========================================================================

  function diag_TermD(ncid, nx, ny, nz, nt) result(d)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: d
    
    real, dimension(nx,ny,nz,nt)   :: f4d
    real, dimension(nx,ny)         :: mfm, f
    real :: rdx
    
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, "RDX", rdx)    
    
    f4d = spread(spread(f,3,nz),4,nt)
    d = spread(spread(mfm,3,nz),4,nt) * (  &
         diag_MS1Unstag(ncid, nx, ny, nz, nt) * cen_diff(f4d, rdx, 2) -  &
         diag_MS2Unstag(ncid, nx, ny, nz, nt) * cen_diff(f4d, rdx, 1) )
   
    print *, maxval(abs(d(:,:,29,5))), maxloc(abs(d(:,:,29,5)))
    print *, mean(abs(d(:,:,29,5)))
  end function diag_TermD

  ! ===========================================================================

  function diag_WaterInv(ncid, nx, ny, nz, nt) result(chi)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: chi
    
    real, dimension(nx,ny,nz,nt) :: qv, qc, qr

    call read_values(ncid, wrfVar(MixR)%wrfName, qv)
    call read_values(ncid, wrfVar(MixRCloud)%wrfName, qc)
    call read_values(ncid, wrfVar(MixRRain)%wrfName, qr)

    chi = One / (qr + qc + qv + One)
  end function diag_WaterInv

  ! ===========================================================================

  function diag_VapPres(ncid, nx, ny, nz, nt) result(e)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: e
   
    real, dimension(nx,ny,nz,nt) :: p, qv

    p = diag_Pressure(ncid, nx, ny, nz, nt, si=.true.)
    call read_values(ncid, wrfVar(MixR)%wrfName, qv)

    e = p / (One + Eps/qv)
  end function diag_VapPres

  ! ===========================================================================

  function diag_GeopS(ncid, nx, ny, nz, nt) result(phi)
    integer, intent(in)            :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz+1,nt) :: phi
    
    real, dimension(nx,ny,nz+1,nt) :: gpb, gpp
    
    call read_values(ncid, wrfVar(GeopBase)%wrfName, gpb)
    call read_values(ncid, wrfVar(GeopPert)%wrfName, gpp)
    phi = gpb + gpp
  end function diag_GeopS

  ! ===========================================================================

  function diag_TermE(ncid, nx, ny, nz, nt) result(e)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: e
    
    real, dimension(nx,ny,nz+1,nt) :: phis
    real, dimension(nx,ny,nz,nt)   :: chi, pv, phi, eta4d, mur4d, m4d,  &
                                      dphideta, mfm4d
    real, dimension(nx,ny,nt)      :: mu, m
    real, dimension(nx,ny)         :: mfm
    real, dimension(nz+1)          :: etaz
    real, dimension(nz)            :: etam
    real    :: rdx, pt
    integer :: i, j, t

    print *, 1
    phis = diag_GeopS(ncid, nx, ny, nz, nt)
    print *, 2
    chi = One  ! For ideal case
!    chi = diag_WaterInv(ncid, nx, ny, nz, nt)
    print *, 3
    pv = diag_VapPres(ncid, nx, ny, nz, nt)
    print *, 4
    mu = diag_DryAir(ncid, nx, ny, nt)
    print *, 5
    m = diag_MCap(ncid, nx, ny, nt)
    print *, 6
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "RDX", rdx)
    call read_values_const(ncid, "P_TOP", pt)

    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             phi(i,j,:,t) =  &
                  unstag_hght(phis(i,j,:,t), etaz, etam, mu(i,j,t), pt)
          end do
       end do
    end do
    eta4d = spread(spread(spread(etam,1,nx),2,ny),4,nt)
    mur4d = spread(One/mu, 3, nz)
    m4d = spread(m,3,nz)
    dphideta = (phis(:,:,2:,:) - phis(:,:,:nz,:)) /  &
         spread(spread(spread(etaz(2:) - etaz(:nz),1,nx),2,ny),4,nt)
    mfm4d = spread(spread(mfm,3,nz),4,nt)
    e = -mfm4d * (  &
         cen_diff(chi * mfm4d * ((One + mur4d * uneven_deriv(pv,etam,3)) *  &
         m4d * cen_diff(phi,rdx,1) - (eta4d * cen_diff(m4d,rdx,1) + m4d *  &
         mur4d * cen_diff(pv,rdx,1)) * dphideta), rdx, 1) +  &
         cen_diff(chi * mfm4d * ((One + mur4d * uneven_deriv(pv,etam,3)) *  &
         m4d * cen_diff(phi,rdx,2) - (eta4d * cen_diff(m4d,rdx,2) + m4d *  &
         mur4d * cen_diff(pv,rdx,2)) * dphideta), rdx, 2) )

    print *, maxval(abs(e(:,:,29,5))), maxloc(abs(e(:,:,29,5)))
    print *, mean(abs(e(:,:,29,5)))
  end function diag_TermE

  ! ===========================================================================
  
  function diag_B1B2CD(ncid, nx, ny, nz, nt) result(total)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: total
    
    total = diag_TermB1(ncid, nx, ny, nz, nt) +  &
         diag_TermB2(ncid, nx, ny, nz, nt) +  &
         diag_TermC(ncid, nx, ny, nz, nt) +  &
         diag_TermD(ncid, nx, ny, nz, nt)
  end function diag_B1B2CD

  ! ===========================================================================

  function diag_ABCD(ncid, nx, ny, nz, nt) result(total)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: total
  
    total = diag_DDDt(ncid,nx,ny,nz,nt) + diag_TermA2(ncid,nx,ny,nz,nt) +  &
         diag_TermA3(ncid,nx,ny,nz,nt) + diag_TermA4(ncid,nx,ny,nz,nt) +  &
         diag_TermA5(ncid,nx,ny,nz,nt) + diag_TermB1(ncid,nx,ny,nz,nt) +  &
         diag_TermB2(ncid,nx,ny,nz,nt) + diag_TermB3(ncid,nx,ny,nz,nt) +  &
         diag_TermC(ncid,nx,ny,nz,nt) + diag_TermD(ncid,nx,ny,nz,nt)
  end function diag_ABCD

  ! ===========================================================================

  function diag_MVort(ncid, nx, ny, nz, nt) result(vort)
    integer, intent(in)              :: ncid, nx, ny, nz, nt
    real, dimension(nx-1,ny-1,nz,nt) :: vort

    real, dimension(0:nx,ny,nz,nt)   :: m1s
    real, dimension(nx,0:ny,nz,nt)   :: m2s
    real, dimension(0:nx,ny)         :: mfu
    real, dimension(nx,0:ny)         :: mfv    
    real, dimension(nx-1,ny-1)       :: mfpsi
    real                             :: rdx

    integer :: i, j, k, t

    m1s = diag_MScript1(ncid, nx, ny, nz, nt)
    m2s = diag_MScript2(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    mfpsi = calc_mfpsi(mfu(1:nx-1,:), mfv(:,1:ny-1))
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
       do k = 1, nz
          vort(:,:,k,t) = mfpsi**2 * (  &
               cen_diff_stag(m2s(:,1:ny-1,k,t)/mfv(:,1:ny-1), rdx, 1) -  &
               cen_diff_stag(m1s(1:nx-1,:,k,t)/mfu(1:nx-1,:), rdx, 2) )
       end do
    end do
!!$    print *
!!$!    print *, m2s(1:3,0:2,10,13)
!!$!    print *, mfv(1:3,0:2)
!!$    print *
!!$    print *, m2s(1:3,1:1,10,13)/mfv(1:3,1:1), rdx
!!$    print *
!!$    print *, cen_diff_stag(m2s(:,1:1,10,13)/mfv(:,1:1), rdx, 1)
!!$    print *, vort(1,1,10,13)
  end function diag_MVort

  ! ===========================================================================

  function diag_Stream(ncid, nx, ny, nz, nt) result(psi)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: psi

!    integer, parameter :: k = 18, t = 40

    real, dimension(0:nx,ny,nz,nt)   :: us  ! Really the x-component of MScr
    real, dimension(nx,0:ny,nz,nt)   :: vs  ! Really the y-component of MScr
    real, dimension(nx,ny,0:nz,nt)   :: ph
    real, dimension(nx,ny,nz,nt)     :: div
    real, dimension(nx-1,ny-1,nz,nt) :: vort
    real(cp), dimension(0:nx+1,0:ny+1) :: ch
    real(cp), dimension(0:nx,0:ny)   :: ps
    real, dimension(0:nx,ny)         :: mfu
    real, dimension(nx,0:ny)         :: mfv
    real, dimension(nx,ny)           :: mfm
    real, dimension(nx-1,ny-1)       :: mfpsi
    real(cp) :: phill
    real    :: rdx
    integer :: i, j, k, t

    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    us = diag_MScript1(ncid, nx, ny, nz, nt)
    vs = diag_MScript2(ncid, nx, ny, nz, nt)
    div = diag_DScript(ncid, nx, ny, nz, nt)
    vort = diag_MVort(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)
    mfpsi = calc_mfpsi(mfu(1:nx-1,:), mfv(:,1:ny-1))
    call read_values_const(ncid, "RDX", rdx)    

    open(33, file="ps.dat", status="new", form="unformatted", action="write")
    write (33) nx, ny, nz, nt
!    t = 13
    do t = 1, nt
       do k = 1, nz
          write (*,"(2(a,i0))") "Time: ", t, "  Level: ", k
          phill = Quarter * (ph(1,1,k-1,t) + ph(1,2,k-1,t) + ph(1,1,k,t) + ph(1,2,k,t))
          call partition(real(us(:,:,k,t),cp), real(vs(:,:,k,t),cp),  &
               real(div(:,:,k,t),cp), real(vort(:,:,k,t),cp), real(mfm,cp),  &
               real(mfpsi,cp), real(mfu,cp), real(mfv,cp), phill, real(rdx,cp), ch,ps)
          write (33) ps
          psi(:,:,k,t) = .25 *  &
!          psi(:,:,k,:) = spread(.25 *  &
               (ps(:nx-1,:ny-1) + ps(:nx-1,1:) + ps(1:,:ny-1) + ps(1:,1:))!, 3, nt)
       end do
    end do
    close(33)
!    do j = 1, ny
!       do i = 1, nx
!          psi(i,j,:,:) = .25 * (ps(i-1,j-1) + ps(i-1,j) + ps(i,j-1) + ps(i,j))
!       end do
!    end do
    ! Correction for 0 in corners
    psi(1,1,:,:) = psi(1,1,:,:) * 1.3333333
    psi(nx,1,:,:) = psi(nx,1,:,:) * 1.3333333
!    psi(1,ny,:,:) = psi(1,ny,:,:) * 1.3333333
!    psi(nx,ny,:,:) = psi(nx,ny,:,:) * 1.3333333
!    print *, maxval(abs(psi(:,:,10,13)))
!    print *, mean(abs(psi(:,:,10,13)))
  end function diag_Stream

  ! ===========================================================================
  
  function diag_VelocPot(ncid, nx, ny, nz, nt) result(chi)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: chi
    
    integer, parameter :: k = 18, t = 40
    
    real, dimension(0:nx,ny,nz,nt)   :: us  ! Really the x-component of MScr
    real, dimension(nx,0:ny,nz,nt)   :: vs  ! Really the y-component of MScr
    real, dimension(nx,ny,0:nz,nt)   :: ph
    real, dimension(nx,ny,nz,nt)     :: div
    real, dimension(nx-1,ny-1,nz,nt) :: vort
    real(cp), dimension(0:nx+1,0:ny+1) :: ch
    real(cp), dimension(0:nx,0:ny)     :: ps
    real, dimension(0:nx,ny)         :: mfu
    real, dimension(nx,0:ny)         :: mfv
    real, dimension(nx,ny)           :: mfm
    real, dimension(nx-1,ny-1)       :: mfpsi
    real(cp) :: phill
    real    :: rdx
    integer :: i, j
    
    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    us = diag_MScript1(ncid, nx, ny, nz, nt)
    vs = diag_MScript2(ncid, nx, ny, nz, nt)
    div = diag_DScript(ncid, nx, ny, nz, nt)
    vort = diag_MVort(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)
    mfpsi = calc_mfpsi(mfu(1:nx-1,:), mfv(:,1:ny-1))
    call read_values_const(ncid, "RDX", rdx)    
    
    phill = Quarter * (ph(1,1,k-1,t) + ph(1,2,k-1,t) + ph(1,1,k,t) + ph(1,2,k,t))
    call partition(real(us(:,:,k,t),cp), real(vs(:,:,k,t),cp),  &
         real(div(:,:,k,t),cp), real(vort(:,:,k,t),cp), real(mfm,cp),  &
         real(mfpsi,cp), real(mfu,cp), real(mfv,cp), phill, real(rdx,cp), ch, ps)

    do j = ny+1, ny-1, -1
       print *, ch(nx-1:nx+1,j)
    end do
    chi = spread(spread(ch(1:nx,1:ny),3,nz),4,nt)
    print *, maxval(abs(chi(:,:,10,13)))
    print *, mean(abs(chi(:,:,10,13)))
  end function diag_VelocPot

  ! ===========================================================================

  function diag_MS1Psi(ncid, nx, ny, nz, nt) result(ms1)
    integer, intent(in)      :: ncid, nx, ny, nz, nt
    real, dimension(0:nx,ny) :: ms1
    
    integer, parameter :: k = 18, t = 40
    
    real, dimension(0:nx,ny,nz,nt)   :: us  ! Really the x-component of MScr
    real, dimension(nx,0:ny,nz,nt)   :: vs  ! Really the y-component of MScr
    real, dimension(nx,ny,0:nz,nt)   :: ph
    real, dimension(nx,ny,nz,nt)     :: div
    real, dimension(nx-1,ny-1,nz,nt) :: vort
    real(cp), dimension(0:nx+1,0:ny+1) :: ch
    real(cp), dimension(0:nx,0:ny)     :: ps
    real, dimension(0:nx,ny)         :: mfu
    real, dimension(nx,0:ny)         :: mfv
    real, dimension(nx,ny)           :: mfm
    real, dimension(nx-1,ny-1)       :: mfpsi
    real(cp) :: phill
    real    :: rdx
    integer :: i, j
    
    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    us = diag_MScript1(ncid, nx, ny, nz, nt)
    vs = diag_MScript2(ncid, nx, ny, nz, nt)
    div = diag_DScript(ncid, nx, ny, nz, nt)
    vort = diag_MVort(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)
    mfpsi = calc_mfpsi(mfu(1:nx-1,:), mfv(:,1:ny-1))
    call read_values_const(ncid, "RDX", rdx)    
    
    phill = Quarter * (ph(1,1,k-1,t) + ph(1,2,k-1,t) + ph(1,1,k,t) + ph(1,2,k,t))
    call partition(real(us(:,:,k,t),cp), real(vs(:,:,k,t),cp),  &
         real(div(:,:,k,t),cp), real(vort(:,:,k,t),cp), real(mfm,cp),  &
         real(mfpsi,cp), real(mfu,cp), real(mfv,cp), phill, real(rdx,cp), ch, ps)
    ms1 = -mfu * rdx * (ps(:,1:) - ps(:,:ny-1))
  end function diag_MS1Psi

  ! ===========================================================================

  function diag_MS2Psi(ncid, nx, ny, nz, nt) result(ms2)
    integer, intent(in)      :: ncid, nx, ny, nz, nt
    real, dimension(nx,0:ny) :: ms2
    
    integer, parameter :: k = 18, t = 40
    
    real, dimension(0:nx,ny,nz,nt)   :: us  ! Really the x-component of MScr
    real, dimension(nx,0:ny,nz,nt)   :: vs  ! Really the y-component of MScr
    real, dimension(nx,ny,0:nz,nt)   :: ph
    real, dimension(nx,ny,nz,nt)     :: div
    real, dimension(nx-1,ny-1,nz,nt) :: vort
    real(cp), dimension(0:nx+1,0:ny+1) :: ch
    real(cp), dimension(0:nx,0:ny)     :: ps
    real, dimension(0:nx,ny)         :: mfu
    real, dimension(nx,0:ny)         :: mfv
    real, dimension(nx,ny)           :: mfm
    real, dimension(nx-1,ny-1)       :: mfpsi
    real(cp) :: phill
    real    :: rdx
    integer :: i, j
    
    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    us = diag_MScript1(ncid, nx, ny, nz, nt)
    vs = diag_MScript2(ncid, nx, ny, nz, nt)
    div = diag_DScript(ncid, nx, ny, nz, nt)
    vort = diag_MVort(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)
    mfpsi = calc_mfpsi(mfu(1:nx-1,:), mfv(:,1:ny-1))
    call read_values_const(ncid, "RDX", rdx)    
    
    phill = Quarter * (ph(1,1,k-1,t) + ph(1,2,k-1,t) + ph(1,1,k,t) + ph(1,2,k,t))
    call partition(real(us(:,:,k,t),cp), real(vs(:,:,k,t),cp),  &
         real(div(:,:,k,t),cp), real(vort(:,:,k,t),cp), real(mfm,cp),  &
         real(mfpsi,cp), real(mfu,cp), real(mfv,cp), phill, real(rdx,cp), ch, ps)
    do j = 0, ny
       do i = 1, nx
          ms2(i,j) = mfv(i,j) * rdx * (ps(i,j) - ps(i-1,j))
       end do
    end do
  end function diag_MS2Psi

  ! ===========================================================================

  function diag_MS1PsiUnstag(ncid, nx, ny, nz, nt) result (ms1u)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: ms1u
   
    real, dimension(nx+1,ny,nz,nt) :: ms1
    
    ms1 = spread(spread(diag_MS1Psi(ncid, nx, ny, nz, nt),3,nz),4,nt)

    ms1u = Half * (ms1(:nx,:,:,:) + ms1(2:,:,:,:))
  end function diag_MS1PsiUnstag

  ! ===========================================================================

  function diag_MS2PsiUnstag(ncid, nx, ny, nz, nt) result (ms2u)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: ms2u
    
    real, dimension(nx,ny+1,nz,nt) :: ms2
    
    ms2 = spread(spread(diag_MS2Psi(ncid, nx, ny, nz, nt),3,nz),4,nt)
    
    ms2u = Half * (ms2(:,:ny,:,:) + ms2(:,2:,:,:))
  end function diag_MS2PsiUnstag

  ! ===========================================================================
  
  function diag_DivPsi(ncid, nx, ny, nz, nt) result(div)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: div
    
    integer, parameter :: k = 18, t = 40
    
    real, dimension(0:nx,ny) :: mfu, ms1
    real, dimension(nx,0:ny) :: mfv, ms2
    real, dimension(nx,ny)   :: mfm
    real    :: rdx
    integer :: i, j
    
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    ms1 = diag_MS1Psi(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    ms2 = diag_MS2Psi(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)
    call read_values_const(ncid, "RDX", rdx)    
    
!!$    div = spread(spread(mfm**2 * rdx * (  &
!!$               ms1(:nx,:)/mfu(:nx,:) - ms1(0:,:)/mfu(0:,:) +  &
!!$               ms2(:,:ny)/mfv(:,:ny) - ms2(:,0:)/mfv(:,0:) ),3,nz),4,nt)
    do j = 1, ny
       do i = 1, nx
          div(i,j,:,:) = mfm(i,j)**2 * rdx * (  &
               ms1(i,j)/mfu(i,j) - ms1(i-1,j)/mfu(i-1,j) +  &
               ms2(i,j)/mfv(i,j) - ms2(i,j-1)/mfv(i,j-1) )
       end do
    end do
    print *, maxval(abs(div(:,:,k,t)))
    print *, mean(abs(div(:,:,k,t)))
  end function diag_DivPsi
  
  ! ===========================================================================

  function diag_MS1Chi(ncid, nx, ny, nz, nt) result(ms1)
    integer, intent(in)      :: ncid, nx, ny, nz, nt
    real, dimension(0:nx,ny) :: ms1
    
    integer, parameter :: k = 18, t = 40
    
    real, dimension(0:nx,ny,nz,nt)   :: us  ! Really the x-component of MScr
    real, dimension(nx,0:ny,nz,nt)   :: vs  ! Really the y-component of MScr
    real, dimension(nx,ny,0:nz,nt)   :: ph
    real, dimension(nx,ny,nz,nt)     :: div
    real, dimension(nx-1,ny-1,nz,nt) :: vort
    real(cp), dimension(0:nx+1,0:ny+1) :: ch
    real(cp), dimension(0:nx,0:ny)     :: ps
    real, dimension(0:nx,ny)         :: mfu
    real, dimension(nx,0:ny)         :: mfv
    real, dimension(nx,ny)           :: mfm
    real, dimension(nx-1,ny-1)       :: mfpsi
    real(cp) :: phill
    real    :: rdx
    integer :: i, j
    
    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    us = diag_MScript1(ncid, nx, ny, nz, nt)
    vs = diag_MScript2(ncid, nx, ny, nz, nt)
    div = diag_DScript(ncid, nx, ny, nz, nt)
    vort = diag_MVort(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)
    mfpsi = calc_mfpsi(mfu(1:nx-1,:), mfv(:,1:ny-1))
    call read_values_const(ncid, "RDX", rdx)    
    
    phill = Quarter * (ph(1,1,k-1,t) + ph(1,2,k-1,t) + ph(1,1,k,t) + ph(1,2,k,t))
    call partition(real(us(:,:,k,t),cp), real(vs(:,:,k,t),cp),  &
         real(div(:,:,k,t),cp), real(vort(:,:,k,t),cp), real(mfm,cp),  &
         real(mfpsi,cp), real(mfu,cp), real(mfv,cp), phill, real(rdx,cp), ch, ps)
    do j = 1, ny
       do i = 0, nx
          ms1(i,j) = mfu(i,j) * rdx * (ch(i+1,j) - ch(i,j))
       end do
    end do
  end function diag_MS1Chi

  ! ===========================================================================

  function diag_MS2Chi(ncid, nx, ny, nz, nt) result(ms2)
    integer, intent(in)      :: ncid, nx, ny, nz, nt
    real, dimension(nx,0:ny) :: ms2
    
    integer, parameter :: k = 18, t = 40
    
    real, dimension(0:nx,ny,nz,nt)   :: us  ! Really the x-component of MScr
    real, dimension(nx,0:ny,nz,nt)   :: vs  ! Really the y-component of MScr
    real, dimension(nx,ny,0:nz,nt)   :: ph
    real, dimension(nx,ny,nz,nt)     :: div
    real, dimension(nx-1,ny-1,nz,nt) :: vort
    real(cp), dimension(0:nx+1,0:ny+1) :: ch
    real(cp), dimension(0:nx,0:ny)     :: ps
    real, dimension(0:nx,ny)         :: mfu
    real, dimension(nx,0:ny)         :: mfv
    real, dimension(nx,ny)           :: mfm
    real, dimension(nx-1,ny-1)       :: mfpsi
    real(cp) :: phill
    real    :: rdx
    integer :: i, j
    
    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    us = diag_MScript1(ncid, nx, ny, nz, nt)
    vs = diag_MScript2(ncid, nx, ny, nz, nt)
    div = diag_DScript(ncid, nx, ny, nz, nt)
    vort = diag_MVort(ncid, nx, ny, nz, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)
    mfpsi = calc_mfpsi(mfu(1:nx-1,:), mfv(:,1:ny-1))
    call read_values_const(ncid, "RDX", rdx)    
    
    phill = Quarter * (ph(1,1,k-1,t) + ph(1,2,k-1,t) + ph(1,1,k,t) + ph(1,2,k,t))
    call partition(real(us(:,:,k,t),cp), real(vs(:,:,k,t),cp),  &
         real(div(:,:,k,t),cp), real(vort(:,:,k,t),cp), real(mfm,cp),  &
         real(mfpsi,cp), real(mfu,cp), real(mfv,cp), phill, real(rdx,cp), ch, ps)
    do j = 0, ny
       do i = 1, nx
          ms2(i,j) = mfv(i,j) * rdx * (ch(i,j+1) - ch(i,j))
       end do
    end do
  end function diag_MS2Chi

  ! ===========================================================================

  function diag_MS1ChiUnstag(ncid, nx, ny, nz, nt) result (ms1u)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: ms1u
   
    real, dimension(nx+1,ny,nz,nt) :: ms1
    
    ms1 = spread(spread(diag_MS1Chi(ncid, nx, ny, nz, nt),3,nz),4,nt)

    ms1u = Half * (ms1(:nx,:,:,:) + ms1(2:,:,:,:))
  end function diag_MS1ChiUnstag

  ! ===========================================================================

  function diag_MS2ChiUnstag(ncid, nx, ny, nz, nt) result (ms2u)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: ms2u
    
    real, dimension(nx,ny+1,nz,nt) :: ms2
    
    ms2 = spread(spread(diag_MS2Chi(ncid, nx, ny, nz, nt),3,nz),4,nt)
    
    ms2u = Half * (ms2(:,:ny,:,:) + ms2(:,2:,:,:))
  end function diag_MS2ChiUnstag

  ! ===========================================================================

  function diag_balLHS(ncid, nx, ny, nz, nt) result (lhs)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: lhs

    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: mfm, f
    real    :: rdx
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)
    mCap = diag_MCap(ncid, nx, ny, nt)   
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
       print *, t
       call balance_lhs(ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(mfu,cp), real(mfv,cp), real(rdx,cp),  &
            lhs(:,:,:,t))
    end do
  end function diag_balLHS

  ! ===========================================================================

  function diag_balLHSExp(ncid, nx, ny, nz, nt) result (lhs)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: lhs
    
    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: mfm, f
    real    :: rdx
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)
    mCap = diag_MCap(ncid, nx, ny, nt)   
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
       print *, t
       call balance_lhs_exp(ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(mfu,cp), real(mfv,cp), real(rdx,cp),  &
            lhs(:,:,:,t))
    end do
    print *, lhs(10,11,9,9)
  end function diag_balLHSExp

  ! ===========================================================================

  function diag_balLHSPert(ncid, nx, ny, nz, nt) result (lhs)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: lhs
    
    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: mfm, f
    real    :: rdx
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)
    ps(5,6,:,:) = ps(5,6,:,:) + 5e6
    mCap = diag_MCap(ncid, nx, ny, nt)   
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, "RDX", rdx)
    
    do t = 1, nt
!       print *, t
       call balance_lhs_exp(ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(mfu,cp), real(mfv,cp), real(rdx,cp),  &
            lhs(:,:,:,t))
    end do
    print *, lhs(10,11,9,9)
  end function diag_balLHSPert

  ! ===========================================================================

  function diag_balLHSTLM(ncid, nx, ny, nz, nt) result (lhs)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: lhs
    
    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps, dps
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: mfm, f
    real    :: rdx
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)

    dps = 0.
    dps(5,6,:,:) = 5e6

    mCap = diag_MCap(ncid, nx, ny, nt)   
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
!       print *, t
       call balance_lhs_tlm(dps(:,:,:,t), ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(mfu,cp), real(mfv,cp), real(rdx,cp),  &
            lhs(:,:,:,t))
    end do
    print *, lhs(10,11,9,9)
  end function diag_balLHSTLM

  ! ===========================================================================

  function diag_balLHSAdj(ncid, nx, ny, nz, nt) result (adj)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: adj
    
    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps, dps, ataq
    real, dimension(nx,ny,nz,nt) :: dlhs
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: mfm, f
    real    :: rdx
    integer :: i, j, k, t

    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)

    dps = 0.
    dps(5,6,:,:) = 5e6

    mCap = diag_MCap(ncid, nx, ny, nt)   
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
!       print *, t
       call balance_lhs_tlm(dps(:,:,:,t), ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(mfu,cp), real(mfv,cp), real(rdx,cp),  &
            dlhs(:,:,:,t))
    end do
    t = 13
    print *, real(sum(real(dlhs,cp)**2))
    do t = 1, nt
       call balance_lhs_adj(real(dlhs(:,:,:,t),cp), ps(:,:,:,t), real(f,cp),  &
            real(mCap(:,:,t),cp), real(mfm,cp), real(mfu,cp), real(mfv,cp),  &
            real(rdx,cp), ataq(:,:,:,t))
    end do
    print *, real(sum(dps*real(ataq,cp)))
    adj = 0
  end function diag_balLHSAdj

  ! ===========================================================================

  function diag_balRHSSimp(ncid, nx, ny, nz, nt) result (rhs)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: rhs
    
    real, dimension(nx,ny,0:nz,nt) :: ph
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(nx,ny)    :: mfm
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real    :: rdx
    integer :: t

    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    mCap = diag_MCap(ncid, nx, ny, nt)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
       print *, t
       call balance_rhs_simp(real(ph(:,:,:,t),cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(etaz,cp), real(etam,cp), real(rdx,cp),  &
            rhs(:,:,:,t))
    end do
  end function diag_balRHSSimp

  ! ===========================================================================

  function diag_balRHSSimpExp(ncid, nx, ny, nz, nt) result (rhs)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: rhs
    
    real, dimension(nx,ny,0:nz,nt) :: ph
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(nx,ny)    :: mfm
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real    :: rdx
    integer :: t
    
    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    mCap = diag_MCap(ncid, nx, ny, nt)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "RDX", rdx)
    
    do t = 1, nt
       print *, t
       call balance_rhs_simp_exp(real(ph(:,:,:,t),cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(etaz,cp), real(etam,cp), real(rdx,cp),  &
            rhs(:,:,:,t))
    end do
  end function diag_balRHSSimpExp

  ! ===========================================================================

  function diag_balRHSSimpPert(ncid, nx, ny, nz, nt) result (rhs)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: rhs
    
    real, dimension(nx,ny,0:nz,nt) :: ph
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(nx,ny)    :: mfm
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real    :: rdx
    integer :: t
    
    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    ph(5,6,:,:) = ph(5,6,:,:) + 50.
    mCap = diag_MCap(ncid, nx, ny, nt)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "RDX", rdx)
    
    do t = 1, nt
       print *, t
       call balance_rhs_simp_exp(real(ph(:,:,:,t),cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(etaz,cp), real(etam,cp), real(rdx,cp),  &
            rhs(:,:,:,t))
    end do
  end function diag_balRHSSimpPert

  ! ===========================================================================
  
  function diag_balRHSSimpTLM(ncid, nx, ny, nz, nt) result (rhs)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: rhs
    
    real(cp), dimension(nx,ny,0:nz,nt) :: dph
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(nx,ny)    :: mfm
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real    :: rdx
    integer :: i, j, k, t
    
    dph = 0.
    dph(2,6,:,:) = 50.

    mCap = diag_MCap(ncid, nx, ny, nt)   
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
!       print *, t
       call balance_rhs_simp_tlm(dph(:,:,:,t), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(etaz,cp), real(etam,cp), real(rdx,cp),  &
            rhs(:,:,:,t))
    end do
    print *, rhs(10,11,9,9)
  end function diag_balRHSSimpTLM

  ! ===========================================================================

  function diag_balRHSSimpAdj(ncid, nx, ny, nz, nt) result (adj)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: adj

    real(cp), dimension(nx,ny,0:nz,nt) :: dph, ataq
    real, dimension(nx,ny,nz,nt) :: drhs
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(nx,ny)    :: mfm
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real    :: rdx
    integer :: i, j, k, t

    dph = 0.
    dph(5,6,:,:) = 50.

    mCap = diag_MCap(ncid, nx, ny, nt)   
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "RDX", rdx)
    
    do t = 1, nt
       call balance_rhs_simp_tlm(dph(:,:,:,t), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(etaz,cp), real(etam,cp), real(rdx,cp),  &
            drhs(:,:,:,t))
    end do
    print *, drhs(10,11,9,9)
    print *, real(sum(real(drhs,cp)**2))
    do t = 1, nt
!       print *, t
       call balance_rhs_simp_adj(real(drhs(:,:,:,t),cp),  &
            real(mCap(:,:,t),cp), real(mfm,cp), real(etaz,cp), real(etam,cp),  &
            real(rdx,cp), ataq(:,:,:,t))
    end do
    print *, real(sum(dph*real(ataq,cp)))
    adj = 0
  end function diag_balRHSSimpAdj

  ! ===========================================================================

  function diag_ImbalSimp(ncid, nx, ny, nz, nt) result (imb)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: imb

    real, dimension(nx,ny,nz,nt) :: lhs, rhs
    real    :: norm
    integer :: t

    lhs = diag_balLHS(ncid, nx, ny, nz, nt)
    rhs = diag_balRHSSimp(ncid, nx, ny, nz, nt)
    imb = lhs - rhs
    
!!$    t = 40
!!$    norm = sum(imb(:,:,:,40)**2)
!!$    print *, norm
!!$    do t = 30, 41
!!$       print *, sum((real(lhs(:,:,:,40),cp)-real(rhs(:,:,:,t),cp))**2)
!!$    end do
  end function diag_ImbalSimp

  ! ===========================================================================

  function diag_ImbalSimpPert(ncid, nx, ny, nz, nt) result (imb)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: imb
    
    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps
    real, dimension(nx,ny,0:nz,nt) :: ph
    real, dimension(nx,ny,nz,nt) :: lhs, rhs
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: mfm, f
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real    :: rdx, norm
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)
    ps(5,6,:,:) = ps(5,6,:,:) + 5e6
    mCap = diag_MCap(ncid, nx, ny, nt)   
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, "RDX", rdx)
    
    do t = 1, nt
       call balance_lhs_exp(ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(mfu,cp), real(mfv,cp), real(rdx,cp),  &
            lhs(:,:,:,t))
    end do

    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    ph(5,6,:,:) = ph(5,6,:,:) + 50.
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    
    do t = 1, nt
       call balance_rhs_simp_exp(real(ph(:,:,:,t),cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(etaz,cp), real(etam,cp), real(rdx,cp),  &
            rhs(:,:,:,t))
    end do

    imb = lhs - rhs
    
!!$    t = 40
!!$    norm = sum(imb(:,:,:,40)**2)
!!$    print *, norm
!!$    do t = 30, 41
!!$       print *, sum((real(lhs(:,:,:,40),cp)-real(rhs(:,:,:,t),cp))**2)
!!$    end do
  end function diag_ImbalSimpPert
  
  ! ===========================================================================

  function diag_ImbalSimpTLM(ncid, nx, ny, nz, nt) result (dimb)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: dimb
    
    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps, dps
    real(cp), dimension(nx,ny,0:nz,nt) :: dph
    real, dimension(nx,ny,nz,nt) :: dlhs, drhs
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: mfm, f
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real    :: rdx, norm
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)

    dps = 0.
    dps(5,6,:,:) = 5e6

    mCap = diag_MCap(ncid, nx, ny, nt)   
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, "RDX", rdx)
    
    do t = 1, nt
       call balance_lhs_tlm(dps(:,:,:,t), ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(mfu,cp), real(mfv,cp), real(rdx,cp), dlhs(:,:,:,t))
    end do

    dph = 0.
    dph(5,6,:,:) = 50.

    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    
    do t = 1, nt
       call balance_rhs_simp_tlm(dph(:,:,:,t), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(etaz,cp), real(etam,cp), real(rdx,cp), drhs(:,:,:,t))
    end do

    dimb = dlhs - drhs
  end function diag_ImbalSimpTLM

  ! ===========================================================================

  function diag_ImbalSimpAdj(ncid, nx, ny, nz, nt) result (adj)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: adj

    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps, dps, adjpsi
    real(cp), dimension(nx,ny,0:nz,nt) :: dph, adjphi
    real(cp), dimension(nx,ny,nz,nt) :: dimb
    real, dimension(nx,ny,nz,nt) :: dlhs, drhs
    real, dimension(nx,ny,nt) :: mCap
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: mfm, f
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real    :: rdx, norm
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)

    dps = 0.
    dps(5,6,:,:) = 5e6

    mCap = diag_MCap(ncid, nx, ny, nt)   
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, "RDX", rdx)
    
    do t = 1, nt
       call balance_lhs_tlm(dps(:,:,:,t), ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(mfu,cp), real(mfv,cp), real(rdx,cp), dlhs(:,:,:,t))
    end do

    dph = 0.
    dph(5,6,:,:) = 50.

    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    
    do t = 1, nt
       call balance_rhs_simp_tlm(dph(:,:,:,t), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(etaz,cp), real(etam,cp), real(rdx,cp), drhs(:,:,:,t))
    end do

    dimb = dlhs - drhs
    print *, real(sum(real(dimb,cp)**2))
    
    do t = 1, nt
       call imbalance_simp_adj(dimb(:,:,:,t), ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mfm,cp), real(mfu,cp), real(mfv,cp), real(etaz,cp), real(etam,cp),  &
            real(rdx,cp), adjpsi(:,:,:,t), adjphi(:,:,:,t))
    end do

    print *, real(sum(dps*real(adjpsi,cp))+sum(dph*real(adjphi,cp)))
    adj = 0
  end function diag_ImbalSimpAdj

  ! ==============================================================
  
  function diag_thetaFromPhi(ncid, nx, ny, nz, nt) result (th)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: th

    real, dimension(nx,ny,0:nz,nt) :: phi
    real, dimension(nx,ny,nz,nt) :: potTemp
    real, dimension(nx,ny,nt) :: mu
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real :: pt
    integer :: i, j, k, t

    potTemp = diag_Theta(ncid, nx, ny, nz, nt)
    phi = diag_GeopS(ncid, nx, ny, nz, nt)
    mu = diag_DryAir(ncid, nx, ny, nt)
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz) 
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "P_TOP", pt)
    
    call theta_phi(phi, mu, etaz, etam, pt, th)
    
    i = 5
    j = 6
    k = 11
    t = 13

    print *, phi(i,j,k,t), phi(i,j,k-1,t)
    print *, etaz(k-1), etam(k), etaz(k)
    print *, mu(i,j,t)
    print *, potTemp(i,j,k,t), th(i,j,k,t)

  end function diag_thetaFromPhi

  ! ==============================================================

  function diag_PotVort(ncid, nx, ny, nz, nt) result (q)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: q

    real, dimension(nx+1,ny,nz,nt) :: u
    real, dimension(nx,ny+1,nz,nt) :: v
    real, dimension(nx,ny,nz+1,nt) :: ht
    real, dimension(nx,ny,nz,nt) :: p, th, qh2o
    real, dimension(nx,ny) :: mf, f
    real, dimension(nz+1) :: etaz
    real, dimension(nz) :: etam
    real :: rdx, dx
    
    call read_values(ncid, wrfVar(UStag)%wrfName, u)
    call read_values(ncid, wrfVar(VStag)%wrfName, v)    
    ht = diag_GeopHghtS(ncid, nx, ny, nz, nt)
    p = diag_Pressure(ncid, nx, ny, nz, nt, si=.true.)
    th = diag_Theta(ncid, nx, ny, nz, nt)
    call read_values(ncid, wrfVar(MixR)%wrfName, qh2o)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mf)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)    
    call read_values_const(ncid, "RDX", rdx)
    
    dx = 1./rdx
    call pvor(u, v, ht, th, density(p, temp_from_theta_p(th, p), qh2o), f, mf, etam, etaz, dx, q)
  end function diag_PotVort
  
  ! ===========================================================

  function diag_PVExp(ncid, nx, ny, nz, nt) result (pv)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: pv
    
    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps
    real, dimension(nx,ny,0:nz,nt) :: ph
    real, dimension(nx,ny,nt) :: mCap, mu
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: f
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real    :: pt, rdx
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)
    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    mCap = diag_MCap(ncid, nx, ny, nt)
    mu = diag_DryAir(ncid, nx, ny, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "P_TOP", pt)
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
       print *, t
       call pv_exp(ps(:,:,:,t), real(ph(:,:,:,t),cp), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mu(:,:,t),cp), real(mfu,cp), real(mfv,cp), real(etaz,cp), real(etam,cp),  &
            real(pt,cp), real(rdx,cp), pv(:,:,:,t))
    end do
  end function diag_PVExp

  ! ===========================================================

  function diag_PVPert(ncid, nx, ny, nz, nt) result (pv)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: pv

    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps
    real, dimension(nx,ny,0:nz,nt) :: ph
    real, dimension(nx,ny,nt) :: mCap, mu
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: f
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real    :: pt, rdx
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)
    ps(5,6,10,:) = ps(5,6,10,:) + 5e6

    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    ph(5,6,10,:) = ph(5,6,10,:) + 50.

    mCap = diag_MCap(ncid, nx, ny, nt)
    mu = diag_DryAir(ncid, nx, ny, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "P_TOP", pt)
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
       call pv_exp(ps(:,:,:,t), real(ph(:,:,:,t),cp), real(f,cp), real(mCap(:,:,t),cp),  &
            real(mu(:,:,t),cp), real(mfu,cp), real(mfv,cp), real(etaz,cp), real(etam,cp),  &
            real(pt,cp), real(rdx,cp), pv(:,:,:,t))
    end do
  end function diag_PVPert

  ! ===========================================================================

  function diag_PVTLM(ncid, nx, ny, nz, nt) result (dpv)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: dpv

    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps, dps
    real(cp), dimension(nx,ny,0:nz,nt) :: dph
    real, dimension(nx,ny,0:nz,nt) :: ph
    real, dimension(nx,ny,nt) :: mCap, mu
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: f
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real :: pt, rdx
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)
    dps = 0.
    dps(5,6,10,:) = 5e6

    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    dph = 0.
    dph(5,6,10,:) = 50.

    mCap = diag_MCap(ncid, nx, ny, nt)
    mu = diag_DryAir(ncid, nx, ny, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "P_TOP", pt)
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
       call pv_tlm(dps(:,:,:,t), dph(:,:,:,t), ps(:,:,:,t), real(ph(:,:,:,t),cp), real(f,cp),  &
            real(mCap(:,:,t),cp), real(mu(:,:,t),cp), real(mfu,cp), real(mfv,cp), real(etaz,cp),  &
            real(etam,cp), real(pt,cp), real(rdx,cp), dpv(:,:,:,t))
    end do
  end function diag_PVTLM

  ! ===========================================================================

  function diag_PVAdj(ncid, nx, ny, nz, nt) result (adj)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: adj

    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps, dps, adjpsi
    real(cp), dimension(nx,ny,0:nz,nt) :: dph, adjphi
    real, dimension(nx,ny,0:nz,nt) :: ph
    real, dimension(nx,ny,nz,nt) :: dpv
    real, dimension(nx,ny,nt) :: mCap, mu
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: f
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real :: pt, rdx
    integer :: i, j, k, t
    
    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)
    dps = 0.
    dps(5,6,10,:) = 5e6

    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    dph = 0.
    dph(5,6,10,:) = 50.

    mCap = diag_MCap(ncid, nx, ny, nt)
    mu = diag_DryAir(ncid, nx, ny, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "P_TOP", pt)
    call read_values_const(ncid, "RDX", rdx)

    do t = 1, nt
       call pv_tlm(dps(:,:,:,t), dph(:,:,:,t), ps(:,:,:,t), real(ph(:,:,:,t),cp), real(f,cp),  &
            real(mCap(:,:,t),cp), real(mu(:,:,t),cp), real(mfu,cp), real(mfv,cp), real(etaz,cp),  &
            real(etam,cp), real(pt,cp), real(rdx,cp), dpv(:,:,:,t))
    end do
    print *, real(sum(real(dpv,cp)**2))

    do t = 1, nt
       call pv_adj(real(dpv(:,:,:,t),cp), ps(:,:,:,t), real(ph(:,:,:,t),cp), real(f,cp),  &
            real(mCap(:,:,t),cp), real(mu(:,:,t),cp), real(mfu,cp), real(mfv,cp), real(etaz,cp),  &
            real(etam,cp), real(pt,cp), real(rdx,cp), adjpsi(:,:,:,t), adjphi(:,:,:,t))
    end do
    print *, real(sum(dps*real(adjpsi,cp))+sum(dph*real(adjphi,cp)))
    adj = 0
  end function diag_PVAdj

  ! ===========================================================================

  function diag_CostGradPhi(ncid, nx, ny, nz, nt) result (grad)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,0:nz,nt) :: gradgem

    real(cp), parameter :: e = 5000.

    real(cp), dimension(0:nx,0:ny,nz,nt) :: ps, gradpsi1, gradpsi2
    real(cp), dimension(nx,ny,0:nz,nt) :: gradphi1, gradphi2, grad
    real(cp), dimension(nx,ny,nz,nt) :: imb
    real(cp), dimension(0:nx,0:ny,nz) :: ps3
    real(cp), dimension(nx,ny,0:nz) :: ph3
    real, dimension(nx,ny,0:nz,nt) :: ph
    real, dimension(nx,ny,nz,nt) :: pv, lhs, rhs, pvGiven
    real, dimension(nx,ny,nt) :: mCap, mu
    real, dimension(0:nx,ny)  :: mfu
    real, dimension(nx,0:ny)  :: mfv
    real, dimension(nx,ny)    :: mfm, f
    real, dimension(0:nz)     :: etaz
    real, dimension(nz)       :: etam
    real(cp) :: costa, costb, cost
    real :: pt, rdx
    integer :: i, j, k, t

    open(33, file="ps.dat", status="old", form="unformatted", action="read")
    read (33) i, j, k, t
    print *, i, j, k, t
    do t = 1, nt
       do k = 1, nz
          read (33) ps(0:nx,0:ny,k,t)
       end do
    end do
    close(33)
    ps(0,0,:,:) = Half * (ps(1,0,:,:) + ps(0,1,:,:))
    ps(nx,0,:,:) = Half * (ps(nx-1,0,:,:) + ps(nx,1,:,:))
    ps(0,ny,:,:) = Half * (ps(0,ny-1,:,:) + ps(1,ny,:,:))
    ps(nx,ny,:,:) = Half * (ps(nx-1,ny,:,:) + ps(nx,ny-1,:,:))

    ph = diag_GeopS(ncid, nx, ny, nz, nt)
    pvGiven = diag_PVExp(ncid, nx, ny, nz, nt)

    mCap = diag_MCap(ncid, nx, ny, nt)
    mu = diag_DryAir(ncid, nx, ny, nt)
    call read_values_const(ncid, wrfVar(MapFacU)%wrfName, mfu)
    call read_values_const(ncid, wrfVar(MapFacV)%wrfName, mfv)
    call read_values_const(ncid, wrfVar(MapFacM)%wrfName, mfm)    
    call read_values_const(ncid, wrfVar(Coriolis)%wrfName, f)
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, etaz)    
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, etam)
    call read_values_const(ncid, "P_TOP", pt)
    call read_values_const(ncid, "RDX", rdx)

    ! Add some perturbations to ps, ph, and pv
    do t = 1, nt
       do k = 2, nz-1
          do j = 1, ny-1
             do i = 1, nx-1
                ps(i,j,k,t) = (i*ps(nx,j,k,t) + (nx-i)*ps(0,j,k,t) + j*ps(i,ny,k,t)  &
                     + (ny-j)*ps(i,0,k,t)) / real(nx+ny)
!                ps(i,j,k,t) = Half * ( (i*ps(nx,j,k,t)+(nx-i)*ps(0,j,k,t))/real(nx)  &
!                     + (j*ps(i,ny,k,t)+(ny-j)*ps(i,0,k,t))/real(ny))
             end do
          end do
       end do
    end do
!    print *, ps(:,10,13,13)
!    print *
!    print *, ps(4,:,13,13)
    !ps(1:nx-1,1:ny-1,2:nz-1,:) = .9 * ps(1:nx-1,1:ny-1,2:nz-1,:)

!    stop
    t = 13

!!$    call balance_lhs_exp(ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp), real(mfm,cp),  &
!!$         real(mfu,cp), real(mfv,cp), real(rdx,cp), lhs(:,:,:,t))
!!$
!!$    call balance_rhs_simp_exp(real(ph(:,:,:,t),cp), real(mCap(:,:,t),cp), real(mfm,cp),  &
!!$         real(etaz,cp), real(etam,cp), real(rdx,cp), rhs(:,:,:,t))

    ! imb = lhs - rhs
!    print *, sum((real(lhs(:,:,:,t),cp)-real(rhs(:,:,:,t),cp))**2)
    call imbalance_exp(ps(:,:,:,t), real(ph(:,:,:,t),cp), real(f,cp), real(mCap(:,:,t),cp),  &
         real(mfm,cp), real(mfu,cp), real(mfv,cp), real(etaz,cp), real(etam,cp), real(rdx,cp),  &
         imb(:,:,:,t))
!    print *, sum(imb**2)
!    costb = Half * e * sum((real(lhs(:,:,:,t),cp)-real(rhs(:,:,:,t),cp))**2)
    costb = Half * e * sum(imb**2)

    call pv_exp(ps(:,:,:,t), real(ph(:,:,:,t),cp), real(f,cp), real(mCap(:,:,t),cp),  &
         real(mu(:,:,t),cp), real(mfu,cp), real(mfv,cp), real(etaz,cp), real(etam,cp),  &
         real(pt,cp), real(rdx,cp), pv(:,:,:,t))

    costa = Half * sum((real(pv(:,:,:,t),cp)-real(pvGiven(:,:,:,t),cp))**2)

    cost = costa + costb
    write(*,*) "Cost A & B:", costa, costb
    write(*,*) "Total cost:", cost

    call pv_adj(real(pv(:,:,:,t),cp)-real(pvGiven(:,:,:,t),cp), ps(:,:,:,t),  &
         real(ph(:,:,:,t),cp), real(f,cp), real(mCap(:,:,t),cp), real(mu(:,:,t),cp),  &
         real(mfu,cp), real(mfv,cp), real(etaz,cp), real(etam,cp), real(pt,cp), real(rdx,cp),  &
         gradpsi1(:,:,:,t), gradphi1(:,:,:,t))
!!$    call imbalance_simp_adj(real(lhs(:,:,:,t),cp)-real(rhs(:,:,:,t),cp), ps(:,:,:,t), real(f,cp), &
!!$         real(mCap(:,:,t),cp), real(mfm,cp), real(mfu,cp), real(mfv,cp), real(etaz,cp),  &
!!$         real(etam,cp), real(rdx,cp), gradpsi2(:,:,:,t), gradphi2(:,:,:,t))
    call imbalance_simp_adj(imb(:,:,:,t), ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp),  &
         real(mfm,cp), real(mfu,cp), real(mfv,cp), real(etaz,cp), real(etam,cp), real(rdx,cp),  &
         gradpsi2(:,:,:,t), gradphi2(:,:,:,t))

    grad(:,:,:,t) = gradphi1(:,:,:,t) + e * gradphi2(:,:,:,t)

    ! Check the gradient
    print *, "Gradient phi 1 & 2: ", gradphi1(5,6,10,t), e * gradphi2(5,6,10,t)
    print *, "Total phi gradient: ", grad(5,6,10,t)
!    ph(5,6,10,t) = ph(5,6,10,t) + .1

!!$    call balance_lhs_exp(ps(:,:,:,t), real(f,cp), real(mCap(:,:,t),cp), real(mfm,cp),  &
!!$         real(mfu,cp), real(mfv,cp), real(rdx,cp), lhs(:,:,:,t))
!!$    call balance_rhs_simp_exp(real(ph(:,:,:,t),cp), real(mCap(:,:,t),cp), real(mfm,cp),  &
!!$         real(etaz,cp), real(etam,cp), real(rdx,cp), rhs(:,:,:,t))
    call imbalance_exp(ps(:,:,:,t), real(ph(:,:,:,t),cp), real(f,cp), real(mCap(:,:,t),cp),  &
         real(mfm,cp), real(mfu,cp), real(mfv,cp), real(etaz,cp), real(etam,cp), real(rdx,cp),  &
         imb(:,:,:,t))
    ! imb = lhs - rhs
!    costb = Half * e * sum((real(lhs(:,:,:,t),cp)-real(rhs(:,:,:,t),cp))**2)
    costb = Half * e * sum(imb**2)

    call pv_exp(ps(:,:,:,t), real(ph(:,:,:,t),cp), real(f,cp), real(mCap(:,:,t),cp),  &
         real(mu(:,:,t),cp), real(mfu,cp), real(mfv,cp), real(etaz,cp), real(etam,cp),  &
         real(pt,cp), real(rdx,cp), pv(:,:,:,t))

    costa = Half * sum((real(pv(:,:,:,t),cp)-real(pvGiven(:,:,:,t),cp))**2)

    cost = costa + costb
    print *, "After perturbation:"
    write(*,*) "Cost A & B:", costa, costb
    write(*,*) "Total cost:", cost

    gradgem = grad

    ps3 = ps(:,:,:,t)
    ph3 = ph(:,:,:,t)
    print *, real(ps3(nx-1,ny-1,11)), real(ph3(nx-1,ny-1,11))   
    call solve_phi_psi(real(pvGiven(:,:,:,t),cp), ps3, ph3, real(f,cp),  &
         real(mCap(:,:,t),cp), real(mu(:,:,t),cp), real(mfm,cp), real(mfu,cp), real(mfv,cp),  &
         real(etaz,cp), real(etam,cp), real(pt,cp), real(rdx,cp))
    print *, real(ps3(nx-1,ny-1,11)), real(ph3(nx-1,ny-1,11))
  end function diag_CostGradPhi

  ! ===========================================================================
  
  function diag_InvertPhi(ncid, nx, ny, nz, nt) result (phiout)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(nx,ny,nz,nt) :: phiout

    integer :: i, j, k, t

    real(cp), dimension(0:nx,0:ny,nz) :: psi
    real(cp), dimension(nx,ny,0:nz) :: phi    
    real(cp), dimension(nx,ny,0:nz,nt) :: phi4    
    real, dimension(nx,ny,nt) :: mu
    real, dimension(nz+1) :: znw
    real, dimension(nz) :: znu
    real :: pt

    open(34, file="psphiinvert.dat", status="old", form="unformatted", action="read")
    read (34) i, j, k
    print *, i, j, k
    read (34) psi
    read (34) phi
    close(34)

    mu = diag_DryAir(ncid, nx, ny, nt)
    call read_values_const(ncid, wrfVar(EtaW)%wrfName, znw)
    call read_values_const(ncid, wrfVar(EtaMass)%wrfName, znu)
    call read_values_const(ncid, "P_TOP", pt)

    phi4 = spread(phi, 4, nt)
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             phiout(i,j,:,t) =  &
                  unstag_hght(real(phi4(i,j,:,t)), znw, znu, mu(i,j,t), pt)
          end do
       end do
    end do
  end function diag_InvertPhi

  ! ===========================================================================
  
  function diag_InvertPsi(ncid, nx, ny, nz, nt) result (psiout)
    integer, intent(in)          :: ncid, nx, ny, nz, nt
    real, dimension(0:nx,0:ny,nz,nt) :: psiout

    integer :: i, j, k
    
    real(cp), dimension(0:nx,0:ny,nz,nt) :: psi4   
    real(cp), dimension(0:nx,0:ny,nz) :: psi
    real(cp), dimension(nx,ny,0:nz) :: phi    

    open(34, file="psphiinvert.dat", status="old", form="unformatted", action="read")
    read (34) i, j, k
    print *, i, j, k
    read (34) psi
    read (34) phi
    close(34)
    
    psi4 = spread(psi, 4, nt)
    psiout = Quarter * (psi4(:nx-1,:ny-1,:,:) + psi4(1:,:ny-1,:,:) + psi4(:nx-1,1:,:,:)  &
         + psi4(1:,1:,:,:))
  end function diag_InvertPsi
end module registry
