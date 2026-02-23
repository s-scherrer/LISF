!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_getdailyfaparpred
! \label{noahmp401_getdailyfaparpred}
!
! !REVISION HISTORY:
! 13 Feb 2020: Sujay Kumar; Initial Specification
! 25 Oct 2025: Samuel Scherrer; adapted from da_LAI
!
! !INTERFACE:
subroutine noahmp401_getdailyfaparpred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use noahmp401_lsmMod
  use noahmp401_dasoilm_Mod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the daily integrated total FAPAR obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  real                   :: obs_tmp
  integer                :: i,t,m,gid,kk
  real                   :: inputs_tp(6)
  character*50           :: units_tp(6)
  real                   :: fapar(LIS_rc%npatch(n,LIS_rc%lsm_index))


  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     fapar(t) = noahmp401_struc(n)%noahmp401(t)%daily_fapar
  enddo
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       fapar,&
       obs_pred)
  
end subroutine noahmp401_getdailyfaparpred


subroutine noahmp401_calculate_instantaneous_fapar(n, k, fpartype, fpar)
    ! !USES:
    use ESMF
    use LIS_constantsMod
    use LIS_coreMod
    use LIS_dataAssimMod
    use LIS_DAobservationsMod
    use noahmp401_lsmMod
    use noahmp401_dasoilm_Mod

    implicit none
    ! !ARGUMENTS: 
    integer, intent(in)    :: n
    integer, intent(in)    :: k
    integer, intent(in)    :: fpartype  ! 1: black-sky, 2: white-sky, else: total
    real, intent(out)      :: wsfpar(LIS_rc%npatch(n,LIS_rc%lsm_index))

    real                   :: cosz
    real                   :: fage, elai, esai, fsno
    real                   :: albold, tauss
    real                   :: bgap, wgap
    real, dimension(1:2)   :: albgrd, albgri, albd, albi
    real, dimension(1:2)   :: fabd, fabi
    real, dimension(1:2)   :: ftdd, ftid, ftii
    real, dimension(1:2)   :: fsun
    real, dimension(1:2)   :: frevd, frevi, fregd, fregi
    integer                :: iloc, jloc, ist, ice
    real                   :: lat, lon
    type(noahmp_parameters) :: param

    ice = 0  ! no ice
    ist = 1  ! land
    iloc = 0 ! dummy - not used
    jloc = 0 ! dummy - not used

    do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)

        param = noahmp401_struc(n)%noahmp401(t)%param

        call calc_elai_esai(param, &
             noahmp401_struc(n)%noahmp401(t)%lai,&
             noahmp401_struc(n)%noahmp401(t)%sai,&
             noahmp401_struc(n)%noahmp401(t)%snowh,&
             elai, esai)

        call calc_fsno(param,&
             noahmp401_struc(n)%noahmp401(t)%sneqv,&
             noahmp401_struc(n)%noahmp401(t)%snowh,&
             fsno)

        albold = noahmp401_struc(n)%noahmp401(t)%albold
        tauss = noahmp401_struc(n)%noahmp401(t)%tauss

        call albedo(&
             param,&                                     ! parameters
             noahmp401_struc(n)%noahmp401(t)%vegetype,&  ! VEGTYP: not used
             ist,&                                       ! IST: surface type = land
             ice,&                                       ! ICE: no ice
             noahmp401_struc(n)%nsoil,&                  ! NSOIL
             noahmp401_struc(n)%ts,&                     ! DT
             noahmp401_struc(n)%noahmp401(t)%cosz,&      ! COSZ
             fage,&                                      ! FAGE (dummy value)
             elai,&                                      ! ELAI
             esai,&                                      ! ESAI
             noahmp401_struc(n)%noahmp401(t)%tg,&        ! TG
             noahmp401_struc(n)%noahmp401(t)%tv,&        ! TV
             noahmp401_struc(n)%noahmp401(t)%snowh,&     ! SNOWH
             fsno,&                                      ! FSNO
             noahmp401_struc(n)%noahmp401(t)%fwet,&      ! FWET
             noahmp401_struc(n)%noahmp401(t)%smc,&       ! SMC
             noahmp401_struc(n)%noahmp401(t)%sneqvo,&    ! SNEQVO
             noahmp401_struc(n)%noahmp401(t)%sneqv,&     ! SNEQV
             noahmp401_struc(n)%noahmp401(t)%qsnow,&     ! QSNOW
             noahmp401_struc(n)%noahmp401(t)%fveg,&      ! FVEG
             iloc, jloc,&
             albold, tauss,&                             ! inout
             albgrd, albgri, albd, albi,fabd,fabi,&      ! out
             ftdd,ftid,ftii,fsun,frevi,frevd,fregd,&     ! out
             fregi,bgap,wgap&                            ! out
             )

        if (fpartype.eq.1) then
            fpar(t) = fabd(1)  ! direct -> black sky
        else if (fpartype.eq.2) then
            fpar(t) = fabi(1)  ! diffuse -> white sky
        else
            ! total instantaneous FAPAR = PSAV / PAR
            ! with
            !  PSAV = CAD(1) + CAI(1)
            !  CAD(1) = SOLAD(1) * FABD(1)
            !  CAI(1) = SOLAI(1) * FABI(1)
            !  PAR = SOLAD(1) + SOLAI(1)    i.e. total incoming visible radiation
            !  SOLAD(1) = SWDOWN * 0.7 * 0.5
            !  SOLAI(1) = SWDOWN * 0.3 * 0.5
            ! therefore:
            !  PAR = SWDOWN * 0.5
            !  PSAV = SWDOWN * 0.5 * (0.7 * FABD(1) + 0.3 * FABI(1))
            !  FAPAR = 0.7 * FABD(1) + 0.3 * FABI(1)
            fpar(t) = 0.7 * fabd(1) + 0.3 * fabi(1)
        endif
    enddo

contains

    subroutine calc_elai_esai(parameters, lai, sai, snowh, elai, esai)
        type (noahmp_parameters), intent(in) :: parameters
        real, intent(in)                     :: lai
        real, intent(in)                     :: sai
        real, intent(in)                     :: snowh
        real, intent(out)                    :: elai
        real, intent(out)                    :: esai

        real :: snowhc, fb
        integer :: croptype

        croptype = 0  ! Noah-MP crop is not supported

        ! copied from noahmp.4.0.1/phys/module_sf_noahmplsm_401.F90
        IF(parameters%HVT> 0. .AND. parameters%HVT <= 1.0) THEN          !MB: change to 1.0 and 0.2 to reflect
            SNOWHC = parameters%HVT*EXP(-SNOWH/0.2)             !      changes to HVT in MPTABLE
            IF(SNOWHC>1.E-06) THEN      !Wanshu: avoid very small SNOWHC induced floating invalid
                FB     = MIN(SNOWH,SNOWHC)/SNOWHC
            ELSE
                !print *,"small snowh in phenology=",snowh
                FB = 1
            END IF
        ENDIF

        ELAI =  LAI*(1.-FB)
        ESAI =  SAI*(1.-FB)
        IF (ESAI < 0.05 .and. CROPTYPE == 0) ESAI = 0.0                   ! MB: ESAI CHECK, change to 0.05 v3.6
        IF ((ELAI < 0.05 .OR. ESAI == 0.0) .and. CROPTYPE == 0) ELAI = 0.0  ! MB: LAI CHECK

    end subroutine calc_elai_esai


    subroutine calc_fsno(parameters, sneqv, snowh, fsno)
        type (noahmp_parameters), intent(in) :: parameters
        real, intent(in)                     :: sneqv
        real, intent(in)                     :: snowh
        real, intent(out)                    :: fsno

        real :: bdsno, fmelt

        REAL, PARAMETER                   :: Z0     = 0.002  ! Bare-soil roughness length (m) (i.e., under the canopy)

        ! copied from noahmp.4.0.1/phys/module_sf_noahmplsm_401.F90
        FSNO = 0.
        IF(SNOWH.GT.0.)  THEN       
            BDSNO    = SNEQV / SNOWH
            FMELT    = (BDSNO/100.)**parameters%MFSNO
            if (FMELT<0.000001) then !Bailing Li, added this for GRACE DA to catch smaller values
                FSNO = 1
                !print *,"small FMELT due to snowh,fmelt,bdsno,para_mfsno",snowh,fmelt,bdsno,parameters%MFSNO
            else
                FSNO     = TANH( SNOWH /(2.5* Z0 * FMELT))
            end if
        ENDIF
    end subroutine calc_fsno

end subroutine noahmp401_calculate_instantaneous_fapar
    


! !ROUTINE: noahmp401_getinstfaparpred
! \label{noahmp401_getinstfaparpred}
!
! !REVISION HISTORY:
! 13 Feb 2020: Sujay Kumar; Initial Specification
! 25 Oct 2025: Samuel Scherrer; adapted from da_LAI
!
! !INTERFACE:
subroutine noahmp401_getinstfaparpred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use noahmp401_lsmMod
  use noahmp401_dasoilm_Mod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the instantaneous total FAPAR obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  real                   :: obs_tmp
  integer                :: i,t,m,gid,kk
  real                   :: inputs_tp(6)
  character*50           :: units_tp(6)
  real                   :: fapar(LIS_rc%npatch(n,LIS_rc%lsm_index))


  call noahmp401_calculate_instantaneous_fapar(n, k, 3, fapar)

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       fapar,&
       obs_pred)
  
end subroutine noahmp401_getinstfaparpred

! !ROUTINE: noahmp401_getinstbsfaparpred
! \label{noahmp401_getinstbsfaparpred}
!
! !REVISION HISTORY:
! 13 Feb 2020: Sujay Kumar; Initial Specification
! 25 Oct 2025: Samuel Scherrer; adapted from da_LAI
!
! !INTERFACE:
subroutine noahmp401_getinstbsfaparpred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use noahmp401_lsmMod
  use noahmp401_dasoilm_Mod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the instantaneous total FAPAR obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  real                   :: obs_tmp
  integer                :: i,t,m,gid,kk
  real                   :: inputs_tp(6)
  character*50           :: units_tp(6)
  real                   :: fapar(LIS_rc%npatch(n,LIS_rc%lsm_index))


  call noahmp401_calculate_instantaneous_fapar(n, k, 1, fapar)

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       fapar,&
       obs_pred)
  
end subroutine noahmp401_getinstbsfaparpred

! !ROUTINE: noahmp401_getinstwsfaparpred
! \label{noahmp401_getinstwsfaparpred}
!
! !REVISION HISTORY:
! 13 Feb 2020: Sujay Kumar; Initial Specification
! 25 Oct 2025: Samuel Scherrer; adapted from da_LAI
!
! !INTERFACE:
subroutine noahmp401_getinstwsfaparpred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use noahmp401_lsmMod
  use noahmp401_dasoilm_Mod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the instantaneous total FAPAR obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  real                   :: obs_tmp
  integer                :: i,t,m,gid,kk
  real                   :: inputs_tp(6)
  character*50           :: units_tp(6)
  real                   :: fapar(LIS_rc%npatch(n,LIS_rc%lsm_index))


  call noahmp401_calculate_instantaneous_fapar(n, k, 2, fapar)

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       fapar,&
       obs_pred)
  
end subroutine noahmp401_getinstwsfaparpred
