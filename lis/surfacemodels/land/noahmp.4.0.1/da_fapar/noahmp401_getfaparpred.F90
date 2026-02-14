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


  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     fapar(t) = noahmp401_struc(n)%noahmp401(t)%fapar
  enddo
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


  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     fapar(t) = noahmp401_struc(n)%noahmp401(t)%bsfpar
  enddo
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


  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     fapar(t) = noahmp401_struc(n)%noahmp401(t)%bsfpar
  enddo
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       fapar,&
       obs_pred)
  
end subroutine noahmp401_getinstwsfaparpred