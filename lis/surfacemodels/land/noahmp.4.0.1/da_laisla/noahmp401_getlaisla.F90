!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_getlaisla
! \label{noahmp401_getlaisla}
!
! !REVISION HISTORY:
! 10 Jan 2023: Samuel Scherrer; Initial Specification
!
! !INTERFACE:
subroutine noahmp401_getlaisla(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
  !EOP

  
  type(ESMF_Field)       :: laiField, slaField

  integer                :: t
  integer                :: status
  real, pointer          :: lai(:)
  real, pointer          :: sla(:)
 
  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_State,"Specific leaf area",slaField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(slaField,localDE=0,farrayPtr=sla,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lai(t) = NOAHMP401_struc(n)%noahmp401(t)%lai
     sla(t) = NOAHMP401_struc(n)%noahmp401(t)%param%sla
  enddo

end subroutine noahmp401_getlaisla

