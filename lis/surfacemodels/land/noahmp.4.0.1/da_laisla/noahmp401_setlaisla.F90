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
! !ROUTINE: noahmp401_setlaisla
!  \label{noahmp401_setlaisla}
!
! !REVISION HISTORY:
! 10 Jan 2023: Samuel Scherrer; Initial Specification based on da_LAI
! 
! !INTERFACE:
subroutine noahmp401_setlaisla(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp401_lsmMod


  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP
  type(ESMF_Field)       :: laiField,slaField
  integer                :: t
  integer                :: status
  real, pointer          :: lai(:),sla(:)
  real                   :: lfmass

  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_State,"Specific leaf area",slaField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(slaField,localDE=0,farrayPtr=sla,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(NOAHMP401_struc(n)%noahmp401(t)%param%sla.ne.0) then 
        NOAHMP401_struc(n)%noahmp401(t)%lai = lai(t)
        NOAHMP401_struc(n)%noahmp401(t)%param%sla = sla(t)
        lfmass = lai(t)*1000.0/sla(t)
        NOAHMP401_struc(n)%noahmp401(t)%lfmass = lfmass
     endif
  enddo

end subroutine noahmp401_setlaisla

