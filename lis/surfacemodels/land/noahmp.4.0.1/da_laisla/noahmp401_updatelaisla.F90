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
! !ROUTINE: noahmp401_updatelaisla
!  \label{noahmp401_updatelaisla}
!
! !REVISION HISTORY:
! 10 Jan 2023: Samuel Scherrer; Initial Specification based on da_LAI
!
! !INTERFACE:
subroutine noahmp401_updatelaisla(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP
  type(ESMF_Field)       :: laiField, laiIncrField, slaField, slaIncrField

  integer                :: t,gid
  integer                :: status
  real, pointer          :: lai(:), laiincr(:), sla(:), slaincr(:)
  real                   :: laitmp,laimax,laimin,slatmp

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: laimean(LIS_rc%ngrid(n))
  real                   :: slamean(LIS_rc%ngrid(n))
  integer                :: nlaimean(LIS_rc%ngrid(n))

 
  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"LAI",laiIncrField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(laiIncrField,localDE=0,farrayPtr=laiincr,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(laiField,"Max Value",laimax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(laiField,"Min Value",laimin,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_State,"Specific leaf area",slaField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Specific leaf area",slaIncrField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(slaField,localDE=0,farrayPtr=sla,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(slaIncrField,localDE=0,farrayPtr=slaincr,rc=status)
  call LIS_verify(status)

  update_flag    = .true.
  perc_violation = 0.0
  laimean       = 0.0
  nlaimean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     laitmp = lai(t) + laiincr(t)
     slatmp = sla(t) + slaincr(t)


     if(laitmp.lt.laimin.or.laitmp.gt.laimax.or.slatmp.le.0.0) then
        update_flag(gid) = .false.
        perc_violation(gid) = perc_violation(gid) +1
     endif

  enddo

  do gid=1,LIS_rc%ngrid(n)
     perc_violation(gid) = perc_violation(gid)/LIS_rc%nensem(n)
  enddo

! For ensembles that are unphysical, compute the
! ensemble average after excluding them. This
! is done only if the majority of the ensemble
! members are good (>60%)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
     laitmp = lai(t) + laiincr(t)
     slatmp = sla(t) + slaincr(t)
     if(.not.update_flag(gid)) then
        if(perc_violation(gid).lt.0.8) then
            if((laitmp.gt.laimin).and.(laitmp.lt.laimax)&
                 .and.slatmp.gt.0.0) then 
              laimean(gid) = laimean(gid) + lai(t) + laiincr(t)
              slamean(gid) = slamean(gid) + sla(t) + slaincr(t)
              nlaimean(gid) = nlaimean(gid) + 1
           endif
        endif
     endif
  enddo

 do gid=1,LIS_rc%ngrid(n)
     if(nlaimean(gid).gt.0) then
        laimean(gid) = laimean(gid)/nlaimean(gid)
        slamean(gid) = slamean(gid)/nlaimean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     laitmp = lai(t) + laiincr(t)
     slatmp = sla(t) + slaincr(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        lai(t) = laitmp
     elseif(perc_violation(gid).lt.0.8) then
        if(laitmp.lt.laimin.or.laitmp.gt.laimax) then
           lai(t) = laimean(gid)
           sla(t) = slamean(gid)
        else
           lai(t) = laitmp
           sla(t) = slatmp
        endif
     endif
  enddo

end subroutine noahmp401_updatelaisla

