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
! !ROUTINE: noahmp401_update_lai_maxlai
!  \label{noahmp401_update_lai_maxlai}
!
! !REVISION HISTORY:
! 13 Feb 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp401_update_lai_maxlai(n, LSM_State, LSM_Incr_State)
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
  type(ESMF_Field)       :: laiField, laiIncrField
  type(ESMF_Field)       :: maxlaiField, maxlaiIncrField

  integer                :: t,gid
  integer                :: status
  real, pointer          :: lai(:), laiincr(:)
  real, pointer          :: maxlai(:), maxlaiincr(:)
  real                   :: laitmp,laimax,laimin
  real                   :: maxlaitmp,maxlaimax,maxlaimin

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: laimean(LIS_rc%ngrid(n))
  integer                :: nlaimean(LIS_rc%ngrid(n))
  real                   :: maxlaimean(LIS_rc%ngrid(n))
  integer                :: nmaxlaimean(LIS_rc%ngrid(n))

 
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

  call ESMF_StateGet(LSM_State,"MAXLAI",maxlaiField,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"MAXLAI",maxlaiIncrField,rc=status)
  call LIS_verify(status)

 
  call ESMF_FieldGet(maxlaiField,localDE=0,farrayPtr=maxlai,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(maxlaiIncrField,localDE=0,farrayPtr=maxlaiincr,rc=status)
  call LIS_verify(status)


  call ESMF_AttributeGet(maxlaiField,"Max Value",maxlaimax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(maxlaiField,"Min Value",maxlaimin,rc=status)
  call LIS_verify(status)

  update_flag    = .true.
  perc_violation = 0.0
  laimean       = 0.0
  nlaimean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     laitmp =  lai(t) + laiincr(t)
     maxlaitmp = maxlai(t) + maxlaiincr(t)


     if(laitmp.lt.laimin.or.laitmp.gt.laimax.or.maxlaitmp.lt.maxlaimin.or.maxlaitmp.gt.maxlaimax) then
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
     if(.not.update_flag(gid)) then
        if(perc_violation(gid).lt.0.8) then
           if((lai(t)+laiincr(t).gt.laimin).and.&
                (lai(t)+laiincr(t).lt.laimax)) then 
              laimean(gid) = laimean(gid) + &
                   lai(t) + laiincr(t)
              nlaimean(gid) = nlaimean(gid) + 1
           endif
           if((maxlai(t)+maxlaiincr(t).gt.maxlaimin).and.&
                (maxlai(t)+maxlaiincr(t).lt.maxlaimax)) then 
              maxlaimean(gid) = maxlaimean(gid) + &
                   maxlai(t) + maxlaiincr(t)
              nmaxlaimean(gid) = nmaxlaimean(gid) + 1
           endif
        endif
     endif
  enddo

 do gid=1,LIS_rc%ngrid(n)
     if(nlaimean(gid).gt.0) then
        laimean(gid) = laimean(gid)/nlaimean(gid)
     endif
     if(nmaxlaimean(gid).gt.0) then
        maxlaimean(gid) = maxlaimean(gid)/nmaxlaimean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     laitmp =  lai(t) + laiincr(t)
     maxlaitmp =  maxlai(t) + maxlaiincr(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        lai(t) = laitmp
        maxlai(t) = maxlaitmp
     elseif(perc_violation(gid).lt.0.8) then
        if(laitmp.lt.laimin.or.laitmp.gt.laimax) then
           lai(t) = laimean(gid)
        else
           lai(t) = lai(t) + laiincr(t)
        endif
        if(maxlaitmp.lt.maxlaimin.or.maxlaitmp.gt.maxlaimax) then
           maxlai(t) = maxlaimean(gid)
        else
           maxlai(t) = maxlai(t) + maxlaiincr(t)
        endif
     endif
  enddo

end subroutine noahmp401_update_lai_maxlai

