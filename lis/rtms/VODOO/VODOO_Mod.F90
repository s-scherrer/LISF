!-----------------------BEGIN---------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module VODOO_Mod
    !BOP
    !
    ! !MODULE: VODOO_Mod
    !
    ! !DESCRIPTION:
    !    This modules implements an observation operator for VOD.
    !
    !    Although it's not a real RTM it programmatically works similarly and is
    !    therefore implemented analogously to a RTM.
    !
    !    The module provides the following options for lis.config:
    !
    !    VODOO parameter file:
    !        Path to the parameter file
    !
    !    The model calculates VOD as linear function LAI and soil moisture;
    !
    !        VOD = a(LAI, RZSM) * SM' + b(LAI, RZSM)
    !
    !    with
    !
    !        SM' = relSM1 - max(0, min(1, relRZSM))
    !        relSM1 = (SM1 - WP) / (FC - WP)
    !        relRZSM = (RZSM - WP) / (FC - WP)
    !        a(LAI, RZSM) = a0 * LAI + a1 * RZSM + a2 * LAI * RZSM + a3
    !        b(LAI, RZSM) = b0 * LAI + b1 * RZSM + b2 * LAI * RZSM + b3
    !
    !    To use it, you need to provide a netCDF file with a single variable
    !    "coefs", with 2 dimensions, (ngrid, ncoef). The coef dimension should
    !    contain the coefficients in this order: a0, ..., a3, b0, ..., b3.
    !
    ! !HISTORY:
    ! 02 Mar 2022: Samuel Scherrer; initial contribution based on WCMRTM
    ! 25 Aug 2022: Samuel Scherrer; support for multiple models
    ! 17 Feb 2023: Samuel Scherrer; complete rewrite with new model structure
    !
    ! !USES:

#if (defined RTMS)

    use ESMF
    use LIS_coreMod
    use LIS_RTMMod
    use LIS_logMod

    implicit none

    PRIVATE

    !-----------------------------------------------------------------------------
    ! !PUBLIC MEMBER FUNCTIONS:
    !-----------------------------------------------------------------------------
    public :: VODOO_initialize
    public :: VODOO_f2t
    public :: VODOO_run
    public :: VODOO_output
    public :: VODOO_geometry
    !-----------------------------------------------------------------------------
    ! !PUBLIC TYPES:
    !-----------------------------------------------------------------------------
    public :: vodoo_struc
    !EOP

    ! The actual implementation of the model equations is done in the
    ! subclasses below, this is just to provide a common interface.
    ! Using instances of this type will raise an error.
    type, public ::  vodoo_type_dec
        character*256 :: parameter_fname
        !-------output------------!
        real, allocatable :: VOD(:)
        contains
        procedure, pass(self) :: initialize => VODOO_initialize_default
        procedure, pass(self) :: run => VODOO_run_default
    end type vodoo_type_dec

    type, extends(vodoo_type_dec) :: vodoo_linear_model
        integer, parameter :: ncoef = 8
        real, allocatable :: coefs(:, :)
        contains
            procedure, pass(self) :: initialize => VODOO_initialize_linear_model
            procedure, pass(self) :: run => VODOO_run_linear_model
    end type vodoo_linear_model

    class(vodoo_type_dec), allocatable :: vodoo_struc(:)

    SAVE

contains
    !BOP
    !
    ! !ROUTINE: VODOO_initialize
    ! \label{VODOO_initialize}
    !
    ! !INTERFACE:
    subroutine VODOO_initialize()
        ! !DESCRIPTION:
        !
        !  This routine creates the datatypes and allocates memory for noahMP3.6-specific
        !  variables. It also invokes the routine to read the runtime specific options
        !  for noahMP3.6 from the configuration file.
        !
        !  The routines invoked are:
        !  \begin{description}
        !   \item[readVODOOcrd](\ref{readVODOOcrd}) \newline
        !
        !EOP
        implicit none

        integer :: rc, ios
        integer :: n, nid, ngrid, ngridId
        character*100 :: modeltype(LIS_rc%nnest)
        character*10 :: predictors(LIS_rc%nnest)
        character*256 :: parameter_fname(LIS_rc%nnest)

        write(LIS_logunit,*) "[INFO] Starting VODOO setup"


        ! read config from file
        call ESMF_ConfigFindLabel(LIS_config, "VODOO parameter file:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, parameter_fname(n), rc=rc)
            call LIS_verify(rc, "VODOO parameter file: not defined")
        enddo

        ! remainder from when there where multiple models available, leaving it
        ! here to make it easier to extend it again
        if (.true.) then
            allocate(vodoo_linear_model :: vodoo_struc(LIS_rc%nnest))
        else
            write(LIS_logunit, *)&
                 "[ERR] VODOO model type must be 'linear' or 'SVR'"
            call LIS_endrun
        endif

        do n=1,LIS_rc%nnest
            vodoo_struc(n)%parameter_fname = parameter_fname(n)
            allocate(vodoo_struc(n)%VOD(LIS_rc%npatch(n,LIS_rc%lsm_index)))

            call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")
            call add_sfc_fields(n,LIS_sfcState(n), "Relative Soil Moisture Layer 1")
            call add_sfc_fields(n,LIS_sfcState(n), "Relative Root Zone Soil Moisture")
            call add_sfc_fields(n,LIS_forwardState(n),"VODOO_VOD")
        enddo


        ! read parameters/initialize model
        do n=1, LIS_rc%nnest
            call vodoo_struc(n)%initialize(n)
        enddo

        write(LIS_logunit,*) '[INFO] Finished VODOO setup'
    end subroutine VODOO_initialize

    subroutine VODOO_initialize_default(self, n)
        class(vodoo_type_dec), intent(inout) :: self
        integer, intent(in) :: n
        write(LIS_logunit,*) "[ERR] VODOO should use linear or svr model"
        call LIS_endrun
    end subroutine VODOO_initialize_default

    subroutine VODOO_initialize_linear_model(self, n)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        class(vodoo_linear_model), intent(inout) :: self
        integer, intent(in) :: n

        integer :: ios, nid, npatch
        integer :: ngridId, ngrid

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODOO requires NETCDF"
        call LIS_endrun
#else
        ! try opening the parameter file
        write(LIS_logunit,*) '[INFO] Reading ',&
            trim(self%parameter_fname)
        ios = nf90_open(path=trim(self%parameter_fname),&
            mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '&
            //trim(vodoo_struc(n)%parameter_fname))

        ! check if ngrid is as expected
        ios = nf90_inq_dimid(nid, "ngrid", ngridId)
        call LIS_verify(ios, "Error nf90_inq_varid: ngrid")
        ios = nf90_inquire_dimension(nid, ngridId, len=ngrid)
        call LIS_verify(ios, "Error nf90_inquire_dimension: ngrid")
        if (ngrid /= LIS_rc%glbngrid_red(n)) then
            write(LIS_logunit, *) "[ERR] ngrid in "//trim(self%parameter_fname)&
                 //" not consistent with expected ngrid: ", ngrid,&
                 " instead of ",LIS_rc%glbngrid_red(n)
            call LIS_endrun
        endif

        ! now that we have ntimes, we can allocate the coefficient arrays
        ! for the nest
        npatch = LIS_rc%npatch(n, LIS_rc%lsm_index)
        allocate(self%coefs(self%ncoef, npatch))

        call read_2d_coef_from_file(n, nid, self%ncoef, ngrid, "coefs", &
            self%coefs, ngrid_first=.false.)
#endif
    end subroutine VODOO_initialize_linear_model

    subroutine read_1d_coef_from_file(n, nid, ngrid, varname, coef)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, ngrid
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:)

        real, allocatable :: coef_file(:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODOO requires NETCDF"
        call LIS_endrun
#else

        allocate(coef_file(ngrid))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))
        allocate(coef_grid(LIS_rc%ngrid(n)))
        call LIS_convertVarToLocalSpace(n, coef_file, coef_grid)
        call gridvar_to_patchvar(&
             n, LIS_rc%lsm_index, coef_grid, coef)
        deallocate(coef_grid)
        deallocate(coef_file)
#endif
    end subroutine read_1d_coef_from_file

    subroutine read_2d_coef_from_file(n, nid, n1, n2, varname, coef, ngrid_first)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, n1, n2
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:,:)
        logical, intent(in), optional :: ngrid_first

        logical :: ngrid_is_first
        real, allocatable :: coef_file(:,:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODOO requires NETCDF"
        call LIS_endrun
#else
        if (present(ngrid_first)) then
            ngrid_is_first = ngrid_first
        else
            ngrid_is_first = .true.
        endif

        allocate(coef_file(n1, n2))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))

        allocate(coef_grid(LIS_rc%ngrid(n)))
        if (ngrid_is_first) then
            do j=1,n2
                call LIS_convertVarToLocalSpace(n, coef_file(:,j), coef_grid)
                call gridvar_to_patchvar(&
                     n, LIS_rc%lsm_index, coef_grid, coef(:, j))
            enddo
        else ! ngrid is the second dimension, we have to loop over n1
            do j=1,n1
                call LIS_convertVarToLocalSpace(n, coef_file(j, :), coef_grid)
                call gridvar_to_patchvar(&
                     n, LIS_rc%lsm_index, coef_grid, coef(j, :))
            enddo
        endif

        deallocate(coef_grid)
        deallocate(coef_file)
#endif
    end subroutine read_2d_coef_from_file

    subroutine read_3d_coef_from_file(n, nid, n1, n2, n3, varname, coef, &
            ngrid_first)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, n1, n2, n3
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:,:,:)
        logical, intent(in), optional :: ngrid_first

        logical :: ngrid_is_first
        real, allocatable :: coef_file(:,:,:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j2, j3

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODOO requires NETCDF"
        call LIS_endrun
#else
        if (present(ngrid_first)) then
            ngrid_is_first = ngrid_first
        else
            ngrid_is_first = .true.
        endif

        allocate(coef_file(n1, n2, n3))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))
        allocate(coef_grid(LIS_rc%ngrid(n)))
        if (ngrid_is_first) then
            do j2=1,n2
                do j3=1,n3
                    call LIS_convertVarToLocalSpace(n, coef_file(:,j2,j3), coef_grid)
                    call gridvar_to_patchvar(&
                         n, LIS_rc%lsm_index, coef_grid, coef(:, j2,j3))
                 enddo
            enddo
        else  ! ngrid is last
            do j2=1,n1
                do j3=1,n2
                    call LIS_convertVarToLocalSpace(n, coef_file(j2,j3, :), coef_grid)
                    call gridvar_to_patchvar(&
                         n, LIS_rc%lsm_index, coef_grid, coef(j2,j3, :))
                 enddo
            enddo
        endif
        deallocate(coef_grid)
        deallocate(coef_file)
#endif
    end subroutine read_3d_coef_from_file

    subroutine gridvar_to_patchvar(n,m,gvar,tvar)
        ! Converts a variable in local gridspace (length LIS_rc%ngrid(n))
        ! to local patch space (length LIS_rc%npatch(n,m))
        ! patch space = ensembles * ngrid

        implicit none

        integer, intent(in) :: n 
        integer, intent(in) :: m
        real, intent(in)    :: gvar(LIS_rc%ngrid(n))
        real, intent(inout) :: tvar(LIS_rc%npatch(n,m))
        integer             :: t,r,c

        do t=1,LIS_rc%npatch(n,m)
            r = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            c = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            if (LIS_domain(n)%gindex(c,r).ge.0) then
                tvar(t) = gvar(LIS_domain(n)%gindex(c,r))
            else
                tvar(t) = LIS_rc%udef
            endif
        enddo

    end subroutine gridvar_to_patchvar
    !!--------------------------------------------------------------------------------



    subroutine add_sfc_fields(n, sfcState,varname)

        implicit none

        integer            :: n
        type(ESMF_State)   :: sfcState
        character(len=*)   :: varname

        type(ESMF_Field)     :: varField
        type(ESMF_ArraySpec) :: arrspec
        integer              :: status
        real :: sum
        call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
             rc=status)
        call LIS_verify(status)

        varField = ESMF_FieldCreate(arrayspec=arrSpec, &
             grid=LIS_vecTile(n), name=trim(varname), &
             rc=status)
        call LIS_verify(status, 'Error in field_create of '//trim(varname))

        call ESMF_StateAdd(sfcState, (/varField/), rc=status)
        call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

    end subroutine add_sfc_fields


    subroutine VODOO_f2t(n)

        implicit none

        integer, intent(in)    :: n

    end subroutine VODOO_f2t


    subroutine VODOO_geometry(n)
        implicit none
        integer, intent(in)    :: n

    end subroutine VODOO_geometry

    subroutine VODOO_run(n)
        use LIS_histDataMod
        ! !USES:
        implicit none

        integer, intent(in) :: n

        integer             :: t
        integer             :: status
        integer             :: col,row
        real, pointer       :: lai(:), relsm1(:), relrzsm(:)
        real, pointer       :: vodval(:)
        real                :: smanom
        real                :: pred(vodoo_struc(n)%ncoef)

        !   map surface properties to SFC
        call getsfcvar(LIS_sfcState(n), "Leaf Area Index", lai)
        call getsfcvar(LIS_sfcState(n), "Relative Soil Moisture Layer 1", relsm1)
        call getsfcvar(LIS_sfcState(n), "Relative Root Zone Soil Moisture", relrzsm)

        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)

            smanom = relsm1(t) - max(0, min(1, relrzsm(t)))

            pred(1) = lai(t) * smanom
            pred(2) = rzsm(t) * smanom
            pred(3) = lai(t) * rzsm(t) * smanom
            pred(4) = smanom
            pred(5) = lai(t)
            pred(6) = rzsm(t)
            pred(7) = lai(t) * rzsm(t) 
            pred(8) = 1

            call vodoo_struc(n)%run(n, t, pred)

            if (vodoo_struc(n)%VOD(t).ne.LIS_rc%udef.and.vodoo_struc(n)%VOD(t).lt.-10) then
                write(LIS_logunit, *) "[WARN] VOD lower than -10"
                vodoo_struc(n)%VOD(t) = LIS_rc%udef
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=vodoo_struc(n)%VOD(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="-")
        enddo

        call getsfcvar(LIS_forwardState(n),"VODOO_VOD", vodval)
        vodval = vodoo_struc(n)%VOD


    end subroutine VODOO_run

    subroutine VODOO_run_default(self, n, t, pred)
        class(vodoo_type_dec), intent(inout) :: self
        integer, intent(in) :: n, t
        real, intent(in) :: pred(self%npred)
        write(LIS_logunit,*) "[ERR] VODOO should use linear or svr model"
        call LIS_endrun
    end subroutine VODOO_run_default

    subroutine VODOO_run_linear_model(self, n, t, pred)
        use LIS_histDataMod
        ! !USES:
        implicit none

        class(vodoo_linear_model), intent(inout) :: self
        integer, intent(in) :: n, t
        real, intent(in) :: pred(self%ncoef)

        real                :: coefs(self%ncoef)
        logical             :: coefs_valid

        coefs = self%coefs(:, t)

        ! For some pixels no VOD was available and therefore no forward
        ! model was fitted. No assimilation will take place over these
        ! pixels anyways, so it's no problem to not predict anything here
        coefs_valid = (.not.isnan(coefs(1)).and.coefs(1).ne.LIS_rc%udef)
        if (coefs_valid) then
            self%VOD(t) = sum(coefs * pred)
        else
            self%VOD(t)=LIS_rc%udef
        endif
    end subroutine VODOO_run_linear_model


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine getsfcvar(sfcState, varname, var)
        ! !USES:

        implicit none

        type(ESMF_State)      :: sfcState
        type(ESMF_Field)      :: varField
        character(len=*)      :: varname
        real, pointer         :: var(:)
        integer               :: status

        call ESMF_StateGet(sfcState, trim(varname), varField, rc=status)
        call LIS_verify(status, "Error in StateGet: VODOO_getsfcvar "//trim(varname))
        call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
        call LIS_verify(status, "Error in FieldGet: VODOO_getsfcvar "//trim(varname))

    end subroutine getsfcvar

!!!!BOP
!!!! !ROUTINE: VODOO_output
!!!! \label{VODOO_output}
!!!!
!!!! !INTERFACE:
    subroutine VODOO_output(n)
        integer, intent(in) :: n
    end subroutine VODOO_output
#endif
end module VODOO_Mod
