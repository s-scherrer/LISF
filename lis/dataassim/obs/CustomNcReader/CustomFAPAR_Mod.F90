! This defines the CustomFAPAR reader, based on the CustomNcReader
! It is basically just an instance of CustomNcReader_dec, and the associated
! subroutines call the respective ones of CustomNcReader_Mod.
! Normally it should be enough to copy, replace FAPAR by your new variable name,
! and potentially adapt the settings in CustomFAPAR_setup.
module CustomFAPAR_Mod
    use CustomNcReader_Mod, only: CustomNcReader_dec

    implicit none

    public :: CustomInstFAPAR_setup, read_CustomFAPAR, write_CustomFAPAR
    public :: CustomFAPAR_struc

    ! declare public reader array
    type(CustomNcReader_dec), allocatable :: CustomFAPAR_struc(:)

contains

    subroutine CustomInstFAPAR_setup(k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: CustomNcReader_setup

        implicit none

        ! !ARGUMENTS:
        integer                   :: k
        type(ESMF_State)          :: OBS_State(LIS_rc%nnest)
        type(ESMF_State)          :: OBS_Pert_State(LIS_rc%nnest)

        integer :: n

        allocate(CustomFAPAR_struc(LIS_rc%nnest))
        do n=1,LIS_rc%nnest
            CustomFAPAR_struc(n)%obsid = "Custom InstFAPAR"
            CustomFAPAR_struc(n)%varname = "FAPAR"
            CustomFAPAR_struc(n)%min_value = 0.0001
            CustomFAPAR_struc(n)%max_value = 10.0
            CustomFAPAR_struc(n)%qcmin_value = 0.0
            CustomFAPAR_struc(n)%qcmax_value = 100.0
        enddo

        call CustomNcReader_setup(CustomFAPAR_struc, k, OBS_State, OBS_Pert_State)
    end subroutine CustomInstFAPAR_setup

    subroutine CustomInstBsFAPAR_setup(k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: CustomNcReader_setup

        implicit none

        ! !ARGUMENTS:
        integer                   :: k
        type(ESMF_State)          :: OBS_State(LIS_rc%nnest)
        type(ESMF_State)          :: OBS_Pert_State(LIS_rc%nnest)

        integer :: n

        allocate(CustomFAPAR_struc(LIS_rc%nnest))
        do n=1,LIS_rc%nnest
            CustomFAPAR_struc(n)%obsid = "Custom InstBlackSkyFAPAR"
            CustomFAPAR_struc(n)%varname = "FAPAR"
            CustomFAPAR_struc(n)%min_value = 0.0001
            CustomFAPAR_struc(n)%max_value = 10.0
            CustomFAPAR_struc(n)%qcmin_value = 0.0
            CustomFAPAR_struc(n)%qcmax_value = 100.0
        enddo

        call CustomNcReader_setup(CustomFAPAR_struc, k, OBS_State, OBS_Pert_State)
    end subroutine CustomInstBsFAPAR_setup


    subroutine CustomInstWsFAPAR_setup(k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: CustomNcReader_setup

        implicit none

        ! !ARGUMENTS:
        integer                   :: k
        type(ESMF_State)          :: OBS_State(LIS_rc%nnest)
        type(ESMF_State)          :: OBS_Pert_State(LIS_rc%nnest)

        integer :: n

        allocate(CustomFAPAR_struc(LIS_rc%nnest))
        do n=1,LIS_rc%nnest
            CustomFAPAR_struc(n)%obsid = "Custom InstWhiteSkyFAPAR"
            CustomFAPAR_struc(n)%varname = "FAPAR"
            CustomFAPAR_struc(n)%min_value = 0.0001
            CustomFAPAR_struc(n)%max_value = 10.0
            CustomFAPAR_struc(n)%qcmin_value = 0.0
            CustomFAPAR_struc(n)%qcmax_value = 100.0
        enddo

        call CustomNcReader_setup(CustomFAPAR_struc, k, OBS_State, OBS_Pert_State)
    end subroutine CustomInstWsFAPAR_setup


    subroutine read_CustomFAPAR(n, k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: read_CustomNetCDF
        integer, intent(in)       :: n, k
        type(ESMF_State)          :: OBS_State
        type(ESMF_State)          :: OBS_Pert_State
        call read_CustomNetCDF(CustomFAPAR_struc, n, k, OBS_State, OBS_Pert_State)
    end subroutine read_CustomFAPAR

    subroutine write_CustomFAPAR(n, k, OBS_State)
        use ESMF
        use CustomNcReader_Mod, only: write_CustomNetCDF
        integer,     intent(in)  :: n, k
        type(ESMF_State)         :: OBS_State
        call write_CustomNetCDF(n, k, OBS_State)
    end subroutine write_CustomFAPAR


end module CustomFAPAR_Mod
