! Routines & modules to calculate temporal averages
!
! Written by Samuel Gotterbarm, 23.08.2025 (assisted by ChatGPT)
module LIS_timeAvgMod
    implicit none

    ! public types
    public :: LIS_TemporalAverage

    type :: LIS_TemporalAverage
        private
        integer :: n = 0                 ! number of values to keep
        real, allocatable :: buffer(:)   ! circular buffer
        integer :: idx = 0               ! current position
        integer :: count = 0             ! how many values stored so far
        real :: total = 0.0              ! running sum
    contains
        procedure :: init_tmp_avg_from_size
        procedure :: init_tmp_avg_from_timestep
        generic   :: init => init_tmp_avg_from_size, init_tmp_avg_from_timestep
        procedure :: add_value    => add_value_for_tmp_avg
        procedure :: get          => get_value_of_tmp_avg
    end type LIS_TemporalAverage

contains
    
    subroutine init_tmp_avg_from_size(this, n)
        class(LIS_TemporalAverage), intent(inout) :: this
        integer, intent(in) :: n

        write(LIS_logunit,*) "[DEBUG] Temporal average initalised with buffer size ", n
        if (allocated(this%buffer)) deallocate(this%buffer)
        allocate(this%buffer(n))
        this%buffer = 0.0
        this%n = n
        this%idx = 0
        this%count = 0
        this%total = 0.0
    end subroutine init_tmp_avg_from_size

    subroutine init_tmp_avg_from_timestep(this, timestep, totaltime)
        class(LIS_TemporalAverage), intent(inout) :: this
        real, intent(in) :: timestep
        real, intent(in) :: totaltime
        integer :: n
        n = max(1, int(totaltime/timestep))   ! ensure at least 1 sample
        call this%init_tmp_avg_from_size(n)
    end subroutine init_tmp_avg_from_timestep

    subroutine add_value_for_tmp_avg(this, value)
        class(LIS_TemporalAverage), intent(inout) :: this
        real, intent(in) :: value

        ! Advance circular index
        this%idx = mod(this%idx, this%n) + 1

        ! Subtract old value if buffer is full
        if (this%count >= this%n) then
            this%total = this%total - this%buffer(this%idx)
        else
            this%count = this%count + 1
        end if

        ! Store new value
        this%buffer(this%idx) = value
        this%total = this%total + value
    end subroutine add_value_for_tmp_avg

    function get_value_of_tmp_avg(this) result(avgval)
        class(LIS_TemporalAverage), intent(in) :: this
        real :: avgval

        if (this%count .gt. 0) then
            avgval = this%total / real(this%count) 
        else
            avgval = 0.0
        end if
    end function get_value_of_tmp_avg

end module
