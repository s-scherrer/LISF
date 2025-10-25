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
        integer :: zerocount = 0         ! number of stored values that are zero
        real :: total = 0.0              ! running sum
    contains
        procedure :: init_tmp_avg_from_size
        procedure :: init_tmp_avg_from_timestep
        generic   :: init => init_tmp_avg_from_size, init_tmp_avg_from_timestep
        procedure :: add_value    => add_value_for_tmp_avg
        procedure :: get          => get_value_of_tmp_avg
        procedure :: print_state  => print_state_of_tmp_avg
    end type LIS_TemporalAverage

contains
    
    subroutine init_tmp_avg_from_size(this, n)
        use LIS_logmod, only: LIS_logunit

        class(LIS_TemporalAverage), intent(inout) :: this
        integer, intent(in) :: n

        ! write(LIS_logunit,*) "[DEBUG] Temporal average initalised with buffer size ", n
        if (allocated(this%buffer)) deallocate(this%buffer)
        allocate(this%buffer(n))
        this%buffer = 0.0
        this%n = n
        this%idx = 0
        this%count = 0
        this%zerocount = 0
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

        real :: oldval

        ! Advance circular index
        this%idx = mod(this%idx, this%n) + 1

        ! Subtract old value if buffer is full
        if (this%count >= this%n) then
            oldval = this%buffer(this%idx)
            if (oldval == 0.0) then
                this%zerocount = this%zerocount - 1
            endif
            this%total = this%total - this%buffer(this%idx)
        else
            this%count = this%count + 1
        end if

        ! Store new value
        this%buffer(this%idx) = value
        this%total = this%total + value
        if (value == 0.0) then
            this%zerocount = this%zerocount + 1
        endif

        ! Calculating the total by subtracting the last
        ! value in the buffer and adding the new value
        ! has the risk of accumulating numerical roundoff
        ! error. I encountered the case where all values
        ! were zero, but the total slightly nonzero. To
        ! capture these cases, an additional zero-value
        ! counter is used, and if all values are zero,
        ! the total sum is set to zero.
        if (this%zerocount == this%count) then
          this%total = 0.0
        endif

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

    subroutine print_state_of_tmp_avg(this)
        use LIS_logmod, only: LIS_logunit
        class(LIS_TemporalAverage), intent(in) :: this

        integer :: i

        write(LIS_logunit,*) "[DEBUG] Size: ", this%n
        write(LIS_logunit,*) "[DEBUG] Count: ", this%count
        write(LIS_logunit,*) "[DEBUG] Total: ", this%total
        write(LIS_logunit,*) "[DEBUG] Idx: ", this%idx
        do i=1, this%n
          write(LIS_logunit,*) "[DEBUG] Buffer: ", i, this%buffer(i)
        enddo
    end subroutine print_state_of_tmp_avg

end module
