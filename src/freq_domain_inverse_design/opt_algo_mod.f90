module opt_algo_mod

    use constants_mod
    use lbfgsb_all_mod

    implicit none

    type TOptAlgo

        character(len=60) :: task
        character(len=60) :: csave
        integer           :: n
        integer           :: m
        integer           :: iprint
        real(dp)          :: factr
        real(dp)          :: pgtol
        real(dp)          :: f
        
        integer , allocatable :: nbd(:)
        integer , allocatable :: iwa(:)
        integer , allocatable :: isave(:)
        logical , allocatable :: lsave(:)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: l(:)
        real(dp), allocatable :: u(:)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: wa(:)
        real(dp), allocatable :: dsave(:)

        contains
            procedure :: init_opt_algo
            procedure :: kill_opt_algo
            procedure :: lbfgsb_optimize_1D
            procedure :: lbfgsb_optimize_2D
            procedure :: lbfgsb_optimize_3D

    end type TOptAlgo

contains
!###################################################################################################
subroutine init_opt_algo(this, n_opt, m_opt, iprint, factr, pgtol, sigma, dim)
    
    class(TOptAlgo), intent(inout) :: this
    integer        , intent(in)    :: n_opt
    integer        , intent(in)    :: m_opt
    integer        , intent(in)    :: iprint
    integer        , intent(in)    :: dim
    real(dp)       , intent(in)    :: factr
    real(dp)       , intent(in)    :: pgtol
    real(dp)       , intent(in)    :: sigma

    this%n = n_opt
    this%m = m_opt
    this%iprint = iprint
    this%factr = factr
    this%pgtol = pgtol

    if (.not. allocated(this%nbd))   allocate(this%nbd(n_opt))
    if (.not. allocated(this%iwa))   allocate(this%iwa(3*n_opt))
    if (.not. allocated(this%isave)) allocate(this%isave(44))
    if (.not. allocated(this%lsave)) allocate(this%lsave(4))
    if (.not. allocated(this%x))     allocate(this%x(n_opt))
    if (.not. allocated(this%l))     allocate(this%l(n_opt))
    if (.not. allocated(this%u))     allocate(this%u(n_opt))
    if (.not. allocated(this%g))     allocate(this%g(n_opt))
    if (.not. allocated(this%wa))    allocate(this%wa(2*m_opt*n_opt + 5*n_opt + 11*m_opt*m_opt + 8*m_opt))
    if (.not. allocated(this%dsave)) allocate(this%dsave(29))

    this%l     = 0.0_dp
    this%u     = 1.0_dp * DSQRT(2 * pi * sigma**2)**dim
    this%nbd   = 2
    this%task  = 'START'
    this%f     = 0.0_dp

end subroutine init_opt_algo

!###################################################################################################

subroutine kill_opt_algo(this)
    class(TOptAlgo), intent(inout) :: this

    if (allocated(this%nbd))   deallocate(this%nbd)
    if (allocated(this%iwa))   deallocate(this%iwa)
    if (allocated(this%isave)) deallocate(this%isave)
    if (allocated(this%lsave)) deallocate(this%lsave)
    if (allocated(this%x))     deallocate(this%x)
    if (allocated(this%l))     deallocate(this%l)
    if (allocated(this%u))     deallocate(this%u)
    if (allocated(this%g))     deallocate(this%g)
    if (allocated(this%wa))    deallocate(this%wa)
    if (allocated(this%dsave)) deallocate(this%dsave)

end subroutine kill_opt_algo

!###################################################################################################

subroutine lbfgsb_optimize_1D(this, rho, grad, nx, opt_region)

    class(TOptAlgo), intent(inout) :: this
    real(dp)       , intent(inout) :: rho(nx)
    real(dp)       , intent(inout) :: grad(nx)
    integer        , intent(in)    :: nx
    logical        , intent(in)    :: opt_region(nx)

    integer :: i
    integer :: ii

    if (this%task(1:2) == 'FG' .or. this%task(1:5) == 'START') then
        ii = 1
        do i = 1, nx
            if (opt_region(i)) then
                this%x(ii) = rho(i)
                this%g(ii) = grad(i)
                ii = ii + 1
            end if
        end do
    end if

    call setulb(this%n, this%m, this%x, this%l, this%u, this%nbd, this%f, this%g, this%factr, &
                this%pgtol, this%wa, this%iwa, this%task, this%iprint, this%csave,       &
                this%lsave, this%isave, this%dsave)

    if (this%task(1:2) == 'FG') then
        ii = 1
        do i = 1, nx
            if (opt_region(i)) then
                rho(i) = this%x(ii)
                ii = ii + 1
            end if
        end do
    end if 

end subroutine lbfgsb_optimize_1D

!###################################################################################################

subroutine lbfgsb_optimize_2D(this, rho, grad, nx, ny, opt_region)

    class(TOptAlgo), intent(inout) :: this
    real(dp)       , intent(inout) :: rho(nx, ny)
    real(dp)       , intent(inout) :: grad(nx, ny)
    integer        , intent(in)    :: nx, ny
    logical        , intent(in)    :: opt_region(nx, ny)

    integer :: i, j
    integer :: ii

    if (this%task(1:2) == 'FG' .or. this%task(1:5) == 'START') then
        ii = 1
        do j = 1, ny
        do i = 1, nx
            if (opt_region(i, j)) then
                this%x(ii) = rho(i, j)
                this%g(ii) = grad(i, j)
                ii = ii + 1
            end if
        end do
        end do
    end if

    call setulb(this%n, this%m, this%x, this%l, this%u, this%nbd, this%f, this%g, this%factr, &
                this%pgtol, this%wa, this%iwa, this%task, this%iprint, this%csave,       &
                this%lsave, this%isave, this%dsave)

    if (this%task(1:2) == 'FG') then
        ii = 1
        do j = 1, ny
        do i = 1, nx
            if (opt_region(i, j)) then
                rho(i, j) = this%x(ii)
                ii = ii + 1
            end if
        end do
        end do
    end if

end subroutine lbfgsb_optimize_2D

!###################################################################################################

subroutine lbfgsb_optimize_3D(this, rho, grad, nx, ny, nz, opt_region)
    class(TOptAlgo), intent(inout) :: this

    real(dp)       , intent(inout) :: rho(nx, ny, nz)
    real(dp)       , intent(inout) :: grad(nx, ny, nz)
    integer        , intent(in)    :: nx, ny, nz
    logical        , intent(in)    :: opt_region(nx, ny, nz)

    integer :: i, j, k
    integer :: ii


    if (this%task(1:2) == 'FG' .or. this%task(1:5) == 'START') then
        ii = 1
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            if (opt_region(i, j, k)) then
                this%x(ii) = rho(i, j, k)
                this%g(ii) = grad(i, j, k)
                ii = ii + 1
            end if
        end do
        end do
        end do
    end if

    call setulb(this%n, this%m, this%x, this%l, this%u, this%nbd, this%f, this%g, this%factr, &
                this%pgtol, this%wa, this%iwa, this%task, this%iprint, this%csave,       &
                this%lsave, this%isave, this%dsave)

    if (this%task(1:2) == 'FG') then

        ii = 1
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            if (opt_region(i, j, k)) then
                rho(i, j, k) = this%x(ii)
                ii = ii + 1
            end if
        end do
        end do
        end do

    end if

end subroutine lbfgsb_optimize_3D

!###################################################################################################

end module opt_algo_mod