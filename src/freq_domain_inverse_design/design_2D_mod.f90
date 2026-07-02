module design_2D_mod

#ifdef USE_MPI
    use mpi
#endif

    use design_base_mod

    implicit none

    type, extends(TDesign) :: TDesign_2D

        real(dp) , allocatable :: ker_mat(:, :)
        real(dp) , allocatable :: rho(:, :)
        real(dp) , allocatable :: rho_conv(:, :)
        real(dp) , allocatable :: rho_old(:, :)
        real(dp) , allocatable :: grad(:, :)
        real(dp) , allocatable :: grad_old(:, :)
        real(dp) , allocatable :: grad_conv(:, :)
        logical  , allocatable :: opt_region(:, :)

        contains
            procedure :: init_design             => init_2D_design
            procedure :: kill_design             => kill_2D_design
            procedure :: collect_opt_regions     => collect_2D_opt_regions
            procedure :: set_opt_algo            => set_2D_opt_algo
            procedure :: collect_FOM             => collect_2D_FOM
            procedure :: collect_gradients       => collect_2D_gradients
            procedure :: apply_kernel_on_rho     => apply_kernel_on_rho_2D
            procedure :: apply_kernel_on_grad    => apply_kernel_on_grad_2D
            procedure :: opt_step                => opt_step_2D
            procedure :: reset_rho_one_step_back => reset_rho_one_step_back_2D
            procedure :: reset_grad              => reset_grad_2D
            procedure :: update_rho              => update_rho_2D
            procedure :: update_grad             => update_grad_2D

    end type TDesign_2D

contains
!###################################################################################################

subroutine init_2D_design(this, dimensions, dr, sigma, grid_Ndims, apply_grad_rho_init, delta_rho)

    class(TDesign_2D), intent(inout) :: this
    logical          , intent(in)     :: apply_grad_rho_init
    integer          , intent(in)     :: dimensions
    integer          , intent(in)     :: grid_Ndims(3)
    real(dp)         , intent(in)     :: dr
    real(dp)         , intent(in)     :: sigma
    real(dp)         , intent(in)     :: delta_rho

    integer       :: i, j
    integer       :: n_ker, nx, ny
    real(dp)      :: ker_sum
    real(dp)      :: x, y

    this%dimensions = dimensions
    this%nx         = grid_Ndims(1)
    this%ny         = grid_Ndims(2)
    this%apply_grad_rho_init = apply_grad_rho_init
    this%drho = delta_rho

    if (sigma <= 0.0_dp) then
        this%n_ker = 0
    else
        this%n_ker = int(3.0_dp*sigma/dr)
    end if

    nx    = this%nx
    ny    = this%ny
    n_ker = this%n_ker

    if (.not. allocated(this%rho))       allocate(this%rho(-n_ker+1:nx+n_ker, -n_ker+1:ny+n_ker))
    if (.not. allocated(this%rho_conv))  allocate(this%rho_conv(nx, ny))
    if (.not. allocated(this%rho_old))   allocate(this%rho_old(nx, ny))
    if (.not. allocated(this%grad))      allocate(this%grad(-n_ker+1:nx+n_ker, -n_ker+1:ny+n_ker))
    if (.not. allocated(this%grad_conv)) allocate(this%grad_conv(nx, ny))
    if (.not. allocated(this%grad_old))  allocate(this%grad_old(nx, ny))
    if (.not. allocated(this%opt_region)) allocate(this%opt_region(-n_ker+1:nx+n_ker, &
                                                                   -n_ker+1:ny+n_ker))
    if (.not. allocated(this%ker_mat))   allocate(this%ker_mat(-n_ker:n_ker, -n_ker:n_ker))

    if (sigma <= 0.0_dp) then
        this%ker_mat(0, 0) = 1.0_dp
    else
        ker_sum = 0.0_dp
        do j = -n_ker, n_ker
        do i = -n_ker, n_ker
            x = real(i, dp) * dr
            y = real(j, dp) * dr
            this%ker_mat(i, j) = EXP(-(x**2 + y**2)/(2.0_dp*sigma**2))
            ker_sum = ker_sum + this%ker_mat(i, j)
        end do
        end do
        this%ker_mat = this%ker_mat / ker_sum
    end if


    this%rho        = 0.0_dp
    this%rho_conv   = 0.0_dp
    this%rho_old    = 0.0_dp
    this%opt_region = .false.
    this%grad_max   = 0.0_dp

    this%continue_opt  = .true.
    this%new_rho_set   = .true.
    this%change_beta   = .false.

end subroutine init_2D_design

!###################################################################################################

subroutine kill_2D_design(this)

    class(TDesign_2D), intent(inout) :: this

    if (allocated(this%rho))       deallocate(this%rho)
    if (allocated(this%rho_conv))  deallocate(this%rho_conv)
    if (allocated(this%rho_old))   deallocate(this%rho_old)
    if (allocated(this%grad))      deallocate(this%grad)
    if (allocated(this%grad_conv)) deallocate(this%grad_conv)
    if (allocated(this%grad_old))  deallocate(this%grad_old)
    if (allocated(this%opt_region)) deallocate(this%opt_region)
    if (allocated(this%ker_mat))   deallocate(this%ker_mat)

    call this%opt_algo%kill_opt_algo()

end subroutine kill_2D_design

!###################################################################################################

subroutine collect_2D_opt_regions(this, opt_region_i, rho_init)

    class(TDesign_2D), intent(inout) :: this
    logical          , intent(in)    :: opt_region_i(:,:,:)
    real(dp)         , intent(in)    :: rho_init

    integer :: i, j
    real(dp) :: x

    call random_seed() 

    do j = 1, this%ny
    do i = 1, this%nx
        if (opt_region_i(i, j, 1)) then
            this%opt_region(i, j) = .true.
            call random_number(x)
            this%rho(i, j) = rho_init + 0.001_dp * x !Random perturbation
        else
            this%opt_region(i, j) = .false.
            this%rho(i, j) = 0.0_dp
        end if
    end do
    end do

end subroutine collect_2D_opt_regions
!###################################################################################################

subroutine set_2D_opt_algo(this, m_opt, iprint, factr, pgtol)

    class(TDesign_2D), intent(inout) :: this
    integer          , intent(in)    :: m_opt
    integer          , intent(in)    :: iprint
    real(dp)         , intent(in)    :: factr
    real(dp)         , intent(in)    :: pgtol

    integer :: i, j
    integer :: n_opt
    integer :: nx
    integer :: ny

    nx = this%nx
    ny = this%ny

    n_opt=0
    do j = 1, ny
    do i = 1, nx
        if (this%opt_region(i,j)) then
            n_opt = n_opt + 1
        end if
    end do
    end do

    call this%opt_algo%init_opt_algo(n_opt, m_opt, iprint, factr, pgtol, this%sigma, &
                                     this%dimensions)

end subroutine set_2D_opt_algo

!###################################################################################################

subroutine collect_2D_FOM(this, w_p, p, n_opt_problems)

    class(TDesign_2D), intent(inout) :: this
    real(dp)         , intent(in)    :: w_p
    integer          , intent(in)    :: p
    integer          , intent(in)    :: n_opt_problems

    integer  :: ierr
    real(dp) :: fom_loc = 0.0_dp
    real(dp) :: fom_sum = 0.0_dp

    if (p == 1) this%fom = 0.0_dp

    this%fom = this%fom + w_p

    !MPI is already used to share w_p across ranks.

end subroutine collect_2D_FOM

!###################################################################################################

subroutine collect_2D_gradients(this, grad_in, p, n_opt_problems)

    class(TDesign_2D), intent(inout) :: this
    real(dp)         , intent(in)    :: grad_in(:,:,:)
    integer          , intent(in)    :: p
    integer          , intent(in)    :: n_opt_problems

    if  (p == 1) this%grad = 0.0_dp

    this%grad(1:this%nx,1:this%ny) = this%grad(1:this%nx,1:this%ny) + &
                                           grad_in(1:this%nx,1:this%ny,1)

    if (p == n_opt_problems) then
        this%grad = 2.0_dp * this%fom * this%grad
    end if

end subroutine collect_2D_gradients

!###################################################################################################

subroutine apply_kernel_on_rho_2D(this)

    class(TDesign_2D), intent(inout) :: this

    integer :: i, j, ii, jj

    this%rho_conv = 0.0_dp

    do j = 1, this%ny
    do i = 1, this%nx
        if (this%opt_region(i,j)) then
            do ii = -this%n_ker, this%n_ker
            do jj = -this%n_ker, this%n_ker
                if (this%opt_region(i+ii,j+jj)) then
                    this%rho_conv(i,j) = this%rho_conv(i,j) + &
                        this%ker_mat(ii,jj)*this%rho(i+ii,j+jj)
                end if
            end do
            end do
        end if
    end do
    end do

end subroutine apply_kernel_on_rho_2D

!###################################################################################################

subroutine apply_kernel_on_grad_2D(this)

    class(TDesign_2D), intent(inout) :: this

    integer :: i, j, ii, jj

    this%grad_conv = 0.0_dp

    do i = 1, this%nx
    do j = 1, this%ny
        if (this%opt_region(i,j)) then
            do ii = -this%n_ker, this%n_ker
            do jj = -this%n_ker, this%n_ker
                if (this%opt_region(i+ii,j+jj)) then
                    this%grad_conv(i,j) = this%grad_conv(i,j) + &
                        this%ker_mat(ii,jj)*this%grad(i+ii,j+jj)
                end if
            end do
            end do
        end if
    end do
    end do

end subroutine apply_kernel_on_grad_2D

!###################################################################################################

subroutine opt_step_2D(this)

    class(TDesign_2D), intent(inout) :: this

    integer  :: nx, ny
    integer  :: i, j
    real(dp) :: norm_loc
    real(dp) :: norm_global

    nx = this%nx
    ny = this%ny

    if (this%first_iter .and. this%apply_grad_rho_init) then

        norm_loc = MAXVAL(ABS(this%grad_conv))

#ifdef USE_MPI
        call MPI_Allreduce(norm_loc, norm_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                        MPI_COMM_WORLD, ierr)
#else
        norm_global = norm_loc
#endif

        if (norm_global <= 0.0_dp) norm_global = 1.0_dp

        this%grad_max = norm_global

        this%rho(1:this%nx,1:this%ny) = this%rho_old(1:this%nx,1:this%ny) - &
                                        this%grad_conv(1:this%nx,1:this%ny) * this%drho / norm_global

        do j = 1, this%ny
        do i = 1, this%nx
            if (this%rho(i,j) < 0.0_dp) this%rho(i,j) = 0.0_dp
        end do
        end do

        this%first_iter = .false.

    else

        if (this%opt_algo%task(1:2) == 'FG' .or. this%opt_algo%task(1:5) == 'START') then
            this%opt_algo%f = -this%fom**2
        end if

        call this%opt_algo%lbfgsb_optimize_2D(this%rho(1:nx,1:ny), this%grad_conv(1:nx,1:ny), &
                                            nx, ny, this%opt_region(1:nx,1:ny))

        if (this%opt_algo%task(1:5) == 'NEW_X') this%new_rho_set = .false.
        if (this%opt_algo%task(1:2) == 'FG')    this%new_rho_set = .true.
        this%change_beta  = .false.

        if (.not.(this%opt_algo%task(1:2) == 'FG' .or. this%opt_algo%task(1:5) == 'NEW_X' .or. &
                this%opt_algo%task(1:5) == 'START')) then

            this%change_beta  = .true.
            this%opt_algo%task = 'START'
            this%new_rho_set  = .true.

        end if
    end if

end subroutine opt_step_2D

!###################################################################################################

subroutine reset_rho_one_step_back_2D(this)

    class(TDesign_2D), intent(inout) :: this

    this%rho(1:this%nx,1:this%ny) = this%rho_old(1:this%nx,1:this%ny)

end subroutine reset_rho_one_step_back_2D

!###################################################################################################

subroutine reset_grad_2D(this)

    class(TDesign_2D), intent(inout) :: this

    this%grad(1:this%nx,1:this%ny) = 0.0_dp

end subroutine reset_grad_2D

!###################################################################################################

subroutine update_rho_2D(this)

    class(TDesign_2D), intent(inout) :: this

    this%rho_old(1:this%nx,1:this%ny) = this%rho(1:this%nx,1:this%ny)

end subroutine update_rho_2D

!###################################################################################################

subroutine update_grad_2D(this)

    class(TDesign_2D), intent(inout) :: this

    this%grad_old(1:this%nx,1:this%ny) = this%grad_conv(1:this%nx,1:this%ny)

end subroutine update_grad_2D

!###################################################################################################

end module design_2D_mod