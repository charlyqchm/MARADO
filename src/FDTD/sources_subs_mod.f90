module sources_subs_mod

    use constants_mod
    use mxll_base_mod
    use mxll_1D_mod

    implicit none

    type :: TSources_list

        integer :: n_p_src
        integer :: n_pw_src
        integer :: n_gb_src
        integer :: n_srcs

        type(TPointSrc)    , allocatable :: points(:)
        type(TPlaneWaveSrc), allocatable :: plane_waves(:)
        type(TGaussbeamSrc), allocatable :: gauss_beams(:)

    contains
        procedure :: read_init_sources
        procedure :: propagate_pw_srcs
        procedure :: propagate_p_srcs
        procedure :: kill_sources

    end type TSources_list

    type :: TPointSrc

        character(len=2)    :: polarization
        integer             :: dim
        integer             :: n_ker
        real(dp)            :: j_amp
        real(dp)            :: w0
        real(dp)            :: t0
        real(dp)            :: tau
        real(dp)            :: t_init
        real(dp)            :: t_final
        real(dp)            :: phase
        real(dp)            :: r0(3)
        real(dp)            :: rad

        real(dp), allocatable  :: ker_mat(:,:,:)
        real(dp), allocatable  :: J_mat(:,:,:)

        integer , allocatable  :: ind_i(:, :, :)
        integer , allocatable  :: ind_j(:, :, :)
        integer , allocatable  :: ind_k(:, :, :)
        logical , allocatable  :: in_this_rank(:,:,:)

    contains

        procedure :: init_point_src
        procedure :: kill_point_src
        procedure :: propagate_point_src

    end type TPointSrc

    type :: TPlaneWaveSrc

        type(TMxll_1D)      :: mxll_inc
        integer             :: dim
        !Logical to indicate whether the limits of the plane wave source in each direction are 
        !delimited by the user or they extend to the whole box.
        !The component 1,2 or 3 become .false. if the option "plane_wave_kx", "plane_wave_ky" &
        !or "plane_wave_kz" is selected, respectively.
        logical             :: limited_axis(3) = .true. 
        real(dp)            :: E_amp
        !Assuming the k vector is initially parallel to z, the E field to x and the H field to y,
        !phi is the angle respect to x, theta is the angle respect to z and psi is the angle
        !of rotation of the E field respect to the plane defined by k and z.  
        real(dp)            :: phi
        real(dp)            :: theta
        real(dp)            :: psi
        real(dp)            :: v_vec(3)
        real(dp)            :: A_vec(3)
        real(dp)            :: w0              
        real(dp)            :: t0
        real(dp)            :: tau
        real(dp)            :: t_init
        real(dp)            :: t_final
        real(dp)            :: phase
        integer             :: i_min, i_max
        integer             :: j_min, j_max
        integer             :: k_min, k_max
        integer             :: i_min_loc, i_max_loc
        integer             :: j_min_loc, j_max_loc
        integer             :: k_min_loc, k_max_loc
        logical             :: i_min_in_this_rank, i_max_in_this_rank
        logical             :: j_min_in_this_rank, j_max_in_this_rank
        logical             :: k_min_in_this_rank, k_max_in_this_rank

    contains

        procedure :: init_plane_wave_src
        procedure :: kill_plane_wave_src
        procedure :: propagate_plane_wave_src

    end type TPlaneWaveSrc

    type :: TGaussbeamSrc

        complex(dp)         :: E_rzt
        real(dp)            :: E_amp
        integer             :: dim
        !Assuming the k vector is initially parallel to z, the E field to x and the H field to y,
        !phi is the angle respect to x, theta is the angle respect to z and psi is the angle
        !of rotation of the E field respect to the plane defined by k and z.  
        real(dp)            :: phi
        real(dp)            :: theta
        real(dp)            :: psi
    
        real(dp)            :: v_vec(3)
        real(dp)            :: freq
        real(dp)            :: t0
        real(dp)            :: dz_ramp
        real(dp)            :: z0_ramp
        real(dp)            :: phase
        real(dp)            :: w0
        real(dp)            :: w
        real(dp)            :: z_R
        real(dp)            :: r0(3) !focus position
        real(dp)            :: lambda
        real(dp)            :: k
        real(dp)            :: r_min, r_max
        integer             :: i_min, i_max
        integer             :: j_min, j_max
        integer             :: k_min, k_max
        integer             :: i_min_loc, i_max_loc
        integer             :: j_min_loc, j_max_loc
        integer             :: k_min_loc, k_max_loc
        logical             :: i_min_in_this_rank, i_max_in_this_rank
        logical             :: j_min_in_this_rank, j_max_in_this_rank
        logical             :: k_min_in_this_rank, k_max_in_this_rank

    contains

        procedure :: init_gaussbeam_src
        procedure :: kill_gaussbeam_src
        procedure :: compute_time_space_profile

    end type TGaussbeamSrc

contains

!###################################################################################################

subroutine read_init_sources(this, dimensions, dt, dr, grid_Ndims, mpi_coords, mpi_dims)
    
    class(TSources_list), intent(inout) :: this
    integer             , intent(in)  :: dimensions
    integer             , intent(in)  :: grid_Ndims(3)
    integer             , intent(in)  :: mpi_coords(3)
    integer             , intent(in)  :: mpi_dims(3)
    real(dp)            , intent(in)  :: dr
    real(dp)            , intent(in)  :: dt

    character(len = 1000) :: input_ch
    character(len = 15)   :: src_type_ch

    character(len=20) :: src_file = "sources.in"
    logical           :: pw_source
    logical           :: p_source
    logical           :: gb_source
    integer           :: ierr, funit
    integer           :: i
    integer           :: n_p_src, n_pw_src, n_gb_src

    ! Check whether file exists.
    inquire (file=src_file, iostat=ierr)
    if (ierr /= 0) then
        write (*, '("Warning: There is no source input file ", A)') trim(src_file)
        return
    end if

    ! Open and read source input file.
    open (action='read', file=src_file, iostat=ierr, newunit=funit)
    if (ierr /= 0) then
        write (*, '("Error: could not open source input file ", A)') trim(src_file)
        error stop
    end if

    n_p_src = 0
    n_pw_src = 0
    n_gb_src = 0

    do

        pw_source = .false.
        p_source  = .false.
        gb_source = .false.

        read (unit=funit, fmt=*, iostat=ierr) src_type_ch

        if (ierr /= 0) exit

        p_source  = (trim(src_type_ch) == "point")
        pw_source = any(trim(src_type_ch) == [character(len=15) :: "plane_wave", &
                    "plane_wave_kx", "plane_wave_ky", "plane_wave_kz"])
        gb_source = (trim(src_type_ch) == "gaussbeam")

        if (p_source) then
            n_p_src = n_p_src + 1
        else if (pw_source) then
            n_pw_src = n_pw_src + 1
        else if (gb_source) then
            n_gb_src = n_gb_src + 1
        else
            write (*, '("Error: unknown source type ", A)') trim(src_type_ch)
            error stop
        end if

    end do

    rewind(funit)

    this%n_p_src = n_p_src
    this%n_pw_src = n_pw_src
    this%n_gb_src = n_gb_src

    if (n_p_src > 0) then
        if (.not. allocated(this%points)) allocate(this%points(this%n_p_src))
    end if

    if (n_pw_src > 0) then
        if (.not. allocated(this%plane_waves)) &
            allocate(this%plane_waves(this%n_pw_src))
    end if

    if (n_gb_src > 0) then
        if (.not. allocated(this%gauss_beams)) &
            allocate(this%gauss_beams(this%n_gb_src))
    end if

    n_p_src  = 1
    n_pw_src = 1
    n_gb_src = 1

    do

        p_source  = .false.
        pw_source = .false.
        gb_source = .false.

        read (unit=funit, fmt='(A)', iostat=ierr) input_ch

        if (ierr /= 0) exit

        read (input_ch, *) src_type_ch

        p_source  = (trim(src_type_ch) == "point")
        pw_source = any(trim(src_type_ch) == [character(len=15) :: "plane_wave", &
                    "plane_wave_kx", "plane_wave_ky", "plane_wave_kz"])
        gb_source = (trim(src_type_ch) == "gaussbeam")

        if (p_source) then
            call this%points(n_p_src)%init_point_src(input_ch, dimensions, dt, dr, &
                                                            grid_Ndims, mpi_coords, mpi_dims)
            n_p_src = n_p_src + 1
        else if (pw_source) then
            call this%plane_waves(n_pw_src)%init_plane_wave_src(input_ch, dimensions, &
                                                        dt , dr, grid_Ndims, mpi_coords, mpi_dims)
            n_pw_src = n_pw_src + 1
        else if (gb_source) then
            call this%gauss_beams(n_gb_src)%init_gaussbeam_src(input_ch, dimensions, dt, dr, &
                                                        grid_Ndims, mpi_coords, mpi_dims)
            n_gb_src = n_gb_src + 1
        else
            write (*, '("Error: unknown source type ", A)') trim(src_type_ch)
            error stop
        end if

    end do

    close(funit)

end subroutine read_init_sources

!###################################################################################################

subroutine kill_sources(this)

    class(TSources_list), intent(inout) :: this
    integer             :: i

    do i = 1, this%n_p_src
        call this%points(i)%kill_point_src()
    end do

    do i = 1, this%n_pw_src
        call this%plane_waves(i)%kill_plane_wave_src()
    end do

    do i = 1, this%n_gb_src
        call this%gauss_beams(i)%kill_gaussbeam_src()
    end do

end subroutine kill_sources

!###################################################################################################

subroutine propagate_p_srcs(this, time)

    class(TSources_list), intent(inout) :: this
    real(dp)            , intent(in)  :: time

    integer :: i

    do i = 1, this%n_p_src
        call this%points(i)%propagate_point_src(time)
    end do

end subroutine propagate_p_srcs

!###################################################################################################

subroutine propagate_pw_srcs(this, time)

    class(TSources_list), intent(inout) :: this
    real(dp)            , intent(in)  :: time

    integer :: i

    do i = 1, this%n_pw_src
        call this%plane_waves(i)%propagate_plane_wave_src(time)
    end do

end subroutine propagate_pw_srcs

!###################################################################################################

subroutine init_point_src(this, input_ch, dim, dt, dr, grid_Ndims, mpi_coords, mpi_dims)

        class(TPointSrc)   ,intent(inout) :: this
        character(len=1000),intent(in)    :: input_ch
        integer            ,intent(in)    :: grid_Ndims(3)
        integer            ,intent(in)    :: mpi_coords(3)
        integer            ,intent(in)    :: mpi_dims(3)
        integer            ,intent(in)    :: dim
        real(dp)           ,intent(in)    :: dr
        real(dp)           ,intent(in)    :: dt


        character(len=2)   :: polarization
        character(len=50)  :: type_src_ch
        real(dp)           :: j_amp
        real(dp)           :: freq
        real(dp)           :: t0
        real(dp)           :: tau
        real(dp)           :: t_init
        real(dp)           :: t_final
        real(dp)           :: phase
        real(dp)           :: r0(3)
        real(dp)           :: radius
        
        integer :: n_ker, i, j

        read(input_ch, *) type_src_ch, polarization, j_amp, freq, t0, tau, &
                          r0(1), r0(2), r0(3), radius, t_init, t_final, phase

        this%dim          = dim
        this%polarization = polarization
        this%j_amp        = j_amp
        this%w0           = freq*ev_to_au
        this%t0           = t0*fs_to_au
        this%tau          = tau*fs_to_au
        this%r0           = r0*nm_to_au
        this%rad          = radius*nm_to_au
        this%t_init       = t_init*fs_to_au
        this%t_final      = t_final*fs_to_au
        this%phase        = phase

        n_ker = int(this%rad/dr)

        this%n_ker = n_ker

        select case (this%dim)
        case (1)
            if (.not. allocated(this%ind_i))        allocate(this%ind_i(-n_ker:n_ker, 1, 1))
            if (.not. allocated(this%ker_mat))      allocate(this%ker_mat(-n_ker:n_ker, 1, 1))
            if (.not. allocated(this%J_mat))        allocate(this%J_mat(-n_ker:n_ker, 1, 1))
            if (.not. allocated(this%in_this_rank)) allocate(this%in_this_rank(-n_ker:n_ker, 1, 1))
        case (2)
            if (.not. allocated(this%ind_i))        allocate(this%ind_i(-n_ker:n_ker, -n_ker:n_ker, 1))
            if (.not. allocated(this%ind_j))        allocate(this%ind_j(-n_ker:n_ker, -n_ker:n_ker, 1))
            if (.not. allocated(this%ker_mat))      allocate(this%ker_mat(-n_ker:n_ker, -n_ker:n_ker, 1))
            if (.not. allocated(this%J_mat))        allocate(this%J_mat(-n_ker:n_ker, -n_ker:n_ker, 1))
            if (.not. allocated(this%in_this_rank)) allocate(this%in_this_rank(-n_ker:n_ker, -n_ker:n_ker, 1))
        case (3)
            if (.not. allocated(this%ind_i))        allocate(this%ind_i(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))
            if (.not. allocated(this%ind_j))        allocate(this%ind_j(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))
            if (.not. allocated(this%ind_k))        allocate(this%ind_k(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))   
            if (.not. allocated(this%ker_mat))      allocate(this%ker_mat(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))
            if (.not. allocated(this%J_mat))        allocate(this%J_mat(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))
            if (.not. allocated(this%in_this_rank)) allocate(this%in_this_rank(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))
        end select

        call compute_kernel(this, dr)


        this%in_this_rank=.true.
        call determine_indx_and_ranks(this, dr, grid_Ndims, mpi_coords, mpi_dims)

        ! do j = -n_ker, n_ker
        ! do i = -n_ker, n_ker
        !     write(*,*) i, j, this%J_mat(i,j,1), this%in_this_rank(i,j,1)
        ! end do
        ! end do

        ! stop

end subroutine init_point_src

!###################################################################################################

subroutine kill_point_src(this)

    class(TPointSrc), intent(inout) :: this

    if (allocated(this%ind_i))        deallocate(this%ind_i)
    if (allocated(this%ind_j))        deallocate(this%ind_j)
    if (allocated(this%ind_k))        deallocate(this%ind_k)
    if (allocated(this%ker_mat))      deallocate(this%ker_mat)
    if (allocated(this%J_mat))        deallocate(this%J_mat)
    if (allocated(this%in_this_rank)) deallocate(this%in_this_rank)

end subroutine kill_point_src

!###################################################################################################

subroutine propagate_point_src(this, time)

    class(TPointSrc), intent(inout) :: this
    real(dp)     , intent(in)    :: time

    real(dp) :: envelope
    real(dp) :: cos_t

    if (time >= this%t_init .and. time <= this%t_final) then

        envelope = DEXP( -((time - this%t0)/this%tau)**2 )
        cos_t    = DCOS( this%w0*(time - this%t0) + this%phase)

        this%J_mat = this%j_amp * envelope * cos_t * this%ker_mat

    else

        this%J_mat = 0.0d0

    end if

end subroutine propagate_point_src

!###################################################################################################

subroutine init_plane_wave_src(this, input_ch, dim, dt, dr, grid_Ndims, mpi_coords, mpi_dims)

    class(TPlaneWaveSrc) ,intent(inout) :: this
    character(len=1000)  ,intent(in)    :: input_ch
    integer              ,intent(in)    :: grid_Ndims(3)
    integer              ,intent(in)    :: mpi_coords(3)
    integer              ,intent(in)    :: mpi_dims(3)
    integer              ,intent(in)    :: dim
    real(dp)             ,intent(in)    :: dr
    real(dp)             ,intent(in)    :: dt

    character(len=50)  :: type_src_ch
    integer            :: i_min, i_max, j_min, j_max, k_min, k_max
    integer            :: rank_x, rank_y, rank_z
    integer            :: aux_grid_Ndim(3)
    integer            :: aux_vec_mpi_coords(3)
    integer            :: aux_boundaries(3)
    real(dp)           :: E_amp
    real(dp)           :: dr_1D
    real(dp)           :: phi
    real(dp)           :: theta
    real(dp)           :: psi
    real(dp)           :: freq
    real(dp)           :: t0
    real(dp)           :: tau
    real(dp)           :: t_init
    real(dp)           :: t_final
    real(dp)           :: phase
    real(dp)           :: x_min, x_max
    real(dp)           :: y_min, y_max
    real(dp)           :: z_min, z_max

    read(input_ch, *) type_src_ch, E_amp, phi, theta, psi, freq, t0, tau, &
                      x_min, x_max, y_min, y_max, z_min, z_max, &
                      t_init, t_final, phase

    select case (trim(type_src_ch))
    case ("plane_wave_kx")
        this%limited_axis(2) = .false.
        this%limited_axis(3) = .false.
        phi                  =  0.0d0
        theta                = 90.0d0
    case ("plane_wave_ky")
        this%limited_axis(1) = .false.
        this%limited_axis(3) = .false.
        phi                  = 90.0d0
        theta                = 90.0d0
    case ("plane_wave_kz")
        this%limited_axis(1) = .false.
        this%limited_axis(2) = .false.
        phi                  = 0.0d0
        theta                = 0.0d0
    end select

    this%dim   = dim
    this%E_amp = E_amp

    this%phi    = phi/180.0d0*pi0
    this%theta  = theta/180.0d0*pi0
    this%psi    = psi/180.0d0*pi0

    if (this%phi < 0.0d0)   this%phi   = this%phi + 2*pi0
    if (this%theta < 0.0d0) this%theta = this%theta + 2*pi0
    if (this%psi < 0.0d0)   this%psi   = this%psi + 2*pi0

    
    this%w0      = freq * ev_to_au
    this%t0      = t0   * fs_to_au
    this%tau     = tau  * fs_to_au
    this%t_init  = t_init  * fs_to_au
    this%t_final = t_final * fs_to_au
    this%phase   = phase

    x_min  = x_min * nm_to_au
    x_max  = x_max * nm_to_au
    y_min  = y_min * nm_to_au
    y_max  = y_max * nm_to_au
    z_min  = z_min * nm_to_au
    z_max  = z_max * nm_to_au

    if (x_min < (1-int(grid_Ndims(1)*mpi_dims(1)/2))*dr .or. x_max > (int(grid_Ndims(1)*mpi_dims(1)/2))*dr) then
        write (*, '("Error: the plane wave source extends beyond the simulation box in x direction &
                & [", F10.4, ", ", F10.4, "].")') (1-int(grid_Ndims(1)*mpi_dims(1)/2))*dr/nm_to_au, &
                                                  (int(grid_Ndims(1)*mpi_dims(1)/2))*dr/nm_to_au
        error stop
    end if

    if (dim > 1) then

        if (y_min < (1-int(grid_Ndims(2)*mpi_dims(2)/2))*dr .or. y_max > (int(grid_Ndims(2)*mpi_dims(2)/2))*dr) then
            write (*, '("Error: the plane wave source extends beyond the simulation box in y direction &
                    & [", F10.4, ", ", F10.4, "].")') (1-int(grid_Ndims(2)*mpi_dims(2)/2))*dr/nm_to_au, &
                                                    (int(grid_Ndims(2)*mpi_dims(2)/2))*dr/nm_to_au
            error stop
        end if

    end if

    if (dim == 3) then

        if (z_min < (1-int(grid_Ndims(3)*mpi_dims(3)/2))*dr .or. z_max > (int(grid_Ndims(3)*mpi_dims(3)/2))*dr) then
            write (*, '("Error: the plane wave source extends beyond the simulation box in z direction &
                    & [", F10.4, ", ", F10.4, "].")') (1-int(grid_Ndims(3)*mpi_dims(3)/2))*dr/nm_to_au, &
                                                    (int(grid_Ndims(3)*mpi_dims(3)/2))*dr/nm_to_au
            error stop
        end if
    end if

    i_min  = FLOOR(x_min/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
    i_max  = FLOOR(x_max/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
    j_min  = FLOOR(y_min/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
    j_max  = FLOOR(y_max/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
    k_min  = FLOOR(z_min/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)
    k_max  = FLOOR(z_max/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)

    if (.not. this%limited_axis(1)) then
        i_min = 1
        i_max = grid_Ndims(1)*mpi_dims(1)
    end if
    
    if (.not. this%limited_axis(2) .and. this%dim > 1) then
        j_min = 1
        j_max = grid_Ndims(2)*mpi_dims(2)
    end if

    if (.not. this%limited_axis(3) .and. this%dim == 3) then
        k_min = 1
        k_max = grid_Ndims(3)*mpi_dims(3)
    end if

    this%i_min = i_min
    this%i_max = i_max
    this%j_min = j_min
    this%j_max = j_max
    this%k_min = k_min
    this%k_max = k_max
    
    this%i_min_in_this_rank = .false.
    this%i_max_in_this_rank = .false.
    this%j_min_in_this_rank = .false.
    this%j_max_in_this_rank = .false.
    this%k_min_in_this_rank = .false.
    this%k_max_in_this_rank = .false.
    
    rank_x = int((i_min-1)/grid_Ndims(1))

    if (rank_x == mpi_coords(1)) then
        this%i_min_in_this_rank = .true.
        this%i_min_loc          = i_min - rank_x*grid_Ndims(1)
    end if

    rank_x = int((i_max-1)/grid_Ndims(1))

    if (rank_x == mpi_coords(1)) then
        this%i_max_in_this_rank = .true.
        this%i_max_loc          = i_max - rank_x*grid_Ndims(1)
    end if

    if (this%dim > 1) then

        rank_y = int((j_min-1)/grid_Ndims(2))

        if (rank_y == mpi_coords(2)) then
            this%j_min_in_this_rank = .true.
            this%j_min_loc          = j_min - rank_y*grid_Ndims(2)
        end if

        rank_y = int((j_max-1)/grid_Ndims(2))

        if (rank_y == mpi_coords(2)) then
            this%j_max_in_this_rank = .true.
            this%j_max_loc          = j_max - rank_y*grid_Ndims(2)
        end if

    end if

    if (this%dim == 3) then

        rank_z = int((k_min-1)/grid_Ndims(3))

        if (rank_z == mpi_coords(3)) then
            this%k_min_in_this_rank = .true.
            this%k_min_loc          = k_min - rank_z*grid_Ndims(3)
        end if

        rank_z = int((k_max-1)/grid_Ndims(3))

        if (rank_z == mpi_coords(3)) then
            this%k_max_in_this_rank = .true.
            this%k_max_loc          = k_max - rank_z*grid_Ndims(3)
        end if

    end if


    !Approximation to the exact matched numerical dispersion method for plane waves
    !to adjust the grid spacing of the auxiliary 1D grid.

    select case (this%dim)
    case (1)
        dr_1D = dr
    case (2)
        this%theta = pi0/2.0d0 
        dr_1D      = dr * SQRT( DCOS(this%phi)**4 + DSIN(this%phi)**4 )
        this%v_vec =  (/DCOS(this%phi), DSIN(this%phi), 0.0d0/)
        this%A_vec = 0.0d0 

        if (phi >= 0.0_dp .and. phi <= 90.0_dp) then
            this%A_vec(1) = (i_min - int(grid_Ndims(1)*mpi_dims(1)/2))*dr
            this%A_vec(2) = (j_min - int(grid_Ndims(2)*mpi_dims(2)/2))*dr
        else if (phi > 90.0_dp .and. phi <= 180.0_dp) then
            this%A_vec(1) = (i_max - int(grid_Ndims(1)*mpi_dims(1)/2))*dr
            this%A_vec(2) = (j_min - int(grid_Ndims(2)*mpi_dims(2)/2))*dr
        else if (phi > 180.0_dp .and. phi <= 270.0_dp) then
            this%A_vec(1) = (i_max - int(grid_Ndims(1)*mpi_dims(1)/2))*dr
            this%A_vec(2) = (j_max - int(grid_Ndims(2)*mpi_dims(2)/2))*dr
        else if (phi > 270.0_dp .and. phi <= 360.0_dp) then
            this%A_vec(1) = (i_min - int(grid_Ndims(1)*mpi_dims(1)/2))*dr
            this%A_vec(2) = (j_max - int(grid_Ndims(2)*mpi_dims(2)/2))*dr
        end if
        
    case (3)
        dr_1D = dr * SQRT( DCOS(this%theta)**4 + &
                           DSIN(this%theta)**4 * (DCOS(this%phi)**4 + DSIN(this%phi)**4) )

        this%v_vec =  (/DCOS(this%phi)*DSIN(this%theta), &
                       DSIN(this%phi)*DSIN(this%theta), &
                       DCOS(this%theta)/)

        this%A_vec = 0.0d0

        if (phi >= 0.0_dp .and. phi <= 90.0_dp) then
            this%A_vec(1) = (i_min - int(grid_Ndims(1)*mpi_dims(1)/2))*dr
            this%A_vec(2) = (j_min - int(grid_Ndims(2)*mpi_dims(2)/2))*dr
        else if (phi > 90.0_dp .and. phi <= 180.0_dp) then
            this%A_vec(1) = (i_max - int(grid_Ndims(1)*mpi_dims(1)/2))*dr
            this%A_vec(2) = (j_min - int(grid_Ndims(2)*mpi_dims(2)/2))*dr
        else if (phi > 180.0_dp .and. phi <= 270.0_dp) then
            this%A_vec(1) = (i_max - int(grid_Ndims(1)*mpi_dims(1)/2))*dr
            this%A_vec(2) = (j_max - int(grid_Ndims(2)*mpi_dims(2)/2))*dr
        else if (phi > 270.0_dp .and. phi <= 360.0_dp) then
            this%A_vec(1) = (i_min - int(grid_Ndims(1)*mpi_dims(1)/2))*dr
            this%A_vec(2) = (j_max - int(grid_Ndims(2)*mpi_dims(2)/2))*dr
        end if

        if (theta >= 0.0_dp .and. theta <= 90.0_dp) then
            this%A_vec(3) = (k_min - int(grid_Ndims(3)*mpi_dims(3)/2))*dr
        else if (theta > 90.0_dp .and. theta <= 180.0_dp) then
            this%A_vec(3) = (k_max - int(grid_Ndims(3)*mpi_dims(3)/2))*dr
        end if

    end select 

    !Estimating the number of points of the auxiliary 1D grid. We consider two-times more points
    !than the maximum distance that the plane wave can propagate in the simulation box. 

    aux_grid_Ndim    = 0
    aux_grid_Ndim(1) = int(2*SQRT((x_max-x_min)**2 + (y_max-y_min)**2 + (z_max-z_min)**2)/dr_1D)

    !TO-DO: the next vector is used to force the 1D case to run in in every rank, but it should
    !be changed in the future to be more general.
    aux_vec_mpi_coords = (/0, 0, 0/)


    aux_boundaries = (/CPML_BOUNDARIES, CLOSE_BOUNDARIES, CLOSE_BOUNDARIES/)
    
    call this%mxll_inc%init(grid_Ndims = aux_grid_Ndim, npml=20, boundaries=aux_boundaries, &
                            dt = dt, dr = dr_1D, mode = AUX_GRID_MODE, n_media = 0, &
                            mpi_coords = aux_vec_mpi_coords, mpi_dims = mpi_dims)

end subroutine init_plane_wave_src

!###################################################################################################

subroutine kill_plane_wave_src(this)

    class(TPlaneWaveSrc), intent(inout) :: this

    call this%mxll_inc%kill()

end subroutine kill_plane_wave_src

!###################################################################################################

subroutine propagate_plane_wave_src(this, time)

    class(TPlaneWaveSrc), intent(inout) :: this
    real(dp)            , intent(in)    :: time

    real(dp) :: envelope
    real(dp) :: cos_t

    if (time >= this%t_init .and. time <= this%t_final) then

        envelope = DEXP( -((time - this%t0)/this%tau)**2 )
        cos_t    = DCOS( this%w0*(time - this%t0) + this%phase)

        call this%mxll_inc%td_propagate_H_field()

        this%mxll_inc%Ex(1) = this%E_amp * envelope * cos_t

        call this%mxll_inc%td_propagate_E_field(0)

    else

        call this%mxll_inc%td_propagate_H_field()
        
        this%mxll_inc%Ex(1) = 0.0d0
        
        call this%mxll_inc%td_propagate_E_field(0)

    end if

end subroutine propagate_plane_wave_src

!###################################################################################################

subroutine init_gaussbeam_src(this, input_ch, dim, dt, dr, grid_Ndims, mpi_coords, mpi_dims)

    class(TGaussbeamSrc) ,intent(inout) :: this
    character(len=1000)  ,intent(in)    :: input_ch
    integer              ,intent(in)    :: grid_Ndims(3)
    integer              ,intent(in)    :: mpi_coords(3)
    integer              ,intent(in)    :: mpi_dims(3)
    integer              ,intent(in)    :: dim
    real(dp)             ,intent(in)    :: dr
    real(dp)             ,intent(in)    :: dt

    character(len=50)  :: type_src_ch
    integer            :: i_min, i_max, j_min, j_max, k_min, k_max
    integer            :: rank_x, rank_y, rank_z

    real(dp)           :: E_amp
    integer            :: aux_grid_Ndim(3)
    integer            :: aux_vec_mpi_coords(3)
    integer            :: aux_boundaries(3)
    real(dp)           :: dr_1D
    real(dp)           :: phi
    real(dp)           :: theta
    real(dp)           :: psi
    real(dp)           :: freq
    real(dp)           :: dz_ramp
    real(dp)           :: z0_ramp
    real(dp)           :: phase
    real(dp)           :: r0(3)
    real(dp)           :: w0
    real(dp)           :: x_min, x_max
    real(dp)           :: y_min, y_max
    real(dp)           :: z_min, z_max

    read(input_ch, *) type_src_ch, E_amp, phi, theta, psi, freq, &
                      x_min, x_max, y_min, y_max, z_min, z_max, &
                      phase, r0(1), r0(2), r0(3), w0, dz_ramp, z0_ramp

    this%dim   = dim
    this%E_amp = E_amp

    this%phi    = phi/180.0d0*pi0
    this%theta  = theta/180.0d0*pi0
    this%psi    = psi/180.0d0*pi0

    if (this%phi < 0.0d0)   this%phi   = this%phi + 2*pi0
    if (this%theta < 0.0d0) this%theta = this%theta + 2*pi0
    if (this%psi < 0.0d0)   this%psi   = this%psi + 2*pi0

    
    this%freq    = freq * ev_to_au
    this%phase   = phase
    this%r0      = r0 * nm_to_au
    this%w0      = w0 * nm_to_au
    this%lambda  = 2*pi0*c0/this%freq
    this%k       = 2*pi0/this%lambda
    this%z_R     = pi0*this%w0**2/this%lambda
    this%dz_ramp = dz_ramp * nm_to_au
    this%z0_ramp = z0_ramp * nm_to_au

    x_min  = x_min * nm_to_au
    x_max  = x_max * nm_to_au
    y_min  = y_min * nm_to_au
    y_max  = y_max * nm_to_au
    z_min  = z_min * nm_to_au
    z_max  = z_max * nm_to_au

    if (x_min < (1-int(grid_Ndims(1)*mpi_dims(1)/2))*dr .or. x_max > (int(grid_Ndims(1)*mpi_dims(1)/2))*dr) then
        write (*, '("Error: the gauss beam source extends beyond the simulation box in x direction &
                & [", F10.4, ", ", F10.4, "].")') (1-int(grid_Ndims(1)*mpi_dims(1)/2))*dr/nm_to_au, &
                                                  (int(grid_Ndims(1)*mpi_dims(1)/2))*dr/nm_to_au
        error stop
    end if

    if (dim > 1) then

        if (y_min < (1-int(grid_Ndims(2)*mpi_dims(2)/2))*dr .or. y_max > (int(grid_Ndims(2)*mpi_dims(2)/2))*dr) then
            write (*, '("Error: the gauss beam source extends beyond the simulation box in y direction &
                    & [", F10.4, ", ", F10.4, "].")') (1-int(grid_Ndims(2)*mpi_dims(2)/2))*dr/nm_to_au, &
                                                    (int(grid_Ndims(2)*mpi_dims(2)/2))*dr/nm_to_au
            error stop
        end if

    end if

    if (dim == 3) then

        if (z_min < (1-int(grid_Ndims(3)*mpi_dims(3)/2))*dr .or. z_max > (int(grid_Ndims(3)*mpi_dims(3)/2))*dr) then
            write (*, '("Error: the gauss beam source extends beyond the simulation box in z direction &
                    & [", F10.4, ", ", F10.4, "].")') (1-int(grid_Ndims(3)*mpi_dims(3)/2))*dr/nm_to_au, &
                                                    (int(grid_Ndims(3)*mpi_dims(3)/2))*dr/nm_to_au
            error stop
        end if
    end if

    i_min  = FLOOR(x_min/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
    i_max  = FLOOR(x_max/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
    j_min  = FLOOR(y_min/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
    j_max  = FLOOR(y_max/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
    k_min  = FLOOR(z_min/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)
    k_max  = FLOOR(z_max/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)
    
    this%i_min = i_min
    this%i_max = i_max
    this%j_min = j_min
    this%j_max = j_max
    this%k_min = k_min
    this%k_max = k_max
    
    this%i_min_in_this_rank = .false.
    this%i_max_in_this_rank = .false.
    this%j_min_in_this_rank = .false.
    this%j_max_in_this_rank = .false.
    this%k_min_in_this_rank = .false.
    this%k_max_in_this_rank = .false.
    
    rank_x = int((i_min-1)/grid_Ndims(1))

    if (rank_x == mpi_coords(1)) then
        this%i_min_in_this_rank = .true.
        this%i_min_loc          = i_min - rank_x*grid_Ndims(1)
    end if

    rank_x = int((i_max-1)/grid_Ndims(1))

    if (rank_x == mpi_coords(1)) then
        this%i_max_in_this_rank = .true.
        this%i_max_loc          = i_max - rank_x*grid_Ndims(1)
    end if

    if (this%dim > 1) then

        rank_y = int((j_min-1)/grid_Ndims(2))

        if (rank_y == mpi_coords(2)) then
            this%j_min_in_this_rank = .true.
            this%j_min_loc          = j_min - rank_y*grid_Ndims(2)
        end if

        rank_y = int((j_max-1)/grid_Ndims(2))

        if (rank_y == mpi_coords(2)) then
            this%j_max_in_this_rank = .true.
            this%j_max_loc          = j_max - rank_y*grid_Ndims(2)
        end if

    end if

    if (this%dim == 3) then

        rank_z = int((k_min-1)/grid_Ndims(3))

        if (rank_z == mpi_coords(3)) then
            this%k_min_in_this_rank = .true.
            this%k_min_loc          = k_min - rank_z*grid_Ndims(3)
        end if

        rank_z = int((k_max-1)/grid_Ndims(3))

        if (rank_z == mpi_coords(3)) then
            this%k_max_in_this_rank = .true.
            this%k_max_loc          = k_max - rank_z*grid_Ndims(3)
        end if

    end if

    !Approximation to the exact matched numerical dispersion method for plane waves
    !to adjust the grid spacing of the auxiliary 1D grid.

    select case (this%dim)
    case (1)
        dr_1D = dr
    case (2)
        this%theta = pi0/2.0d0 
        this%v_vec =  (/DCOS(this%phi), DSIN(this%phi), 0.0d0/)

    case (3)

        this%v_vec =  (/DCOS(this%phi)*DSIN(this%theta), &
                       DSIN(this%phi)*DSIN(this%theta), &
                       DCOS(this%theta)/)

    end select 

end subroutine init_gaussbeam_src

!###################################################################################################

subroutine kill_gaussbeam_src(this)

    class(TGaussbeamSrc), intent(inout) :: this

    !Currently, there are no dynamic resources to free.

end subroutine kill_gaussbeam_src
!###################################################################################################

subroutine compute_time_space_profile(this, r, z, time)

    class(TGaussbeamSrc), intent(inout) :: this
    real(dp)            , intent(in)    :: r
    real(dp)            , intent(in)    :: z
    real(dp)            , intent(in)    :: time

    complex(dp) :: E_t
    complex(dp) :: E_rz
    real(dp)    :: envelope
    real(dp)    :: cos_t
    real(dp)    :: sin_t
    real(dp)    :: z_ramp
    real(dp)    :: z0_ramp
    real(dp)    :: z_min
    real(dp)    :: z_max
    real(dp)    :: dz_ramp
    real(dp) :: w_z
    real(dp) :: inv_R_z
    real(dp) :: psi_z

    z0_ramp = this%z0_ramp
    dz_ramp = this%dz_ramp

    z_ramp = time*c0 + z0_ramp

    z_min = z_ramp - 0.5_dp*dz_ramp
    z_max = z_ramp + 0.5_dp*dz_ramp


    envelope = EXP(-(z-z_ramp)**2/(2*dz_ramp**2))

    E_t   = this%E_amp * envelope * (DCOS(this%freq*(time)) + Z_I*DSIN(this%freq*(time)))

    w_z   = this%w0*SQRT(1.0d0 + (z/this%z_R)**2)
    inv_R_z = z/(z**2 + this%z_R**2)
    psi_z = ATAN(z/this%z_R)

    E_rz = (this%w0/w_z)*EXP(-r**2/w_z**2)*(DCOS(this%k*z + (this%k*r**2)*(0.5d0*inv_R_z)-psi_z) &
                                      - Z_I*DSIN(this%k*z + (this%k*r**2)*(0.5d0*inv_R_z)-psi_z))

    this%E_rzt = E_t * E_rz

end subroutine compute_time_space_profile


!###################################################################################################
!######################################GENERAL SUBS#################################################
!###################################################################################################

subroutine compute_kernel(this, dr)

    class(TPointSrc), intent(inout) :: this
    real(dp)        , intent(in)    :: dr

    integer  :: i, j, k
    integer  :: n_ker
    real(dp) :: x, y, z
    real(dp) :: r_max, r, r0
    real(dp) :: norm

    n_ker = this%n_ker

    r_max = (this%n_ker*2 + 2)*dr
    r0    = r_max / 2.0d0

    norm = 0.0d0
    select case (this%dim)
    case (1)
        do i = -n_ker, n_ker

            x = i*dr

            this%ker_mat(i,1,1) = aBH(1)+ &
                          aBH(2)*DCOS(2.0*pi0*(x-r0)/r_max)+ &
                          aBH(3)*DCOS(2.0*pi0*2.0*(x-r0)/r_max)+ &
                          aBH(4)*DCOS(2.0*pi0*3.0*(x-r0)/r_max)

            norm = norm + this%ker_mat(i,1,1)
        end do
        
    case (2)
        do j = -n_ker, n_ker
        do i = -n_ker, n_ker
                
            x = i*dr
            y = j*dr
            
            r = SQRT(x**2 + y**2)
            
            this%ker_mat(i,j,1) = aBH(1)+ &
            aBH(2)*DCOS(2.0*pi0*(r-r0)/r_max)+ &
            aBH(3)*DCOS(2.0*pi0*2.0*(r-r0)/r_max)+ &
            aBH(4)*DCOS(2.0*pi0*3.0*(r-r0)/r_max)
            
            if (this%ker_mat(i,j,1) < 0.0d0) this%ker_mat(i,j,1) = 0.0d0
            
            norm = norm + this%ker_mat(i,j,1)
        end do
        end do
    case (3)
        do k = -n_ker, n_ker
        do j = -n_ker, n_ker
        do i = -n_ker, n_ker

            x = i*dr
            y = j*dr
            z = k*dr

            r = SQRT(x**2 + y**2 + z**2)

            this%ker_mat(i,j,k) = aBH(1)+ &
                          aBH(2)*DCOS(2.0*pi0*(r-r0)/r_max)+ &
                          aBH(3)*DCOS(2.0*pi0*2.0*(r-r0)/r_max)+ &
                          aBH(4)*DCOS(2.0*pi0*3.0*(r-r0)/r_max)

            if (this%ker_mat(i,j,k) < 0.0d0) this%ker_mat(i,j,k) = 0.0d0

            norm = norm + this%ker_mat(i,j,k)
        end do
        end do
        end do
    end select

    this%ker_mat = this%ker_mat / norm

end subroutine compute_kernel

!###################################################################################################

subroutine determine_indx_and_ranks(this, dr, grid_Ndims, mpi_coords, mpi_dims)

    class(TPointSrc) , intent(inout) :: this
    real(dp)         , intent(in)    :: dr
    integer          , intent(in)    :: grid_Ndims(3)
    integer          , intent(in)    :: mpi_coords(3)
    integer          , intent(in)    :: mpi_dims(3)

    integer :: i, j, k
    integer :: rank_x, rank_y, rank_z

    !By default, we consider that the origin is where global_i, global_j and/or global_k
    !reach the half of the global grid size.

    select case (this%dim)
    case (1)
        do i = -this%n_ker, this%n_ker
            this%ind_i(i,1,1) = i+int(this%r0(1)/dr) + int(grid_Ndims(1)/2)
        end do
    case (2)
        do j = -this%n_ker, this%n_ker
        do i = -this%n_ker, this%n_ker
            this%ind_i(i,j,1) = i+int(this%r0(1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
            this%ind_j(i,j,1) = j+int(this%r0(2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
        end do
        end do
    case (3)
        do k = -this%n_ker, this%n_ker
        do j = -this%n_ker, this%n_ker
        do i = -this%n_ker, this%n_ker
            this%ind_i(i,j,k) = i+int(this%r0(1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
            this%ind_j(i,j,k) = j+int(this%r0(2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
            this%ind_k(i,j,k) = k+int(this%r0(3)/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)
        end do
        end do
        end do
    end select

#ifdef USE_MPI

    select case (this%dim)
    case (1)
    case (2)
        do j = -this%n_ker, this%n_ker
        do i = -this%n_ker, this%n_ker
            rank_x = int((this%ind_i(i,j,1)-1)/grid_Ndims(1))
            rank_y = int((this%ind_j(i,j,1)-1)/grid_Ndims(2))

            if (rank_x == mpi_coords(1) .and. rank_y == mpi_coords(2)) then
                this%in_this_rank(i,j,1) = .true.
                this%ind_i(i,j,1) = this%ind_i(i,j,1) - rank_x*grid_Ndims(1)
                this%ind_j(i,j,1) = this%ind_j(i,j,1) - rank_y*grid_Ndims(2)
            else
                this%in_this_rank(i,j,1) = .false.
            end if

        end do
        end do
    case (3)
        do k = -this%n_ker, this%n_ker
        do j = -this%n_ker, this%n_ker
        do i = -this%n_ker, this%n_ker
            rank_x = int((this%ind_i(i,j,k)-1)/grid_Ndims(1))
            rank_y = int((this%ind_j(i,j,k)-1)/grid_Ndims(2))
            rank_z = int((this%ind_k(i,j,k)-1)/grid_Ndims(3))

            if (rank_x == mpi_coords(1) .and. rank_y == mpi_coords(2) .and. rank_z == mpi_coords(3)) then
                this%in_this_rank(i,j,k) = .true.
                this%ind_i(i,j,k) = this%ind_i(i,j,k) - rank_x*grid_Ndims(1)
                this%ind_j(i,j,k) = this%ind_j(i,j,k) - rank_y*grid_Ndims(2)
                this%ind_k(i,j,k) = this%ind_k(i,j,k) - rank_z*grid_Ndims(3)
            else
                this%in_this_rank(i,j,k) = .false.
            end if

        end do
        end do
        end do
    end select

#endif

end subroutine determine_indx_and_ranks

!###################################################################################################



!###################################################################################################

end module sources_subs_mod