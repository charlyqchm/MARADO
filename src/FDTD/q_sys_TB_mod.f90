module q_sys_TB_mod

    use constants_mod
    use q_sys_base_mod

    implicit none

    type, extends(TQ_sys_base) :: TQ_sys_TB

        integer   :: n_basis
        integer   :: n_particles
        logical   :: basis_mo
        logical   :: print_occ

        complex(dp), allocatable :: rho_at(:,:)
        complex(dp), allocatable :: rho_mo(:,:)
        complex(dp), allocatable :: rho(:,:)
        complex(dp), allocatable :: rho_new(:,:)
        complex(dp), allocatable :: rho_old(:,:)
        real(dp)   , allocatable :: H0_at(:,:)
        real(dp)   , allocatable :: H0_mo(:,:)
        real(dp)   , allocatable :: mu_x_at(:,:)
        real(dp)   , allocatable :: mu_y_at(:,:)
        real(dp)   , allocatable :: mu_z_at(:,:)
        real(dp)   , allocatable :: mu_x_mo(:,:)
        real(dp)   , allocatable :: mu_y_mo(:,:)
        real(dp)   , allocatable :: mu_z_mo(:,:)
        real(dp)   , allocatable :: c_vecs(:,:)
        real(dp)   , allocatable :: cT_vecs(:,:)
        real(dp)   , allocatable :: Ht(:,:)
        real(dp)   , allocatable :: H_tot(:,:)

        contains

            procedure :: init => init_TB
            procedure :: kill => kill_TB
            procedure :: gs_calculate => gs_calculate_TB
            procedure :: td_propagate => td_propagate_TB
            procedure :: write_output => write_output_TB

    end type TQ_sys_TB

contains

!###################################################################################################

subroutine init_TB(this, id, id_file, dt, t_steps, rank, print_on)

    class(TQ_sys_TB), intent(inout) :: this
    integer         , intent(in)    :: id
    integer         , intent(in)    :: id_file
    integer         , intent(in)    :: t_steps
    integer         , intent(in)    :: rank
    real(dp)        , intent(in)    :: dt
    logical         , intent(in)    :: print_on


    logical            :: set_initial_occ = .false.
    character(len=20)  :: file_name = "molecule_"
    character(len=20)  :: file_exten = ".in"
    character(len=20)  :: file_number
    character(len=20)  :: input_name
    character(len=30)  :: print_occ_str
    character(len=30)  :: external_src_on_str
    character(len=30)  :: set_initial_occ_str
    character(len=100) :: basis_type_str
    character(len=100) :: external_src_file
    character(len=200) :: iomsg
    integer            :: n_basis
    integer            :: n_particles
    integer            :: ierr
    integer            :: funit
    integer            :: i, j
    real(dp)           :: H_ij, mux_ij, muy_ij, muz_ij
    real(dp)           :: Re_rho_ij, Im_rho_ij

    this%id       = id
    this%id_file  = id_file
    this%dt       = dt
    this%t_steps  = t_steps
    this%rank     = rank
    this%print_on = print_on

    write(file_number,'(I7.7)') id_file
    input_name = trim(file_name)//trim(file_number)//trim(file_exten)

    open (action='read', file=input_name, iostat=ierr, newunit=funit)

    if (ierr /= 0) then
        write(*,*) "Error: Could not open molecule file ", trim(input_name)
        stop
    end if

    read (unit=funit, fmt=*, iostat=ierr) n_basis, n_particles
    read (unit=funit, fmt=*, iostat=ierr) basis_type_str, print_occ_str, set_initial_occ_str
    
    this%n_basis = n_basis
    this%n_particles = n_particles
    
    this%print_occ  = .false.
    set_initial_occ = .false.
    
    
    if (ierr /= 0) then
        write(*,*) "Error: Could not read basis type, print occupancy flag and"
        write(*,*) "initial occupancy flag from file ", trim(input_name)
        stop
    end if
    
    read (unit=funit, fmt=*, iostat=ierr) external_src_on_str
    
    if (ierr /= 0) then
        write(*,*) "Error: Could not read external source flag from file ", trim(input_name)
        stop
    end if
    
    if (external_src_on_str == "external_src_on") then
        this%external_src = .true.
        read (unit=funit, fmt=*, iostat=ierr) external_src_file

        if (ierr /= 0) then
            write(*,*) "Error: Could not read external source file name from file ", &
                        trim(input_name)
            stop
        end if

    else if (external_src_on_str == "external_src_off") then
        this%external_src = .false.
    else
        write(*,*) "Error: Unknown external source flag ", trim(external_src_on_str)
        stop
    end if


    if (this%print_on) then
        if (print_occ_str == "print_occ_on") then
            this%print_occ = .true.
        else if (print_occ_str == "print_occ_off") then
            this%print_occ = .false.
        else
            write(*,*) "Error: Unknown print occupancy flag ", trim(print_occ_str)
            stop
        end if
    end if

    select case (trim(set_initial_occ_str))
    case ('set_initial_occ_on')
        set_initial_occ = .true.
    case ('set_initial_occ_off')
        set_initial_occ = .false.
    case default
        write(*,*) "Error: Unknown initial occupancy flag ", trim(set_initial_occ_str)
        stop
    end select

    select case (trim(basis_type_str))
    case ('atomic_basis')
        this%basis_mo = .false.
    case ('molecular_orbital_basis')
        this%basis_mo = .true.
    case default
        write(*,*) "Error: Unknown TB basis type ", trim(basis_type_str)
        stop
    end select

    if (.not. allocated(this%H0_mo))      allocate(this%H0_mo(this%n_basis, this%n_basis))
    if (.not. allocated(this%rho_mo))     allocate(this%rho_mo(this%n_basis, this%n_basis))
    if (.not. allocated(this%rho))        allocate(this%rho(this%n_basis, this%n_basis))
    if (.not. allocated(this%rho_old))    allocate(this%rho_old(this%n_basis, this%n_basis))
    if (.not. allocated(this%rho_new))    allocate(this%rho_new(this%n_basis, this%n_basis))
    if (.not. allocated(this%mu_x_mo))    allocate(this%mu_x_mo(this%n_basis, this%n_basis))
    if (.not. allocated(this%mu_y_mo))    allocate(this%mu_y_mo(this%n_basis, this%n_basis))
    if (.not. allocated(this%mu_z_mo))    allocate(this%mu_z_mo(this%n_basis, this%n_basis))
    if (.not. allocated(this%Ht))         allocate(this%Ht(this%n_basis, this%n_basis))
    if (.not. allocated(this%H_tot))      allocate(this%H_tot(this%n_basis, this%n_basis))
    if (.not. allocated(this%dipole))     allocate(this%dipole(3))
    if (.not. allocated(this%dip_old))    allocate(this%dip_old(3))
    if (.not. allocated(this%dPt_dt))     allocate(this%dPt_dt(3))
    if (.not. allocated(this%dPt_dt_old)) allocate(this%dPt_dt_old(3))

    this%H0_mo   = M_ZERO
    this%rho_mo  = Z_0
    this%rho     = Z_0
    this%rho_old = Z_0
    this%rho_new = Z_0
    this%mu_x_mo = M_ZERO
    this%mu_y_mo = M_ZERO
    this%mu_z_mo = M_ZERO
    this%Ht      = M_ZERO    
    this%H_tot   = M_ZERO

    this%dipole     = M_ZERO
    this%dip_old    = M_ZERO
    this%dPt_dt     = M_ZERO
    this%dPt_dt_old = M_ZERO  

    if (.not. this%basis_mo) then

        if (.not. allocated(this%H0_at))   allocate(this%H0_at(this%n_basis, this%n_basis))
        if (.not. allocated(this%rho_at))  allocate(this%rho_at(this%n_basis, this%n_basis))
        if (.not. allocated(this%mu_x_at)) allocate(this%mu_x_at(this%n_basis, this%n_basis))
        if (.not. allocated(this%mu_y_at)) allocate(this%mu_y_at(this%n_basis, this%n_basis))
        if (.not. allocated(this%mu_z_at)) allocate(this%mu_z_at(this%n_basis, this%n_basis))
        if (.not. allocated(this%c_vecs))  allocate(this%c_vecs(this%n_basis, this%n_basis))
        if (.not. allocated(this%cT_vecs)) allocate(this%cT_vecs(this%n_basis, this%n_basis))

        this%H0_at   = M_ZERO
        this%rho_at  = Z_0
        this%mu_x_at = M_ZERO
        this%mu_y_at = M_ZERO
        this%mu_z_at = M_ZERO
        this%c_vecs  = M_ZERO
        this%cT_vecs = M_ZERO

    end if

    do

        if (set_initial_occ) then
            read (unit=funit, fmt=*, iostat=ierr) i, j, H_ij, mux_ij, muy_ij, muz_ij, Re_rho_ij, Im_rho_ij
            if (ierr /= 0) exit
        else
            read (unit=funit, fmt=*, iostat=ierr) i, j, H_ij, mux_ij, muy_ij, muz_ij
            if (ierr /= 0) exit
        end if

        if (this%basis_mo) then
            
            this%H0_mo(i,j)   = H_ij*ev_to_au
            this%mu_x_mo(i,j) = mux_ij
            this%mu_y_mo(i,j) = muy_ij
            this%mu_z_mo(i,j) = muz_ij

        else

            this%H0_at(i,j)   = H_ij*ev_to_au
            this%mu_x_at(i,j) = mux_ij
            this%mu_y_at(i,j) = muy_ij
            this%mu_z_at(i,j) = muz_ij

        end if

        if (set_initial_occ) this%rho_mo(i,j) = Re_rho_ij + Z_I*Im_rho_ij

    end do

    close(unit=funit)

    if (.not. set_initial_occ) then

        do i = 1, this%n_particles
            this%rho_mo(i,i) = Z_ONE
        end do

    end if

    if (this%external_src) then

        this%n_steps_src = t_steps
        this%src_step    = 1

        if (.not. allocated(this%E_ext_field)) allocate(this%E_ext_field(3, this%n_steps_src))

        open(newunit=funit, file=trim(external_src_file), status='old', action='read', iostat=ierr, iomsg=iomsg)

        if (ierr /= 0) then
            write(*,*) "Error opening external source file: '"//trim(external_src_file)//"'"
            write(*,*) "iostat = ", ierr, " iomsg = ", trim(iomsg)
            stop
        end if

        do i = 1, this%n_steps_src
            read(unit=funit, fmt=*, iostat=ierr) this%E_ext_field(1, i), this%E_ext_field(2, i), &
                                                 this%E_ext_field(3, i)
            if (ierr /= 0) then
                write(*,*) "Error reading external source file: "//trim(external_src_file)
                stop
            end if
        end do

        close(funit)

    end if

end subroutine init_TB

!###################################################################################################

subroutine kill_TB(this)

    class(TQ_sys_TB), intent(inout) :: this

    if (allocated(this%H0_mo))    deallocate(this%H0_mo)
    if (allocated(this%mu_x_mo))  deallocate(this%mu_x_mo)
    if (allocated(this%mu_y_mo))  deallocate(this%mu_y_mo)
    if (allocated(this%mu_z_mo))  deallocate(this%mu_z_mo)
    if (allocated(this%rho_mo))   deallocate(this%rho_mo)
    if (allocated(this%c_vecs))   deallocate(this%c_vecs)
    if (allocated(this%cT_vecs))  deallocate(this%cT_vecs)
    if (allocated(this%Ht))       deallocate(this%Ht)
    if (allocated(this%H_tot))    deallocate(this%H_tot)

    if (allocated(this%H0_at))    deallocate(this%H0_at)
    if (allocated(this%mu_x_at))  deallocate(this%mu_x_at)
    if (allocated(this%mu_y_at))  deallocate(this%mu_y_at)
    if (allocated(this%mu_z_at))  deallocate(this%mu_z_at)
    if (allocated(this%rho_at))   deallocate(this%rho_at)

    if (allocated(this%E_ext_field)) deallocate(this%E_ext_field)
    if (allocated(this%dipole))      deallocate(this%dipole)
    if (allocated(this%dip_old))     deallocate(this%dip_old)
    if (allocated(this%dPt_dt))      deallocate(this%dPt_dt)
    if (allocated(this%dPt_dt_old))  deallocate(this%dPt_dt_old)

end subroutine

!###################################################################################################

subroutine gs_calculate_TB(this)
    class(TQ_sys_TB), intent(inout) :: this

    complex(dp), allocatable :: aux_Im_mat(:,:)
    real(dp)   , allocatable :: aux_Re_mat(:,:)
    real(dp)   , allocatable :: e_vals(:)
    real(dp)   , allocatable :: work(:)
    real(dp)              :: work_query(1)
    integer               :: lwork, info
    integer               :: i, j

    allocate(aux_Re_mat(this%n_basis, this%n_basis))
    allocate(aux_Im_mat(this%n_basis, this%n_basis))

    if (.not. this%basis_mo) then

        if (.not. allocated(this%H0_at) .or. .not. allocated(this%c_vecs)) then
            write(*,*) "Error in gs_calculate_TB: H0_at or c_vecs not allocated."
            stop
        end if

        allocate(e_vals(this%n_basis))

        ! dsyev overwrites the input matrix with the eigenvectors, so work on a copy of H0_at.
        this%c_vecs = this%H0_at

        lwork = -1
        call dsyev('V', 'U', this%n_basis, this%c_vecs, this%n_basis, e_vals, work_query, lwork, info)

        lwork = max(1, int(work_query(1)))
        allocate(work(lwork))

        call dsyev('V', 'U', this%n_basis, this%c_vecs, this%n_basis, e_vals, work, lwork, info)

        if (info /= 0) then
            write(*,*) "Error diagonalizing H0_at with LAPACK dsyev. info = ", info
            stop
        end if

        do i = 1, this%n_basis
            this%H0_mo(i,i) = e_vals(i)
        end do

        this%cT_vecs = TRANSPOSE(this%c_vecs)

        deallocate(work)
        deallocate(e_vals)

        aux_Re_mat = MATMUL(this%mu_x_at, this%c_vecs)
        this%mu_x_mo = MATMUL(this%cT_vecs, aux_Re_mat)

        aux_Re_mat = MATMUL(this%mu_y_at, this%c_vecs)
        this%mu_y_mo = MATMUL(this%cT_vecs, aux_Re_mat)

        aux_Re_mat = MATMUL(this%mu_z_at, this%c_vecs)
        this%mu_z_mo = MATMUL(this%cT_vecs, aux_Re_mat)

    end if

    aux_Im_mat = MATMUL(this%rho_mo, this%H0_mo)

    this%E0 = M_ZERO

    do i = 1, this%n_basis
        this%E0 = this%E0 + REAL(aux_Im_mat(i,i))
    end do

    this%dipole     = M_ZERO

    aux_Im_mat = MATMUL(this%rho_mo, this%mu_x_mo)
    do i = 1, this%n_basis
        this%dipole(1) = this%dipole(1) + REAL(aux_Im_mat(i,i))
    end do

    aux_Im_mat = MATMUL(this%rho_mo, this%mu_y_mo)
    do i = 1, this%n_basis
        this%dipole(2) = this%dipole(2) + REAL(aux_Im_mat(i,i))
    end do

    aux_Im_mat = MATMUL(this%rho_mo, this%mu_z_mo)
    do i = 1, this%n_basis
        this%dipole(3) = this%dipole(3) + REAL(aux_Im_mat(i,i))
    end do

    this%dip_old    = this%dipole

    this%rho     = this%rho_mo
    this%rho_old = this%rho_mo

    deallocate(aux_Re_mat)
    deallocate(aux_Im_mat)

end subroutine gs_calculate_TB

!###################################################################################################

subroutine td_propagate_TB(this, tq_step, E_field_in)
    class(TQ_sys_TB), intent(inout) :: this
    integer           , intent(in)  :: tq_step
    real(dp)          , intent(in)  :: E_field_in(3)
    
    integer  :: i
    real(dp) :: E_field(3)
    
    complex(dp), allocatable :: commut_rhoH(:,:)
    complex(dp), allocatable :: aux_Im_mat(:,:)

    allocate(commut_rhoH(this%n_basis, this%n_basis))
    allocate(aux_Im_mat(this%n_basis, this%n_basis))
    
    E_field = E_field_in

    this%dPt_dt_old = this%dPt_dt
    
    if (this%external_src) then
        E_field = E_field + this%E_ext_field(:, this%src_step)
        this%src_step = this%src_step + 1
    end if
    
    this%Ht    = E_field(1) * this%mu_x_mo + E_field(2) * this%mu_y_mo + E_field(3) * this%mu_z_mo
    this%H_tot = this%H0_mo - this%Ht
    
    commut_rhoH = -2.0d0*Z_I*(MATMUL(this%H_tot, this%rho) - MATMUL(this%rho, this%H_tot))
    
    this%rho_new = this%rho_old + commut_rhoH * this%dt
    
    this%dipole = M_ZERO
    
    aux_Im_mat = MATMUL(this%rho_new, this%mu_x_mo)
    do i = 1, this%n_basis
        this%dipole(1) = this%dipole(1) + REAL(aux_Im_mat(i,i))
    end do
    
    aux_Im_mat = MATMUL(this%rho_new, this%mu_y_mo)
    do i = 1, this%n_basis
        this%dipole(2) = this%dipole(2) + REAL(aux_Im_mat(i,i))
    end do
    
    aux_Im_mat = MATMUL(this%rho_new, this%mu_z_mo)
    do i = 1, this%n_basis
        this%dipole(3) = this%dipole(3) + REAL(aux_Im_mat(i,i))
    end do
    
    this%Et = M_ZERO
    
    aux_Im_mat = MATMUL(this%rho_new, this%H_tot)
    do i = 1, this%n_basis
        this%Et = this%Et + REAL(aux_Im_mat(i,i))
    end do
    
    this%dPt_dt = (this%dipole - this%dip_old) / this%dt
    
    this%dip_old    = this%dipole
    this%rho_old    = this%rho
    this%rho        = this%rho_new

    deallocate(commut_rhoH)
    deallocate(aux_Im_mat)

end subroutine td_propagate_TB 

!###################################################################################################

subroutine write_output_TB(this, time, i_group, print_q_step)
    class(TQ_sys_TB), intent(in) :: this
    real(dp), intent(in) :: time
    integer, intent(in) :: i_group
    integer, intent(in) :: print_q_step

    character(len=20) :: dirname = "./output_q_group_"
    character(len=4)  :: group_number
    character(len=20) :: mol_number
    character(len=20) :: file_name = "occupations_mol_"
    character(len=99) :: full_dirname
    integer :: i
    integer :: ierr
    integer :: funit

    if (this%print_on .and. this%print_occ) then

        write(group_number,'(I4.4)') i_group
        write(mol_number,'(I7.7)') this%id

        full_dirname = trim(dirname)//trim(group_number)//"/"//trim(file_name)//trim(mol_number)//".dat"

        if (print_q_step == 0) then
            open(newunit=funit, file=full_dirname, status='replace', action="write", iostat=ierr)
            write(funit, *) "# Time (a.u.) ", "Occupations from 1 to N"
        else 
            open(newunit=funit, file=full_dirname, status='old', position='append', action="write", iostat=ierr)
        end if


        write(funit, *) time, (REAL(this%rho(i,i)), i = 1, this%n_basis)

        close(funit)

    end if

end subroutine write_output_TB

!###################################################################################################

end module q_sys_TB_mod