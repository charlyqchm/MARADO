module q_sys_dftb_mod

    use constants_mod
    use, intrinsic :: iso_fortran_env, only : output_unit
    use dftbp_common_constants, only : AA__Bohr, V_m__au, eV__Hartree, fs__au, imag
    use dftbplus
    ! Only needed for the internal test system
    !    use testhelpers, only : writeAutotestTag
    use q_sys_base_mod

    implicit none

    type, extends(TQ_sys_base) :: TQ_sys_dftb

        character(len=20)    :: temp_file
        logical              :: print_coordinates
        logical              :: init_velocities
        integer              :: dyn_type
        integer              :: n_atoms
        integer              :: n_at_typ
        integer              :: euler_steps
        integer              :: std_unit=100
        real(dp)             :: Kinetic_Energy 
        type(TDftbPlus)      :: dftbp
        type(TDftbPlusInput) :: dftbp_input

        character(len=2)    , allocatable  :: atom_names(:)
        real(dp)            , allocatable  :: coor(:,:)
        real(dp)            , allocatable  :: coor_old(:,:)
        real(dp)            , allocatable  :: coor_new(:,:)
        real(dp)            , allocatable  :: forces(:,:)
        real(dp)            , allocatable  :: vel(:,:)
        real(dp)            , allocatable  :: at_masses(:)
        real(dp)            , allocatable  :: at_charges(:)

        type(fnode), pointer :: pRoot
        type(fnode), pointer :: pGeo 
        type(fnode), pointer :: pHam
        type(fnode), pointer :: pDftb
        type(fnode), pointer :: pMaxAng
        type(fnode), pointer :: pSlakos
        type(fnode), pointer :: pType2Files
        type(fnode), pointer :: pElecDyn
        type(fnode), pointer :: pElecStatic
        type(fnode), pointer :: pExternal
        type(fnode), pointer :: pPerturb
        type(fnode), pointer :: pLaser
        type(fnode), pointer :: pAnalysis

        contains
            procedure :: init => init_dftb
            procedure :: kill => kill_dftb
            procedure :: gs_calculate => gs_calculate_dftb
            procedure :: td_propagate => td_propagate_dftb
            procedure :: write_output => write_output_dftb
    end type TQ_sys_dftb

    !Node of a parsed DFTB+-like (HSD) input. The nodes are stored in the order they appear in
    !the file, and "depth" gives the nesting level (0 = top level).
    type :: hsd_node
        character(len=64)  :: name     = ""
        character(len=256) :: value    = ""
        integer            :: depth    = 0
        logical            :: is_block = .false.
    end type hsd_node

    type :: fnode_ptr
        type(fnode), pointer :: p => null()
    end type fnode_ptr

contains

!###################################################################################################

subroutine init_dftb(this, id, id_file, dt, t_steps, rank, print_on)

    class(TQ_sys_dftb), intent(inout) :: this
    integer           , intent(in)    :: id
    integer           , intent(in)    :: id_file
    integer           , intent(in)    :: t_steps
    integer           , intent(in)    :: rank
    real(dp)          , intent(in)    :: dt
    logical           , intent(in)    :: print_on

    character(len=20)  :: file_name = "molecule_"
    character(len=20)  :: file_exten = ".in"
    character(len=20)  :: file_number
    character(len=20)  :: input_name
    character(len=30)  :: print_coord_ch
    character(len=100) :: dyn_type_str
    character(len=100) :: external_src_file
    integer            :: n_atoms
    integer            :: n_at_typ
    integer            :: euler_steps
    integer            :: i, ii, j
    integer            :: ierr
    integer            :: funit
    integer            :: ic
    character(len=200) :: iomsg
    logical            :: atom_type_exists
    logical            :: periodic = .false. !This could be read from the file in the future
    logical            :: ion_dyn = .false.
    logical            :: scc = .true. !This could be read from the file in the future
    real(dp)           :: scc_tol
    integer            :: n_hsd
    integer            :: i_hsd

    type(hsd_node)    , allocatable :: hsd_nodes(:)
    character(len = 2), allocatable :: atom_type_list(:)
    character(len = 2), allocatable :: max_ang_orb(:)
    real(dp)          , allocatable :: atomNetCarges(:,:)
    real(dp)          , allocatable :: coor_aux(:,:)
    integer           , allocatable :: atom_type(:)

    this%id       = id
    this%id_file  = id_file
    this%dt       = dt
    this%t_steps  = t_steps
    this%rank     = rank
    this%print_on = print_on

    this%std_unit = this%std_unit + this%rank

    write(file_number,'(I7.7)') id_file
    input_name = trim(file_name)//trim(file_number)//trim(file_exten)

    write(file_number,'(I7.7)') this%rank
    this%temp_file = "temp_output_"//trim(file_number)

    open (action='read', file=input_name, iostat=ierr, newunit=funit)

    if (ierr /= 0) then
        write(*,*) "Error: Could not open molecule file ", trim(input_name)
        stop
    end if

    read (unit=funit, fmt=*, iostat=ierr) n_atoms, n_at_typ
    read (unit=funit, fmt=*, iostat=ierr) dyn_type_str, print_coord_ch, euler_steps
    
    this%print_coordinates = .false.

    if (this%print_on) then
        
        if (ierr /= 0) then
            write(*,*) "Error: Could not read dynamics type, euler steps and"
            write(*,*) "print coordinates flag from file ", trim(input_name)
            stop
        end if

        if (print_coord_ch == "print_coordinates_on") then
            this%print_coordinates = .true.
        else if (print_coord_ch == "print_coordinates_off") then
            this%print_coordinates = .false.
        else
            write(*,*) "Error: Unknown print coordinates flag ", trim(print_coord_ch)
            stop
        end if
    end if

    this%n_atoms = n_atoms
    this%n_at_typ = n_at_typ
    this%euler_steps = euler_steps

    select case (trim(dyn_type_str))
    case ('electrons')
        this%dyn_type = DFTB_ELEC_DYN
        if (this%print_coordinates) then
            write(*,*) "Warning: print_coordinates flag is on, "
            write(*,*) "but it will be ignored for electron dynamics."
            this%print_coordinates = .false.
        end if
    case ('ehrenfest')
        this%dyn_type = DFTB_EHREN_DYN
    case ('born-oppenheimer')
        this%dyn_type = DFTB_BO_DYN
    case default
        write(*,*) "Error: Unknown DFTB dynamics type ", trim(dyn_type_str)
        stop
    end select

    this%init_velocities = .false.
    if (this%dyn_type == DFTB_BO_DYN .or. this%dyn_type == DFTB_EHREN_DYN) then
        this%init_velocities = .true.
    end if

    if (.not. allocated(this%atom_names)) allocate(this%atom_names(n_atoms))
    if (.not. allocated(this%coor))       allocate(this%coor(3, n_atoms))
    if (.not. allocated(this%at_charges)) allocate(this%at_charges(n_atoms))
    if (.not. allocated(this%dipole))     allocate(this%dipole(3))
    if (.not. allocated(this%dip_old))    allocate(this%dip_old(3))
    if (.not. allocated(this%dPt_dt))     allocate(this%dPt_dt(3))
    if (.not. allocated(this%dPt_dt_old)) allocate(this%dPt_dt_old(3))

    this%coor       = M_ZERO
    this%at_charges = M_ZERO
    this%dipole     = M_ZERO
    this%dip_old    = M_ZERO
    this%dPt_dt     = M_ZERO
    this%dPt_dt_old = M_ZERO    

    if (.not. allocated(atom_type_list))  allocate(atom_type_list(n_at_typ))
    if (.not. allocated(max_ang_orb))     allocate(max_ang_orb(n_at_typ))
    if (.not. allocated(atom_type))       allocate(atom_type(n_atoms))
    if (.not. allocated(coor_aux))        allocate(coor_aux(3, n_atoms))

    if (this%dyn_type == DFTB_BO_DYN .or. this%dyn_type == DFTB_EHREN_DYN) then
        if (this%init_velocities) then
            if (.not. allocated(this%vel)) allocate(this%vel(3, n_atoms))
            this%vel = M_ZERO
        end if
    end if

    if (this%dyn_type == DFTB_BO_DYN) then
        if (.not. allocated(this%coor_old))   allocate(this%coor_old(3, n_atoms))
        if (.not. allocated(this%coor_new))   allocate(this%coor_new(3, n_atoms))
        if (.not. allocated(this%forces))     allocate(this%forces(3, n_atoms))
        if (.not. allocated(this%at_masses))  allocate(this%at_masses(n_atoms))

        this%at_masses  = M_ZERO
        this%forces     = M_ZERO
        this%coor_old   = M_ZERO
        this%coor_new   = M_ZERO

    end if 

    read (unit=funit, fmt=*, iostat=ierr) (atom_type_list(i), i=1,n_at_typ)
    if (ierr /= 0) then
        write(*,*) "Error: Could not read atom type list from file ", trim(input_name)
        stop
    end if

    read (unit=funit, fmt=*, iostat=ierr) (max_ang_orb(i), i=1,n_at_typ)
    if (ierr /= 0) then
        write(*,*) "Error: Could not read max angular momenta from file ", trim(input_name)
        stop
    end if

    read (unit=funit, fmt=*, iostat=ierr) scc_tol
    if (ierr /= 0) then
        write(*,*) "Error: Could not read SCC tolerance from file ", trim(input_name)
        stop
    end if

    read (unit=funit, fmt=*, iostat=ierr) 

    do i=1, n_atoms

        if (this%init_velocities) then

            read (unit=funit, fmt=*, iostat=ierr) this%atom_names(i), &
                                                this%coor(1,i), &
                                                this%coor(2,i), &
                                                this%coor(3,i), &
                                                this%vel(1,i), &
                                                this%vel(2,i), &
                                                this%vel(3,i)
            if (ierr /= 0) then
                write(*,*) "Error: Could not read atom name, coordinates and velocity for &
                           &atom", i, "from file ", trim(input_name)
                write(*,*) "Check that the file has 7 columns per atom line &
                           &(name, x, y, z, vx, vy, vz) for this dynamics type."
                stop
            end if
            this%vel(:, i) = this%vel(:, i) * AA_to_au / ps_to_au
        else

            read (unit=funit, fmt=*, iostat=ierr) this%atom_names(i), &
                                                this%coor(1,i), &
                                                this%coor(2,i), &
                                                this%coor(3,i)
            if (ierr /= 0) then
                write(*,*) "Error: Could not read atom name and coordinates for atom", i, &
                           "from file ", trim(input_name)
                write(*,*) "Check that the file has 4 columns per atom line &
                           &(name, x, y, z) for this dynamics type."
                stop
            end if

        end if

        atom_type_exists = .false.
        do ii=1, n_at_typ
            if (this%atom_names(i) == atom_type_list(ii)) then
                atom_type(i) = ii
                atom_type_exists = .true.
            end if
        end do
        if (.not. atom_type_exists) then
            write(*,*) "ERROR.", this%atom_names(i), "is not an atom type."
            stop
        end if
    end do

    !This file should save data printed by DFTB+. The program will overwrite it
    !in different stages to not accumulate too much data.

    
    open(newunit=this%std_unit, file=trim(this%temp_file), status='replace', action="write")

    call TDftbPlus_init(this%dftbp, outputUnit=this%std_unit)

    call this%dftbp%getEmptyInput(this%dftbp_input)
    call this%dftbp_input%getRootNode(this%pRoot)


    call setChild(this%pRoot, "Geometry", this%pGeo)
    call setChildValue(this%pGeo, "Periodic", periodic)
    call setChildValue(this%pGeo, "TypeNames", atom_type_list)

    coor_aux = M_ZERO
    call setChildValue(this%pGeo, "TypesAndCoordinates", reshape(atom_type, &
                       [1, size(atom_type)]), coor_aux)

    call setChild(this%pRoot, "Hamiltonian", this%pHam)
    call setChild(this%pHam, "Dftb", this%pDftb)
    call setChildValue(this%pDftb, "Scc", scc)
    call setChildValue(this%pDftb, "SccTolerance", scc_tol)

    call setChild(this%pDftb, "MaxAngularMomentum", this%pMaxAng)
    do j=1, n_at_typ
        call setChildValue(this%pMaxAng, atom_type_list(j), max_ang_orb(j))
    end do

    !Extra DFTB+ options (HubbardDerivs, Dispersion, ThirdOrderFull, ...) are read from the
    !DFTB+-like block "Hamiltonian = DFTB { ... }" of the molecule file.
    call hsd_read_file(input_name, hsd_nodes, n_hsd)
    i_hsd = hsd_find_node(hsd_nodes, n_hsd, "Hamiltonian/Dftb")
    if (i_hsd > 0) call hsd_apply_children(hsd_nodes, n_hsd, i_hsd, this%pDftb)

    call setChild(this%pDftb, "SlaterKosterFiles", this%pSlakos)
    call setChild(this%pSlakos, "Type2FileNames", this%pType2Files)
    call setChildValue(this%pType2Files, "Prefix", "./")
    call setChildValue(this%pType2Files, "Separator", "-")
    call setChildValue(this%pType2Files, "Suffix", ".skf")

    if (this%dyn_type == DFTB_BO_DYN) then
        call setChild(this%pRoot, "Analysis", this%pAnalysis)
        call setChildValue(this%pAnalysis, "PrintForces", .true.)
        
        call setChild(this%pDftb, "ElectricField", this%pElecStatic)
        call setChild(this%pElecStatic, "External", this%pExternal)
        call setChildValue(this%pExternal, "Strength", 0.0_dp)
        call setChildValue(this%pExternal, "Direction", [1.0_dp, 1.0_dp, 1.0_dp])
    else
        !  set up electron dynamics options
        call setChild(this%pRoot, "ElectronDynamics", this%pElecDyn)
        call setChildValue(this%pElecDyn, "Steps", this%t_steps)
        call setChildValue(this%pElecDyn, "TimeStep", this%dt)
        call setChildValue(this%pElecDyn, "FieldStrength", 1.0_dp)
     
        if (this%dyn_type == DFTB_EHREN_DYN) ion_dyn = .true.
     
        call setChildValue(this%pElecDyn, "IonDynamics", ion_dyn)
        call setChildValue(this%pElecDyn, "InitialTemperature", 0.0_dp)
        call setChildValue(this%pElecDyn, "EulerFrequency", this%euler_steps)
        call setChildValue(this%pElecDyn, "VerboseDynamics", .false.)

        if (this%init_velocities) then
            call setChildValue(this%pElecDyn, "Velocities", this%vel)
        end if

        call setChild(this%pElecDyn, "Perturbation", this%pPerturb)
        call setChild(this%pPerturb, "Laser", this%pLaser)
    ! these twovalues will be overriden
        call setChildValue(this%pLaser, "PolarisationDirection", [1.0_dp, 1.0_dp, 1.0_dp])
        call setChildValue(this%pLaser, "LaserEnergy", 1.0_dp)
    end if

    call this%dftbp%setupCalculator(this%dftbp_input)


!Checking the presence of the name of the external source file
    read(unit=funit, fmt="(A)", iostat=ierr) external_src_file

    !Strip stray CR (Windows line endings) and other control characters that
    !survive trim() and would otherwise make the filename invalid
    do ic = 1, len(external_src_file)
        if (iachar(external_src_file(ic:ic)) < 32) external_src_file(ic:ic) = ' '
    end do
    external_src_file = adjustl(external_src_file)

    !The DFTB+-like block can follow directly, in which case there is no external source
    if (index(external_src_file, '=') > 0 .or. index(external_src_file, '{') > 0) then
        external_src_file = ' '
    end if

    if (ierr /= 0 .or. len_trim(external_src_file) == 0) then
        this%external_src = .false.
    else
        this%external_src = .true.

        this%n_steps_src = t_steps

        this%src_step = 1

        if (.not. allocated(this%E_ext_field)) allocate(this%E_ext_field(3, this%n_steps_src))

    end if

    close(funit)

!Trying to read the external source file into the E_ext_field array

    if (this%external_src) then

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

    if (allocated(atom_type_list)) deallocate(atom_type_list)
    if (allocated(max_ang_orb))    deallocate(max_ang_orb)
    if (allocated(atom_type))      deallocate(atom_type)
    if (allocated(atomNetCarges))  deallocate(atomNetCarges)

end subroutine init_dftb

!###################################################################################################

subroutine kill_dftb(this)

    class(TQ_sys_dftb), intent(inout) :: this

    call TDftbPlus_destruct(this%dftbp)

    if (allocated(this%atom_names)) deallocate(this%atom_names)
    if (allocated(this%coor))       deallocate(this%coor)
    if (allocated(this%coor_old))   deallocate(this%coor_old)
    if (allocated(this%coor_new))   deallocate(this%coor_new)
    if (allocated(this%forces))     deallocate(this%forces)
    if (allocated(this%vel))        deallocate(this%vel)
    if (allocated(this%at_masses))  deallocate(this%at_masses)
    if (allocated(this%at_charges)) deallocate(this%at_charges)
    if (allocated(this%dipole))     deallocate(this%dipole)
    if (allocated(this%dip_old))    deallocate(this%dip_old)
    if (allocated(this%dPt_dt))     deallocate(this%dPt_dt)
    if (allocated(this%dPt_dt_old)) deallocate(this%dPt_dt_old)

    close(this%std_unit)

end subroutine kill_dftb

!###################################################################################################

subroutine gs_calculate_dftb(this)

    class(TQ_sys_dftb), intent(inout) :: this

    integer  :: i, kk

    this%E0      = M_ZERO
    this%Et      = M_ZERO
    this%dipole  = M_ZERO
    this%dip_old = M_ZERO

    do kk = 1, this%n_atoms
        this%coor(:, kk) = this%coor(:, kk) * AA__Bohr
    end do

    call this%dftbp%setGeometry(this%coor)
    call this%dftbp%getEnergy(this%E0)

    call this%dftbp%getGrossCharges(this%at_charges)

    do i = 1, this%n_atoms
        this%dipole(1) = this%dipole(1) + this%at_charges(i) * this%coor(1,i)
        this%dipole(2) = this%dipole(2) + this%at_charges(i) * this%coor(2,i)
        this%dipole(3) = this%dipole(3) + this%at_charges(i) * this%coor(3,i)
    end do

    if (this%dyn_type == DFTB_BO_DYN) then
        call this%dftbp%getAtomicMasses(this%at_masses)
        this%Kinetic_Energy = M_ZERO
    else
    !This nitialize time propagation for Ehrenfest dynamics as well
    !TO-DO: maybe this should be somewhere else.
        call this%dftbp%initializeTimeProp(this%dt, .true., .false.)
    end if

    if (allocated(this%coor_old)) this%coor_old = this%coor
    this%dip_old = this%dipole

end subroutine gs_calculate_dftb

!###################################################################################################

subroutine td_propagate_dftb(this, tq_step, E_field_in)

    class(TQ_sys_dftb), intent(inout) :: this
    integer           , intent(in)    :: tq_step
    real(dp)          , intent(in)    :: E_field_in(3)
    
    integer  :: i
    real(dp) :: E_amp
    real(dp) :: vv
    real(dp) :: E_field(3)
    real(dp), allocatable :: aux_at_charges(:,:)
    real(dp), allocatable :: dip_aux(:,:)


    if(.not. allocated(aux_at_charges)) allocate(aux_at_charges(this%n_atoms, 1))
    if(.not. allocated(dip_aux)) allocate(dip_aux(3, 1))

    
    E_field = E_field_in
    
    if (this%external_src) then
        E_field = E_field + this%E_ext_field(:, this%src_step)
        this%src_step = this%src_step + 1
    end if

    E_amp = dsqrt(E_field(1)**2 + E_field(2)**2 + E_field(3)**2)

    this%dPt_dt_old = this%dPt_dt

    if (this%dyn_type == DFTB_BO_DYN) then
        this%Kinetic_Energy = M_ZERO

        if (E_amp /= M_ZERO) E_field = E_field / E_amp
        
        call this%dftbp%setGeometry(this%coor)
        call this%dftbp%setExternalEfield(E_amp, E_field)
        call this%dftbp%getGradients(this%forces) !The force has the opposite sign of the gradient.
        ! call this%dftbp%getEnergy(this%Et)
        call this%dftbp%getGrossCharges(this%at_charges)

        this%dip_old = this%dipole
        this%dipole  = M_ZERO

        do i = 1, this%n_atoms

            if (this%init_velocities .and. tq_step == 1 ) then

                this%coor_new(:, i) = this%coor(:, i) + this%vel(:, i) * this%dt - &
                                    0.5*(this%forces(:, i) / this%at_masses(i)) * this%dt**2
            
            else
                
                this%coor_new(:, i) = 2.0d0*this%coor(:, i) - this%coor_old(:, i) - &
                                    (this%forces(:, i) / this%at_masses(i)) * this%dt**2
            end if

            this%dipole(:) = this%dipole(:) + this%at_charges(i) * this%coor(:,i)

            this%vel(:, i) = (this%coor_new(:, i) - this%coor_old(:, i)) / (2.0d0 * this%dt)
            vv = this%vel(1,i)**2 + this%vel(2,i)**2 + this%vel(3,i)**2

            this%Kinetic_Energy = this%Kinetic_Energy + 0.5d0 * this%at_masses(i) * vv

        end do

        this%Et = this%Kinetic_Energy

        this%coor_old = this%coor
        this%coor     = this%coor_new

        this%dPt_dt = (this%dipole - this%dip_old) / this%dt

    else ! electron dynamics and Ehrenfest dynamics.

        this%dip_old = this%dipole
        call this%dftbp%setTdElectricField(E_field)

        call this%dftbp%doOneTdStep(tq_step, dipole=dip_aux, energy=this%Et,  &
                                   atomNetCharges=aux_at_charges, coord=this%coor, &
                                   force=this%forces)

        this%dipole     = dip_aux(:,1)
        this%at_charges = aux_at_charges(:,1)

        this%dPt_dt = (this%dipole - this%dip_old) / this%dt

    end if

    if (allocated(aux_at_charges)) deallocate(aux_at_charges)
    if (allocated(dip_aux))        deallocate(dip_aux)

end subroutine td_propagate_dftb

!###################################################################################################

subroutine write_output_dftb(this, time, i_group, print_q_step)

    class(TQ_sys_dftb), intent(in) :: this
    real(dp), intent(in) :: time
    integer, intent(in) :: i_group
    integer, intent(in) :: print_q_step

    character(len=20) :: dirname = "./output_q_group_"
    character(len=4)  :: group_number
    character(len=20) :: mol_number
    character(len=20) :: file_name = "molecule_"
    character(len=99) :: full_dirname
    integer :: i
    integer :: n
    integer :: ierr
    integer :: funit

    if (this%print_on .and. this%print_coordinates) then

        write(group_number,'(I4.4)') i_group
        write(mol_number,'(I7.7)') this%id

        full_dirname = trim(dirname)//trim(group_number)//"/"//trim(file_name)//trim(mol_number)//".xyz"

        if (print_q_step == 0) then
            open(newunit=funit, file=full_dirname, status='replace', action="write", iostat=ierr)
        else 
            open(newunit=funit, file=full_dirname, status='old', position='append', action="write", iostat=ierr)
        end if

        write(funit, *) this%n_atoms
        write(funit, '("# Time = ", ES18.8, " (a.u.)")') time   

        do i = 1, this%n_atoms
            write(funit, *) this%atom_names(i), &
                            this%coor(1,i)/AA__Bohr, &
                            this%coor(2,i)/AA__Bohr, &
                            this%coor(3,i)/AA__Bohr
        end do

        close(funit)

    end if

end subroutine write_output_dftb

!###################################################################################################
! Reader of a DFTB+-like (HSD) input block inside the molecule_*.in file.
!
! Supported syntax (one item per line, '#' starts a comment):
!
!   Hamiltonian = DFTB {
!       ThirdOrderFull = Yes
!       HubbardDerivs {
!           C = -0.1492
!       }
!       Dispersion = DftD3 {
!           Damping = BeckeJohnson {
!               a1 = 0.5719
!           }
!       }
!   }
!
! "Name = Method {" is equivalent to "Name { Method {". Not supported: unit modifiers
! ("Name [Angstrom] = ..."), data blocks ("Name = {" ... "}"), and braces sharing a line with values.
! The block starts at the first line containing '=' or '{', so it must come after the
! legacy lines (and after the external source line, which may be left empty).
!###################################################################################################

subroutine hsd_read_file(file_name, nodes, n_nodes)

    character(len=*)           , intent(in)  :: file_name
    type(hsd_node), allocatable, intent(out) :: nodes(:)
    integer                    , intent(out) :: n_nodes

    integer, parameter :: MAX_OPEN = 64

    character(len=512) :: line
    character(len=512) :: head
    character(len=64)  :: name
    character(len=64)  :: method
    integer            :: funit, ierr
    integer            :: i_line, i_eq, n_line
    integer            :: depth, n_open, n_levels
    integer            :: levels(MAX_OPEN)
    logical            :: started

    allocate(nodes(32))
    n_nodes = 0
    depth   = 0
    n_open  = 0
    i_line  = 0
    started = .false.

    open(newunit=funit, file=file_name, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write(*,*) "Error: Could not open molecule file ", trim(file_name)
        stop
    end if

    do
        read(funit, '(A)', iostat=ierr) line
        if (ierr /= 0) exit
        i_line = i_line + 1

        call hsd_clean_line(line)
        if (len_trim(line) == 0) cycle

        if (.not. started) then
            if (index(line, '=') == 0 .and. index(line, '{') == 0) cycle
            started = .true.
        end if

        do while (line(1:1) == '}')
            if (n_open == 0) call hsd_error(file_name, i_line, "unmatched '}'")
            depth  = depth - levels(n_open)
            n_open = n_open - 1
            line   = adjustl(line(2:))
        end do
        if (len_trim(line) == 0) cycle

        n_line = len_trim(line)

        if (line(n_line:n_line) == '{') then

            head = line(:n_line-1)
            i_eq = index(head, '=')

            if (i_eq > 0) then
                name   = adjustl(head(:i_eq-1))
                method = adjustl(head(i_eq+1:))
                if (len_trim(method) == 0) call hsd_error(file_name, i_line, &
                    "data blocks 'Name = {' are not supported")
                call hsd_check_name(file_name, i_line, name)
                call hsd_check_name(file_name, i_line, method)
                call hsd_add_node(nodes, n_nodes, name,   "", depth,     .true.)
                call hsd_add_node(nodes, n_nodes, method, "", depth + 1, .true.)
                n_levels = 2
            else
                name = adjustl(head)
                call hsd_check_name(file_name, i_line, name)
                call hsd_add_node(nodes, n_nodes, name, "", depth, .true.)
                n_levels = 1
            end if

            if (n_open == MAX_OPEN) call hsd_error(file_name, i_line, "blocks nested too deeply")
            n_open = n_open + 1
            levels(n_open) = n_levels
            depth = depth + n_levels

        else

            i_eq = index(line, '=')
            if (i_eq == 0) call hsd_error(file_name, i_line, "expected 'Name = value' or 'Name {'")

            name = adjustl(line(:i_eq-1))
            call hsd_check_name(file_name, i_line, name)
            if (len_trim(line(i_eq+1:)) == 0) call hsd_error(file_name, i_line, "empty value")
            call hsd_add_node(nodes, n_nodes, name, adjustl(line(i_eq+1:)), depth, .false.)

        end if
    end do

    close(funit)

    if (n_open /= 0) call hsd_error(file_name, i_line, "missing '}' at the end of the file")

end subroutine hsd_read_file

!###################################################################################################

subroutine hsd_error(file_name, i_line, msg)

    character(len=*), intent(in) :: file_name
    integer         , intent(in) :: i_line
    character(len=*), intent(in) :: msg

    write(*,*) "Error reading DFTB+-like input from ", trim(file_name), " (line", i_line, "): ", msg
    stop

end subroutine hsd_error

!###################################################################################################

subroutine hsd_check_name(file_name, i_line, name)

    character(len=*), intent(in) :: file_name
    integer         , intent(in) :: i_line
    character(len=*), intent(in) :: name

    if (len_trim(name) == 0) call hsd_error(file_name, i_line, "missing name")
    if (index(trim(name), ' ') > 0) call hsd_error(file_name, i_line, "invalid name '"//trim(name)//"'")

end subroutine hsd_check_name

!###################################################################################################

subroutine hsd_clean_line(line)

    character(len=*), intent(inout) :: line

    integer :: i

    do i = 1, len(line)
        if (line(i:i) == char(9) .or. line(i:i) == char(13)) line(i:i) = ' '
    end do

    i = index(line, '#')
    if (i > 0) line(i:) = ' '

    line = adjustl(line)

end subroutine hsd_clean_line

!###################################################################################################

subroutine hsd_add_node(nodes, n_nodes, name, value, depth, is_block)

    type(hsd_node), allocatable, intent(inout) :: nodes(:)
    integer                    , intent(inout) :: n_nodes
    character(len=*)           , intent(in)    :: name
    character(len=*)           , intent(in)    :: value
    integer                    , intent(in)    :: depth
    logical                    , intent(in)    :: is_block

    type(hsd_node), allocatable :: tmp(:)

    if (n_nodes == size(nodes)) then
        allocate(tmp(2*size(nodes)))
        tmp(1:n_nodes) = nodes(1:n_nodes)
        call move_alloc(tmp, nodes)
    end if

    n_nodes = n_nodes + 1
    nodes(n_nodes)%name     = trim(name)
    nodes(n_nodes)%value    = trim(value)
    nodes(n_nodes)%depth    = depth
    nodes(n_nodes)%is_block = is_block

end subroutine hsd_add_node

!###################################################################################################

function hsd_to_lower(str) result(low)

    character(len=*), intent(in) :: str
    character(len=len(str))      :: low

    integer :: i, ic

    low = str
    do i = 1, len(str)
        ic = iachar(str(i:i))
        if (ic >= iachar('A') .and. ic <= iachar('Z')) low(i:i) = achar(ic + 32)
    end do

end function hsd_to_lower

!###################################################################################################

!Returns the index of the block given by a "/"-separated path (case insensitive), or 0 if absent.
function hsd_find_node(nodes, n_nodes, path) result(idx)

    type(hsd_node), intent(in) :: nodes(:)
    integer       , intent(in) :: n_nodes
    character(len=*), intent(in) :: path
    integer                    :: idx

    character(len=len(path)) :: rest
    character(len=64)        :: comp
    integer :: i_sl, j, k, lo, hi, d

    idx  = 0
    lo   = 1
    hi   = n_nodes
    d    = 0
    rest = adjustl(path)

    do while (len_trim(rest) > 0)

        i_sl = index(rest, '/')
        if (i_sl == 0) then
            comp = rest
            rest = ' '
        else
            comp = rest(:i_sl-1)
            rest = adjustl(rest(i_sl+1:))
        end if

        idx = 0
        do j = lo, hi
            if (nodes(j)%depth == d .and. nodes(j)%is_block .and. &
                hsd_to_lower(trim(nodes(j)%name)) == hsd_to_lower(trim(comp))) then
                idx = j
                exit
            end if
        end do
        if (idx == 0) return

        k = idx + 1
        do while (k <= n_nodes)
            if (nodes(k)%depth <= d) exit
            k = k + 1
        end do
        lo = idx + 1
        hi = k - 1
        d  = d + 1

    end do

end function hsd_find_node

!###################################################################################################

!Creates in the DFTB+ tree "target" all the descendants of the block nodes(idx).
subroutine hsd_apply_children(nodes, n_nodes, idx, target)

    type(hsd_node), intent(in) :: nodes(:)
    integer       , intent(in) :: n_nodes
    integer       , intent(in) :: idx
    type(fnode), pointer       :: target

    type(fnode_ptr), allocatable :: ptrs(:)
    integer :: j, d, d0

    d0 = nodes(idx)%depth
    allocate(ptrs(0:maxval(nodes(1:n_nodes)%depth) + 1))
    ptrs(d0)%p => target

    do j = idx + 1, n_nodes
        d = nodes(j)%depth
        if (d <= d0) exit

        if (nodes(j)%is_block) then
            call setChild(ptrs(d-1)%p, trim(nodes(j)%name), ptrs(d)%p)
        else
            call hsd_set_value(ptrs(d-1)%p, trim(nodes(j)%name), trim(nodes(j)%value))
        end if
    end do

    deallocate(ptrs)

end subroutine hsd_apply_children

!###################################################################################################

!Sets "value" as a logical, integer, real, numeric array or string, according to its content.
subroutine hsd_set_value(parent, name, value)

    type(fnode), pointer         :: parent
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value

    character(len=256), allocatable :: tok(:)
    real(dp)          , allocatable :: r_arr(:)
    integer           , allocatable :: i_arr(:)
    character(len=256) :: low
    integer :: n_tok, i
    logical :: all_int, all_real

    if (value(1:1) == '"' .or. value(1:1) == "'") then
        call setChildValue(parent, name, value(2:len(value)-1))
        return
    end if

    call hsd_split_tokens(value, tok, n_tok)

    all_int  = .true.
    all_real = .true.
    do i = 1, n_tok
        if (.not. hsd_is_integer(tok(i))) all_int = .false.
        if (.not. hsd_is_real(tok(i)))    all_real = .false.
    end do

    if (n_tok == 1) then
        low = hsd_to_lower(tok(1))
        if (trim(low) == 'yes' .or. trim(low) == 'true') then
            call setChildValue(parent, name, .true.)
        else if (trim(low) == 'no' .or. trim(low) == 'false') then
            call setChildValue(parent, name, .false.)
        else if (all_int) then
            allocate(i_arr(1))
            read(tok(1), *) i_arr(1)
            call setChildValue(parent, name, i_arr(1))
        else if (all_real) then
            allocate(r_arr(1))
            read(tok(1), *) r_arr(1)
            call setChildValue(parent, name, r_arr(1))
        else
            call setChildValue(parent, name, trim(value))
        end if
    else if (all_int) then
        allocate(i_arr(n_tok))
        read(value, *) i_arr
        call setChildValue(parent, name, i_arr)
    else if (all_real) then
        allocate(r_arr(n_tok))
        read(value, *) r_arr
        call setChildValue(parent, name, r_arr)
    else
        call setChildValue(parent, name, trim(value))
    end if

end subroutine hsd_set_value

!###################################################################################################

subroutine hsd_split_tokens(str, tok, n_tok)

    character(len=*)               , intent(in)  :: str
    character(len=256), allocatable, intent(out) :: tok(:)
    integer                        , intent(out) :: n_tok

    integer :: i, i0, n
    logical :: in_tok

    n = len_trim(str)

    n_tok  = 0
    in_tok = .false.
    do i = 1, n
        if (str(i:i) == ' ' .or. str(i:i) == ',') then
            in_tok = .false.
        else if (.not. in_tok) then
            in_tok = .true.
            n_tok  = n_tok + 1
        end if
    end do

    allocate(tok(max(n_tok, 1)))
    tok = ' '

    n_tok  = 0
    in_tok = .false.
    i0     = 1
    do i = 1, n + 1
        if (i > n) then
            if (in_tok) tok(n_tok) = str(i0:n)
        else if (str(i:i) == ' ' .or. str(i:i) == ',') then
            if (in_tok) tok(n_tok) = str(i0:i-1)
            in_tok = .false.
        else if (.not. in_tok) then
            in_tok = .true.
            n_tok  = n_tok + 1
            i0     = i
        end if
    end do

end subroutine hsd_split_tokens

!###################################################################################################

function hsd_is_integer(tok) result(ok)

    character(len=*), intent(in) :: tok
    logical                      :: ok

    integer :: n, i0

    n  = len_trim(tok)
    i0 = 1
    if (n > 0) then
        if (tok(1:1) == '+' .or. tok(1:1) == '-') i0 = 2
    end if

    ok = n >= i0
    if (ok) ok = verify(tok(i0:n), '0123456789') == 0

end function hsd_is_integer

!###################################################################################################

function hsd_is_real(tok) result(ok)

    character(len=*), intent(in) :: tok
    logical                      :: ok

    real(dp) :: val
    integer  :: n, ierr

    n  = len_trim(tok)
    ok = n > 0
    if (ok) ok = verify(tok(:n), '0123456789+-.eEdD') == 0
    if (ok) ok = scan(tok(:n), '0123456789') > 0
    if (ok) then
        read(tok(:n), *, iostat=ierr) val
        ok = ierr == 0
    end if

end function hsd_is_real

!###################################################################################################
end module q_sys_dftb_mod