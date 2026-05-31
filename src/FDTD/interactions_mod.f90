module interactions_mod

#ifdef USE_MPI
    use mpi
#endif

    use constants_mod
    use sources_subs_mod
    use mxll_base_mod
    use mxll_1D_mod
    use mxll_2D_mod
    use mxll_3D_mod
    use q_group_mod

    implicit none

contains

!###################################################################################################

subroutine point_source_interactions(mxll, sources)
    class(TMxll)       ,intent(inout) :: mxll
    type(TSources_list),intent(inout) :: sources

    integer  :: i, j, k, s
    integer  :: i_id, j_id, k_id
    integer  :: n_ker
    real(dp) :: J_av
    real(dp) :: c_src

    c_src = mxll%dt_eps0/mxll%dr/c0/2.0d0

    select type(mxll)
    class is(TMxll_1D)

#ifdef USE_MPI
        if (.not. allocated(mxll%Ex)) return
#endif

        do s=1, sources%n_p_src
            n_ker = sources%points(s)%n_ker
            select case (sources%points(s)%polarization)
            case ('x')
                do i=-n_ker, n_ker
                    i_id = sources%points(s)%ind_i(i,1,1)
                    J_av = sources%points(s)%J_mat(i,1,1)
                    mxll%Ex(i_id) = mxll%Ex(i_id) + c_src * J_av
                end do
            case default
                print *, "Error: Unknown source direction in 1D mxll."
                stop
            end select
        end do

    class is(TMxll_2D)
    
        do s=1, sources%n_p_src
            n_ker = sources%points(s)%n_ker
            select case (sources%points(s)%polarization)
            case ('x')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                    if (sources%points(s)%in_this_rank(i,j,1)) then
                        i_id = sources%points(s)%ind_i(i,j,1)
                        j_id = sources%points(s)%ind_j(i,j,1)
                        J_av = 0.5d0*(sources%points(s)%J_mat(i,j,1) + &
                                        sources%points(s)%J_mat(i+1,j,1))
                        mxll%Ex(i_id,j_id) = mxll%Ex(i_id,j_id) + c_src * J_av
                    end if
                end do
                end do
            case ('y')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                    if (sources%points(s)%in_this_rank(i,j,1)) then
                        i_id = sources%points(s)%ind_i(i,j,1)
                        j_id = sources%points(s)%ind_j(i,j,1)
                        J_av = 0.5d0*(sources%points(s)%J_mat(i,j,1) + &
                                        sources%points(s)%J_mat(i,j+1,1))
                        mxll%Ey(i_id,j_id) = mxll%Ey(i_id,j_id) + c_src * J_av
                    end if
                end do
                end do
            case ('z')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                    if (sources%points(s)%in_this_rank(i,j,1)) then
                        i_id = sources%points(s)%ind_i(i,j,1)
                        j_id = sources%points(s)%ind_j(i,j,1)
                        J_av = sources%points(s)%J_mat(i,j,1)
                        mxll%Ez(i_id,j_id) = mxll%Ez(i_id,j_id) + c_src * J_av
                    end if
                end do
                end do
            case default
                print *, "Error: Unknown source direction in 2D mxll."
                stop
            end select
        end do

    class is(TMxll_3D)
    
        do s=1, sources%n_p_src
            n_ker = sources%points(s)%n_ker
            select case (sources%points(s)%polarization)
            case ('x')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                do k=-n_ker, n_ker
                    if (sources%points(s)%in_this_rank(i,j,k)) then
                        i_id = sources%points(s)%ind_i(i,j,k)
                        j_id = sources%points(s)%ind_j(i,j,k)
                        k_id = sources%points(s)%ind_k(i,j,k)
                        J_av = 0.5d0*(sources%points(s)%J_mat(i,j,k) + &
                                        sources%points(s)%J_mat(i+1,j,k))
                        mxll%Ex(i_id,j_id,k_id) = mxll%Ex(i_id,j_id,k_id) + c_src * J_av
                    end if
                end do
                end do
                end do
            case ('y')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                do k=-n_ker, n_ker
                    if (sources%points(s)%in_this_rank(i,j,k)) then
                        i_id = sources%points(s)%ind_i(i,j,k)
                        j_id = sources%points(s)%ind_j(i,j,k)
                        k_id = sources%points(s)%ind_k(i,j,k)
                        J_av = 0.5d0*(sources%points(s)%J_mat(i,j,k) + &
                                        sources%points(s)%J_mat(i,j+1,k))
                        mxll%Ey(i_id,j_id,k_id) = mxll%Ey(i_id,j_id,k_id) + c_src * J_av
                    end if
                end do
                end do
                end do
            case ('z')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                do k=-n_ker, n_ker
                    if (sources%points(s)%in_this_rank(i,j,k)) then
                        i_id = sources%points(s)%ind_i(i,j,k)
                        j_id = sources%points(s)%ind_j(i,j,k)
                        k_id = sources%points(s)%ind_k(i,j,k)
                        J_av = 0.5d0*(sources%points(s)%J_mat(i,j,k) + &
                                        sources%points(s)%J_mat(i,j,k+1))
                        mxll%Ez(i_id,j_id,k_id) = mxll%Ez(i_id,j_id,k_id) + c_src * J_av
                    end if
                end do
                end do
                end do
            case default
                print *, "Error: Unknown source direction in 3D mxll."
                stop
            end select
        end do

    class default
        print *, "Error: Unknown mxll type in point_source_interactions."
        stop
    end select

end subroutine point_source_interactions

!###################################################################################################

subroutine plane_waves_E_interactions(mxll, sources, mpi_coords, mpi_dims, time)
    class(TMxll)       ,intent(inout) :: mxll
    type(TSources_list),intent(inout) :: sources
    integer            , intent(in)   :: mpi_coords(3)
    integer            , intent(in)   :: mpi_dims(3)
    real(dp)           , intent(in)   :: time

    integer  :: s
    integer  :: i_min, j_min, k_min
    integer  :: i_max, j_max, k_max
    integer  :: i, j ,k
    integer  :: i0, j0, k0
    integer  :: i_min_loc, j_min_loc, k_min_loc
    integer  :: i_max_loc, j_max_loc, k_max_loc
    integer  :: nx, ny, nz
    integer  :: d_int
    integer  :: m0 = 10
    real(dp) :: sin_psi
    real(dp) :: cos_psi
    real(dp) :: cos_phi
    real(dp) :: sin_phi
    real(dp) :: d, d_p, d_pp
    real(dp) :: dx_aux
    real(dp) :: dr_main
    real(dp) :: A_vec(3)
    real(dp) :: P_vec(3)
    real(dp) :: v_vec(3)
    real(dp) :: w_vec(3)
    real(dp) :: uz_vec(3) = (/0.0d0, 0.0d0, 1.0d0/)
    real(dp) :: v1_vec(3) = 0.0d0
    real(dp) :: v3_vec(3) = 0.0d0
    real(dp) :: E_vec(3)
    real(dp) :: E_inc
    real(dp) :: Ex_inc, Ey_inc, Ez_inc
    real(dp) :: dt_mu

    select type(mxll)
    class is(TMxll_1D)
    class is(TMxll_2D)

        nx       = mxll%nx
        ny       = mxll%ny
        dt_mu    = mxll%dt/mu0/mxll%dr
        dr_main  = mxll%dr

        do s=1, sources%n_pw_src

            i_min   = sources%plane_waves(s)%i_min
            i_max   = sources%plane_waves(s)%i_max
            j_min   = sources%plane_waves(s)%j_min
            j_max   = sources%plane_waves(s)%j_max

            P_vec    = 0.0d0
            w_vec    = 0.0d0
            A_vec    = sources%plane_waves(s)%A_vec
            v_vec    = sources%plane_waves(s)%v_vec
            dx_aux   = sources%plane_waves(s)%mxll_inc%dr
            
            cos_psi  = DCOS(sources%plane_waves(s)%psi)
            sin_psi  = DSIN(sources%plane_waves(s)%psi)
            cos_phi  = DCOS(sources%plane_waves(s)%phi - pi0/2)
            sin_phi  = DSIN(sources%plane_waves(s)%phi - pi0/2)

            if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then

                if (mxll%mode == TMZ_2D_MODE) sin_psi = 1.0d0

                if (sources%plane_waves(s)%i_min_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(1)) then

                    i_min_loc = sources%plane_waves(s)%i_min_loc

                    P_vec(1) = (i_min - INT(mpi_dims(1)*nx/2))*dr_main
                    
                    j0 = ny*mpi_coords(2)

                    do j = 1, ny 
                        if (j0+j <= j_max .and. j0+j >= j_min) then

                            P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_int      = FLOOR(d/dx_aux)
                            d_p        = (d - d_int*dx_aux)/dx_aux
                            
                            Ez_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                     d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                            mxll%Hy(i_min_loc-1,j) = mxll%Hy(i_min_loc-1,j) - &
                                                     dt_mu * sin_psi * Ez_inc
                    
                        end if
                    end do

                end if

                if (sources%plane_waves(s)%i_max_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(1)) then

                    i_max_loc = sources%plane_waves(s)%i_max_loc

                    P_vec(1) = (i_max - INT(mpi_dims(1)*nx/2))*dr_main
                    
                    j0 = ny*mpi_coords(2)

                    do j = 1, ny 
                        if (j0+j <= j_max .and. j0+j >= j_min) then

                            P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_int      = FLOOR(d/dx_aux)
                            d_p        = (d - d_int*dx_aux)/dx_aux
                            
                            Ez_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                     d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                            mxll%Hy(i_max_loc,j) = mxll%Hy(i_max_loc, j) + &
                                                   dt_mu * sin_psi * Ez_inc
                    
                        end if
                    end do

                end if

                if (sources%plane_waves(s)%j_min_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(2)) then

                    j_min_loc = sources%plane_waves(s)%j_min_loc

                    P_vec(2) = (j_min - INT(mpi_dims(2)*ny/2))*dr_main
                    
                    i0 = nx*mpi_coords(1)

                    do i = 1, nx 
                        if (i0+i <= i_max .and. i0+i >= i_min) then

                            P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_int      = FLOOR(d/dx_aux)
                            d_p        = (d - d_int*dx_aux)/dx_aux
                            
                            Ez_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                        d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                            mxll%Hx(i,j_min_loc-1) = mxll%Hx(i,j_min_loc-1) + &
                                                     dt_mu * sin_psi * Ez_inc
                        end if
                    end do

                end if

                if (sources%plane_waves(s)%j_max_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(2)) then
                    j_max_loc = sources%plane_waves(s)%j_max_loc

                    P_vec(2) = (j_max - INT(mpi_dims(2)*ny/2))*dr_main
                    
                    i0 = nx*mpi_coords(1)

                    do i = 1, nx 
                        if (i0+i <= i_max .and. i0+i >= i_min) then

                            P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_int      = FLOOR(d/dx_aux)
                            d_p        = (d - d_int*dx_aux)/dx_aux
                            
                            Ez_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                        d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                            mxll%Hx(i,j_max_loc) = mxll%Hx(i,j_max_loc) - &
                                                     dt_mu * sin_psi * Ez_inc
                        end if
                    end do

                end if
                
            end if

            if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then

                if (mxll%mode == TEZ_2D_MODE) cos_psi = 1.0d0

                if (sources%plane_waves(s)%i_min_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(1)) then

                    i_min_loc = sources%plane_waves(s)%i_min_loc

                    P_vec(1) = (i_min - INT(mpi_dims(1)*nx/2))*dr_main
                    
                    j0 = ny*mpi_coords(2)

                    do j = 1, ny 
                        if (j0+j <= j_max-1 .and. j0+j >= j_min) then

                            P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_int      = FLOOR(d/dx_aux)
                            d_p        = (d - d_int*dx_aux)/dx_aux

                            Ey_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                        d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                            mxll%Hz(i_min_loc-1,j) = mxll%Hz(i_min_loc-1,j) + &
                                                     dt_mu * sin_phi * cos_psi * Ey_inc

                        end if
                    end do
                end if

                if (sources%plane_waves(s)%i_max_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(1)) then

                    i_max_loc = sources%plane_waves(s)%i_max_loc

                    P_vec(1) = (i_max - INT(mpi_dims(1)*nx/2))*dr_main
                    
                    j0 = ny*mpi_coords(2)

                    do j = 1, ny 
                        if (j0+j <= j_max-1 .and. j0+j >= j_min) then

                            P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_int      = FLOOR(d/dx_aux)
                            d_p        = (d - d_int*dx_aux)/dx_aux

                            Ey_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                        d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                            mxll%Hz(i_max_loc,j) = mxll%Hz(i_max_loc,j) - &
                                                    dt_mu * sin_phi * cos_psi * Ey_inc
                        end if
                    end do
                end if

                if (sources%plane_waves(s)%j_min_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(2)) then

                    j_min_loc = sources%plane_waves(s)%j_min_loc

                    P_vec(2) = (j_min - INT(mpi_dims(2)*ny/2))*dr_main
                    
                    i0 = nx*mpi_coords(1)

                    do i = 1, nx 
                        if (i0+i <= i_max-1 .and. i0+i >= i_min) then

                            P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_int      = FLOOR(d/dx_aux)
                            d_p        = (d - d_int*dx_aux)/dx_aux

                            Ex_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                        d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)
                        
                            mxll%Hz(i,j_min_loc-1) = mxll%Hz(i,j_min_loc-1) - &
                                                     dt_mu * cos_phi * cos_psi * Ex_inc
                        end if
                    end do

                end if

                if (sources%plane_waves(s)%j_max_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(2)) then
                    j_max_loc = sources%plane_waves(s)%j_max_loc

                    P_vec(2) = (j_max - INT(mpi_dims(2)*ny/2))*dr_main
                    
                    i0 = nx*mpi_coords(1)

                    do i = 1, nx 
                        if (i0+i <= i_max-1 .and. i0+i >= i_min) then

                            P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_int      = FLOOR(d/dx_aux)
                            d_p        = (d - d_int*dx_aux)/dx_aux

                            Ex_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                        d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)
                        
                            mxll%Hz(i,j_max_loc) = mxll%Hz(i,j_max_loc) + &
                                                     dt_mu * cos_phi * cos_psi * Ex_inc
                        end if
                    end do
                end if

            end if

        end do
    class is(TMxll_3D)
    
        nx     = mxll%nx
        ny     = mxll%ny
        nz     = mxll%nz
        dt_mu  = mxll%dt/mu0/mxll%dr
        dr_main = mxll%dr
    
        do s=1, sources%n_pw_src

            i_min   = sources%plane_waves(s)%i_min
            i_max   = sources%plane_waves(s)%i_max
            j_min   = sources%plane_waves(s)%j_min
            j_max   = sources%plane_waves(s)%j_max
            k_min   = sources%plane_waves(s)%k_min
            k_max   = sources%plane_waves(s)%k_max

            P_vec    = 0.0d0
            w_vec    = 0.0d0
            A_vec    = sources%plane_waves(s)%A_vec
            v_vec    = sources%plane_waves(s)%v_vec
            dx_aux   = sources%plane_waves(s)%mxll_inc%dr

            v1_vec   = CROSS_PRODUCT(v_vec, uz_vec)
            v3_vec   = CROSS_PRODUCT(v1_vec, v_vec)

            cos_psi  = DCOS(sources%plane_waves(s)%psi)
            sin_psi  = DSIN(sources%plane_waves(s)%psi)

            if (sources%plane_waves(s)%i_min_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(1)) then

                i_min_loc = sources%plane_waves(s)%i_min_loc

                P_vec(1) = (i_min - INT(mpi_dims(1)*nx/2))*dr_main

                j0 = ny*mpi_coords(2)
                k0 = nz*mpi_coords(3)

                do k = 1, nz
                do j = 1, ny

                    if ((j0+j >= j_min .and. k0+k >= k_min) .and. &
                        (j0+j <= j_max-1 .and. k0+k <= k_max)) then

                        P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                        P_vec(3)   = (k0 + k - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux

                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ey_inc = E_vec(2)
                        mxll%Hz(i_min_loc-1,j,k) = mxll%Hz(i_min_loc-1,j,k) + dt_mu * Ey_inc

                    end if
                        
                    if ((j0+j >= j_min .and. k0+k >= k_min) .and. &
                        (j0+j <= j_max .and. k0+k <= k_max-1)) then
                        
                        P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                        P_vec(3)   = (k0 + k + 0.5d0- INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux
                        
                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                        d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)
                        
                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ez_inc = E_vec(3)
                        mxll%Hy(i_min_loc-1,j,k) = mxll%Hy(i_min_loc-1,j,k) - dt_mu * Ez_inc

                    end if 
                end do
                end do

            end if

            if (sources%plane_waves(s)%i_max_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(1)) then

                i_max_loc = sources%plane_waves(s)%i_max_loc
                P_vec(1) = (i_max - INT(mpi_dims(1)*nx/2))*dr_main

                j0 = ny*mpi_coords(2)
                k0 = nz*mpi_coords(3)

                do k = 1, nz
                do j = 1, ny
                    if ((j0+j >= j_min .and. k0+k >= k_min) .and. &
                        (j0+j <= j_max-1 .and. k0+k <= k_max)) then

                        P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                        P_vec(3)   = (k0 + k - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux

                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ey_inc = E_vec(2)
                        mxll%Hz(i_max_loc,j,k) = mxll%Hz(i_max_loc,j,k) - dt_mu * Ey_inc

                    end if
                        
                    if ((j0+j >= j_min .and. k0+k >= k_min) .and. &
                        (j0+j <= j_max .and. k0+k <= k_max-1)) then
                        
                        P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                        P_vec(3)   = (k0 + k + 0.5d0- INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux
                        
                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ez_inc = E_vec(3)
                        mxll%Hy(i_max_loc,j,k) = mxll%Hy(i_max_loc,j,k) + dt_mu * Ez_inc

                    end if

                end do
                end do

            end if

            if (sources%plane_waves(s)%j_min_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(2)) then

                j_min_loc = sources%plane_waves(s)%j_min_loc
                P_vec(2) = (j_min - INT(mpi_dims(2)*ny/2))*dr_main

                i0 = nx*mpi_coords(1)
                k0 = nz*mpi_coords(3)

                do k = 1, nz
                do i = 1, nx
                    if ((i0+i >= i_min .and. k0+k >= k_min) .and. &
                        (i0+i <= i_max-1 .and. k0+k <= k_max)) then

                        P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(3)   = (k0 + k - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux

                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ex_inc = E_vec(1)
                        mxll%Hz(i,j_min_loc-1,k) = mxll%Hz(i,j_min_loc-1,k) - dt_mu * Ex_inc

                    end if
                        
                    if ((i0+i >= i_min .and. k0+k >= k_min) .and. &
                        (i0+i <= i_max .and. k0+k <= k_max-1)) then
                        
                        P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(3)   = (k0 + k + 0.5d0 - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux
                        
                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ez_inc = E_vec(3)
                        mxll%Hx(i,j_min_loc-1,k) = mxll%Hx(i,j_min_loc-1,k) + dt_mu * Ez_inc

                    end if

                end do
                end do

            end if

            if (sources%plane_waves(s)%j_max_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(2)) then

                j_max_loc = sources%plane_waves(s)%j_max_loc
                P_vec(2) = (j_max - INT(mpi_dims(2)*ny/2))*dr_main

                i0 = nx*mpi_coords(1)
                k0 = nz*mpi_coords(3)

                do k = 1, nz
                do i = 1, nx
                    if ((i0+i >= i_min .and. k0+k >= k_min) .and. &
                        (i0+i <= i_max-1 .and. k0+k <= k_max)) then

                        P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(3)   = (k0 + k - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux

                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ex_inc = E_vec(1)
                        mxll%Hz(i,j_max_loc,k) = mxll%Hz(i,j_max_loc,k) + dt_mu * Ex_inc

                    end if
                        
                    if ((i0+i >= i_min .and. k0+k >= k_min) .and. &
                        (i0+i <= i_max .and. k0+k <= k_max-1)) then
                        
                        P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(3)   = (k0 + k + 0.5d0 - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux

                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ez_inc = E_vec(3)
                        mxll%Hx(i,j_max_loc,k) = mxll%Hx(i,j_max_loc,k) - dt_mu * Ez_inc

                    end if

                end do
                end do
            end if

            if (sources%plane_waves(s)%k_min_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(3)) then

                k_min_loc = sources%plane_waves(s)%k_min_loc
                P_vec(3) = (k_min - INT(mpi_dims(3)*nz/2))*dr_main

                i0 = nx*mpi_coords(1)
                j0 = ny*mpi_coords(2)

                do j = 1, ny
                do i = 1, nx
                    if ((i0+i >= i_min .and. j0+j >= j_min) .and. &
                        (i0+i <= i_max-1 .and. j0+j <= j_max)) then

                        P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux

                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ex_inc = E_vec(1)

                        mxll%Hy(i,j,k_min_loc-1) = mxll%Hy(i,j,k_min_loc-1) + dt_mu * Ex_inc

                    end if

                    if ((i0+i >= i_min .and. j0+j >= j_min) .and. &
                        (i0+i <= i_max .and. j0+j <= j_max-1)) then

                        P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux

                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ey_inc = E_vec(2)

                        mxll%Hx(i,j,k_min_loc-1) = mxll%Hx(i,j,k_min_loc-1) - dt_mu * Ey_inc

                    end if

                end do
                end do

            end if

            if (sources%plane_waves(s)%k_max_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(3)) then

                k_max_loc = sources%plane_waves(s)%k_max_loc
                P_vec(3) = (k_max - INT(mpi_dims(3)*nz/2))*dr_main

                i0 = nx*mpi_coords(1)
                j0 = ny*mpi_coords(2)

                do j = 1, ny
                do i = 1, nx
                    if ((i0+i >= i_min .and. j0+j >= j_min) .and. &
                        (i0+i <= i_max-1 .and. j0+j <= j_max)) then

                        P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux

                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ex_inc = E_vec(1)

                        mxll%Hy(i,j,k_max_loc) = mxll%Hy(i,j,k_max_loc) - dt_mu * Ex_inc

                    end if

                    if ((i0+i >= i_min .and. j0+j >= j_min) .and. &
                        (i0+i <= i_max .and. j0+j <= j_max-1)) then

                        P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_int      = FLOOR(d/dx_aux)
                        d_p        = (d - d_int*dx_aux)/dx_aux
                        
                        E_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0) + &
                                  d_p   * sources%plane_waves(s)%mxll_inc%Ex(d_int+m0+1)

                        E_vec = E_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Ey_inc = E_vec(2)

                        mxll%Hx(i,j,k_max_loc) = mxll%Hx(i,j,k_max_loc) + dt_mu * Ey_inc

                    end if

                end do
                end do

            end if

        end do

    end select

end subroutine plane_waves_E_interactions
!###################################################################################################

subroutine plane_waves_H_interactions(mxll, sources, mpi_coords, mpi_dims, time)
    class(TMxll)       ,intent(inout) :: mxll
    type(TSources_list),intent(inout) :: sources
    integer            , intent(in)   :: mpi_coords(3)
    integer            , intent(in)   :: mpi_dims(3)
    real(dp)           , intent(in)   :: time

    integer  :: s
    integer  :: i_min, j_min, k_min
    integer  :: i_max, j_max, k_max
    integer  :: i, j ,k
    integer  :: i0, j0, k0
    integer  :: i_min_loc, j_min_loc, k_min_loc
    integer  :: i_max_loc, j_max_loc, k_max_loc
    integer  :: nx, ny, nz
    integer  :: d_int
    integer  :: m0 = 10
    real(dp) :: sin_psi
    real(dp) :: cos_psi
    real(dp) :: cos_phi
    real(dp) :: sin_phi
    real(dp) :: d, d_p, d_pp
    real(dp) :: dx_aux
    real(dp) :: dr_main
    real(dp) :: A_vec(3)
    real(dp) :: P_vec(3)
    real(dp) :: v_vec(3)
    real(dp) :: w_vec(3)
    real(dp) :: uz_vec(3) = (/0.0d0, 0.0d0, 1.0d0/)
    real(dp) :: v1_vec(3) = 0.0d0
    real(dp) :: v3_vec(3) = 0.0d0
    real(dp) :: H_vec(3)
    real(dp) :: H_inc
    real(dp) :: Hx_inc, Hy_inc, Hz_inc
    real(dp) :: dt_eps

    select type(mxll)
    class is(TMxll_1D)
    class is(TMxll_2D)

        nx       = mxll%nx
        ny       = mxll%ny
        dt_eps   = mxll%dt/eps0/mxll%dr
        dr_main  = mxll%dr

        do s = 1, sources%n_pw_src

            i_min   = sources%plane_waves(s)%i_min
            i_max   = sources%plane_waves(s)%i_max
            j_min   = sources%plane_waves(s)%j_min
            j_max   = sources%plane_waves(s)%j_max

            P_vec    = 0.0d0
            w_vec    = 0.0d0
            A_vec    = sources%plane_waves(s)%A_vec
            v_vec    = sources%plane_waves(s)%v_vec
            dx_aux   = sources%plane_waves(s)%mxll_inc%dr
            
            cos_psi  = DCOS(sources%plane_waves(s)%psi - pi0/2)
            sin_psi  = DSIN(sources%plane_waves(s)%psi - pi0/2)
            cos_phi  = DCOS(sources%plane_waves(s)%phi - pi0/2)
            sin_phi  = DSIN(sources%plane_waves(s)%phi - pi0/2)

            if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then

                if (mxll%mode == TMZ_2D_MODE) cos_psi = 1.0d0

                if (sources%plane_waves(s)%i_min_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(1)) then

                    i_min_loc = sources%plane_waves(s)%i_min_loc

                    P_vec(1) = (i_min - 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                    
                    j0 = ny*mpi_coords(2)

                    do j = 1, ny 
                        if (j0+j <= j_max .and. j0+j >= j_min) then

                            P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_pp       = d + 0.5d0*dx_aux
                            d_int      = FLOOR(d_pp/dx_aux)
                            d_p        = (d_pp - d_int*dx_aux)/dx_aux
                            
                            Hy_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                      d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                            mxll%Ez(i_min_loc, j) = mxll%Ez(i_min_loc, j) - &
                                                     dt_eps * sin_phi * cos_psi * Hy_inc

                        end if
                    end do

                end if

                if (sources%plane_waves(s)%i_max_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(1)) then

                    i_max_loc = sources%plane_waves(s)%i_max_loc

                    P_vec(1) = (i_max + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                    
                    j0 = ny*mpi_coords(2)

                    do j = 1, ny 
                        if (j0+j <= j_max .and. j0+j >= j_min) then

                            P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_pp       = d + 0.5d0*dx_aux
                            d_int      = FLOOR(d_pp/dx_aux)
                            d_p        = (d_pp - d_int*dx_aux)/dx_aux
                            
                            Hy_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                      d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                            mxll%Ez(i_max_loc, j) = mxll%Ez(i_max_loc, j) + &
                                                       dt_eps * sin_phi * cos_psi * Hy_inc
                                      
                        end if
                    end do

                end if

                if (sources%plane_waves(s)%j_min_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(2)) then

                    j_min_loc = sources%plane_waves(s)%j_min_loc

                    P_vec(2) = (j_min - 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                    
                    i0 = nx*mpi_coords(1)

                    do i = 1, nx 
                        if (i0+i <= i_max .and. i0+i >= i_min) then

                            P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_pp       = d + 0.5d0*dx_aux
                            d_int      = FLOOR(d_pp/dx_aux)
                            d_p        = (d_pp - d_int*dx_aux)/dx_aux

                            Hx_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                      d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)
                        
                            mxll%Ez(i,j_min_loc) = mxll%Ez(i,j_min_loc) + &
                                                     dt_eps * cos_phi * cos_psi * Hx_inc
                        end if
                    end do

                end if

                if (sources%plane_waves(s)%j_max_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(2)) then
                    j_max_loc = sources%plane_waves(s)%j_max_loc

                    P_vec(2) = (j_max + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                    
                    i0 = nx*mpi_coords(1)

                    do i = 1, nx 
                        if (i0+i <= i_max .and. i0+i >= i_min) then

                            P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_pp       = d + 0.5d0*dx_aux
                            d_int      = FLOOR(d_pp/dx_aux)
                            d_p        = (d_pp - d_int*dx_aux)/dx_aux

                            Hx_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                      d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)
                        
                            mxll%Ez(i,j_max_loc) = mxll%Ez(i,j_max_loc) - &
                                                     dt_eps * cos_phi * cos_psi * Hx_inc
                        end if
                    end do

                end if

            end if

            if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then

                if (mxll%mode == TEZ_2D_MODE) sin_psi = -1.0d0

                if (sources%plane_waves(s)%i_min_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(1)) then

                    i_min_loc = sources%plane_waves(s)%i_min_loc

                    P_vec(1) = (i_min - 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                    
                    j0 = ny*mpi_coords(2)

                    do j = 1, ny 
                        if (j0+j <= j_max-1 .and. j0+j >= j_min) then

                            P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_pp       = d + 0.5d0*dx_aux
                            d_int      = FLOOR(d_pp/dx_aux)
                            d_p        = (d_pp - d_int*dx_aux)/dx_aux

                            Hz_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                            mxll%Ey(i_min_loc, j) = mxll%Ey(i_min_loc, j) + &
                                                    dt_eps * sin_psi * Hz_inc
                            
                        end if
                    end do
                end if
    
                if (sources%plane_waves(s)%i_max_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(1)) then

                    i_max_loc = sources%plane_waves(s)%i_max_loc

                    P_vec(1) = (i_max + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main

                    j0 = ny*mpi_coords(2)

                    do j = 1, ny 
                        if (j0+j <= j_max-1 .and. j0+j >= j_min) then

                            P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_pp       = d + 0.5d0*dx_aux
                            d_int      = FLOOR(d_pp/dx_aux)
                            d_p        = (d_pp - d_int*dx_aux)/dx_aux

                            Hz_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                            mxll%Ey(i_max_loc, j) = mxll%Ey(i_max_loc, j) - &
                                                    dt_eps * sin_psi * Hz_inc
                            
                        end if
                    end do
                end if
    
                if (sources%plane_waves(s)%j_min_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(2)) then

                    j_min_loc = sources%plane_waves(s)%j_min_loc

                    P_vec(2) = (j_min - 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                    
                    i0 = nx*mpi_coords(1)

                    do i = 1, nx 
                        if (i0+i <= i_max-1 .and. i0+i >= i_min) then

                            P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_pp       = d + 0.5d0*dx_aux
                            d_int      = FLOOR(d_pp/dx_aux)
                            d_p        = (d_pp - d_int*dx_aux)/dx_aux

                            Hz_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                            mxll%Ex(i,j_min_loc) = mxll%Ex(i,j_min_loc) - &
                                                    dt_eps * sin_psi * Hz_inc
                        end if
                    end do
                end if

                if (sources%plane_waves(s)%j_max_in_this_rank .and. &
                    sources%plane_waves(s)%limited_axis(2)) then
                    j_max_loc = sources%plane_waves(s)%j_max_loc

                    P_vec(2) = (j_max + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                    
                    i0 = nx*mpi_coords(1)

                    do i = 1, nx 
                        if (i0+i <= i_max-1 .and. i0+i >= i_min) then

                            P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                            w_vec(1:2) = P_vec(1:2) - A_vec(1:2)  
                            d          = DOT_PRODUCT(w_vec(1:2),v_vec(1:2))
                            d_pp       = d + 0.5d0*dx_aux
                            d_int      = FLOOR(d_pp/dx_aux)
                            d_p        = (d_pp - d_int*dx_aux)/dx_aux

                            Hz_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                            mxll%Ex(i,j_max_loc) = mxll%Ex(i,j_max_loc) + &
                                                    dt_eps * sin_psi * Hz_inc
                        end if
                    end do
                end if

            end if

        end do

    class is(TMxll_3D)

        nx      = mxll%nx
        ny      = mxll%ny
        nz      = mxll%nz
        dt_eps  = mxll%dt/eps0/mxll%dr
        dr_main = mxll%dr
    
        do s=1, sources%n_pw_src

            i_min   = sources%plane_waves(s)%i_min
            i_max   = sources%plane_waves(s)%i_max
            j_min   = sources%plane_waves(s)%j_min
            j_max   = sources%plane_waves(s)%j_max
            k_min   = sources%plane_waves(s)%k_min
            k_max   = sources%plane_waves(s)%k_max

            P_vec    = 0.0d0
            w_vec    = 0.0d0
            A_vec    = sources%plane_waves(s)%A_vec
            v_vec    = sources%plane_waves(s)%v_vec
            dx_aux   = sources%plane_waves(s)%mxll_inc%dr

            v1_vec   = CROSS_PRODUCT(v_vec, uz_vec)
            v3_vec   = CROSS_PRODUCT(v1_vec, v_vec)

            cos_psi  = DCOS(sources%plane_waves(s)%psi-pi0/2)
            sin_psi  = DSIN(sources%plane_waves(s)%psi-pi0/2)

            if (sources%plane_waves(s)%i_min_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(1)) then

                i_min_loc = sources%plane_waves(s)%i_min_loc

                P_vec(1) = (i_min - 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                
                j0 = ny*mpi_coords(2)
                k0 = nz*mpi_coords(3)

                i_min_loc = sources%plane_waves(s)%i_min_loc

                P_vec(1) = (i_min - 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                
                j0 = ny*mpi_coords(2)
                k0 = nz*mpi_coords(3)

                do k = 1, nz
                do j = 1, ny 
                    if ((j0+j >= j_min .and. k0+k >= k_min) .and. &
                        (j0+j <= j_max-1 .and. k0+k <= k_max)) then

                        P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                        P_vec(3)   = (k0 + k - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux
                        
                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hz_inc = H_vec(3)

                        mxll%Ey(i_min_loc, j, k) = mxll%Ez(i_min_loc, j, k) + dt_eps * Hz_inc

                    end if

                    if ((j0+j >= j_min .and. k0+k >= k_min) .and. &
                        (j0+j <= j_max .and. k0+k <= k_max-1)) then
                        
                        P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                        P_vec(3)   = (k0 + k + 0.5d0 - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux

                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hy_inc = H_vec(2)

                        mxll%Ex(i_min_loc, j, k) = mxll%Ex(i_min_loc, j, k) - dt_eps * Hy_inc
                    
                    end if

                end do
                end do

            end if

            if (sources%plane_waves(s)%i_max_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(1)) then

                i_max_loc = sources%plane_waves(s)%i_max_loc

                P_vec(1) = (i_max + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                
                j0 = ny*mpi_coords(2)
                k0 = nz*mpi_coords(3)

                do k = 1, nz
                do j = 1, ny 
                    if ((j0+j >= j_min .and. k0+k >= k_min) .and. &
                        (j0+j <= j_max-1 .and. k0+k <= k_max)) then

                        P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                        P_vec(3)   = (k0 + k - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux
                        
                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hz_inc = H_vec(3)

                        mxll%Ey(i_max_loc, j, k) = mxll%Ey(i_max_loc, j, k) - dt_eps * Hz_inc

                    end if

                    if ((j0+j >= j_min .and. k0+k >= k_min) .and. &
                        (j0+j <= j_max .and. k0+k <= k_max-1)) then
                        
                        P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                        P_vec(3)   = (k0 + k + 0.5d0 - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux

                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hy_inc = H_vec(2)

                        mxll%Ez(i_max_loc, j, k) = mxll%Ez(i_max_loc, j, k) + dt_eps * Hy_inc
                    end if

                end do
                end do
            
            end if

            if (sources%plane_waves(s)%j_min_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(2)) then

                j_min_loc = sources%plane_waves(s)%j_min_loc

                P_vec(2) = (j_min - 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                
                i0 = nx*mpi_coords(1)
                k0 = nz*mpi_coords(3)

                do k = 1, nz
                do i = 1, nx 
                    if ((i0+i >= i_min .and. k0+k >= k_min) .and. &
                        (i0+i <= i_max-1 .and. k0+k <= k_max)) then

                        P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(3)   = (k0 + k - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux

                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hz_inc = H_vec(3)

                        mxll%Ex(i, j_min_loc, k) = mxll%Ex(i, j_min_loc, k) - dt_eps * Hz_inc
                        
                    end if  

                    if ((i0+i >= i_min .and. k0+k >= k_min) .and. &
                        (i0+i <= i_max .and. k0+k <= k_max-1)) then
                        
                        P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(3)   = (k0 + k + 0.5d0 - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux

                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hx_inc = H_vec(1)

                        mxll%Ez(i, j_min_loc, k) = mxll%Ez(i, j_min_loc, k) + dt_eps * Hx_inc
                    
                    end if

                end do
                end do

            end if

            if (sources%plane_waves(s)%j_max_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(2)) then

                j_max_loc = sources%plane_waves(s)%j_max_loc

                P_vec(2) = (j_max + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                
                i0 = nx*mpi_coords(1)
                k0 = nz*mpi_coords(3)

                do k = 1, nz
                do i = 1, nx 
                    if ((i0+i >= i_min .and. k0+k >= k_min) .and. &
                        (i0+i <= i_max-1 .and. k0+k <= k_max)) then

                        P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(3)   = (k0 + k - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux

                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hz_inc = H_vec(3)

                        mxll%Ex(i, j_max_loc, k) = mxll%Ex(i, j_max_loc, k) + dt_eps * Hz_inc
                    
                    end if  

                    if ((i0+i >= i_min .and. k0+k >= k_min) .and. &
                        (i0+i <= i_max .and. k0+k <= k_max-1)) then
                        
                        P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(3)   = (k0 + k + 0.5d0 - INT(mpi_dims(3)*nz/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux

                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hx_inc = H_vec(1)

                        mxll%Ez(i, j_max_loc, k) = mxll%Ez(i, j_max_loc, k) - dt_eps * Hx_inc

                    end if

                end do
                end do

            end if

            if (sources%plane_waves(s)%k_min_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(3)) then

                k_min_loc = sources%plane_waves(s)%k_min_loc

                P_vec(3) = (k_min - 0.5d0 - INT(mpi_dims(3)*nz/2))*dr_main
                
                i0 = nx*mpi_coords(1)
                j0 = ny*mpi_coords(2)

                do j = 1, ny
                do i = 1, nx 
                    if ((i0+i >= i_min .and. j0+j >= j_min) .and. &
                        (i0+i <= i_max-1 .and. j0+j <= j_max)) then

                        P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux

                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hy_inc = H_vec(2)

                        mxll%Ex(i, j, k_min_loc) = mxll%Ex(i, j, k_min_loc) + dt_eps * Hy_inc
                    
                    end if  

                    if ((i0+i >= i_min .and. j0+j >= j_min) .and. &
                        (i0+i <= i_max .and. j0+j <= j_max-1)) then
                        
                        P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux

                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hx_inc = H_vec(1)

                        mxll%Ey(i, j, k_min_loc) = mxll%Ey(i, j, k_min_loc) - dt_eps * Hx_inc

                    end if

                end do
                end do

            end if

            if (sources%plane_waves(s)%k_max_in_this_rank .and. &
                sources%plane_waves(s)%limited_axis(3)) then

                k_max_loc = sources%plane_waves(s)%k_max_loc

                P_vec(3) = (k_max + 0.5d0 - INT(mpi_dims(3)*nz/2))*dr_main
                
                i0 = nx*mpi_coords(1)
                j0 = ny*mpi_coords(2)

                do j = 1, ny
                do i = 1, nx 
                    if ((i0+i >= i_min .and. j0+j >= j_min) .and. &
                        (i0+i <= i_max-1 .and. j0+j <= j_max)) then

                        P_vec(1)   = (i0 + i + 0.5d0 - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(2)   = (j0 + j - INT(mpi_dims(2)*ny/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux

                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hy_inc = H_vec(2)

                        mxll%Ex(i, j, k_max_loc) = mxll%Ex(i, j, k_max_loc) - dt_eps * Hy_inc
                    
                    end if  

                    if ((i0+i >= i_min .and. j0+j >= j_min) .and. &
                        (i0+i <= i_max .and. j0+j <= j_max-1)) then
                        
                        P_vec(1)   = (i0 + i - INT(mpi_dims(1)*nx/2))*dr_main
                        P_vec(2)   = (j0 + j + 0.5d0 - INT(mpi_dims(2)*ny/2))*dr_main
                        w_vec(1:3) = P_vec(1:3) - A_vec(1:3)  
                        d          = DOT_PRODUCT(w_vec(1:3),v_vec(1:3))
                        d_pp       = d + 0.5d0*dx_aux
                        d_int      = FLOOR(d_pp/dx_aux)
                        d_p        = (d_pp - d_int*dx_aux)/dx_aux

                        H_inc = (1-d_p) * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0-1) + &
                                    d_p   * sources%plane_waves(s)%mxll_inc%Hy(d_int+m0)

                        H_vec = H_inc * (cos_psi * v1_vec + sin_psi * v3_vec)
                        Hx_inc = H_vec(1)

                        mxll%Ey(i, j, k_max_loc) = mxll%Ey(i, j, k_max_loc) + dt_eps * Hx_inc

                    end if

                end do
                end do
            
            end if

        end do

    end select

end subroutine plane_waves_H_interactions

!###################################################################################################    
subroutine send_E_to_J_ranks(mxll, q_group, move_q_system, myrank)

    class(TMxll)  , intent(inout) :: mxll
    type(TQ_Group), intent(inout) :: q_group
    logical       , intent(in)    :: move_q_system
    integer       , intent(in)    :: myrank
    integer :: ierr

    if (.not. move_q_system) return

#ifdef USE_MPI
    ! call MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif

    select type(mxll)
    class is(TMxll_1D)

        call send_E_1D_to_J_ranks(mxll, q_group, myrank)

    class is(TMxll_2D)

        call send_E_2D_to_J_ranks(mxll, q_group, myrank)
    class is(TMxll_3D)

        call send_E_3D_to_J_ranks(mxll, q_group, myrank)
    end select

#ifdef USE_MPI
    ! call MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif

end subroutine send_E_to_J_ranks

!###################################################################################################

subroutine send_J_to_E_ranks(mxll, q_group, move_q_system, myrank)

    class(TMxll),  intent(inout)  :: mxll
    type(TQ_Group), intent(inout) :: q_group
    logical       , intent(in)    :: move_q_system
    integer       , intent(in)    :: myrank
    integer :: ierr

    if (.not. move_q_system) return

#ifdef USE_MPI
    ! call MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif

    select type(mxll)
    class is(TMxll_1D)

        call send_J_to_E_1D_ranks(mxll, q_group, myrank)

    class is(TMxll_2D)

        call send_J_to_E_2D_ranks(mxll, q_group, myrank)
    class is(TMxll_3D)

        call send_J_to_E_3D_ranks(mxll, q_group, myrank)
    end select

#ifdef USE_MPI
    ! call MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif

end subroutine send_J_to_E_ranks

!###################################################################################################

subroutine send_E_1D_to_J_ranks(mxll, q_group, myrank)
    class(TMxll_1D),  intent(inout)  :: mxll
    type(TQ_Group),   intent(inout)  :: q_group
    integer       , intent(in)       :: myrank

    integer :: n
    integer :: i_idx
    integer :: n_mol
    real(dp) :: E_field_send
    real(dp) :: E_field_get

#ifdef USE_MPI

    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
    
    select case(q_group%group_type)

    case(Q_MATERIAL)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)

            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then
                
                q_group%E_field_list(n_mol, 1) = mxll%Ex(i_idx)
                n_mol = n_mol + 1

            else

                if (q_group%map(n,2) == myrank) then

                    E_field_send = mxll%Ex(i_idx)
                    call mpi_send(E_field_send,1,mpi_double_precision, &
                                  q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)


                else if (q_group%map(n,1) == myrank) then

                    call mpi_recv(E_field_get,1,mpi_double_precision, &
                                  q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get
                    n_mol = n_mol + 1

                end if
            end if
        end do

    case(Q_SINGLE)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%kernel_map(n,0,0,0,3)

            if(q_group%kernel_map(n,0,0,0,1) == myrank .and. q_group%kernel_map(n,0,0,0,2) == myrank) then
            
                q_group%E_field_list(n_mol, 1) = mxll%Ex(i_idx)
                n_mol = n_mol + 1

            else

                if (q_group%kernel_map(n,0,0,0,2) == myrank) then

                    E_field_send = mxll%Ex(i_idx)

                    call mpi_send(E_field_send,1,mpi_double_precision, &
                                  q_group%kernel_map(n,0,0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)


                else if (q_group%kernel_map(n,0,0,0,1) == myrank) then

                    call mpi_recv(E_field_get,1,mpi_double_precision, &
                                  q_group%kernel_map(n,0,0,0,2), MPI_GOOD_TAG, MPI_COMM_WORLD, &
                                  istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get
                    n_mol = n_mol + 1

                end if
            end if
        end do

    end select

#else

    select case(q_group%group_type)

    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            n_mol = n
            i_idx = q_group%map(n,3)
            q_group%E_field_list(n_mol, 1) = mxll%Ex(i_idx)
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
            n_mol = n
            i_idx = q_group%kernel_map(n,0,0,0,3)
            q_group%E_field_list(n_mol, 1) = mxll%Ex(i_idx)
        end do
    end select

#endif

end subroutine send_E_1D_to_J_ranks

!###################################################################################################

subroutine send_E_2D_to_J_ranks(mxll, q_group, myrank)
    class(TMxll_2D),  intent(inout) :: mxll
    type(TQ_Group) ,  intent(inout) :: q_group
    integer        , intent(in)     :: myrank
 
    integer :: n
    integer :: i_idx
    integer :: j_idx
    integer :: n_mol
    real(dp) :: E_field_send(3)
    real(dp) :: E_field_get(3)

#ifdef USE_MPI

    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    E_field_send = M_ZERO
    E_field_get  = M_ZERO

    select case(q_group%group_type)
    case(Q_MATERIAL)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            
            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then
            
                if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                
                    q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                            mxll%Ex(i_idx, j_idx))
                    q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                            mxll%Ey(i_idx, j_idx))
                else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                    q_group%E_field_list(n_mol, 3) = mxll%Ez(i_idx, j_idx)

                end if
                
                n_mol = n_mol + 1

            else
                if (q_group%map(n,2) == myrank) then

                    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        E_field_send(1) =  0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                    mxll%Ex(i_idx, j_idx))
                        E_field_send(2) =  0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                    mxll%Ey(i_idx, j_idx))
                    else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        E_field_send(3) = mxll%Ez(i_idx, j_idx)
                    end if

                    call mpi_send(E_field_send,3,mpi_double_precision, &
                                    q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                else if (q_group%map(n,1) == myrank) then

                    call mpi_recv(E_field_get,3,mpi_double_precision, &
                                    q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get(1)
                    q_group%E_field_list(n_mol, 2) = E_field_get(2)
                    q_group%E_field_list(n_mol, 3) = E_field_get(3) 
                    n_mol = n_mol + 1

                end if
            end if
        end do

    case(Q_SINGLE)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%kernel_map(n,0,0,0,3)
            j_idx = q_group%kernel_map(n,0,0,0,4)
            
                if(q_group%kernel_map(n,0,0,0,1) == myrank .and. &
                    q_group%kernel_map(n,0,0,0,2) == myrank) then
            
                if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                
                    q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                            mxll%Ex(i_idx, j_idx))
                    q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                            mxll%Ey(i_idx, j_idx))

                else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                    q_group%E_field_list(n_mol, 3) = mxll%Ez(i_idx, j_idx)
                end if
                
                n_mol = n_mol + 1

            else
                if (q_group%kernel_map(n,0,0,0,2) == myrank) then

                    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        E_field_send(1) =  0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                    mxll%Ex(i_idx, j_idx))
                        E_field_send(2) =  0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                    mxll%Ey(i_idx, j_idx))
                    else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        E_field_send(3) = mxll%Ez(i_idx, j_idx)
                    end if
                
                    call mpi_send(E_field_send,3,mpi_double_precision, &
                                    q_group%kernel_map(n,0,0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)
                
                else if (q_group%kernel_map(n,0,0,0,1) == myrank) then
                    
                    call mpi_recv(E_field_get,3,mpi_double_precision, &
                                    q_group%kernel_map(n,0,0,0,2), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get(1)
                    q_group%E_field_list(n_mol, 2) = E_field_get(2)
                    q_group%E_field_list(n_mol, 3) = E_field_get(3) 
                    n_mol = n_mol + 1

                end if
            end if
        end do

    end select

#else

    select case(q_group%group_type)
    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            
            n_mol = n
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            
            if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                
                q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                        mxll%Ex(i_idx, j_idx))
                q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                        mxll%Ey(i_idx, j_idx))
            else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                q_group%E_field_list(n_mol, 3) = mxll%Ez(i_idx, j_idx)

            end if
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
            
            n_mol = n
            i_idx = q_group%kernel_map(n,0,0,0,3)
            j_idx = q_group%kernel_map(n,0,0,0,4)
            
            if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                
                q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                        mxll%Ex(i_idx, j_idx))
                q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                        mxll%Ey(i_idx, j_idx))
            else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                q_group%E_field_list(n_mol, 3) = mxll%Ez(i_idx, j_idx)
            end if
        end do
    
    end select

#endif

end subroutine send_E_2D_to_J_ranks

!###################################################################################################

subroutine send_E_3D_to_J_ranks(mxll, q_group, myrank)
    class(TMxll_3D),  intent(inout) :: mxll
    type(TQ_Group) ,  intent(inout) :: q_group
    integer        , intent(in)     :: myrank
 
    integer :: n
    integer :: i_idx
    integer :: j_idx
    integer :: k_idx
    integer :: n_mol
    real(dp) :: E_field_send(3)
    real(dp) :: E_field_get(3)
    
#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
    

    E_field_send = M_ZERO
    E_field_get  = M_ZERO

    select case(q_group%group_type)
    case(Q_MATERIAL)

        n_mol = 1
        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            k_idx = q_group%map(n,5)
            
            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then
            
                q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                                        mxll%Ex(i_idx, j_idx, k_idx))
                q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                                        mxll%Ey(i_idx, j_idx, k_idx))
                q_group%E_field_list(n_mol, 3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                                        mxll%Ez(i_idx, j_idx, k_idx))
                n_mol = n_mol + 1

            else
                
                if (q_group%map(n,2) == myrank) then

                    E_field_send(1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                             mxll%Ex(i_idx, j_idx, k_idx))
                    E_field_send(2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                             mxll%Ey(i_idx, j_idx, k_idx))
                    E_field_send(3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                             mxll%Ez(i_idx, j_idx, k_idx))

                    call mpi_send(E_field_send,3,mpi_double_precision, &
                                  q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                else if (q_group%map(n,1) == myrank) then

                    call mpi_recv(E_field_get,3,mpi_double_precision, &
                                  q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get(1)
                    q_group%E_field_list(n_mol, 2) = E_field_get(2)
                    q_group%E_field_list(n_mol, 3) = E_field_get(3) 
                    n_mol = n_mol + 1

                end if

            end if
        end do

    case(Q_SINGLE)

        n_mol = 1
        do n = 1, q_group%n_systems
            
            i_idx = q_group%kernel_map(n,0,1,1,3)
            j_idx = q_group%kernel_map(n,0,1,1,4)
            k_idx = q_group%kernel_map(n,0,1,1,5)
            
                if(q_group%kernel_map(n,0,0,0,1) == myrank .and. &
                    q_group%kernel_map(n,0,0,0,2) == myrank) then
            
                q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                                        mxll%Ex(i_idx, j_idx, k_idx))
                q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                                        mxll%Ey(i_idx, j_idx, k_idx))
                q_group%E_field_list(n_mol, 3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                                        mxll%Ez(i_idx, j_idx, k_idx))
                n_mol = n_mol + 1

            else
                
                if (q_group%kernel_map(n,0,0,0,2) == myrank) then

                    E_field_send(1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                             mxll%Ex(i_idx, j_idx, k_idx))
                    E_field_send(2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                             mxll%Ey(i_idx, j_idx, k_idx))
                    E_field_send(3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                             mxll%Ez(i_idx, j_idx, k_idx))

                    call mpi_send(E_field_send,3,mpi_double_precision, &
                                  q_group%kernel_map(n,0,0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                else if (q_group%kernel_map(n,0,0,0,1) == myrank) then

                    call mpi_recv(E_field_get,3,mpi_double_precision, &
                                  q_group%kernel_map(n,0,0,0,2), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get(1)
                    q_group%E_field_list(n_mol, 2) = E_field_get(2)
                    q_group%E_field_list(n_mol, 3) = E_field_get(3) 
                    n_mol = n_mol + 1

                end if
            end if
        end do

    end select  

#else

    select case(q_group%group_type)
    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            
            n_mol = n
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            k_idx = q_group%map(n,5)
            
            q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                                    mxll%Ex(i_idx, j_idx, k_idx))
            q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                                    mxll%Ey(i_idx, j_idx, k_idx))
            q_group%E_field_list(n_mol, 3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                                    mxll%Ez(i_idx, j_idx, k_idx))
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
            
            n_mol = n
            i_idx = q_group%kernel_map(n,0,0,0,3)
            j_idx = q_group%kernel_map(n,0,0,0,4)
            k_idx = q_group%kernel_map(n,0,0,0,5)
            
            q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                                    mxll%Ex(i_idx, j_idx, k_idx))
            q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                                    mxll%Ey(i_idx, j_idx, k_idx))
            q_group%E_field_list(n_mol, 3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                                    mxll%Ez(i_idx, j_idx, k_idx))
        end do
    
    end select

#endif

end subroutine send_E_3D_to_J_ranks

!###################################################################################################

subroutine send_J_to_E_1D_ranks(mxll, q_group, myrank)
    class(TMxll_1D),  intent(inout)  :: mxll
    type(TQ_Group),   intent(inout)  :: q_group
    integer       , intent(in)       :: myrank

    integer :: n, i
    integer :: i_idx
    integer :: n_mol
    logical  :: move_mol
    real(dp) :: J_field_send
    real(dp) :: J_field_get

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)


    if (myrank == 0) then
        mxll%Jx_old = mxll%Jx
    end if

#else
    mxll%Jx_old = mxll%Jx
#endif

    
#ifdef USE_MPI    
    

    select case(q_group%group_type)

    case(Q_MATERIAL)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)


            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then

                mxll%Jx(i_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                n_mol = n_mol + 1

            else

                if (q_group%map(n,1) == myrank) then
                    J_field_send = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)

                    call mpi_send(J_field_send,1,mpi_double_precision, &
                                  q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                    n_mol = n_mol + 1

                else if (q_group%map(n,2) == myrank) then
                    call mpi_recv(J_field_get,1,mpi_double_precision, &
                                  q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    mxll%Jx(i_idx) = J_field_get

                end if
            end if
        end do

    case(Q_SINGLE)

        n_mol = 1

        do n = 1, q_group%n_systems
            move_mol = .false.
            do i = -q_group%n_ker, q_group%n_ker
                
                i_idx = q_group%kernel_map(n,i,0,0,3)

                if(q_group%kernel_map(n,i,0,0,1) == myrank .and. q_group%kernel_map(n,i,0,0,2) == myrank) then
                
                    mxll%Jx(i_idx) = q_group%density*q_group%mat_kernel(i,1,1)*  &
                                     q_group%q_sys(n_mol)%dPt_dt(1)
                    move_mol = .true.
                else

                    if (q_group%kernel_map(n,i,0,0,1) == myrank) then

                        J_field_send = q_group%density*q_group%mat_kernel(i,1,1)*  &
                                       q_group%q_sys(n_mol)%dPt_dt(1)

                        call mpi_send(J_field_send,1,mpi_double_precision, &
                                    q_group%kernel_map(n,i,0,0,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)
                        move_mol = .true.

                    else if (q_group%kernel_map(n,i,0,0,2) == myrank) then
                        
                        call mpi_recv(J_field_get,1,mpi_double_precision, &
                                    q_group%kernel_map(n,i,0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)
                        mxll%Jx(i_idx) = J_field_get
                        
                    end if
                end if
            end do
            if (move_mol) n_mol = n_mol + 1
        end do
    end select

    if (myrank == 0) then
        mxll%dJx = (mxll%Jx - mxll%Jx_old)/q_group%dt
    end if

#else

    select case(q_group%group_type)

    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            n_mol = n
            i_idx = q_group%map(n,3)
            mxll%Jx(i_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
        do i = -q_group%n_ker, q_group%n_ker
            n_mol = n
            i_idx = q_group%kernel_map(n,i,0,0,3)
            mxll%Jx(i_idx) = q_group%density*q_group%mat_kernel(i,1,1)*  &
                             q_group%q_sys(n_mol)%dPt_dt(1)
        end do
        end do
    end select
    
    mxll%dJx = (mxll%Jx - mxll%Jx_old)/q_group%dt

#endif


end subroutine send_J_to_E_1D_ranks

!###################################################################################################

subroutine send_J_to_E_2D_ranks(mxll, q_group, myrank)
    class(TMxll_2D),  intent(inout) :: mxll
    type(TQ_Group) ,  intent(inout) :: q_group
    integer        , intent(in)     :: myrank
 
    integer :: n
    integer :: i, j
    integer :: i_idx
    integer :: j_idx
    integer :: n_mol
    logical  :: move_mol
    real(dp) :: J_field_send(3)
    real(dp) :: J_field_get(3)

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
#endif

    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        mxll%Jx_old = mxll%Jx
        mxll%Jy_old = mxll%Jy
    end if
    if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        mxll%Jz_old = mxll%Jz
    end if

#ifdef USE_MPI

    select case(q_group%group_type)
    case(Q_MATERIAL)
        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            
            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then
            
                if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                    mxll%Jx(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                    mxll%Jy(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
                end if
                if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                    mxll%Jz(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)
                end if
                
                n_mol = n_mol + 1

            else
            
                if (q_group%map(n,1) == myrank) then

                    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        J_field_send(1) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                        J_field_send(2) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
                    end if
                    if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        J_field_send(3) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)
                    end if

                    call mpi_send(J_field_send,3,mpi_double_precision, &
                                    q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                    n_mol = n_mol + 1

                else if (q_group%map(n,2) == myrank) then

                    call mpi_recv(J_field_get,3,mpi_double_precision, &
                                    q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        mxll%Jx(i_idx, j_idx) = J_field_get(1)
                        mxll%Jy(i_idx, j_idx) = J_field_get(2)
                    end if
                    if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        mxll%Jz(i_idx, j_idx) = J_field_get(3)
                    end if
                end if

            end if
        end do

    case(Q_SINGLE)

        n_mol = 1

        do n = 1, q_group%n_systems
            move_mol = .false.
            do j = -q_group%n_ker, q_group%n_ker
            do i = -q_group%n_ker, q_group%n_ker
                
                i_idx = q_group%kernel_map(n,i,j,0,3)
                j_idx = q_group%kernel_map(n,i,j,0,4)

                     if(q_group%kernel_map(n,i,j,0,1) == myrank .and. &
                         q_group%kernel_map(n,i,j,0,2) == myrank) then
                
                    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        mxll%Jx(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                                q_group%q_sys(n_mol)%dPt_dt(1)
                        mxll%Jy(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                                q_group%q_sys(n_mol)%dPt_dt(2)
                    end if
                    if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        mxll%Jz(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                                q_group%q_sys(n_mol)%dPt_dt(3)
                    end if
                    move_mol = .true.
                else

                    if (q_group%kernel_map(n,i,j,0,1) == myrank) then

                        if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                            J_field_send(1) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                              q_group%q_sys(n_mol)%dPt_dt(1)
                            J_field_send(2) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                              q_group%q_sys(n_mol)%dPt_dt(2)
                        end if
                        if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                            J_field_send(3) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                              q_group%q_sys(n_mol)%dPt_dt(3)
                        end if
                        call mpi_send(J_field_send,3,mpi_double_precision, &
                                    q_group%kernel_map(n,i,j,0,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)
                        move_mol = .true.

                    else if (q_group%kernel_map(n,i,j,0,2) == myrank) then
                        
                         call mpi_recv(J_field_get,3,mpi_double_precision, &
                                    q_group%kernel_map(n,i,j,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                        if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                            mxll%Jx(i_idx, j_idx) = J_field_get(1)
                            mxll%Jy(i_idx, j_idx) = J_field_get(2)
                        end if
                        if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                            mxll%Jz(i_idx, j_idx) = J_field_get(3)
                        end if
                        
                    end if
                end if
            end do
            end do
            if (move_mol) n_mol = n_mol + 1
        end do
    end select

#else

    select case(q_group%group_type)
    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            n_mol = n
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            
            if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                mxll%Jx(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                mxll%Jy(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
            end if
            if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                mxll%Jz(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)
            end if
                
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
        do j = -q_group%n_ker, q_group%n_ker
        do i = -q_group%n_ker, q_group%n_ker
            n_mol = n
            i_idx = q_group%kernel_map(n,i,j,0,3)
            j_idx = q_group%kernel_map(n,i,j,0,4)

            if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                mxll%Jx(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)* &
                                        q_group%q_sys(n_mol)%dPt_dt(1)
                mxll%Jy(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)* &
                                        q_group%q_sys(n_mol)%dPt_dt(2)
            end if
            if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                mxll%Jz(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                        q_group%q_sys(n_mol)%dPt_dt(3)
            end if
        end do
        end do
        end do
    end select

#endif

    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        mxll%dJx = (mxll%Jx - mxll%Jx_old)/q_group%dt
        mxll%dJy = (mxll%Jy - mxll%Jy_old)/q_group%dt
    else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        mxll%dJz = (mxll%Jz - mxll%Jz_old)/q_group%dt
    end if

end subroutine send_J_to_E_2D_ranks

!###################################################################################################

subroutine send_J_to_E_3D_ranks(mxll, q_group, myrank)
    class(TMxll_3D),  intent(inout) :: mxll
    type(TQ_Group) ,  intent(inout) :: q_group
    integer        , intent(in)     :: myrank
 
    integer :: n
    integer :: i, j, k
    integer :: i_idx
    integer :: j_idx
    integer :: k_idx
    integer :: n_mol
    logical  :: move_mol
    real(dp) :: J_field_send(3)
    real(dp) :: J_field_get(3)

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
#endif


    mxll%Jx_old = mxll%Jx
    mxll%Jy_old = mxll%Jy
    mxll%Jz_old = mxll%Jz

#ifdef USE_MPI

    select case(q_group%group_type)
    case(Q_MATERIAL)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            k_idx = q_group%map(n,5)
            
            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then
            
                mxll%Jx(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                mxll%Jy(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
                mxll%Jz(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)
                
                n_mol = n_mol + 1

            else
            
                if (q_group%map(n,1) == myrank) then

                    J_field_send(1) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                    J_field_send(2) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
                    J_field_send(3) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)

                    call mpi_send(J_field_send,3,mpi_double_precision, &
                                  q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                    n_mol = n_mol + 1

                else if (q_group%map(n,2) == myrank) then

                    call mpi_recv(J_field_get,3,mpi_double_precision, &
                                  q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    mxll%Jx(i_idx, j_idx, k_idx) = J_field_get(1)
                    mxll%Jy(i_idx, j_idx, k_idx) = J_field_get(2)
                    mxll%Jz(i_idx, j_idx, k_idx) = J_field_get(3)

                end if

            end if
        end do

    case(Q_SINGLE)

        n_mol = 1

        do n = 1, q_group%n_systems
            move_mol = .false.
            do k = -q_group%n_ker, q_group%n_ker
            do j = -q_group%n_ker, q_group%n_ker
            do i = -q_group%n_ker, q_group%n_ker
                
                i_idx = q_group%kernel_map(n,i,j,k,3)
                j_idx = q_group%kernel_map(n,i,j,k,4)
                k_idx = q_group%kernel_map(n,i,j,k,5)

                     if(q_group%kernel_map(n,i,j,k,1) == myrank .and. &
                         q_group%kernel_map(n,i,j,k,2) == myrank) then
                
                    mxll%Jx(i_idx, j_idx, k_idx) = q_group%density*&
                                                   q_group%mat_kernel(i,j,k)*&
                                                   q_group%q_sys(n_mol)%dPt_dt(1)
                    mxll%Jy(i_idx, j_idx, k_idx) = q_group%density*&
                                                   q_group%mat_kernel(i,j,k)*&
                                                   q_group%q_sys(n_mol)%dPt_dt(2)
                    mxll%Jz(i_idx, j_idx, k_idx) = q_group%density*&
                                                   q_group%mat_kernel(i,j,k)*&
                                                   q_group%q_sys(n_mol)%dPt_dt(3)
                    move_mol = .true.
                else

                    if (q_group%kernel_map(n,i,j,k,1) == myrank) then

                        J_field_send(1) = q_group%density*q_group%mat_kernel(i,j,k)* &
                                          q_group%q_sys(n_mol)%dPt_dt(1)
                        J_field_send(2) = q_group%density*q_group%mat_kernel(i,j,k)* &
                                          q_group%q_sys(n_mol)%dPt_dt(2)
                        J_field_send(3) = q_group%density*q_group%mat_kernel(i,j,k)* &
                                          q_group%q_sys(n_mol)%dPt_dt(3)

                        call mpi_send(J_field_send,3,mpi_double_precision, &
                                    q_group%kernel_map(n,i,j,k,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)
                        move_mol = .true.

                    else if (q_group%kernel_map(n,i,j,k,2) == myrank) then
                        
                         call mpi_recv(J_field_get,3,mpi_double_precision, &
                                    q_group%kernel_map(n,i,j,k,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                        mxll%Jx(i_idx, j_idx, k_idx) = J_field_get(1)
                        mxll%Jy(i_idx, j_idx, k_idx) = J_field_get(2)
                        mxll%Jz(i_idx, j_idx, k_idx) = J_field_get(3)
                    end if
                end if
            end do
            end do
            end do
            if (move_mol) n_mol = n_mol + 1
        end do
    end select
#else
    select case(q_group%group_type)
    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            n_mol = n
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            k_idx = q_group%map(n,5)
            
            mxll%Jx(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
            mxll%Jy(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
            mxll%Jz(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
        do k = -q_group%n_ker, q_group%n_ker
        do j = -q_group%n_ker, q_group%n_ker
        do i = -q_group%n_ker, q_group%n_ker
            n_mol = n
            i_idx = q_group%kernel_map(n,i,j,k,3)
            j_idx = q_group%kernel_map(n,i,j,k,4)
            k_idx = q_group%kernel_map(n,i,j,k,5)

            mxll%Jx(i_idx, j_idx, k_idx) = q_group%density*q_group%mat_kernel(i,j,k)*  &
                                           q_group%q_sys(n_mol)%dPt_dt(1)
            mxll%Jy(i_idx, j_idx, k_idx) = q_group%density*q_group%mat_kernel(i,j,k)*  &
                                           q_group%q_sys(n_mol)%dPt_dt(2)
            mxll%Jz(i_idx, j_idx, k_idx) = q_group%density*q_group%mat_kernel(i,j,k)*  &
                                           q_group%q_sys(n_mol)%dPt_dt(3)
        end do
        end do
        end do
        end do
    end select
#endif

    mxll%dJx = (mxll%Jx - mxll%Jx_old)/q_group%dt
    mxll%dJy = (mxll%Jy - mxll%Jy_old)/q_group%dt
    mxll%dJz = (mxll%Jz - mxll%Jz_old)/q_group%dt

end subroutine send_J_to_E_3D_ranks

!###################################################################################################

function CROSS_PRODUCT(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

end function CROSS_PRODUCT

!###################################################################################################

end module interactions_mod