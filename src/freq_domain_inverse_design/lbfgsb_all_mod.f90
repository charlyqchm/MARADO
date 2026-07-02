! The credit of this modules and subroutines corresponds to the developers of the 
! L-BFGS-B Fortran code which can be found in the following github repository:
! https://github.com/stephenbeckr/L-BFGS-B-C.git


module lbfgsb_all_mod
  implicit none
contains

  subroutine active(n, l, u, nbd, x, iwhere, iprint, &
                    prjctd, cnstnd, boxed)
  !> \brief This subroutine initializes iwhere and projects the initial x to
  !>        the feasible set if necessary.
  !>
  !> This subroutine initializes iwhere and projects the initial x to
  !> the feasible set if necessary.
  !>
  !> @param n number of parameters
  !> @param l lower bounds on parameters
  !> @param u upper bounds on parameters
  !> @param nbd indicates which bounds are present
  !> @param x position
  !> @param iwhere On entry iwhere is unspecified.<br/>
  !>               On exit: iwhere(i)=<ul><li>-1  if x(i) has no bounds</li>
  !>                                      <li> 3   if l(i)=u(i),</li>
  !>                                      <li> 0   otherwise.</li></ul>
  !>               In cauchy, iwhere is given finer gradations.
  !> @param iprint console output flag
  !> @param prjctd On exit .true. if any input x(i) had to be clipped onto its
  !>                bound (the user supplied an infeasible starting point).
  !> @param cnstnd On exit .true. if any variable has at least one bound
  !>               (any nbd(i) /= 0); .false. for purely unconstrained problems.
  !> @param boxed On exit .true. iff every variable has both lower and upper
  !>              bounds (every nbd(i) == 2); .false. otherwise.
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

  logical          prjctd, cnstnd, boxed
  integer          n, iprint, nbd(n), iwhere(n)
  double precision x(n), l(n), u(n)

  !     ************

  integer          nbdd,i
  double precision zero
  parameter        (zero=0.0d0)

  !     Initialize nbdd, prjctd, cnstnd and boxed.

  nbdd = 0
  prjctd = .false.
  cnstnd = .false.
  boxed = .true.

  !     Project the initial x to the feasible set if necessary.

  do 10 i = 1, n
     if (nbd(i) .gt. 0) then
        if (nbd(i) .le. 2 .and. x(i) .le. l(i)) then
           if (x(i) .lt. l(i)) then
              prjctd = .true.
              x(i) = l(i)
           endif
           nbdd = nbdd + 1
        else if (nbd(i) .ge. 2 .and. x(i) .ge. u(i)) then
           if (x(i) .gt. u(i)) then
              prjctd = .true.
              x(i) = u(i)
           endif
           nbdd = nbdd + 1
        endif
     endif
  10 continue

  !     Initialize iwhere and assign values to cnstnd and boxed.

  do 20 i = 1, n
     if (nbd(i) .ne. 2) boxed = .false.
     if (nbd(i) .eq. 0) then
  !                                this variable is always free
        iwhere(i) = -1

  !           otherwise set x(i)=mid(x(i), u(i), l(i)).
     else
        cnstnd = .true.
        if (nbd(i) .eq. 2 .and. u(i) - l(i) .le. zero) then
  !                   this variable is always fixed
           iwhere(i) = 3
        else
           iwhere(i) = 0
        endif
     endif
  20 continue

  if (iprint .ge. 0) then
     if (prjctd) write (6,*) &
  & 'The initial X is infeasible.  Restart with its projection.'
     if (.not. cnstnd) &
  & write (6,*) 'This problem is unconstrained.'
  endif

  if (iprint .gt. 0) write (6,1001) nbdd

  1001 format (/,'At X0 ',i9,' variables are exactly at the bounds')

  return

  end subroutine active

  subroutine bmv(m, sy, wt, col, v, p)

  !> \brief This subroutine computes the product of the 2m x 2m middle matrix
  !>        in the compact L-BFGS formula of B and a 2m vector v.
  !>
  !> This subroutine computes the product of the 2m x 2m middle matrix
  !> in the compact L-BFGS formula of B and a 2m vector v;
  !> it returns the product in p.
  !>
  !> @param m On entry m is the maximum number of variable metric corrections
  !>             used to define the limited memory matrix.<br/>
  !>          On exit m is unchanged.
  !>
  !> @param sy On entry sy specifies the matrix S'Y.<br/>
  !>           On exit sy is unchanged.
  !>
  !> @param wt On entry wt specifies the upper triangular matrix J' which is
  !>              the Cholesky factor of (thetaS'S+LD^(-1)L').<br/>
  !>           On exit wt is unchanged.
  !>
  !> @param col On entry col specifies the number of s-vectors (or y-vectors)
  !>               stored in the compact L-BFGS formula.<br/>
  !>            On exit col is unchanged.
  !>
  !> @param v On entry v specifies vector v.<br/>
  !>          On exit v is unchanged.
  !>
  !> @param p On entry p is unspecified.<br/>
  !>          On exit p is the product Mv.
  !>
  !> Historical note: this routine used to take an `info` output parameter
  !> for the LINPACK `dtrsl` triangular-solve return status. Since LAPACK's
  !> `dtrsm` replacement cannot fail on a non-singular triangular factor
  !> (and `formt` ensures `wt` is non-singular), the parameter was always
  !> 0 on exit and has been removed.


  integer m, col
  double precision sy(m, m), wt(m, m), v(2*col), p(2*col)

  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer          i,k,i2
  double precision sum,one
  parameter        (one=1.0d0)

  if (col .eq. 0) return

  !     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
  !                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].

  !       solve Jp2=v2+LD^(-1)v1.
  p(col + 1) = v(col + 1)
  do 20 i = 2, col
     i2 = col + i
     sum = 0.0d0
     do 10 k = 1, i - 1
        sum = sum + sy(i,k)*v(k)/sy(k,k)
  10 continue
     p(i2) = v(i2) + sum
  20 continue
  !     Solve the triangular system
  call dtrsm('l','u','t','n',col,1,one,wt,m,p(col+1),col)

  !       solve D^(1/2)p1=v1.
  do 30 i = 1, col
     p(i) = v(i)/sqrt(sy(i,i))
  30 continue

  !     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
  !                    [  0         J'           ] [ p2 ]   [ p2 ].

  !       solve J^Tp2=p2.
  call dtrsm('l','u','n','n',col,1,one,wt,m,p(col+1),col)

  !       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
  !                 =-D^(-1/2)p1+D^(-1)L'p2.
  do 40 i = 1, col
     p(i) = -p(i)/sqrt(sy(i,i))
  40 continue
  do 60 i = 1, col
     sum = 0.d0
     do 50 k = i + 1, col
        sum = sum + sy(k,i)*p(col+k)/sy(i,i)
  50 continue
     p(i) = p(i) + sum
  60 continue

  return

  end subroutine bmv

  subroutine cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, &
                    m, wy, ws, sy, wt, theta, col, head, p, c, wbp, &
                    v, nseg, iprint, sbgnrm, epsmch)      

  !> \brief Compute the Generalized Cauchy Point along the projected gradient direction.
  !>
  !> For given x, l, u, g (with sbgnrm > 0), and a limited memory
  !> BFGS matrix B defined in terms of matrices WY, WS, WT, and
  !> scalars head, col, and theta, this subroutine computes the
  !> generalized Cauchy point (GCP), defined as the first local
  !> minimizer of the quadratic
  !>
  !>            Q(x + s) = g's + 1/2 s'Bs
  !>
  !> along the projected gradient direction P(x-tg,l,u).
  !> The routine returns the GCP in xcp.
  !>
  !> @param n On entry n is the dimension of the problem.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param x On entry x is the starting point for the GCP computation.<br/>
  !>          On exit x is unchanged.
  !>
  !> @param l On entry l is the lower bound of x.<br/>
  !>          On exit l is unchanged.
  !>
  !> @param u On entry u is the upper bound of x.<br/>
  !>          On exit u is unchanged.
  !>
  !> @param nbd On entry nbd represents the type of bounds imposed on the
  !>               variables, and must be specified as follows:
  !>               nbd(i)=<ul><li>0 if x(i) is unbounded,</li>
  !>                          <li>1 if x(i) has only a lower bound,</li>
  !>                          <li>2 if x(i) has both lower and upper bounds, and</li>
  !>                          <li>3 if x(i) has only an upper bound.</li></ul>
  !>            On exit nbd is unchanged.
  !>
  !> @param g On entry g is the gradient of f(x). g must be a nonzero vector.<br/>
  !>          On exit g is unchanged.
  !>
  !> @param iorder iorder will be used to store the breakpoints in the piecewise
  !>               linear path and free variables encountered.<br/>
  !>               On exit,<ul><li>iorder(1),...,iorder(nleft)
  !>                               are indices of breakpoints
  !>                               which have not been encountered;</li>
  !>                           <li>iorder(nleft+1),...,iorder(nbreak)
  !>                               are indices of
  !>                               encountered breakpoints; and</li>
  !>                           <li>iorder(nfree),...,iorder(n)
  !>                               are indices of variables which
  !>                               have no bound constraits along the search direction.</li></ul>
  !>
  !> @param iwhere On entry iwhere indicates only the permanently fixed (iwhere=3)
  !>                  or free (iwhere= -1) components of x.<br/>
  !>               On exit iwhere records the status of the current x variables.
  !>                  iwhere(i)=<ul><li>-3  if x(i) is free and has bounds, but is not moved</li>
  !>                                <li> 0   if x(i) is free and has bounds, and is moved</li>
  !>                                <li> 1   if x(i) is fixed at l(i), and l(i) .ne. u(i)</li>
  !>                                <li> 2   if x(i) is fixed at u(i), and u(i) .ne. l(i)</li>
  !>                                <li> 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)</li>
  !>                                <li>-1  if x(i) is always free, i.e., it has no bounds.</li></ul>
  !>
  !> @param t working array; will be used to store the break points.
  !>
  !> @param d the Cauchy direction P(x-tg)-x
  !>
  !> @param xcp returns the GCP on exit
  !>
  !> @param m On entry m is the maximum number of variable metric corrections
  !>             used to define the limited memory matrix.<br/>
  !>          On exit m is unchanged.
  !>
  !> @param ws On entry this stores S, a set of s-vectors, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param wy On entry this stores Y, a set of y-vectors, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param sy On entry this stores S'Y, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param wt On entry this stores the
  !>              Cholesky factorization of (theta*S'S+LD^(-1)L'), that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param theta On entry theta is the scaling factor specifying B_0 = theta I.<br/>
  !>              On exit theta is unchanged.
  !>
  !> @param col On entry col is the actual number of variable metric
  !>               corrections stored so far.<br/>
  !>            On exit col is unchanged.
  !>
  !> @param head On entry head is the location of the first s-vector (or y-vector)
  !>                in S (or Y).<br/>
  !>             On exit col is unchanged.
  !>
  !> @param p will be used to store the vector p = W^(T)d.
  !>
  !> @param c will be used to store the vector c = W^(T)(xcp-x).
  !>
  !> @param wbp will be used to store the row of W corresponding
  !>              to a breakpoint.
  !>
  !> @param v working array
  !>
  !> @param nseg On exit nseg records the number of quadratic segments explored
  !>                in searching for the GCP.
  !>
  !> @param iprint variable that must be set by the user.<br/>
  !>       It controls the frequency and type of output generated:
  !>       <ul><li>iprint<0    no output is generated;</li>
  !>           <li>iprint=0    print only one line at the last iteration;</li>
  !>           <li>0<iprint<99 print also f and |proj g| every iprint iterations;</li>
  !>           <li>iprint=99   print details of every iteration except n-vectors;</li>
  !>           <li>iprint=100  print also the changes of active set and final x;</li>
  !>           <li>iprint>100  print details of every iteration including x and g;</li></ul>
  !>       When iprint > 0, the file iterate.dat will be created to
  !>                        summarize the iteration.
  !>
  !> @param sbgnrm On entry sbgnrm is the norm of the projected gradient at x.<br/>
  !>               On exit sbgnrm is unchanged.
  !>
  !> @param epsmch machine precision epsilon
  !>
  !> Historical note: this routine used to take an `info` output parameter
  !> to forward errors from the embedded `bmv` calls. Since `bmv` cannot
  !> fail under LAPACK `dtrsm`, the parameter was always 0 on exit and
  !> has been removed.

   implicit none
   integer          n, m, head, col, nseg, iprint, &
   & nbd(n), iorder(n), iwhere(n)
   double precision theta, epsmch, &
   & x(n), l(n), u(n), g(n), t(n), d(n), xcp(n), &
   & wy(n, col), ws(n, col), sy(m, m), &
   & wt(m, m), p(2*m), c(2*m), wbp(2*m), v(2*m)
  !
  !     References:
  !
  !       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
  !       memory algorithm for bound constrained optimization'',
  !       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
  !
  !       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
  !       Subroutines for Large Scale Bound Constrained Optimization''
  !       Tech. Report, NAM-11, EECS Department, Northwestern University,
  !       1994.
  !
  !       (Postscript files of these papers are available via anonymous
  !        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

   logical          xlower,xupper,bnded,skip_segment_move
   integer          i,j,col2,nfree,nbreak,pointr, &
   & ibp,nleft,ibkmin,iter
   double precision f1,f2,dt,dtm,tsum,dibp,zibp,dibp2,bkmin, &
   & tu,tl,wmc,wmp,wmw,ddot,tj,tj0,neggi,sbgnrm, &
   & f2_org
   double precision one,zero
   parameter        (one=1.0d0,zero=0.0d0)

  !     Check the status of the variables, reset iwhere(i) if necessary;
  !       compute the Cauchy direction d and the breakpoints t; initialize
  !       the derivative f1 and the vector p = W'd (for theta = 1).

  if (sbgnrm .le. zero) then
     if (iprint .ge. 0) write (6,*) 'Subgnorm = 0.  GCP = X.'
     call dcopy(n,x,1,xcp,1)
     return
  endif
  bnded = .true.
  nfree = n + 1
  nbreak = 0
  ibkmin = 0
  bkmin = zero
  col2 = 2*col
  f1 = zero
  if (iprint .ge. 99) write (6,3010)

  !     We set p to zero and build it up as we determine d.

  do 20 i = 1, col2
     p(i) = zero
  20 continue

  !     In the following loop we determine for each variable its bound
  !        status and its breakpoint, and update p accordingly.
  !        Smallest breakpoint is identified.

  do 50 i = 1, n
     neggi = -g(i)
     if (iwhere(i) .ne. 3 .and. iwhere(i) .ne. -1) then
  !             if x(i) is not a constant and has bounds,
  !             compute the difference between x(i) and its bounds.
        if (nbd(i) .le. 2) tl = x(i) - l(i)
        if (nbd(i) .ge. 2) tu = u(i) - x(i)

  !           If a variable is close enough to a bound
  !             we treat it as at bound.
        xlower = nbd(i) .le. 2 .and. tl .le. zero
        xupper = nbd(i) .ge. 2 .and. tu .le. zero

  !              reset iwhere(i).
        iwhere(i) = 0
        if (xlower) then
           if (neggi .le. zero) iwhere(i) = 1
        else if (xupper) then
           if (neggi .ge. zero) iwhere(i) = 2
        else
           if (abs(neggi) .le. zero) iwhere(i) = -3
        endif
     endif
     pointr = head
     if (iwhere(i) .ne. 0 .and. iwhere(i) .ne. -1) then
        d(i) = zero
     else
        d(i) = neggi
        f1 = f1 - neggi*neggi
  !             calculate p := p - W'e_i* (g_i).
        do 40 j = 1, col
           p(j) = p(j) +  wy(i,pointr)* neggi
           p(col + j) = p(col + j) + ws(i,pointr)*neggi
           pointr = mod(pointr,m) + 1
  40 continue
        if (nbd(i) .le. 2 .and. nbd(i) .ne. 0 &
  & .and. neggi .lt. zero) then
  !                                 x(i) + d(i) is bounded; compute t(i).
           nbreak = nbreak + 1
           iorder(nbreak) = i
           t(nbreak) = tl/(-neggi)
           if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
              bkmin = t(nbreak)
              ibkmin = nbreak
           endif
        else if (nbd(i) .ge. 2 .and. neggi .gt. zero) then
  !                                 x(i) + d(i) is bounded; compute t(i).
           nbreak = nbreak + 1
           iorder(nbreak) = i
           t(nbreak) = tu/neggi
           if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
              bkmin = t(nbreak)
              ibkmin = nbreak
           endif
        else
  !                x(i) + d(i) is not bounded.
           nfree = nfree - 1
           iorder(nfree) = i
           if (abs(neggi) .gt. zero) bnded = .false.
        endif
     endif
  50 continue

  !     The indices of the nonzero components of d are now stored
  !       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
  !       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.

  if (theta .ne. one) then
  !                   complete the initialization of p for theta not= one.
     call dscal(col,theta,p(col+1),1)
  endif

  !     Initialize GCP xcp = x.

  call dcopy(n,x,1,xcp,1)

  if (nbreak .eq. 0 .and. nfree .eq. n + 1) then
  !                  is a zero vector, return with the initial xcp as GCP.
     if (iprint .gt. 100) write (6,1010) (xcp(i), i = 1, n)
     return
  endif

  !     Initialize c = W'(xcp - x) = 0.

  do 60 j = 1, col2
     c(j) = zero
  60 continue

  !     Initialize derivative f2.

  f2 =  -theta*f1
  f2_org  =  f2
  if (col .gt. 0) then
     call bmv(m,sy,wt,col,p,v)
     f2 = f2 - ddot(col2,v,1,p,1)
  endif
  dtm = -f1/f2
  tsum = zero
  nseg = 1
  if (iprint .ge. 99) &
  & write (6,*) 'There are ',nbreak,'  breakpoints '

  !     If there are no breakpoints, locate the GCP and return.

  skip_segment_move = .false.
  if (nbreak .gt. 0) then
     nleft = nbreak
     iter = 1
     tj = zero

     !------------------- the beginning of the loop ----------------------
     do

  !        Find the next smallest breakpoint;
  !          compute dt = t(nleft) - t(nleft + 1).

        tj0 = tj
        if (iter .eq. 1) then
  !            Since we already have the smallest breakpoint we need not do
  !            heapsort yet. Often only one breakpoint is used and the
  !            cost of heapsort is avoided.
           tj = bkmin
           ibp = iorder(ibkmin)
        else
           if (iter .eq. 2) then
  !               Replace the already used smallest breakpoint with the
  !               breakpoint numbered nbreak > nlast, before heapsort call.
              if (ibkmin .ne. nbreak) then
                 t(ibkmin) = t(nbreak)
                 iorder(ibkmin) = iorder(nbreak)
              endif
  !             Update heap structure of breakpoints
  !               (if iter=2, initialize heap).
           endif
           call hpsolb(nleft,t,iorder,iter-2)
           tj = t(nleft)
           ibp = iorder(nleft)
        endif

        dt = tj - tj0

        if (dt .ne. zero .and. iprint .ge. 100) then
           write (6,4011) nseg,f1,f2
           write (6,5010) dt
           write (6,6010) dtm
        endif

  !        If a minimizer is within this interval, locate the GCP and return.
        if (dtm .lt. dt) exit

  !        Otherwise fix one variable and
  !          reset the corresponding component of d to zero.

        tsum = tsum + dt
        nleft = nleft - 1
        iter = iter + 1
        dibp = d(ibp)
        d(ibp) = zero
        if (dibp .gt. zero) then
           zibp = u(ibp) - x(ibp)
           xcp(ibp) = u(ibp)
           iwhere(ibp) = 2
        else
           zibp = l(ibp) - x(ibp)
           xcp(ibp) = l(ibp)
           iwhere(ibp) = 1
        endif
        if (iprint .ge. 100) write (6,*) 'Variable  ',ibp,'  is fixed.'
        if (nleft .eq. 0 .and. nbreak .eq. n) then
  !                                               all n variables are fixed,
  !                                                  return with xcp as GCP.
           dtm = dt
           skip_segment_move = .true.
           exit
        endif

  !        Update the derivative information.

        nseg = nseg + 1
        dibp2 = dibp**2

  !        Update f1 and f2.

  !           temporarily set f1 and f2 for col=0.
        f1 = f1 + dt*f2 + dibp2 - theta*dibp*zibp
        f2 = f2 - theta*dibp2

        if (col .gt. 0) then
  !                             update c = c + dt*p.
           call daxpy(col2,dt,p,1,c,1)

  !              choose wbp,
  !              the row of W corresponding to the breakpoint encountered.
           pointr = head
           do 70 j = 1,col
              wbp(j) = wy(ibp,pointr)
              wbp(col + j) = theta*ws(ibp,pointr)
              pointr = mod(pointr,m) + 1
  70       continue

  !              compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
           call bmv(m,sy,wt,col,wbp,v)
           wmc = ddot(col2,c,1,v,1)
           wmp = ddot(col2,p,1,v,1)
           wmw = ddot(col2,wbp,1,v,1)

  !              update p = p - dibp*wbp.
           call daxpy(col2,-dibp,wbp,1,p,1)

  !              complete updating f1 and f2 while col > 0.
           f1 = f1 + dibp*wmc
           f2 = f2 + 2.0d0*dibp*wmp - dibp2*wmw
        endif

        f2 = max(epsmch*f2_org,f2)
        if (nleft .gt. 0) then
  !                     repeat for unsearched intervals.
           dtm = -f1/f2
           cycle
        else if (bnded) then
           f1 = zero
           f2 = zero
           dtm = zero
        else
           dtm = -f1/f2
        endif

        exit
     end do
     !------------------- the end of the loop ----------------------------
  endif

  if (.not. skip_segment_move) then
     if (iprint .ge. 99) then
        write (6,*)
        write (6,*) 'GCP found in this segment'
        write (6,4010) nseg,f1,f2
        write (6,6010) dtm
     endif
     if (dtm .le. zero) dtm = zero
     tsum = tsum + dtm

  !     Move free variables (i.e., the ones w/o breakpoints) and
  !       the variables whose breakpoints haven't been reached.
     call daxpy(n,tsum,d,1,xcp,1)
  endif

  !     Update c = c + dtm*p = W'(x^c - x)
  !       which will be used in computing r = Z'(B(x^c - x) + g).

  if (col .gt. 0) call daxpy(col2,dtm,p,1,c,1)
  if (iprint .gt. 100) write (6,1010) (xcp(i),i = 1,n)
  if (iprint .ge. 99) write (6,2010)

  1010 format ('Cauchy X =  ',/,(4x,1p,6(1x,d11.4)))
  2010 format (/,'---------------- exit CAUCHY----------------------',/)
  3010 format (/,'---------------- CAUCHY entered-------------------')
  4010 format ('Piece    ',i3,' --f1, f2 at start point ',1p,2(1x,d11.4))
  4011 format (/,'Piece    ',i3,' --f1, f2 at start point ', &
  & 1p,2(1x,d11.4))
  5010 format ('Distance to the next break point =  ',1p,d11.4)
  6010 format ('Distance to the stationary point =  ',1p,d11.4)

  return

  end subroutine cauchy

  subroutine cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, &
                    theta, col, head, nfree, cnstnd)

  !> \brief This subroutine computes r=-Z'B(xcp-xk)-Z'g by using
  !>        wa(2m+1)=W'(xcp-x) from subroutine cauchy.
  !>
  !> This subroutine computes r=-Z'B(xcp-xk)-Z'g by using
  !> wa(2m+1)=W'(xcp-x) from subroutine cauchy.
  !>
  !> @param n number of parameters
  !> @param m history size of Hessian approximation
  !> @param x position
  !> @param g gradient
  !> @param ws part of L-BFGS matrix
  !> @param wy part of L-BFGS matrix
  !> @param sy part of L-BFGS matrix
  !> @param wt part of L-BFGS matrix
  !> @param z The generalized Cauchy point xcp computed by cauchy.
  !>          Used here as the linearisation point: r encodes -B(z-x) - g
  !>          restricted to the free variables.
  !> @param r On exit r(1:nfree) contains -theta*(z(k)-x(k)) - g(k) plus the
  !>          W*M^{-1}*W' correction, where k = index(i). Caller uses this as
  !>          the residual for the subspace minimisation problem.
  !>          When cnstnd=.false. and col>0, the shortcut path sets
  !>          r(1:n) = -g(:) without using the index array.
  !> @param wa Length-4m workspace shared with cauchy. On entry, the segment
  !>           wa(2m+1 : 2m+2col) holds W'(z-x) (filled by cauchy). On exit
  !>           wa(1 : 2col) holds M^{-1}*W'(z-x) from the bmv call.
  !> @param index Permutation of (1..n): index(1..nfree) lists the indices of
  !>              variables that are free at the GCP and are the active
  !>              optimisation variables here.
  !> @param theta Scaling factor specifying the initial Hessian B_0 = theta*I.
  !> @param col Number of stored (s,y) correction pairs (0 on the first
  !>            iteration; up to m thereafter).
  !> @param head Index in the cyclic WS/WY buffer of the oldest stored
  !>             correction. Used to walk the columns in chronological order.
  !> @param nfree Number of free variables; size of the subspace problem.
  !> @param cnstnd .true. if the problem has bounds; controls the shortcut
  !>               path described under @param r.
  !>
  !> Historical note: this routine used to take an `info` output parameter
  !> to forward errors from the embedded `bmv` call. Since `bmv` cannot
  !> fail under LAPACK `dtrsm`, the parameter was always 0 on exit and
  !> has been removed.

  logical          cnstnd
  integer          n, m, col, head, nfree, index(n)
  double precision theta, &
  & x(n), g(n), z(n), r(n), wa(4*m), &
  & ws(n, m), wy(n, m), sy(m, m), wt(m, m)
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer          i,j,k,pointr
  double precision a1,a2

  if (.not. cnstnd .and. col .gt. 0) then
     do 26 i = 1, n
        r(i) = -g(i)
  26 continue
  else
     do 30 i = 1, nfree
        k = index(i)
        r(i) = -theta*(z(k) - x(k)) - g(k)
  30 continue
     call bmv(m,sy,wt,col,wa(2*m+1),wa(1))
     pointr = head
     do 34 j = 1, col
        a1 = wa(j)
        a2 = theta*wa(col + j)
        do 32 i = 1, nfree
           k = index(i)
           r(i) = r(i) + wy(k,pointr)*a1 + ws(k,pointr)*a2
  32 continue
        pointr = mod(pointr,m) + 1
  34 continue
  endif

  return

  end subroutine cmprlb

  !>             function at stp.<br/>
  !>          On exit f is the value of the function at stp.
  !> @param g On initial entry g is the derivative of the function at 0.<br/>
  !>          On subsequent entries g is the derivative of the
  !>             function at stp.<br/>
  !>          On exit g is the derivative of the function at stp.
  !> @param stp On entry stp is the current estimate of a satisfactory
  !>               step. On initial entry, a positive initial estimate
  !>               must be provided.<br/>
  !>            On exit stp is the current estimate of a satisfactory step
  !>               if task = 'FG'. If task = 'CONV' then stp satisfies
  !>               the sufficient decrease and curvature condition.
  !> @param ftol On entry ftol specifies a nonnegative tolerance for the
  !>                sufficient decrease condition.<br/>
  !>             On exit ftol is unchanged.
  !> @param gtol On entry gtol specifies a nonnegative tolerance for the
  !>                curvature condition.<br/>
  !>             On exit gtol is unchanged.
  !> @param xtol On entry xtol specifies a nonnegative relative tolerance
  !>                for an acceptable step. The subroutine exits with a
  !>                warning if the relative difference between sty and stx
  !>                is less than xtol.<br/>
  !>             On exit xtol is unchanged.
  !> @param stpmin On entry stpmin is a nonnegative lower bound for the step.<br/>
  !>               On exit stpmin is unchanged.
  !> @param stpmax On entry stpmax is a nonnegative upper bound for the step.<br/>
  !>               On exit stpmax is unchanged.
  !> @param task On initial entry task must be set to 'START'.<br/>
  !>             On exit task indicates the required action:
  !>             <ul>
  !>             <li>If task(1:2) = 'FG' then evaluate the function and
  !>                 derivative at stp and call dcsrch again.</li>
  !>             <li>If task(1:4) = 'CONV' then the search is successful.</li>
  !>             <li>If task(1:4) = 'WARN' then the subroutine is not able
  !>                 to satisfy the convergence conditions. The exit value of
  !>                 stp contains the best point found during the search.</li>
  !>             <li>If task(1:5) = 'ERROR' then there is an error in the
  !>                 input arguments.</li>
  !>             </ul>
  !>             On exit with convergence, a warning or an error, the
  !>                variable task contains additional information.
  !> @param isave work array
  !> @param dsave work array

  subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax, &
                    task,isave,dsave)

  !> \brief This subroutine finds a step that satisfies a sufficient
  !>        decrease condition and a curvature condition.
  !>
  !> This subroutine finds a step that satisfies a sufficient
  !> decrease condition and a curvature condition.
  !>
  !> Each call of the subroutine updates an interval with
  !> endpoints stx and sty. The interval is initially chosen
  !> so that it contains a minimizer of the modified function
  !>
  !>       psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
  !>
  !> If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
  !> interval is chosen so that it contains a minimizer of f.
  !>
  !> The algorithm is designed to find a step that satisfies
  !> the sufficient decrease condition
  !>
  !>       f(stp) <= f(0) + ftol*stp*f'(0),
  !>
  !> and the curvature condition
  !>
  !>       abs(f'(stp)) <= gtol*abs(f'(0)).
  !>
  !> If ftol is less than gtol and if, for example, the function
  !> is bounded below, then there is always a step which satisfies
  !> both conditions.
  !>
  !> If no step can be found that satisfies both conditions, then
  !> the algorithm stops with a warning. In this case stp only
  !> satisfies the sufficient decrease condition.
  !>
  !> A typical invocation of dcsrch has the following outline:
  !>
  !> ```Fortran
  !>     task = 'START'
  !>  10 continue
  !>     call dcsrch( ... )
  !>     if (task .eq. 'FG') then
  !>       Evaluate the function and the gradient at stp
  !>     goto 10
  !>     end if
  !> ```
  !>
  !> NOTE: The user must no alter work arrays between calls.
  !>
  !> @param f On initial entry f is the value of the function at 0.<br/>
  !>          On subsequent entries f is the value of the

  character*(*) task
  integer isave(2)
  double precision f,g,stp,ftol,gtol,xtol,stpmin,stpmax
  double precision dsave(13)
  !
  !     MINPACK-1 Project. June 1983.
  !     Argonne National Laboratory.
  !     Jorge J. More' and David J. Thuente.
  !
  !     MINPACK-2 Project. October 1993.
  !     Argonne National Laboratory and University of Minnesota.
  !     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
  !
  !     **********
  double precision zero,p5,p66
  parameter(zero=0.0d0,p5=0.5d0,p66=0.66d0)
  double precision xtrapl,xtrapu
  parameter(xtrapl=1.1d0,xtrapu=4.0d0)

  logical brackt
  integer stage
  double precision finit,ftest,fm,fx,fxm,fy,fym,ginit,gtest, &
  & gm,gx,gxm,gy,gym,stx,sty,stmin,stmax,width,width1

  !     Initialization block.

   if (task(1:5) .eq. 'START') then

  !        Check the input arguments for errors.

     if (stp .lt. stpmin) task = 'ERROR: STP .LT. STPMIN'
     if (stp .gt. stpmax) task = 'ERROR: STP .GT. STPMAX'
     if (g .ge. zero) task = 'ERROR: INITIAL G .GE. ZERO'
     if (ftol .lt. zero) task = 'ERROR: FTOL .LT. ZERO'
     if (gtol .lt. zero) task = 'ERROR: GTOL .LT. ZERO'
     if (xtol .lt. zero) task = 'ERROR: XTOL .LT. ZERO'
     if (stpmin .lt. zero) task = 'ERROR: STPMIN .LT. ZERO'
     if (stpmax .lt. stpmin) task = 'ERROR: STPMAX .LT. STPMIN'

  !        Exit if there are errors on input.

     if (task(1:5) .eq. 'ERROR') return

  !        Initialize local variables.

     brackt = .false.
     stage = 1
     finit = f
     ginit = g
     gtest = ftol*ginit
     width = stpmax - stpmin
     width1 = width/p5

  !        The variables stx, fx, gx contain the values of the step,
  !        function, and derivative at the best step.
  !        The variables sty, fy, gy contain the value of the step,
  !        function, and derivative at sty.
  !        The variables stp, f, g contain the values of the step,
  !        function, and derivative at stp.

     stx = zero
     fx = finit
     gx = ginit
     sty = zero
     fy = finit
     gy = ginit
     stmin = zero
     stmax = stp + xtrapu*stp
   task = 'FG'

  else

  !        Restore local variables.

     if (isave(1) .eq. 1) then
        brackt = .true.
     else
        brackt = .false.
     endif
     stage = isave(2)
     ginit = dsave(1)
     gtest = dsave(2)
     gx = dsave(3)
     gy = dsave(4)
     finit = dsave(5)
     fx = dsave(6)
     fy = dsave(7)
     stx = dsave(8)
     sty = dsave(9)
     stmin = dsave(10)
     stmax = dsave(11)
     width = dsave(12)
     width1 = dsave(13)

  endif

  !     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
  !     algorithm enters the second stage.

  ftest = finit + stp*gtest
  if (stage .eq. 1 .and. f .le. ftest .and. g .ge. zero) &
  & stage = 2

  !     Test for warnings.

  if (brackt .and. (stp .le. stmin .or. stp .ge. stmax)) &
  & task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
  if (brackt .and. stmax - stmin .le. xtol*stmax) &
  & task = 'WARNING: XTOL TEST SATISFIED'
  if (stp .eq. stpmax .and. f .le. ftest .and. g .le. gtest) &
  & task = 'WARNING: STP = STPMAX'
  if (stp .eq. stpmin .and. (f .gt. ftest .or. g .ge. gtest)) &
  & task = 'WARNING: STP = STPMIN'

  !     Test for convergence.

  if (f .le. ftest .and. abs(g) .le. gtol*(-ginit)) &
  & task = 'CONVERGENCE'

  !     Test for termination.

  if (task(1:4) .ne. 'WARN' .and. task(1:4) .ne. 'CONV') then

  !        A modified function is used to predict the step during the
  !        first stage if a lower function value has been obtained but
  !        the decrease is not sufficient.

     if (stage .eq. 1 .and. f .le. fx .and. f .gt. ftest) then

  !           Define the modified function and derivative values.

        fm = f - stp*gtest
        fxm = fx - stx*gtest
        fym = fy - sty*gtest
        gm = g - gtest
        gxm = gx - gtest
        gym = gy - gtest

  !           Call dcstep to update stx, sty, and to compute the new step.

        call dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm, &
  & brackt,stmin,stmax)

  !           Reset the function and derivative values for f.

        fx = fxm + stx*gtest
        fy = fym + sty*gtest
        gx = gxm + gtest
        gy = gym + gtest

     else

  !          Call dcstep to update stx, sty, and to compute the new step.

       call dcstep(stx,fx,gx,sty,fy,gy,stp,f,g, &
  & brackt,stmin,stmax)

     endif

  !        Decide if a bisection step is needed.

     if (brackt) then
        if (abs(sty-stx) .ge. p66*width1) stp = stx + p5*(sty - stx)
        width1 = width
        width = abs(sty-stx)
     endif

  !        Set the minimum and maximum steps allowed for stp.

     if (brackt) then
        stmin = min(stx,sty)
        stmax = max(stx,sty)
     else
        stmin = stp + xtrapl*(stp - stx)
        stmax = stp + xtrapu*(stp - stx)
     endif

  !        Force the step to be within the bounds stpmax and stpmin.

     stp = max(stp,stpmin)
     stp = min(stp,stpmax)

  !        If further progress is not possible, let stp be the best
  !        point obtained during the search.

     if (brackt .and. (stp .le. stmin .or. stp .ge. stmax) &
  & .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx

  !        Obtain another function and derivative.

     task = 'FG'
  endif

  !     Save local variables.

  if (brackt) then
     isave(1) = 1
  else
     isave(1) = 0
  endif
  isave(2) = stage
  dsave(1) =  ginit
  dsave(2) =  gtest
  dsave(3) =  gx
  dsave(4) =  gy
  dsave(5) =  finit
  dsave(6) =  fx
  dsave(7) =  fy
  dsave(8) =  stx
  dsave(9) =  sty
  dsave(10) = stmin
  dsave(11) = stmax
  dsave(12) = width
  dsave(13) = width1

  return
  end subroutine dcsrch

  !> The subroutine assumes that if brackt is set to .true. then
  !>
  !>       min(stx,sty) < stp < max(stx,sty),
  !>
  !> and that the derivative at stx is negative in the direction
  !> of the step.
  !>
  !> @param stx On entry stx is the best step obtained so far and is an
  !>               endpoint of the interval that contains the minimizer.<br/>
  !>            On exit stx is the updated best step.
  !>
  !> @param fx On entry fx is the function at stx.<br/>
  !>           On exit fx is the function at stx.
  !>
  !> @param dx On entry dx is the derivative of the function at
  !>              stx. The derivative must be negative in the direction of
  !>              the step, that is, dx and stp - stx must have opposite
  !>              signs.<br/>
  !>           On exit dx is the derivative of the function at stx.
  !>
  !> @param sty On entry sty is the second endpoint of the interval that
  !>               contains the minimizer.<br/>
  !>            On exit sty is the updated endpoint of the interval that
  !>               contains the minimizer.
  !>
  !> @param fy On entry fy is the function at sty.<br/>
  !>           On exit fy is the function at sty.
  !>
  !> @param dy On entry dy is the derivative of the function at sty.<br/>
  !>           On exit dy is the derivative of the function at the exit sty.
  !>
  !> @param stp On entry stp is the current step. If brackt is set to .true.
  !>               then on input stp must be between stx and sty.<br/>
  !>            On exit stp is a new trial step.
  !>
  !> @param fp On entry fp is the function at stp.<br/>
  !>           On exit fp is unchanged.
  !>
  !> @param dp On entry dp is the the derivative of the function at stp.<br/>
  !>           On exit dp is unchanged.
  !>
  !> @param brackt On entry brackt specifies if a minimizer has been bracketed.
  !>                  Initially brackt must be set to .false.<br/>
  !>               On exit brackt specifies if a minimizer has been bracketed.
  !>                  When a minimizer is bracketed brackt is set to .true.
  !>
  !> @param stpmin On entry stpmin is a lower bound for the step.<br/>
  !>               On exit stpmin is unchanged.
  !>
  !> @param stpmax On entry stpmax is an upper bound for the step.<br/>
  !>               On exit stpmax is unchanged.
  subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt, &
                    stpmin,stpmax)

  !> \brief This subroutine computes a safeguarded step for a search
  !>        procedure and updates an interval that contains a step that
  !>        satisfies a sufficient decrease and a curvature condition.
  !>
  !> This subroutine computes a safeguarded step for a search
  !> procedure and updates an interval that contains a step that
  !> satisfies a sufficient decrease and a curvature condition.
  !>
  !> The parameter stx contains the step with the least function
  !> value. If brackt is set to .true. then a minimizer has
  !> been bracketed in an interval with endpoints stx and sty.
  !> The parameter stp contains the current step.

  logical brackt
  double precision stx,fx,dx,sty,fy,dy,stp,fp,dp,stpmin,stpmax
  !
  !     MINPACK-1 Project. June 1983
  !     Argonne National Laboratory.
  !     Jorge J. More' and David J. Thuente.
  !
  !     MINPACK-2 Project. October 1993.
  !     Argonne National Laboratory and University of Minnesota.
  !     Brett M. Averick and Jorge J. More'.
  !
  !     **********
  double precision zero,p66,two,three
  parameter(zero=0.0d0,p66=0.66d0,two=2.0d0,three=3.0d0)

  double precision gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta

  sgnd = dp*(dx/abs(dx))

  !     First case: A higher function value. The minimum is bracketed.
  !     If the cubic step is closer to stx than the quadratic step, the
  !     cubic step is taken, otherwise the average of the cubic and
  !     quadratic steps is taken.

  if (fp .gt. fx) then
     theta = three*(fx - fp)/(stp - stx) + dx + dp
     s = max(abs(theta),abs(dx),abs(dp))
     gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
     if (stp .lt. stx) gamma = -gamma
     p = (gamma - dx) + theta
     q = ((gamma - dx) + gamma) + dp
     r = p/q
     stpc = stx + r*(stp - stx)
     stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/two)* &
  & (stp - stx)
     if (abs(stpc-stx) .lt. abs(stpq-stx)) then
        stpf = stpc
     else
        stpf = stpc + (stpq - stpc)/two
     endif
     brackt = .true.

  !     Second case: A lower function value and derivatives of opposite
  !     sign. The minimum is bracketed. If the cubic step is farther from
  !     stp than the secant step, the cubic step is taken, otherwise the
  !     secant step is taken.

  else if (sgnd .lt. zero) then
     theta = three*(fx - fp)/(stp - stx) + dx + dp
     s = max(abs(theta),abs(dx),abs(dp))
     gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
     if (stp .gt. stx) gamma = -gamma
     p = (gamma - dp) + theta
     q = ((gamma - dp) + gamma) + dx
     r = p/q
     stpc = stp + r*(stx - stp)
     stpq = stp + (dp/(dp - dx))*(stx - stp)
     if (abs(stpc-stp) .gt. abs(stpq-stp)) then
        stpf = stpc
     else
        stpf = stpq
     endif
     brackt = .true.

  !     Third case: A lower function value, derivatives of the same sign,
  !     and the magnitude of the derivative decreases.

  else if (abs(dp) .lt. abs(dx)) then

  !        The cubic step is computed only if the cubic tends to infinity
  !        in the direction of the step or if the minimum of the cubic
  !        is beyond stp. Otherwise the cubic step is defined to be the
  !        secant step.

     theta = three*(fx - fp)/(stp - stx) + dx + dp
     s = max(abs(theta),abs(dx),abs(dp))

  !        The case gamma = 0 only arises if the cubic does not tend
  !        to infinity in the direction of the step.

     gamma = s*sqrt(max(zero,(theta/s)**2-(dx/s)*(dp/s)))
     if (stp .gt. stx) gamma = -gamma
     p = (gamma - dp) + theta
     q = (gamma + (dx - dp)) + gamma
     r = p/q
     if (r .lt. zero .and. gamma .ne. zero) then
        stpc = stp + r*(stx - stp)
     else if (stp .gt. stx) then
        stpc = stpmax
     else
        stpc = stpmin
     endif
     stpq = stp + (dp/(dp - dx))*(stx - stp)

     if (brackt) then

  !           A minimizer has been bracketed. If the cubic step is
  !           closer to stp than the secant step, the cubic step is
  !           taken, otherwise the secant step is taken.

        if (abs(stpc-stp) .lt. abs(stpq-stp)) then
           stpf = stpc
        else
           stpf = stpq
        endif
        if (stp .gt. stx) then
           stpf = min(stp+p66*(sty-stp),stpf)
        else
           stpf = max(stp+p66*(sty-stp),stpf)
        endif
     else

  !           A minimizer has not been bracketed. If the cubic step is
  !           farther from stp than the secant step, the cubic step is
  !           taken, otherwise the secant step is taken.

        if (abs(stpc-stp) .gt. abs(stpq-stp)) then
           stpf = stpc
        else
           stpf = stpq
        endif
        stpf = min(stpmax,stpf)
        stpf = max(stpmin,stpf)
     endif

  !     Fourth case: A lower function value, derivatives of the same sign,
  !     and the magnitude of the derivative does not decrease. If the
  !     minimum is not bracketed, the step is either stpmin or stpmax,
  !     otherwise the cubic step is taken.

  else
     if (brackt) then
        theta = three*(fp - fy)/(sty - stp) + dy + dp
        s = max(abs(theta),abs(dy),abs(dp))
        gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
        if (stp .gt. sty) gamma = -gamma
        p = (gamma - dp) + theta
        q = ((gamma - dp) + gamma) + dy
        r = p/q
        stpc = stp + r*(sty - stp)
        stpf = stpc
     else if (stp .gt. stx) then
        stpf = stpmax
     else
        stpf = stpmin
     endif
  endif

  !     Update the interval which contains a minimizer.

  if (fp .gt. fx) then
     sty = stp
     fy = fp
     dy = dp
  else
     if (sgnd .lt. zero) then
        sty = stx
        fy = fx
        dy = dx
     endif
     stx = stp
     fx = fp
     dx = dp
  endif

  !     Compute the new step.

  stp = stpf

  return
  end subroutine dcstep

  subroutine errclb(n, m, factr, l, u, nbd, task, info, k)

  !> \brief This subroutine checks the validity of the input data.
  !>
  !> This subroutine checks the validity of the input data.
  !>
  !> @param n number of parameters
  !> @param m history size of approximated Hessian
  !> @param factr convergence criterion on function value
  !> @param l lower bounds for parameters
  !> @param u upper bounds for parameters
  !> @param nbd indicates which bounds are present
  !> @param task if an error occurs, contains a human-readable error message
  !> @param info =0 on success; =-6 if nbd(k) was invalid; =-7 if both limits are given but l(k) > u(k)
  !> @param k index of last errournous parameter


  character*60     task
  integer          n, m, info, k, nbd(n)
  double precision factr, l(n), u(n)
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer          i
  double precision zero
  parameter        (zero=0.0d0)

  !     Check the input arguments for errors.

  if (n .le. 0) task = 'ERROR: N .LE. 0'
  if (m .le. 0) task = 'ERROR: M .LE. 0'
  if (factr .lt. zero) task = 'ERROR: FACTR .LT. 0'

  !     Check the validity of the arrays nbd(i), u(i), and l(i).

  do 10 i = 1, n
     if (nbd(i) .lt. 0 .or. nbd(i) .gt. 3) then
  !                                                   return
        task = 'ERROR: INVALID NBD'
        info = -6
        k = i
     endif
     if (nbd(i) .eq. 2) then
        if (l(i) .gt. u(i)) then
  !                                    return
           task = 'ERROR: NO FEASIBLE SOLUTION'
           info = -7
           k = i
        endif
     endif
  10 continue

  return

  end subroutine errclb

  subroutine formk(n, nsub, ind, nenter, ileave, indx2, iupdat, &
                   updatd, wn, wn1, m, ws, wy, sy, theta, col, &
                   head, info)

  !> \brief Forms the LEL^T factorization of the indefinite matrix K.
  !>
  !> This subroutine forms the LEL^T factorization of the indefinite matrix
  !>
  !> K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  !>     [L_a -R_z           theta*S'AA'S ]
  !>                                    where E = [-I  0]
  !>                                              [ 0  I]
  !>
  !> The matrix K can be shown to be equal to the matrix M^[-1]N
  !>   occurring in section 5.1 of [1], as well as to the matrix
  !>   Mbar^[-1] Nbar in section 5.3.
  !>
  !> @param n On entry n is the dimension of the problem.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param nsub On entry nsub is the number of subspace variables in free set.<br/>
  !>             On exit nsub is not changed.
  !>
  !> @param ind On entry ind specifies the indices of subspace variables.<br/>
  !>            On exit ind is unchanged.
  !>
  !> @param nenter On entry nenter is the number of variables entering the
  !>                  free set.<br/>
  !>               On exit nenter is unchanged.
  !>
  !> @param ileave On entry indx2(ileave),...,indx2(n) are the variables leaving
  !>                  the free set.<br/>
  !>               On exit ileave is unchanged.
  !>
  !> @param indx2 On entry indx2(1),...,indx2(nenter) are the variables entering
  !>                 the free set, while indx2(ileave),...,indx2(n) are the
  !>                 variables leaving the free set.<br/>
  !>              On exit indx2 is unchanged.
  !>
  !> @param iupdat On entry iupdat is the total number of BFGS updates made so far.<br/>
  !>               On exit iupdat is unchanged.
  !>
  !> @param updatd On entry 'updatd' is true if the L-BFGS matrix is updated.<br/>
  !>               On exit 'updatd' is unchanged.
  !>
  !> @param wn On entry wn is unspecified.<br/>
  !>           On exit the upper triangle of wn stores the LEL^T factorization
  !>              of the 2*col x 2*col indefinite matrix
  !>                         [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  !>                         [L_a -R_z           theta*S'AA'S ]
  !>
  !> @param wn1 On entry wn1 stores the lower triangular part of
  !>                          [Y' ZZ'Y   L_a'+R_z']
  !>                          [L_a+R_z   S'AA'S   ]
  !>               in the previous iteration.<br/>
  !>            On exit wn1 stores the corresponding updated matrices.<br/>
  !>            The purpose of wn1 is just to store these inner products
  !>            so they can be easily updated and inserted into wn.
  !>
  !> @param m On entry m is the maximum number of variable metric corrections
  !>             used to define the limited memory matrix.<br/>
  !>          On exit m is unchanged.
  !>
  !> @param ws On entry this stores S, a set of s-vectors, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param wy On entry this stores Y, a set of y-vectors, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param sy On entry this stores S'Y, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param theta On entry theta is the scaling factor specifying B_0 = theta I.<br/>
  !>              On exit theta is unchanged.
  !>
  !> @param col On entry col is the actual number of variable metric
  !>               corrections stored so far.<br/>
  !>            On exit col is unchanged.
  !>
  !> @param head On entry head is the location of the first s-vector (or y-vector)
  !>                in S (or Y).<br/>
  !>             On exit col is unchanged.
  !>
  !> @param info On entry info is unspecified.<br/>
  !>             On exit info<ul><li>=  0 for normal return;</li>
  !>                             <li>= -1 when the 1st Cholesky factorization failed;</li>
  !>                             <li>= -2 when the 2st Cholesky factorization failed.</li></ul>

  integer          n, nsub, m, col, head, nenter, ileave, iupdat, &
  & info, ind(n), indx2(n)
  double precision theta, wn(2*m, 2*m), wn1(2*m, 2*m), &
  & ws(n, m), wy(n, m), sy(m, m)
  logical          updatd
  !
  !     References:
  !       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
  !       memory algorithm for bound constrained optimization'',
  !       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
  !
  !       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
  !       limited memory FORTRAN code for solving bound constrained
  !       optimization problems'', Tech. Report, NAM-11, EECS Department,
  !       Northwestern University, 1994.
  !
  !       (Postscript files of these papers are available via anonymous
  !        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer          m2,ipntr,jpntr,iy,is,jy,js,is1,js1,k1,i,k, &
  & col2,pbegin,pend,dbegin,dend,upcl
  double precision ddot,temp1,temp2,temp3,temp4
  double precision one,zero
  parameter        (one=1.0d0,zero=0.0d0)

  !     Form the lower triangular part of
  !               WN1 = [Y' ZZ'Y   L_a'+R_z']
  !                     [L_a+R_z   S'AA'S   ]
  !        where L_a is the strictly lower triangular part of S'AA'Y
  !              R_z is the upper triangular part of S'ZZ'Y.

  if (updatd) then
     if (iupdat .gt. m) then
  !                                 shift old part of WN1.
        do 10 jy = 1, m - 1
           js = m + jy
           call dcopy(m-jy,wn1(jy+1,jy+1),1,wn1(jy,jy),1)
           call dcopy(m-jy,wn1(js+1,js+1),1,wn1(js,js),1)
           call dcopy(m-1,wn1(m+2,jy+1),1,wn1(m+1,jy),1)
  10 continue
     endif

  !          put new rows in blocks (1,1), (2,1) and (2,2).
     pbegin = 1
     pend = nsub
     dbegin = nsub + 1
     dend = n
     iy = col
     is = m + col
     ipntr = head + col - 1
     if (ipntr .gt. m) ipntr = ipntr - m
     jpntr = head
     do 20 jy = 1, col
        js = m + jy
        temp1 = zero
        temp2 = zero
        temp3 = zero
  !             compute element jy of row 'col' of Y'ZZ'Y
        do 15 k = pbegin, pend
           k1 = ind(k)
           temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
  15 continue
  !             compute elements jy of row 'col' of L_a and S'AA'S
        do 16 k = dbegin, dend
           k1 = ind(k)
           temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
           temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  16 continue
        wn1(iy,jy) = temp1
        wn1(is,js) = temp2
        wn1(is,jy) = temp3
        jpntr = mod(jpntr,m) + 1
  20 continue

  !          put new column in block (2,1).
     jy = col
     jpntr = head + col - 1
     if (jpntr .gt. m) jpntr = jpntr - m
     ipntr = head
     do 30 i = 1, col
        is = m + i
        temp3 = zero
  !             compute element i of column 'col' of R_z
        do 25 k = pbegin, pend
           k1 = ind(k)
           temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  25 continue
        ipntr = mod(ipntr,m) + 1
        wn1(is,jy) = temp3
  30 continue
     upcl = col - 1
  else
     upcl = col
  endif

  !       modify the old parts in blocks (1,1) and (2,2) due to changes
  !       in the set of free variables.
  ipntr = head
  do 45 iy = 1, upcl
     is = m + iy
     jpntr = head
     do 40 jy = 1, iy
        js = m + jy
        temp1 = zero
        temp2 = zero
        temp3 = zero
        temp4 = zero
        do 35 k = 1, nenter
           k1 = indx2(k)
           temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
           temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
  35 continue
        do 36 k = ileave, n
           k1 = indx2(k)
           temp3 = temp3 + wy(k1,ipntr)*wy(k1,jpntr)
           temp4 = temp4 + ws(k1,ipntr)*ws(k1,jpntr)
  36 continue
        wn1(iy,jy) = wn1(iy,jy) + temp1 - temp3
        wn1(is,js) = wn1(is,js) - temp2 + temp4
        jpntr = mod(jpntr,m) + 1
  40 continue
     ipntr = mod(ipntr,m) + 1
  45 continue

  !       modify the old parts in block (2,1).
  ipntr = head
  do 60 is = m + 1, m + upcl
     jpntr = head
     do 55 jy = 1, upcl
        temp1 = zero
        temp3 = zero
        do 50 k = 1, nenter
           k1 = indx2(k)
           temp1 = temp1 + ws(k1,ipntr)*wy(k1,jpntr)
  50 continue
        do 51 k = ileave, n
           k1 = indx2(k)
           temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  51 continue
     if (is .le. jy + m) then
           wn1(is,jy) = wn1(is,jy) + temp1 - temp3
        else
           wn1(is,jy) = wn1(is,jy) - temp1 + temp3
        endif
        jpntr = mod(jpntr,m) + 1
  55 continue
     ipntr = mod(ipntr,m) + 1
  60 continue

  !     Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ]
  !                                     [-L_a +R_z        S'AA'S*theta]

  m2 = 2*m
  do 70 iy = 1, col
     is = col + iy
     is1 = m + iy
     do 65 jy = 1, iy
        js = col + jy
        js1 = m + jy
        wn(jy,iy) = wn1(iy,jy)/theta
        wn(js,is) = wn1(is1,js1)*theta
  65 continue
     do 66 jy = 1, iy - 1
        wn(jy,is) = -wn1(is1,jy)
  66 continue
     do 67 jy = iy, col
        wn(jy,is) = wn1(is1,jy)
  67 continue
     wn(iy,iy) = wn(iy,iy) + sy(iy,iy)
  70 continue

  !     Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')]
  !                                    [(-L_a +R_z)L'^-1   S'AA'S*theta  ]

  !        first Cholesky factor (1,1) block of wn to get LL'
  !                          with L' stored in the upper triangle of wn.
  call dpotrf('U',col,wn,m2,info)

  if (info .ne. 0) then
     info = -1
     return
  endif
  !        then form L^-1(-L_a'+R_z') in the (1,2) block.
  col2 = 2*col
  !     TODO: can combine this loop into a single dtrsm call?
  do 71 js = col+1 ,col2
     call dtrsm('l','u','t','n',col,1,one,wn,m2,wn(1,js),col)
  71 continue

  !     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
  !        upper triangle of (2,2) block of wn.


  do 72 is = col+1, col2
     do 74 js = is, col2
           wn(is,js) = wn(is,js) + ddot(col,wn(1,is),1,wn(1,js),1)
  74 continue
  72 continue

  !     Cholesky factorization of (2,2) block of wn.

  call dpotrf('U',col,wn(col+1,col+1),m2,info)

  if (info .ne. 0) then
     info = -2
     return
  endif

  return

  end subroutine formk

  subroutine formt(m, wt, sy, ss, col, theta, info)

  !> \brief Forms the upper half of the pos. def. and symm. T.
  !>
  !> This subroutine forms the upper half of the pos. def. and symm.
  !> T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
  !> of the array wt, and performs the Cholesky factorization of T
  !> to produce J*J', with J' stored in the upper triangle of wt.
  !>
  !> @param m history size of approximated Hessian
  !> @param wt part of L-BFGS matrix
  !> @param sy part of L-BFGS matrix
  !> @param ss part of L-BFGS matrix
  !> @param col On entry col is the actual number of variable metric
  !>               corrections stored so far.<br/>
  !>            On exit col is unchanged.
  !> @param theta On entry theta is the scaling factor specifying B_0 = theta I.<br/>
  !>              On exit theta is unchanged.
  !>
  !> @param info error/success indicator


  integer          m, col, info
  double precision theta, wt(m, m), sy(m, m), ss(m, m)
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer          i,j,k,k1
  double precision ddum
  double precision zero
  parameter        (zero=0.0d0)


  !     Form the upper half of  T = theta*SS + L*D^(-1)*L',
  !        store T in the upper triangle of the array wt.

  do 52 j = 1, col
     wt(1,j) = theta*ss(1,j)
  52 continue
  do 55 i = 2, col
     do 54 j = i, col
        k1 = min(i,j) - 1
        ddum  = zero
        do 53 k = 1, k1
           ddum  = ddum + sy(i,k)*sy(j,k)/sy(k,k)
  53 continue
        wt(i,j) = ddum + theta*ss(i,j)
  54 continue
  55 continue

  !     Cholesky factorize T to J*J' with
  !        J' stored in the upper triangle of wt.

  call dpotrf('U',col,wt,m,info)
  if (info .ne. 0) then
     info = -3
  endif

  return

  end subroutine formt

  subroutine freev(n, nfree, index, nenter, ileave, indx2, &
                    iwhere, wrk, updatd, cnstnd, iprint, iter)

  !> \brief This subroutine counts the entering and leaving variables when
  !>        iter > 0, and finds the index set of free and active variables
  !>        at the GCP.
  !>
  !> This subroutine counts the entering and leaving variables when
  !> iter > 0, and finds the index set of free and active variables
  !> at the GCP.
  !>
  !> @param n number of parameters
  !> @param nfree number of free parameters, i.e., those not at their bounds
  !> @param index for i=1,...,nfree, index(i) are the indices of free variables<br/>
  !>              for i=nfree+1,...,n, index(i) are the indices of bound variables<br/>
  !>              On entry after the first iteration, index gives
  !>                the free variables at the previous iteration.<br/>
  !>              On exit it gives the free variables based on the determination
  !>                in cauchy using the array iwhere.
  !> @param nenter On exit nenter is the number of variables that entered the
  !>                free set this iteration (were active, now free at the GCP).
  !> @param ileave On exit indx2(ileave),...,indx2(n) list the variables that
  !>               left the free set this iteration. ileave starts at n+1 and
  !>               is decremented each time a leaving variable is recorded.
  !> @param indx2 On entry indx2 is unspecified.<br/>
  !>              On exit with iter>0, indx2 indicates which variables
  !>                 have changed status since the previous iteration.<br/>
  !>              For i= 1,...,nenter, indx2(i) have changed from bound to free.<br/>
  !>              For i= ileave+1,...,n, indx2(i) have changed from free to bound.
  !> @param iwhere On entry iwhere(i) classifies each variable's bound status
  !>               (set by cauchy): <=0 means free at GCP, >0 means at-bound.
  !>               Used here to compare against the previous index/nfree to
  !>               detect leaving and entering variables.
  !> @param wrk On exit .true. if the active-set or L-BFGS bookkeeping has
  !>            changed enough that the workspace WN needs to be rebuilt
  !>            (some variable entered/left, or updatd is .true.).
  !> @param updatd On entry .true. if the L-BFGS matrix was updated in the
  !>               previous iteration. Combined with the entering/leaving
  !>               counts to set wrk.
  !> @param cnstnd Whether bounds are present (true if at least one variable
  !>               is bounded). When false, the entering/leaving counting
  !>               loop is skipped.
  !> @param iprint Console output flag (>=99 prints summary, >=100 prints
  !>               per-variable change records).
  !> @param iter Current outer iteration number. The entering/leaving
  !>             counting loop only runs when iter > 0 (the first iteration
  !>             has no "previous" set to compare against).

  integer n, nfree, nenter, ileave, iprint, iter, &
  & index(n), indx2(n), iwhere(n)
  logical wrk, updatd, cnstnd
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer iact,i,k

  nenter = 0
  ileave = n + 1
  if (iter .gt. 0 .and. cnstnd) then
  !                           count the entering and leaving variables.
     do 20 i = 1, nfree
        k = index(i)

  !            write(6,*) ' k  = index(i) ', k
  !            write(6,*) ' index = ', i

        if (iwhere(k) .gt. 0) then
           ileave = ileave - 1
           indx2(ileave) = k
           if (iprint .ge. 100) write (6,*) &
  & 'Variable ',k,' leaves the set of free variables'
        endif
  20 continue
     do 22 i = 1 + nfree, n
        k = index(i)
        if (iwhere(k) .le. 0) then
           nenter = nenter + 1
           indx2(nenter) = k
           if (iprint .ge. 100) write (6,*) &
  & 'Variable ',k,' enters the set of free variables'
        endif
  22 continue
     if (iprint .ge. 99) write (6,*) &
  & n+1-ileave,' variables leave; ',nenter,' variables enter'
  endif
  wrk = (ileave .lt. n+1) .or. (nenter .gt. 0) .or. updatd

  !     Find the index set of free and active variables at the GCP.

  nfree = 0
  iact = n + 1
  do 24 i = 1, n
     if (iwhere(i) .le. 0) then
        nfree = nfree + 1
        index(nfree) = i
     else
        iact = iact - 1
        index(iact) = i
     endif
  24 continue
  if (iprint .ge. 99) write (6,*) &
  & nfree,' variables are free at GCP ',iter + 1

  return

  end subroutine freev

  subroutine hpsolb(n, t, iorder, iheap)
  !> \file hpsolb.f

  !> This subroutine sorts out the least element of t, and puts the
  !>   remaining elements of t in a heap.
  !>
  !> @param n On entry n is the dimension of the arrays t and iorder.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param t On entry t stores the elements to be sorted.<br/>
  !>          On exit t(n) stores the least elements of t, and t(1) to t(n-1)
  !>             stores the remaining elements in the form of a heap.
  !>
  !> @param iorder On entry iorder(i) is the index of t(i).<br/>
  !>               On exit iorder(i) is still the index of t(i), but iorder may be
  !>                  permuted in accordance with t.
  !>
  !> @param iheap On entry iheap should be set as follows:<ul>
  !>                 <li>iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,</li>
  !>                 <li>iheap .ne. 0 if otherwise.</li></ul>
  !>              On exit iheap is unchanged.

  integer          iheap, n, iorder(n)
  double precision t(n)
  !
  !     References:
  !       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !     ************

  integer          i,j,k,indxin,indxou
  double precision ddum,out

  if (iheap .eq. 0) then

  !        Rearrange the elements t(1) to t(n) to form a heap.

     do 20 k = 2, n
        ddum  = t(k)
        indxin = iorder(k)

  !           Add ddum to the heap.
        i = k
        do while (i .gt. 1)
           j = i/2
           if (ddum .lt. t(j)) then
              t(i) = t(j)
              iorder(i) = iorder(j)
              i = j
           else
              exit
           endif
        end do
        t(i) = ddum
        iorder(i) = indxin
  20 continue
  endif

  !     Assign to 'out' the value of t(1), the least member of the heap,
  !        and rearrange the remaining members to form a heap as
  !        elements 1 to n-1 of t.

  if (n .gt. 1) then
     i = 1
     out = t(1)
     indxou = iorder(1)
     ddum  = t(n)
     indxin  = iorder(n)

  !        Restore the heap
     do
        j = i+i
        if (j .le. n-1) then
           if (t(j+1) .lt. t(j)) j = j+1
           if (t(j) .lt. ddum ) then
              t(i) = t(j)
              iorder(i) = iorder(j)
              i = j
              cycle
           endif
        endif
        exit
     end do
     t(i) = ddum
     iorder(i) = indxin

  !     Put the least member in t(n).

     t(n) = out
     iorder(n) = indxou
  endif

  return

  end subroutine hpsolb

  subroutine lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, &
                    z, stp, dnorm, dtd, xstep, stpmx, iter, ifun, &
                    iback, nfgv, info, task, boxed, cnstnd, csave, &
                    isave, dsave)
  !> \file lnsrlb.f

  !> \brief This subroutine calls subroutine dcsrch from the Minpack2 library
  !>        to perform the line search.  Subroutine dscrch is safeguarded so
  !>        that all trial points lie within the feasible region.
  !>
  !> @param n number of parameters
  !> @param l lower bounds of parameters
  !> @param u upper bounds of parameters
  !> @param nbd On entry nbd represents the type of bounds imposed on the
  !>               variables, and must be specified as follows:
  !>               nbd(i)=<ul><li>0 if x(i) is unbounded,</li>
  !>                          <li>1 if x(i) has only a lower bound,</li>
  !>                          <li>2 if x(i) has both lower and upper bounds, and</li>
  !>                          <li>3 if x(i) has only an upper bound.</li></ul>
  !>            On exit nbd is unchanged.
  !> @param x position
  !> @param f function value at x
  !> @param fold Function value at the start of this line search (i.e. the
  !>              accepted value from the previous iteration). Saved at entry
  !>              so the caller can restore x, g, f if the line search fails.
  !> @param gd Directional derivative g'd at the current trial step. Computed
  !>           on every entry and passed to dcsrch as its `g` argument.
  !> @param gdold Directional derivative at stp=0 (i.e. the initial g'd before
  !>              any line-search progress this iteration). Saved on the first
  !>              call and used by mainlb to test the curvature condition
  !>              after the line search returns.
  !> @param g Gradient of f at x.
  !> @param d Search direction (z - x_current). Length-n vector; the candidate
  !>          step is t + stp*d.
  !> @param r Workspace: copy of g at the start of this line search. Used
  !>          alongside fold to restore the previous iterate on abnormal
  !>          line-search termination (mainlb does the restore from r and t).
  !> @param t Workspace: copy of x at the start of this line search.
  !> @param z Pre-projected candidate (cauchy/subsm output). When stp=1
  !>          exactly, lnsrlb sets x := z directly; otherwise it computes
  !>          x := t + stp*d.
  !> @param stp Current trial step length. On the first entry of a line
  !>            search lnsrlb initialises stp; subsequent dcsrch calls
  !>            update it.
  !> @param dnorm 2-norm of d (||d||).
  !> @param dtd Squared 2-norm of d (d'd).
  !> @param xstep On exit stp * ||d||, the actual length of the step in x-space.
  !> @param stpmx Maximum allowed step. For unconstrained problems set to a
  !>              large constant (1e10); for bounded problems lnsrlb scans
  !>              the active bounds and tightens stpmx so x + stpmx*d stays
  !>              feasible.
  !> @param iter Outer iteration number from mainlb.
  !> @param ifun On exit number of f/g evaluations performed in this line
  !>             search; reset to 0 on each new line search.
  !> @param iback On exit number of "backtracks" (ifun - 1). mainlb aborts
  !>              the line search if iback >= 20.
  !> @param nfgv Cumulative count of f/g evaluations across all iterations;
  !>             incremented by 1 per evaluation requested.
  !> @param info On exit 0 on success; -4 if the projected directional
  !>             derivative gd is non-negative on the first call (no descent
  !>             possible).
  !> @param task Reverse-comm task. Initial entry: 'START'. While the line
  !>             search is running, lnsrlb returns 'FG_LNSRCH' (the user
  !>             evaluates f, g at the new x and re-enters with task starting
  !>             with 'FG_LN'). On line-search success lnsrlb returns 'NEW_X'.
  !> @param boxed .true. if every variable has both lower and upper bounds.
  !>              When true, the initial trial step is unit (stp=1) regardless
  !>              of d's magnitude.
  !> @param cnstnd .true. if the problem has at least one bound. Controls the
  !>               stpmx-from-bounds scan.
  !> @param csave working array
  !> @param isave working array
  !> @param dsave working array

  character*60     task, csave
  logical          boxed, cnstnd
  integer          n, iter, ifun, iback, nfgv, info, &
  & nbd(n), isave(2)
  double precision f, fold, gd, gdold, stp, dnorm, dtd, xstep, &
  & stpmx, x(n), l(n), u(n), g(n), d(n), r(n), t(n), &
  & z(n), dsave(13)
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     **********

  integer          i
  double           precision ddot,a1,a2
  double precision one,zero,big
  parameter        (one=1.0d0,zero=0.0d0,big=1.0d+10)
  !     Wolfe-condition tolerances passed to dcsrch:
  !       ftol -- sufficient-decrease constant (alpha in eq (2.5) of the
  !               algorithm tech report; 1.0d-3 here, vs 1.0d-4 stated in
  !               that report). The looser value 1.0d-3 matches the
  !               implementation that ships with Algorithm 778; neither
  !               the ACM paper (docs/code.pdf) nor the 2011 remark
  !               documents the change explicitly.
  !       gtol -- curvature constant (beta in eq (2.6) of the tech report,
  !               0.9 -- matches the paper).
  !       xtol -- relative-step tolerance for dcsrch's bracketing safeguard.
  double precision ftol,gtol,xtol
  parameter        (ftol=1.0d-3,gtol=0.9d0,xtol=0.1d0)

  if (task(1:5) .ne. 'FG_LN') then
     dtd = ddot(n,d,1,d,1)
     dnorm = sqrt(dtd)

  !     Determine the maximum step length.

     stpmx = big
     if (cnstnd) then
        if (iter .eq. 0) then
           stpmx = one
        else
           do 43 i = 1, n
              a1 = d(i)
              if (nbd(i) .ne. 0) then
                 if (a1 .lt. zero .and. nbd(i) .le. 2) then
                    a2 = l(i) - x(i)
                    if (a2 .ge. zero) then
                       stpmx = zero
                    else if (a1*stpmx .lt. a2) then
                       stpmx = a2/a1
                    endif
                 else if (a1 .gt. zero .and. nbd(i) .ge. 2) then
                    a2 = u(i) - x(i)
                    if (a2 .le. zero) then
                       stpmx = zero
                    else if (a1*stpmx .gt. a2) then
                       stpmx = a2/a1
                    endif
                 endif
              endif
  43       continue
        endif
     endif

     if (iter .eq. 0 .and. .not. boxed) then
        stp = min(one/dnorm, stpmx)
     else
        stp = one
     endif

     call dcopy(n,x,1,t,1)
     call dcopy(n,g,1,r,1)
     fold = f
     ifun = 0
     iback = 0
     csave = 'START'
  endif
  gd = ddot(n,g,1,d,1)
  if (ifun .eq. 0) then
     gdold=gd
     if (gd .ge. zero) then
  !                               the directional derivative >=0.
  !                               Line search is impossible.
        write(6,*)' ascent direction in projection gd = ', gd
        info = -4
        return
     endif
  endif

  call dcsrch(f,gd,stp,ftol,gtol,xtol,zero,stpmx,csave,isave,dsave)

  xstep = stp*dnorm
  if (csave(1:4) .ne. 'CONV' .and. csave(1:4) .ne. 'WARN') then
     task = 'FG_LNSRCH'
     ifun = ifun + 1
     nfgv = nfgv + 1
     iback = ifun - 1
     if (stp .eq. one) then
        call dcopy(n,z,1,x,1)
     else
        do 41 i = 1, n
           x(i) = stp*d(i) + t(i)
  41 continue
     endif
  else
     task = 'NEW_X'
  endif

  return

  end subroutine lnsrlb

  !>              subroutine formk with this information.
  !>
  !> @param task working string indicating
  !>             the current job when entering and leaving this subroutine.
  !>
  !> @param iprint It controls the frequency and type of output generated:<ul>
  !>               <li>iprint<0    no output is generated;</li>
  !>               <li>iprint=0    print only one line at the last iteration;</li>
  !>               <li>0<iprint<99 print also f and |proj g| every iprint iterations;</li>
  !>               <li>iprint=99   print details of every iteration except n-vectors;</li>
  !>               <li>iprint=100  print also the changes of active set and final x;</li>
  !>               <li>iprint>100  print details of every iteration including x and g;</li></ul>
  !>               When iprint > 0, the file iterate.dat will be created to
  !>                                summarize the iteration.
  !>
  !> @param csave working string
  !>
  !> @param lsave working array
  !>
  !> @param isave working array
  !>
  !> @param dsave working array
  subroutine mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy, &
  !> \file mainlb.f

  !> \brief This subroutine solves bound constrained optimization problems by
  !>        using the compact formula of the limited memory BFGS updates.
  !>
  !> This subroutine solves bound constrained optimization problems by
  !> using the compact formula of the limited memory BFGS updates.
  !>
  !> @param n On entry n is the number of variables.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param m On entry m is the maximum number of variable metric
  !>             corrections allowed in the limited memory matrix.<br/>
  !>          On exit m is unchanged.
  !>
  !> @param x On entry x is an approximation to the solution.<br/>
  !>          On exit x is the current approximation.
  !>
  !> @param l On entry l is the lower bound of x.<br/>
  !>          On exit l is unchanged.
  !>
  !> @param u On entry u is the upper bound of x.<br/>
  !>          On exit u is unchanged.
  !>
  !> @param nbd On entry nbd represents the type of bounds imposed on the
  !>               variables, and must be specified as follows:
  !>               nbd(i)=<ul><li>0 if x(i) is unbounded,</li>
  !>                          <li>1 if x(i) has only a lower bound,</li>
  !>                          <li>2 if x(i) has both lower and upper bounds,</li>
  !>                          <li>3 if x(i) has only an upper bound.</li></ul>
  !>            On exit nbd is unchanged.
  !>
  !> @param f On first entry f is unspecified.<br/>
  !>          On final exit f is the value of the function at x.
  !>
  !> @param g On first entry g is unspecified.<br/>
  !>          On final exit g is the value of the gradient at x.
  !>
  !> @param factr On entry factr >= 0 is specified by the user.  The iteration
  !>                 will stop when<br/>
  !>                        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch<br/>
  !>                 where epsmch is the machine precision, which is automatically
  !>                 generated by the code.<br/>
  !>              On exit factr is unchanged.
  !>
  !> @param pgtol On entry pgtol >= 0 is specified by the user.  The iteration
  !>                 will stop when<br/>
  !>                        max{|proj g_i | i = 1, ..., n} <= pgtol<br/>
  !>                 where pg_i is the ith component of the projected gradient.<br/>
  !>              On exit pgtol is unchanged.
  !>
  !> @param ws On entry this stores S, a set of s-vectors, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param wy On entry this stores Y, a set of y-vectors, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param sy On entry this stores S'Y, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param ss On entry this stores S'S, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param wt On entry this stores the
  !>              Cholesky factorization of (theta*S'S+LD^(-1)L'), that defines the
  !>              limited memory BFGS matrix. See eq. (2.26) in [3].<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param wn working array used to store the LEL^T factorization of the indefinite matrix
  !>                     K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  !>                         [L_a -R_z           theta*S'AA'S ]<br/>
  !>           where     E = [-I  0]
  !>                         [ 0  I]
  !>
  !> @param snd working array used to store the lower triangular part of
  !>                      N = [Y' ZZ'Y   L_a'+R_z']
  !>                          [L_a +R_z  S'AA'S   ]
  !>
  !> @param z working array used at different times to store the Cauchy point and
  !>          the Newton point.
  !> @param r working array
  !> @param d working array
  !> @param t working array
  !> @param xp working array used to safeguard the projected Newton direction
  !> @param wa working array
  !>
  !> @param index In subroutine freev, index is used to store the free and fixed
  !>              variables at the Generalized Cauchy Point (GCP).
  !>
  !> @param iwhere working array used to record
  !>               the status of the vector x for GCP computation.<br/>
  !>               iwhere(i)=<ul><li>0 or -3 if x(i) is free and has bounds,</li>
  !>                             <li> 1       if x(i) is fixed at l(i), and l(i) .ne. u(i)</li>
  !>                             <li> 2       if x(i) is fixed at u(i), and u(i) .ne. l(i)</li>
  !>                             <li> 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)</li>
  !>                             <li>-1       if x(i) is always free, i.e., no bounds on it.</li></ul>
  !>
  !> @param indx2 working array<br/>
  !>              Within subroutine cauchy, indx2 corresponds to the array iorder.<br/>
  !>              In subroutine freev, a list of variables entering and leaving
  !>              the free set is stored in indx2, and it is passed on to

  & sy, ss, wt, wn, snd, z, r, d, t, xp, wa, &
  & index, iwhere, indx2, task, &
  & iprint, csave, lsave, isave, dsave)
  implicit none
  character*60     task, csave
  logical          lsave(4)
  integer          n, m, iprint, nbd(n), index(n), &
  & iwhere(n), indx2(n), isave(23)
  double precision f, factr, pgtol, &
  & x(n), l(n), u(n), g(n), z(n), r(n), d(n), t(n), &
  !-jlm-jn
  & xp(n), &
  & wa(8*m), &
  & ws(n, m), wy(n, m), sy(m, m), ss(m, m), &
  & wt(m, m), wn(2*m, 2*m), snd(2*m, 2*m), dsave(29)
  !
  !     References:
  !
  !       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
  !       memory algorithm for bound constrained optimization'',
  !       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
  !
  !       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
  !       Subroutines for Large Scale Bound Constrained Optimization''
  !       Tech. Report, NAM-11, EECS Department, Northwestern University,
  !       1994.
  !
  !       [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
  !       Quasi-Newton Matrices and their use in Limited Memory Methods'',
  !       Mathematical Programming 63 (1994), no. 4, pp. 129-156.
  !
  !       (Postscript files of these papers are available via anonymous
  !        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  logical          prjctd,cnstnd,boxed,updatd,wrk
  character*3      word
  integer          i,k,nintol,itfile,iback,nskip, &
  & head,col,iter,itail,iupdat, &
  & nseg,nfgv,info,ifun, &
  & iword,nfree,nact,ileave,nenter
  double precision theta,fold,ddot,dr,rr,tol, &
  & xstep,sbgnrm,ddum,dnorm,dtd,epsmch, &
  & cpu1,cpu2,cachyt,sbtime,lnscht,time1,time2, &
  & gd,gdold,stp,stpmx,time
  double precision one,zero
  parameter        (one=1.0d0,zero=0.0d0)

  if (task .eq. 'START') then

     epsmch = epsilon(one)

     call timer(time1)

  !        Initialize counters and scalars when task='START'.

  !           for the limited memory BFGS matrices:
     col    = 0
     head   = 1
     theta  = one
     iupdat = 0
     updatd = .false.
     iback  = 0
     itail  = 0
     iword  = 0
     nact   = 0
     ileave = 0
     nenter = 0
     fold   = zero
     dnorm  = zero
     cpu1   = zero
     gd     = zero
     stpmx  = zero
     sbgnrm = zero
     stp    = zero
     gdold  = zero
     dtd    = zero

  !           for operation counts:
     iter   = 0
     nfgv   = 0
     nseg   = 0
     nintol = 0
     nskip  = 0
     nfree  = n
     ifun   = 0
  !           for stopping tolerance:
     tol = factr*epsmch

  !           for measuring running time:
     cachyt = 0
     sbtime = 0
     lnscht = 0

  !           'word' records the status of subspace solutions.
     word = '---'

  !           'info' records the termination information.
     info = 0

     itfile = 8
     if (iprint .ge. 1) then
  !                                open a summary file 'iterate.dat'
        open (itfile, file = 'iterate.dat', status = 'unknown')
     endif

  !        Check the input arguments for errors.

     call errclb(n,m,factr,l,u,nbd,task,info,k)
     if (task(1:5) .eq. 'ERROR') then
        call prn3lb(n,x,f,task,iprint,info,itfile, &
  & iter,nfgv,nintol,nskip,nact,sbgnrm, &
  & zero,nseg,word,iback,stp,xstep,k, &
  & cachyt,sbtime,lnscht)
        return
     endif

     call prn1lb(n,m,l,u,x,iprint,itfile,epsmch)

  !        Initialize iwhere & project x onto the feasible set.

     call active(n,l,u,nbd,x,iwhere,iprint,prjctd,cnstnd,boxed)

  !        The end of the initialization.

  else
  !          restore local variables.

     prjctd = lsave(1)
     cnstnd = lsave(2)
     boxed  = lsave(3)
     updatd = lsave(4)

     nintol = isave(1)
     itfile = isave(3)
     iback  = isave(4)
     nskip  = isave(5)
     head   = isave(6)
     col    = isave(7)
     itail  = isave(8)
     iter   = isave(9)
     iupdat = isave(10)
     nseg   = isave(12)
     nfgv   = isave(13)
     info   = isave(14)
     ifun   = isave(15)
     iword  = isave(16)
     nfree  = isave(17)
     nact   = isave(18)
     ileave = isave(19)
     nenter = isave(20)

     theta  = dsave(1)
     fold   = dsave(2)
     tol    = dsave(3)
     dnorm  = dsave(4)
     epsmch = dsave(5)
     cpu1   = dsave(6)
     cachyt = dsave(7)
     sbtime = dsave(8)
     lnscht = dsave(9)
     time1  = dsave(10)
     gd     = dsave(11)
     stpmx  = dsave(12)
     sbgnrm = dsave(13)
     stp    = dsave(14)
     gdold  = dsave(15)
     dtd    = dsave(16)

  !        After returning from the driver go to the point where execution
  !        is to resume.

     if (task(1:5) .eq. 'FG_LN') goto 666
     if (task(1:5) .eq. 'NEW_X') goto 777
     if (task(1:5) .eq. 'FG_ST') goto 111
     if (task(1:4) .eq. 'STOP') then
        if (task(7:9) .eq. 'CPU') then
  !                                          restore the previous iterate.
           call dcopy(n,t,1,x,1)
           call dcopy(n,r,1,g,1)
           f = fold
        endif
        goto 999
     endif
  endif

  !     Compute f0 and g0.

  task = 'FG_START'
  !          return to the driver to calculate f and g; reenter at 111.
  goto 1000
  111 continue
  nfgv = 1

  !     Compute the infinity norm of the (-) projected gradient.

  call projgr(n,l,u,nbd,x,g,sbgnrm)

  if (iprint .ge. 1) then
     write (6,1002) iter,f,sbgnrm
     write (itfile,1003) iter,nfgv,sbgnrm,f
  endif
  if (sbgnrm .le. pgtol) then
  !                                terminate the algorithm.
     task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
     goto 999
  endif

  ! ----------------- the beginning of the loop --------------------------

  222 continue
  if (iprint .ge. 99) write (6,1001) iter + 1
  iword = -1
  !
  if (.not. cnstnd .and. col .gt. 0) then
  !                                            skip the search for GCP.
     call dcopy(n,x,1,z,1)
     wrk = updatd
     nseg = 0
     goto 333
  endif

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     Compute the Generalized Cauchy Point (GCP).
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  call timer(cpu1)
  call cauchy(n,x,l,u,nbd,g,indx2,iwhere,t,d,z, &
  & m,wy,ws,sy,wt,theta,col,head, &
  & wa(1),wa(2*m+1),wa(4*m+1),wa(6*m+1),nseg, &
  & iprint, sbgnrm, epsmch)
  call timer(cpu2)
  cachyt = cachyt + cpu2 - cpu1
  nintol = nintol + nseg

  !     Count the entering and leaving variables for iter > 0;
  !     find the index set of free and active variables at the GCP.

  call freev(n,nfree,index,nenter,ileave,indx2, &
  & iwhere,wrk,updatd,cnstnd,iprint,iter)
  nact = n - nfree

  333 continue

  !     If there are no free variables or B=theta*I, then
  !                                        skip the subspace minimization.

  if (nfree .eq. 0 .or. col .eq. 0) goto 555

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     Subspace minimization.
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  call timer(cpu1)

  !     Form  the LEL^T factorization of the indefinite
  !       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  !                     [L_a -R_z           theta*S'AA'S ]
  !       where     E = [-I  0]
  !                     [ 0  I]

  if (wrk) call formk(n,nfree,index,nenter,ileave,indx2,iupdat, &
  & updatd,wn,snd,m,ws,wy,sy,theta,col,head,info)
  if (info .ne. 0) then
  !          nonpositive definiteness in Cholesky factorization;
  !          refresh the lbfgs memory and restart the iteration.
     if(iprint .ge. 1) write (6, 1006)
     info   = 0
     col    = 0
     head   = 1
     theta  = one
     iupdat = 0
     updatd = .false.
     call timer(cpu2)
     sbtime = sbtime + cpu2 - cpu1
     goto 222
  endif

  !        compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x)
  !                                                   from 'cauchy').
  call cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa,index, &
  & theta,col,head,nfree,cnstnd)

  !-jlm-jn   call the direct method.

  call subsm( n, m, nfree, index, l, u, nbd, z, r, xp, ws, wy, &
  & theta, x, g, col, head, iword, wa, wn, iprint)

  call timer(cpu2)
  sbtime = sbtime + cpu2 - cpu1
  555 continue

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     Line search and optimality tests.
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  !     Generate the search direction d:=z-x.

  do 40 i = 1, n
     d(i) = z(i) - x(i)
  40 continue
  call timer(cpu1)
  666 continue
  call lnsrlb(n,l,u,nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm, &
  & dtd,xstep,stpmx,iter,ifun,iback,nfgv,info,task, &
  & boxed,cnstnd,csave,isave(22),dsave(17))
  if (info .ne. 0 .or. iback .ge. 20) then
  !          restore the previous iterate.
     call dcopy(n,t,1,x,1)
     call dcopy(n,r,1,g,1)
     f = fold
     if (col .eq. 0) then
  !             abnormal termination.
        if (info .eq. 0) then
           info = -9
  !                restore the actual number of f and g evaluations etc.
           nfgv = nfgv - 1
           ifun = ifun - 1
           iback = iback - 1
        endif
        task = 'ABNORMAL_TERMINATION_IN_LNSRCH'
        iter = iter + 1
        goto 999
     else
  !             refresh the lbfgs memory and restart the iteration.
        if(iprint .ge. 1) write (6, 1008)
        if (info .eq. 0) nfgv = nfgv - 1
        info   = 0
        col    = 0
        head   = 1
        theta  = one
        iupdat = 0
        updatd = .false.
        task   = 'RESTART_FROM_LNSRCH'
        call timer(cpu2)
        lnscht = lnscht + cpu2 - cpu1
        goto 222
     endif
  else if (task(1:5) .eq. 'FG_LN') then
  !          return to the driver for calculating f and g; reenter at 666.
     goto 1000
  else
  !          calculate and print out the quantities related to the new X.
     call timer(cpu2)
     lnscht = lnscht + cpu2 - cpu1
     iter = iter + 1

  !        Compute the infinity norm of the projected (-)gradient.

     call projgr(n,l,u,nbd,x,g,sbgnrm)

  !        Print iteration information.

     call prn2lb(n,x,f,g,iprint,itfile,iter,nfgv,nact, &
  & sbgnrm,nseg,word,iword,iback,stp,xstep)
     goto 1000
  endif
  777 continue

  !     Test for termination.

  if (sbgnrm .le. pgtol) then
  !                                terminate the algorithm.
     task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
     goto 999
  endif

  ddum = max(abs(fold), abs(f), one)
  if ((fold - f) .le. tol*ddum) then
  !                                        terminate the algorithm.
     task = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
     if (iback .ge. 10) info = -5
  !           i.e., to issue a warning if iback>10 in the line search.
     goto 999
  endif

  !     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.

  do 42 i = 1, n
     r(i) = g(i) - r(i)
  42 continue
  rr = ddot(n,r,1,r,1)
  if (stp .eq. one) then
     dr = gd - gdold
     ddum = -gdold
  else
     dr = (gd - gdold)*stp
     call dscal(n,stp,d,1)
     ddum = -gdold*stp
  endif

  if (dr .le. epsmch*ddum) then
  !                            skip the L-BFGS update.
     nskip = nskip + 1
     updatd = .false.
     if (iprint .ge. 1) write (6,1004) dr, ddum
     goto 888
  endif

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     Update the L-BFGS matrix.
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  updatd = .true.
  iupdat = iupdat + 1

  !     Update matrices WS and WY and form the middle matrix in B.

  call matupd(n,m,ws,wy,sy,ss,d,r,itail, &
  & iupdat,col,head,theta,rr,dr,stp,dtd)

  !     Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
  !        Store T in the upper triangular of the array wt;
  !        Cholesky factorize T to J*J' with
  !           J' stored in the upper triangular of wt.

  call formt(m,wt,sy,ss,col,theta,info)

  if (info .ne. 0) then
  !          nonpositive definiteness in Cholesky factorization;
  !          refresh the lbfgs memory and restart the iteration.
     if(iprint .ge. 1) write (6, 1007)
     info = 0
     col = 0
     head = 1
     theta = one
     iupdat = 0
     updatd = .false.
     goto 222
  endif

  !     Now the inverse of the middle matrix in B is

  !       [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
  !       [ -L*D^(-1/2)   J ] [  0        J'          ]

  888 continue

  ! -------------------- the end of the loop -----------------------------

  goto 222
  999 continue
  call timer(time2)
  time = time2 - time1
  call prn3lb(n,x,f,task,iprint,info,itfile, &
  & iter,nfgv,nintol,nskip,nact,sbgnrm, &
  & time,nseg,word,iback,stp,xstep,k, &
  & cachyt,sbtime,lnscht)
  1000 continue

  !     Save local variables.

  lsave(1)  = prjctd
  lsave(2)  = cnstnd
  lsave(3)  = boxed
  lsave(4)  = updatd

  isave(1)  = nintol
  isave(3)  = itfile
  isave(4)  = iback
  isave(5)  = nskip
  isave(6)  = head
  isave(7)  = col
  isave(8)  = itail
  isave(9)  = iter
  isave(10) = iupdat
  isave(12) = nseg
  isave(13) = nfgv
  isave(14) = info
  isave(15) = ifun
  isave(16) = iword
  isave(17) = nfree
  isave(18) = nact
  isave(19) = ileave
  isave(20) = nenter

  dsave(1)  = theta
  dsave(2)  = fold
  dsave(3)  = tol
  dsave(4)  = dnorm
  dsave(5)  = epsmch
  dsave(6)  = cpu1
  dsave(7)  = cachyt
  dsave(8)  = sbtime
  dsave(9)  = lnscht
  dsave(10) = time1
  dsave(11) = gd
  dsave(12) = stpmx
  dsave(13) = sbgnrm
  dsave(14) = stp
  dsave(15) = gdold
  dsave(16) = dtd

  1001 format (//,'ITERATION ',i5)
  1002 format &
  & (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
  1003 format (2(1x,i4),5x,'-',5x,'-',3x,'-',5x,'-',5x,'-',8x,'-',3x, &
  & 1p,2(1x,d10.3))
  1004 format ('  ys=',1p,e10.3,'  -gs=',1p,e10.3,' BFGS update SKIPPED')
  1006 format (/, &
  & ' Nonpositive definiteness in Cholesky factorization in formk;',/, &
  & '   refresh the lbfgs memory and restart the iteration.')
  1007 format (/, &
  & ' Nonpositive definiteness in Cholesky factorization in formt;',/, &
  & '   refresh the lbfgs memory and restart the iteration.')
  1008 format (/, &
  & ' Bad direction in the line search;',/, &
  & '   refresh the lbfgs memory and restart the iteration.')

  return

  end subroutine mainlb

  subroutine matupd(n, m, ws, wy, sy, ss, d, r, itail, &
  !> \file matupd.f

  !> \brief This subroutine updates matrices WS and WY, and forms the
  !>        middle matrix in B.
  !>
  !> This subroutine updates matrices WS and WY, and forms the
  !> middle matrix in B.
  !>
  !> @param n On entry n is the number of variables.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param m On entry m is the maximum number of variable metric
  !>             corrections allowed in the limited memory matrix.<br/>
  !>          On exit m is unchanged.
  !>
  !> @param ws On entry this stores S, a set of s-vectors, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param wy On entry this stores Y, a set of y-vectors, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param sy On entry this stores S'Y, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param ss On entry this stores S'S, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !> @param d Search direction at the current iteration. After the line
  !>           search, the new s-vector is stp*d; matupd stores it as a
  !>           column of WS (writing d directly, since the stp scaling is
  !>           folded into the stored ss diagonal entry below).
  !> @param r The accepted gradient difference y = g_{k+1} - g_k. Stored as
  !>          a column of WY.
  !> @param itail On entry the column index in WS/WY that the previous
  !>              update wrote. On exit the column index this update wrote;
  !>              advances cyclically modulo m.
  !> @param iupdat Total number of L-BFGS updates performed so far
  !>               (incremented by mainlb before each matupd call). When
  !>               iupdat <= m, col simply grows; when iupdat > m, the
  !>               history wraps and the oldest column is discarded.
  !> @param col On entry col is the actual number of variable metric
  !>               corrections stored so far.<br/>
  !>            On exit col is unchanged.
  !>
  !> @param head On entry head is the location of the first s-vector (or y-vector)
  !>                in S (or Y).<br/>
  !>             On exit col is unchanged.
  !>
  !> @param theta On entry theta is the scaling factor specifying B_0 = theta I.<br/>
  !>              On exit theta is unchanged.
  !> @param rr Squared 2-norm of r (i.e. ||y||^2 = y'y). Used to set
  !>           theta := rr/dr = y'y / (s'y), the standard initial Hessian
  !>           scaling for the next L-BFGS update.
  !> @param dr Inner product d'r = s'y / stp (the curvature condition; must
  !>           be positive for the update to be safely accepted -- mainlb
  !>           checks this and skips matupd if dr <= eps*rr).
  !> @param stp Line-search step length. Used to recover the s-vector from
  !>            d (s = stp*d) when computing ss(col,col) = stp^2 * d'd.
  !> @param dtd Squared 2-norm of d (i.e. d'd). Combined with stp to form
  !>            ss(col,col) = stp^2 * dtd = ||s||^2 for the new column.

  & iupdat, col, head, theta, rr, dr, stp, dtd)

  integer          n, m, itail, iupdat, col, head
  double precision theta, rr, dr, stp, dtd, d(n), r(n), &
  & ws(n, m), wy(n, m), sy(m, m), ss(m, m)

  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer          j,pointr
  double precision ddot
  double precision one
  parameter        (one=1.0d0)

  !     Set pointers for matrices WS and WY.

  if (iupdat .le. m) then
     col = iupdat
     itail = mod(head+iupdat-2,m) + 1
  else
     itail = mod(itail,m) + 1
     head = mod(head,m) + 1
  endif

  !     Update matrices WS and WY.

  call dcopy(n,d,1,ws(1,itail),1)
  call dcopy(n,r,1,wy(1,itail),1)

  !     Set theta=yy/ys.

  theta = rr/dr

  !     Form the middle matrix in B.

  !        update the upper triangle of SS,
  !                                         and the lower triangle of SY:
  if (iupdat .gt. m) then
  !                              move old information
     do 50 j = 1, col - 1
        call dcopy(j,ss(2,j+1),1,ss(1,j),1)
        call dcopy(col-j,sy(j+1,j+1),1,sy(j,j),1)
  50 continue
  endif
  !        add new information: the last row of SY
  !                                             and the last column of SS:
  pointr = head
  do 51 j = 1, col - 1
     sy(col,j) = ddot(n,d,1,wy(1,pointr),1)
     ss(j,col) = ddot(n,ws(1,pointr),1,d,1)
     pointr = mod(pointr,m) + 1
  51 continue
  if (stp .eq. one) then
     ss(col,col) = dtd
  else
     ss(col,col) = stp*stp*dtd
  endif
  sy(col,col) = dr

  return

  end subroutine matupd

  subroutine prn1lb(n, m, l, u, x, iprint, itfile, epsmch)
  !> \file prn1lb.f

  !> \brief This subroutine prints the input data, initial point, upper and
  !>        lower bounds of each variable, machine precision, as well as
  !>        the headings of the output.
  !>
  !> This subroutine prints the input data, initial point, upper and
  !>        lower bounds of each variable, machine precision, as well as
  !>        the headings of the output.
  !>
  !> @param n On entry n is the number of variables.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param m On entry m is the maximum number of variable metric
  !>             corrections allowed in the limited memory matrix.<br/>
  !>          On exit m is unchanged.
  !>
  !> @param l On entry l is the lower bound of x.<br/>
  !>          On exit l is unchanged.
  !>
  !> @param u On entry u is the upper bound of x.<br/>
  !>          On exit u is unchanged.
  !>
  !> @param x On entry x is an approximation to the solution.<br/>
  !>          On exit x is the current approximation.
  !>
  !> @param iprint It controls the frequency and type of output generated:<ul>
  !>               <li>iprint<0    no output is generated;</li>
  !>               <li>iprint=0    print only one line at the last iteration;</li>
  !>               <li>0<iprint<99 print also f and |proj g| every iprint iterations;</li>
  !>               <li>iprint=99   print details of every iteration except n-vectors;</li>
  !>               <li>iprint=100  print also the changes of active set and final x;</li>
  !>               <li>iprint>100  print details of every iteration including x and g;</li></ul>
  !>               When iprint > 0, the file iterate.dat will be created to
  !>                                summarize the iteration.
  !>
  !> @param itfile unit number of iterate.dat file
  !> @param epsmch machine precision epsilon


  integer n, m, iprint, itfile
  double precision epsmch, x(n), l(n), u(n)

  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer i

  if (iprint .ge. 0) then
     write (6,7001) epsmch
     write (6,*) 'N = ',n,'    M = ',m
     if (iprint .ge. 1) then
        write (itfile,2001) epsmch
        write (itfile,*)'N = ',n,'    M = ',m
        write (itfile,9001)
        if (iprint .gt. 100) then
           write (6,1004) 'L =',(l(i),i = 1,n)
           write (6,1004) 'X0 =',(x(i),i = 1,n)
           write (6,1004) 'U =',(u(i),i = 1,n)
        endif
     endif
  endif

  1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
  2001 format ('RUNNING THE L-BFGS-B CODE',/,/, &
  & 'it    = iteration number',/, &
  & 'nf    = number of function evaluations',/, &
  & 'nseg  = number of segments explored during the Cauchy search',/, &
  & 'nact  = number of active bounds at the generalized Cauchy point' &
  & ,/, &
  & 'sub   = manner in which the subspace minimization terminated:' &
  & ,/,'        con = converged, bnd = a bound was reached',/, &
  & 'itls  = number of iterations performed in the line search',/, &
  & 'stepl = step length used',/, &
  & 'tstep = norm of the displacement (total step)',/, &
  & 'projg = norm of the projected gradient',/, &
  & 'f     = function value',/,/, &
  & '           * * *',/,/, &
  & 'Machine precision =',1p,d10.3)
  7001 format ('RUNNING THE L-BFGS-B CODE',/,/, &
  & '           * * *',/,/, &
  & 'Machine precision =',1p,d10.3)
  9001 format (/,3x,'it',3x,'nf',2x,'nseg',2x,'nact',2x,'sub',2x,'itls', &
  & 2x,'stepl',4x,'tstep',5x,'projg',8x,'f')

  return

  end subroutine prn1lb

  subroutine prn2lb(n, x, f, g, iprint, itfile, iter, nfgv, nact, &
  !> \file prn2lb.f

  !> \brief This subroutine prints out new information after a successful
  !>        line search.
  !>
  !> This subroutine prints out new information after a successful
  !> line search.
  !>
  !> @param n On entry n is the number of variables.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param x On entry x is an approximation to the solution.<br/>
  !>          On exit x is the current approximation.
  !>
  !> @param f On first entry f is unspecified.<br/>
  !>          On final exit f is the value of the function at x.
  !>
  !> @param g On first entry g is unspecified.<br/>
  !>          On final exit g is the value of the gradient at x.
  !>
  !> @param iprint It controls the frequency and type of output generated:<ul>
  !>               <li>iprint<0    no output is generated;</li>
  !>               <li>iprint=0    print only one line at the last iteration;</li>
  !>               <li>0<iprint<99 print also f and |proj g| every iprint iterations;</li>
  !>               <li>iprint=99   print details of every iteration except n-vectors;</li>
  !>               <li>iprint=100  print also the changes of active set and final x;</li>
  !>               <li>iprint>100  print details of every iteration including x and g;</li></ul>
  !>               When iprint > 0, the file iterate.dat will be created to
  !>                                summarize the iteration.
  !>
  !> @param itfile unit number of iterate.dat file
  !>
  !> @param iter Current outer iteration number; printed on the per-iterate
  !>              summary line.
  !> @param nfgv Cumulative number of f/g evaluations across the whole run;
  !>             logged into the iterate.dat record.
  !> @param nact Number of variables active at their bounds at the current
  !>             iterate.
  !> @param sbgnrm Infinity norm of the projected gradient at x; the natural
  !>               first-order optimality measure printed alongside f.
  !> @param nseg Number of breakpoint segments traversed by cauchy at the
  !>             current iteration; printed for diagnostic accounting.
  !> @param word On exit a 3-character status code summarising the subspace
  !>             solution: 'con' (converged), 'bnd' (hit a bound), 'TNT'
  !>             (truncated-Newton step used), or '---' otherwise. Determined
  !>             from iword.
  !> @param iword Status code from subsm: 0 = subspace minimisation
  !>              converged, 1 = stopped at a bound, 5 = truncated-Newton
  !>              step used. Maps to word.
  !> @param iback Number of backtracks the line search performed.
  !> @param stp Final step length accepted by the line search.
  !> @param xstep stp * ||d|| -- the actual length of the step in x-space.

  & sbgnrm, nseg, word, iword, iback, stp, xstep)

  character*3      word
  integer          n, iprint, itfile, iter, nfgv, nact, nseg, &
  & iword, iback
  double precision f, sbgnrm, stp, xstep, x(n), g(n)

  !     ************
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer i,imod

  !           'word' records the status of subspace solutions.
  if (iword .eq. 0) then
  !                            the subspace minimization converged.
     word = 'con'
  else if (iword .eq. 1) then
  !                          the subspace minimization stopped at a bound.
     word = 'bnd'
  else if (iword .eq. 5) then
  !                             the truncated Newton step has been used.
     word = 'TNT'
  else
     word = '---'
  endif
  if (iprint .ge. 99) then
     write (6,*) 'LINE SEARCH',iback,' times; norm of step = ',xstep
     write (6,2001) iter,f,sbgnrm
     if (iprint .gt. 100) then
        write (6,1004) 'X =',(x(i), i = 1, n)
        write (6,1004) 'G =',(g(i), i = 1, n)
     endif
  else if (iprint .gt. 0) then
     imod = mod(iter,iprint)
     if (imod .eq. 0) write (6,2001) iter,f,sbgnrm
  endif
  if (iprint .ge. 1) write (itfile,3001) &
  & iter,nfgv,nseg,nact,word,iback,stp,xstep,sbgnrm,f

  1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
  2001 format &
  & (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
  3001 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),1p,2(1x,d10.3))

  return

  end subroutine prn2lb

  subroutine prn3lb(n, x, f, task, iprint, info, itfile, &
                    iter, nfgv, nintol, nskip, nact, sbgnrm, &
                    time, nseg, word, iback, stp, xstep, k, &
                    cachyt, sbtime, lnscht)
  !> \file prn3lb.f

  !> \brief This subroutine prints out information when either a built-in
  !>        convergence test is satisfied or when an error message is
  !>        generated.
  !>
  !> This subroutine prints out information when either a built-in
  !> convergence test is satisfied or when an error message is
  !> generated.
  !>
  !> @param n On entry n is the number of variables.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param x On entry x is an approximation to the solution.<br/>
  !>          On exit x is the current approximation.
  !>
  !> @param f On first entry f is unspecified.<br/>
  !>          On final exit f is the value of the function at x.
  !>
  !> @param task working string indicating
  !>             the current job when entering and leaving this subroutine.
  !> @param iprint Console output flag; same convention as in setulb.
  !>                Negative suppresses all output. >=0 prints the per-run
  !>                summary. >=1 also writes to iterate.dat. >=100 also
  !>                prints the final x vector.
  !> @param info Termination status. 0 on normal convergence. Negative
  !>             values map to specific error messages: -1, -2 (Cholesky
  !>             not pos.def. in formk's 1st/2nd factor), -3 (formt's
  !>             Cholesky), -4 (non-descent direction caught by lnsrlb),
  !>             -5 (>10 evals in a line search), -6 (invalid nbd(k)),
  !>             -7 (l(k) > u(k)), -8 (singular triangular system), -9
  !>             (>20 evals in the last line search).
  !> @param itfile unit number of iterate.dat file
  !> @param iter Total number of outer iterations.
  !> @param nfgv Total number of f/g evaluations across the run.
  !> @param nintol Total number of breakpoint segments traversed by cauchy.
  !> @param nskip Number of L-BFGS updates that were skipped because the
  !>              curvature condition failed (s'y too small).
  !> @param nact Number of active bounds at the final generalised Cauchy point.
  !> @param sbgnrm Final infinity norm of the projected gradient at x.
  !> @param time Total wall-clock time of the optimisation run, in seconds.
  !> @param nseg Number of breakpoint segments traversed at the LAST cauchy
  !>             call (used in the iterate.dat record on abnormal exits).
  !> @param word 3-character status code from prn2lb: 'con', 'bnd', 'TNT',
  !>             or '---' (see prn2lb.f). Logged in iterate.dat on info=-4
  !>             or info=-9 termination.
  !> @param iback Number of backtracks in the LAST line search.
  !> @param stp Final step length from the LAST line search.
  !> @param xstep stp * ||d|| from the LAST line search (i.e. the actual
  !>              step length in x-space).
  !> @param k Index of the offending parameter when info is -6 (invalid
  !>          nbd) or -7 (infeasible bound). Set by errclb.
  !> @param cachyt Cumulative wall-clock time spent in cauchy, in seconds.
  !> @param sbtime Cumulative wall-clock time spent in subsm, in seconds.
  !> @param lnscht Cumulative wall-clock time spent in lnsrlb, in seconds.

  character*60     task
  character*3      word
  integer          n, iprint, info, itfile, iter, nfgv, nintol, &
  & nskip, nact, nseg, iback, k
  double precision f, sbgnrm, time, stp, xstep, cachyt, sbtime, &
  & lnscht, x(n)

  !     ************
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer i

  if (task(1:5) .ne. 'ERROR') then

     if (iprint .ge. 0) then
        write (6,3003)
        write (6,3004)
        write(6,3005) n,iter,nfgv,nintol,nskip,nact,sbgnrm,f
        if (iprint .ge. 100) then
           write (6,1004) 'X =',(x(i),i = 1,n)
        endif
        if (iprint .ge. 1) write (6,*) ' F =',f
     endif
  endif
  if (iprint .ge. 0) then
     write (6,3009) task
     if (info .ne. 0) then
        if (info .eq. -1) write (6,9011)
        if (info .eq. -2) write (6,9012)
        if (info .eq. -3) write (6,9013)
        if (info .eq. -4) write (6,9014)
        if (info .eq. -5) write (6,9015)
        if (info .eq. -6) write (6,*)' Input nbd(',k,') is invalid.'
        if (info .eq. -7) &
  & write (6,*)' l(',k,') > u(',k,').  No feasible solution.'
        if (info .eq. -8) write (6,9018)
        if (info .eq. -9) write (6,9019)
     endif
     if (iprint .ge. 1) write (6,3007) cachyt,sbtime,lnscht
     write (6,3008) time
     if (iprint .ge. 1) then
        if (info .eq. -4 .or. info .eq. -9) then
           write (itfile,3002) &
  & iter,nfgv,nseg,nact,word,iback,stp,xstep
        endif
        write (itfile,3009) task
        if (info .ne. 0) then
           if (info .eq. -1) write (itfile,9011)
           if (info .eq. -2) write (itfile,9012)
           if (info .eq. -3) write (itfile,9013)
           if (info .eq. -4) write (itfile,9014)
           if (info .eq. -5) write (itfile,9015)
           if (info .eq. -8) write (itfile,9018)
           if (info .eq. -9) write (itfile,9019)
        endif
        write (itfile,3008) time
     endif
  endif

  1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
  3002 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),6x,'-',10x,'-')
  3003 format (/, &
  & '           * * *',/,/, &
  & 'Tit   = total number of iterations',/, &
  & 'Tnf   = total number of function evaluations',/, &
  & 'Tnint = total number of segments explored during', &
  & ' Cauchy searches',/, &
  & 'Skip  = number of BFGS updates skipped',/, &
  & 'Nact  = number of active bounds at final generalized', &
  & ' Cauchy point',/, &
  & 'Projg = norm of the final projected gradient',/, &
  & 'F     = final function value',/,/, &
  & '           * * *')
  3004 format (/,3x,'N',4x,'Tit',5x,'Tnf',2x,'Tnint',2x, &
  & 'Skip',2x,'Nact',5x,'Projg',8x,'F')
  3005 format (i5,2(1x,i6),(1x,i6),(2x,i4),(1x,i5),1p,2(2x,d10.3))
  3007 format (/,' Cauchy                time',1p,e10.3,' seconds.',/ &
  & ' Subspace minimization time',1p,e10.3,' seconds.',/ &
  & ' Line search           time',1p,e10.3,' seconds.')
  3008 format (/,' Total User time',1p,e10.3,' seconds.',/)
  3009 format (/,a60)
  9011 format (/, &
  & ' Matrix in 1st Cholesky factorization in formk is not Pos. Def.')
  9012 format (/, &
  & ' Matrix in 2st Cholesky factorization in formk is not Pos. Def.')
  9013 format (/, &
  & ' Matrix in the Cholesky factorization in formt is not Pos. Def.')
  9014 format (/, &
  & ' Derivative >= 0, backtracking line search impossible.',/, &
  & '   Previous x, f and g restored.',/, &
  & ' Possible causes: 1 error in function or gradient evaluation;',/, &
  & '                  2 rounding errors dominate computation.')
  9015 format (/, &
  & ' Warning:  more than 10 function and gradient',/, &
  & '   evaluations in the last line search.  Termination',/, &
  & '   may possibly be caused by a bad search direction.')
  9018 format (/,' The triangular system is singular.')
  9019 format (/, &
  & ' Line search cannot locate an adequate point after 20 function',/ &
  & ,'  and gradient evaluations.  Previous x, f and g restored.',/, &
  & ' Possible causes: 1 error in function or gradient evaluation;',/, &
  & '                  2 rounding errors dominate computation.')

  return

  end subroutine prn3lb

  subroutine projgr(n, l, u, nbd, x, g, sbgnrm)
  !> \file projgr.f

  !> \brief This subroutine computes the infinity norm of the projected
  !>        gradient.
  !>
  !> This subroutine computes the infinity norm of the projected
  !> gradient.
  !>
  !> @param n On entry n is the number of variables.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param l On entry l is the lower bound of x.<br/>
  !>          On exit l is unchanged.
  !>
  !> @param u On entry u is the upper bound of x.<br/>
  !>          On exit u is unchanged.
  !>
  !> @param nbd On entry nbd represents the type of bounds imposed on the
  !>               variables, and must be specified as follows:
  !>               nbd(i)=<ul><li>0 if x(i) is unbounded,</li>
  !>                          <li>1 if x(i) has only a lower bound,</li>
  !>                          <li>2 if x(i) has both lower and upper bounds,</li>
  !>                          <li>3 if x(i) has only an upper bound.</li></ul>
  !>            On exit nbd is unchanged.
  !>
  !> @param x On entry x is an approximation to the solution.<br/>
  !>          On exit x is unchanged.
  !>
  !> @param g On entry g is the gradient.<br/>
  !>          On exit g is unchanged.
  !>
  !> @param sbgnrm infinity norm of projected gradient


  integer          n, nbd(n)
  double precision sbgnrm, x(n), l(n), u(n), g(n)

  !     ************
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer i
  double precision gi
  double precision zero
  parameter        (zero=0.0d0)

  sbgnrm = zero
  do 15 i = 1, n
    gi = g(i)
    if (nbd(i) .ne. 0) then
       if (gi .lt. zero) then
          if (nbd(i) .ge. 2) gi = max((x(i)-u(i)),gi)
       else
          if (nbd(i) .le. 2) gi = min((x(i)-l(i)),gi)
       endif
    endif
    sbgnrm = max(sbgnrm,abs(gi))
  15 continue

  return

  end subroutine projgr

  subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
                    task, iprint, csave, lsave, isave, dsave)

  !> \brief This subroutine partitions the working arrays wa and iwa, and
  !>        then uses the limited memory BFGS method to solve the bound
  !>        constrained optimization problem by calling mainlb.
  !>
  !> This subroutine partitions the working arrays wa and iwa, and
  !> then uses the limited memory BFGS method to solve the bound
  !> constrained optimization problem by calling mainlb.
  !> (The direct method will be used in the subspace minimization.)
  !>
  !> @param n On entry n is the dimension of the problem.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param m On entry m is the maximum number of variable metric corrections
  !>             used to define the limited memory matrix.<br/>
  !>          On exit m is unchanged.
  !>
  !> @param x On entry x is an approximation to the solution.<br/>
  !>          On exit x is the current approximation.
  !>
  !> @param l On entry l is the lower bound on x.<br/>
  !>          On exit l is unchanged.
  !>
  !> @param u On entry u is the upper bound on x.<br/>
  !>          On exit u is unchanged.
  !>
  !> @param nbd On entry nbd represents the type of bounds imposed on the
  !>               variables, and must be specified as follows:
  !>               nbd(i)=<ul><li>0 if x(i) is unbounded,</li>
  !>                          <li>1 if x(i) has only a lower bound,</li>
  !>                          <li>2 if x(i) has both lower and upper bounds, and</li>
  !>                          <li>3 if x(i) has only an upper bound.</li></ul>
  !>            On exit nbd is unchanged.
  !>
  !> @param f On first entry f is unspecified.<br/>
  !>          On final exit f is the value of the function at x.
  !>
  !> @param g On first entry g is unspecified.<br/>
  !>          On final exit g is the value of the gradient at x.
  !>
  !> @param factr On entry factr >= 0 is specified by the user.  The iteration
  !>                 will stop when<br/>
  !>                 (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch<br/>
  !>                 where epsmch is the machine precision, which is automatically
  !>                 generated by the code.<br/>
  !>                 Typical values for factr:<ul>
  !>                    <li>1.d+12 for low accuracy;</li>
  !>                    <li>1.d+7  for moderate accuracy;</li>
  !>                    <li>1.d+1 for extremely high accuracy.</li></ul>
  !>              On exit factr is unchanged.
  !>
  !> @param pgtol On entry pgtol >= 0 is specified by the user.  The iteration
  !>                 will stop when<br/>
  !>                        max{|proj g_i | i = 1, ..., n} <= pgtol<br/>
  !>                 where pg_i is the ith component of the projected gradient.<br/>
  !>              On exit pgtol is unchanged.
  !>
  !> @param wa working array
  !>
  !> @param iwa working array
  !>
  !> @param task working string indicating
  !>             the current job when entering and quitting this subroutine.
  !>
  !> @param iprint Must be set by the user.
  !>               It controls the frequency and type of output generated:<ul>
  !>               <li>iprint<0    no output is generated;</li>
  !>               <li>iprint=0    print only one line at the last iteration;</li>
  !>               <li>0<iprint<99 print also f and |proj g| every iprint iterations;</li>
  !>               <li>iprint=99   print details of every iteration except n-vectors;</li>
  !>               <li>iprint=100  print also the changes of active set and final x;</li>
  !>               <li>iprint>100  print details of every iteration including x and g;</li></ul>
  !>               When iprint > 0, the file iterate.dat will be created to
  !>                                summarize the iteration.
  !>
  !> @param csave working string
  !>
  !> @param lsave working array;
  !>              On exit with 'task' = NEW_X, the following information is
  !>                                                                    available:<ul>
  !>                <li>If lsave(1) = .true.  then  the initial X has been replaced by
  !>                                                its projection in the feasible set;</li>
  !>                <li>If lsave(2) = .true.  then  the problem is constrained;</li>
  !>                <li>If lsave(3) = .true.  then  each variable has upper and lower
  !>                                                bounds;</li></ul>
  !>
  !> @param isave working array;
  !>              On exit with 'task' = NEW_X, the following information is
  !>                                                                    available:<ul>
  !>                <li>isave(22) = the total number of intervals explored in the
  !>                                    search of Cauchy points;</li>
  !>                <li>isave(26) = the total number of skipped BFGS updates before
  !>                                    the current iteration;</li>
  !>                <li>isave(30) = the number of current iteration;</li>
  !>                <li>isave(31) = the total number of BFGS updates prior the current
  !>                                    iteration;</li>
  !>                <li>isave(33) = the number of intervals explored in the search of
  !>                                    Cauchy point in the current iteration;</li>
  !>                <li>isave(34) = the total number of function and gradient
  !>                                    evaluations;</li>
  !>                <li>isave(36) = the number of function value or gradient
  !>                                              evaluations in the current iteration;</li>
  !>                <li>if isave(37) = 0  then the subspace argmin is within the box;</li>
  !>                <li>if isave(37) = 1  then the subspace argmin is beyond the box;</li>
  !>                <li>isave(38) = the number of free variables in the current
  !>                                    iteration;</li>
  !>                <li>isave(39) = the number of active constraints in the current
  !>                                    iteration;</li>
  !>                <li>n + 1 - isave(40) = the number of variables leaving the set of
  !>                                        active constraints in the current iteration;</li>
  !>                <li>isave(41) = the number of variables entering the set of active
  !>                                    constraints in the current iteration.</li></ul>
  !>
  !> @param dsave working array;
  !>              On exit with 'task' = NEW_X, the following information is
  !>                                                                    available:<ul>
  !>                <li>dsave(1) = current 'theta' in the BFGS matrix;</li>
  !>                <li>dsave(2) = f(x) in the previous iteration;</li>
  !>                <li>dsave(3) = factr*epsmch;</li>
  !>                <li>dsave(4) = 2-norm of the line search direction vector;</li>
  !>                <li>dsave(5) = the machine precision epsmch generated by the code;</li>
  !>                <li>dsave(7) = the accumulated time spent on searching for
  !>                                                             Cauchy points;</li>
  !>                <li>dsave(8) = the accumulated time spent on
  !>                                                        subspace minimization;</li>
  !>                <li>dsave(9) = the accumulated time spent on line search;</li>
  !>                <li>dsave(11) = the slope of the line search function at
  !>                                         the current point of line search;</li>
  !>                <li>dsave(12) = the maximum relative step length imposed in
  !>                                                                  line search;</li>
  !>                <li>dsave(13) = the infinity norm of the projected gradient;</li>
  !>                <li>dsave(14) = the relative step length in the line search;</li>
  !>                <li>dsave(15) = the slope of the line search function at
  !>                                        the starting point of the line search;</li>
  !>                <li>dsave(16) = the square of the 2-norm of the line search
  !>                                                             direction vector.</li></ul>

  character*60     task, csave
  logical          lsave(4)
  integer          n, m, iprint, &
  & nbd(n), iwa(3*n), isave(44)
  double precision f, factr, pgtol, x(n), l(n), u(n), g(n), &
  !
  !-jlm-jn
  & wa(2*m*n + 5*n + 11*m*m + 8*m), dsave(29)

  !     ************
  !
  !     References:
  !
  !       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
  !       memory algorithm for bound constrained optimization'',
  !       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
  !
  !       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
  !       limited memory FORTRAN code for solving bound constrained
  !       optimization problems'', Tech. Report, NAM-11, EECS Department,
  !       Northwestern University, 1994.
  !
  !       (Postscript files of these papers are available via anonymous
  !        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************
  !-jlm-jn
  integer   lws,lr,lz,lt,ld,lxp,lwa, &
  & lwy,lsy,lss,lwt,lwn,lsnd

  if (task .eq. 'START') then
     isave(1)  = m*n
     isave(2)  = m**2
     isave(3)  = 4*m**2
     isave(4)  = 1                      ! ws      m*n
     isave(5)  = isave(4)  + isave(1)   ! wy      m*n
     isave(6)  = isave(5)  + isave(1)   ! wsy     m**2
     isave(7)  = isave(6)  + isave(2)   ! wss     m**2
     isave(8)  = isave(7)  + isave(2)   ! wt      m**2
     isave(9)  = isave(8)  + isave(2)   ! wn      4*m**2
     isave(10) = isave(9)  + isave(3)   ! wsnd    4*m**2
     isave(11) = isave(10) + isave(3)   ! wz      n
     isave(12) = isave(11) + n          ! wr      n
     isave(13) = isave(12) + n          ! wd      n
     isave(14) = isave(13) + n          ! wt      n
     isave(15) = isave(14) + n          ! wxp     n
     isave(16) = isave(15) + n          ! wa      8*m
  endif
  lws  = isave(4)
  lwy  = isave(5)
  lsy  = isave(6)
  lss  = isave(7)
  lwt  = isave(8)
  lwn  = isave(9)
  lsnd = isave(10)
  lz   = isave(11)
  lr   = isave(12)
  ld   = isave(13)
  lt   = isave(14)
  lxp  = isave(15)
  lwa  = isave(16)

  call mainlb(n,m,x,l,u,nbd,f,g,factr,pgtol, &
  & wa(lws),wa(lwy),wa(lsy),wa(lss), wa(lwt), &
  & wa(lwn),wa(lsnd),wa(lz),wa(lr),wa(ld),wa(lt),wa(lxp), &
  & wa(lwa), &
  & iwa(1),iwa(n+1),iwa(2*n+1),task,iprint, &
  & csave,lsave,isave(22),dsave)

  return

  end subroutine setulb



  subroutine subsm ( n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy, &
                    theta, xx, gg, col, head, iword, wv, wn, iprint )

  !> \brief Performs the subspace minimization.
  !>
  !> Given xcp, l, u, r, an index set that specifies
  !> the active set at xcp, and an l-BFGS matrix B
  !> (in terms of WY, WS, SY, WT, head, col, and theta),
  !> this subroutine computes an approximate solution
  !> of the subspace problem
  !>
  !> (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)
  !>
  !> subject to l<=x<=u
  !>           x_i=xcp_i for all i in A(xcp)
  !>
  !> along the subspace unconstrained Newton direction
  !>
  !>    d = -(Z'BZ)^(-1) r.
  !>
  !> The formula for the Newton direction, given the L-BFGS matrix
  !> and the Sherman-Morrison formula, is
  !>
  !>    d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.
  !>
  !> where
  !>           K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  !>               [L_a -R_z           theta*S'AA'S ]
  !>
  !> Note that this procedure for computing d differs
  !> from that described in [1]. One can show that the matrix K is
  !> equal to the matrix M^[-1]N in that paper.
  !>
  !> @param n On entry n is the dimension of the problem.<br/>
  !>          On exit n is unchanged.
  !>
  !> @param m On entry m is the maximum number of variable metric corrections
  !>             used to define the limited memory matrix.<br/>
  !>          On exit m is unchanged.
  !>
  !> @param nsub On entry nsub is the number of free variables.<br/>
  !>             On exit nsub is unchanged.
  !>
  !> @param ind On entry ind specifies the coordinate indices of free variables.<br/>
  !>            On exit ind is unchanged.
  !>
  !> @param l On entry l is the lower bound of x.<br/>
  !>          On exit l is unchanged.
  !>
  !> @param u On entry u is the upper bound of x.<br/>
  !>          On exit u is unchanged.
  !>
  !> @param nbd On entry nbd represents the type of bounds imposed on the
  !>               variables, and must be specified as follows:
  !>               nbd(i)=<ul><li>0 if x(i) is unbounded,</li>
  !>                          <li>1 if x(i) has only a lower bound,</li>
  !>                          <li>2 if x(i) has both lower and upper bounds, and</li>
  !>                          <li>3 if x(i) has only an upper bound.</li></ul>
  !>            On exit nbd is unchanged.
  !>
  !> @param x On entry x specifies the Cauchy point xcp.<br/>
  !>          On exit x(i) is the minimizer of Q over the subspace of
  !>                                                        free variables.
  !>
  !> @param d On entry d is the reduced gradient of Q at xcp.<br/>
  !>          On exit d is the Newton direction of Q.
  !>
  !> @param xp used to safeguard the projected Newton direction<br/>
  !>
  !> @param xx On entry it holds the current iterate.<br/>
  !>           On output it is unchanged.
  !>
  !> @param gg On entry it holds the gradient at the current iterate.<br/>
  !>           On output it is unchanged.
  !>
  !> @param ws On entry this stores S, a set of s-vectors, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param wy On entry this stores Y, a set of y-vectors, that defines the
  !>              limited memory BFGS matrix.<br/>
  !>           On exit this array is unchanged.
  !>
  !> @param theta On entry theta is the scaling factor specifying B_0 = theta I.<br/>
  !>              On exit theta is unchanged.
  !>
  !> @param col On entry col is the actual number of variable metric
  !>               corrections stored so far.<br/>
  !>            On exit col is unchanged.
  !>
  !> @param head On entry head is the location of the first s-vector (or y-vector)
  !>                in S (or Y).<br/>
  !>             On exit col is unchanged.
  !>
  !> @param iword On entry iword is unspecified.<br/>
  !>              On exit iword specifies the status of the subspace solution.
  !>                 iword = <ul><li>0 if the solution is in the box,</li>
  !>                             <li>1 if some bound is encountered.</li></ul>
  !>
  !> @param wv working array
  !>
  !> @param wn On entry the upper triangle of wn stores the LEL^T factorization
  !>              of the indefinite matrix<br/>
  !>                   K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
  !>                       [L_a -R_z           theta*S'AA'S ]<br/>
  !>              where E = [-I  0]
  !>                        [ 0  I]<br/>
  !>           On exit wn is unchanged.
  !>
  !> @param iprint must be set by the user;
  !>               It controls the frequency and type of output generated:<ul>
  !>               <li>iprint<0    no output is generated;</li>
  !>               <li>iprint=0    print only one line at the last iteration;</li>
  !>               <li>0<iprint<99 print also f and |proj g| every iprint iterations;</li>
  !>               <li>iprint=99   print details of every iteration except n-vectors;</li>
  !>               <li>iprint=100  print also the changes of active set and final x;</li>
  !>               <li>iprint>100  print details of every iteration including x and g;</li></ul>
  !>               When iprint > 0, the file iterate.dat will be created to
  !>                                summarize the iteration.
  !>
  !>
  !> Historical note: this routine used to take an `info` output parameter
  !> for an "ill-conditioned K" status. Since K's conditioning is checked
  !> in `formt`/`formk` (which fail loudly with their own `info`) and the
  !> two `dtrsm` calls inside `subsm` cannot fail on a non-singular factor,
  !> the parameter was always 0 on exit and has been removed.

  implicit none
  integer          n, m, nsub, col, head, iword, iprint, &
  & ind(nsub), nbd(n)
  double precision theta, &
  & l(n), u(n), x(n), d(n), xp(n), xx(n), gg(n), &
  & ws(n, m), wy(n, m), &
  & wv(2*m), wn(2*m, 2*m)

  !     **********************************************************************
  !
  !     This routine contains the major changes in the updated version.
  !     The changes are described in the accompanying paper
  !
  !      Jose Luis Morales, Jorge Nocedal
  !      "Remark On Algorithm 788: L-BFGS-B: Fortran Subroutines for Large-Scale
  !       Bound Constrained Optimization". Decemmber 27, 2010.
  !
  !             J.L. Morales  Departamento de Matematicas,
  !                           Instituto Tecnologico Autonomo de Mexico
  !                           Mexico D.F.
  !
  !             J, Nocedal    Department of Electrical Engineering and
  !                           Computer Science.
  !                           Northwestern University. Evanston, IL. USA
  !
  !                           January 17, 2011
  !
  !      **********************************************************************
  !
  !     References:
  !
  !       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
  !       memory algorithm for bound constrained optimization'',
  !       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
  !
  !                           *  *  *
  !
  !     NEOS, November 1994. (Latest revision June 1996.)
  !     Optimization Technology Center.
  !     Argonne National Laboratory and Northwestern University.
  !     Written by
  !                        Ciyou Zhu
  !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
  !
  !
  !     ************

  integer          pointr,m2,col2,ibd,jy,js,i,j,k
  double precision alpha, xk, dk, temp1, temp2
  double precision one,zero
  parameter        (one=1.0d0,zero=0.0d0)
  !
  double precision dd_p

  if (nsub .le. 0) return
  if (iprint .ge. 99) write (6,1001)

  !     Compute wv = W'Zd.

  pointr = head
  do 20 i = 1, col
     temp1 = zero
     temp2 = zero
     do 10 j = 1, nsub
        k = ind(j)
        temp1 = temp1 + wy(k,pointr)*d(j)
        temp2 = temp2 + ws(k,pointr)*d(j)
  10 continue
     wv(i) = temp1
     wv(col + i) = theta*temp2
     pointr = mod(pointr,m) + 1
  20 continue

  !     Compute wv:=K^(-1)wv.

  m2 = 2*m
  col2 = 2*col

  call dtrsm('l','u','t','n',col2,1,one,wn,m2,wv,col2)

  do 25 i = 1, col
     wv(i) = -wv(i)
  25 continue

  call dtrsm('l','u','n','n',col2,1,one,wn,m2,wv,col2)

  !     Compute d = (1/theta)d + (1/theta**2)Z'W wv.

  pointr = head
  do 40 jy = 1, col
     js = col + jy
     do 30 i = 1, nsub
        k = ind(i)
        d(i) = d(i) + wy(k,pointr)*wv(jy)/theta &
  & + ws(k,pointr)*wv(js)
  30 continue
     pointr = mod(pointr,m) + 1
  40 continue

  call dscal( nsub, one/theta, d, 1 )
  !
  !-----------------------------------------------------------------
  !     Let us try the projection, d is the Newton direction

  iword = 0

  call dcopy ( n, x, 1, xp, 1 )
  !
  do 50 i=1, nsub
     k  = ind(i)
     dk = d(i)
     xk = x(k)
     if ( nbd(k) .ne. 0 ) then
  !
        if ( nbd(k).eq.1 ) then          ! lower bounds only
           x(k) = max( l(k), xk + dk )
           if ( x(k).eq.l(k) ) iword = 1
        else
  !
           if ( nbd(k).eq.2 ) then       ! upper and lower bounds
              xk   = max( l(k), xk + dk )
              x(k) = min( u(k), xk )
              if ( x(k).eq.l(k) .or. x(k).eq.u(k) ) iword = 1
           else
  !
              if ( nbd(k).eq.3 ) then    ! upper bounds only
                 x(k) = min( u(k), xk + dk )
                 if ( x(k).eq.u(k) ) iword = 1
              end if
           end if
        end if
  !
     else                                ! free variables
        x(k) = xk + dk
     end if
  50 continue
  !
  if ( iword.ne.0 ) then
  !
  !     check sign of the directional derivative
  !
     dd_p = zero
     do 55 i=1, n
        dd_p  = dd_p + (x(i) - xx(i))*gg(i)
  55 continue
     if ( dd_p .gt.zero ) then
        call dcopy( n, xp, 1, x, 1 )
        write(6,*) ' Positive dir derivative in projection '
        write(6,*) ' Using the backtracking step '

  !
  !-----------------------------------------------------------------
  !
        alpha = one
        temp1 = alpha
        ibd   = 0
        do 60 i = 1, nsub
           k = ind(i)
           dk = d(i)
           if (nbd(k) .ne. 0) then
              if (dk .lt. zero .and. nbd(k) .le. 2) then
                 temp2 = l(k) - x(k)
                 if (temp2 .ge. zero) then
                    temp1 = zero
                 else if (dk*alpha .lt. temp2) then
                    temp1 = temp2/dk
                 endif
              else if (dk .gt. zero .and. nbd(k) .ge. 2) then
                 temp2 = u(k) - x(k)
                 if (temp2 .le. zero) then
                    temp1 = zero
                 else if (dk*alpha .gt. temp2) then
                    temp1 = temp2/dk
                 endif
              endif
              if (temp1 .lt. alpha) then
                 alpha = temp1
                 ibd = i
              endif
           endif
  60    continue

        if (alpha .lt. one) then
           dk = d(ibd)
           k = ind(ibd)
           if (dk .gt. zero) then
              x(k) = u(k)
              d(ibd) = zero
           else if (dk .lt. zero) then
              x(k) = l(k)
              d(ibd) = zero
           endif
        endif
        do 70 i = 1, nsub
           k    = ind(i)
           x(k) = x(k) + alpha*d(i)
  70    continue
     endif
  end if

  if (iprint .ge. 99) write (6,1004)

  1001 format (/,'----------------SUBSM entered-----------------',/)
  1004 format (/,'----------------exit SUBSM --------------------',/)

  return

  end subroutine subsm

  subroutine timer(ttime)
  !> \file timer.f

  !> \brief This routine computes cpu time in double precision.
  !>
  !> This routine computes cpu time in double precision; it makes use of
  !> the intrinsic f90 cpu_time therefore a conversion type is
  !> needed.
  !>
  !> @param ttime CPU time in double precision

  double precision ttime
  !
  real temp
  !
  !           J.L Morales  Departamento de Matematicas,
  !                        Instituto Tecnologico Autonomo de Mexico
  !                        Mexico D.F.
  !
  !           J.L Nocedal  Department of Electrical Engineering and
  !                        Computer Science.
  !                        Northwestern University. Evanston, IL. USA
  !
  !                        January 21, 2011
  !
  temp = sngl(ttime)
  call cpu_time(temp)
  ttime = dble(temp)

  return

  end subroutine timer


end module lbfgsb_all_mod
