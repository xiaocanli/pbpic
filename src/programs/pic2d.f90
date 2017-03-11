!******************************************************************************
! modules for program pic2d
!******************************************************************************
module parameters
    implicit none
  ! numerics
    integer, parameter :: singlep=selected_real_kind(6)
    integer, parameter :: doublep=selected_real_kind(12)
    real(kind=doublep), parameter :: pi=3.14159265358970, small=1.0e-10

  ! simulation domain sizes: # of grid points, grid length
    integer, parameter :: nx=256, ny=256
!    integer, parameter :: nx=2560, ny=2560
!    integer, parameter :: nx=2048, ny=2048
    real, parameter :: dx=0.04, dy=0.04
    real , parameter :: idx = 1.0/dx, idy = 1.0/dy
  
  ! particle information
    real, parameter :: mime=25.0 ! mass ratio of ion and electron
    real, parameter :: cva=15.0  ! ratio light speed and alfven speed
    real, parameter :: cva2= cva*cva
    real, parameter :: vthec=0.092 ! ratio of electron thermal speed and c
    real, parameter :: vdrift=1.0/cva ! drifting speed along z
    real, parameter :: dcs = 0.5   ! half width of the current sheets
    ! ratio of thermal speed of current sheet electron and light speed c
    real, parameter :: tiecs=5.0 ! ion-to-electron temperature ratio for cs
    real, parameter :: tiebg=1.0 ! ion-to-electron temperature ratio for bg
    real, parameter :: tbgcs=0.1 ! bg to cs temperature ratio
    real, parameter :: nbgcs=0.2 ! bg to cs center density ratio
    integer, parameter :: ppg=64 ! # of particles per grid
    real, parameter :: in0=1.0/291.0 ! 1/(n_0*dx*dy)
!    real, parameter :: in0=(nx*dx*0.2+2.0)/(ppg*nx*dx)
!    real, parameter :: in0=1.0/161.897 ! 1/(n_0*dx*dy)
    integer, save :: nptl        ! total # of particles for each species in each
                                 ! computing thread
    real, parameter :: totact=2.0! ratio of total particles to actual particles
    integer, save :: nptltot     ! the size of the particle information array
    integer, save :: nptlact     ! actual number of particles in a sub-domain

  ! time array: # of time steps, interval of steps, modes for diagnostics
    integer, parameter :: ntime=30
    real, parameter :: dt=0.045/mime
    real, parameter :: dth = dt/2.0
    integer :: itime

  ! structures for particle and fields
    type particle
        real dx, dy, ux, uy, uz
        integer grid, q, pid
    end type

    type fields
        real ex, dex_dy, ey, dey_dx
        real ez, dez_dx, dez_dy, dive
        real bx, dbx_dy, by, dby_dx
        real bz, dbz_dx, dbz_dy, pad1
    end type

    ! fields on grid points
    type fields_grid
        real ex, ey, ez, pad1
        real bx, by, bz, pad2
    end type

    type density
        real jx, jy, jz, rho
    end type

  ! particle and fields array, current density and charge density
    type(particle), dimension(:), allocatable, save :: ptl
    type(fields), dimension(:,:), allocatable, save :: emf
    type(fields_grid), dimension(:,:), allocatable, save :: emfgrid
    type(density), dimension(:,:), allocatable, save :: denjr
    !! charge density array for separate treatment of charge density
    !real, dimension(:,:), allocatable, save :: rhoq

  ! particle  number density
    integer, dimension(:,:), allocatable, save :: numdens

  ! particle array and particle allocation for counting sort
    type(particle), dimension(:), allocatable, save :: ptlsort
    integer, dimension(:), allocatable, save :: ptlalloc

  ! buffers for communication of em fields between neighbours
    real, dimension(:), allocatable, save :: bufsendt
    real, dimension(:), allocatable, save :: bufsendb
    real, dimension(:), allocatable, save :: bufrecvt
    real, dimension(:), allocatable, save :: bufrecvb
    integer, save :: bufsize

  ! buffers for communication of current and particle density between neighbours
    real, dimension(:), allocatable, save :: densendt
    real, dimension(:), allocatable, save :: densendb
    real, dimension(:), allocatable, save :: denrecvt
    real, dimension(:), allocatable, save :: denrecvb
    integer, save :: densize

  ! buffers for exchange of particles between neighbours
    type(particle), dimension(:), allocatable :: ptlsendt
    type(particle), dimension(:), allocatable :: ptlsendb
    type(particle), dimension(:), allocatable :: ptlrecvt
    type(particle), dimension(:), allocatable :: ptlrecvb
    integer :: nptlbuf, nptlcrosst, nptlcrossb, nptlcross
    integer :: nptlrecvt, nptlrecvb

  ! array saving the index of crossing particles
    integer, dimension(:), allocatable :: crossindex

  ! id and # of computation processes
    integer, save :: taskid, numtasks, ierr
    integer, save :: mx, my

  ! mpi left(below) and right(above) neighbours
    integer, save :: left, right
  
  ! derived particle data type for mpi
    integer :: particletype, oldtypes(0:1), blockcounts(0:1), &
        offsets(0:1), extent

  ! array saving the magnetic, electric and ion and electron kinetic energies.
  ! the energies are normalized by m_e c^2.
    real, dimension(:,:), allocatable :: totene
    real, dimension(:,:), allocatable :: totene2
  ! array saving particle energy spectra
    integer, dimension(:,:), allocatable :: espectra
    integer, dimension(:,:), allocatable :: espectra_tot
    real, parameter :: minene = 0.0
    real, parameter :: maxenee = 5.0
    real, parameter :: maxenei = 5.0*tiecs
    real, parameter :: dee = 0.005
    real, parameter :: dei = 0.005*tiecs
    integer, parameter :: nbins = (maxenee-minene)/dee
  ! spectral information of fields fluctuations
    integer, parameter :: ndiagy = 2   ! number of diagnostic points along x
    integer, parameter :: ndiagx = 2   ! number of diagnostic points along y
    integer, parameter :: ntdiag = 10  ! read out the data every 1000 steps
    integer, dimension(ndiagx) :: diagindexx
    integer, dimension(:), allocatable :: ndiagxmpi ! allocation of ndiagx in each mpi process
    integer, dimension(ndiagy) :: diagindexy
    type(fields_grid), dimension(:,:,:), allocatable :: fieldsx
    type(fields_grid), dimension(:,:,:), allocatable :: fieldsy
end module parameters

!******************************************************************************
! pic2d main program for 2-d em particle simulation of vlasov plasmas
! xiaocan li 03/24/2014
!******************************************************************************
program pic2d
    use parameters
    implicit none
    include 'mpif.h'
    integer :: itime0

    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, taskid, ierr)
    call mpi_comm_size(mpi_comm_world, numtasks, ierr)

    itime0 = 0
    itime = itime0
    call init
    call countingsort
    call fieldsinfo(1)
!    call particlesinfo(1)
!    call fieldsinfo(0)
!    call particlesinfo(0)
!    call fieldsinfo(1)
!    call particlesinfo(1)
    call energyinfo
!    call dispersioninfo
    itime0 = itime0 + 1
    do itime = itime0, ntime
        print*, itime
        denjr%jx = 0.0
        denjr%jy = 0.0
        denjr%jz = 0.0
        denjr%rho = 0.0
        call pusher
        call fieldsolver
        call dispersioninfo
        if (mod(itime,50) .eq. 0) then
            call countingsort
        endif
        if (mod(itime,ntdiag) .eq. 0) then
            call fieldsinfo(1)
            call energyinfo
            call particlesinfo(1)
        endif
    enddo
    call releasememory
    ! call fftfields
    call mpi_finalize(ierr)
end program pic2d

!******************************************************************************
! initialization of particles and fields
!******************************************************************************
subroutine init
    use parameters
    implicit none
    include 'mpif.h'
    real(kind=doublep) :: lx, xl, xr
    integer :: i, j
    lx = nx*dx
    xl = 0.25*lx
    xr = 0.75*lx
  ! 1d-decomposition along y
    mx = nx
    my = ny/numtasks
    j = mod(ny, numtasks)
    if (taskid .lt. j) then
        my = my + 1
    endif

  ! array saving the magnetic, electric and ion and electron kinetic energies
    allocate(totene(4,0:ntime))
    totene = 0.0

    call initfields

  ! current density and charge density
    allocate(denjr(0:mx+2,0:my+2))
    !allocate(rhoq(0:mx+2,0:my+2))
    denjr%jx = 0.0
    denjr%jy = 0.0
    denjr%jz = 0.0
    denjr%rho = 0.0
    !rhoq = 0.0

  ! buffers for communication for fields between neighbours
    bufsize = (mx+2)*2
    allocate(bufsendt(bufsize))
    allocate(bufsendb(bufsize))
    allocate(bufrecvt(bufsize))
    allocate(bufrecvb(bufsize))
    bufsendt = 0.0
    bufsendb = 0.0
    bufrecvt = 0.0
    bufrecvb = 0.0

  ! buffers for communication for density fields between neighbours
    densize = (mx+2)*5
    allocate(densendt(densize))
    allocate(densendb(densize))
    allocate(denrecvt(densize))
    allocate(denrecvb(densize))
    densendt = 0.0
    densendb = 0.0
    denrecvt = 0.0
    denrecvb = 0.0

    allocate(ptlalloc(mx*my))
    ptlalloc = 0

    call initparticles

    allocate(numdens(0:mx+1,0:my+1))
    numdens = 0
    allocate(ptlsort(nptltot))
    ptlsort = ptl

  ! initialization of buffers for exchange particle information
    nptlbuf = nptlact*0.1
    nptlcrosst = 0 ! # of particles crossing the top
    nptlcrossb = 0 ! # of particles crossing the bottom
    nptlcross = 0  ! total # of particles crossing the domain boundary
    allocate(ptlsendt(nptlbuf))
    allocate(ptlsendb(nptlbuf))
    allocate(ptlrecvt(nptlbuf))
    allocate(ptlrecvb(nptlbuf))
    allocate(crossindex(nptlbuf*2))
    do i = 1, nptlbuf
        ptlsendt(i) = particle(0.0,0.0,0.0,0.0,0.0,1,0,0)
    enddo
    ptlsendb = ptlsendt
    ptlrecvt = ptlsendt
    ptlrecvb = ptlsendt
    crossindex = 0

  ! setup description of the 5 mpi_real data dx, dy, ux, uy, uz
    offsets(0) = 0
    oldtypes(0) = mpi_real
    blockcounts(0) = 5
  ! setup description of the 3 mpi_integer grid, q, test
  ! need to first figure offset by getting size of mpi_real
    call mpi_type_extent(mpi_real, extent, ierr)
    offsets(1) = 5*extent
    oldtypes(1) = mpi_integer
    blockcounts(1) = 3
  ! now define structured type and commit it
    call mpi_type_struct(2, blockcounts, offsets, oldtypes, &
                         particletype, ierr)
    call mpi_type_commit(particletype, ierr)

  ! get the mpi left(below) and right(above) neighbours
    if (taskid .eq. 0) then
        left = numtasks-1
        right = 1
    else if (taskid .eq. numtasks-1) then
        left = taskid-1
        right = 0
    else
        left = taskid-1
        right = taskid+1
    endif

  ! spectral density information for fields fluctuations
    call initspectdens
end subroutine init

!******************************************************************************
! initialization of spectral density information for fields fluctuations
!******************************************************************************
subroutine initspectdens
    use parameters
    implicit none
    include 'mpif.h'
    integer :: i, j
    allocate(ndiagxmpi(-1:numtasks-1))
    ndiagxmpi = 0
  ! the diagnostics indices need to specified where is interested in
    do i = 1, ndiagx
        diagindexx(i) = ny*i/2 - ny/4 - 2 ! index in the whole simulation box
        j = (diagindexx(i)-1)/my
        diagindexx(i) = mod(diagindexx(i)-1,my) + 1
        ndiagxmpi(j) = ndiagxmpi(j) + 1
    enddo
    do i = 0, numtasks-1
        ndiagxmpi(i) = ndiagxmpi(i-1)+ndiagxmpi(i)
    enddo
    j = ndiagxmpi(taskid) - ndiagxmpi(taskid-1)
    allocate(fieldsx(j,mx,ntdiag))
    fieldsx%ex = 0.0
    fieldsx%ey = 0.0
    fieldsx%ez = 0.0
    fieldsx%bx = 0.0
    fieldsx%by = 0.0
    fieldsx%bz = 0.0
    fieldsx%pad1 = 0.0
    fieldsx%pad2 = 0.0
    allocate(fieldsy(ndiagy,my,ntdiag))
    fieldsy%ex = 0.0
    fieldsy%ey = 0.0
    fieldsy%ez = 0.0
    fieldsy%bx = 0.0
    fieldsy%by = 0.0
    fieldsy%bz = 0.0
    fieldsy%pad1 = 0.0
    fieldsy%pad2 = 0.0
    do i = 1, ndiagy
        diagindexy(i) = mx*i/2 - mx/4
    enddo
end subroutine initspectdens
!******************************************************************************
! initialization of particles in a double current sheet and background
! parameters are from oka et al. 2010 apj
!******************************************************************************
subroutine initparticles
    use parameters
    implicit none
    include 'mpif.h'
    real(kind=doublep) :: lx, x, y, xl, xr, ran
    real(kind=singlep) :: ran1, ran2
    real(kind=singlep) :: vthecs, vthics, vthebg, vthibg
    real(kind=singlep) :: uthecs, uthics, uthebg, uthibg, udrift
    real(kind=doublep) :: ux, uy, uz, gama, ms
    integer :: nptlcs, nptlbg ! # of cs and bg particles
    integer :: i, j, i1, i2, ix, iy
    lx = nx*dx
    xl = 0.25*lx
    xr = 0.75*lx
  ! 1d-decomposition along y
    mx = nx
    my = ny/numtasks

    nptl = ppg*mx*my
    nptlcs = ceiling(2.0*nptl/(2.0+0.2*lx)) ! see the document
    nptlcs = (nptlcs/2)*2         ! even number for convenience
    nptlbg = nptl - nptlcs
    vthecs = vthec                ! normalized electron thermal speed for cs
    vthebg = vthecs*sqrt(tbgcs)   ! normalized electron thermal speed for bg
    vthics = vthecs*sqrt(tiecs/mime) ! for ion
    vthibg = vthebg*sqrt(tiebg/mime) ! for ion
    uthecs = vthecs/sqrt(1.0-vthecs**2)
    uthics = vthics/sqrt(1.0-vthics**2)
    uthebg = vthebg/sqrt(1.0-vthebg**2)
    uthibg = vthibg/sqrt(1.0-vthibg**2)
    udrift = vdrift/sqrt(1.0-vdrift**2)
    nptlact = nptl*2              ! actual number of particles
    nptltot = nptl*2*totact       ! the actual array size, making room for
                                  ! particles crossing the sub-domain boundaries
    allocate(ptl(nptltot))
    call init_random_seed()
    ! -------- current sheet particles --------!
    csptl: do i = 1, nptlcs/2
        ! # of ions and electrons are initially same in each cell
        i1 = 2*i-1 ! for ion
        i2 = 2*i   ! for electron
        ix = 0
        do
            call random_number(ran)
          ! avoiding extreme small numbers
            if (ran > 0.5) then
                ran = min(ran,1.0-small)
            else
                ran = max(ran,small)
            endif
            x = dcs*log(ran/(1.0-ran))/2.0
            ix = ceiling((x+xl)/dx)
            if ((ix .gt. 0) .and. (ix .le. mx)) then
                exit
            endif
        enddo
        call random_number(ran)
        y = ran*my*dy
        iy = ceiling(y/dy)
        if (iy .gt. my) then
            iy = my
        end if

        ptl(i1)%dx = (x+xl)/dx - (ix-0.5)
        ptl(i1)%dy = y/dy - (iy-0.5)
        ptl(i1)%grid = my*(ix-1) + iy
        ptl(i1)%q = 1.0
        ptl(i1)%pid = i1  ! for diagnostics test particles
        ptl(i2)%grid = ptl(i1)%grid
        ptl(i2)%dx = ptl(i1)%dx
        call random_number(ran)
        ptl(i2)%dy = ran-0.5
        ptl(i2)%q = -1.0
        ptl(i2)%pid = i2

      ! for counting sort
        ptlalloc(ptl(i1)%grid) = ptlalloc(ptl(i1)%grid) + 1
        ptlalloc(ptl(i2)%grid) = ptlalloc(ptl(i2)%grid) + 1

        call gaussianrand(ran1, ran2)
        ptl(i1)%ux = ran1*uthics
        ptl(i1)%uy = ran2*uthics
        !call gaussianrand(ran1, ran2)
        !ptl(i1)%uz = ran1*uthics + udrift
        ptl(i1)%uz = udrift

        call gaussianrand(ran1, ran2)
        ptl(i2)%ux = ran1*uthecs
        ptl(i2)%uy = ran2*uthecs
        !call gaussianrand(ran1, ran2)
        !ptl(i2)%uz = ran1*uthecs - udrift
        ptl(i2)%uz = -udrift

      ! the other half of the simulation domain
        i1 = i1 + nptlcs
        i2 = i2 + nptlcs
        ix = 0
        do
            call random_number(ran)
            if (ran > 0.5) then
                ran = min(ran,1.0-small)
            else
                ran = max(ran,small)
            endif
            x = dcs*log(ran/(1.0-ran))/2.0
            ix = ceiling((x+xr)/dx)
            if ((ix .gt. 0) .and. (ix .le. mx)) then
                exit
            endif
        enddo
        call random_number(ran)
        y = ran*my*dy
        iy = ceiling(y/dy)
        if (iy .gt. my) then
            iy = my
        end if

        ptl(i1)%dx = (x+xr)/dx - (ix-0.5)
        ptl(i1)%dy = y/dy - (iy-0.5)
        ptl(i1)%grid = my*(ix-1) + iy
        ptl(i1)%q = 1.0
        ptl(i1)%pid = i1  ! for diagnostics test particles
        ptl(i2)%grid = ptl(i1)%grid
        ptl(i2)%dx = ptl(i1)%dx
        call random_number(ran)
        ptl(i2)%dy = ran-0.5
        ptl(i2)%q = -1.0
        ptl(i2)%pid = i2

      ! for counting sort
        ptlalloc(ptl(i1)%grid) = ptlalloc(ptl(i1)%grid) + 1
        ptlalloc(ptl(i2)%grid) = ptlalloc(ptl(i2)%grid) + 1

        call gaussianrand(ran1, ran2)
        ptl(i1)%ux = ran1*uthics
        ptl(i1)%uy = ran2*uthics
        !call gaussianrand(ran1, ran2)
        !ptl(i1)%uz = ran1*uthics + udrift
        ptl(i1)%uz = -udrift

        call gaussianrand(ran1, ran2)
        ptl(i2)%ux = ran1*uthecs
        ptl(i2)%uy = ran2*uthecs
        !call gaussianrand(ran1, ran2)
        !ptl(i2)%uz = ran1*uthecs - udrift
        ptl(i2)%uz = udrift
    end do csptl

    ! -------- background particles --------!
    bgptl: do i = 1, nptlbg
        i1 = nptlcs*2 + 2*i - 1 ! ion
        i2 = nptlcs*2 + 2*i     ! electron
        call random_number(ran)
        x = lx*ran
        ix = ceiling(x/dx)
        if (ix .gt. mx) then
            ix = mx
        end if
        call random_number(ran)
        y = ran*my*dy
        iy = ceiling(y/dy)
        if (iy .gt. my) then
            iy = my
        end if

        ptl(i1)%dx = x/dx - (ix-0.5)
        ptl(i1)%dy = y/dy - (iy-0.5) 
        ptl(i1)%grid = my*(ix-1) + iy
        ptl(i1)%q = 1.0
        ptl(i1)%pid = i1

        ptl(i2)%grid = ptl(i1)%grid
        ptl(i2)%dx = ptl(i1)%dx
        call random_number(ran)
        ptl(i2)%dy = ran-0.5
        ptl(i2)%q = -1.0
        ptl(i2)%pid = i2

      ! for counting sort
        ptlalloc(ptl(i1)%grid) = ptlalloc(ptl(i1)%grid) + 1
        ptlalloc(ptl(i2)%grid) = ptlalloc(ptl(i2)%grid) + 1

        call gaussianrand(ran1, ran2)
        ptl(i1)%ux = ran1*uthibg
        ptl(i1)%uy = ran2*uthibg
        !call gaussianrand(ran1, ran2)
        !ptl(i1)%uz = ran1*uthics + udrift
        ptl(i1)%uz = 0.0

        call gaussianrand(ran1, ran2)
        ptl(i2)%ux = ran1*uthebg
        ptl(i2)%uy = ran2*uthebg
        !call gaussianrand(ran1, ran2)
        !ptl(i2)%uz = ran1*uthecs - udrift
        ptl(i2)%uz = 0.0
    end do bgptl

  ! total kinetic energy for ions and electrons
    do i = 1, nptlact
        ux = ptl(i)%ux
        uy = ptl(i)%uy
        uz = ptl(i)%uz
        gama = sqrt(1.0+ux*ux+uy*uy+uz*uz)
        if (ptl(i)%q .eq. 1) then
            ms = mime
            j = 2
        else
            ms = 1
            j = 1
        endif
        totene(j,itime) = totene(j,itime) + (gama-1.0)*ms
    enddo
  ! empty memory place for more particles
    ptl(nptlact+1:nptltot)%dx = 0.0
    ptl(nptlact+1:nptltot)%dy = 0.0
    ptl(nptlact+1:nptltot)%grid = mx*my + 1
    ptl(nptlact+1:nptltot)%pid = nptltot + 1
    ptl(nptlact+1:nptltot)%ux = 0.0
    ptl(nptlact+1:nptltot)%uy = 0.0
    ptl(nptlact+1:nptltot)%uz = 0.0
    ptl(nptlact+1:nptltot)%q = 0
end subroutine initparticles

!******************************************************************************
! initialization of fields and curl of fields
!******************************************************************************
subroutine initfields
    use parameters
    implicit none
    include 'mpif.h'
    real(kind=doublep) :: lx, x, xl, xr
    integer :: i

    lx = nx*dx
    xl = 0.25*lx
    xr = 0.75*lx
    allocate(emf(0:mx+1,0:my+1))
    emf%ex = 0.0
    emf%ey = 0.0
    emf%ez = 0.0
    emf%dex_dy = 0.0
    emf%dey_dx = 0.0
    emf%dez_dx = 0.0
    emf%dez_dy = 0.0
    emf%bx = 0.0
    emf%bz = 0.0
    emf%dbx_dy = 0.0
    emf%dbz_dx = 0.0
    emf%dbz_dy = 0.0
  ! equation from oka et al. 2010
    do i = 1, mx+1
        x = (i-1)*dx
        emf(i,:)%by = tanh((x-xl)/dcs)-tanh((x-xr)/dcs)-1.0
        ! dby_dx is in the middle of a grid
        x = x+dx/2.0
        emf(i,:)%dby_dx = 2.0/cosh((x-xl)/dcs)**2 - &
                          2.0/cosh((x-xr)/dcs)**2
    enddo
  ! keep the structure as a 4-vector
    emf%dive = 0.0
    emf%pad1 = 0.0

  ! periodic along x
    emf(1,:)%by = emf(mx+1,:)%by
    emf(0,:)%by = emf(mx,:)%by
    emf(0,:)%dby_dx = emf(mx,:)%dby_dx
    emf(mx+1,:)%dby_dx = emf(1,:)%dby_dx

  ! fields on grid points
    allocate(emfgrid(0:mx+1,0:my+1))
    emfgrid%ex = 0.0
    emfgrid%ey = 0.0
    emfgrid%ez = 0.0
    emfgrid%pad1 = 0.0
    emfgrid%bx = 0.0
    emfgrid%by = 0.0
    emfgrid%bz = 0.0
    emfgrid%pad2 = 0.0

    emfgrid(1:mx+1,:)%ex = (emf(1:mx+1,:)%ex+emf(0:mx,:)%ex)*0.5
    emfgrid(:,1:my+1)%ey = (emf(:,1:my+1)%ey+emf(:,0:my)%ey)*0.5
    emfgrid(1:mx+1,1:my+1)%ez = (emf(1:mx+1,1:my+1)%ez+emf(0:mx,1:my+1)%ez+ &
        emf(1:mx+1,0:my)%ez+emf(0:mx,0:my)%ez)*0.25
    emfgrid(1:mx+1,:)%bx = (emf(1:mx+1,:)%bx+emf(0:mx,:)%bx)*0.5
    emfgrid(:,1:my+1)%by = (emf(:,1:my+1)%by+emf(:,0:my)%by)*0.5
    emfgrid(1:mx+1,1:my+1)%bz = (emf(1:mx+1,1:my+1)%bz+emf(0:mx,1:my+1)%bz+ &
        emf(1:mx+1,0:my)%bz+emf(0:mx,0:my)%bz)*0.25

  ! calculate the total fields energy
    call fieldsenergy

end subroutine initfields

!******************************************************************************
! solving maxwell equations
!******************************************************************************
subroutine fieldsenergy
    use parameters
    implicit none
    include 'mpif.h'
    integer :: i, j
    real :: ene
    do j = 1, my
        do i  = 1, mx
          ! electric energy
            ene = emfgrid(i,j)%ex**2 + emfgrid(i,j)%ey**2 + &
                emfgrid(i,j)%ez**2
            totene(3,itime) = totene(3,itime) + ene
          ! magnetic energy
            ene = emfgrid(i,j)%bx**2 + emfgrid(i,j)%by**2 + &
                emfgrid(i,j)%bz**2
            totene(4,itime) = totene(4,itime) + ene
        enddo
    enddo
    totene(3,itime) = 0.5*totene(3,itime)*mime/(cva2*cva2*in0)
    totene(4,itime) = 0.5*totene(4,itime)*mime/(cva2*in0)
end subroutine fieldsenergy

!******************************************************************************
! solving maxwell equations
!******************************************************************************
subroutine fieldsolver
    use parameters
    implicit none
    include 'mpif.h'
    integer :: tag1, tag2
    integer :: stat(mpi_status_size)
    tag1 = 1
    tag2 = 2

  ! half step advance of magnetic field
    emf%bx = emf%bx - dth*emf%dez_dy
    emf%by = emf%by + dth*emf%dez_dx
    emf%bz = emf%bz + dth*(emf%dex_dy-emf%dey_dx)

  ! periodic boundary condition along x
    emf(0,:)%bx = emf(mx,:)%bx
    emf(1,:)%bx = emf(mx+1,:)%bx
    emf(0,:)%by = emf(mx,:)%by
    emf(1,:)%by = emf(mx+1,:)%by
    emf(0,:)%bz = emf(mx,:)%bz
    emf(1,:)%bz = emf(mx+1,:)%bz

  ! exchange the information of by, bz with neighbours
    bufsendt(1:mx+2) = emf(:,my)%by
    bufsendt(mx+3:bufsize) = emf(:,my)%bz
    bufsendb(1:mx+2) = emf(:,1)%by
    bufsendb(mx+3:bufsize) = emf(:,1)%bz
    if (mod(taskid,2) .eq. 0) then
        call mpi_send(bufsendt, bufsize, mpi_real, &
            right, tag1, mpi_comm_world, ierr)
        call mpi_send(bufsendb, bufsize, mpi_real, &
            left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(bufrecvb, bufsize, mpi_real, &
            left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(bufrecvt, bufsize, mpi_real, &
            right, tag2, mpi_comm_world, stat, ierr)
    endif
    if (mod(taskid,2) .eq. 1) then
        call mpi_send(bufsendt, bufsize, mpi_real, &
            right, tag1, mpi_comm_world, ierr)
        call mpi_send(bufsendb, bufsize, mpi_real, &
            left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(bufrecvb, bufsize, mpi_real, &
            left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(bufrecvt, bufsize, mpi_real, &
            right, tag2, mpi_comm_world, stat, ierr)
    endif
    emf(:,my+1)%by = bufrecvt(1:mx+2)
    emf(:,my+1)%bz = bufrecvt(mx+3:bufsize)
    emf(:,0)%by = bufrecvb(1:mx+2)
    emf(:,0)%bz = bufrecvb(mx+3:bufsize)

  ! update terms of curl b
    emf(:,1:my)%dbx_dy = (emf(:,2:my+1)%bx-emf(:,1:my)%bx)*idy
    emf(0:mx,:)%dby_dx = (emf(1:mx+1,:)%by-emf(0:mx,:)%by)*idx
    emf(1:mx+1,:)%dbz_dx = (emf(1:mx+1,:)%bz-emf(0:mx,:)%bz)*idx
    emf(:,1:my+1)%dbz_dy = (emf(:,1:my+1)%bz-emf(:,0:my)%bz)*idy

  ! periodic along x
    emf(mx+1,:)%dby_dx = emf(1,:)%dby_dx

  ! full step advance of electric field
    emf%ex = emf%ex + dt*cva2*(emf%dbz_dy-denjr%jx)
    emf%ey = emf%ey - dt*cva2*(emf%dbz_dx+denjr%jy)
    emf%ez = emf%ez + dt*cva2*(emf%dby_dx-emf%dbx_dy-denjr%jz)

  ! periodic along x
    emf(0,:)%ex = emf(mx,:)%ex
    emf(1,:)%ex = emf(mx+1,:)%ex
    emf(0,:)%ey = emf(mx,:)%ey
    emf(1,:)%ey = emf(mx+1,:)%ey
    emf(0,:)%ez = emf(mx,:)%ez
    emf(1,:)%ez = emf(mx+1,:)%ez

!    call divergenceclean
!    call suppresscerenkov

  ! exchange the information of ey, ez with neighbours
    bufsendt(1:mx+2) = emf(:,my)%ey
    bufsendt(mx+3:bufsize) = emf(:,my)%ez
    bufsendb(1:mx+2) = emf(:,1)%ey
    bufsendb(mx+3:bufsize) = emf(:,1)%ez
    if (mod(taskid,2) .eq. 0) then
        call mpi_send(bufsendt, bufsize, mpi_real, &
            right, tag1, mpi_comm_world, ierr)
        call mpi_send(bufsendb, bufsize, mpi_real, &
            left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(bufrecvb, bufsize, mpi_real, &
            left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(bufrecvt, bufsize, mpi_real, &
            right, tag2, mpi_comm_world, stat, ierr)
    endif
    if (mod(taskid,2) .eq. 1) then
        call mpi_send(bufsendt, bufsize, mpi_real, &
            right, tag1, mpi_comm_world, ierr)
        call mpi_send(bufsendb, bufsize, mpi_real, &
            left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(bufrecvb, bufsize, mpi_real, &
            left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(bufrecvt, bufsize, mpi_real, &
            right, tag2, mpi_comm_world, stat, ierr)
    endif
    emf(:,my+1)%ey = bufrecvt(1:mx+2)
    emf(:,my+1)%ez = bufrecvt(mx+3:bufsize)
    emf(:,0)%ey = bufrecvb(1:mx+2)
    emf(:,0)%ez = bufrecvb(mx+3:bufsize)

  ! update terms of curl and divergence of e
    emf(:,1:my)%dex_dy = (emf(:,2:my+1)%ex-emf(:,1:my)%ex)*idy
    emf(0:mx,:)%dey_dx = (emf(1:mx+1,:)%ey-emf(0:mx,:)%ey)*idx
    emf(1:mx+1,:)%dez_dx = (emf(1:mx+1,:)%ez-emf(0:mx,:)%ez)*idx
    emf(:,1:my+1)%dez_dy = (emf(:,1:my+1)%ez-emf(:,0:my)%ez)*idy
    emf(1:mx+1,1:my+1)%dive = (emf(1:mx+1,1:my+1)%ex-emf(0:mx,1:my+1)%ex)*idx + &
                              (emf(1:mx+1,1:my+1)%ey-emf(1:mx+1,0:mx)%ey)*idy

  ! periodic along x
    emf(mx+1,:)%dey_dx = emf(1,:)%dey_dx

  ! half step advance of magnetic field
    emf%bx = emf%bx - dth*emf%dez_dy
    emf%by = emf%by + dth*emf%dez_dx
    emf%bz = emf%bz + dth*(emf%dex_dy-emf%dey_dx)

  ! exchange the information of by, bz with neighbours
    bufsendt(1:mx+2) = emf(:,my)%by
    bufsendt(mx+3:bufsize) = emf(:,my)%bz
    bufsendb(1:mx+2) = emf(:,1)%by
    bufsendb(mx+3:bufsize) = emf(:,1)%bz
    if (mod(taskid,2) .eq. 0) then
        call mpi_send(bufsendt, bufsize, mpi_real, &
            right, tag1, mpi_comm_world, ierr)
        call mpi_send(bufsendb, bufsize, mpi_real, &
            left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(bufrecvb, bufsize, mpi_real, &
            left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(bufrecvt, bufsize, mpi_real, &
            right, tag2, mpi_comm_world, stat, ierr)
    endif
    if (mod(taskid,2) .eq. 1) then
        call mpi_send(bufsendt, bufsize, mpi_real, &
            right, tag1, mpi_comm_world, ierr)
        call mpi_send(bufsendb, bufsize, mpi_real, &
            left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(bufrecvb, bufsize, mpi_real, &
            left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(bufrecvt, bufsize, mpi_real, &
            right, tag2, mpi_comm_world, stat, ierr)
    endif
    emf(:,my+1)%by = bufrecvt(1:mx+2)
    emf(:,my+1)%bz = bufrecvt(mx+3:bufsize)
    emf(:,0)%by = bufrecvb(1:mx+2)
    emf(:,0)%bz = bufrecvb(mx+3:bufsize)

  ! fields on grid points
    emfgrid(1:mx+1,:)%ex = (emf(1:mx+1,:)%ex+emf(0:mx,:)%ex)*0.5
    emfgrid(:,1:my+1)%ey = (emf(:,1:my+1)%ey+emf(:,0:my)%ey)*0.5
    emfgrid(1:mx+1,1:my+1)%ez = (emf(1:mx+1,1:my+1)%ez+emf(0:mx,1:my+1)%ez+ &
        emf(1:mx+1,0:my)%ez+emf(0:mx,0:my)%ez)*0.25
    emfgrid(1:mx+1,:)%bx = (emf(1:mx+1,:)%bx+emf(0:mx,:)%bx)*0.5
    emfgrid(:,1:my+1)%by = (emf(:,1:my+1)%by+emf(:,0:my)%by)*0.5
    emfgrid(1:mx+1,1:my+1)%bz = (emf(1:mx+1,1:my+1)%bz+emf(0:mx,1:my+1)%bz+ &
        emf(1:mx+1,0:my)%bz+emf(0:mx,0:my)%bz)*0.25

  ! calculate the total fields energy
    call fieldsenergy

!  ! perfect electric conductor (pec) along x
!    emf(1,:)%ey = 0.0
!    emf(mx+1,:)%ey = 0.0
!    emf(0,:)%ez = -emf(1,:)%ez
!    emf(mx+1,:)%ez = -emf(mx,:)%ez
!  ! perfect magnetic conductor (pmc) along x
!    emf(1,:)%by = 0.0
!    emf(mx+1,:)%by = 0.0
!    emf(0,:)%bz = -emf(1,:)%bz
!    emf(mx+1,:)%bz = -emf(mx,:)%bz
end subroutine fieldsolver

!******************************************************************************
! marder passes for divergence clean.
! marder, barry. "a method for incorporating gauss' law into electromagnetic pic
! codes." journal of computational physics 68.1 (1987): 48-55.
!******************************************************************************
subroutine divergenceclean
    use parameters
    implicit none
    include 'mpif.h'
    real, dimension(:,:), allocatable :: marderf
    real, save :: dconst, dtdx, dtdy
    dconst = 0.001
    dtdx = dt*dconst/dx
    dtdy = dt*dconst/dy
    allocate(marderf(0:mx+1,0:my+1))
    marderf = 0.0
    marderf(1:mx+1,1:my+1) = emf(1:mx+1,1:my+1)%dive - &
        cva2*denjr(1:mx+1,1:my+1)%rho
    emf(1:mx,1:my+1)%ex = emf(1:mx,1:my+1)%ex + &
        dtdx*(marderf(2:mx+1,1:my+1)-marderf(1:mx,1:my+1))
    emf(1:mx+1,1:my)%ey = emf(1:mx+1,1:my)%ey + &
        dtdy*(marderf(1:mx+1,2:my+1)-marderf(1:mx+1,1:my))
    deallocate(marderf)
end subroutine divergenceclean

!******************************************************************************
! non-physical cerenkov radiation suppression
!******************************************************************************
subroutine suppresscerenkov
    use parameters
    implicit none
    include 'mpif.h'
    real, save :: tau, dtdx, dtdy
    tau = 0.125*dt
    dtdx = dt*cva2*tau/dx
    dtdy = dt*cva2*tau/dy
    emf(1:mx,1:my+1)%ex = emf(1:mx,1:my+1)%ex - &
        dtdy*(emf(1:mx,1:my+1)%dey_dx-emf(1:mx,0:my)%dey_dx-&
        emf(1:mx,1:my+1)%dex_dy+emf(1:mx,0:my)%dex_dy)
    emf(1:mx+1,1:my)%ey = emf(1:mx+1,1:my)%ey - &
        dtdx*(-emf(1:mx+1,1:my)%dey_dx+emf(0:mx,1:my)%dey_dx+&
        emf(1:mx+1,1:my)%dex_dy-emf(0:mx,1:my)%dex_dy)
    emf(1:mx,1:my)%ez = emf(1:mx,1:my)%ez + &
        dtdx*(emf(2:mx+1,1:my)%dez_dx-emf(1:mx,1:my)%dez_dx) + &
        dtdy*(emf(1:mx,2:my+1)%dez_dy-emf(1:mx,1:my)%dez_dy)
end subroutine suppresscerenkov

!******************************************************************************
! push particles
!******************************************************************************
subroutine pusher
    use parameters
    implicit none
    include 'mpif.h'
    real :: ux, uy, uz, u0x, u0y, u0z, vz
    real :: umx, umy, umz, upx, upy, upz
    real :: tx, ty, tz, t2, sx, sy, sz
    real :: ex, ey, ez, bx, by, bz
    real :: s1, s2, s3, s4
    real :: gama, c0, c1, c2, c3, ms
    real :: dx0, dy0, dx1, dy1
    integer :: i, j, ix, iy, ix1, iy1 
    integer :: dcellx, dcelly ! particle movement in cells
    integer :: tag1, tag2, istart, iend, mx1
    integer :: stat(mpi_status_size)
    tag1 = 1
    tag2 = 2

    mx1 = mx + 2
    nptlcrosst = 0
    nptlcrossb = 0
    nptlcross = 0
    do i = 1, nptlact
        ux = ptl(i)%ux
        uy = ptl(i)%uy
        uz = ptl(i)%uz
        ix = (ptl(i)%grid-1)/my + 1
        iy = mod(ptl(i)%grid-1,my) + 1

!        if ((ix .gt. mx) .or. (iy .gt. my)) then
!            print*, 1, taskid, ix, iy, ptl(i)%grid
!        endif
!        if ((ix .lt. 0) .or. (iy .lt. 0)) then
!            print*, 2, taskid, ix, iy, ptl(i)%grid 
!        endif

        s1 = (0.5+ptl(i)%dx)*(0.5+ptl(i)%dy)
        s2 = 0.5+ptl(i)%dx - s1
        s3 = 0.5+ptl(i)%dy - s1
        s4 = 1.0-s1-s2-s3
        ex = emfgrid(ix,iy)%ex*s4 + emfgrid(ix+1,iy)%ex*s2 + &
             emfgrid(ix,iy+1)%ex*s3 + emfgrid(ix+1,iy+1)%ex*s1
        ey = emfgrid(ix,iy)%ey*s4 + emfgrid(ix+1,iy)%ey*s2 + &
             emfgrid(ix,iy+1)%ey*s3 + emfgrid(ix+1,iy+1)%ey*s1
        ez = emfgrid(ix,iy)%ez*s4 + emfgrid(ix+1,iy)%ez*s2 + &
             emfgrid(ix,iy+1)%ez*s3 + emfgrid(ix+1,iy+1)%ez*s1
        bx = emfgrid(ix,iy)%bx*s4 + emfgrid(ix+1,iy)%bx*s2 + &
             emfgrid(ix,iy+1)%bx*s3 + emfgrid(ix+1,iy+1)%bx*s1
        by = emfgrid(ix,iy)%by*s4 + emfgrid(ix+1,iy)%by*s2 + &
             emfgrid(ix,iy+1)%by*s3 + emfgrid(ix+1,iy+1)%by*s1
        bz = emfgrid(ix,iy)%bz*s4 + emfgrid(ix+1,iy)%bz*s2 + &
             emfgrid(ix,iy+1)%bz*s3 + emfgrid(ix+1,iy+1)%bz*s1
        if(ptl(i)%q .gt. 0) then
            ms = mime
        else
            ms = 1.0
        endif
        c0 = ptl(i)%q*dt*mime/ms/2.0
        c1 = c0/cva
        umx = ux + c1*ex
        umy = uy + c1*ey
        umz = uz + c1*ez
        gama = sqrt(1.0+umx*umx+umy*umy+umz*umz)

        c2 = c0/gama
        tx = c2*bx
        ty = c2*by
        tz = c2*bz
        t2 = 2.0/(1.0+tx*tx+ty*ty+tz*tz)
        sx = t2*tx
        sy = t2*ty
        sz = t2*tz

        u0x = umx + (uy*tz-uz*ty)
        u0y = umy + (uz*tx-ux*tz)
        u0z = umz + (ux*ty-uy*tx)

        upx = umx + (u0y*sz-u0z*sy)
        upy = umy + (u0z*sx-u0x*sz)
        upz = umz + (u0x*sy-u0y*sx)

        ux = upx + c1*ex
        uy = upy + c1*ey
        uz = upz + c1*ez
        gama = sqrt(1.0+ux*ux+uy*uy+uz*uz)

        ptl(i)%ux = ux
        ptl(i)%uy = uy
        ptl(i)%uz = uz

        dx0 = ptl(i)%dx
        dy0 = ptl(i)%dy
        c3 = dt*cva/gama
        dx1 = dx0 + c3*ux*idx
        dy1 = dy0 + c3*uy*idy

        if(abs(dx1) .gt. 0.5) then
            dcellx = nint(sign(1.0,dx1))
            !ptl(i)%dx = sign(abs(dx1)-1.0,-dx1)
            ptl(i)%dx = dx1 + sign(1.0,-dx1)
        else
            dcellx = 0
            ptl(i)%dx = dx1
        endif
        if(abs(dy1) .gt. 0.5) then
            dcelly = nint(sign(1.0,dy1))
            !ptl(i)%dy = sign(abs(dy1)-1.0,-dy1)
            ptl(i)%dy = dy1 + sign(1.0,-dy1)
        else
            dcelly = 0
            ptl(i)%dy = dy1
        endif

        ix1 = ix + dcellx
        iy1 = iy + dcelly

        vz = uz*cva/gama
        call currentdeposit(ix, ix1, iy, iy1, &
            dx0, ptl(i)%dx, dy0, ptl(i)%dy, ptl(i)%q, vz)

        ix = mod(ix1+mx-1,mx) + 1   ! periodic along x
        if (iy1 .gt. my) then
            iy = iy1 - my
            ptl(i)%grid = my*(ix-1) + iy 
            nptlcrosst = nptlcrosst + 1
            nptlcross = nptlcross + 1
            ptlsendt(nptlcrosst) = ptl(i)
            !print*, nptlcrosst, ptlsendt(nptlcrosst)%test
            ptl(i)%pid = nptltot+1 ! flag for crossing particle
            crossindex(nptlcross) = i
        else if( iy1 .lt. 1) then
            iy = iy1 + my
            ptl(i)%grid = my*(ix-1) + iy 
            nptlcrossb = nptlcrossb + 1
            nptlcross = nptlcross + 1
            ptlsendb(nptlcrossb) = ptl(i)
            ptl(i)%pid = nptltot+1
            crossindex(nptlcross) = i
        else
            iy = iy1
            ptl(i)%grid = my*(ix-1) + iy 
        endif

      ! total kinetic energy for ions and electrons
!        if (ptl(i)%q .eq. 1) then
!            j = 2
!        else
!            j = 1
!        endif
        j = ceiling(ptl(i)%q*0.5+1.0)
        totene(j,itime) = totene(j,itime) + (gama-1.0)*ms
    enddo

  ! normalized current density and charge density
    denjr%jx = denjr%jx * in0
    denjr%jy = denjr%jy * in0
    denjr%jz = denjr%jz * in0
    denjr%rho = denjr%rho * in0

  ! reuse the memory of the particles crossing the domain boundaries
    j = nptlact
    do i = 1, nptlcross
        do while (ptl(j)%pid .gt. nptltot)
            j = j - 1
        enddo
        if(j .lt. crossindex(i)) then
            exit
        endif
        ptl(crossindex(i)) = ptl(j)
        j = j - 1
    enddo
    nptlact = nptlact - nptlcross

  ! send and receive the particles numbers first
    if (mod(taskid,2) .eq. 0) then
        call mpi_send(nptlcrosst, 1, mpi_integer, &
            right, tag1, mpi_comm_world, ierr)
        call mpi_send(nptlcrossb, 1, mpi_integer, &
            left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(nptlrecvb, 1, mpi_integer, &
            left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(nptlrecvt, 1, mpi_integer, &
            right, tag2, mpi_comm_world, stat, ierr)
    endif
    if (mod(taskid,2) .eq. 1) then
        call mpi_send(nptlcrosst, 1, mpi_integer, &
            right, tag1, mpi_comm_world, ierr)
        call mpi_send(nptlcrossb, 1, mpi_integer, &
            left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(nptlrecvb, 1, mpi_integer, &
            left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(nptlrecvt, 1, mpi_integer, &
            right, tag2, mpi_comm_world, stat, ierr)
    endif

    !print*, taskid, itime, nptlcrosst, nptlcrossb, nptlrecvt, nptlrecvb

  ! exchange particles with neighbours
    if (mod(taskid,2) .eq. 0) then
        call mpi_send(ptlsendt(1:nptlcrosst), nptlcrosst, &
            particletype, right, tag1, mpi_comm_world, ierr)
        call mpi_send(ptlsendb(1:nptlcrossb), nptlcrossb, &
            particletype, left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(ptlrecvb(1:nptlrecvb), nptlrecvb, &
            particletype, left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(ptlrecvt(1:nptlrecvt), nptlrecvt, &
            particletype, right, tag2, mpi_comm_world, stat, ierr)
    endif

    if (mod(taskid,2) .eq. 1) then
        call mpi_send(ptlsendt(1:nptlcrosst), nptlcrosst, &
            particletype, right, tag1, mpi_comm_world, ierr)
        call mpi_send(ptlsendb(1:nptlcrossb), nptlcrossb, &
            particletype, left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(ptlrecvb(1:nptlrecvb), nptlrecvb, &
            particletype, left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(ptlrecvt(1:nptlrecvt), nptlrecvt, &
            particletype, right, tag2, mpi_comm_world, stat, ierr)
    endif

    !print*, taskid, itime, nptlcrosst, nptlcrossb, nptlrecvt, nptlrecvb

    istart = nptlact + 1
    iend = nptlact + nptlrecvt
    ptl(istart:iend) = ptlrecvt(1:nptlrecvt)
    nptlact = nptlact + nptlrecvt

    istart = nptlact + 1
    iend = nptlact + nptlrecvb
    ptl(istart:iend) = ptlrecvb(1:nptlrecvb)
    nptlact = nptlact + nptlrecvb

  ! electric current density and charge density information exchange
    ! periodic along x
    denjr(0,:)%jx = denjr(0,:)%jx + denjr(mx,:)%jx
    denjr(1,:)%jx = denjr(1,:)%jx + denjr(mx+1,:)%jx
    denjr(mx,:)%jx = denjr(0,:)%jx
    denjr(mx+1,:)%jx = denjr(1,:)%jx
    denjr(0,:)%jy = denjr(0,:)%jy + denjr(mx,:)%jy
    denjr(1,:)%jy = denjr(1,:)%jy + denjr(mx+1,:)%jy
    !denjr(2,:)%jy = denjr(2,:)%jy + denjr(mx+2,:)%jy
    denjr(mx,:)%jy = denjr(0,:)%jy
    denjr(mx+1,:)%jy = denjr(1,:)%jy
    !denjr(mx+2,:)%jy = denjr(2,:)%jy
    denjr(0,:)%jz = denjr(0,:)%jz + denjr(mx,:)%jz
    denjr(1,:)%jz = denjr(1,:)%jz + denjr(mx+1,:)%jz
    !denjr(2,:)%jz = denjr(2,:)%jz + denjr(mx+2,:)%jz
    denjr(mx,:)%jz = denjr(0,:)%jz
    denjr(mx+1,:)%jz = denjr(1,:)%jz
    !denjr(mx+2,:)%jz = denjr(2,:)%jz

    denjr(1,:)%rho = denjr(1,:)%rho + denjr(mx+1,:)%rho
    denjr(mx+1,:)%rho = denjr(1,:)%rho
    denjr(0,:)%rho = denjr(mx,:)%rho

    densendt(1:mx1) = denjr(0:mx+1,my+1)%jx
    densendt(mx1+1:mx1*2) = denjr(0:mx+1,my+2)%jx
    densendt(mx1*2+1:mx1*3) = denjr(0:mx+1,my+1)%jy
    densendt(mx1*3+1:mx1*4) = denjr(0:mx+1,my+1)%jz
    densendt(mx1*4+1:mx1*5) = denjr(0:mx+1,my+1)%rho
    densendb(1:mx1) = denjr(0:mx+1,0)%jx
    densendb(mx1+1:mx1*2) = denjr(0:mx+1,1)%jx
    densendb(mx1*2+1:mx1*3) = denjr(0:mx+1,0)%jy
    densendb(mx1*3+1:mx1*4) = denjr(0:mx+1,0)%jz
    densendb(mx1*4+1:mx1*5) = denjr(0:mx+1,1)%rho

    if (mod(taskid,2) .eq. 0) then
        call mpi_send(densendt, densize, mpi_real, &
            right, tag1, mpi_comm_world, ierr)
        call mpi_send(densendb, densize, mpi_real, &
            left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(denrecvb, densize, mpi_real, &
            left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(denrecvt, densize, mpi_real, &
            right, tag2, mpi_comm_world, stat, ierr)
    endif
    if (mod(taskid,2) .eq. 1) then
        call mpi_send(densendt, densize, mpi_real, &
            right, tag1, mpi_comm_world, ierr)
        call mpi_send(densendb, densize, mpi_real, &
            left, tag2, mpi_comm_world, ierr)
    else
        call mpi_recv(denrecvb, densize, mpi_real, &
            left, tag1, mpi_comm_world, stat, ierr)
        call mpi_recv(denrecvt, densize, mpi_real, &
            right, tag2, mpi_comm_world, stat, ierr)
    endif

    !print*, taskid, itime, nptlcrosst, nptlcrossb, nptlrecvt, nptlrecvb

    denjr(0:mx+1,1)%jx = denjr(0:mx+1,1)%jx + denrecvb(1:mx1)
    denjr(0:mx+1,2)%jx = denjr(0:mx+1,2)%jx + denrecvb(mx1+1:mx1*2)
    denjr(0:mx+1,1)%jy = denjr(0:mx+1,1)%jy + denrecvb(mx1*2+1:mx1*3)
    denjr(0:mx+1,1)%jz = denjr(0:mx+1,1)%jz + denrecvb(mx1*3+1:mx1*4)
    denjr(0:mx+1,1)%rho = denjr(0:mx+1,1)%rho + denrecvb(mx1*4+1:mx1*5) 
    denjr(0:mx+1,my)%jx = denjr(0:mx+1,my)%jx + denrecvt(1:mx1)
    denjr(0:mx+1,my+1)%jx = denjr(0:mx+1,my+1)%jx + denrecvt(mx1+1:mx1*2)
    denjr(0:mx+1,my)%jy = denjr(0:mx+1,my)%jy + denrecvt(mx1*2+1:mx1*3)
    denjr(0:mx+1,my)%jz = denjr(0:mx+1,my)%jz + denrecvt(mx1*3+1:mx1*4)
    denjr(0:mx+1,my+1)%rho = denjr(0:mx+1,my+1)%rho + denrecvt(mx1*4+1:mx1*5) 
end subroutine pusher

!******************************************************************************
! current deposit to get current density after pushing particles
! zigzag scheme from
! umeda, takayuki, et al. "a new charge conservation method in electromagnetic
! particle-in-cell simulations." computer physics communications 156.1 (2003):
! 73-85.
!******************************************************************************
subroutine currentdeposit(i1,i2,j1,j2,dx1,dx2,dy1,dy2,q0, uz)
    use parameters
    implicit none
    include 'mpif.h'
    integer, intent(in) :: i1, i2, j1, j2, q0
    real, intent(in) :: dx1, dx2, dy1, dy2, uz
    real, save :: fx1, fx2, fy1, fy2, wx1, wx2, wy1, wy2
    real, save :: ir, jr, im, jm, dimid, djmid, qx, qy
    real, save :: s1, s2, s3, s4, quz
    integer, save :: di, dj, ix1, iy1
  
    qx = q0*dx/dt
    qy = q0*dy/dt
    im = (i1+i2+dx1+dx2+1)*0.5
    jm = (j1+j2+dy1+dy2+1)*0.5
    ir = min(min(i1,i2)+1.0,max(max(i1,i2)+0.0,im))
    jr = min(min(j1,j2)+1.0,max(max(j1,j2)+0.0,jm))

    fx1 = qx*(ir-i1-dx1-0.5)
    fy1 = qy*(jr-j1-dy1-0.5)
    fx2 = qx*(i2+dx2+0.5-ir)
    fy2 = qy*(j2+dy2+0.5-jr)

    wx1 = (ir-i1+dx1+0.5)*0.5
    wy1 = (jr-j1+dy1+0.5)*0.5
    wx2 = (ir-i2+dx2+0.5)*0.5
    wy2 = (jr-j2+dy2+0.5)*0.5

    denjr(i1,j1)%jx = denjr(i1,j1)%jx + fx1*(1.0-wy1)
    denjr(i1,j1+1)%jx = denjr(i1,j1+1)%jx + fx1*wy1
    denjr(i1,j1)%jy = denjr(i1,j1)%jy + fy1*(1.0-wx1)
    denjr(i1+1,j1)%jy = denjr(i1+1,j1)%jy + fy1*wx1
    denjr(i2,j2)%jx = denjr(i2,j2)%jx + fx2*(1.0-wy2)
    denjr(i2,j2+1)%jx = denjr(i2,j2+1)%jx + fx2*wy2
    denjr(i2,j2)%jy = denjr(i2,j2)%jy + fy2*(1.0-wx2)
    denjr(i2+1,j2)%jy = denjr(i2+1,j2)%jy + fy2*wx2

    dimid = im - i1 - 0.5
    djmid = jm - j1 - 0.5

    di = nint(sign(1.0,dimid))
    dj = nint(sign(1.0,djmid))

    dimid = abs(dimid)
    djmid = abs(djmid)

    s1 = dimid*djmid
    s2 = dimid - s1
    s3 = djmid - s1
    s4 = 1.0-s1-s2-s3

    quz = q0*uz
    denjr(i1,j1)%jz = denjr(i1,j1)%jz + quz*s4
    denjr(i1+di,j1)%jz = denjr(i1+di,j1)%jz + quz*s2
    denjr(i1,j1+dj)%jz = denjr(i1,j1+dj)%jz + quz*s3
    denjr(i1+di,j1+dj)%jz = denjr(i1+di,j1+dj)%jz + quz*s1

  ! charge density accumulation with cic particle shape
  ! the density is accumulated at time n
    s1 = (0.5-dx1)*(0.5-dy1)
    s2 = 0.5 - dx1 - s1
    s3 = 0.5 - dy1 - s1
    s4 = 1.0-s1-s2-s3
    
    ix1 = i1 + 1
    iy1 = j1 + 1
    denjr(i1,j1)%rho = denjr(i1,j1)%rho + s4*q0 
    denjr(ix1,j1)%rho = denjr(ix1,j1)%rho + s2*q0
    denjr(i1,iy1)%rho = denjr(i1,iy1)%rho + s3*q0
    denjr(ix1,iy1)%rho = denjr(ix1,iy1)%rho + s1*q0
end subroutine currentdeposit

!******************************************************************************
! counting sort for particles to accelerate the simulation
! reference: bowers, k. j. "accelerating a particle-in-cell simulation using a
! hybrid counting sort." journal of computational physics 173.2 (2001): 393-411.
!******************************************************************************
subroutine countingsort
    use parameters
    implicit none
    include 'mpif.h'
    integer :: i, j, k
    ptlalloc = 0
  ! count the number of particles in each cell
    do i = 1, nptlact
        j = ptl(i)%grid
        ptlalloc(j) = ptlalloc(j) + 1
    enddo
  ! convert ptlalloc to an allocation
    k = 0
    do i = 1, mx*my
        j = ptlalloc(i)
        ptlalloc(i) = k
        k = k + j
    enddo
  ! sort particle information into ptlsort
    do i = 1, nptlact
        j = ptl(i)%grid
        k = ptlalloc(j) + 1
        ptlalloc(j) = k
        ptlsort(k) = ptl(i)
    enddo
    ptl(1:nptlact) = ptlsort(1:nptlact)
end subroutine countingsort

!******************************************************************************
! diagnostics for particle, fields and energy
!******************************************************************************
subroutine diagnostics
    use parameters
    implicit none
    include 'mpif.h'
end subroutine diagnostics

!******************************************************************************
! particle number density
!******************************************************************************
subroutine numdensity
    use parameters
    implicit none
    integer :: ix, iy, n
    numdens = 0
    do n = 1, nptlact
        ix = (ptl(n)%grid-1)/my + 1
        iy = mod(ptl(n)%grid-1,my) + 1
        numdens(ix,iy) = numdens(ix,iy) + 1
    enddo
    !numdens = numdens * in0
end subroutine numdensity

!******************************************************************************
! diagnostics for fields (e, b, j, rho)
!******************************************************************************
subroutine particlesinfo(isave)
    use parameters
    use hdf5
    implicit none
    include 'mpif.h'
    integer, intent(in) :: isave
    character(len=14) :: fname ! file name
    character(len=2), parameter :: dx_id = "dx"
    character(len=2), parameter :: dy_id = "dy"
    character(len=2), parameter :: ux_id = "ux"
    character(len=2), parameter :: uy_id = "uy"
    character(len=2), parameter :: uz_id = "uz"
    character(len=4), parameter :: grid_id = "grid"
    character(len=1), parameter :: q_id = "q"
    character(len=4), parameter :: pid_id = "pid"
    character(len=7), parameter :: dsetnptl="nptlact"
    character(len=11), parameter :: dsetnpre="nptltot_pre"

    integer(hid_t) :: file_id   ! file identifier
    integer(hid_t) :: dset_dx, dset_dy
    integer(hid_t) :: dset_ux, dset_uy, dset_uz
    integer(hid_t) :: dset_grid, dset_q, dset_pid
    integer(hid_t) :: filespace ! dataspace identifier in file
    integer(hid_t) :: memspace  ! dataspace identifier in memory
    integer(hid_t) :: plist_id  ! property list identifier

    integer(hid_t) :: dset_nptlact
    integer(hid_t) :: dspace_nptlact
    integer(hid_t) :: mem_nptlact

    integer(hsize_t), dimension(1) :: dimsf, dimsfi, dimsna

    integer(hsize_t), dimension(1) :: count
    integer(hsize_t), dimension(1) :: offset
    integer :: rank = 1 ! dataset rank

    integer, dimension(:), allocatable :: nptlarray
    integer :: error, info, comm, i, tag1
    integer :: stat(mpi_status_size)
    comm = mpi_comm_world
!    info = mpi_info_null
    call mpi_info_create(info, ierr)
  ! disable romio's data-sieving
    call mpi_info_set(info, "romio_ds_read", "disable", ierr)
    call mpi_info_set(info, "romio_ds_write", "disable", ierr)
  ! enable romio's collective buffering
    call mpi_info_set(info, "romio_cb_read", "enable", ierr)
    call mpi_info_set(info, "romio_cb_write", "enable", ierr)

    tag1 = 1

    dimsf = (/ppg*nx*ny*2/)
    dimsfi = (/ppg*nx*ny*2/)
    dimsna = (/numtasks/)

    write(fname, "(a12)") "particles.h5"
  ! initialize fortran predefined data types
    call h5open_f(error)

  ! setup file access property list with parallel i/o access
    call h5pcreate_f(h5p_file_access_f, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    call mpi_info_free(info, ierr)

    if (isave .eq. 1) then
      ! create the file collectively
        call h5fcreate_f(fname, h5f_acc_trunc_f, file_id, error, access_prp=plist_id)
    else if (isave .eq. 0) then
      ! read the file collectively
        call h5fopen_f (fname, h5f_acc_rdonly_f, file_id, error, access_prp=plist_id)
    endif
    call h5pclose_f(plist_id, error)

    if (isave .eq. 1) then
      ! create the data space for the dataset
        call h5screate_simple_f(rank, dimsf, filespace, error)
      ! create the dataset with default properties
        call h5dcreate_f(file_id, dx_id, h5t_native_real, filespace, &
                         dset_dx, error)
        call h5dcreate_f(file_id, dy_id, h5t_native_real, filespace, &
                         dset_dy, error)
        call h5dcreate_f(file_id, ux_id, h5t_native_real, filespace, &
                         dset_ux, error)
        call h5dcreate_f(file_id, uy_id, h5t_native_real, filespace, &
                         dset_uy, error)
        call h5dcreate_f(file_id, uz_id, h5t_native_real, filespace, &
                         dset_uz, error)
        call h5dcreate_f(file_id, grid_id, h5t_native_integer, filespace, &
                         dset_grid, error)
        call h5dcreate_f(file_id, q_id, h5t_native_integer, filespace, &
                         dset_q, error)
        call h5dcreate_f(file_id, pid_id, h5t_native_integer, filespace, &
                         dset_pid, error)
    else if (isave .eq. 0) then
        call h5dopen_f(file_id, dx_id, dset_dx, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, dy_id, dset_dy, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, ux_id, dset_ux, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, uy_id, dset_uy, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, uz_id, dset_uz, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, grid_id, dset_grid, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, q_id, dset_q, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, pid_id, dset_pid, error, dapl_id=h5p_default_f)
    endif

    if (isave .eq. 1) then
      ! create the data space for the dataset of actual particle numbers
        call h5screate_simple_f(rank, dimsna, dspace_nptlact, error)
      ! create the dataset
        call h5dcreate_f(file_id, dsetnptl, h5t_native_integer, &
            dspace_nptlact, dset_nptlact, error)
    else if (isave .eq. 0) then
        call h5dopen_f(file_id, dsetnptl, dset_nptlact, error, dapl_id=h5p_default_f)
    endif

  ! create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)

  ! each process defines dataset in memory and writes it to the 
  ! hyperslab in the file
    count(1) = 1
    offset(1) = taskid
    call h5screate_simple_f(rank, count, mem_nptlact, error)

  ! select hyperslab in the file
    call h5dget_space_f(dset_nptlact, dspace_nptlact, error)
    call h5sselect_hyperslab_f(dspace_nptlact, h5s_select_set_f, offset, &
        count, error)
    
    if (isave .eq. 1) then
        call h5dwrite_f(dset_nptlact, h5t_native_integer, nptlact, dimsna, &
            error, file_space_id=dspace_nptlact, mem_space_id=mem_nptlact, &
            xfer_prp=plist_id)
    else if(isave .eq. 0) then
        call h5dread_f(dset_nptlact, h5t_native_integer, nptlact, dimsna, &
            error, file_space_id=dspace_nptlact, mem_space_id=mem_nptlact, &
            xfer_prp=plist_id)
    endif

  ! each process defines dataset in memory and writes it to the hyperslab
  ! in the file

    if (taskid .eq. 0) then
        allocate(nptlarray(0:numtasks-1))
        nptlarray(0) = nptlact
    endif
    
    if (taskid .ne. 0) then
        call mpi_send(nptlact, 1, mpi_integer, 0, tag1, comm, ierr)
    else
        do i = 1, numtasks-1
            call mpi_recv(nptlarray(i), 1, mpi_integer, &
                i, tag1, comm, stat, ierr)
            nptlarray(i) = nptlarray(i) + nptlarray(i-1)
        enddo
    endif

    if (taskid .eq. 0) then
        do i = 1, numtasks-1
            call mpi_send(nptlarray(i-1), 1, mpi_integer, i, tag1, comm, ierr)
        enddo
        offset(1) = 0
    else
        call mpi_recv(offset(1), 1, mpi_integer, 0, tag1, comm, stat, ierr)
    endif
   
    call mpi_barrier(comm, ierr)

    if (taskid .eq. 0) then
        deallocate(nptlarray)
    endif

    count(1) = nptlact
    call h5screate_simple_f(rank, count, memspace, error)

  ! select hyperslab in the file
    call h5dget_space_f(dset_dx, filespace, error)
    call h5dget_space_f(dset_dy, filespace, error)
    call h5dget_space_f(dset_ux, filespace, error)
    call h5dget_space_f(dset_uy, filespace, error)
    call h5dget_space_f(dset_uz, filespace, error)
    call h5dget_space_f(dset_grid, filespace, error)
    call h5dget_space_f(dset_q, filespace, error)
    call h5dget_space_f(dset_pid, filespace, error)
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, offset, &
        count, error)

  ! write the dataset collectively
    if (isave .eq. 1) then
        call h5dwrite_f(dset_dx, h5t_native_real, ptl(1:nptlact)%dx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_dy, h5t_native_real, ptl(1:nptlact)%dy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_ux, h5t_native_real, ptl(1:nptlact)%ux, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_uy, h5t_native_real, ptl(1:nptlact)%uy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_uz, h5t_native_real, ptl(1:nptlact)%uz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_grid, h5t_native_integer, ptl(1:nptlact)%grid, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_q, h5t_native_integer, ptl(1:nptlact)%q, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_pid, h5t_native_integer, ptl(1:nptlact)%pid, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
    else if(isave .eq. 0) then
        call h5dread_f(dset_dx, h5t_native_real, ptl(1:nptlact)%dx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_dy, h5t_native_real, ptl(1:nptlact)%dy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_ux, h5t_native_real, ptl(1:nptlact)%ux, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_uy, h5t_native_real, ptl(1:nptlact)%uy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_uz, h5t_native_real, ptl(1:nptlact)%uz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_grid, h5t_native_integer, ptl(1:nptlact)%grid, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_q, h5t_native_integer, ptl(1:nptlact)%q, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_pid, h5t_native_integer, ptl(1:nptlact)%pid, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
    endif

  ! close dataspaces
    call h5sclose_f(filespace, error)
    call h5sclose_f(dspace_nptlact, error)
    call h5sclose_f(memspace, error)
    call h5sclose_f(mem_nptlact, error)

  ! close the dataset and property list
    call h5dclose_f(dset_dx, error)
    call h5dclose_f(dset_dy, error)
    call h5dclose_f(dset_ux, error)
    call h5dclose_f(dset_uy, error)
    call h5dclose_f(dset_uz, error)
    call h5dclose_f(dset_grid, error)
    call h5dclose_f(dset_q, error)
    call h5dclose_f(dset_pid, error)
    call h5pclose_f(plist_id, error)
    call h5dclose_f(dset_nptlact, error)

  ! close the file
    call h5fclose_f(file_id, error)
  ! close fortran interface
    call h5close_f(error)
end subroutine particlesinfo

!******************************************************************************
! diagnostics for fields (e, b, j, rho)
!******************************************************************************
subroutine fieldsinfo(isave)
    use parameters
    use hdf5
    implicit none
    include 'mpif.h'
    integer, intent(in) :: isave
    character(len=14) :: fname ! file name
    character(len=2), parameter :: ex = "ex"
    character(len=2), parameter :: ey = "ey"
    character(len=2), parameter :: ez = "ez"
    character(len=2), parameter :: bx = "bx"
    character(len=2), parameter :: by = "by"
    character(len=2), parameter :: bz = "bz"
    character(len=2), parameter :: jx = "jx"
    character(len=2), parameter :: jy = "jy"
    character(len=2), parameter :: jz = "jz"
    character(len=3), parameter :: rho = "rho"
    character(len=3), parameter :: num = "num"
    character(len=4), parameter :: dive = "dive"
    integer(hid_t) :: file_id   ! file identifier
    integer(hid_t) :: dset_ex, dset_ey, dset_ez
    integer(hid_t) :: dset_bx, dset_by, dset_bz
    integer(hid_t) :: dset_jx, dset_jy, dset_jz
    integer(hid_t) :: dset_rho, dset_num, dset_dive
    integer(hid_t) :: filespace ! dataspace identifier in file
    integer(hid_t) :: memspace  ! dataspace identifier in memory
    integer(hid_t) :: plist_id  ! property list identifier

    integer(hsize_t), dimension(2) :: dimsf = (/nx+2,ny+1/) ! dataset dimensions.
    integer(hsize_t), dimension(2) :: dimsfi = (/nx+2,ny+1/)

    integer(hsize_t), dimension(2) :: count
    integer(hsize_t), dimension(2) :: offset
    integer :: rank = 2 ! dataset rank

    integer :: error, info, comm, jend, jstart
    integer :: tag1, tag2
    integer :: stat(mpi_status_size)
    comm = mpi_comm_world
!    info = mpi_info_null
    call mpi_info_create(info, ierr)
  ! disable romio's data-sieving
    call mpi_info_set(info, "romio_ds_read", "disable", ierr)
    call mpi_info_set(info, "romio_ds_write", "disable", ierr)
  ! enable romio's collective buffering
    call mpi_info_set(info, "romio_cb_read", "enable", ierr)
    call mpi_info_set(info, "romio_cb_write", "enable", ierr)

    tag1 = 1
    tag2 = 2

    jend = 0
    jstart = 0

    write(fname, "(a6,i5.5,a3)") "fields", itime, ".h5"
  ! initialize fortran predefined data types
    call h5open_f(error)

  ! setup file access property list with parallel i/o access
    call h5pcreate_f(h5p_file_access_f, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    call mpi_info_free(info, ierr)

    if (isave .eq. 1) then
      ! create the file collectively
        call h5fcreate_f(fname, h5f_acc_trunc_f, file_id, error, access_prp=plist_id)
    else if (isave .eq. 0) then
      ! read the file collectively
        call h5fopen_f (fname, h5f_acc_rdonly_f, file_id, error, access_prp=plist_id)
    endif
    call h5pclose_f(plist_id, error)

    if (isave .eq. 1) then
      ! create the data space for the dataset
        call h5screate_simple_f(rank, dimsf, filespace, error)
      ! create the dataset with default properties
        call h5dcreate_f(file_id, ex, h5t_native_real, filespace, &
                         dset_ex, error)
        call h5dcreate_f(file_id, ey, h5t_native_real, filespace, &
                         dset_ey, error)
        call h5dcreate_f(file_id, ez, h5t_native_real, filespace, &
                         dset_ez, error)
        call h5dcreate_f(file_id, bx, h5t_native_real, filespace, &
                         dset_bx, error)
        call h5dcreate_f(file_id, by, h5t_native_real, filespace, &
                         dset_by, error)
        call h5dcreate_f(file_id, bz, h5t_native_real, filespace, &
                         dset_bz, error)
        call h5dcreate_f(file_id, jx, h5t_native_real, filespace, &
                         dset_jx, error)
        call h5dcreate_f(file_id, jy, h5t_native_real, filespace, &
                         dset_jy, error)
        call h5dcreate_f(file_id, jz, h5t_native_real, filespace, &
                         dset_jz, error)
        call h5dcreate_f(file_id, rho, h5t_native_real, filespace, &
                         dset_rho, error)
        call h5dcreate_f(file_id, num, h5t_native_integer, filespace, &
                         dset_num, error)
        call h5dcreate_f(file_id, dive, h5t_native_real, filespace, &
                         dset_dive, error)
        call h5sclose_f(filespace, error)
    else if (isave .eq. 0) then
        call h5dopen_f(file_id, ex, dset_ex, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, ey, dset_ey, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, ez, dset_ez, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, bx, dset_bx, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, by, dset_by, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, bz, dset_bz, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, jx, dset_jx, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, jy, dset_jy, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, jz, dset_jz, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, rho, dset_rho, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, num, dset_num, error, dapl_id=h5p_default_f)
        call h5dopen_f(file_id, dive, dset_dive, error, dapl_id=h5p_default_f)
    endif

    if (isave .eq. 1) then
      ! each process defines dataset in memory and writes it to the hyperslab
      ! in the file
        count(1) = mx+2
        if (taskid .eq. numtasks-1) then
            jend = my+1
        else
            jend = my
        endif
        count(2) = jend
        offset(1) = 0
        offset(2) = taskid * my
    else if (isave .eq. 0) then
        count(1) = mx+2
        offset(1) = 0
        if (taskid .eq. 0) then
            jstart = 1
            offset(2) = taskid*my
        else
            jstart = 0
            offset(2) = taskid*my-1
        endif
        count(2) = my+2-jstart
    endif
    call h5screate_simple_f(rank, count, memspace, error)

  ! select hyperslab in the file
    call h5dget_space_f(dset_ex, filespace, error)
    call h5dget_space_f(dset_ey, filespace, error)
    call h5dget_space_f(dset_ez, filespace, error)
    call h5dget_space_f(dset_bx, filespace, error)
    call h5dget_space_f(dset_by, filespace, error)
    call h5dget_space_f(dset_bz, filespace, error)
    call h5dget_space_f(dset_jx, filespace, error)
    call h5dget_space_f(dset_jy, filespace, error)
    call h5dget_space_f(dset_jz, filespace, error)
    call h5dget_space_f(dset_rho, filespace, error)
    call h5dget_space_f(dset_num, filespace, error)
    call h5dget_space_f(dset_dive, filespace, error)
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, offset, &
        count, error)
  
  ! create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)

  ! write the dataset collectively
    if (isave .eq. 1) then
        call numdensity ! get particle number density

        call h5dwrite_f(dset_ex, h5t_native_real, emf(0:mx+1,1:jend)%ex, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_ey, h5t_native_real, emf(0:mx+1,1:jend)%ey, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_ez, h5t_native_real, emf(0:mx+1,1:jend)%ez, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_bx, h5t_native_real, emf(0:mx+1,1:jend)%bx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_by, h5t_native_real, emf(0:mx+1,1:jend)%by, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_bz, h5t_native_real, emf(0:mx+1,1:jend)%bz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_dive, h5t_native_real, emf(0:mx+1,1:jend)%dive, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_jx, h5t_native_real, denjr(0:mx+1,1:jend)%jx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_jy, h5t_native_real, denjr(0:mx+1,1:jend)%jy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_jz, h5t_native_real, denjr(0:mx+1,1:jend)%jz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_rho, h5t_native_real, denjr(0:mx+1,1:jend)%rho, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dwrite_f(dset_num, h5t_native_integer, numdens(0:mx+1,1:jend), &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
    else if(isave .eq. 0) then
        call h5dread_f(dset_ex, h5t_native_real, emf(0:mx+1,jstart:my+1)%ex, &
            dimsfi, error, mem_space_id=memspace, file_space_id=filespace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_ey, h5t_native_real, emf(0:mx+1,jstart:my+1)%ey, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_ez, h5t_native_real, emf(0:mx+1,jstart:my+1)%ez, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_bx, h5t_native_real, emf(0:mx+1,jstart:my+1)%bx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_by, h5t_native_real, emf(0:mx+1,jstart:my+1)%by, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_bz, h5t_native_real, emf(0:mx+1,jstart:my+1)%bz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_dive, h5t_native_real, emf(0:mx+1,jstart:my+1)%dive, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_jx, h5t_native_real, denjr(0:mx+1,jstart:my+1)%jx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_jy, h5t_native_real, denjr(0:mx+1,jstart:my+1)%jy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_jz, h5t_native_real, denjr(0:mx+1,jstart:my+1)%jz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_rho, h5t_native_real, denjr(0:mx+1,jstart:my+1)%rho, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        call h5dread_f(dset_num, h5t_native_integer, numdens(0:mx+1,jstart:my+1), &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)

      ! exchange the information of by, bz, ey, ez for taskid=0 and numtasks-1
        if (taskid .eq. 0) then
            bufsendb(1:mx+2) = emf(:,1)%by
            bufsendb(mx+3:bufsize) = emf(:,1)%bz
            call mpi_send(bufsendb, bufsize, mpi_real, &
                left, tag2, mpi_comm_world, ierr)
            call mpi_recv(bufrecvb, bufsize, mpi_real, &
                left, tag1, mpi_comm_world, stat, ierr)
            emf(:,0)%by = bufrecvb(1:mx+2)
            emf(:,0)%bz = bufrecvb(mx+3:bufsize)
        else if (taskid .eq. numtasks-1) then
            bufsendt(1:mx+2) = emf(:,my)%by
            bufsendt(mx+3:bufsize) = emf(:,my)%bz
            call mpi_recv(bufrecvt, bufsize, mpi_real, &
                right, tag2, mpi_comm_world, stat, ierr)
            call mpi_send(bufsendt, bufsize, mpi_real, &
                right, tag1, mpi_comm_world, ierr)
            emf(:,my+1)%by = bufrecvt(1:mx+2)
            emf(:,my+1)%bz = bufrecvt(mx+3:bufsize)
        endif
        if (taskid .eq. 0) then
            bufsendb(1:mx+2) = emf(:,1)%ey
            bufsendb(mx+3:bufsize) = emf(:,1)%ez
            call mpi_send(bufsendb, bufsize, mpi_real, &
                left, tag2, mpi_comm_world, ierr)
            call mpi_recv(bufrecvb, bufsize, mpi_real, &
                left, tag1, mpi_comm_world, stat, ierr)
            emf(:,0)%ey = bufrecvb(1:mx+2)
            emf(:,0)%ez = bufrecvb(mx+3:bufsize)
        else if (taskid .eq. numtasks-1) then
            bufsendt(1:mx+2) = emf(:,my)%ey
            bufsendt(mx+3:bufsize) = emf(:,my)%ez
            call mpi_recv(bufrecvt, bufsize, mpi_real, &
                right, tag2, mpi_comm_world, stat, ierr)
            call mpi_send(bufsendt, bufsize, mpi_real, &
                right, tag1, mpi_comm_world, ierr)
            emf(:,my+1)%ey = bufrecvt(1:mx+2)
            emf(:,my+1)%ez = bufrecvt(mx+3:bufsize)
        endif

        emf(:,1:my)%dbx_dy = (emf(:,2:my+1)%bx-emf(:,1:my)%bx)*idy
        emf(0:mx,:)%dby_dx = (emf(1:mx+1,:)%by-emf(0:mx,:)%by)*idx
        emf(1:mx+1,:)%dbz_dx = (emf(1:mx+1,:)%bz-emf(0:mx,:)%bz)*idx
        emf(:,1:my+1)%dbz_dy = (emf(:,1:my+1)%bz-emf(:,0:my)%bz)*idy
        emf(:,1:my)%dex_dy = (emf(:,2:my+1)%ex-emf(:,1:my)%ex)*idy
        emf(0:mx,:)%dey_dx = (emf(1:mx+1,:)%ey-emf(0:mx,:)%ey)*idx
        emf(1:mx+1,:)%dez_dx = (emf(1:mx+1,:)%ez-emf(0:mx,:)%ez)*idx
        emf(:,1:my+1)%dez_dy = (emf(:,1:my+1)%ez-emf(:,0:my)%ez)*idy
        emf(1:mx+1,1:my+1)%dive = (emf(1:mx+1,1:my+1)%ex-emf(0:mx,1:my+1)%ex)*idx + &
                                  (emf(1:mx+1,1:my+1)%ey-emf(1:mx+1,0:mx)%ey)*idy
      ! periodic along x
        emf(mx+1,:)%dby_dx = emf(1,:)%dby_dx
        emf(mx+1,:)%dey_dx = emf(1,:)%dey_dx

        emfgrid(1:mx+1,:)%ex = (emf(1:mx+1,:)%ex+emf(0:mx,:)%ex)*0.5
        emfgrid(:,1:my+1)%ey = (emf(:,1:my+1)%ey+emf(:,0:my)%ey)*0.5
        emfgrid(1:mx+1,1:my+1)%ez = (emf(1:mx+1,1:my+1)%ez+emf(0:mx,1:my+1)%ez+ &
            emf(1:mx+1,0:my)%ez+emf(0:mx,0:my)%ez)*0.25
        emfgrid(1:mx+1,:)%bx = (emf(1:mx+1,:)%bx+emf(0:mx,:)%bx)*0.5
        emfgrid(:,1:my+1)%by = (emf(:,1:my+1)%by+emf(:,0:my)%by)*0.5
        emfgrid(1:mx+1,1:my+1)%bz = (emf(1:mx+1,1:my+1)%bz+emf(0:mx,1:my+1)%bz+ &
            emf(1:mx+1,0:my)%bz+emf(0:mx,0:my)%bz)*0.25
    endif
  ! close dataspaces
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)

  ! close the dataset and property list
    call h5dclose_f(dset_ex, error)
    call h5dclose_f(dset_ey, error)
    call h5dclose_f(dset_ez, error)
    call h5dclose_f(dset_bx, error)
    call h5dclose_f(dset_by, error)
    call h5dclose_f(dset_bz, error)
    call h5dclose_f(dset_jx, error)
    call h5dclose_f(dset_jy, error)
    call h5dclose_f(dset_jz, error)
    call h5dclose_f(dset_rho, error)
    call h5dclose_f(dset_num, error)
    call h5dclose_f(dset_dive, error)
    call h5pclose_f(plist_id, error)

  ! close the file
    call h5fclose_f(file_id, error)
  ! close fortran interface
    call h5close_f(error)
end subroutine fieldsinfo

!******************************************************************************
! accumulate to get the energy spectra for ions and electrons.
!******************************************************************************
subroutine energyspectra
    use parameters
    implicit none
    include 'mpif.h'
    real(kind=doublep) :: ux, uy, uz, ms, gama, de, maxe, ene
    integer :: i, ibin, q, j
    espectra = 0
    do i = 1, nptlact
        ux = ptl(i)%ux
        uy = ptl(i)%uy
        uz = ptl(i)%uz
        q = ptl(i)%q
        gama = sqrt(1.0+ux*ux+uy*uy+uz*uz)
        if (q .eq. 1) then
          ! ion
            ms = mime
            de = dei
            maxe = maxenei
            j = 2
        else
          ! electron
            ms = 1.0
            de = dee
            maxe = maxenee
            j = 1
        endif
        
        ene = (gama-1.0)*ms
        if (ene .lt. maxe) then
            ibin = ceiling((ene-minene)/de)
            espectra(ibin, j) = espectra(ibin, j) + 1
        endif
    enddo
end subroutine energyspectra

!******************************************************************************
! diagnostics for energy spectra, including the total magnetic energy,
! electric energy, particle kinetic energy and particle energy spectra.
!******************************************************************************
subroutine energyinfo
    use parameters
    use hdf5
    implicit none
    include 'mpif.h'
    character(len=15) :: fname1
    character(len=15) :: fname2
    character(len=5) :: gname
    character(len=15) :: dsetname1
    character(len=15) :: dsetname2

    integer(hid_t) :: file_id         ! file identifier
    integer(hid_t) :: group_id         ! group identifier
    integer(hid_t) :: dset_tote        ! total energy
    integer(hid_t) :: dset_spect       ! particle energy spectra
    integer(hid_t) :: dset_ebins       ! energy spectra bins
    integer(hid_t) :: filespace        ! dataspace identifier in file
!    integer(hid_t) :: plist_id         ! property list identifier
    integer :: rank = 2                ! datasets rank
    integer :: i

    integer(hsize_t), dimension(2) :: dims1 = (/2,nbins/)
    integer(hsize_t), dimension(2) :: dims2 = (/4,ntime/)

    real, dimension(:,:), allocatable, save :: enebins

    integer :: error, comm
    comm = mpi_comm_world
!!    info = mpi_info_null
!    call mpi_info_create(info, ierr)
!  ! disable romio's data-sieving
!    call mpi_info_set(info, "romio_ds_read", "disable", ierr)
!    call mpi_info_set(info, "romio_ds_write", "disable", ierr)
!  ! enable romio's collective buffering
!    call mpi_info_set(info, "romio_cb_read", "enable", ierr)
!    call mpi_info_set(info, "romio_cb_write", "enable", ierr)


    write(fname1, "(a10)") "spectra.h5"
    write(fname2, "(a13)") "energy_tot.h5"
    write(gname, "(i5.5)") itime
    dsetname1 = 'spectra'
    dsetname2 = "tote"

    allocate(espectra(2,nbins))
    allocate(espectra_tot(2,nbins))
    allocate(totene2(4,0:ntime))
    if (taskid .eq. 0) then
        allocate(enebins(2,nbins))
        do i = 1, nbins
            enebins(1,i) = i*dee
            enebins(2,i) = i*dei
        enddo
    endif
    call energyspectra
    call mpi_reduce(espectra, espectra_tot, 2*nbins, &
        mpi_integer, mpi_sum, 0, comm, ierr)
    call mpi_reduce(totene, totene2, 4*(ntime+1), &
        mpi_real, mpi_sum, 0, comm, ierr)

    if (taskid .eq. 0) then
      ! initialize fortran interface
        call h5open_f(error)
      ! create or open a file using default properties
        if (itime .eq. 0) then
            call h5fcreate_f(fname1, h5f_acc_trunc_f, file_id, error)
        else
            call h5fopen_f(fname1, h5f_acc_rdwr_f, file_id, error)
        endif

        call h5gcreate_f(file_id, gname, group_id, error)
      ! create the data space for the datasets
        call h5screate_simple_f(rank, dims1, filespace, error)
      ! create the dataset
        call h5dcreate_f(group_id, dsetname1, h5t_native_integer, &
            filespace, dset_spect, error)
      ! write the dataset to the group
        call h5dwrite_f(dset_spect, h5t_native_integer, &
            espectra_tot, dims1, error)
      ! close the dataset
        call h5dclose_f(dset_spect, error)
      ! close the data space for the dataset
        call h5sclose_f(filespace, error)
      ! close the groups
        call h5gclose_f(group_id, error)

        if (itime .eq. 0) then
            call h5screate_simple_f(rank, dims1, filespace, error)
            call h5dcreate_f(file_id, "ebins", h5t_native_real, &
                filespace, dset_ebins, error)
            call h5dwrite_f(dset_ebins, h5t_native_real, &
                enebins, dims1, error)
            call h5dclose_f(dset_ebins, error)
            call h5sclose_f(filespace, error)
        endif
      ! close the files
        call h5fclose_f(file_id, error)

        deallocate(enebins)

        call h5screate_simple_f(rank, dims2, filespace, error)
        if (itime .eq. 0) then
            call h5fcreate_f(fname2, h5f_acc_trunc_f, file_id, error)
            call h5dcreate_f(file_id, dsetname2, h5t_native_real, &
                filespace, dset_tote, error)
        else
            call h5fopen_f(fname2, h5f_acc_rdwr_f, file_id, error)
            call h5dopen_f(file_id, dsetname2, dset_tote, error)
        endif
        call h5dwrite_f(dset_tote, h5t_native_real, totene2, dims2, error)

        call h5dclose_f(dset_tote, error)
        call h5sclose_f(filespace, error)
        call h5fclose_f(file_id, error)
        call h5close_f(error)
    endif
    deallocate(espectra_tot)
    deallocate(espectra)
    deallocate(totene2)
end subroutine energyinfo

!******************************************************************************
! diagnostics for phase space information of fields fluctuation
!******************************************************************************
subroutine dispersioninfo
    use parameters
    implicit none
    include 'mpif.h'

    integer :: comm
    integer :: i, j, k, itdiag
    comm = mpi_comm_world

    itdiag = mod(itime-1, ntdiag) + 1
    do i = 1, ndiagy
        fieldsy(i,1:my,itdiag) = emfgrid(diagindexy(i),1:my)
    enddo
    j = ndiagxmpi(taskid) - ndiagxmpi(taskid-1)
    do i = 1, j
        k = ndiagxmpi(taskid-1)+i
        k = diagindexx(k)
        fieldsx(i,:,itdiag) = emfgrid(1:mx,k)
    enddo
    if (mod(itime, ntdiag) .eq. 0) then
        call fieldsdiag
    endif
end subroutine dispersioninfo

! !******************************************************************************
! ! fft transform of fields components
! !******************************************************************************
! subroutine fftfields
!     use parameters
!     use, intrinsic :: iso_c_binding
!     use hdf5
!     implicit none
!     include 'mpif.h'
!     include 'fftw3-mpi.f03'
!     character(len=15) :: fname
!     character(len=15), dimension(2) :: gnamelist
!     character(len=5), dimension(6) :: dnamelist, dnamelist1
!     integer(hid_t) :: file_id
!     integer(hid_t) :: group_id
!     integer(hid_t) :: dset_id
!     integer(hid_t) :: filespace
!     integer(hid_t) :: memspace
!     integer(hid_t) :: plist_id
! !    integer(hid_t) :: crp_list    ! dataset creation property
!     real, dimension(:,:,:), allocatable :: fields_tot
!     real(kind=doublep), dimension(:,:,:), allocatable :: fieldsfft
! 
!     integer(hsize_t), dimension(3) :: dimsr, maxdimsr
!     integer(hsize_t), dimension(3) :: offset1, count1
!     integer(hsize_t), dimension(3) :: dimsc
! 
!     integer :: rank = 3
!     integer :: error, comm, info, ndiag
!     integer :: i, j, k
! 
!     integer(c_intptr_t) :: l, m
!     type(c_ptr) :: plan, cdata
!     complex(c_double_complex), pointer :: data(:,:)
!     integer(c_intptr_t) :: i1, alloc_local, local_m, local_j_offset
! 
!     comm = mpi_comm_world
! !    info = mpi_info_null
!     call mpi_info_create(info, ierr)
!   ! disable romio's data-sieving
!     call mpi_info_set(info, "romio_ds_read", "disable", ierr)
!     call mpi_info_set(info, "romio_ds_write", "disable", ierr)
!   ! enable romio's collective buffering
!     call mpi_info_set(info, "romio_cb_read", "enable", ierr)
!     call mpi_info_set(info, "romio_cb_write", "enable", ierr)
! 
! 
!     fname = "fields.h5"
!     gnamelist = (/'fieldsalongx', 'fieldsalongy'/)
!     dnamelist = (/'ex', 'ey', 'ez', 'bx', 'by', 'bz'/)
!     dnamelist1 = (/'fftex', 'fftey', 'fftez', 'fftbx', 'fftby', 'fftbz'/)
! 
!     call fftw_mpi_init
!     call h5open_f(error)
!     call h5pcreate_f(h5p_file_access_f, plist_id, error)
!     call h5pset_fapl_mpio_f(plist_id, comm, info, error)
!     call mpi_info_free(info, ierr)
! 
!     call h5fopen_f (fname, h5f_acc_rdwr_f, file_id, &
!         error, access_prp=plist_id)
!     call h5pclose_f(plist_id, error)
! 
!   ! create property list for independent dataset write
!     call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
!     call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)
! 
!     l = ntime
!     dimsc = (/1,1,1/)
!     do i = 1, 2
!         if(i .eq. 1) then
!             m = nx
!             local_m = m/numtasks
!             ndiag = ndiagx
!         else
!             m = ny
!             local_m = my
!             ndiag = ndiagy
!         endif
!         local_j_offset = taskid*local_m
!         alloc_local = fftw_mpi_local_size_2d(m, l, comm, local_m, local_j_offset)
!         cdata = fftw_alloc_complex(alloc_local)
!         call c_f_pointer(cdata, data, [l,local_m])
!       ! create mpi plan for in-place forward dft (note dimension reversal)
!         plan = fftw_mpi_plan_dft_2d(m, l, data, data, comm, fftw_forward, fftw_measure)
!         call h5gopen_f(file_id, gnamelist(i), group_id, error)
!         allocate(fields_tot(ndiag,local_m,l))
!         allocate(fieldsfft(ndiag,local_m,l))
!         do j = 1, 6
!             call h5dopen_f(group_id, dnamelist(j), dset_id, error)
!             call h5dget_space_f(dset_id, filespace, error)
!             call h5sget_simple_extent_dims_f(filespace, dimsr, maxdimsr, error)
! 
!             count1(1) = ndiag
!             count1(2) = local_m
!             count1(3) = l
!             offset1(1) = 0
!             offset1(2) = local_j_offset
!             offset1(3) = 0
! 
!             call h5screate_simple_f(rank, count1, memspace, error)
!             call h5sselect_hyperslab_f(filespace, h5s_select_set_f, offset1, &
!                 count1, error)
!             call h5dread_f(dset_id, h5t_native_real, fields_tot, &
!                 dimsr, error, file_space_id=filespace, mem_space_id=memspace, &
!                 xfer_prp=plist_id)
!             call h5sclose_f(memspace, error)
!             call h5dclose_f(dset_id, error)
!             call h5sclose_f(filespace, error)
! 
!             do k = 1, ndiag
!                 do i1 = 1, local_m
!                     data(:,i1) = dble(fields_tot(k,i1,:))
!                 enddo
!                 call fftw_mpi_execute_dft(plan, data, data)
!                 do i1 = 1, local_m
!                     fieldsfft(k,i1,:) = abs(data(:,i1))
!                 enddo
!             enddo
! 
! !            maxdims = (/h5s_unlimited_f, h5s_unlimited_f, h5s_unlimited_f/)
! !            call h5pcreate_f(h5p_dataset_create_f, crp_list, error)
! !            call h5pset_chunk_f(crp_list, rank, dimsc, error)
! 
!             call h5screate_simple_f(rank, dimsr, filespace, error)
!             call h5dcreate_f(group_id, dnamelist1(j), h5t_native_double, &
!                 filespace, dset_id, error)
!             call h5screate_simple_f(rank, count1, memspace, error)
!             call h5sselect_hyperslab_f(filespace, h5s_select_set_f, offset1, &
!                 count1, error)
!             call h5dwrite_f(dset_id, h5t_native_double, fieldsfft, dimsr, error, &
!                 file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
!             call h5sclose_f(memspace, error)
! 
!             call h5dclose_f(dset_id, error)
!             call h5sclose_f(filespace, error)
! !            call h5pclose_f(crp_list, error)
!         enddo
!         deallocate(fields_tot, fieldsfft)
!         call h5gclose_f(group_id, error)
!         call fftw_destroy_plan(plan)
!         call fftw_free(cdata)
!     enddo
!     call h5pclose_f(plist_id, error)
!     call h5fclose_f(file_id,error)
!     call h5close_f(error)
! end subroutine fftfields
!******************************************************************************
! fields diagnostics for fft transform
!******************************************************************************
subroutine fieldsdiag
    use parameters
    use hdf5
    implicit none
    include 'mpif.h'
    character(len=2), dimension(6) :: dnamelist
    character(len=15) :: fname
    character(len=15) :: groupname
    integer(hid_t) :: file_id
    integer(hid_t) :: group_id
    integer(hid_t) :: dset_ex, dset_ey, dset_ez
    integer(hid_t) :: dset_bx, dset_by, dset_bz
    integer(hid_t) :: filespace
    integer(hid_t) :: memspace
    integer(hid_t) :: plist_id
!    integer(hid_t) :: crp_list    ! dataset creation property

    integer(hsize_t), dimension(3) :: dimsx, dimsy
!    integer(hsize_t), dimension(3) :: maxdims
    integer(hsize_t), dimension(3) :: dimsr, maxdimsr
    integer(hsize_t), dimension(3) :: size1, offset1, count1

    integer :: rank = 3
    integer :: j, error, comm, info
    comm = mpi_comm_world
!    info = mpi_info_null
    call mpi_info_create(info, ierr)
  ! disable romio's data-sieving
    call mpi_info_set(info, "romio_ds_read", "disable", ierr)
    call mpi_info_set(info, "romio_ds_write", "disable", ierr)
  ! enable romio's collective buffering
    call mpi_info_set(info, "romio_cb_read", "enable", ierr)
    call mpi_info_set(info, "romio_cb_write", "enable", ierr)

    dimsx = (/ndiagxmpi(numtasks-1),nx,ntime/)
    dimsy = (/ndiagy,ny,ntime/)

    dnamelist = (/'ex', 'ey', 'ez', 'bx', 'by', 'bz'/)
    write(fname, '(a9)') 'fields.h5'

    call h5open_f(error)
    call h5pcreate_f(h5p_file_access_f, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    call mpi_info_free(info, ierr)

    if (itime .eq. ntdiag) then
        call h5fcreate_f(fname, h5f_acc_trunc_f, file_id, &
            error, access_prp=plist_id)
    else
        call h5fopen_f (fname, h5f_acc_rdwr_f, file_id, &
            error, access_prp=plist_id)
    endif
    call h5pclose_f(plist_id, error)

    groupname = "fieldsalongx"
    if (itime .eq. ntdiag) then
        call h5gcreate_f(file_id, groupname, group_id, error)
        call h5screate_simple_f(rank, dimsx, filespace, error)
        call h5dcreate_f(group_id, dnamelist(1), h5t_native_real, &
            filespace, dset_ex, error)
        call h5dcreate_f(group_id, dnamelist(2), h5t_native_real, &
            filespace, dset_ey, error)
        call h5dcreate_f(group_id, dnamelist(3), h5t_native_real, &
            filespace, dset_ez, error)
        call h5dcreate_f(group_id, dnamelist(4), h5t_native_real, &
            filespace, dset_bx, error)
        call h5dcreate_f(group_id, dnamelist(5), h5t_native_real, &
            filespace, dset_by, error)
        call h5dcreate_f(group_id, dnamelist(6), h5t_native_real, &
            filespace, dset_bz, error)
    else
        call h5gopen_f(file_id, groupname, group_id, error)
        call h5dopen_f(group_id, dnamelist(1), dset_ex, error)
        call h5dopen_f(group_id, dnamelist(2), dset_ey, error)
        call h5dopen_f(group_id, dnamelist(3), dset_ez, error)
        call h5dopen_f(group_id, dnamelist(4), dset_bx, error)
        call h5dopen_f(group_id, dnamelist(5), dset_by, error)
        call h5dopen_f(group_id, dnamelist(6), dset_bz, error)
        call h5dget_space_f(dset_ex, filespace, error)
    endif

  ! create property list for independent dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)
    
    j = ndiagxmpi(taskid) - ndiagxmpi(taskid-1)
    count1(1) = j
    count1(2) = nx
    count1(3) = ntdiag
    offset1(1) = ndiagxmpi(taskid-1)
    offset1(2) = 0
    offset1(3) = itime-ntdiag

    call h5screate_simple_f(rank, count1, memspace, error)
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, offset1, &
        count1, error)

    call h5dwrite_f(dset_ex, h5t_native_real, fieldsx%ex, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    call h5dwrite_f(dset_ey, h5t_native_real, fieldsx%ey, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    call h5dwrite_f(dset_ez, h5t_native_real, fieldsx%ez, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    call h5dwrite_f(dset_bx, h5t_native_real, fieldsx%bx, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    call h5dwrite_f(dset_by, h5t_native_real, fieldsx%by, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    call h5dwrite_f(dset_bz, h5t_native_real, fieldsx%bz, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)

    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_ex, error)
    call h5dclose_f(dset_ey, error)
    call h5dclose_f(dset_ez, error)
    call h5dclose_f(dset_bx, error)
    call h5dclose_f(dset_by, error)
    call h5dclose_f(dset_bz, error)
    call h5sclose_f(filespace, error)
    call h5gclose_f(group_id, error)
    call h5pclose_f(plist_id, error)

  ! create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)
    groupname = "fieldsalongy"

    count1(1) = ndiagy
    count1(2) = my
    count1(3) = ntdiag
    offset1(1) = 0
    offset1(2) = my*taskid
    offset1(3) = itime-ntdiag

    if (itime .eq. ntdiag) then
        call h5gcreate_f(file_id, groupname, group_id, error)
        call h5screate_simple_f(rank, dimsy, filespace, error)
        call h5dcreate_f(group_id, dnamelist(1), h5t_native_real, &
            filespace, dset_ex, error)
        call h5dcreate_f(group_id, dnamelist(2), h5t_native_real, &
            filespace, dset_ey, error)
        call h5dcreate_f(group_id, dnamelist(3), h5t_native_real, &
            filespace, dset_ez, error)
        call h5dcreate_f(group_id, dnamelist(4), h5t_native_real, &
            filespace, dset_bx, error)
        call h5dcreate_f(group_id, dnamelist(5), h5t_native_real, &
            filespace, dset_by, error)
        call h5dcreate_f(group_id, dnamelist(6), h5t_native_real, &
            filespace, dset_bz, error)
    else
        call h5gopen_f(file_id, groupname, group_id, error)
        call h5dopen_f(group_id, dnamelist(1), dset_ex, error)
        call h5dopen_f(group_id, dnamelist(2), dset_ey, error)
        call h5dopen_f(group_id, dnamelist(3), dset_ez, error)
        call h5dopen_f(group_id, dnamelist(4), dset_bx, error)
        call h5dopen_f(group_id, dnamelist(5), dset_by, error)
        call h5dopen_f(group_id, dnamelist(6), dset_bz, error)
        call h5dget_space_f(dset_ex, filespace, error)
    endif
    call h5screate_simple_f(rank, count1, memspace, error)
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, offset1, &
        count1, error)

    call h5dwrite_f(dset_ex, h5t_native_real, fieldsy%ex, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    call h5dwrite_f(dset_ey, h5t_native_real, fieldsy%ey, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    call h5dwrite_f(dset_ez, h5t_native_real, fieldsy%ez, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    call h5dwrite_f(dset_bx, h5t_native_real, fieldsy%bx, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    call h5dwrite_f(dset_by, h5t_native_real, fieldsy%by, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    call h5dwrite_f(dset_bz, h5t_native_real, fieldsy%bz, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)

    call h5dclose_f(dset_ex, error)
    call h5dclose_f(dset_ey, error)
    call h5dclose_f(dset_ez, error)
    call h5dclose_f(dset_bx, error)
    call h5dclose_f(dset_by, error)
    call h5dclose_f(dset_bz, error)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5gclose_f(group_id, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
end subroutine fieldsdiag
!******************************************************************************
! initialize the random seed with a varying seed in order to ensure a different
! random number sequence for each invocation of the program.
! source: http://gcc.gnu.org/onlinedocs/gfortran/random_005fseed.html
!******************************************************************************
subroutine init_random_seed()
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size=n)
    allocate(seed(n))
    ! first try if the os provides a random number generator
    open(unit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
        read(un) seed
        close(un)
    else
        ! fallback to xor:ing the current time and pid. the pid is
        ! useful in case one launches multiple instances of the same
        ! program in parallel.
        call system_clock(count)
        if (count /= 0) then
            t = transfer(count, t)
        else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
            t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279 ! add a prime
        s = ieor(s, pid)
        if (n >= 3) then
            seed(1) = t(1) + 36269
            seed(2) = t(2) + 72551
            seed(3) = pid
            if (n > 3) then
                seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
        else
            seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
    end if
    call random_seed(put=seed)
end subroutine init_random_seed

!******************************************************************************
! generating 2 independent gaussian random numbers using box-muller method
! source: section 7.3.4 of numerical recipes 2007
!******************************************************************************
subroutine gaussianrand(ran1, ran2)
    real, intent(out) :: ran1, ran2
    real :: v1, v2, rsq, ran, fac
    rsq = 2.0
    do while ((rsq .ge. 1.0) .or. (rsq .eq. 0.0))
        call random_number(ran)
        v1 = 2.0*ran-1.0
        call random_number(ran)
        v2 = 2.0*ran-1.0
        rsq = v1*v1 + v2*v2
    end do
    fac = sqrt(-2.0*log(rsq)/rsq)
    ran1 = v1*fac
    ran2 = v2*fac
end subroutine gaussianrand

!******************************************************************************
! free used memory
!******************************************************************************
subroutine releasememory 
    use parameters
    implicit none
    include 'mpif.h'
    integer :: j
    deallocate(emf, emfgrid)
    deallocate(denjr)
    !deallocate(rhoq)
    deallocate(bufsendt, bufsendb)
    deallocate(bufrecvt, bufrecvb)
    deallocate(densendt, densendb)
    deallocate(denrecvt, denrecvb)
    deallocate(ptl)
    deallocate(numdens)
    deallocate(ptlsendt, ptlsendb)
    deallocate(ptlrecvt, ptlrecvb)
    deallocate(ptlsort, ptlalloc)
    deallocate(crossindex)
    deallocate(totene)

    j = ndiagxmpi(taskid) - ndiagxmpi(taskid-1)
    deallocate(fieldsx)
    deallocate(fieldsy)
    deallocate(ndiagxmpi)
    call mpi_type_free(particletype, ierr)
eND SUBROUTINE ReleaseMemory 
