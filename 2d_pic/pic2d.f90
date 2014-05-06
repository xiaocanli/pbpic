!******************************************************************************
! Modules for program PIC2D
!******************************************************************************
MODULE parameters
    IMPLICIT NONE
  ! Numerics
    INTEGER, PARAMETER :: singlep=selected_real_kind(6)
    INTEGER, PARAMETER :: doublep=selected_real_kind(12)
    REAL, PARAMETER :: pi=3.14159265358970, small=1.0e-10

  ! Simulation domain sizes: # of grid points, grid length
    INTEGER, PARAMETER :: nx=256, ny=256
!    INTEGER, PARAMETER :: nx=2560, ny=2560
!    INTEGER, PARAMETER :: nx=2048, ny=2048
    REAL, PARAMETER :: dx=0.04, dy=0.04
    REAL , PARAMETER :: idx = 1.0/dx, idy = 1.0/dy
  
  ! Particle information
    REAL, PARAMETER :: mime=25.0 ! Mass ratio of ion and electron
    REAL, PARAMETER :: cva=15.0  ! Ratio light speed and Alfven speed
    REAL, PARAMETER :: cva2= cva*cva
    REAL, PARAMETER :: vthec=0.092 ! Ratio of electron thermal speed and c
    REAL, PARAMETER :: vdrift=1.0/cva ! Drifting speed along z
    REAL, PARAMETER :: dcs = 0.5   ! Half width of the current sheets
    ! Ratio of thermal speed of current sheet electron and light speed c
    REAL, PARAMETER :: tiecs=5.0 ! Ion-to-electron temperature ratio for CS
    REAL, PARAMETER :: tiebg=1.0 ! Ion-to-electron temperature ratio for BG
    REAL, PARAMETER :: tbgcs=0.1 ! BG to CS temperature ratio
    REAL, PARAMETER :: nbgcs=0.2 ! BG to CS center density ratio
    INTEGER, PARAMETER :: ppg=64 ! # of particles per grid
    REAL, PARAMETER :: in0=1.0/291.0 ! 1/(n_0*dx*dy)
!    REAL, PARAMETER :: in0=(nx*dx*0.2+2.0)/(ppg*nx*dx)
!    REAL, PARAMETER :: in0=1.0/161.897 ! 1/(n_0*dx*dy)
    INTEGER, SAVE :: nptl        ! Total # of particles for each species in each
                                 ! Computing thread
    REAL, PARAMETER :: totact=2.0! Ratio of total particles to actual particles
    INTEGER, SAVE :: nptltot     ! The size of the particle information array
    INTEGER, SAVE :: nptlact     ! Actual number of particles in a sub-domain

  ! Time array: # of time steps, interval of steps, modes for diagnostics
    INTEGER, PARAMETER :: ntime=30
    REAL, PARAMETER :: dt=0.045/mime
    REAL, PARAMETER :: dth = dt/2.0
    INTEGER :: itime

  ! Structures for particle and fields
    TYPE particle
        REAL dx, dy, ux, uy, uz
        INTEGER grid, q, pid
    END TYPE

    TYPE fields
        REAL ex, dex_dy, ey, dey_dx
        REAL ez, dez_dx, dez_dy, divE
        REAL bx, dbx_dy, by, dby_dx
        REAL bz, dbz_dx, dbz_dy, pad1
    END TYPE

    ! Fields on grid points
    TYPE fields_grid
        REAL ex, ey, ez, pad1
        REAL bx, by, bz, pad2
    END TYPE

    TYPE density
        REAL jx, jy, jz, rho
    END TYPE

  ! Particle and fields array, current density and charge density
    TYPE(particle), DIMENSION(:), ALLOCATABLE, SAVE :: ptl
    TYPE(fields), DIMENSION(:,:), ALLOCATABLE, SAVE :: emf
    TYPE(fields_grid), DIMENSION(:,:), ALLOCATABLE, SAVE :: emfGrid
    TYPE(density), DIMENSION(:,:), ALLOCATABLE, SAVE :: denjr
    !! Charge density array for separate treatment of charge density
    !REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: rhoq

  ! Particle  number density
    INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: numdens

  ! Particle array and particle allocation for counting sort
    TYPE(particle), DIMENSION(:), ALLOCATABLE, SAVE :: ptlSort
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ptlAlloc

  ! Buffers for communication of em fields between neighbours
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: bufsendt
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: bufsendb
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: bufrecvt
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: bufrecvb
    INTEGER, SAVE :: bufsize

  ! Buffers for communication of current and particle density between neighbours
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: densendt
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: densendb
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: denrecvt
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: denrecvb
    INTEGER, SAVE :: densize

  ! Buffers for exchange of particles between neighbours
    TYPE(particle), DIMENSION(:), ALLOCATABLE :: ptlsendt
    TYPE(particle), DIMENSION(:), ALLOCATABLE :: ptlsendb
    TYPE(particle), DIMENSION(:), ALLOCATABLE :: ptlrecvt
    TYPE(particle), DIMENSION(:), ALLOCATABLE :: ptlrecvb
    INTEGER :: nptlBuf, nptlCrosst, nptlCrossb, nptlCross
    INTEGER :: nptlRecvt, nptlRecvb

  ! Array saving the index of crossing particles
    INTEGER, DIMENSION(:), ALLOCATABLE :: CrossIndex

  ! Id and # of computation processes
    INTEGER, SAVE :: taskid, numtasks, ierr
    INTEGER, SAVE :: mx, my

  ! MPI left(below) and right(above) neighbours
    INTEGER, SAVE :: left, right
  
  ! Derived particle data type for MPI
    INTEGER :: particletype, oldtypes(0:1), blockcounts(0:1), &
        offsets(0:1), extent

  ! Array saving the magnetic, electric and ion and electron kinetic energies.
  ! The energies are normalized by m_e c^2.
    REAL, DIMENSION(:,:), ALLOCATABLE :: totene
    REAL, DIMENSION(:,:), ALLOCATABLE :: totene2
  ! Array saving particle energy spectra
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: espectra
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: espectra_tot
    REAL, PARAMETER :: minene = 0.0
    REAL, PARAMETER :: maxenee = 5.0
    REAL, PARAMETER :: maxenei = 5.0*tiecs
    REAL, PARAMETER :: dee = 0.005
    REAL, PARAMETER :: dei = 0.005*tiecs
    INTEGER, PARAMETER :: nbins = (maxenee-minene)/dee
  ! Spectral information of fields fluctuations
    INTEGER, PARAMETER :: ndiagy = 2   ! Number of diagnostic points along x
    INTEGER, PARAMETER :: ndiagx = 2   ! Number of diagnostic points along y
    INTEGER, PARAMETER :: ntdiag = 10  ! Read out the data every 1000 steps
    INTEGER, DIMENSION(ndiagx) :: diagIndexx
    INTEGER, DIMENSION(:), ALLOCATABLE :: ndiagxmpi ! Allocation of ndiagx in each mpi process
    INTEGER, DIMENSION(ndiagy) :: diagIndexy
    TYPE(fields_grid), DIMENSION(:,:,:), ALLOCATABLE :: fieldsx
    TYPE(fields_grid), DIMENSION(:,:,:), ALLOCATABLE :: fieldsy
END MODULE parameters

!******************************************************************************
! PIC2D main program for 2-D EM particle simulation of Vlasov plasmas
! Xiaocan Li 03/24/2014
!******************************************************************************
PROGRAM PIC2D
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER :: itime0

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

    itime0 = 0
    itime = itime0
    CALL Init
    CALL CountingSort
    CALL FieldsInfo(1)
!    CALL ParticlesInfo(1)
!    CALL FieldsInfo(0)
!    CALL ParticlesInfo(0)
!    CALL FieldsInfo(1)
!    CALL ParticlesInfo(1)
    CALL EnergyInfo
!    CALL DispersionInfo
    itime0 = itime0 + 1
    DO itime = itime0, ntime
        PRINT*, itime
        denjr%jx = 0.0
        denjr%jy = 0.0
        denjr%jz = 0.0
        denjr%rho = 0.0
        CALL pusher
        CALL FieldSolver
        CALL DispersionInfo
        IF (MOD(itime,50) .EQ. 0) THEN
            CALL CountingSort
        ENDIF
        IF (MOD(itime,ntdiag) .EQ. 0) THEN
            CALL FieldsInfo(1)
            CALL EnergyInfo
            CALL ParticlesInfo(1)
        ENDIF
    ENDDO
    CALL ReleaseMemory
    CALL FFTFields
    CALL MPI_FINALIZE(ierr)
END PROGRAM PIC2D

!******************************************************************************
! Initialization of particles and fields
!******************************************************************************
SUBROUTINE Init
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=doublep) :: lx, xl, xr
    INTEGER :: i, j
    lx = nx*dx
    xl = 0.25*lx
    xr = 0.75*lx
  ! 1D-decomposition along y
    mx = nx
    my = ny/numtasks
    j = MOD(ny, numtasks)
    IF (taskid .LT. j) THEN
        my = my + 1
    ENDIF

  ! Array saving the magnetic, electric and ion and electron kinetic energies
    ALLOCATE(totene(4,0:ntime))
    totene = 0.0

    CALL InitFields

  ! Current density and charge density
    ALLOCATE(denjr(0:mx+2,0:my+2))
    !ALLOCATE(rhoq(0:mx+2,0:my+2))
    denjr%jx = 0.0
    denjr%jy = 0.0
    denjr%jz = 0.0
    denjr%rho = 0.0
    !rhoq = 0.0

  ! Buffers for communication for fields between neighbours
    bufsize = (mx+2)*2
    ALLOCATE(bufsendt(bufsize))
    ALLOCATE(bufsendb(bufsize))
    ALLOCATE(bufrecvt(bufsize))
    ALLOCATE(bufrecvb(bufsize))
    bufsendt = 0.0
    bufsendb = 0.0
    bufrecvt = 0.0
    bufrecvb = 0.0

  ! Buffers for communication for density fields between neighbours
    densize = (mx+2)*5
    ALLOCATE(densendt(densize))
    ALLOCATE(densendb(densize))
    ALLOCATE(denrecvt(densize))
    ALLOCATE(denrecvb(densize))
    densendt = 0.0
    densendb = 0.0
    denrecvt = 0.0
    denrecvb = 0.0

    ALLOCATE(ptlAlloc(mx*my))
    ptlAlloc = 0

    CALL InitParticles

    ALLOCATE(numdens(0:mx+1,0:my+1))
    numdens = 0
    ALLOCATE(ptlSort(nptltot))
    ptlSort = ptl

  ! Initialization of buffers for exchange particle information
    nptlBuf = nptlact*0.1
    nptlCrosst = 0 ! # of particles crossing the top
    nptlCrossb = 0 ! # of particles crossing the bottom
    nptlCross = 0  ! Total # of particles crossing the domain boundary
    ALLOCATE(ptlsendt(nptlBuf))
    ALLOCATE(ptlsendb(nptlBuf))
    ALLOCATE(ptlrecvt(nptlBuf))
    ALLOCATE(ptlrecvb(nptlBuf))
    ALLOCATE(CrossIndex(nptlBuf*2))
    DO i = 1, nptlBuf
        ptlsendt(i) = particle(0.0,0.0,0.0,0.0,0.0,1,0,0)
    ENDDO
    ptlsendb = ptlsendt
    ptlrecvt = ptlsendt
    ptlrecvb = ptlsendt
    CrossIndex = 0

  ! Setup description of the 5 MPI_REAL data dx, dy, ux, uy, uz
    offsets(0) = 0
    oldtypes(0) = MPI_REAL
    blockcounts(0) = 5
  ! Setup description of the 3 MPI_INTEGER grid, q, test
  ! Need to first figure offset by getting size of MPI_REAL
    CALL MPI_TYPE_EXTENT(MPI_REAL, extent, ierr)
    offsets(1) = 5*extent
    oldtypes(1) = MPI_INTEGER
    blockcounts(1) = 3
  ! Now define structured type and commit it
    CALL MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, &
                         particletype, ierr)
    CALL MPI_TYPE_COMMIT(particletype, ierr)

  ! Get the MPI left(below) and right(above) neighbours
    IF (taskid .EQ. 0) THEN
        left = numtasks-1
        right = 1
    ELSE IF (taskid .EQ. numtasks-1) THEN
        left = taskid-1
        right = 0
    ELSE
        left = taskid-1
        right = taskid+1
    ENDIF

  ! Spectral density information for fields fluctuations
    CALL InitSpectDens
END SUBROUTINE Init

!******************************************************************************
! Initialization of spectral density information for fields fluctuations
!******************************************************************************
SUBROUTINE InitSpectDens
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER :: i, j
    ALLOCATE(ndiagxmpi(-1:numtasks-1))
    ndiagxmpi = 0
  ! The diagnostics indices need to specified where is interested in
    DO i = 1, ndiagx
        diagIndexx(i) = ny*i/2 - ny/4 - 2 ! Index in the whole simulation box
        j = (diagIndexx(i)-1)/my
        diagIndexx(i) = MOD(diagIndexx(i)-1,my) + 1
        ndiagxmpi(j) = ndiagxmpi(j) + 1
    ENDDO
    DO i = 0, numtasks-1
        ndiagxmpi(i) = ndiagxmpi(i-1)+ndiagxmpi(i)
    ENDDO
    j = ndiagxmpi(taskid) - ndiagxmpi(taskid-1)
    ALLOCATE(fieldsx(j,mx,ntdiag))
    fieldsx%ex = 0.0
    fieldsx%ey = 0.0
    fieldsx%ez = 0.0
    fieldsx%bx = 0.0
    fieldsx%by = 0.0
    fieldsx%bz = 0.0
    fieldsx%pad1 = 0.0
    fieldsx%pad2 = 0.0
    ALLOCATE(fieldsy(ndiagy,my,ntdiag))
    fieldsy%ex = 0.0
    fieldsy%ey = 0.0
    fieldsy%ez = 0.0
    fieldsy%bx = 0.0
    fieldsy%by = 0.0
    fieldsy%bz = 0.0
    fieldsy%pad1 = 0.0
    fieldsy%pad2 = 0.0
    DO i = 1, ndiagy
        diagIndexy(i) = mx*i/2 - mx/4
    ENDDO
END SUBROUTINE InitSpectDens
!******************************************************************************
! Initialization of particles in a double current sheet and background
! Parameters are from Oka et al. 2010 ApJ
!******************************************************************************
SUBROUTINE InitParticles
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=doublep) :: lx, x, y, xl, xr, ran
    REAL(KIND=singlep) :: ran1, ran2
    REAL(KIND=singlep) :: vthecs, vthics, vthebg, vthibg
    REAL(KIND=singlep) :: uthecs, uthics, uthebg, uthibg, udrift
    REAL(KIND=doublep) :: ux, uy, uz, gama, ms
    INTEGER :: nptlcs, nptlbg ! # of CS and BG particles
    INTEGER :: i, j, i1, i2, ix, iy
    lx = nx*dx
    xl = 0.25*lx
    xr = 0.75*lx
  ! 1D-decomposition along y
    mx = nx
    my = ny/numtasks

    nptl = ppg*mx*my
    nptlcs = CEILING(2.0*nptl/(2.0+0.2*lx)) ! See the document
    nptlcs = (nptlcs/2)*2         ! Even number for convenience
    nptlbg = nptl - nptlcs
    vthecs = vthec                ! Normalized electron thermal speed for CS
    vthebg = vthecs*SQRT(tbgcs)   ! Normalized electron thermal speed for BG
    vthics = vthecs*SQRT(tiecs/mime) ! For ion
    vthibg = vthebg*SQRT(tiebg/mime) ! For ion
    uthecs = vthecs/SQRT(1.0-vthecs**2)
    uthics = vthics/SQRT(1.0-vthics**2)
    uthebg = vthebg/SQRT(1.0-vthebg**2)
    uthibg = vthibg/SQRT(1.0-vthibg**2)
    udrift = vdrift/SQRT(1.0-vdrift**2)
    nptlact = nptl*2              ! Actual number of particles
    nptltot = nptl*2*totact       ! The actual array size, making room for
                                  ! particles crossing the sub-domain boundaries
    ALLOCATE(ptl(nptltot))
    CALL init_random_seed()
    ! -------- current sheet particles --------!
    CSPTL: DO i = 1, nptlcs/2
        ! # of ions and electrons are initially same in each cell
        i1 = 2*i-1 ! For ion
        i2 = 2*i   ! For electron
        ix = 0
        DO
            CALL RANDOM_NUMBER(ran)
          ! Avoiding extreme small numbers
            IF (ran > 0.5) THEN
                ran = MIN(ran,1.0-small)
            ELSE
                ran = MAX(ran,small)
            ENDIF
            x = dcs*LOG(ran/(1.0-ran))/2.0
            ix = CEILING((x+xl)/dx)
            IF ((ix .GT. 0) .AND. (ix .LE. mx)) THEN
                EXIT
            ENDIF
        ENDDO
        CALL RANDOM_NUMBER(ran)
        y = ran*my*dy
        iy = CEILING(y/dy)
        IF (iy .GT. my) THEN
            iy = my
        END IF

        ptl(i1)%dx = (x+xl)/dx - (ix-0.5)
        ptl(i1)%dy = y/dy - (iy-0.5)
        ptl(i1)%grid = my*(ix-1) + iy
        ptl(i1)%q = 1.0
        ptl(i1)%pid = i1  ! For diagnostics test particles
        ptl(i2)%grid = ptl(i1)%grid
        ptl(i2)%dx = ptl(i1)%dx
        CALL RANDOM_NUMBER(ran)
        ptl(i2)%dy = ran-0.5
        ptl(i2)%q = -1.0
        ptl(i2)%pid = i2

      ! For counting sort
        ptlAlloc(ptl(i1)%grid) = ptlAlloc(ptl(i1)%grid) + 1
        ptlAlloc(ptl(i2)%grid) = ptlAlloc(ptl(i2)%grid) + 1

        CALL GaussianRand(ran1, ran2)
        ptl(i1)%ux = ran1*uthics
        ptl(i1)%uy = ran2*uthics
        !CALL GaussianRand(ran1, ran2)
        !ptl(i1)%uz = ran1*uthics + udrift
        ptl(i1)%uz = udrift

        CALL GaussianRand(ran1, ran2)
        ptl(i2)%ux = ran1*uthecs
        ptl(i2)%uy = ran2*uthecs
        !CALL GaussianRand(ran1, ran2)
        !ptl(i2)%uz = ran1*uthecs - udrift
        ptl(i2)%uz = -udrift

      ! The other half of the simulation domain
        i1 = i1 + nptlcs
        i2 = i2 + nptlcs
        ix = 0
        DO
            CALL RANDOM_NUMBER(ran)
            IF (ran > 0.5) THEN
                ran = MIN(ran,1.0-small)
            ELSE
                ran = MAX(ran,small)
            ENDIF
            x = dcs*LOG(ran/(1.0-ran))/2.0
            ix = CEILING((x+xr)/dx)
            IF ((ix .GT. 0) .AND. (ix .LE. mx)) THEN
                EXIT
            ENDIF
        ENDDO
        CALL RANDOM_NUMBER(ran)
        y = ran*my*dy
        iy = CEILING(y/dy)
        IF (iy .GT. my) THEN
            iy = my
        END IF

        ptl(i1)%dx = (x+xr)/dx - (ix-0.5)
        ptl(i1)%dy = y/dy - (iy-0.5)
        ptl(i1)%grid = my*(ix-1) + iy
        ptl(i1)%q = 1.0
        ptl(i1)%pid = i1  ! For diagnostics test particles
        ptl(i2)%grid = ptl(i1)%grid
        ptl(i2)%dx = ptl(i1)%dx
        CALL RANDOM_NUMBER(ran)
        ptl(i2)%dy = ran-0.5
        ptl(i2)%q = -1.0
        ptl(i2)%pid = i2

      ! For counting sort
        ptlAlloc(ptl(i1)%grid) = ptlAlloc(ptl(i1)%grid) + 1
        ptlAlloc(ptl(i2)%grid) = ptlAlloc(ptl(i2)%grid) + 1

        CALL GaussianRand(ran1, ran2)
        ptl(i1)%ux = ran1*uthics
        ptl(i1)%uy = ran2*uthics
        !CALL GaussianRand(ran1, ran2)
        !ptl(i1)%uz = ran1*uthics + udrift
        ptl(i1)%uz = -udrift

        CALL GaussianRand(ran1, ran2)
        ptl(i2)%ux = ran1*uthecs
        ptl(i2)%uy = ran2*uthecs
        !CALL GaussianRand(ran1, ran2)
        !ptl(i2)%uz = ran1*uthecs - udrift
        ptl(i2)%uz = udrift
    END DO CSPTL

    ! -------- background particles --------!
    BGPTL: DO i = 1, nptlbg
        i1 = nptlcs*2 + 2*i - 1 ! Ion
        i2 = nptlcs*2 + 2*i     ! Electron
        CALL RANDOM_NUMBER(ran)
        x = lx*ran
        ix = CEILING(x/dx)
        IF (ix .GT. mx) THEN
            ix = mx
        END IF
        CALL RANDOM_NUMBER(ran)
        y = ran*my*dy
        iy = CEILING(y/dy)
        IF (iy .GT. my) THEN
            iy = my
        END IF

        ptl(i1)%dx = x/dx - (ix-0.5)
        ptl(i1)%dy = y/dy - (iy-0.5) 
        ptl(i1)%grid = my*(ix-1) + iy
        ptl(i1)%q = 1.0
        ptl(i1)%pid = i1

        ptl(i2)%grid = ptl(i1)%grid
        ptl(i2)%dx = ptl(i1)%dx
        CALL RANDOM_NUMBER(ran)
        ptl(i2)%dy = ran-0.5
        ptl(i2)%q = -1.0
        ptl(i2)%pid = i2

      ! For counting sort
        ptlAlloc(ptl(i1)%grid) = ptlAlloc(ptl(i1)%grid) + 1
        ptlAlloc(ptl(i2)%grid) = ptlAlloc(ptl(i2)%grid) + 1

        CALL GaussianRand(ran1, ran2)
        ptl(i1)%ux = ran1*uthibg
        ptl(i1)%uy = ran2*uthibg
        !CALL GaussianRand(ran1, ran2)
        !ptl(i1)%uz = ran1*uthics + udrift
        ptl(i1)%uz = 0.0

        CALL GaussianRand(ran1, ran2)
        ptl(i2)%ux = ran1*uthebg
        ptl(i2)%uy = ran2*uthebg
        !CALL GaussianRand(ran1, ran2)
        !ptl(i2)%uz = ran1*uthecs - udrift
        ptl(i2)%uz = 0.0
    END DO BGPTL

  ! Total kinetic energy for ions and electrons
    DO i = 1, nptlact
        ux = ptl(i)%ux
        uy = ptl(i)%uy
        uz = ptl(i)%uz
        gama = SQRT(1.0+ux*ux+uy*uy+uz*uz)
        IF (ptl(i)%q .EQ. 1) THEN
            ms = mime
            j = 2
        ELSE
            ms = 1
            j = 1
        ENDIF
        totene(j,itime) = totene(j,itime) + (gama-1.0)*ms
    ENDDO
  ! Empty memory place for more particles
    ptl(nptlact+1:nptltot)%dx = 0.0
    ptl(nptlact+1:nptltot)%dy = 0.0
    ptl(nptlact+1:nptltot)%grid = mx*my + 1
    ptl(nptlact+1:nptltot)%pid = nptltot + 1
    ptl(nptlact+1:nptltot)%ux = 0.0
    ptl(nptlact+1:nptltot)%uy = 0.0
    ptl(nptlact+1:nptltot)%uz = 0.0
    ptl(nptlact+1:nptltot)%q = 0
END SUBROUTINE InitParticles

!******************************************************************************
! Initialization of fields and curl of fields
!******************************************************************************
SUBROUTINE InitFields
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=doublep) :: lx, x, xl, xr
    INTEGER :: i

    lx = nx*dx
    xl = 0.25*lx
    xr = 0.75*lx
    ALLOCATE(emf(0:mx+1,0:my+1))
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
  ! Equation from Oka et al. 2010
    DO i = 1, mx+1
        x = (i-1)*dx
        emf(i,:)%by = TANH((x-xl)/dcs)-TANH((x-xr)/dcs)-1.0
        ! dby_dx is in the middle of a grid
        x = x+dx/2.0
        emf(i,:)%dby_dx = 2.0/COSH((x-xl)/dcs)**2 - &
                          2.0/COSH((x-xr)/dcs)**2
    ENDDO
  ! Keep the structure as a 4-vector
    emf%divE = 0.0
    emf%pad1 = 0.0

  ! Periodic along x
    emf(1,:)%by = emf(mx+1,:)%by
    emf(0,:)%by = emf(mx,:)%by
    emf(0,:)%dby_dx = emf(mx,:)%dby_dx
    emf(mx+1,:)%dby_dx = emf(1,:)%dby_dx

  ! Fields on grid points
    ALLOCATE(emfGrid(0:mx+1,0:my+1))
    emfGrid%ex = 0.0
    emfGrid%ey = 0.0
    emfGrid%ez = 0.0
    emfGrid%pad1 = 0.0
    emfGrid%bx = 0.0
    emfGrid%by = 0.0
    emfGrid%bz = 0.0
    emfGrid%pad2 = 0.0

    emfGrid(1:mx+1,:)%ex = (emf(1:mx+1,:)%ex+emf(0:mx,:)%ex)*0.5
    emfGrid(:,1:my+1)%ey = (emf(:,1:my+1)%ey+emf(:,0:my)%ey)*0.5
    emfGrid(1:mx+1,1:my+1)%ez = (emf(1:mx+1,1:my+1)%ez+emf(0:mx,1:my+1)%ez+ &
        emf(1:mx+1,0:my)%ez+emf(0:mx,0:my)%ez)*0.25
    emfGrid(1:mx+1,:)%bx = (emf(1:mx+1,:)%bx+emf(0:mx,:)%bx)*0.5
    emfGrid(:,1:my+1)%by = (emf(:,1:my+1)%by+emf(:,0:my)%by)*0.5
    emfGrid(1:mx+1,1:my+1)%bz = (emf(1:mx+1,1:my+1)%bz+emf(0:mx,1:my+1)%bz+ &
        emf(1:mx+1,0:my)%bz+emf(0:mx,0:my)%bz)*0.25

  ! Calculate the total fields energy
    CALL FieldsEnergy

END SUBROUTINE InitFields

!******************************************************************************
! Solving Maxwell equations
!******************************************************************************
SUBROUTINE FieldsEnergy
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER :: i, j
    REAL :: ene
    DO j = 1, my
        DO i  = 1, mx
          ! Electric energy
            ene = emfGrid(i,j)%ex**2 + emfGrid(i,j)%ey**2 + &
                emfGrid(i,j)%ez**2
            totene(3,itime) = totene(3,itime) + ene
          ! Magnetic energy
            ene = emfGrid(i,j)%bx**2 + emfGrid(i,j)%by**2 + &
                emfGrid(i,j)%bz**2
            totene(4,itime) = totene(4,itime) + ene
        ENDDO
    ENDDO
    totene(3,itime) = 0.5*totene(3,itime)*mime/(cva2*cva2*in0)
    totene(4,itime) = 0.5*totene(4,itime)*mime/(cva2*in0)
END SUBROUTINE FieldsEnergy

!******************************************************************************
! Solving Maxwell equations
!******************************************************************************
SUBROUTINE FieldSolver
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER :: tag1, tag2
    INTEGER :: stat(MPI_STATUS_SIZE)
    tag1 = 1
    tag2 = 2

  ! Half step advance of magnetic field
    emf%bx = emf%bx - dth*emf%dez_dy
    emf%by = emf%by + dth*emf%dez_dx
    emf%bz = emf%bz + dth*(emf%dex_dy-emf%dey_dx)

  ! Periodic boundary condition along x
    emf(0,:)%bx = emf(mx,:)%bx
    emf(1,:)%bx = emf(mx+1,:)%bx
    emf(0,:)%by = emf(mx,:)%by
    emf(1,:)%by = emf(mx+1,:)%by
    emf(0,:)%bz = emf(mx,:)%bz
    emf(1,:)%bz = emf(mx+1,:)%bz

  ! Exchange the information of By, Bz with neighbours
    bufsendt(1:mx+2) = emf(:,my)%by
    bufsendt(mx+3:bufsize) = emf(:,my)%bz
    bufsendb(1:mx+2) = emf(:,1)%by
    bufsendb(mx+3:bufsize) = emf(:,1)%bz
    IF (MOD(taskid,2) .EQ. 0) THEN
        CALL MPI_SEND(bufsendt, bufsize, MPI_REAL, &
            right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(bufsendb, bufsize, MPI_REAL, &
            left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(bufrecvb, bufsize, MPI_REAL, &
            left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(bufrecvt, bufsize, MPI_REAL, &
            right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF
    IF (MOD(taskid,2) .EQ. 1) THEN
        CALL MPI_SEND(bufsendt, bufsize, MPI_REAL, &
            right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(bufsendb, bufsize, MPI_REAL, &
            left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(bufrecvb, bufsize, MPI_REAL, &
            left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(bufrecvt, bufsize, MPI_REAL, &
            right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF
    emf(:,my+1)%by = bufrecvt(1:mx+2)
    emf(:,my+1)%bz = bufrecvt(mx+3:bufsize)
    emf(:,0)%by = bufrecvb(1:mx+2)
    emf(:,0)%bz = bufrecvb(mx+3:bufsize)

  ! Update terms of curl B
    emf(:,1:my)%dbx_dy = (emf(:,2:my+1)%bx-emf(:,1:my)%bx)*idy
    emf(0:mx,:)%dby_dx = (emf(1:mx+1,:)%by-emf(0:mx,:)%by)*idx
    emf(1:mx+1,:)%dbz_dx = (emf(1:mx+1,:)%bz-emf(0:mx,:)%bz)*idx
    emf(:,1:my+1)%dbz_dy = (emf(:,1:my+1)%bz-emf(:,0:my)%bz)*idy

  ! Periodic along x
    emf(mx+1,:)%dby_dx = emf(1,:)%dby_dx

  ! Full step advance of electric field
    emf%ex = emf%ex + dt*cva2*(emf%dbz_dy-denjr%jx)
    emf%ey = emf%ey - dt*cva2*(emf%dbz_dx+denjr%jy)
    emf%ez = emf%ez + dt*cva2*(emf%dby_dx-emf%dbx_dy-denjr%jz)

  ! Periodic along x
    emf(0,:)%ex = emf(mx,:)%ex
    emf(1,:)%ex = emf(mx+1,:)%ex
    emf(0,:)%ey = emf(mx,:)%ey
    emf(1,:)%ey = emf(mx+1,:)%ey
    emf(0,:)%ez = emf(mx,:)%ez
    emf(1,:)%ez = emf(mx+1,:)%ez

!    CALL DivergenceClean
!    CALL SuppressCerenkov

  ! Exchange the information of Ey, Ez with neighbours
    bufsendt(1:mx+2) = emf(:,my)%ey
    bufsendt(mx+3:bufsize) = emf(:,my)%ez
    bufsendb(1:mx+2) = emf(:,1)%ey
    bufsendb(mx+3:bufsize) = emf(:,1)%ez
    IF (MOD(taskid,2) .EQ. 0) THEN
        CALL MPI_SEND(bufsendt, bufsize, MPI_REAL, &
            right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(bufsendb, bufsize, MPI_REAL, &
            left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(bufrecvb, bufsize, MPI_REAL, &
            left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(bufrecvt, bufsize, MPI_REAL, &
            right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF
    IF (MOD(taskid,2) .EQ. 1) THEN
        CALL MPI_SEND(bufsendt, bufsize, MPI_REAL, &
            right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(bufsendb, bufsize, MPI_REAL, &
            left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(bufrecvb, bufsize, MPI_REAL, &
            left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(bufrecvt, bufsize, MPI_REAL, &
            right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF
    emf(:,my+1)%ey = bufrecvt(1:mx+2)
    emf(:,my+1)%ez = bufrecvt(mx+3:bufsize)
    emf(:,0)%ey = bufrecvb(1:mx+2)
    emf(:,0)%ez = bufrecvb(mx+3:bufsize)

  ! Update terms of curl and divergence of E
    emf(:,1:my)%dex_dy = (emf(:,2:my+1)%ex-emf(:,1:my)%ex)*idy
    emf(0:mx,:)%dey_dx = (emf(1:mx+1,:)%ey-emf(0:mx,:)%ey)*idx
    emf(1:mx+1,:)%dez_dx = (emf(1:mx+1,:)%ez-emf(0:mx,:)%ez)*idx
    emf(:,1:my+1)%dez_dy = (emf(:,1:my+1)%ez-emf(:,0:my)%ez)*idy
    emf(1:mx+1,1:my+1)%divE = (emf(1:mx+1,1:my+1)%ex-emf(0:mx,1:my+1)%ex)*idx + &
                              (emf(1:mx+1,1:my+1)%ey-emf(1:mx+1,0:mx)%ey)*idy

  ! Periodic along x
    emf(mx+1,:)%dey_dx = emf(1,:)%dey_dx

  ! Half step advance of magnetic field
    emf%bx = emf%bx - dth*emf%dez_dy
    emf%by = emf%by + dth*emf%dez_dx
    emf%bz = emf%bz + dth*(emf%dex_dy-emf%dey_dx)

  ! Exchange the information of By, Bz with neighbours
    bufsendt(1:mx+2) = emf(:,my)%by
    bufsendt(mx+3:bufsize) = emf(:,my)%bz
    bufsendb(1:mx+2) = emf(:,1)%by
    bufsendb(mx+3:bufsize) = emf(:,1)%bz
    IF (MOD(taskid,2) .EQ. 0) THEN
        CALL MPI_SEND(bufsendt, bufsize, MPI_REAL, &
            right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(bufsendb, bufsize, MPI_REAL, &
            left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(bufrecvb, bufsize, MPI_REAL, &
            left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(bufrecvt, bufsize, MPI_REAL, &
            right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF
    IF (MOD(taskid,2) .EQ. 1) THEN
        CALL MPI_SEND(bufsendt, bufsize, MPI_REAL, &
            right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(bufsendb, bufsize, MPI_REAL, &
            left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(bufrecvb, bufsize, MPI_REAL, &
            left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(bufrecvt, bufsize, MPI_REAL, &
            right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF
    emf(:,my+1)%by = bufrecvt(1:mx+2)
    emf(:,my+1)%bz = bufrecvt(mx+3:bufsize)
    emf(:,0)%by = bufrecvb(1:mx+2)
    emf(:,0)%bz = bufrecvb(mx+3:bufsize)

  ! Fields on grid points
    emfGrid(1:mx+1,:)%ex = (emf(1:mx+1,:)%ex+emf(0:mx,:)%ex)*0.5
    emfGrid(:,1:my+1)%ey = (emf(:,1:my+1)%ey+emf(:,0:my)%ey)*0.5
    emfGrid(1:mx+1,1:my+1)%ez = (emf(1:mx+1,1:my+1)%ez+emf(0:mx,1:my+1)%ez+ &
        emf(1:mx+1,0:my)%ez+emf(0:mx,0:my)%ez)*0.25
    emfGrid(1:mx+1,:)%bx = (emf(1:mx+1,:)%bx+emf(0:mx,:)%bx)*0.5
    emfGrid(:,1:my+1)%by = (emf(:,1:my+1)%by+emf(:,0:my)%by)*0.5
    emfGrid(1:mx+1,1:my+1)%bz = (emf(1:mx+1,1:my+1)%bz+emf(0:mx,1:my+1)%bz+ &
        emf(1:mx+1,0:my)%bz+emf(0:mx,0:my)%bz)*0.25

  ! Calculate the total fields energy
    CALL FieldsEnergy

!  ! Perfect electric conductor (PEC) along x
!    emf(1,:)%ey = 0.0
!    emf(mx+1,:)%ey = 0.0
!    emf(0,:)%ez = -emf(1,:)%ez
!    emf(mx+1,:)%ez = -emf(mx,:)%ez
!  ! Perfect magnetic conductor (PMC) along x
!    emf(1,:)%by = 0.0
!    emf(mx+1,:)%by = 0.0
!    emf(0,:)%bz = -emf(1,:)%bz
!    emf(mx+1,:)%bz = -emf(mx,:)%bz
END SUBROUTINE FieldSolver

!******************************************************************************
! Marder passes for divergence clean.
! Marder, Barry. "A method for incorporating Gauss' law into electromagnetic PIC
! codes." Journal of Computational Physics 68.1 (1987): 48-55.
!******************************************************************************
SUBROUTINE DivergenceClean
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL, DIMENSION(:,:), ALLOCATABLE :: MarderF
    REAL, SAVE :: dconst, dtdx, dtdy
    dconst = 0.001
    dtdx = dt*dconst/dx
    dtdy = dt*dconst/dy
    ALLOCATE(MarderF(0:mx+1,0:my+1))
    MarderF = 0.0
    MarderF(1:mx+1,1:my+1) = emf(1:mx+1,1:my+1)%divE - &
        cva2*denjr(1:mx+1,1:my+1)%rho
    emf(1:mx,1:my+1)%ex = emf(1:mx,1:my+1)%ex + &
        dtdx*(MarderF(2:mx+1,1:my+1)-MarderF(1:mx,1:my+1))
    emf(1:mx+1,1:my)%ey = emf(1:mx+1,1:my)%ey + &
        dtdy*(MarderF(1:mx+1,2:my+1)-MarderF(1:mx+1,1:my))
    DEALLOCATE(MarderF)
END SUBROUTINE DivergenceClean

!******************************************************************************
! Non-physical Cerenkov radiation suppression
!******************************************************************************
SUBROUTINE SuppressCerenkov
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL, SAVE :: tau, dtdx, dtdy
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
END SUBROUTINE SuppressCerenkov

!******************************************************************************
! Push particles
!******************************************************************************
SUBROUTINE pusher
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL :: ux, uy, uz, u0x, u0y, u0z, vz
    REAL :: umx, umy, umz, upx, upy, upz
    REAL :: tx, ty, tz, t2, sx, sy, sz
    REAL :: ex, ey, ez, bx, by, bz
    REAL :: s1, s2, s3, s4
    REAL :: gama, c0, c1, c2, c3, ms
    REAL :: dx0, dy0, dx1, dy1
    INTEGER :: i, j, ix, iy, ix1, iy1 
    INTEGER :: dCellx, dCelly ! Particle movement in cells
    INTEGER :: tag1, tag2, istart, iend, mx1
    INTEGER :: stat(MPI_STATUS_SIZE)
    tag1 = 1
    tag2 = 2

    mx1 = mx + 2
    nptlCrosst = 0
    nptlCrossb = 0
    nptlCross = 0
    DO i = 1, nptlact
        ux = ptl(i)%ux
        uy = ptl(i)%uy
        uz = ptl(i)%uz
        ix = (ptl(i)%grid-1)/my + 1
        iy = MOD(ptl(i)%grid-1,my) + 1

!        IF ((ix .GT. mx) .OR. (iy .GT. my)) THEN
!            PRINT*, 1, taskid, ix, iy, ptl(i)%grid
!        ENDIF
!        IF ((ix .LT. 0) .OR. (iy .LT. 0)) THEN
!            PRINT*, 2, taskid, ix, iy, ptl(i)%grid 
!        ENDIF

        s1 = (0.5+ptl(i)%dx)*(0.5+ptl(i)%dy)
        s2 = 0.5+ptl(i)%dx - s1
        s3 = 0.5+ptl(i)%dy - s1
        s4 = 1.0-s1-s2-s3
        ex = emfGrid(ix,iy)%ex*s4 + emfGrid(ix+1,iy)%ex*s2 + &
             emfGrid(ix,iy+1)%ex*s3 + emfGrid(ix+1,iy+1)%ex*s1
        ey = emfGrid(ix,iy)%ey*s4 + emfGrid(ix+1,iy)%ey*s2 + &
             emfGrid(ix,iy+1)%ey*s3 + emfGrid(ix+1,iy+1)%ey*s1
        ez = emfGrid(ix,iy)%ez*s4 + emfGrid(ix+1,iy)%ez*s2 + &
             emfGrid(ix,iy+1)%ez*s3 + emfGrid(ix+1,iy+1)%ez*s1
        bx = emfGrid(ix,iy)%bx*s4 + emfGrid(ix+1,iy)%bx*s2 + &
             emfGrid(ix,iy+1)%bx*s3 + emfGrid(ix+1,iy+1)%bx*s1
        by = emfGrid(ix,iy)%by*s4 + emfGrid(ix+1,iy)%by*s2 + &
             emfGrid(ix,iy+1)%by*s3 + emfGrid(ix+1,iy+1)%by*s1
        bz = emfGrid(ix,iy)%bz*s4 + emfGrid(ix+1,iy)%bz*s2 + &
             emfGrid(ix,iy+1)%bz*s3 + emfGrid(ix+1,iy+1)%bz*s1
        IF(ptl(i)%q .GT. 0) THEN
            ms = mime
        ELSE
            ms = 1.0
        ENDIF
        c0 = ptl(i)%q*dt*mime/ms/2.0
        c1 = c0/cva
        umx = ux + c1*ex
        umy = uy + c1*ey
        umz = uz + c1*ez
        gama = SQRT(1.0+umx*umx+umy*umy+umz*umz)

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
        gama = SQRT(1.0+ux*ux+uy*uy+uz*uz)

        ptl(i)%ux = ux
        ptl(i)%uy = uy
        ptl(i)%uz = uz

        dx0 = ptl(i)%dx
        dy0 = ptl(i)%dy
        c3 = dt*cva/gama
        dx1 = dx0 + c3*ux*idx
        dy1 = dy0 + c3*uy*idy

        IF(ABS(dx1) .GT. 0.5) THEN
            dCellx = NINT(SIGN(1.0,dx1))
            !ptl(i)%dx = SIGN(ABS(dx1)-1.0,-dx1)
            ptl(i)%dx = dx1 + SIGN(1.0,-dx1)
        ELSE
            dCellx = 0
            ptl(i)%dx = dx1
        ENDIF
        IF(ABS(dy1) .GT. 0.5) THEN
            dCelly = NINT(SIGN(1.0,dy1))
            !ptl(i)%dy = SIGN(ABS(dy1)-1.0,-dy1)
            ptl(i)%dy = dy1 + SIGN(1.0,-dy1)
        ELSE
            dCelly = 0
            ptl(i)%dy = dy1
        ENDIF

        ix1 = ix + dCellx
        iy1 = iy + dCelly

        vz = uz*cva/gama
        CALL CurrentDeposit(ix, ix1, iy, iy1, &
            dx0, ptl(i)%dx, dy0, ptl(i)%dy, ptl(i)%q, vz)

        ix = MOD(ix1+mx-1,mx) + 1   ! Periodic along x
        IF (iy1 .GT. my) THEN
            iy = iy1 - my
            ptl(i)%grid = my*(ix-1) + iy 
            nptlCrosst = nptlCrosst + 1
            nptlCross = nptlCross + 1
            ptlsendt(nptlCrosst) = ptl(i)
            !PRINT*, nptlCrosst, ptlsendt(nptlCrosst)%test
            ptl(i)%pid = nptltot+1 ! Flag for crossing particle
            CrossIndex(nptlCross) = i
        ELSE IF( iy1 .LT. 1) THEN
            iy = iy1 + my
            ptl(i)%grid = my*(ix-1) + iy 
            nptlCrossb = nptlCrossb + 1
            nptlCross = nptlCross + 1
            ptlsendb(nptlCrossb) = ptl(i)
            ptl(i)%pid = nptltot+1
            CrossIndex(nptlCross) = i
        ELSE
            iy = iy1
            ptl(i)%grid = my*(ix-1) + iy 
        ENDIF

      ! Total kinetic energy for ions and electrons
!        IF (ptl(i)%q .EQ. 1) THEN
!            j = 2
!        ELSE
!            j = 1
!        ENDIF
        j = CEILING(ptl(i)%q*0.5+1.0)
        totene(j,itime) = totene(j,itime) + (gama-1.0)*ms
    ENDDO

  ! Normalized current density and charge density
    denjr%jx = denjr%jx * in0
    denjr%jy = denjr%jy * in0
    denjr%jz = denjr%jz * in0
    denjr%rho = denjr%rho * in0

  ! Reuse the memory of the particles crossing the domain boundaries
    j = nptlact
    DO i = 1, nptlCross
        DO WHILE (ptl(j)%pid .GT. nptltot)
            j = j - 1
        ENDDO
        IF(j .LT. CrossIndex(i)) THEN
            EXIT
        ENDIF
        ptl(CrossIndex(i)) = ptl(j)
        j = j - 1
    ENDDO
    nptlact = nptlact - nptlCross

  ! Send and receive the particles numbers first
    IF (MOD(taskid,2) .EQ. 0) THEN
        CALL MPI_SEND(nptlCrosst, 1, MPI_INTEGER, &
            right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(nptlCrossb, 1, MPI_INTEGER, &
            left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(nptlRecvb, 1, MPI_INTEGER, &
            left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(nptlRecvt, 1, MPI_INTEGER, &
            right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF
    IF (MOD(taskid,2) .EQ. 1) THEN
        CALL MPI_SEND(nptlCrosst, 1, MPI_INTEGER, &
            right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(nptlCrossb, 1, MPI_INTEGER, &
            left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(nptlRecvb, 1, MPI_INTEGER, &
            left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(nptlRecvt, 1, MPI_INTEGER, &
            right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF

    !PRINT*, taskid, itime, nptlCrosst, nptlCrossb, nptlRecvt, nptlRecvb

  ! Exchange particles with neighbours
    IF (MOD(taskid,2) .EQ. 0) THEN
        CALL MPI_SEND(ptlsendt(1:nptlCrosst), nptlCrosst, &
            particletype, right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(ptlsendb(1:nptlCrossb), nptlCrossb, &
            particletype, left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(ptlrecvb(1:nptlRecvb), nptlRecvb, &
            particletype, left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(ptlrecvt(1:nptlRecvt), nptlRecvt, &
            particletype, right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF

    IF (MOD(taskid,2) .EQ. 1) THEN
        CALL MPI_SEND(ptlsendt(1:nptlCrosst), nptlCrosst, &
            particletype, right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(ptlsendb(1:nptlCrossb), nptlCrossb, &
            particletype, left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(ptlrecvb(1:nptlRecvb), nptlRecvb, &
            particletype, left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(ptlrecvt(1:nptlRecvt), nptlRecvt, &
            particletype, right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF

    !PRINT*, taskid, itime, nptlCrosst, nptlCrossb, nptlRecvt, nptlRecvb

    istart = nptlact + 1
    iend = nptlact + nptlRecvt
    ptl(istart:iend) = ptlrecvt(1:nptlRecvt)
    nptlact = nptlact + nptlRecvt

    istart = nptlact + 1
    iend = nptlact + nptlRecvb
    ptl(istart:iend) = ptlrecvb(1:nptlRecvb)
    nptlact = nptlact + nptlRecvb

  ! Electric current density and charge density information exchange
    ! Periodic along x
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

    IF (MOD(taskid,2) .EQ. 0) THEN
        CALL MPI_SEND(densendt, densize, MPI_REAL, &
            right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(densendb, densize, MPI_REAL, &
            left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(denrecvb, densize, MPI_REAL, &
            left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(denrecvt, densize, MPI_REAL, &
            right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF
    IF (MOD(taskid,2) .EQ. 1) THEN
        CALL MPI_SEND(densendt, densize, MPI_REAL, &
            right, tag1, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(densendb, densize, MPI_REAL, &
            left, tag2, MPI_COMM_WORLD, ierr)
    ELSE
        CALL MPI_RECV(denrecvb, densize, MPI_REAL, &
            left, tag1, MPI_COMM_WORLD, stat, ierr)
        CALL MPI_RECV(denrecvt, densize, MPI_REAL, &
            right, tag2, MPI_COMM_WORLD, stat, ierr)
    ENDIF

    !PRINT*, taskid, itime, nptlCrosst, nptlCrossb, nptlRecvt, nptlRecvb

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
END SUBROUTINE pusher

!******************************************************************************
! Current deposit to get current density after pushing particles
! Zigzag scheme from
! Umeda, Takayuki, et al. "A new charge conservation method in electromagnetic
! particle-in-cell simulations." Computer Physics Communications 156.1 (2003):
! 73-85.
!******************************************************************************
SUBROUTINE CurrentDeposit(i1,i2,j1,j2,dx1,dx2,dy1,dy2,q0, uz)
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER, INTENT(IN) :: i1, i2, j1, j2, q0
    REAL, INTENT(IN) :: dx1, dx2, dy1, dy2, uz
    REAL, SAVE :: fx1, fx2, fy1, fy2, wx1, wx2, wy1, wy2
    REAL, SAVE :: ir, jr, im, jm, diMid, djMid, qx, qy
    REAL, SAVE :: s1, s2, s3, s4, quz
    INTEGER, SAVE :: di, dj, ix1, iy1
  
    qx = q0*dx/dt
    qy = q0*dy/dt
    im = (i1+i2+dx1+dx2+1)*0.5
    jm = (j1+j2+dy1+dy2+1)*0.5
    ir = MIN(MIN(i1,i2)+1.0,MAX(MAX(i1,i2)+0.0,im))
    jr = MIN(MIN(j1,j2)+1.0,MAX(MAX(j1,j2)+0.0,jm))

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

    diMid = im - i1 - 0.5
    djMid = jm - j1 - 0.5

    di = NINT(SIGN(1.0,diMid))
    dj = NINT(SIGN(1.0,djMid))

    diMid = ABS(diMid)
    djMid = ABS(djMid)

    s1 = diMid*djMid
    s2 = diMid - s1
    s3 = djMid - s1
    s4 = 1.0-s1-s2-s3

    quz = q0*uz
    denjr(i1,j1)%jz = denjr(i1,j1)%jz + quz*s4
    denjr(i1+di,j1)%jz = denjr(i1+di,j1)%jz + quz*s2
    denjr(i1,j1+dj)%jz = denjr(i1,j1+dj)%jz + quz*s3
    denjr(i1+di,j1+dj)%jz = denjr(i1+di,j1+dj)%jz + quz*s1

  ! Charge density accumulation with CIC particle shape
  ! The density is accumulated at time n
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
END SUBROUTINE CurrentDeposit

!******************************************************************************
! Counting sort for particles to accelerate the simulation
! Reference: Bowers, K. J. "Accelerating a particle-in-cell simulation using a
! hybrid counting sort." Journal of Computational Physics 173.2 (2001): 393-411.
!******************************************************************************
SUBROUTINE CountingSort
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER :: i, j, k
    ptlAlloc = 0
  ! Count the number of particles in each cell
    DO i = 1, nptlact
        j = ptl(i)%grid
        ptlAlloc(j) = ptlAlloc(j) + 1
    ENDDO
  ! Convert ptlAlloc to an allocation
    k = 0
    DO i = 1, mx*my
        j = ptlAlloc(i)
        ptlAlloc(i) = k
        k = k + j
    ENDDO
  ! Sort particle information into ptlSort
    Do i = 1, nptlact
        j = ptl(i)%grid
        k = ptlAlloc(j) + 1
        ptlAlloc(j) = k
        ptlSort(k) = ptl(i)
    ENDDO
    ptl(1:nptlact) = ptlSort(1:nptlact)
END SUBROUTINE CountingSort

!******************************************************************************
! Diagnostics for particle, fields and energy
!******************************************************************************
SUBROUTINE diagnostics
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
END SUBROUTINE diagnostics

!******************************************************************************
! Particle number density
!******************************************************************************
SUBROUTINE NumDensity
    USE parameters
    IMPLICIT NONE
    INTEGER :: ix, iy, n
    numdens = 0
    DO n = 1, nptlact
        ix = (ptl(n)%grid-1)/my + 1
        iy = MOD(ptl(n)%grid-1,my) + 1
        numdens(ix,iy) = numdens(ix,iy) + 1
    ENDDO
    !numdens = numdens * in0
END SUBROUTINE NumDensity

!******************************************************************************
! Diagnostics for fields (E, B, J, rho)
!******************************************************************************
SUBROUTINE ParticlesInfo(isave)
    USE parameters
    USE hdf5
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER, INTENT(IN) :: isave
    CHARACTER(LEN=14) :: fname ! File name
    CHARACTER(LEN=2), PARAMETER :: dx_id = "dx"
    CHARACTER(LEN=2), PARAMETER :: dy_id = "dy"
    CHARACTER(LEN=2), PARAMETER :: ux_id = "ux"
    CHARACTER(LEN=2), PARAMETER :: uy_id = "uy"
    CHARACTER(LEN=2), PARAMETER :: uz_id = "uz"
    CHARACTER(LEN=4), PARAMETER :: grid_id = "grid"
    CHARACTER(LEN=1), PARAMETER :: q_id = "q"
    CHARACTER(LEN=4), PARAMETER :: pid_id = "pid"
    CHARACTER(LEN=7), PARAMETER :: dsetnptl="nptlact"
    CHARACTER(LEN=11), PARAMETER :: dsetnpre="nptltot_pre"

    INTEGER(HID_T) :: file_id   ! File identifier
    INTEGER(HID_T) :: dset_dx, dset_dy
    INTEGER(HID_T) :: dset_ux, dset_uy, dset_uz
    INTEGER(HID_T) :: dset_grid, dset_q, dset_pid
    INTEGER(HID_T) :: filespace ! Dataspace identifier in file
    INTEGER(HID_T) :: memspace  ! Dataspace identifier in memory
    INTEGER(HID_T) :: plist_id  ! Property list identifier

    INTEGER(HID_T) :: dset_nptlact
    INTEGER(HID_T) :: dspace_nptlact
    INTEGER(HID_T) :: mem_nptlact

    INTEGER(HSIZE_T), DIMENSION(1) :: dimsf, dimsfi, dimsna

    INTEGER(HSIZE_T), DIMENSION(1) :: count
    INTEGER(HSIZE_T), DIMENSION(1) :: offset
    INTEGER :: rank = 1 ! Dataset rank

    INTEGER, DIMENSION(:), ALLOCATABLE :: nptlarray
    INTEGER :: error, info, comm, i, tag1
    INTEGER :: stat(MPI_STATUS_SIZE)
    comm = MPI_COMM_WORLD
!    info = MPI_INFO_NULL
    CALL MPI_INFO_CREATE(info, ierr)
  ! Disable ROMIO's data-sieving
    CALL MPI_INFO_SET(info, "romio_ds_read", "disable", ierr)
    CALL MPI_INFO_SET(info, "romio_ds_write", "disable", ierr)
  ! Enable ROMIO's collective buffering
    CALL MPI_INFO_SET(info, "romio_cb_read", "enable", ierr)
    CALL MPI_INFO_SET(info, "romio_cb_write", "enable", ierr)

    tag1 = 1

    dimsf = (/ppg*nx*ny*2/)
    dimsfi = (/ppg*nx*ny*2/)
    dimsna = (/numtasks/)

    WRITE(fname, "(A12)") "particles.h5"
  ! Initialize FORTRAN predefined data types
    CALL h5open_f(error)

  ! Setup file access property list with parallel I/O access
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
    CALL MPI_INFO_FREE(info, ierr)

    IF (isave .EQ. 1) THEN
      ! Create the file collectively
        CALL h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
    ELSE IF (isave .EQ. 0) THEN
      ! Read the file collectively
        CALL h5fopen_f (fname, H5F_ACC_RDONLY_F, file_id, error, access_prp=plist_id)
    ENDIF
    CALL h5pclose_f(plist_id, error)

    IF (isave .EQ. 1) THEN
      ! Create the data space for the dataset
        CALL h5screate_simple_f(rank, dimsf, filespace, error)
      ! Create the dataset with default properties
        CALL h5dcreate_f(file_id, dx_id, H5T_NATIVE_REAL, filespace, &
                         dset_dx, error)
        CALL h5dcreate_f(file_id, dy_id, H5T_NATIVE_REAL, filespace, &
                         dset_dy, error)
        CALL h5dcreate_f(file_id, ux_id, H5T_NATIVE_REAL, filespace, &
                         dset_ux, error)
        CALL h5dcreate_f(file_id, uy_id, H5T_NATIVE_REAL, filespace, &
                         dset_uy, error)
        CALL h5dcreate_f(file_id, uz_id, H5T_NATIVE_REAL, filespace, &
                         dset_uz, error)
        CALL h5dcreate_f(file_id, grid_id, H5T_NATIVE_INTEGER, filespace, &
                         dset_grid, error)
        CALL h5dcreate_f(file_id, q_id, H5T_NATIVE_INTEGER, filespace, &
                         dset_q, error)
        CALL h5dcreate_f(file_id, pid_id, H5T_NATIVE_INTEGER, filespace, &
                         dset_pid, error)
    ELSE IF (isave .EQ. 0) THEN
        CALL h5dopen_f(file_id, dx_id, dset_dx, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, dy_id, dset_dy, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, ux_id, dset_ux, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, uy_id, dset_uy, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, uz_id, dset_uz, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, grid_id, dset_grid, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, q_id, dset_q, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, pid_id, dset_pid, error, dapl_id=H5P_DEFAULT_F)
    ENDIF

    IF (isave .EQ. 1) THEN
      ! Create the data space for the dataset of actual particle numbers
        CALL h5screate_simple_f(rank, dimsna, dspace_nptlact, error)
      ! Create the dataset
        CALL h5dcreate_f(file_id, dsetnptl, H5T_NATIVE_INTEGER, &
            dspace_nptlact, dset_nptlact, error)
    ELSE IF (isave .EQ. 0) THEN
        CALL h5dopen_f(file_id, dsetnptl, dset_nptlact, error, dapl_id=H5P_DEFAULT_F)
    ENDIF

  ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  ! Each process defines dataset in memory and writes it to the 
  ! hyperslab in the file
    count(1) = 1
    offset(1) = taskid
    CALL h5screate_simple_f(rank, count, mem_nptlact, error)

  ! Select hyperslab in the file
    CALL h5dget_space_f(dset_nptlact, dspace_nptlact, error)
    CALL h5sselect_hyperslab_f(dspace_nptlact, H5S_SELECT_SET_F, offset, &
        count, error)
    
    IF (isave .EQ. 1) THEN
        CALL h5dwrite_f(dset_nptlact, H5T_NATIVE_INTEGER, nptlact, dimsna, &
            error, file_space_id=dspace_nptlact, mem_space_id=mem_nptlact, &
            xfer_prp=plist_id)
    ELSE IF(isave .EQ. 0) THEN
        CALL h5dread_f(dset_nptlact, H5T_NATIVE_INTEGER, nptlact, dimsna, &
            error, file_space_id=dspace_nptlact, mem_space_id=mem_nptlact, &
            xfer_prp=plist_id)
    ENDIF

  ! Each process defines dataset in memory and writes it to the hyperslab
  ! in the file

    IF (taskid .EQ. 0) THEN
        ALLOCATE(nptlarray(0:numtasks-1))
        nptlarray(0) = nptlact
    ENDIF
    
    IF (taskid .NE. 0) THEN
        CALL MPI_SEND(nptlact, 1, MPI_INTEGER, 0, tag1, comm, ierr)
    ELSE
        DO i = 1, numtasks-1
            CALL MPI_RECV(nptlarray(i), 1, MPI_INTEGER, &
                i, tag1, comm, stat, ierr)
            nptlarray(i) = nptlarray(i) + nptlarray(i-1)
        ENDDO
    ENDIF

    IF (taskid .EQ. 0) THEN
        DO i = 1, numtasks-1
            CALL MPI_SEND(nptlarray(i-1), 1, MPI_INTEGER, i, tag1, comm, ierr)
        ENDDO
        offset(1) = 0
    ELSE
        CALL MPI_RECV(offset(1), 1, MPI_INTEGER, 0, tag1, comm, stat, ierr)
    ENDIF
   
    CALL MPI_Barrier(comm, ierr)

    IF (taskid .EQ. 0) THEN
        DEALLOCATE(nptlarray)
    ENDIF

    count(1) = nptlact
    CALL h5screate_simple_f(rank, count, memspace, error)

  ! Select hyperslab in the file
    CALL h5dget_space_f(dset_dx, filespace, error)
    CALL h5dget_space_f(dset_dy, filespace, error)
    CALL h5dget_space_f(dset_ux, filespace, error)
    CALL h5dget_space_f(dset_uy, filespace, error)
    CALL h5dget_space_f(dset_uz, filespace, error)
    CALL h5dget_space_f(dset_grid, filespace, error)
    CALL h5dget_space_f(dset_q, filespace, error)
    CALL h5dget_space_f(dset_pid, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
        count, error)

  ! Write the dataset collectively
    IF (isave .EQ. 1) THEN
        CALL h5dwrite_f(dset_dx, H5T_NATIVE_REAL, ptl(1:nptlact)%dx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_dy, H5T_NATIVE_REAL, ptl(1:nptlact)%dy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_ux, H5T_NATIVE_REAL, ptl(1:nptlact)%ux, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_uy, H5T_NATIVE_REAL, ptl(1:nptlact)%uy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_uz, H5T_NATIVE_REAL, ptl(1:nptlact)%uz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_grid, H5T_NATIVE_INTEGER, ptl(1:nptlact)%grid, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_q, H5T_NATIVE_INTEGER, ptl(1:nptlact)%q, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_pid, H5T_NATIVE_INTEGER, ptl(1:nptlact)%pid, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
    ELSE IF(isave .EQ. 0) THEN
        CALL h5dread_f(dset_dx, H5T_NATIVE_REAL, ptl(1:nptlact)%dx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_dy, H5T_NATIVE_REAL, ptl(1:nptlact)%dy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_ux, H5T_NATIVE_REAL, ptl(1:nptlact)%ux, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_uy, H5T_NATIVE_REAL, ptl(1:nptlact)%uy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_uz, H5T_NATIVE_REAL, ptl(1:nptlact)%uz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_grid, H5T_NATIVE_INTEGER, ptl(1:nptlact)%grid, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_q, H5T_NATIVE_INTEGER, ptl(1:nptlact)%q, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_pid, H5T_NATIVE_INTEGER, ptl(1:nptlact)%pid, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
    ENDIF

  ! Close dataspaces
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(dspace_nptlact, error)
    CALL h5sclose_f(memspace, error)
    CALL h5sclose_f(mem_nptlact, error)

  ! Close the dataset and property list
    CALL h5dclose_f(dset_dx, error)
    CALL h5dclose_f(dset_dy, error)
    CALL h5dclose_f(dset_ux, error)
    CALL h5dclose_f(dset_uy, error)
    CALL h5dclose_f(dset_uz, error)
    CALL h5dclose_f(dset_grid, error)
    CALL h5dclose_f(dset_q, error)
    CALL h5dclose_f(dset_pid, error)
    CALL h5pclose_f(plist_id, error)
    CALL h5dclose_f(dset_nptlact, error)

  ! Close the file
    CALL h5fclose_f(file_id, error)
  ! Close FORTRAN interface
    CALL h5close_f(error)
END SUBROUTINE ParticlesInfo

!******************************************************************************
! Diagnostics for fields (E, B, J, rho)
!******************************************************************************
SUBROUTINE FieldsInfo(isave)
    USE parameters
    USE hdf5
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER, INTENT(IN) :: isave
    CHARACTER(LEN=14) :: fname ! File name
    CHARACTER(LEN=2), PARAMETER :: Ex = "Ex"
    CHARACTER(LEN=2), PARAMETER :: Ey = "Ey"
    CHARACTER(LEN=2), PARAMETER :: Ez = "Ez"
    CHARACTER(LEN=2), PARAMETER :: Bx = "Bx"
    CHARACTER(LEN=2), PARAMETER :: By = "By"
    CHARACTER(LEN=2), PARAMETER :: Bz = "Bz"
    CHARACTER(LEN=2), PARAMETER :: Jx = "Jx"
    CHARACTER(LEN=2), PARAMETER :: Jy = "Jy"
    CHARACTER(LEN=2), PARAMETER :: Jz = "Jz"
    CHARACTER(LEN=3), PARAMETER :: Rho = "Rho"
    CHARACTER(LEN=3), PARAMETER :: Num = "Num"
    CHARACTER(LEN=4), PARAMETER :: divE = "divE"
    INTEGER(HID_T) :: file_id   ! File identifier
    INTEGER(HID_T) :: dset_ex, dset_ey, dset_ez
    INTEGER(HID_T) :: dset_bx, dset_by, dset_bz
    INTEGER(HID_T) :: dset_jx, dset_jy, dset_jz
    INTEGER(HID_T) :: dset_rho, dset_num, dset_divE
    INTEGER(HID_T) :: filespace ! Dataspace identifier in file
    INTEGER(HID_T) :: memspace  ! Dataspace identifier in memory
    INTEGER(HID_T) :: plist_id  ! Property list identifier

    INTEGER(HSIZE_T), DIMENSION(2) :: dimsf = (/nx+2,ny+1/) ! Dataset dimensions.
    INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi = (/nx+2,ny+1/)

    INTEGER(HSIZE_T), DIMENSION(2) :: count
    INTEGER(HSIZE_T), DIMENSION(2) :: offset
    INTEGER :: rank = 2 ! Dataset rank

    INTEGER :: error, info, comm, jend, jstart
    INTEGER :: tag1, tag2
    INTEGER :: stat(MPI_STATUS_SIZE)
    comm = MPI_COMM_WORLD
!    info = MPI_INFO_NULL
    CALL MPI_INFO_CREATE(info, ierr)
  ! Disable ROMIO's data-sieving
    CALL MPI_INFO_SET(info, "romio_ds_read", "disable", ierr)
    CALL MPI_INFO_SET(info, "romio_ds_write", "disable", ierr)
  ! Enable ROMIO's collective buffering
    CALL MPI_INFO_SET(info, "romio_cb_read", "enable", ierr)
    CALL MPI_INFO_SET(info, "romio_cb_write", "enable", ierr)

    tag1 = 1
    tag2 = 2

    jend = 0
    jstart = 0

    WRITE(fname, "(A6,I5.5,A3)") "fields", itime, ".h5"
  ! Initialize FORTRAN predefined data types
    CALL h5open_f(error)

  ! Setup file access property list with parallel I/O access
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
    CALL MPI_INFO_FREE(info, ierr)

    IF (isave .EQ. 1) THEN
      ! Create the file collectively
        CALL h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
    ELSE IF (isave .EQ. 0) THEN
      ! Read the file collectively
        CALL h5fopen_f (fname, H5F_ACC_RDONLY_F, file_id, error, access_prp=plist_id)
    ENDIF
    CALL h5pclose_f(plist_id, error)

    IF (isave .EQ. 1) THEN
      ! Create the data space for the dataset
        CALL h5screate_simple_f(rank, dimsf, filespace, error)
      ! Create the dataset with default properties
        CALL h5dcreate_f(file_id, Ex, H5T_NATIVE_REAL, filespace, &
                         dset_ex, error)
        CALL h5dcreate_f(file_id, Ey, H5T_NATIVE_REAL, filespace, &
                         dset_ey, error)
        CALL h5dcreate_f(file_id, Ez, H5T_NATIVE_REAL, filespace, &
                         dset_ez, error)
        CALL h5dcreate_f(file_id, Bx, H5T_NATIVE_REAL, filespace, &
                         dset_bx, error)
        CALL h5dcreate_f(file_id, By, H5T_NATIVE_REAL, filespace, &
                         dset_by, error)
        CALL h5dcreate_f(file_id, Bz, H5T_NATIVE_REAL, filespace, &
                         dset_bz, error)
        CALL h5dcreate_f(file_id, Jx, H5T_NATIVE_REAL, filespace, &
                         dset_jx, error)
        CALL h5dcreate_f(file_id, Jy, H5T_NATIVE_REAL, filespace, &
                         dset_jy, error)
        CALL h5dcreate_f(file_id, Jz, H5T_NATIVE_REAL, filespace, &
                         dset_jz, error)
        CALL h5dcreate_f(file_id, Rho, H5T_NATIVE_REAL, filespace, &
                         dset_rho, error)
        CALL h5dcreate_f(file_id, Num, H5T_NATIVE_INTEGER, filespace, &
                         dset_num, error)
        CALL h5dcreate_f(file_id, divE, H5T_NATIVE_REAL, filespace, &
                         dset_divE, error)
        CALL h5sclose_f(filespace, error)
    ELSE IF (isave .EQ. 0) THEN
        CALL h5dopen_f(file_id, Ex, dset_ex, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, Ey, dset_ey, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, Ez, dset_ez, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, Bx, dset_bx, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, By, dset_by, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, Bz, dset_bz, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, Jx, dset_jx, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, Jy, dset_jy, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, Jz, dset_jz, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, Rho, dset_rho, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, Num, dset_num, error, dapl_id=H5P_DEFAULT_F)
        CALL h5dopen_f(file_id, divE, dset_divE, error, dapl_id=H5P_DEFAULT_F)
    ENDIF

    IF (isave .EQ. 1) THEN
      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file
        count(1) = mx+2
        IF (taskid .EQ. numtasks-1) THEN
            jend = my+1
        ELSE
            jend = my
        ENDIF
        count(2) = jend
        offset(1) = 0
        offset(2) = taskid * my
    ELSE IF (isave .EQ. 0) THEN
        count(1) = mx+2
        offset(1) = 0
        IF (taskid .EQ. 0) THEN
            jstart = 1
            offset(2) = taskid*my
        ELSE
            jstart = 0
            offset(2) = taskid*my-1
        ENDIF
        count(2) = my+2-jstart
    ENDIF
    CALL h5screate_simple_f(rank, count, memspace, error)

  ! Select hyperslab in the file
    CALL h5dget_space_f(dset_ex, filespace, error)
    CALL h5dget_space_f(dset_ey, filespace, error)
    CALL h5dget_space_f(dset_ez, filespace, error)
    CALL h5dget_space_f(dset_bx, filespace, error)
    CALL h5dget_space_f(dset_by, filespace, error)
    CALL h5dget_space_f(dset_bz, filespace, error)
    CALL h5dget_space_f(dset_jx, filespace, error)
    CALL h5dget_space_f(dset_jy, filespace, error)
    CALL h5dget_space_f(dset_jz, filespace, error)
    CALL h5dget_space_f(dset_rho, filespace, error)
    CALL h5dget_space_f(dset_num, filespace, error)
    CALL h5dget_space_f(dset_divE, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
        count, error)
  
  ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  ! Write the dataset collectively
    IF (isave .EQ. 1) THEN
        CALL NumDensity ! Get particle number density

        CALL h5dwrite_f(dset_ex, H5T_NATIVE_REAL, emf(0:mx+1,1:jend)%ex, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_ey, H5T_NATIVE_REAL, emf(0:mx+1,1:jend)%ey, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_ez, H5T_NATIVE_REAL, emf(0:mx+1,1:jend)%ez, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_bx, H5T_NATIVE_REAL, emf(0:mx+1,1:jend)%bx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_by, H5T_NATIVE_REAL, emf(0:mx+1,1:jend)%by, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_bz, H5T_NATIVE_REAL, emf(0:mx+1,1:jend)%bz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_divE, H5T_NATIVE_REAL, emf(0:mx+1,1:jend)%divE, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_jx, H5T_NATIVE_REAL, denjr(0:mx+1,1:jend)%jx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_jy, H5T_NATIVE_REAL, denjr(0:mx+1,1:jend)%jy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_jz, H5T_NATIVE_REAL, denjr(0:mx+1,1:jend)%jz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_rho, H5T_NATIVE_REAL, denjr(0:mx+1,1:jend)%rho, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dwrite_f(dset_num, H5T_NATIVE_INTEGER, numdens(0:mx+1,1:jend), &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
    ELSE IF(isave .EQ. 0) THEN
        CALL h5dread_f(dset_ex, H5T_NATIVE_REAL, emf(0:mx+1,jstart:my+1)%ex, &
            dimsfi, error, mem_space_id=memspace, file_space_id=filespace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_ey, H5T_NATIVE_REAL, emf(0:mx+1,jstart:my+1)%ey, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_ez, H5T_NATIVE_REAL, emf(0:mx+1,jstart:my+1)%ez, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_bx, H5T_NATIVE_REAL, emf(0:mx+1,jstart:my+1)%bx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_by, H5T_NATIVE_REAL, emf(0:mx+1,jstart:my+1)%by, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_bz, H5T_NATIVE_REAL, emf(0:mx+1,jstart:my+1)%bz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_divE, H5T_NATIVE_REAL, emf(0:mx+1,jstart:my+1)%divE, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_jx, H5T_NATIVE_REAL, denjr(0:mx+1,jstart:my+1)%jx, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_jy, H5T_NATIVE_REAL, denjr(0:mx+1,jstart:my+1)%jy, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_jz, H5T_NATIVE_REAL, denjr(0:mx+1,jstart:my+1)%jz, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_rho, H5T_NATIVE_REAL, denjr(0:mx+1,jstart:my+1)%rho, &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)
        CALL h5dread_f(dset_num, H5T_NATIVE_INTEGER, numdens(0:mx+1,jstart:my+1), &
            dimsfi, error, file_space_id=filespace, mem_space_id=memspace, &
            xfer_prp=plist_id)

      ! Exchange the information of By, Bz, Ey, Ez for taskid=0 and numtasks-1
        IF (taskid .EQ. 0) THEN
            bufsendb(1:mx+2) = emf(:,1)%by
            bufsendb(mx+3:bufsize) = emf(:,1)%bz
            CALL MPI_SEND(bufsendb, bufsize, MPI_REAL, &
                left, tag2, MPI_COMM_WORLD, ierr)
            CALL MPI_RECV(bufrecvb, bufsize, MPI_REAL, &
                left, tag1, MPI_COMM_WORLD, stat, ierr)
            emf(:,0)%by = bufrecvb(1:mx+2)
            emf(:,0)%bz = bufrecvb(mx+3:bufsize)
        ELSE IF (taskid .EQ. numtasks-1) THEN
            bufsendt(1:mx+2) = emf(:,my)%by
            bufsendt(mx+3:bufsize) = emf(:,my)%bz
            CALL MPI_RECV(bufrecvt, bufsize, MPI_REAL, &
                right, tag2, MPI_COMM_WORLD, stat, ierr)
            CALL MPI_SEND(bufsendt, bufsize, MPI_REAL, &
                right, tag1, MPI_COMM_WORLD, ierr)
            emf(:,my+1)%by = bufrecvt(1:mx+2)
            emf(:,my+1)%bz = bufrecvt(mx+3:bufsize)
        ENDIF
        IF (taskid .EQ. 0) THEN
            bufsendb(1:mx+2) = emf(:,1)%ey
            bufsendb(mx+3:bufsize) = emf(:,1)%ez
            CALL MPI_SEND(bufsendb, bufsize, MPI_REAL, &
                left, tag2, MPI_COMM_WORLD, ierr)
            CALL MPI_RECV(bufrecvb, bufsize, MPI_REAL, &
                left, tag1, MPI_COMM_WORLD, stat, ierr)
            emf(:,0)%ey = bufrecvb(1:mx+2)
            emf(:,0)%ez = bufrecvb(mx+3:bufsize)
        ELSE IF (taskid .EQ. numtasks-1) THEN
            bufsendt(1:mx+2) = emf(:,my)%ey
            bufsendt(mx+3:bufsize) = emf(:,my)%ez
            CALL MPI_RECV(bufrecvt, bufsize, MPI_REAL, &
                right, tag2, MPI_COMM_WORLD, stat, ierr)
            CALL MPI_SEND(bufsendt, bufsize, MPI_REAL, &
                right, tag1, MPI_COMM_WORLD, ierr)
            emf(:,my+1)%ey = bufrecvt(1:mx+2)
            emf(:,my+1)%ez = bufrecvt(mx+3:bufsize)
        ENDIF

        emf(:,1:my)%dbx_dy = (emf(:,2:my+1)%bx-emf(:,1:my)%bx)*idy
        emf(0:mx,:)%dby_dx = (emf(1:mx+1,:)%by-emf(0:mx,:)%by)*idx
        emf(1:mx+1,:)%dbz_dx = (emf(1:mx+1,:)%bz-emf(0:mx,:)%bz)*idx
        emf(:,1:my+1)%dbz_dy = (emf(:,1:my+1)%bz-emf(:,0:my)%bz)*idy
        emf(:,1:my)%dex_dy = (emf(:,2:my+1)%ex-emf(:,1:my)%ex)*idy
        emf(0:mx,:)%dey_dx = (emf(1:mx+1,:)%ey-emf(0:mx,:)%ey)*idx
        emf(1:mx+1,:)%dez_dx = (emf(1:mx+1,:)%ez-emf(0:mx,:)%ez)*idx
        emf(:,1:my+1)%dez_dy = (emf(:,1:my+1)%ez-emf(:,0:my)%ez)*idy
        emf(1:mx+1,1:my+1)%divE = (emf(1:mx+1,1:my+1)%ex-emf(0:mx,1:my+1)%ex)*idx + &
                                  (emf(1:mx+1,1:my+1)%ey-emf(1:mx+1,0:mx)%ey)*idy
      ! Periodic along x
        emf(mx+1,:)%dby_dx = emf(1,:)%dby_dx
        emf(mx+1,:)%dey_dx = emf(1,:)%dey_dx

        emfGrid(1:mx+1,:)%ex = (emf(1:mx+1,:)%ex+emf(0:mx,:)%ex)*0.5
        emfGrid(:,1:my+1)%ey = (emf(:,1:my+1)%ey+emf(:,0:my)%ey)*0.5
        emfGrid(1:mx+1,1:my+1)%ez = (emf(1:mx+1,1:my+1)%ez+emf(0:mx,1:my+1)%ez+ &
            emf(1:mx+1,0:my)%ez+emf(0:mx,0:my)%ez)*0.25
        emfGrid(1:mx+1,:)%bx = (emf(1:mx+1,:)%bx+emf(0:mx,:)%bx)*0.5
        emfGrid(:,1:my+1)%by = (emf(:,1:my+1)%by+emf(:,0:my)%by)*0.5
        emfGrid(1:mx+1,1:my+1)%bz = (emf(1:mx+1,1:my+1)%bz+emf(0:mx,1:my+1)%bz+ &
            emf(1:mx+1,0:my)%bz+emf(0:mx,0:my)%bz)*0.25
    ENDIF
  ! Close dataspaces
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)

  ! Close the dataset and property list
    CALL h5dclose_f(dset_ex, error)
    CALL h5dclose_f(dset_ey, error)
    CALL h5dclose_f(dset_ez, error)
    CALL h5dclose_f(dset_bx, error)
    CALL h5dclose_f(dset_by, error)
    CALL h5dclose_f(dset_bz, error)
    CALL h5dclose_f(dset_jx, error)
    CALL h5dclose_f(dset_jy, error)
    CALL h5dclose_f(dset_jz, error)
    CALL h5dclose_f(dset_rho, error)
    CALL h5dclose_f(dset_num, error)
    CALL h5dclose_f(dset_divE, error)
    CALL h5pclose_f(plist_id, error)

  ! Close the file
    CALL h5fclose_f(file_id, error)
  ! Close FORTRAN interface
    CALL h5close_f(error)
END SUBROUTINE FieldsInfo

!******************************************************************************
! Accumulate to get the energy spectra for ions and electrons.
!******************************************************************************
SUBROUTINE EnergySpectra
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=doublep) :: ux, uy, uz, ms, gama, de, maxe, ene
    INTEGER :: i, ibin, q, j
    espectra = 0
    DO i = 1, nptlact
        ux = ptl(i)%ux
        uy = ptl(i)%uy
        uz = ptl(i)%uz
        q = ptl(i)%q
        gama = SQRT(1.0+ux*ux+uy*uy+uz*uz)
        IF (q .EQ. 1) THEN
          ! Ion
            ms = mime
            de = dei
            maxe = maxenei
            j = 2
        ELSE
          ! Electron
            ms = 1.0
            de = dee
            maxe = maxenee
            j = 1
        ENDIF
        
        ene = (gama-1.0)*ms
        IF (ene .LT. maxe) THEN
            ibin = CEILING((ene-minene)/de)
            espectra(ibin, j) = espectra(ibin, j) + 1
        ENDIF
    ENDDO
END SUBROUTINE EnergySpectra

!******************************************************************************
! Diagnostics for energy spectra, including the total magnetic energy,
! electric energy, particle kinetic energy and particle energy spectra.
!******************************************************************************
SUBROUTINE EnergyInfo
    USE parameters
    USE hdf5
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    CHARACTER(LEN=15) :: fname1
    CHARACTER(LEN=15) :: fname2
    CHARACTER(LEN=5) :: gname
    CHARACTER(LEN=15) :: dsetname1
    CHARACTER(LEN=15) :: dsetname2

    INTEGER(HID_T) :: file_id         ! File identifier
    INTEGER(HID_T) :: group_id         ! Group identifier
    INTEGER(HID_T) :: dset_tote        ! Total energy
    INTEGER(HID_T) :: dset_spect       ! Particle energy spectra
    INTEGER(HID_T) :: dset_ebins       ! Energy spectra bins
    INTEGER(HID_T) :: filespace        ! Dataspace identifier in file
!    INTEGER(HID_T) :: plist_id         ! Property list identifier
    INTEGER :: rank = 2                ! Datasets rank
    INTEGER :: i

    INTEGER(HSIZE_T), DIMENSION(2) :: dims1 = (/2,nbins/)
    INTEGER(HSIZE_T), DIMENSION(2) :: dims2 = (/4,ntime/)

    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: enebins

    INTEGER :: error, comm
    comm = MPI_COMM_WORLD
!!    info = MPI_INFO_NULL
!    CALL MPI_INFO_CREATE(info, ierr)
!  ! Disable ROMIO's data-sieving
!    CALL MPI_INFO_SET(info, "romio_ds_read", "disable", ierr)
!    CALL MPI_INFO_SET(info, "romio_ds_write", "disable", ierr)
!  ! Enable ROMIO's collective buffering
!    CALL MPI_INFO_SET(info, "romio_cb_read", "enable", ierr)
!    CALL MPI_INFO_SET(info, "romio_cb_write", "enable", ierr)


    WRITE(fname1, "(A10)") "spectra.h5"
    WRITE(fname2, "(A13)") "energy_tot.h5"
    WRITE(gname, "(I5.5)") itime
    dsetname1 = 'spectra'
    dsetname2 = "TotE"

    ALLOCATE(espectra(2,nbins))
    ALLOCATE(espectra_tot(2,nbins))
    ALLOCATE(totene2(4,0:ntime))
    IF (taskid .EQ. 0) THEN
        ALLOCATE(enebins(2,nbins))
        DO i = 1, nbins
            enebins(1,i) = i*dee
            enebins(2,i) = i*dei
        ENDDO
    ENDIF
    CALL EnergySpectra
    CALL MPI_REDUCE(espectra, espectra_tot, 2*nbins, &
        MPI_INTEGER, MPI_SUM, 0, comm, ierr)
    CALL MPI_REDUCE(totene, totene2, 4*(ntime+1), &
        MPI_REAL, MPI_SUM, 0, comm, ierr)

    IF (taskid .EQ. 0) THEN
      ! Initialize FORTRAN interface
        CALL h5open_f(error)
      ! Create or open a file using default properties
        IF (itime .EQ. 0) THEN
            CALL h5fcreate_f(fname1, H5F_ACC_TRUNC_F, file_id, error)
        ELSE
            CALL h5fopen_f(fname1, H5F_ACC_RDWR_F, file_id, error)
        ENDIF

        CALL h5gcreate_f(file_id, gname, group_id, error)
      ! Create the data space for the datasets
        CALL h5screate_simple_f(rank, dims1, filespace, error)
      ! Create the dataset
        CALL h5dcreate_f(group_id, dsetname1, H5T_NATIVE_INTEGER, &
            filespace, dset_spect, error)
      ! Write the dataset to the group
        CALL h5dwrite_f(dset_spect, H5T_NATIVE_INTEGER, &
            espectra_tot, dims1, error)
      ! Close the dataset
        CALL h5dclose_f(dset_spect, error)
      ! Close the data space for the dataset
        CALL h5sclose_f(filespace, error)
      ! Close the groups
        CALL h5gclose_f(group_id, error)

        IF (itime .EQ. 0) THEN
            CALL h5screate_simple_f(rank, dims1, filespace, error)
            CALL h5dcreate_f(file_id, "eBins", H5T_NATIVE_REAL, &
                filespace, dset_ebins, error)
            CALL h5dwrite_f(dset_ebins, H5T_NATIVE_REAL, &
                enebins, dims1, error)
            CALL h5dclose_f(dset_ebins, error)
            CALL h5sclose_f(filespace, error)
        ENDIF
      ! Close the files
        CALL h5fclose_f(file_id, error)

        DEALLOCATE(enebins)

        CALL h5screate_simple_f(rank, dims2, filespace, error)
        IF (itime .EQ. 0) THEN
            CALL h5fcreate_f(fname2, H5F_ACC_TRUNC_F, file_id, error)
            CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_REAL, &
                filespace, dset_tote, error)
        ELSE
            CALL h5fopen_f(fname2, H5F_ACC_RDWR_F, file_id, error)
            CALL h5dopen_f(file_id, dsetname2, dset_tote, error)
        ENDIF
        CALL h5dwrite_f(dset_tote, H5T_NATIVE_REAL, totene2, dims2, error)

        CALL h5dclose_f(dset_tote, error)
        CALL h5sclose_f(filespace, error)
        CALL h5fclose_f(file_id, error)
        CALL h5close_f(error)
    ENDIF
    DEALLOCATE(espectra_tot)
    DEALLOCATE(espectra)
    DEALLOCATE(totene2)
END SUBROUTINE EnergyInfo

!******************************************************************************
! Diagnostics for phase space information of fields fluctuation
!******************************************************************************
SUBROUTINE DispersionInfo
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: comm
    INTEGER :: i, j, k, itdiag
    comm = MPI_COMM_WORLD

    itdiag = MOD(itime-1, ntdiag) + 1
    DO i = 1, ndiagy
        fieldsy(i,1:my,itdiag) = emfGrid(diagIndexy(i),1:my)
    ENDDO
    j = ndiagxmpi(taskid) - ndiagxmpi(taskid-1)
    DO i = 1, j
        k = ndiagxmpi(taskid-1)+i
        k = diagIndexx(k)
        fieldsx(i,:,itdiag) = emfGrid(1:mx,k)
    ENDDO
    IF (MOD(itime, ntdiag) .EQ. 0) THEN
        CALL FieldsDiag
    ENDIF
END SUBROUTINE DispersionInfo

!******************************************************************************
! FFT transform of fields components
!******************************************************************************
SUBROUTINE FFTFields
    USE parameters
    USE, intrinsic :: iso_c_binding
    USE hdf5
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INCLUDE 'fftw3-mpi.f03'
    CHARACTER(LEN=15) :: fname
    CHARACTER(LEN=15), DIMENSION(2) :: gnamelist
    CHARACTER(LEN=5), DIMENSION(6) :: dnamelist, dnamelist1
    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: group_id
    INTEGER(HID_T) :: dset_id
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: plist_id
!    INTEGER(HID_T) :: crp_list    ! Dataset creation property
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: fields_tot
    REAL(kind=doublep), DIMENSION(:,:,:), ALLOCATABLE :: fieldsfft

    INTEGER(HSIZE_T), DIMENSION(3) :: dimsr, maxdimsr
    INTEGER(HSIZE_T), DIMENSION(3) :: offset1, count1
    INTEGER(HSIZE_T), DIMENSION(3) :: dimsc

    INTEGER :: rank = 3
    INTEGER :: error, comm, info, ndiag
    INTEGER :: i, j, k

    INTEGER(C_INTPTR_T) :: L, M
    TYPE(C_PTR) :: plan, cdata
    COMPLEX(C_DOUBLE_COMPLEX), pointer :: data(:,:)
    INTEGER(C_INTPTR_T) :: i1, alloc_local, local_M, local_j_offset

    comm = MPI_COMM_WORLD
!    info = MPI_INFO_NULL
    CALL MPI_INFO_CREATE(info, ierr)
  ! Disable ROMIO's data-sieving
    CALL MPI_INFO_SET(info, "romio_ds_read", "disable", ierr)
    CALL MPI_INFO_SET(info, "romio_ds_write", "disable", ierr)
  ! Enable ROMIO's collective buffering
    CALL MPI_INFO_SET(info, "romio_cb_read", "enable", ierr)
    CALL MPI_INFO_SET(info, "romio_cb_write", "enable", ierr)


    fname = "fields.h5"
    gnamelist = (/'FieldsAlongx', 'FieldsAlongy'/)
    dnamelist = (/'Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz'/)
    dnamelist1 = (/'FFTEx', 'FFTEy', 'FFTEz', 'FFTBx', 'FFTBy', 'FFTBz'/)

    CALL fftw_mpi_init
    CALL h5open_f(error)
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
    CALL MPI_INFO_FREE(info, ierr)

    CALL h5fopen_f (fname, H5F_ACC_RDWR_F, file_id, &
        error, access_prp=plist_id)
    CALL h5pclose_f(plist_id, error)

  ! Create property list for independent dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    L = ntime
    dimsc = (/1,1,1/)
    DO i = 1, 2
        IF(i .EQ. 1) THEN
            M = nx
            local_M = M/numtasks
            ndiag = ndiagx
        ELSE
            M = ny
            local_M = my
            ndiag = ndiagy
        ENDIF
        local_j_offset = taskid*local_M
        alloc_local = fftw_mpi_local_size_2d(M, L, comm, local_M, local_j_offset)
        cdata = fftw_alloc_complex(alloc_local)
        CALL c_f_pointer(cdata, data, [L,local_M])
      ! Create MPI plan for in-place forward DFT (note dimension reversal)
        plan = fftw_mpi_plan_dft_2d(M, L, data, data, comm, FFTW_FORWARD, FFTW_MEASURE)
        CALL h5gopen_f(file_id, gnamelist(i), group_id, error)
        ALLOCATE(fields_tot(ndiag,local_M,L))
        ALLOCATE(fieldsFFT(ndiag,local_M,L))
        DO j = 1, 6
            CALL h5dopen_f(group_id, dnamelist(j), dset_id, error)
            CALL h5dget_space_f(dset_id, filespace, error)
            CALL h5sget_simple_extent_dims_f(filespace, dimsr, maxdimsr, error)

            count1(1) = ndiag
            count1(2) = local_M
            count1(3) = L
            offset1(1) = 0
            offset1(2) = local_j_offset
            offset1(3) = 0

            CALL h5screate_simple_f(rank, count1, memspace, error)
            CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset1, &
                count1, error)
            CALL h5dread_f(dset_id, H5T_NATIVE_REAL, fields_tot, &
                dimsr, error, file_space_id=filespace, mem_space_id=memspace, &
                xfer_prp=plist_id)
            CALL h5sclose_f(memspace, error)
            CALL h5dclose_f(dset_id, error)
            CALL h5sclose_f(filespace, error)

            DO k = 1, ndiag
                DO i1 = 1, local_M
                    data(:,i1) = DBLE(fields_tot(k,i1,:))
                ENDDO
                CALL fftw_mpi_execute_dft(plan, data, data)
                DO i1 = 1, local_M
                    fieldsFFT(k,i1,:) = ABS(data(:,i1))
                ENDDO
            ENDDO

!            maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
!            CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
!            CALL h5pset_chunk_f(crp_list, rank, dimsc, error)

            CALL h5screate_simple_f(rank, dimsr, filespace, error)
            CALL h5dcreate_f(group_id, dnamelist1(j), H5T_NATIVE_DOUBLE, &
                filespace, dset_id, error)
            CALL h5screate_simple_f(rank, count1, memspace, error)
            CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset1, &
                count1, error)
            CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fieldsFFT, dimsr, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            CALL h5sclose_f(memspace, error)

            CALL h5dclose_f(dset_id, error)
            CALL h5sclose_f(filespace, error)
!            CALL h5pclose_f(crp_list, error)
        ENDDO
        DEALLOCATE(fields_tot, fieldsFFT)
        CALL h5gclose_f(group_id, error)
        CALL fftw_destroy_plan(plan)
        CALL fftw_free(cdata)
    ENDDO
    CALL h5pclose_f(plist_id, error)
    CALL h5fclose_f(file_id,error)
    CALL h5close_f(error)
END SUBROUTINE FFTFields
!******************************************************************************
! Fields diagnostics for FFT transform
!******************************************************************************
SUBROUTINE FieldsDiag
    USE parameters
    USE hdf5
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    CHARACTER(LEN=2), DIMENSION(6) :: dnamelist
    CHARACTER(LEN=15) :: fname
    CHARACTER(LEN=15) :: groupname
    INTEGER(HID_T) :: file_id
    INTEGER(HID_T) :: group_id
    INTEGER(HID_T) :: dset_ex, dset_ey, dset_ez
    INTEGER(HID_T) :: dset_bx, dset_by, dset_bz
    INTEGER(HID_T) :: filespace
    INTEGER(HID_T) :: memspace
    INTEGER(HID_T) :: plist_id
!    INTEGER(HID_T) :: crp_list    ! Dataset creation property

    INTEGER(HSIZE_T), DIMENSION(3) :: dimsx, dimsy
!    INTEGER(HSIZE_T), DIMENSION(3) :: maxdims
    INTEGER(HSIZE_T), DIMENSION(3) :: dimsr, maxdimsr
    INTEGER(HSIZE_T), DIMENSION(3) :: size1, offset1, count1

    INTEGER :: rank = 3
    INTEGER :: j, error, comm, info
    comm = MPI_COMM_WORLD
!    info = MPI_INFO_NULL
    CALL MPI_INFO_CREATE(info, ierr)
  ! Disable ROMIO's data-sieving
    CALL MPI_INFO_SET(info, "romio_ds_read", "disable", ierr)
    CALL MPI_INFO_SET(info, "romio_ds_write", "disable", ierr)
  ! Enable ROMIO's collective buffering
    CALL MPI_INFO_SET(info, "romio_cb_read", "enable", ierr)
    CALL MPI_INFO_SET(info, "romio_cb_write", "enable", ierr)

    dimsx = (/ndiagxmpi(numtasks-1),nx,ntime/)
    dimsy = (/ndiagy,ny,ntime/)

    dnamelist = (/'Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz'/)
    WRITE(fname, '(A9)') 'fields.h5'

    CALL h5open_f(error)
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
    CALL MPI_INFO_FREE(info, ierr)

    IF (itime .EQ. ntdiag) THEN
        CALL h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, &
            error, access_prp=plist_id)
    ELSE
        CALL h5fopen_f (fname, H5F_ACC_RDWR_F, file_id, &
            error, access_prp=plist_id)
    ENDIF
    CALL h5pclose_f(plist_id, error)

    groupname = "FieldsAlongx"
    IF (itime .EQ. ntdiag) THEN
        CALL h5gcreate_f(file_id, groupname, group_id, error)
        CALL h5screate_simple_f(rank, dimsx, filespace, error)
        CALL h5dcreate_f(group_id, dnamelist(1), H5T_NATIVE_REAL, &
            filespace, dset_ex, error)
        CALL h5dcreate_f(group_id, dnamelist(2), H5T_NATIVE_REAL, &
            filespace, dset_ey, error)
        CALL h5dcreate_f(group_id, dnamelist(3), H5T_NATIVE_REAL, &
            filespace, dset_ez, error)
        CALL h5dcreate_f(group_id, dnamelist(4), H5T_NATIVE_REAL, &
            filespace, dset_bx, error)
        CALL h5dcreate_f(group_id, dnamelist(5), H5T_NATIVE_REAL, &
            filespace, dset_by, error)
        CALL h5dcreate_f(group_id, dnamelist(6), H5T_NATIVE_REAL, &
            filespace, dset_bz, error)
    ELSE
        CALL h5gopen_f(file_id, groupname, group_id, error)
        CALL h5dopen_f(group_id, dnamelist(1), dset_ex, error)
        CALL h5dopen_f(group_id, dnamelist(2), dset_ey, error)
        CALL h5dopen_f(group_id, dnamelist(3), dset_ez, error)
        CALL h5dopen_f(group_id, dnamelist(4), dset_bx, error)
        CALL h5dopen_f(group_id, dnamelist(5), dset_by, error)
        CALL h5dopen_f(group_id, dnamelist(6), dset_bz, error)
        CALL h5dget_space_f(dset_ex, filespace, error)
    ENDIF

  ! Create property list for independent dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    j = ndiagxmpi(taskid) - ndiagxmpi(taskid-1)
    count1(1) = j
    count1(2) = nx
    count1(3) = ntdiag
    offset1(1) = ndiagxmpi(taskid-1)
    offset1(2) = 0
    offset1(3) = itime-ntdiag

    CALL h5screate_simple_f(rank, count1, memspace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset1, &
        count1, error)

    CALL h5dwrite_f(dset_ex, H5T_NATIVE_REAL, fieldsx%ex, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    CALL h5dwrite_f(dset_ey, H5T_NATIVE_REAL, fieldsx%ey, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    CALL h5dwrite_f(dset_ez, H5T_NATIVE_REAL, fieldsx%ez, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    CALL h5dwrite_f(dset_bx, H5T_NATIVE_REAL, fieldsx%bx, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    CALL h5dwrite_f(dset_by, H5T_NATIVE_REAL, fieldsx%by, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    CALL h5dwrite_f(dset_bz, H5T_NATIVE_REAL, fieldsx%bz, dimsx, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)

    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dset_ex, error)
    CALL h5dclose_f(dset_ey, error)
    CALL h5dclose_f(dset_ez, error)
    CALL h5dclose_f(dset_bx, error)
    CALL h5dclose_f(dset_by, error)
    CALL h5dclose_f(dset_bz, error)
    CALL h5sclose_f(filespace, error)
    CALL h5gclose_f(group_id, error)
    CALL h5pclose_f(plist_id, error)

  ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    groupname = "FieldsAlongy"

    count1(1) = ndiagy
    count1(2) = my
    count1(3) = ntdiag
    offset1(1) = 0
    offset1(2) = my*taskid
    offset1(3) = itime-ntdiag

    IF (itime .EQ. ntdiag) THEN
        CALL h5gcreate_f(file_id, groupname, group_id, error)
        CALL h5screate_simple_f(rank, dimsy, filespace, error)
        CALL h5dcreate_f(group_id, dnamelist(1), H5T_NATIVE_REAL, &
            filespace, dset_ex, error)
        CALL h5dcreate_f(group_id, dnamelist(2), H5T_NATIVE_REAL, &
            filespace, dset_ey, error)
        CALL h5dcreate_f(group_id, dnamelist(3), H5T_NATIVE_REAL, &
            filespace, dset_ez, error)
        CALL h5dcreate_f(group_id, dnamelist(4), H5T_NATIVE_REAL, &
            filespace, dset_bx, error)
        CALL h5dcreate_f(group_id, dnamelist(5), H5T_NATIVE_REAL, &
            filespace, dset_by, error)
        CALL h5dcreate_f(group_id, dnamelist(6), H5T_NATIVE_REAL, &
            filespace, dset_bz, error)
    ELSE
        CALL h5gopen_f(file_id, groupname, group_id, error)
        CALL h5dopen_f(group_id, dnamelist(1), dset_ex, error)
        CALL h5dopen_f(group_id, dnamelist(2), dset_ey, error)
        CALL h5dopen_f(group_id, dnamelist(3), dset_ez, error)
        CALL h5dopen_f(group_id, dnamelist(4), dset_bx, error)
        CALL h5dopen_f(group_id, dnamelist(5), dset_by, error)
        CALL h5dopen_f(group_id, dnamelist(6), dset_bz, error)
        CALL h5dget_space_f(dset_ex, filespace, error)
    ENDIF
    CALL h5screate_simple_f(rank, count1, memspace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset1, &
        count1, error)

    CALL h5dwrite_f(dset_ex, H5T_NATIVE_REAL, fieldsy%ex, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    CALL h5dwrite_f(dset_ey, H5T_NATIVE_REAL, fieldsy%ey, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    CALL h5dwrite_f(dset_ez, H5T_NATIVE_REAL, fieldsy%ez, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    CALL h5dwrite_f(dset_bx, H5T_NATIVE_REAL, fieldsy%bx, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    CALL h5dwrite_f(dset_by, H5T_NATIVE_REAL, fieldsy%by, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)
    CALL h5dwrite_f(dset_bz, H5T_NATIVE_REAL, fieldsy%bz, dimsy, error, &
        file_space_id=filespace, mem_space_id=memspace,xfer_prp=plist_id)

    CALL h5dclose_f(dset_ex, error)
    CALL h5dclose_f(dset_ey, error)
    CALL h5dclose_f(dset_ez, error)
    CALL h5dclose_f(dset_bx, error)
    CALL h5dclose_f(dset_by, error)
    CALL h5dclose_f(dset_bz, error)
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5gclose_f(group_id, error)
    CALL h5pclose_f(plist_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
END SUBROUTINE FieldsDiag
!******************************************************************************
! Initialize the random seed with a varying seed in order to ensure a different
! random number sequence for each invocation of the program.
! Source: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!******************************************************************************
SUBROUTINE init_random_seed()
    IMPLICIT NONE
    INTEGER, ALLOCATABLE :: seed(:)
    INTEGER :: i, n, un, istat, dt(8), pid, t(2), s
    INTEGER(8) :: count, tms

    CALL RANDOM_SEED(size=n)
    ALLOCATE(seed(n))
    ! First try if the OS provides a random number generator
    OPEN(unit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    IF (istat == 0) THEN
        READ(un) seed
        CLOSE(un)
    ELSE
        ! Fallback to XOR:ing the current time and pid. The PID is
        ! useful in case one launches multiple instances of the same
        ! program in parallel.
        CALL SYSTEM_CLOCK(count)
        IF (count /= 0) THEN
            t = transfer(count, t)
        ELSE
            CALL DATE_AND_TIME(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
            t = TRANSFER(tms, t)
        END IF
        s = IEOR(t(1), t(2))
        pid = getpid() + 1099279 ! Add a prime
        s = IEOR(s, pid)
        IF (n >= 3) THEN
            seed(1) = t(1) + 36269
            seed(2) = t(2) + 72551
            seed(3) = pid
            IF (n > 3) THEN
                seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            END IF
        ELSE
            seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        END IF
    END IF
    CALL RANDOM_SEED(put=seed)
END SUBROUTINE init_random_seed

!******************************************************************************
! Generating 2 independent Gaussian random numbers using Box-Muller method
! Source: section 7.3.4 of Numerical Recipes 2007
!******************************************************************************
SUBROUTINE GaussianRand(ran1, ran2)
    REAL, INTENT(OUT) :: ran1, ran2
    REAL :: v1, v2, rsq, ran, fac
    rsq = 2.0
    DO WHILE ((rsq .GE. 1.0) .OR. (rsq .EQ. 0.0))
        CALL RANDOM_NUMBER(ran)
        v1 = 2.0*ran-1.0
        CALL RANDOM_NUMBER(ran)
        v2 = 2.0*ran-1.0
        rsq = v1*v1 + v2*v2
    END DO
    fac = SQRT(-2.0*LOG(rsq)/rsq)
    ran1 = v1*fac
    ran2 = v2*fac
END SUBROUTINE GaussianRand

!******************************************************************************
! Free used memory
!******************************************************************************
SUBROUTINE ReleaseMemory 
    USE parameters
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER :: j
    DEALLOCATE(emf, emfGrid)
    DEALLOCATE(denjr)
    !DEALLOCATE(rhoq)
    DEALLOCATE(bufsendt, bufsendb)
    DEALLOCATE(bufrecvt, bufrecvb)
    DEALLOCATE(densendt, densendb)
    DEALLOCATE(denrecvt, denrecvb)
    DEALLOCATE(ptl)
    DEALLOCATE(numdens)
    DEALLOCATE(ptlsendt, ptlsendb)
    DEALLOCATE(ptlrecvt, ptlrecvb)
    DEALLOCATE(ptlSort, ptlAlloc)
    DEALLOCATE(CrossIndex)
    DEALLOCATE(totene)

    j = ndiagxmpi(taskid) - ndiagxmpi(taskid-1)
    DEALLOCATE(fieldsx)
    DEALLOCATE(fieldsy)
    DEALLOCATE(ndiagxmpi)
    CALL MPI_TYPE_FREE(particletype, ierr)
END SUBROUTINE ReleaseMemory 
