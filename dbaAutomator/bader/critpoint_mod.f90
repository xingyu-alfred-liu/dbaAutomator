  MODULE critpoints_mod
    USE kind_mod
    USE matrix_mod
    USE bader_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    USE ions_mod
    USE weight_mod
    USE dsyevj3_mod
    IMPLICIT NONE

    PRIVATE 
    PUBLIC :: critpoint_find

    TYPE hessian

      ! du dv dw are derivatives of the three original lattice vectors read from CHGCAR
      REAL(q2),DIMENSION(3) ::  du, dv, dw
      REAL(q2) :: dudu, dvdv, dwdw, dudv, dudw, dvdw
      ! eigval and eigvec are eigenvalues and eigenvectors of hessian matrix
    END TYPE
    CONTAINS

!-----------------------------------------------------------------------------------!
!critpoint_find: find critical points on the boundary of the Bader volumes
!NOTE: this subroutine should be called after refine_edge
!      in order to restrict the calculation to edge points
!-----------------------------------------------------------------------------------!
  SUBROUTINE critpoint_find(bdr,chg,opts,ions)


! These are for screening CP
    TYPE cpc ! stands for critical point candidate
      INTEGER,DIMENSION(3) :: ind, colatint ! these are the indices of the cp
      REAL(q2),DIMENSION(3) :: force
      REAL(q2),DIMENSION(3) :: tempcart, tempind
      REAL(q2),DIMENSION(3,3,3) :: dx, dy, dz ! first derivatives of neighbors
!      REAL(q2),DIMENSION(3,3,3) :: du, dv, dw
      REAL(q2),DIMENSION(3) :: du, dv, dw ! 1 2 3 are backward, present, forward
      REAL(q2),DIMENSION(3) :: eigvals, r, cocart, colat, tempr
      INTEGER(q2),DIMENSION(8,3) :: nnind ! indices of neighbors 
      ! indices from 0 to 8 are zyx 000 001 010 011 100 101 110 111
!      REAL(q2),DIMENSION(8,3) :: nngrad ! gradients of nn mentioned above.
!      REAL(q2),DIMENSION(6,3) :: intnngrad ! gradients of interpolated neighbors
      ! used to find interpolated hessians.
      REAL(q2),DIMENSION(3,3) :: eigvecs
      LOGICAL :: proxy, keep
    END TYPE
    

    TYPE(hessian) :: hes
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(cpc),ALLOCATABLE,DIMENSION(:) :: cpl, cplt ! critical point list and a temporary
! copy
! for points, 1 and 2 are +1, -1
    INTEGER,DIMENSION(3) :: p, pt, ptt, ptx1, ptx2, pty1, pty2, ptz1, ptz2
    INTEGER,DIMENSION(3) :: tempind
    INTEGER :: n1, n2, n3, d1, d2, d3, cptnum, ucptnum, i, j, debugnum
    REAL(q2),DIMENSION(3) :: eigvec1, eigvec2, eigvec3, tempVec
    REAL(q2),DIMENSION(3) :: tempforce 
    ! to be used in newton method in finding unique critical points.
    REAL(q2),DIMENSION(3,3) ::  hessianMatrix, bkhessianMatrix
    ! these are vectors orthogonal to eigenvectors
    REAL(q2),DIMENSION(3) :: tem, tem2a,tem2b,tem2c, force,eigvals,carts
    REAL(q2) :: umag, vmag, wmag, threshhold, minmag, tempreal
    REAL(q2),DIMENSION(3,3) :: eigvecs, A, inverseHessian
    ! linearized approximated derivatives for proxy critical screening
    REAL(q2),DIMENSION(3,3) :: transformationmatrix ! normalized lattice vectors
    LOGICAL :: trilinear ! determines if hessian calculation uses interpolation
    REAL(q2) :: dx0,dx1,dy0,dy1,dz0,dz1 ! outputs from interpolated gradients
    REAL(q2),DIMENSION(6,3) :: intcarts ! positions in cartesian coordinates of 6 interpolated points
    ! row 1 2 are + and 1 x, then + and - y, then + and - z
    REAL(q2),DIMENSION(6,3) :: intgrads ! gradients of interpolated points
    REAL(q2),DIMENSION(6) :: intrhos ! rhos of interpolated points 
    REAL(q2),DIMENSION(6,3) :: intinds ! fraction indicies for interpolated
    REAL(q2) :: rhocur ! rho of current point
    REAL(q2) :: stepsize
    REAL(q2),DIMENSION(3) :: preal
    INTEGER,DIMENSION(8,3) :: nn ! alternative trilinear approx points
    REAL(q2),DIMENSION(8) :: vals
    ! points
    LOGICAL,DIMENSION(3) :: cartcoor ! check if axis are alone cartesian
    ! The followings are for finding unique critical points
    REAL(q2),DIMENSION(8,3) :: nngrad
    REAL(q2),DIMENSION(6,3) :: intnngrad
    REAL(q2),DIMENSION(3,3) :: temphessian
    REAL(q2),DIMENSION(3) :: nexttem
    INTEGER :: stepcount
    WRITE(*,'(A)')  'FINDING CRITICAL POINTS'
    
    PRINT *, chg%car2lat
    trilinear = .FALSE. ! do not interpolate, as it can get mesy
    ucptnum = 0
    cptnum = 0
    PRINT * , "These code requires -vac auto or -vac #"
    PRINT *, '-----------                  WARNING             -----------'
    PRINT *, ' An analysis of valence charges will not yield sensible results'
    PRINT *, ' Use the total charge density for finding CPs   '
    PRINT *, '____________________________________________________________'
    OPEN(97,FILE='CPF.dat',STATUS='REPLACE',ACTION='WRITE')
    OPEN(98,FILE='CPFU.dat',STATUS='REPLACE',ACTION='WRITE')
    debugnum = 0

!    PRINT *, rho_val(chg,36,42,35),rho_val(chg,37,42,35),rho_val(chg,38,42,35)
!    PRINT *, rho_val(chg,36,41,35),rho_val(chg,37,41,35),rho_val(chg,38,41,35)
!    PRINT *, rho_val(chg,36,40,35),rho_val(chg,37,40,35),rho_val(chg,38,40,35)
!    PRINT *, rho_val(chg,36,42,34),rho_val(chg,37,42,34),rho_val(chg,38,42,34)
!    PRINT *, rho_val(chg,36,41,34),rho_val(chg,37,41,34),rho_val(chg,38,41,34)
!    PRINT *, rho_val(chg,36,40,34),rho_val(chg,37,40,34),rho_val(chg,38,40,34)
!    PRINT *, rho_val(chg,36,42,33),rho_val(chg,37,42,33),rho_val(chg,38,42,33)
!    PRINT *, rho_val(chg,36,41,33),rho_val(chg,37,41,33),rho_val(chg,38,41,33)
!    PRINT *, rho_val(chg,36,40,33),rho_val(chg,37,40,33),rho_val(chg,38,40,33)
!    PRINT *, rho_val(chg,36,42,32),rho_val(chg,37,42,32),rho_val(chg,38,42,32)
!    PRINT *, rho_val(chg,36,41,32),rho_val(chg,37,41,32),rho_val(chg,38,41,32)
!    PRINT *, rho_val(chg,36,40,32),rho_val(chg,37,40,32),rho_val(chg,38,40,32)
    umag = SQRT(ions%lattice(1,1)**2 + ions%lattice(1,2)**2 + &
      ions%lattice(1,3)**2) / REAL(chg%npts(1),q2)
    vmag = SQRT(ions%lattice(2,1)**2 + ions%lattice(2,2)**2 + &
      ions%lattice(2,3)**2) / REAL(chg%npts(2),q2)
    wmag = SQRT(ions%lattice(3,1)**2 + ions%lattice(3,2)**2 + &
      ions%lattice(3,3)**2) / REAL(chg%npts(3),q2)
    minmag = MIN(umag,vmag,wmag)
    ! check if axis are cartesian
!    cartcoor = coorcheck(ions%lattice)
    IF ( ALL(cartcoor) ) trilinear = .FALSE.
    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
            ! check to see if this point is in the vacuum
            IF (bdr%volnum(n1,n2,n3) == bdr%bnum + 1) THEN
              debugnum = debugnum + 1
              CYCLE
            END IF
            p = (/n1,n2,n3/)
            preal = p
            trilinear = .FALSE.
            IF (trilinear) THEN
              ! get cartesian coordinates of current point
              !carts = getcart(p,chg%lat2car)
              carts = MATMUL(chg%lat2car,p)
              ! get positions of where interpolated points should be
              intcarts = getintcarts(carts,umag,vmag,wmag)
              ! get indices of the interpolated points
              intinds = getinds(chg%car2lat,intcarts)
              ! get gradients of interpolated points
              DO i = 1, 6
                ! indices 1 to 6 are + - x, + - y, + - z
                intgrads(i,:) = rho_grad(chg,intinds(i,:),intrhos(i))
!                PRINT *, 'old grads', intgrads(1,:)
                ! alternative : find closest neighbors to do trilinear
                ! approximation
!                nn = findnn(p,intcarts(i,:),chg)
!                intgrads(i,:) = nn_grad(chg,intinds(i,:),intrhos(i),nn)
!                PRINT *, 'new grads', intgrads(1,:)
              END DO
              ! get second derivatives
              ! dxdx
              hessianMatrix(1,1) = (intgrads(1,1) - intgrads(2,1)) / (2 * minmag)
              ! dydy
              hessianMatrix(2,2) = (intgrads(3,2) - intgrads(4,2)) / (2 * minmag)
              ! dzdz
              hessianMatrix(3,3) = (intgrads(5,3) - intgrads(6,3)) / (2 * minmag)
              ! dxdy
              hessianMatrix(1,2) = (intgrads(1,2) - intgrads(2,2)) / (2 * minmag)
              hessianMatrix(2,1) = hessianMatrix(1,2)
              ! dxdz
              hessianMatrix(1,3) = (intgrads(1,3) - intgrads(2,3)) / (2 * minmag)
              hessianMatrix(3,1) = hessianMatrix(1,3)
              ! dydz
              hessianMatrix(2,3) = (intgrads(3,3) - intgrads(4,3)) / (2 * minmag)
              hessianMatrix(3,2) = hessianMatrix(3,2)
              ! force
              force = rho_grad(chg,preal,rhocur)
!              DO i = 1 , 6
!                PRINT *, intinds(i,:)
!              END DO
!              PRINT *, '26 48 49'
!              PRINT *, rho_val(chg,26,48,49)
!              PRINT *, '25 26 27'
!              PRINT *,rho_val(chg,25,48,49),rho_val(chg,26,48,49), &
!                rho_val(chg,27,48,49)
!              PRINT *, '47 48 49'
!              PRINT *,rho_val(chg,26,47,49),rho_val(chg,26,48,49), &
!                rho_val(chg,26,49,49)
!              PRINT *, '48 49 50'
!              PRINT *,rho_val(chg,26,48,48),rho_val(chg,26,48,49), &
!                rho_val(chg,26,48,50)
            ELSE
  
              
  
  !-----------------------------------------------------------------------------------!
  ! after finding candidate critical edge points, now find the hessian
  !-----------------------------------------------------------------------------------!
               CALL getgradhes(p,chg,hes,force)
               hessianMatrix(1,1) = hes%dudu
               hessianMatrix(1,2) = hes%dudv
               hessianMatrix(1,3) = hes%dudw
               hessianMatrix(2,1) = hes%dudv
               hessianMatrix(2,2) = hes%dvdv
               hessianMatrix(2,3) = hes%dvdw
               hessianMatrix(3,1) = hes%dudw
               hessianMatrix(3,2) = hes%dvdw
               hessianMatrix(3,3) = hes%dwdw
               hessianMatrix = MATMUL(chg%car2lat,hessianMatrix)
               hessianMatrix = MATMUL(hessianMatrix,TRANSPOSE(chg%car2lat))
               force = MATMUL(chg%car2lat,force)
               ! now everything is cartesian
             END IF
!             tem2 = tem2 * umag 
!             tem2a = (/ ions%lattice(1,1),ions%lattice(2,1), & 
!               ions%lattice(3,1) /) / chg%npts(1)
             IF (n1==22.AND.n2==26.AND.n3==51) THEN
               PRINT *, n1,n2,n3
               PRINT *, force
             ELSE IF (n1==22.AND.n2==26.AND.n3==52)THEN
               PRINT *, n1,n2,n3
               PRINT *, force
             ELSE IF (n1==22.AND.n2==27.AND.n3==51)THEN
               PRINT *, n1,n2,n3
               PRINT *, force
             ELSE IF (n1==22.AND.n2==27.AND.n3==52)THEN
               PRINT *, n1,n2,n3
               PRINT *, force
             ELSE IF (n1==23.AND.n2==26.AND.n3==51)THEN
               PRINT *, n1,n2,n3
               PRINT *, force
             ELSE IF (n1==23.AND.n2==26.AND.n3==52)THEN
               PRINT *, n1,n2,n3
               PRINT *, force
             ELSE IF (n1==23.AND.n2==27.AND.n3==51)THEN
               PRINT *, n1,n2,n3
               PRINT *, force
             ELSE IF (n1==23.AND.n2==27.AND.n3==52)THEN
               PRINT *, n1,n2,n3
               PRINT *, force
             END IF
             inverseHessian = inverse(hessianMatrix)
             tem = MATMUL(inverseHessian,force)
             ! convert from cartesian to lattice
             tem = MATMUL(chg%car2lat,tem)
             IF (ABS(tem(1)) <= 0.5) THEN
               IF (ABS(tem(2)) <= 0.5) THEN
                 IF (ABS(tem(3)) <= 0.5) THEN
                   cptnum = cptnum + 1
                   WRITE(97,*) '*********** A NEW ENTRY *************'
                   bkhessianMatrix = hessianMatrix
                   CALL DSYEVJ3(hessianMatrix,eigvecs,eigvals)
                   WRITE(97,*) 'Critical point number: ', cptnum
                   WRITE(97,*) "Indices are"
                   WRITE(97,*) p(1),p(2),p(3)
                   WRITE(97,*) "Density at this point is" 
                   WRITE(97,*) rho_val(chg,p(1),p(2),p(3))
                   WRITE(97,*) 'Eigenvalues: '
                   WRITE(97,*) eigvals
                   WRITE(97,'(3(1X,E18.11))') 
                   WRITE(97,*) 'Eigenvectors:'
                   WRITE(97,*) eigvecs(1,:)
                   WRITE(97,*) eigvecs(2,:)
                   WRITE(97,*) eigvecs(3,:)
                   WRITE(97,*) 'tem', tem
                   WRITE(97,*) 'force is'
                   WRITE(97,*)  force
                   WRITE(97,*) 'Hessian is'
                   WRITE(97,*)  bkhessianMatrix
                   WRITE(97,*) 'du and dudu is'
                   WRITE(97,*) hes%du, hes%dudu
                   IF (cptnum == 1)  THEN
                     ALLOCATE(cpl(1))
                     cpl(1)%du = hes%du
                     cpl(1)%dv = hes%dv  
                     cpl(1)%dw = hes%dw
                     cpl(1)%ind(1) = n1
                     cpl(1)%ind(2) = n2
                     cpl(1)%ind(3) = n3
                     cpl(1)%force = force
                     cpl(1)%proxy = .FALSE.
                     cpl(1)%keep = .TRUE.
                     cpl(1)%eigvals = eigvals
                     cpl(1)%eigvecs = eigvecs
                     cpl(1)%r = tem
                     cpl(1)%tempcart = MATMUL(tem + p,chg%car2lat)
                   ELSE
                     ALLOCATE(cplt(cptnum))
                     DO i = 1, cptnum -1
                       cplt(i) = cpl(i)
                     END DO
                     DEALLOCATE(cpl)
                     ALLOCATE(cpl(cptnum))
                     DO i = 1, cptnum - 1
                       cpl(i)=cplt(i)
                     END DO
                     DEALLOCATE(cplt)
                     cpl(cptnum)%du = hes%du 
                     cpl(cptnum)%dv = hes%dv
                     cpl(cptnum)%dw = hes%dw
                     cpl(cptnum)%ind(1) = n1
                     cpl(cptnum)%ind(2) = n2
                     cpl(cptnum)%ind(3) = n3
                     cpl(cptnum)%force = force
                     cpl(cptnum)%proxy = .FALSE.
                     cpl(cptnum)%keep = .TRUE.
                     cpl(cptnum)%eigvals = eigvals
                     cpl(cptnum)%eigvecs = eigvecs
                     cpl(cptnum)%r = tem
                     cpl(cptnum)%tempcart = MATMUL(tem + p, chg%car2lat)
                   END IF
                 END IF
               END IF
!               ELSE 
             END IF
        END DO
      END DO
    END DO
    PRINT *,'CRITICAL POINTS INFO WRITEN TO CPF.dat'
    PRINT *, "CRITICAL POINTS FOUND: ", cptnum 
    PRINT *, 'FINDING UNIQUE CRITICAL POINTS...'
!! ************* DONT DELETE THESE CODE **************************
!!    ! now go through candidates, for each candidate, see if there is another
!!    ! within 1 distance on u v w direction.
!!    DO i = 1 , cptnum - 1
!!      DO j = i + 1 , cptnum
!!        IF (i < j) THEN
!!          IF (ABS(cpl(i)%ind(1) - cpl(j)%ind(1)) <= 1 .AND. &
!!             ABS(cpl(i)%ind(2) - cpl(j)%ind(2)) <= 1 .AND. & 
!!             ABS(cpl(i)%ind(3) - cpl(j)%ind(3)) <= 1) THEN
!!             cpl(i)%proxy = .TRUE.
!!             cpl(j)%proxy = .TRUE.
!!             ! we know the cp is being shared. follow r to see
!!             ! if we go out of cell. 
!!             ! first step is to linearize dx
!!             ldu = cpl(i)%du(2) - cpl(j)%du(2)
!!             ldv = cpl(i)%dv(2) - cpl(j)%dv(2)
!!             ldw = cpl(i)%dw(2) - cpl(j)%dw(2) 
!!!             END IF
!!             WRITE(97,*)  i, j   
!!          END IF
!!        END IF
!!      END DO
!!    END DO
!!    PRINT *, 'number of unique cp :',debugnum
!!*******************************************************************
!!*******************************************************************
    ! To find critical points (unique), start with a cell that contains a
    ! critical point and its hessian and force. Use Newton's method to make a
    ! move. Interpolate the force inside the voxel. 
    ! Once moved, get the new force through trilinear interpolation, and
    ! get the new hessian which will be a matrix of constants, make moves until
    ! r is zero. get the coordinates of the new true critical point. If this
    ! point is within half lattice to another, do not record this new point.
    DO i = 1, cptnum
      stepcount = 0
      PRINT *, 'indices of cpt is'
      PRINT *,cpl(i)%ind
!      PRINT *, 'r is'
!      PRINT *,cpl(i)%r
      ! move to r in lattice units
      cpl(i)%colat = cpl(i)%ind + cpl(i)%r 
      ! find the nearest neighbors to this point.
      ! first find a nearby grid point.
      cpl(i)%colatint = (/ &
         FLOOR(cpl(i)%colat(1)),&
         FLOOR(cpl(i)%colat(2)),& 
         FLOOR(cpl(i)%colat(3))/)
      ! get cartesians of the current location
!      PRINT *,'colatint is', cpl(i)%colatint 
      cpl(i)%cocart = MATMUL(chg%lat2car,cpl(i)%colat)
!      PRINT *,'cocart is', cpl(i)%cocart
      ! find nearest neighbors
!      cpl(i)%nnind = findnn(cpl(i)%colatint,cpl(i)%cocart,chg)
      ! try not to find NN first. rather just assume NN's.
      ! x-1, y-1, z-1 000
      cpl(i)%nnind(1,:) = cpl(i)%ind + & 
        (/FLOOR(cpl(i)%r(1)),FLOOR(cpl(i)%r(2)),FLOOR(cpl(i)%r(3))/)
      ! x+1, y-1, z-1 100
      cpl(i)%nnind(2,:) = cpl(i)%ind + &
        (/CEILING(cpl(i)%r(1)),FLOOR(cpl(i)%r(2)),FLOOR(cpl(i)%r(3))/)
      ! x-1, y+1, z-1 010
      cpl(i)%nnind(3,:) = cpl(i)%ind + &
        (/FLOOR(cpl(i)%r(1)),CEILING(cpl(i)%r(2)),FLOOR(cpl(i)%r(3))/)
      ! x+1, y+1, z-1 110
      cpl(i)%nnind(4,:) = cpl(i)%ind + &
        (/CEILING(cpl(i)%r(1)),CEILING(cpl(i)%r(2)),FLOOR(cpl(i)%r(3))/)
      ! x-1, y-1, z+1 001
      cpl(i)%nnind(5,:) = cpl(i)%ind + &
        (/FLOOR(cpl(i)%r(1)),FLOOR(cpl(i)%r(2)),CEILING(cpl(i)%r(3))/)
      ! x+1, y-1, z+1 101
      cpl(i)%nnind(6,:) = cpl(i)%ind + &
        (/CEILING(cpl(i)%r(1)),FLOOR(cpl(i)%r(2)),CEILING(cpl(i)%r(3))/)
      ! x-1, y+1, z+1 011
      cpl(i)%nnind(7,:) = cpl(i)%ind + &
        (/FLOOR(cpl(i)%r(1)),CEILING(cpl(i)%r(2)),CEILING(cpl(i)%r(3))/)
      ! x            111
      cpl(i)%nnind(8,:) = cpl(i)%ind + &
        (/CEILING(cpl(i)%r(1)),CEILING(cpl(i)%r(2)),CEILING(cpl(i)%r(3))/)
      ! get gradients "force" of all nearest neighbors
!      PRINT *, 'nn indices are'
      DO j = 1, 8
        tempind(1) = cpl(i)%nnind(j,1)
        tempind(2) = cpl(i)%nnind(j,2)
        tempind(3) = cpl(i)%nnind(j,3)
        ! getting nn indices are checked. pass. 
        CALL getgradhes(tempind,chg,hes,nngrad(j,:))
        ! this force is not adjusted with lattice
        !cpl(i)%nngrad(j,:) = MATMUL(chg%car2lat,cpl(i)%nngrad(j,:))
        !PRINT *, cpl(i)%nngrad(j,:)
        ! now things are cartesian. forces are checked to be correct.
        ! correct: if I make things cartesian, forces are checked out to be
        ! correct. However I still need to find hessian in lattice first. 
        ! So I do not convert force to cartesian yet.
        ! force is checked to be correct.
      END DO
      ! Now start newton method iterations
      DO stepcount = 1,10
        ! The next big step is to interpolate the force at predicted critical
        ! point.
        tempforce = trilinear_interpol_grad(nngrad,cpl(i)%r) ! val r interpol
        PRINT *, 'force is'
        PRINT *, MATMUL(chg%car2lat,tempforce)
        ! assume the force is calculated correctly. 
        ! The next step is to find a reasonable step size. 
        ! cpl(i)%r is in lattice units so step size should also be in lattice
        ! units. 
        stepsize = findstepsize(cpl(i)%r)
        PRINT *, 'stepsize is'
        PRINT *, stepsize
        ! it seems to be able to find the correct step size
        ! the next step is to find the hessian at the new location.
        ! to do this, forces are needed at the 8 neighbor sites. 
        ! list of intnngrad is ordered as follows from index 1 to 6
        ! +x -x +y -y +z -z 
  !      cpl(i)%tempr = cpl(i)%r + (/stepsize,REAL(0.,q2),REAL(0.,q2)/)
        intnngrad(1,:) = trilinear_interpol_grad(nngrad, &
          cpl(i)%r + (/stepsize,REAL(0.,q2),REAL(0.,q2)/))
        intnngrad(2,:) = trilinear_interpol_grad(nngrad, &
          cpl(i)%r + (/-stepsize,REAL(0.,q2),REAL(0.,q2)/))
        intnngrad(3,:) = trilinear_interpol_grad(nngrad, &
          cpl(i)%r + (/REAL(0.,q2),stepsize,REAL(0.,q2)/))
        intnngrad(4,:) = trilinear_interpol_grad(nngrad, &
          cpl(i)%r + (/REAL(0.,q2),-stepsize,REAL(0.,q2)/))
        intnngrad(5,:) = trilinear_interpol_grad(nngrad, &
          cpl(i)%r + (/REAL(0.,q2),REAL(0.,q2),stepsize/))
        intnngrad(6,:) = trilinear_interpol_grad(nngrad, &
          cpl(i)%r + (/REAL(0.,q2),REAL(0.,q2),-stepsize/))
        ! should expect gradient of neighbors with different diretions
        ! nngrad forces are in lattice. so these forces are also in lattice. 
        temphessian = inthessian(intnngrad,stepsize) 
        PRINT *, 'hessian is'
        PRINT *, MATMUL(MATMUL(chg%car2lat,temphessian),TRANSPOSE(chg%car2lat))
        ! assuming it is correct.
        ! now get force and hessian back in cartesian
        DO j = 1, 6
          intnngrad(j,:) = MATMUL(chg%car2lat,intnngrad(j,:))
        END DO
        temphessian = MATMUL(chg%car2lat,temphessian)
        temphessian = MATMUL(temphessian,TRANSPOSE(chg%car2lat))
        nexttem = newtonstep(chg%car2lat,chg%lat2car,temphessian,tempforce)
        PRINT *, 'nexttem is'
        PRINT *, nexttem
        ! lets not worry about the nearest neighbors yet, for finding gradient 
        ! update critical point location
        cpl(i)%r = cpl(i)%r + nexttem
        PRINT *, 'new location is'
        PRINT *, cpl(i)%r 
        IF (cpl(i)%r(1)>=0.5.OR.cpl(i)%r(2)>=0.5.OR.cpl(i)%r(3)>=0.5) THEN
          PRINT *,stepcount
          EXIT
        END IF 
        STOP
      END DO
    END DO    
    DEALLOCATE (cpl)
    CLOSE(97)
    CLOSE(98)
    END SUBROUTINE critpoint_find


    REAL(q2) FUNCTION projection(r,rp)
      ! r is the vector pointing towards a CP
      ! rp is a cell vector to be projected onto
      REAL(q2),DIMENSION(3) :: r, rp
      projection = DOT_PRODUCT(r,rp)/sqrt((rp(1)**2 + &
      rp(2)**2 + rp(3)**2))
      RETURN
    END FUNCTION

    ! get cartesian coordinates of a point
    FUNCTION getcart(ind,lat2car)
      !ind is indicies of the current point
      !cart is the cartisian coordinates of the current point
      INTEGER,DIMENSION(3),INTENT(IN) :: ind
      REAL(q2),DIMENSION(3) :: getcart
      REAL(q2),DIMENSION(3,3) :: lat2car
      getcart(1) = ind(1) * lat2car(1,1) + & 
        ind(2) * lat2car(1,2) + & 
        ind(3) *  lat2car(1,3)
      getcart(2) = ind(1) * lat2car(2,1) + &
        ind(2) * lat2car(2,2) + &
        ind(3) *  lat2car(2,3)
      getcart(3) = ind(1) * lat2car(3,1) + &
        ind(2) * lat2car(3,2) + &
        ind(3) *  lat2car(3,3)

      RETURN
    END FUNCTION

    ! get cartesian coordinates of interpolated points
    FUNCTION  getintcarts(carts,umag,vmag,wmag)
      REAL(q2) :: minmag,umag,vmag,wmag
      REAL(q2),DIMENSION(3) :: carts
      REAL(q2),DIMENSION(6,3) :: getintcarts
      INTEGER :: i,j,k
      ! first find which one is the smallest mag
      minmag = MIN(umag,vmag,wmag)
      DO i = 1, 6
        getintcarts(i,:) = carts
      END DO
      getintcarts(1,1) = getintcarts(1,1) + minmag
      getintcarts(2,1) = getintcarts(2,1) - minmag
      getintcarts(3,2) = getintcarts(3,2) + minmag
      getintcarts(4,2) = getintcarts(4,2) - minmag
      getintcarts(5,3) = getintcarts(5,3) + minmag
      getintcarts(6,3) = getintcarts(6,3) - minmag
      RETURN
    END FUNCTION
    
    ! get indices from cartesian coordinates
    FUNCTION  getinds(car2lat,intcarts) 
      REAL(q2),DIMENSION(6,3) :: getinds, intcarts
      REAL(q2),DIMENSION(3,3) :: car2lat,W
      REAL(q2),DIMENSION(3) :: Z
      INTEGER :: i,j,k
      DO i = 1, 6 
        getinds(i,:) = MATMUL(car2lat,intcarts(i,:))
      END DO
      RETURN
    END FUNCTION

    ! finds nearest grid point to a interpolated point
    ! used for doing trilinear interpolation
    ! p is the centered on grid point. intcart is a nearby point to be
    ! interpolated in cartesian coordinates.
    FUNCTION findnn(p,intcart,chg) ! THIS FUNCTION IS NOT STABLE
      TYPE(weight_obj),DIMENSION(27) :: nndist ! record distances of a point to all nearest on
      !grid points
      TYPE(charge_obj) :: chg
      INTEGER,DIMENSION(8,3) :: findnn,tfindnn
      REAL(q2),DIMENSION(3) :: intcart,p2cart
      INTEGER,DIMENSION(3) :: p,maxs
      REAL(q2),DIMENSION(27) :: dist
      INTEGER,DIMENSION(27) :: scores
      INTEGER,DIMENSION(27,3) :: p2
      INTEGER :: i,j,k,counter
      ! since only closest neigherbos of current grid point
      ! will be used for interpolation
      counter = 1
      !CALL pbc(p,chg%npts)
      !to calculate distance it is not necessary to run pbc
      !infact pbc should be avoided at this stage
      DO i = -1,1
        DO j = -1,1
          DO k = -1,1
            p2(counter,:) = (/i,j,k/) + p
            !CALL pbc(p2(counter,:),chg%npts)
            PRINT *,'p2 is', p2(counter,:)
            !p2cart = getcart(p2,chg%lat2car)
            p2cart = MATMUL(p2(counter,:),chg%lat2car)
!            PRINT *, 'p2 cart is', p2cart
!            PRINT *,'counter is', counter
            nndist(counter)%rho = SQRT( &
              (intcart(1) - p2cart(1))**2 + & 
              (intcart(2) - p2cart(2))**2 + &
              (intcart(3) - p2cart(3))**2 &
              )
            PRINT *,'dist is',dist(counter)
            counter = counter + 1
          END DO
        END DO
      END DO
      ! now find the top 4 smallest array
      ! compare each element to all
      ! if it is smaller, it gets score
      ! keep the ones with highest scores
      ! there should not be points with equal scores
      ! Each point should have a score ranging from 0 to 26
      DO i = 1, 27
        scores(i) = 0
        DO j = 1, 27
          IF (i == j) THEN
            CYCLE
          END IF
          IF (dist(i) <= dist(j)) THEN
            scores(i) = scores(i) + 1
          END IF
        END DO
        IF (scores(i) == 26) THEN
          tfindnn(1,:) = p2(i,:)
          CYCLE
        END IF
        IF (scores(i) == 25) THEN
          tfindnn(2,:) = p2(i,:)
          CYCLE
        END IF
        IF (scores(i) == 24) THEN
          tfindnn(3,:) = p2(i,:)
          CYCLE
        END IF
        IF (scores(i) == 23) THEN
          tfindnn(4,:) = p2(i,:)
          CYCLE
        END IF
        IF (scores(i) == 22) THEN
          tfindnn(5,:) = p2(i,:)
          CYCLE
        END IF
        IF (scores(i) == 21) THEN
          tfindnn(6,:) = p2(i,:)
          CYCLE
        END IF
        IF (scores(i) == 20) THEN
          tfindnn(7,:) = p2(i,:)
          CYCLE
        END IF
        IF (scores(i) == 19) THEN
          tfindnn(8,:) = p2(i,:)
          CYCLE
        END IF
      END DO     
      PRINT *,'scores', scores
      PRINT *, 'tfindnn'
      DO i = 1, 8
        PRINT *, tfindnn(i,:)
      END DO
      ! then these needs to be rearranged 
      ! so that is follows this format:
      !  1   2   3   4   5   6   7   8
      ! 000 001 010 100 011 101 110 111
      ! indices should have only 2 values
      ! on each axis. min gives 0 on that 
      ! axis and max gives 1. 
      maxs = (/0,0,0/)
      DO i = 1, 8
        DO j = 1,3
          IF (tfindnn(i,j)>= maxs(j) ) THEN
            maxs(j) = tfindnn(i,j)
          END IF
        END DO
      END DO
      PRINT *, 'maxs is', maxs
      DO i = 1, 8
        PRINT *, 'tfindnn here is'
        PRINT * ,tfindnn(i,:)
        IF (ALL(tfindnn(i,:) == maxs)) THEN
          PRINT *, '8'
          findnn(8,:) = tfindnn(i,:)
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,1,1/))) THEN
          findnn(1,:) = tfindnn(i,:)
          PRINT *, '1'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,1,0/))) THEN
          findnn(2,:) = tfindnn(i,:)
          PRINT *, '2'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,0,1/))) THEN
          findnn(3,:) = tfindnn(i,:)
          PRINT *, '3'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/0,1,1/))) THEN
          findnn(4,:) = tfindnn(i,:)
          PRINT *, '4'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/1,0,0/))) THEN
          findnn(5,:) = tfindnn(i,:)
          PRINT *, '5'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/0,1,0/))) THEN
          findnn(6,:) = tfindnn(i,:)
          PRINT *, '6'
          CYCLE
        END IF
        IF (ALL(tfindnn(i,:) == maxs - (/0,0,1/))) THEN
          findnn(7,:) = tfindnn(i,:)
          PRINT *, '7'
          CYCLE
        END IF
      END DO 
      RETURN
    END FUNCTION

    ! This function takes in a list of nearest neighbors
    FUNCTION nn_grad(chg,r,rho,nn)
      TYPE(charge_obj) :: chg
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),INTENT(OUT) :: rho
      REAL(q2),DIMENSION(3) :: nn_grad
      INTEGER :: p1, p2, p3
      REAL(q2),DIMENSION(3) :: rho_grad_lat
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
      REAL(q2) :: rho00_, rho01_, rho10_, rho11_
      REAL(q2) :: rho0__, rho1__, rho_0_, rho_1_, rho__0, rho__1
      REAL(q2) :: rho_00, rho_01, rho_10, rho_11
      INTEGER,DIMENSION(8,3) :: nn
      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
      rho000 = rho_val(chg,nn(1,1),nn(1,2),nn(1,3))
      rho001 = rho_val(chg,nn(2,1),nn(2,2),nn(2,3))
      rho010 = rho_val(chg,nn(3,1),nn(3,2),nn(3,3))
      rho100 = rho_val(chg,nn(4,1),nn(4,2),nn(4,3))
      rho011 = rho_val(chg,nn(5,1),nn(5,2),nn(5,3))
      rho101 = rho_val(chg,nn(6,1),nn(6,2),nn(6,3))
      rho110 = rho_val(chg,nn(7,1),nn(7,2),nn(7,3))
      rho111 = rho_val(chg,nn(8,1),nn(8,2),nn(8,3))
      rho00_ = rho000*g3 + rho001*f3
      rho01_ = rho010*g3 + rho011*f3
      rho10_ = rho100*g3 + rho101*f3
      rho11_ = rho110*g3 + rho111*f3
      rho0__ = rho00_*g2 + rho01_*f2
      rho1__ = rho10_*g2 + rho11_*f2
      rho = rho0__*g1 + rho1__*f1
  ! More work for gradients
      rho_0_ = rho00_*g1 + rho10_*f1
      rho_1_ = rho01_*g1 + rho11_*f1
      rho_00 = rho000*g1 + rho100*f1
      rho_01 = rho001*g1 + rho101*f1
      rho_10 = rho010*g1 + rho110*f1
      rho_11 = rho011*g1 + rho111*f1
      rho__0 = rho_00*g2 + rho_10*f2
      rho__1 = rho_01*g2 + rho_11*f2
      rho_grad_lat(1) = rho1__ - rho0__
      rho_grad_lat(2) = rho_1_ - rho_0_
      rho_grad_lat(3) = rho__1 - rho__0
  !   CALL vector_matrix(rho_grad_lat, chg%car2lat, rho_grad)
      nn_grad = MATMUL(rho_grad_lat, chg%car2lat)
    RETURN
    END FUNCTION nn_grad

    ! this funciton takes in 8 values, return a
    ! trilinear interpolated gradient of the values.
    ! the 8 value list value order is 
    ! 000 001 010 100 011 101 110 111
    ! Note 02042019: the above order is what I wrote previously
    ! Note 02042019: I believe the actuall order is the following
    ! 000 100 010 110 001 101 011 111
    ! r is the indice of the predicted critical point
    FUNCTION trilinear_interpol_grad(vals,r)
      REAL(q2),DIMENSION(3),INTENT(IN) :: r
      REAL(q2),DIMENSION(3) :: trilinear_interpol_grad
      INTEGER :: p1, p2, p3
      REAL(q2) :: f1, f2, f3, g1, g2, g3
      REAL(q2),DIMENSION(8,3) :: vals
      !REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
      REAL(q2),DIMENSION(3) :: val00_, val01_, val10_, val11_
      REAL(q2),DIMENSION(3) :: val0__, val1__, val_0_, val_1_, val__0, val__1
      REAL(q2),DIMENSION(3) :: val_00, val_01, val_10, val_11
      p1 = FLOOR(r(1))
      p2 = FLOOR(r(2))
      p3 = FLOOR(r(3))
      f1 = r(1) - REAL(p1,q2)
      f2 = r(2) - REAL(p2,q2)
      f3 = r(3) - REAL(p3,q2)
      ! f1 f2 f3 are checked to be correct. 
      ! they should equal to tem for the first step
      g1 = 1._q2-f1
      g2 = 1._q2-f2
      g3 = 1._q2-f3
!     rho00_ = rho000*g3 + rho001*f3
!     rho01_ = rho010*g3 + rho011*f3
!     rho10_ = rho100*g3 + rho101*f3
!     rho11_ = rho110*g3 + rho111*f3
!     rho0__ = rho00_*g2 + rho01_*f2
!     rho1__ = rho10_*g2 + rho11_*f2
!     rho = rho0__*g1 + rho1__*f1
      val00_ = vals(1,:)*g3 + vals(5,:)*f3
      val01_ = vals(3,:)*g3 + vals(7,:)*f3
      val10_ = vals(2,:)*g3 + vals(6,:)*f3
      val11_ = vals(4,:)*g3 + vals(8,:)*f3
      val0__ = val00_*g2 + val01_*f2
      val1__ = val10_*g2 + val11_*f2
      trilinear_interpol_grad = val0__*g1 + val1__*f1
!      PRINT *, 'in function, trilinear is'
!      PRINT *, trilinear_interpol_grad 
      ! starting from z or x axis makes no difference 
      val00_ = vals(1,:)*g1 + vals(2,:)*f1
      val01_ = vals(5,:)*g1 + vals(6,:)*f1
      val10_ = vals(3,:)*g1 + vals(4,:)*f1
      val11_ = vals(7,:)*g1 + vals(8,:)*f1
      val0__ = val00_*g2 + val10_*f2
      val1__ = val01_*g2 + val11_*f2
      trilinear_interpol_grad = val0__*g3 + val1__*f3
!      PRINT *, 'in function, trilinear is'
!      PRINT *, trilinear_interpol_grad
  ! More work for gradients 
  ! In this case, gradients are input values. 
  ! it doesn't make sense to find gradients here. 
!      val_0_ = val00_*g1 + val10_*f1
!      val_1_ = val01_*g1 + val11_*f1
!      val_00 = val000*g1 + val100*f1
!      val_01 = val001*g1 + val101*f1
!      val_10 = val010*g1 + val110*f1
!      val_11 = val011*g1 + val111*f1
!      val_00 = vals(1)*g1 + vals(1)*f1
!      val_01 = vals(4)*g1 + vals(6)*f1
!      val_10 = vals(3)*g1 + vals(4)*f1
!      val_11 = vals(7)*g1 + vals(8)*f1
!      val__0 = val_00*g2 + val_10*f2
!      val__1 = val_01*g2 + val_11*f2
!      trilinear_interpol_grad(1) = val1__ - val0__
!      trilinear_interpol_grad(2) = val_1_ - val_0_
!      trilinear_interpol_grad(3) = val__1 - val__0
      ! despite the names, these interpolated forces should be in cartesian.
      ! Because this function should be receiving cartesian inputs
  !    CALL vector_matrix(rho_grad_lat, chg%car2lat, rho_grad)
    RETURN
    END FUNCTION trilinear_interpol_grad
 

    FUNCTION coorcheck(lattice)
      LOGICAL,DIMENSION(3) :: coorcheck
      REAL(q2) :: prod1,prod2,prod3
      REAL(q2),DIMENSION(3,3) :: lattice
      INTEGER :: i
      ! elements 1, 2, 3 are for whether coordinates align with x, y, z axis.
      coorcheck = (/.TRUE.,.TRUE.,.TRUE./)
      IF (lattice(1,2) .NE. 0 .OR. lattice(1,3) .NE. 0) THEN
        coorcheck(1)=.FALSE.
      END IF
      IF (lattice(2,1) .NE. 0 .OR. lattice(2,3) .NE. 0) THEN
        coorcheck(2)=.FALSE.
      END IF
      IF (lattice(3,2) .NE. 0 .OR. lattice(3,1) .NE. 0) THEN
        coorcheck(3)=.FALSE.
      END IF
    END FUNCTION
    
    SUBROUTINE getgradhes(p,chg,hes,force)
    TYPE(hessian) :: hes
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3) :: ptxy1, ptxy2, ptxy3, ptxy4, ptxz1, ptxz2, ptxz3 
    INTEGER,DIMENSION(3) :: ptxz4, ptyz1, ptyz2, ptyz3, ptyz4
    REAL(q2),DIMENSION(3) :: force
    INTEGER,DIMENSION(3) :: p, pt, ptt, ptx1, ptx2, pty1, pty2, ptz1, ptz2
    ! calculate hessian matrix, the second derivative at a point

    ! is calculated in the following way:
    ! first order derivatives are calculated between the point
    ! and its first neighbor point to get a central difference derivative
    ! the forward and backward central derivatives are used to calculate
    ! the second order derivatives
    CALL pbc(p,chg%npts)
    ptx1 = p + (/1,0,0/)
    ptx2 = p + (/-1,0,0/)
    pty1 = p + (/0,1,0/)
    pty2 = p + (/0,-1,0/)
    ptz1 = p + (/0,0,1/)
    ptz2 = p + (/0,0,-1/)
    CALL pbc(ptx1,chg%npts)
    CALL pbc(ptx2,chg%npts)
    CALL pbc(pty1,chg%npts)
    CALL pbc(pty2,chg%npts)
    CALL pbc(ptz1,chg%npts)
    CALL pbc(ptz2,chg%npts)
    ! these are the first order derivatives
    hes%du(1) = &
      (rho_val(chg,p(1),p(2),p(3)) - &
       rho_val(chg,ptx2(1),ptx2(2),ptx2(3)))
    hes%du(2) = 0.5 * & 
      (rho_val(chg,ptx1(1),ptx1(2),ptx1(3)) - &
       rho_val(chg,ptx2(1),ptx2(2),ptx2(3)))
    hes%du(3) = & 
      (rho_val(chg,ptx1(1),ptx1(2),ptx1(3)) - &
       rho_val(chg,p(1),p(2),p(3)))
    hes%dv(1) = & 
      (rho_val(chg,p(1),p(2),p(3)) - &
       rho_val(chg,pty2(1),pty2(2),pty2(3)))
    hes%dv(2) = 0.5 * & 
      (rho_val(chg,pty1(1),pty1(2),pty1(3)) - &
       rho_val(chg,pty2(1),pty2(2),pty2(3)))
    hes%dv(3) = & 
      (rho_val(chg,pty1(1),pty1(2),pty1(3)) - &
       rho_val(chg,p(1),p(2),p(3)))
    hes%dw(1) = & 
      (rho_val(chg,p(1),p(2),p(3)) - &
       rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))
    hes%dw(2) = 0.5 * & 
      (rho_val(chg,ptz1(1),ptz1(2),ptz1(3)) - &
       rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))
    hes%dw(3) = & 
      (rho_val(chg,ptz1(1),ptz1(2),ptz1(3)) - &
       rho_val(chg,p(1),p(2),p(3)))
    ! these are the second order derivatives
    hes%dudu = (hes%du(3) - hes%du(1)) 
    hes%dvdv = (hes%dv(3) - hes%dv(1)) 
    hes%dwdw = (hes%dw(3) - hes%dw(1)) 
    ! for dudv, calculate dv of half points on u direction, using
    ! extrapolated points, such as - 0.5 u, +- 1 v
    ptxy1 = p + (/-1,-1,0/)
    ptxy2 = p + (/-1,+1,0/)
    ptxy3 = p + (/+1,+1,0/)
    ptxy4 = p + (/+1,-1,0/)
    ptxz1 = p + (/-1,0,-1/)
    ptxz2 = p + (/-1,0,+1/)
    ptxz3 = p + (/+1,0,+1/)
    ptxz4 = p + (/+1,0,-1/)
    ptyz1 = p + (/0,-1,-1/)
    ptyz2 = p + (/0,-1,+1/)
    ptyz3 = p + (/0,+1,+1/)
    ptyz4 = p + (/0,+1,-1/)
    CALL pbc(ptxy1,chg%npts)
    CALL pbc(ptxy2,chg%npts)
    CALL pbc(ptxy3,chg%npts)
    CALL pbc(ptxy4,chg%npts)
    CALL pbc(ptxz1,chg%npts)
    CALL pbc(ptxz2,chg%npts)
    CALL pbc(ptxz3,chg%npts)
    CALL pbc(ptxz4,chg%npts)
    CALL pbc(ptyz1,chg%npts)
    CALL pbc(ptyz2,chg%npts)
    CALL pbc(ptyz3,chg%npts)
    CALL pbc(ptyz4,chg%npts)
!    hes%dudv = 1 / umag * & 
    hes%dudv = &
      ( &
      ! this is the backward dv
!      - 0.5_q2 / vmag * &
      - 0.25_q2 * & 
      ((rho_val(chg,ptxy2(1),ptxy2(2),ptxy2(3)) + &
          rho_val(chg,pty1(1),pty1(2),pty1(3)))  - &
        (rho_val(chg,ptxy1(1),ptxy1(2),ptxy1(3)) + &
          rho_val(chg,pty2(1),pty2(2),pty2(3)))  ) & 
      ! this is the forward dv
      + 0.25_q2 * & 
      ((rho_val(chg,ptxy3(1),ptxy3(2),ptxy3(3)) + &
          rho_val(chg,pty1(1),pty1(2),pty1(3)))  - &
        (rho_val(chg,ptxy4(1),ptxy4(2),ptxy4(3)) + &
          rho_val(chg,pty2(1),pty2(2),pty2(3)))  ) &
      )
    hes%dudw = &
      ( &
      ! this is the bacward dw
      - 0.25_q2 * &
      ((rho_val(chg,ptxz2(1),ptxz2(2),ptxz2(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))   - &
        (rho_val(chg,ptxz1(1),ptxz1(2),ptxz1(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))   )  &
      ! this is the forward dw
      + 0.25_q2 * & 
      ((rho_val(chg,ptxz3(1),ptxz3(2),ptxz3(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))   - &
        (rho_val(chg,ptxz4(1),ptxz4(2),ptxz4(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))   ) &
      ) 
  
    hes%dvdw = &
      ( &
      ! this is the bacward dw
      - 0.25_q2 * &
      ((rho_val(chg,ptyz2(1),ptyz2(2),ptyz2(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))  - &
        (rho_val(chg,ptyz1(1),ptyz1(2),ptyz1(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))  ) & 
      ! this is the forward dw
      + 0.25_q2 * & 
      ((rho_val(chg,ptyz3(1),ptyz3(2),ptyz3(3)) + &
          rho_val(chg,ptz1(1),ptz1(2),ptz1(3)))  - &
        (rho_val(chg,ptyz4(1),ptyz4(2),ptyz4(3)) + &
          rho_val(chg,ptz2(1),ptz2(2),ptz2(3)))  ) &
      )
    
      force(1) = hes%du(2)
      force(2) = hes%dv(2)
      force(3) = hes%dw(2)
    END SUBROUTINE
    

    ! This function takes in current position and grid point position in
    ! lattice, finds the distance this point is to the nearest cell, take half
    ! the distance as step size. 
    FUNCTION findstepsize(r)
      REAL(q2) :: findstepsize    
      REAL(q2),DIMENSION(3) :: r
      REAL(q2) :: f1,f2,f3,c1,c2,c3
      
      f1 = MIN(ABS(r(1)),ABS(0.5-ABS(r(1)))) 
      f2 = MIN(ABS(r(2)),ABS(0.5-ABS(r(2))))
      f3 = MIN(ABS(r(3)),ABS(0.5-ABS(r(3))))
      findstepsize = MIN(MIN(f1,f2),MIN(f1,f3))/2
    RETURN
    END FUNCTION findstepsize

    ! this funciton finds hessian of a interpolated point using interpolated
    ! nearest neighbor gradiants. Also gradiants taken in here should be in
    ! lattice.
    FUNCTION inthessian(grad,stepsize)
      REAL(q2),DIMENSION(6,3) :: grad
      REAL(q2),DIMENSION(3,3) :: inthessian
      REAL(q2) :: stepsize
      ! again, intnngrad is following this order:
      ! +x -x +y -y +z -z
      inthessian(1,1) = (grad(1,1)-grad(2,1))*0.5_q2/stepsize
      inthessian(2,2) = (grad(3,2)-grad(4,2))*0.5_q2/stepsize
      inthessian(3,3) = (grad(5,3)-grad(6,3))*0.5_q2/stepsize
      inthessian(1,2) = (grad(3,1)-grad(4,1))*0.5_q2/stepsize
      inthessian(2,1) = inthessian(1,2)
      inthessian(1,3) = (grad(5,1)-grad(6,1))*0.5_q2/stepsize
      inthessian(3,1) = inthessian(1,3)
      inthessian(2,3) = (grad(5,2)-grad(6,2))*0.5_q2/stepsize
      inthessian(3,2) = inthessian(2,3)
    ! assuming that this function is fine
    RETURN
    END FUNCTION inthessian

    FUNCTION newtonstep(car2lat,lat2car,hessian,force)
      REAL(q2),DIMENSION(3,3) :: car2lat,lat2car,hessian
      REAL(q2),DIMENSION(3) :: force,newtonstep
      newtonstep = MATMUL( inverse(hessian),force)
      newtonstep = MATMUL(car2lat,newtonstep)
    RETURN
    END FUNCTION newtonstep

  END MODULE

