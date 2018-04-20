PROGRAM monte_carlo
  USE utilities
  IMPLICIT NONE

  INTEGER, PARAMETER :: sp = 4, dp = 8
  ! REAL(kind = sp), PARAMETER :: l = 4.0
  REAL (KIND= dp), ALLOCATABLE :: H(:), J(:,:), A(:,:), J_buff(:)
  REAL (KIND= dp):: dx, x!checksum , my_sum, tolerance, dx
  REAL (KIND= dp):: time1, time2!, time3
  REAL (KIND= dp):: T
  INTEGER (KIND=dp):: i1, i2, k, N, M, nsweep, i, nthreads
  INTEGER (KIND=sp):: i3
  ! PRINT *,'Sweeps N M Temp'
  ! READ *, nsweep
  ! READ *, nsweep, N, M, T
  nsweep = 1000
  N = 1000
  M = 1
  T = 0.1
  !
  CALL CPU_TIME(time1)
  i3 = NINT(time1)
  CALL random_SEED(i3)
  !

  ALLOCATE(H(m))
  ALLOCATE(J(n,n))
  ALLOCATE(A(m,n))
  ALLOCATE(J_buff(n))

  nthreads = 4
  CALL OMP_SET_NUM_THREADS(nthreads)

  !$OMP PARALLEL DO
  ! OPEN(unit=130,file='A',status='unknown')
  DO i1=1,N,1
     DO i2 = 1, m, 1
        CALL random_NUMBER(x)
        ! IF ( (NINT(X)) .EQ. 1 ) THEN
        A(i2,i1)=1
        ! END IF
     END DO
  END DO
  ! CLOSE(130)
  !$OMP END PARALLEL DO

  dx = 0.1
  OPEN(unit=120,file='mag',status='unknown')
  OPEN(unit=180,file='Avg_energy',status='unknown')
  !$OMP PARALLEL DO
  DO i3 = 1, 40, 1
     T = i3*dx
     CALL corr(N,J, T)
     DO i1 = 1, nsweep, 1
        DO k = 1, N, 1
           J_buff(:) = J(:,k)
           CALL DGEMM('N','N', M, 1, N, -4.0d0, A, M, J_buff, N, 0.0d0, H,M)
           H(:) = EXP(-H(:))
           CALL random_NUMBER(dx)
           DO i = 1, m, 1
              CALL random_NUMBER(dx)
              IF ( dx .LE. H(i) ) THEN
                 A(i,k) = -A(i,k)
              END IF
           END DO
        END DO
        ! CALL energy(J, A, N, M, H)

     END DO
     DO k = 1, m, 1
        H(k) = SUM(A(k,:))/N
     END DO
     WRITE(120,*) SUM(H)/m
  ENDDO

  !$OMP END PARALLEL DO
  CLOSE(120)
  CLOSE(180)





  IF (ALLOCATED(A)) DEALLOCATE(A)
  IF (ALLOCATED(H)) DEALLOCATE(H)
  IF (ALLOCATED(J)) DEALLOCATE(J)
  IF (ALLOCATED(J_buff)) DEALLOCATE(J_buff)

  CALL CPU_TIME(TIME2)
  WRITE(*,*) TIME2-TIME1

END PROGRAM monte_carlo
