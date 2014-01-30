PROGRAM Disc_Fourier_Transform
 IMPLICIT NONE
 INTEGER            :: n_t  = 0
 INTEGER            :: n_k  = 0
 REAL               :: tau  = 1.0
 INTEGER            :: IERR = 0
 INTEGER            :: k,t
 REAL, ALLOCATABLE  :: a(:),alpha(:),w(:),phi(:)
 REAL               :: ar,ai,cte,beta
 REAL, PARAMETER    :: pi=ACOS(-1.0)
 CHARACTER (LEN=80) :: line
 INTEGER            :: CHUNK,NTHREADS,TID
!
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID) 
! TID=OMP_GET_THREAD_NUM()
! IF(TID==0)NTHREADS=OMP_GET_NUM_THREADS()
! NTHREADS=1
!$OMP END PARALLEL
 OPEN(100,FILE='input',STATUS='OLD',ACTION='READ',IOSTAT=IERR )
 OPEN(999,file='out' )
 file_open: IF( IERR == 0) THEN
  DO
    READ(100,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT
    n_t=n_t+1
  ENDDO
 ENDIF file_open
 REWIND(100)
! {{ reservamos
 allocate(a(0:n_t-1))
 allocate(w(0:n_t-1))
 allocate(alpha(0:n_t-1))
 allocate(phi(0:n_t-1))
! }} {{ numero de datos
 read_file: DO k=0,n_t-1
   READ(100,*) t,a(t-1)
 ENDDO read_file
!
 beta=2.0*pi/(n_t*tau)
 cte=beta/sqrt(2.0*pi)
!$OMP PARALLEL DEFAULT(SHARED)
 DO k=1,int(n_t/2.0)-1
    ar=0.0
    ai=0.0
    DO t=0,n_t-1
       ar=ar+cte*a(t)*cos(beta*k*t*tau)
       ai=ai+cte*a(t)*sin(beta*k*t*tau)
    ENDDO
    w(k) = beta*k 
    alpha(k)= sqrt(ar*ar+ai*ai)
    phi(k)  = sign(pi/2.0,ai) - atan(ar/ai)
    WRITE(999,*) w(k),alpha(k),phi(k)
 ENDDO
!$OMP END PARALLEL
 deallocate(a,alpha,w,phi)
 CLOSE(100)
 CLOSE(999)
!
END PROGRAM
