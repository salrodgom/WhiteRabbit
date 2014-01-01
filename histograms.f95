PROGRAM histogram
 IMPLICIT NONE
 INTEGER            :: n_datos = 0
 INTEGER            :: n_boxs  = 200
 INTEGER            :: ierr,i,j,k,n_T,n_atoms
 REAL               :: max_,min_,suma
 CHARACTER (LEN=80) :: line
 REAL, allocatable  :: values(:),delta(:)
 LOGICAL            :: normalise = .false.
!
 OPEN(100,FILE='TIEMPO_ATOMS',STATUS='OLD',ACTION='READ',IOSTAT=IERR)
 OPEN(111,FILE="input",STATUS='OLD',ACTION='READ',IOSTAT=IERR)
 READ(100,*) n_T
 READ(100,*) n_atoms
 CLOSE(100)
 fileopen: IF( IERR == 0) THEN
  read_: DO
    READ (111,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT
    n_datos = n_datos + 1
  END DO read_
  REWIND( 111 )
 ENDIF fileopen
 ALLOCATE(values(1:n_datos) ,STAT=IERR)
 ALLOCATE(delta(0:n_boxs+1) ,STAT=IERR)
 IF(IERR/=0) STOP '[ERROR] variables sin alicatar en memoria.'
 datas: DO i=1,n_datos
   READ(111,*) values(i)
 END DO datas
 max_  = MAXVAL(values)
 min_  = MINVAL(values)
 delta(0) = min_
 DO j=1,n_boxs+1
    delta(j)=delta(j-1) + (max_ - min_ )/ REAL(n_boxs)
 ENDDO
 CALL print_history( values, delta, n_T, n_datos, n_boxs+1 )
 CLOSE(111)
 DEALLOCATE(values)
 DEALLOCATE(delta)
CONTAINS
 SUBROUTINE print_history( data, bound , n_T, j, k )
   IMPLICIT NONE
   INTEGER :: j,k
   INTEGER , intent(in) :: n_T
   REAL,     intent(in) :: data(1:j)
   REAL,     intent(in) :: bound(0:k)
   INTEGER :: i
   do i = 1,k
       WRITE(6,'(f10.4,f10.4)') bound(i), count( data <= bound(i) .and. data > bound(i-1))/REAL(n_T)
   enddo
 END SUBROUTINE print_history
END PROGRAM histogram
