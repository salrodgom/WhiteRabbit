PROGRAM remove_fake_rings
 IMPLICIT NONE
 INTEGER            :: i,j,k
 INTEGER            :: atom_4(1:8),atom_6(1:12),atom_8(1:16)
 CHARACTER (LEN=80) :: line
 LOGICAL            :: FLAG = .false.
! 
 OPEN(666, file="12-ring.txt",STATUS='OLD',ACTION='READ',IOSTAT=IERR )
 OPEN(888, file="12-ring.txt",STATUS='OLD',ACTION='READ',IOSTAT=IERR )
 OPEN(606, file="12-ring.txt.out",STATUS='OLD',ACTION='WRITE',IOSTAT=IERR )
 OPEN(808, file="12-ring.txt.out",STATUS='OLD',ACTION='WRITE',IOSTAT=IERR )
 READ(5,'(A)',,IOSTAT=IERR) line
 IF( IERR == 0) THEN READ(line,*)( atom_4(j) , j=1,8 )
 read_file_12: DO
   READ(666,'(A)') line
   IF( IERR /= 0) EXIT
   READ(line,*)( atom(j) , j=1,12 )
   DO i=1,8
      IF(atom_4(i)/=atom_6(i).and.atom_4(i+1)/=atom_6(i+1).and.atom_4(i+2)/=atom_6(i+2).and.atom_4(i+3)/=atom_6(i+3)) FLAG=.true.
   END DO
 END DO read_file
END PROGRAM remove_fake_rings
