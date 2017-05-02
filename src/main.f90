PROGRAM main
! {{ MAIN }}
! USE OMP_LIB
 USE vector_module
 IMPLICIT NONE
 INTEGER            :: n_atoms = 0
 INTEGER            :: n_T = 0
 REAL               :: r1  = 1.0
 REAL               :: r2  = 2.2
 INTEGER            :: IERR
 INTEGER            :: i,j,k,l,p,m,o,n_ring,n,h
 INTEGER            :: ei_atoms(16),six_atoms(14),fou_atoms(8),ten_atoms(20)
 REAL, PARAMETER    :: pi=ACOS(-1.0)
 REAL               :: atom(3),ouratom(3),dist,distance,a_ring,b_ring,delta,dist_(1:3),a_average
 REAL               :: cell_0(6),q,r,s,rv(3,3),vr(3,3),delta_ring,make_distance_car_cm,volume
 REAL               :: x1,x2,x3,average(3)!,centres(3,8),dist_(1:3)!,pointout_strgline(3)
 REAL               :: image(1:3,1:27)!,image_car(1:3,1:27),medio(1:3,1:27)
 REAL, ALLOCATABLE  :: xcryst(:,:),xinbox(:,:)
 REAL               :: xcrys_in_ring(1:3,0:20),xinbox_in_ring(1:3,0:20),d(1:20),d0(1:20),e(1:20)
 INTEGER, ALLOCATABLE  :: id(:),adj(:,:),O_ident(:)
 CHARACTER (LEN=1), ALLOCATABLE :: adj_char(:,:)
 CHARACTER (LEN=3)  :: input_type = "PDB"
 CHARACTER (LEN=6)  :: atomc
 CHARACTER (LEN=80) :: line,string
 CHARACTER (LEN=4)  :: mol
 CHARACTER (LEN=5)  :: model
 CHARACTER (LEN=6)  :: charct(6)
 CHARACTER (LEN=4)  :: typ(3)
 CHARACTER (LEN=2)  :: t1,t2,t3,t4
 CHARACTER (LEN=4), ALLOCATABLE :: label(:,:)
 CHARACTER (LEN=10) :: spacegroup = "P1"
 LOGICAL            :: compute_distance,FLAG_stop
 LOGICAL            :: FLAG = .true.
! INTEGER            :: CHUNK,NTHREADS,TID
 TYPE (vector)      :: o0,u,ou,v,ov
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID)
! TID=OMP_GET_THREAD_NUM()
! IF(TID==0)THEN
!   NTHREADS=OMP_GET_NUM_THREADS()
!   NTHREADS=1
!   WRITE(*,'(A,I3,A)')'[CPU: ',NTHREADS,' hilo/s]'
! ENDIF
!$OMP END PARALLEL
 READ(5,*) compute_distance
 READ(5,*) input_type
 WRITE(6,'(A,A)')'Input type: .', input_type
! {{ leemos el movie
 IF (input_type=='PDB') THEN
  OPEN(100,FILE='input.pdb',STATUS='OLD',ACTION='READ',IOSTAT=IERR )
  fileopen_PDB: IF( IERR == 0) THEN
   DO
    READ (100,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT
    IF (line(1:5)=='MODEL') THEN
      n_T = n_T + 1
      n_atoms = 0
    ENDIF
    IF (line(1:4)=='ATOM') n_atoms = n_atoms + 1
   END DO
  ENDIF fileopen_PDB
  REWIND( 100 )
 ELSE IF (input_type=='ARC') THEN 
  OPEN(101,FILE='input.arc',STATUS='OLD',ACTION='READ',IOSTAT=IERR )
  fileopen_ARC: IF( IERR == 0 ) THEN
   DO
    READ (101,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT
    IF (line(1:3)=='PBC') THEN
      n_T = n_T + 1
      n_atoms = 0
    ENDIF
    IF (line(1:2)=='Ge'.OR.line(1:2)=='Si'.OR.line(1:2)=='O '.or.&
        line(1:2)=='AL'.OR.line(1:2)=='SI'.OR.line(1:2)=='O '.or.&
        line(1:2)=='GE'.OR.line(1:2)=='NA'.OR.line(1:2)=='O2') n_atoms = n_atoms + 1
   END DO
  END IF fileopen_ARC
  REWIND( 101 )
 END IF
 WRITE(*,*)'Steps: ',n_T
 WRITE(*,*)'Atoms: ',n_atoms
 IF(compute_distance.EQV..false.) THEN
  OPEN(123,FILE='TIEMPO_ATOMS')
   WRITE(123,*)n_T
   WRITE(123,*)n_atoms
  CLOSE(123)
 END IF
! }}
! {{ alicatamos las variables x,y,z,labels
 ALLOCATE(xcryst(0:3,n_atoms) ,STAT=IERR)
 ALLOCATE(xinbox(3,n_atoms)   ,STAT=IERR)
 ALLOCATE(id(n_atoms)           ,STAT=IERR)
 ALLOCATE(LABEL(n_atoms,2)      ,STAT=IERR)
 IF(IERR/=0) STOP '[ERROR] variables sin alicatar en memoria.'
! }}
 IF (input_type=='PDB') OPEN(999,FILE='out.pdb')
 IF (input_type=='ARC') THEN
    OPEN(909,FILE='out.pdb')
    !WRITE(909,'(A)')'!BIOSYM archive 2'
    !WRITE(909,'(A)')'PBC=ON'
 ENDIF
 IF (compute_distance) OPEN(888,file='8-distance.txt')
 IF (compute_distance) OPEN(666,file='6-distance.txt')
 IF (compute_distance) OPEN(444,file='4-distance.txt')
 IF (compute_distance) OPEN(382,file='10-distance.txt')
 if (compute_distance) open(220,file='20-ring.txt')
 IF (compute_distance) OPEN(880,file='16-ring.txt')
 IF (compute_distance) OPEN(660,file='12-ring.txt')
 IF (compute_distance) OPEN(440,file='8-ring.txt' )
 steps: DO p=1,n_T
  type_input_question: IF(input_type=='PDB') THEN
   READ (100,'(A)') line
   WRITE(999,'(A)') line
   READ (100,'(A)') line
   READ (line,'(6x,3f9.3,3f7.2,1x,a10)') &
    cell_0(1),cell_0(2),cell_0(3),cell_0(4),cell_0(5),cell_0(6),spacegroup
   spacegroup = "P1"
   CALL cell(rv,vr,cell_0) ! crea la matriz H para cambios entre coordenadas y su inversa.
   a_average = a_average + ((volume(rv))**(1.0/3.0))/real(n_T)
   WRITE(999,'(A6,3f9.3,3f7.2,1x,a10)') &
    'CRYST1',cell_0(1),cell_0(2),cell_0(3),cell_0(4),cell_0(5),cell_0(6),spacegroup
   average(1)=0.0
   average(2)=0.0
   average(3)=0.0
! {{ comentarios de RASPA sobre la lectura y escritura de archivos PDB }}
! PDB file-format
!  1 -  6        Record name     "ATOM  "
!  7 - 11        Integer         serial        Atom serial number
! 13 - 14        Atom            name          Chemical symbol (right justified)
!      15        Remoteness indicator
!      16        Branch designator    
! 17             Character       altLoc        Alternate location indicator
! 18 - 20        Residue name    resName       Residue name
!      21        Reserved
! 22             Character       chainID       Chain identifier
! 23 - 26        Integer         resSeq        Residue sequence number
! 27             AChar           iCode         Code for insertion of residues
! 31 - 38        Real(8.3)       x             Orthogonal coordinates for X
! 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y
! 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z
! 55 - 60        Real(6.2)       occupancy     Occupancy
! 61 - 66        Real(6.2)       tempFactor    Isotropic B-factor
! 73 - 76        LString(4)      segID         Segment identifier, left-justified, may
!                                              include a space, e.g., CH86, A 1, NASE.
! 77 - 78        LString(2)      element       Element symbol, right-justified
! 79 - 80        LString(2)      charge        Charge on the atom
! Typical Format:  (6A1,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,2A4)
! Cols.  1-6    Record name "CRYST1"
!      7-15    a (Angstrom)
!     16-24    b (Angstrom)
!     25-33    c (Angstrom)
!     34-40    alpha (degrees)
!     41-47    beta  (degrees)
!     48-54    gamma (degrees)
!     56-66    Space group symbol, left justified
!     67-70    Z   Z value is the number of polymeric chains in a unit cell. In the case of heteropolymers,
!                  Z is the number of occurrences of the most populous chain.
! Typical Format:  (6A1,3F9.3,3F7.2,1X,11A1,I4)
! transformation from the orthogonal coordinates contained in the entry to fractional
! crystallographic coordinates
! # Notes:     If the orthogonal Angstroms coordinates are X, Y, Z, and the fractional
! cell coordinates are xfrac, yfrac, zfrac, then:
! xfrac = S11X + S12Y + S13Z + U1
! yfrac = S21X + S22Y + S23Z + U2
! zfrac = S31X + S32Y + S33Z + U3
! Col. 1-6     Record name SCALEn
!  1 -  6       Record name    "SCALEn" (n=1, 2, or 3)
! 11 - 20       Real(10.6)     s[n][1]                         
! 21 - 30       Real(10.6)     s[n][2]                       
! 31 - 40       Real(10.6)     s[n][3]                         
! 46 - 55       Real(10.5)     u[n] 
! Record:   MODEL
! Contains:    the model serial number when a single coordinate entry contains multiple structures
! # Notes:     Models are numbered sequentially beginning with 1.
! # If an entry contains more than 99,999 total atoms,
! then it must be divided among multiple models.
! # Each MODEL must have a corresponding ENDMDL record.
! # In the case of an NMR entry the EXPDTA record states the number of model structures
! that are present in the individual entry.
!  1 -  6       Record name    "MODEL "                                            
! 11 - 14       Integer        Model serial number
! Record: ENDMDL
! Contains:    these records are paired with MODEL records to group individual structures found in a coordinate entry
! # Notes:     MODEL/ENDMDL records are used only when more than one structure
! is presented in the entry, or if there are more than 99,999 atoms.
! # Every MODEL record has an associated ENDMDL record. 
!  1 -  6         Record name      "ENDMDL"
   read_coor_PDB: DO i=1,n_atoms
    READ(100, '(A)', iostat = IERR ) line
    IF(line(1:6)=='ENDMDL') exit read_coor_PDB
    IF (IERR/=0) EXIT read_coor_PDB          ! formatos de lectura para PDB 
       atomc = line(1:6)                     ! leemos
       mol   = line(18:20)
       READ(line(7:11),*)  id(i)
       READ(line(31:38),*) xinbox(1,id(i))
       READ(line(39:46),*) xinbox(2,id(i))
       READ(line(47:54),*) xinbox(3,id(i))
       label(id(i),1)= line(13:14)
       label(id(i),2)= label(id(i),1)
    FORALL ( j = 1:3 )
       average(j)=average(j)+xinbox(j,id(i))/real(n_atoms)
    END FORALL
   ENDDO read_coor_PDB
   print_coor_PDB: DO i=1,n_atoms ! definimos las coordenadas cartesianas
    FORALL ( j=1:3 )
       xinbox(j,id(i)) = xinbox(j,id(i)) - average(j)
    END FORALL
    FORALL ( j=1:3 )
       xcryst(j,id(i)) = vr(j,1)*xinbox(1,id(i))+vr(j,2)*xinbox(2,id(i))+vr(j,3)*xinbox(3,id(i))
    end FORALL
    x1=xinbox(1,id(i))
    x2=xinbox(2,id(i))
    x3=xinbox(3,id(i))
    typ(1)=label(id(i),1)
    typ(2)=label(id(i),2)
!ATOM      1 Ge   MOL  ****       0.331  -6.716 -14.604******  0.00  Ge
    WRITE(999,'(a6,i5,1x,2(a4,1x),i4,4x,3f8.3,2f6.2,2X,a4)') &
     atomc,id(i),typ(1),mol,0,x1,x2,x3,0.0,0.0,typ(2)
   ENDDO print_coor_PDB
   READ (100,'(A)') line
   WRITE(999,'(A)') line
! {{ procedimiento ante archivo ARC }}
  ELSE IF (input_type=='ARC') THEN
   READ(101,'(A)', iostat = IERR ) line
   IF (IERR/=0) EXIT steps
   do1: DO
    IF(line(1:3)=='PBC')THEN
      EXIT do1
    ELSE
      READ (101,'(A)') line
    END IF
   END DO do1
   READ(line,*),mol,(cell_0(i),i=1,6)
   !READ(line(6:13),*) cell_0(1)
   !READ(line(16:23),*)cell_0(2)
   !READ(line(26:33),*)cell_0(3)
   !READ(line(36:43),*)cell_0(4)
   !READ(line(46:53),*)cell_0(5)
   !READ(line(56:63),*)cell_0(6)
   CALL cell(rv,vr,cell_0)
   WRITE(909,'(A5,i5)')'MODEL',p
   WRITE(909,'(A6,3f9.3,3f7.2,1x,a10)') &
    'CRYST1',cell_0(1),cell_0(2),cell_0(3),cell_0(4),cell_0(5),cell_0(6),spacegroup
   read_coor_ARC: DO i=1,n_atoms
     READ(101, '(A)', iostat = IERR ) line
     IF (IERR/=0) EXIT read_coor_ARC
     id(i)=i
!     READ(line,'(a4,2x,3(f13.9,1x),a80)') typ(1),x1,x2,x3,string
     READ(line,*) typ(1),x1,x2,x3,string
     xinbox(1,id(i))=x1
     xinbox(2,id(i))=x2
     xinbox(3,id(i))=x3
     label(id(i),1)=typ(1)
     label(id(i),2)=label(id(i),1)
     FORALL ( j=1:3 )
       xcryst(j,id(i)) = mod(vr(j,1)*xinbox(1,id(i))+vr(j,2)*xinbox(2,id(i))+vr(j,3)*xinbox(3,id(i))+100.0,1.0)
     END FORALL
     typ(2)=typ(1)
     mol=' MOL'
     atomc='ATOM  '
     q=0.0
     r=0.0
     m=0
     WRITE(909,'(a6,i5,1x,2(a4,1x),i4,4x,3f8.3,2f6.2,2X,a4)') &
      atomc,id(i),typ(1),mol,m,xinbox(1,id(i)),xinbox(2,id(i)),xinbox(3,id(i)),q,r,typ(2)
   ENDDO read_coor_ARC
   READ (101,'(A)') line
   READ (101,'(A)') line
   WRITE(909,'(A)')'ENDMDL'
  END IF type_input_question
  call escritura_cif(xcryst,n_atoms,LABEL,cell_0,rv)
! {{ para calcular propiedades geometricas del grafo calculo primero la matriz de adyacencias adj }}
  make_graph: IF( FLAG.EQV..true. ) THEN
    ALLOCATE(adj(1:n_atoms,1:n_atoms)      ,STAT=IERR)
    IF(IERR/=0) STOP '[ERROR] variable adj sin alicatar en memoria.'
    ALLOCATE(adj_char(1:n_atoms,1:n_atoms) ,STAT=IERR)
    IF(IERR/=0) STOP '[ERROR] variable adj_char sin alicatar en memoria.'
    ALLOCATE(O_ident(1:n_atoms)            ,STAT=IERR)
    IF(IERR/=0) STOP '[ERROR] variable O_ident sin alicatar en memoria.'
    pairs: DO i=1,n_atoms
      DO j=i,n_atoms
      adj(i,j)=0.0
      adj_CHAR(i,j)=' '
      IF(j/=i)THEN
       FORALL ( k=1:3 )
          atom(k)= xcryst(k,j)
          ouratom(k) = xcryst(k,i)
       END FORALL
! {{ condicion de enlace en un par T-O 
       call make_distances(.false.,cell_0,ouratom,atom,rv,dist_,r)
       IF(r>=r1.AND.r<r2)THEN
         adj(i,j)=1.0
         adj_CHAR(I,J)='@'
       ENDIF
! }}
      ENDIF
      adj(j,i)=adj(i,j)
      adj_CHAR(J,I)=adj_CHAR(I,J)
     ENDDO
    ENDDO pairs
    CALL write_grafo(adj,n_atoms,id,222)
    CALL ESCRITURA_GRAFO(adj,n_atoms,123)
! {{ el siguiente bloque revida la conectividad de cada tetraedro, debe ser 4
! independiente del grupo espacial SPACE
   conectivity: DO i=1,n_atoms ! conectividad de los ATOM  ! k=4 y k=2    en zeolitas.
     k=0
     DO j=1,n_atoms
        k=k+adj(i,j)
     ENDDO
     xcryst(0,i)=k
!     IF(k/=0.0.and.k/=2.0.and.k/=4.0) WRITE(6,'(3A,I2)') '[OJO] Connectivity of ',label(i,1),'is ',k
   ENDDO conectivity
   WRITE(6,'(A)')'[loops] Connectivity between nodes [Degree]:'
   IF(compute_distance.eqv..false.) WRITE(*,'(1000(I1))')(int(xcryst(0,k)),k=1,n_atoms)
   DEALLOCATE(adj,adj_char,O_ident)
   FLAG=.false.
  ENDIF make_graph
  RING_PROPERTIES: IF (compute_distance) THEN
  n_ring=0
  distancias_4: do
! 440 es 8-ring.txt  
    READ(440,*,IOSTAT=IERR)( fou_atoms(j), j=1,8 )
    IF(IERR/=0) EXIT distancias_4
    d(1:16)=0.0
    d0(1:16)=0.0
    e(1:16)=0.0
    n_ring=n_ring+1
! {{ coloco todos los atomos en el mismo sistema de referencia y calculo
! las distancias sin PBC }}
    elijo_pivote_8: FORALL ( k=1:3 , j=1:8 )
       xcrys_in_ring(k,j) = xcryst(k,fou_atoms(j))
    END FORALL elijo_pivote_8
    DO j=2,8
       pivoteo_8: FORALL ( k=1:3 )
          atom(k)   = xcrys_in_ring(k,j)   ! { atomo que cambia }
          ouratom(k)= xcrys_in_ring(k,j-1) ! { atomo pivote }
       END FORALL pivoteo_8
       CALL make_distances(.true.,cell_0,atom,ouratom,rv,dist_,r)
       FORALL ( k=1:3 )
          xcrys_in_ring(k,j) = dist_(k)
       END FORALL
    END DO
! {{ paso a coordenadas cartesianas
    DO i=1,8
     FORALL ( j=1:3 )
       xinbox_in_ring(j,i) = rv(j,1)*xcrys_in_ring(1,i)+rv(j,2)*xcrys_in_ring(2,i)+rv(j,3)*xcrys_in_ring(3,i)
     END FORALL
    END DO 
! }}
!   {{ calculo el centro geometrico del anillo
    x1=0.0
    DO j=1,3
       x1=0.0
       DO l=1,8
          x1=x1+xinbox_in_ring(j,l)/8.0
       ENDDO
       xinbox_in_ring(j,0)=x1
    ENDDO
    o0%x= xinbox_in_ring(1,0)
    o0%y= xinbox_in_ring(2,0)
    o0%z= xinbox_in_ring(3,0)
!   }}
    x1=0.0 ! {{ partial window area }}
    area_8: DO j=1,8
      l=j+1
      IF (l>8) l=l-8
      u%x= xinbox_in_ring(1,j)
      u%y= xinbox_in_ring(2,j)
      u%z= xinbox_in_ring(3,j)
      v%x= xinbox_in_ring(1,l)
      v%y= xinbox_in_ring(2,l)
      v%z= xinbox_in_ring(3,l)
      ou = vector_sub(u,o0)
      ov = vector_sub(v,o0)
      x1=x1+0.5*absvec(cross(ou,ov))
    ENDDO area_8
! }}
    atom1_8: DO j=1,8
     atom2_8: DO l=1,8
      IF( (j/=l).and. &
       ((label(fou_atoms(j),1)==" O  ".or.label(fou_atoms(j),1)=="O   ").and. &
        (label(fou_atoms(l),1)==" O  ".or.label(fou_atoms(l),1)=="O   "))) THEN
       FORALL ( k=1:3 )
         atom(k)    = xcrys_in_ring(k,j) !xcryst(k,ei_atoms(j))
         ouratom(k) = xcrys_in_ring(k,l) !xcryst(k,ei_atoms(l))
       END FORALL
       d0(l)=DISTANCE(atom,ouratom,rv)
      ELSE
       d0(l)=0.0
      ENDIF
     ENDDO atom2_8
     d(j)=MAXVAL(d0) ! de cada vuelta cojo el diametor mas grande
     e(j)=d(j)
     IF(d(j)<=0.001) d(j)=9999.0
    ENDDO atom1_8
    limpia_8: DO j=9,16 ! valores no validos de d.
       d(j)=9999.0      ! los relleno con infinitos
    ENDDO limpia_8       
    b_ring=MINVAL(d) ! de los diametros calculados cojo el mas pequeño
    a_ring=MAXVAL(e)
    delta_ring = 0.5*abs( b_ring - a_ring )
    ! x1 es el area, en este caso
    WRITE(444,*)p,n_ring,b_ring,delta_ring,x1
   ENDDO distancias_4
! }}
! {{ propiedades de los anillos de 12.
   n_ring=0
   distancias_6: do
    READ(660,*,IOSTAT=IERR)( six_atoms(j), j=1,12 )
    IF(IERR/=0) EXIT distancias_6
    d(1:16)=0.0
    e(1:16)=0.0
    d0(1:16)=0.0
    n_ring=n_ring+1
    elijo_pivote_12: FORALL ( k=1:3 , j=1:12 )
       xcrys_in_ring(k,j) = xcryst(k,six_atoms(j))
    END FORALL elijo_pivote_12
    DO j=2,12
       pivoteo_12: FORALL ( k=1:3 )
          atom(k)   = xcrys_in_ring(k,j)
          ouratom(k)= xcrys_in_ring(k,j-1) ! { atomo pivote }
       END FORALL pivoteo_12
       CALL make_distances(.true.,cell_0,atom,ouratom,rv,dist_,r)
       FORALL ( k=1:3 )
          xcrys_in_ring(k,j) = dist_(k)
       END FORALL
    END DO
! {{ paso a coordenadas cartesianas
     FORALL ( i=1:12, j=1:3 )
       xinbox_in_ring(j,i) = rv(j,1)*xcrys_in_ring(1,i)+rv(j,2)*xcrys_in_ring(2,i)+rv(j,3)*xcrys_in_ring(3,i)
     END FORALL
! }}
!   {{ calculo el centro geometrico del anillo
    x1=0.0
    DO j=1,3
       x1=0.0
       DO l=1,12
          x1=x1+xinbox_in_ring(j,l)/12.0
       ENDDO
       xinbox_in_ring(j,0)=x1
    ENDDO
    o0%x= xinbox_in_ring(1,0)
    o0%y= xinbox_in_ring(2,0)
    o0%z= xinbox_in_ring(3,0)
!   }}
    x1=0.0 ! {{ partial window area }}
    area_12: DO j=1,12
      l=j+1
      IF (l>12) l=l-12
      u%x= xinbox_in_ring(1,j)
      u%y= xinbox_in_ring(2,j)
      u%z= xinbox_in_ring(3,j)
      v%x= xinbox_in_ring(1,l)
      v%y= xinbox_in_ring(2,l)
      v%z= xinbox_in_ring(3,l)
      ou = vector_sub(u,o0)
      ov = vector_sub(v,o0)
      x1=x1+0.5*absvec(cross(ou,ov))
    ENDDO area_12
! }}
! {{ 
    atom1_12: DO j=1,12
      atom2_12: DO l=1,12
       IF( (j/=l).and. &
        ((label(six_atoms(j),1)==" O  ".or.label(six_atoms(j),1)=="O   ").and. &
         (label(six_atoms(l),1)==" O  ".or.label(six_atoms(l),1)=="O   "))) THEN
        FORALL ( k=1:3 )
         atom(k)    = xcrys_in_ring(k,j) !xcryst(k,ei_atoms(j))
         ouratom(k) = xcrys_in_ring(k,l) !xcryst(k,ei_atoms(l))
        END FORALL
        d0(l)=DISTANCE(atom,ouratom,rv)
       ELSE
        d0(l)=0.0
       ENDIF
      ENDDO atom2_12
       d(j)=MAXVAL(d0) ! de cada vuelta cojo el diametor mas grande
       e(j)=d(j)
       IF(d(j)<=0.001) d(j)=9999.0
     ENDDO atom1_12
     limpia_12: DO j=13,16
       d(j) = 9999.0
     END DO limpia_12
     b_ring=MINVAL(d) ! de los diametros calculados cojo el mas pequeño
     a_ring=MAXVAL(e)
     delta_ring = 0.5*abs( a_ring - b_ring )
     WRITE(666,*)p,n_ring,b_ring,delta_ring,x1
   ENDDO distancias_6
! }}
! {{ calculo propiedades de los anillos de 16.
   n_ring=0
   distancias_8: do
    d(1:16)  = 0.0
    e(1:16)  = 0.0
    d0(1:16) = 0.0
    READ(880,*,IOSTAT=IERR)( ei_atoms(j), j=1,16 )
    IF(IERR/=0) EXIT distancias_8
    n_ring=n_ring+1
! {{ BLOQUE: minima distancia entre atomos opuestos
! independiente de si empieza por T o por O.
     forall ( k=1:3 , j=1:16 )
       xcrys_in_ring(k,j) = xcryst(k,ei_atoms(j))
     end forall    
     do j=2,16
       forall ( k=1:3 )
          atom(k)   = xcrys_in_ring(k,j)   
          ouratom(k)= xcrys_in_ring(k,j-1) ! { atomo pivote }
       end forall
       CALL make_distances(.true.,cell_0,atom,ouratom,rv,dist_,r)
       forall ( k=1:3 )
          xcrys_in_ring(k,j) = dist_(k)
       end forall
     end do
! {{ paso a coordenadas cartesianas
     FORALL ( i=1:16, j=1:3 )
       xinbox_in_ring(j,i) = rv(j,1)*xcrys_in_ring(1,i)+rv(j,2)*xcrys_in_ring(2,i)+rv(j,3)*xcrys_in_ring(3,i)
     END FORALL
! }}
!   {{ calculo el centro geometrico del anillo
    x1=0.0 
    DO j=1,3
       x1=0.0
       DO l=1,16
          x1=x1+xinbox_in_ring(j,l)/16.0
       ENDDO
       xinbox_in_ring(j,0)=x1
    ENDDO
    o0%x= xinbox_in_ring(1,0)
    o0%y= xinbox_in_ring(2,0)
    o0%z= xinbox_in_ring(3,0)
! {{ if el CG del anillo está a una distancia de un 
    do j=1, n_atoms
     do i=1,16
      if(id(j)/=ei_atoms(j)) then
       atom(1)=o0%x
       atom(2)=o0%y
       atom(3)=o0%z
       ouratom(1)=xinbox(1,id(j))
       ouratom(2)=xinbox(2,id(j))
       ouratom(3)=xinbox(3,id(j))
       if( DISTANCE(atom,ouratom,rv) <= 2.0 ) cycle distancias_8
      end if
     end do
    end do
!   }}
    x1=0.0 ! {{ partial window area }}
    area_16: DO j=1,16
      l=j+1
      IF (l>16) l=l-16
      u%x= xinbox_in_ring(1,j)
      u%y= xinbox_in_ring(2,j)
      u%z= xinbox_in_ring(3,j)
      v%x= xinbox_in_ring(1,l)
      v%y= xinbox_in_ring(2,l)
      v%z= xinbox_in_ring(3,l)
      ou = vector_sub(u,o0)
      ov = vector_sub(v,o0)
      x1=x1+0.5*absvec(cross(ou,ov))
    ENDDO area_16
! }}
     atom1_16: DO j=1,16
       atom2_16: DO l=1,16
         IF( (j/=l).and. &
          ((label(ei_atoms(j),1)==" O  ".or.label(ei_atoms(j),1)=="O   ").and. &
           (label(ei_atoms(l),1)==" O  ".or.label(ei_atoms(l),1)=="O   "))) THEN
           FORALL (k=1:3)
             atom(k)    = xcrys_in_ring(k,j) !xcryst(k,ei_atoms(j))
             ouratom(k) = xcrys_in_ring(k,l) !xcryst(k,ei_atoms(l))
           END FORALL
           d0(l)=DISTANCE(atom,ouratom,rv) !CALL make_distances(PBC,cell_0,atom,ouratom,rv,dist_,r)
         ELSE
           d0(l)=0.00001
         ENDIF
       ENDDO atom2_16
       d(j)=MAXVAL(d0) ! de cada vuelta cojo el diametor mas grande
       e(j)=d(j)
       IF(d(j)<=0.001) d(j)=9999.0
     ENDDO atom1_16
     b_ring=MINVAL(d) ! de los diametros calculados cojo el mas pequeño
     a_ring=MAXVAL(e)
     delta_ring=0.5*abs(b_ring-a_ring)
     WRITE(888,*)p,n_ring,b_ring,delta_ring,x1,a_ring
   ENDDO distancias_8
!  
! {{ calculo propiedades de los anillos de 20:
   n_ring=0
   distancias_10: do
    d(1:20)  = 0.0
    e(1:20)  = 0.0
    d0(1:20) = 0.0
    read(220,*,IOSTAT=IERR)( ten_atoms(j), j=1,20 )
    IF(IERR/=0) EXIT distancias_10
    n_ring=n_ring+1
! {{ BLOQUE: minima distancia entre atomos opuestos
! independiente de si empieza por T o por O.
     forall ( k=1:3 , j=1:20 )
       xcrys_in_ring(k,j) = xcryst(k,ten_atoms(j))
     end forall
     do j=2,20
       forall ( k=1:3 )
          atom(k)   = xcrys_in_ring(k,j)
          ouratom(k)= xcrys_in_ring(k,j-1) ! { atomo pivote }
       end forall
       CALL make_distances(.true.,cell_0,atom,ouratom,rv,dist_,r)
       forall ( k=1:3 )
        xcrys_in_ring(k,j) = dist_(k)
       end forall
     end do
! {{ paso a coordenadas cartesianas
     FORALL ( i=1:20, j=1:3 )
       xinbox_in_ring(j,i) = rv(j,1)*xcrys_in_ring(1,i)+rv(j,2)*xcrys_in_ring(2,i)+rv(j,3)*xcrys_in_ring(3,i)
     END FORALL
! }}
!   {{ calculo el centro geometrico del anillo
    x1=0.0
    DO j=1,3
       x1=0.0
       DO l=1,20
          x1=x1+xinbox_in_ring(j,l)/20.0
       ENDDO
       xinbox_in_ring(j,0)=x1
    ENDDO
    o0%x= xinbox_in_ring(1,0)
    o0%y= xinbox_in_ring(2,0)
    o0%z= xinbox_in_ring(3,0)
! {{ if el CG del anillo está a una distancia de un 
    do j=1, n_atoms
     do i=1,20
      if(id(j)/=ten_atoms(j)) then
       atom(1)=o0%x
       atom(2)=o0%y
       atom(3)=o0%z
       ouratom(1)=xinbox(1,id(j))
       ouratom(2)=xinbox(2,id(j))
       ouratom(3)=xinbox(3,id(j))
       if( DISTANCE(atom,ouratom,rv) <= 4.0 ) cycle distancias_10
      end if
     end do
    end do
!   }}

!   }}
    x1=0.0 ! {{ partial window area }}
    area_20: DO j=1,20
      l=j+1
      IF (l>20) l=l-20
      u%x= xinbox_in_ring(1,j)
      u%y= xinbox_in_ring(2,j)
      u%z= xinbox_in_ring(3,j)
      v%x= xinbox_in_ring(1,l)
      v%y= xinbox_in_ring(2,l)
      v%z= xinbox_in_ring(3,l)
      ou = vector_sub(u,o0)
      ov = vector_sub(v,o0)
      x1=x1+0.5*absvec(cross(ou,ov))
    ENDDO area_20
! {{
    atom1_20: DO j=1,20
      atom2_20: DO l=1,20
       IF( (j/=l).and. &
        ((label(ten_atoms(j),1)==" O  ".or.label(ten_atoms(j),1)=="O   ").and. &
         (label(ten_atoms(l),1)==" O  ".or.label(ten_atoms(l),1)=="O   "))) THEN
        FORALL ( k=1:3 )
         atom(k)    = xcrys_in_ring(k,j) !xcryst(k,ei_atoms(j))
         ouratom(k) = xcrys_in_ring(k,l) !xcryst(k,ei_atoms(l))
        END FORALL
        d0(l)=DISTANCE(atom,ouratom,rv)
       ELSE
        d0(l)=0.0
       ENDIF
      ENDDO atom2_20
       d(j)=MAXVAL(d0) ! de cada vuelta cojo el diametor mas grande
       e(j)=d(j)
       IF(d(j)<=0.001) d(j)=9999.0
     ENDDO atom1_20
     limpia_20: DO j=17,20
       d(j) = 9999.0
     END DO limpia_20
     b_ring=MINVAL(d) ! de los diametros calculados cojo el mas pequeño
     a_ring=MAXVAL(e)
     delta_ring = 0.5*abs( a_ring - b_ring )
     WRITE(382,*)p,n_ring,b_ring,delta_ring,x1
   ENDDO distancias_10
! }}
! {{ rebobino los ficheros de los anillos }}
   rewind ( 220 )
   REWIND ( 880 )
   REWIND ( 660 )
   REWIND ( 440 )
  ENDIF RING_PROPERTIES
 ENDDO steps
 IF (compute_distance)  CLOSE( 660 )
 if (compute_distance)  close( 382 )
 if (compute_distance)  close( 220 )
 IF (compute_distance)  CLOSE( 666 )
 IF (compute_distance)  CLOSE( 888 )
 IF (compute_distance)  CLOSE( 444 )
 IF (compute_distance)  CLOSE( 440 )
 IF (compute_distance)  CLOSE( 880 )
 CLOSE( 100 )
 CLOSE( 999 )
 CLOSE( 555 )
 WRITE(6,'(A,1x,f10.5)')'cell average: .', a_average
 IF (compute_distance.EQV..false.) STOP "PASO 1"
 IF (compute_distance.EQV..true.) STOP "PASO 2"
END PROGRAM main
!
SUBROUTINE PIKSRT(N,ARR)
! {{ ... }}
 IMPLICIT NONE
 INTEGER :: n,j,i
 REAL    :: a,arr(n)
 DO j=2,n
    a=ARR(j)
    do i=j-1,1,-1
      if (ARR(i)<=a) goto 10
      ARR(i+1)=ARR(i)
    end do
    i=0
10  ARR(i+1)=a
 END DO
 RETURN
END SUBROUTINE
!
SUBROUTINE make_distances(flag,cell_0,r2,r1,rv,r3,dist)
! {{ prepara para el calculo de distancias en una celda triclinica }}
 IMPLICIT NONE
 REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist,r3(1:3)                     ! matriz de distancias N X N
 REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
 REAL                 :: distance,rcm(3),phi
 INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
 REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
 LOGICAL              :: flag
! {{ calculamos la matriz de distancias }}
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
   enddo
  enddo
  dist=MINVAL(d_image)
  if(flag)then
   phi=1000.0
   k=1
   do l=1,27
    if(d_image(l)<=phi)then
 !     PRINT*,d_image(l),( image(m,l), m=1,3 )
      phi=d_image(l) ! seleccionamos el parametro menor
      k=l            ! y el contador correspondiente.
    endif
   enddo
   forall ( l=1:3) 
     r3(l)=image(l,k)
   end forall
  else
   r3(1:3)=0.0
  end if
 RETURN
END SUBROUTINE
!
SUBROUTINE StrgLine2PointDist( Point1, Point2, dist1, pointout_strgline)
 IMPLICIT NONE
 REAL, intent(in)  :: Point1(1:3),Point2(1:3),dist1
 REAL, intent(out) :: pointout_strgline(1:3)
 REAL              :: dist1_2,d(3)
 INTEGER           :: j
 FORALL ( j=1:3 )
   d(j)=Point1(j)-Point2(j)
 END FORALL
 dist1_2=sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
 FORALL ( j=1:3 )
  pointout_strgline(j) = Point1(j) + (Point2(j) - Point1(j))*dist1/dist1_2
 END FORALL
 RETURN
END SUBROUTINE
!
REAL FUNCTION DISTANCE(atom,ouratom,rv)
 IMPLICIT NONE
 INTEGER :: j
 REAL :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
 REAL :: rv(3,3)
 FORALL ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
 END FORALL
 DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
END FUNCTION
!
REAL FUNCTION make_distance_car_cm(atom,ouratom,rv,medio)
 IMPLICIT NONE
 INTEGER :: j
 REAL :: atom(3),ouratom(3),medio(3),dist(3),o_atom(3),o_ouratom(3)
 REAL :: rv(3,3)
 FORALL ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j) - medio(j)
 END FORALL
 make_distance_car_cm = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
END FUNCTION
!
SUBROUTINE write_grafo(A,N,ident,U)
 integer, intent(in) :: N,A(N,N),U,ident(N)
 integer :: X,Y
 OPEN(U,FILE="grafo.red")
 nodos: DO X=1,N
   DO Y=X+1,N
      IF(A(X,Y)>0.0) THEN
        IF(ident(X)>ident(Y))THEN
         WRITE(U,'(I5,1x,I5,1x,A)') ident(X),ident(Y),'+'
        ELSE
         WRITE(U,'(I5,1x,I5,1x,A)') ident(Y),ident(X),'+'
        ENDIF
      ENDIF
   ENDDO
 ENDDO nodos
 WRITE(*,*)'[loops] Escrito grafo.red en formato: i j +'
 CLOSE(U)
 RETURN
END SUBROUTINE write_grafo
!
SUBROUTINE ESCRITURA_GRAFO(A,N,U)
!--------------------------------------------
! ADAPTAMOS LA MATRIS DE ADYACENCIA PARA QUE EL PROGRAMA DE
! REPRESENTACION DE GRAFOS LO ENTIENDA. LA SALIDA LA TIENEN UN ARCHIVO
! LLAMADO: "GRAFO.DOT"
!--------------------------------------------
 INTEGER N,X,Y,A(N,N),U
 OPEN(U,file='grafo.dot')
 write(U,*)'graph G {'
 WRITE(U,*)'node[label="",height=0.1,width=0.1,fontsize=1]'
 do X=1,N
    do Y=X+1,N
       if(A(X,Y)>0.0)then
         write(U,'(A,I4,A,I4,A)')'   ',X,' -- ',Y,' ;'
       endif
    enddo
 enddo
 write(U,*)'}'
 close(u)
 RETURN
END SUBROUTINE
