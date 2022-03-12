! Projet realise dans le cadre du cours de 3eme annee de l'enseirb-matmeca AN304 "Calcul parallel"
! Auteurs : Benjamin Benon - bbenon@enseirb-matmeca.fr
!           Leo Cicerale - lcicerale@enseirb-matmeca.fr
!           Yoan Gorschka - ygorschka@enseirb-matmeca.fr



Program projet
  Use matrice
  Use param
  Use fonctions
  Use gradient
  Use MPI
  Use charge
  Use communications
  Implicit None
  Real(PR)                            :: dx, dy, alpha, beta, gamma, error
  Real(PR), Dimension(:), Allocatable :: u, u0, b, x, recep1, recep2, recep1_old, recep2_old
  Integer                             :: i1, iN, Me, Np, Statinfo, i, max, iter
  Character(len=13)                   :: name

  Call MPI_Init(Statinfo)

  Call MPI_Comm_size(MPI_COMM_WORLD,Np,Statinfo)
  Call MPI_Comm_rank(MPI_COMM_WORLD,Me,Statinfo)

  !On lit les paramètres
  Call parameters()

  !On affecte les charges aux processeurs
  Call charge_(Ny,Np,Me,i1,iN)

  !Calcul des pas
  dx=Lx/(Nx+1)
  dy=Ly/(Ny+1)

  !On crée A
  Call matA(dx,dy,alpha,beta,gamma)

  Allocate(u0(i1:iN),u(i1:iN),b(i1:iN))!,recep1(1:Nx),recep2(1:Nx),recep1_old(1:Nx),recep2_old(1:Nx))
  If(alpha_rob.eq.0) Then
    Allocate(recep1(1:Nx),recep2(1:Nx),recep1_old(1:Nx),recep2_old(1:Nx))
  Else
    Allocate(recep1(1:3*Nx),recep2(1:3*Nx),recep1_old(1:3*Nx),recep2_old(1:3*Nx))
  End If

  u0(i1:iN)=0._PR

  !Boucle en temps
  Do i = 1, Nmax
    iter=0
    error=1._PR
    recep1_old=1._PR
    recep2_old=1._PR
    Do While((error.ge.epsilon_schwarz).and.(iter.le.Nmax_schwarz))
      Call comms(u0,i1,iN,Me,Np,recep1,recep2)

      !On crée le second membre
      Call second_membre(u0,(i+1)*dt,dx,dy,b,i1,iN,Me,Np,recep1,recep2)

      !Résolution du système
      Call bicgstab(alpha,beta,gamma,b,u,i1,iN,dx,dy,Me,Np)

      error=error_schwarz(u,i1,iN,Me,Np,recep1_old,recep2_old)

      u0(i1:iN)=u(i1:iN)

      recep1_old=recep1
      recep2_old=recep2

      iter=iter+1

      If(Np .eq. 1)Then
          exit
      End If

    End Do
    !Mise à jour
    u0(i1:iN)=u(i1:iN)
  End Do

  !Écriture de la solution de chaque processeur
  Call writesol(Me,name,i1,iN)

  Deallocate(u0,u,b,recep1,recep2,recep1_old,recep2_old)

  Call MPI_Finalize(Statinfo)

CONTAINS

  Subroutine writesol(Me,name,i1,iN)
    !Écriture de la solution du processeur dans un fichier
    Integer, Intent(In)            :: Me, i1, iN
    Character(len=13), Intent(Out) :: name
    Character(len=3)               :: tn
    Integer                        :: i0, i2, i3, i
    i0 = Me/100
    i2 =( Me - 100*i0)/10
    i3 = Me - 100*i0 -10*i2
    tn = char(i0+48)//char(i2+48)//char(i3+48)
    name="sol"//tn//".txt"
    Open(unit=10,file=name)
    Do i=i1,iN
      Write(10,*) (modulo(i-1,Nx)+1)*dx, (((i-1)/(Nx))+1)*dy, u(i)
    End Do
    Close(10)
  End Subroutine writesol

End Program projet
