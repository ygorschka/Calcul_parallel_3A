Module param
  Implicit None
  Integer, Parameter  :: PR=8
  Real(PR), Parameter :: Pi=4._PR*Atan(1._PR)
  Integer             :: Nx, Ny, case_choice, kmax, Nmax, recouvrement, Nmax_schwarz
  Real(PR)            :: Lx, Ly, D, dt, epsilon, epsilon_schwarz, alpha_rob, beta_rob

Contains

  Subroutine parameters()
    Open(Unit=10, File='data.txt')
    Read(10,*) Nx ! 1 Nombre de points de discrétisation selon x
    Read(10,*) Ny ! 2 Nombre de points de discrétisation selon y
    Read(10,*) Lx ! 3 Longueur selon x
    Read(10,*) Ly ! 4 Longueur selon x
    Read(10,*) D ! 5 Coefficient de diffusion
    Read(10,*) dt ! 6 Pas de temps
    Read(10,*) case_choice ! 7 Choix du cas
    Read(10,*) kmax ! 8 Nombre maximal d'itérations pour le gradient conjugué
    Read(10,*) epsilon ! 9 Tolérance pour le gradient conjugué
    Read(10,*) Nmax ! 10 Nombre d'itérations en temps
    Read(10,*) recouvrement ! 11 Longueur du recouvrement
    Read(10,*) epsilon_schwarz ! 12 Tolerance pour les iterations de schwarz
    Read(10,*) Nmax_schwarz ! 13 Nombre maximal d'iterations pour schwarz
    Read(10,*) alpha_rob ! 14 Paramètre relatif à la dérivée dans les conditions de Robin
    Read(10,*) beta_rob ! 15 Paramètre relatif à Dirichlet dans les conditions de Robin
    Close(10)

  End Subroutine parameters

End Module param
