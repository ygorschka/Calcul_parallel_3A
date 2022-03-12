Module matrice
  Use param
  Use fonctions
  Use MPI
  Use charge
  Implicit None

CONTAINS

  Subroutine matA(dx,dy,alpha,beta,gamma)
    !Construit les 3 coefficients de A
    Real(PR), Intent(In)  :: dx, dy
    Real(PR), Intent(Out) :: alpha, beta, gamma

    alpha = 1._PR+ 2._PR*D*dt*(1._PR/(dx**2)+1._PR/(dy**2))
    beta = -D*dt/(dx**2)
    gamma = -D*dt/(dy**2)

  End Subroutine matA

  Function MatVect(alpha,beta,gamma,x,i1,iN,dx,dy,Me,Np) Result(prod)
    !Effectue le produit de A par x du processeur
    Real(PR), Intent(In)                      :: alpha, beta, gamma, dx, dy
    Integer, Intent(In)                       :: i1, iN, Me, Np
    Real(PR), Dimension(i1:iN), Intent(In)    :: x
    Real(PR), Dimension(i1:iN)                :: prod
    Integer                                   :: i

    prod=0

    Do i=i1+1,iN-1
          prod(i)=prod(i)+alpha*x(i)
          If (modulo(i,Nx)==1) Then
             prod(i)=prod(i)+beta*x(i+1)
          ElseIf (modulo(i,Nx)==0) Then
             prod(i)=prod(i)+beta*x(i-1)
          Else
             prod(i)=prod(i)+beta*(x(i+1)+x(i-1))
          End If
    End Do

    prod(i1)=prod(i1)+alpha*x(i1)
    If (modulo(i1,Nx)==1) Then
       prod(i1)=prod(i1)+beta*x(i1+1)
    End If

    prod(iN)=prod(iN)+alpha*x(iN)

    If (modulo(iN,Nx)==0) Then
       prod(iN)=prod(iN)+beta*x(iN-1)
    End If

    Do i=1,iN-Nx-i1+1
      prod(i+i1-1)=prod(i+i1-1)+gamma*x(i+i1+Nx-1)
    End Do
    Do i=i1-iN+2*Nx,Nx
      prod(iN-Nx+i)=prod(iN-Nx+i)+gamma*x(iN-2*Nx+i)
    End Do

    If(alpha_rob.eq.0) Then
      If (Me.ne.0) Then
        Do i=i1,i1+Nx-1
          prod(i)=x(i)/(dx**2+dy**2)
        End Do
      End If
      If (Me.ne.(Np-1)) Then
        Do i=iN-Nx+1,iN
          prod(i)=x(i)/(dx**2+dy**2)
        End Do
      End If
    Else
      If(Me.ne.0) Then

        Do i=i1,i1+Nx-1
          prod(i)=(1._PR+D*dt*(2._PR/(dx**2)+2._PR/(dy**2)+2._PR*beta_rob/(alpha_rob*dy)))*x(i)
        End Do

        Do i=i1+1,i1+Nx-2
          prod(i)=prod(i)-D*dt*(x(i-1)+x(i+1))/(dx**2)
        End Do

        prod(i1)=prod(i1)-D*dt*x(i1+1)/(dx**2)
        prod(i1+Nx-1)=prod(i1+Nx-1)-D*dt*x(i1+Nx-2)/(dx**2)

        Do i=i1+Nx,i1+2*Nx-1
          !prod(i)=-2._PR*D*dt/(dy**2)*x(i)
          prod(i-Nx)=prod(i-Nx)-2._PR*D*dt/(dy**2)*x(i)
        End Do

      End If
      If(Me.ne.(Np-1)) Then

        Do i=iN-Nx+1,iN
          prod(i)=(1._PR+D*dt*(2._PR/(dx**2)+2._PR/(dy**2)+2._PR*beta_rob/(alpha_rob*dy)))*x(i)
        End Do

        Do i=iN-Nx+2,iN-1
          prod(i)=prod(i)-D*dt*(x(i-1)+x(i+1))/(dx**2)
        End Do

        prod(iN-Nx+1)=prod(iN-Nx+1)-D*dt*x(iN-Nx+2)/(dx**2)
        prod(iN)=prod(iN)-D*dt*x(iN-1)/(dx**2)

        Do i=iN-2*Nx+1,iN-Nx
          !prod(i)=-2._PR*D*dt/(dy**2)*x(i)
          prod(i+Nx)=prod(i+Nx)-2._PR*D*dt/(dy**2)*x(i)
        End Do

      End If
    End If

  End Function MatVect

  Subroutine second_membre(u0,t,dx,dy,b,i1,iN,Me,Np,recep1,recep2)
    !Crée la partie du second membre du processeur
    Integer, Intent(In)                       :: i1, iN, Me, Np
    Real(PR), Dimension(i1:iN), Intent(In)    :: u0
    Real(PR), Intent(In)                      :: t, dx, dy
    Real(PR), Dimension(i1:iN), Intent(InOut) :: b
    Real(PR), Dimension(:), Intent(In)        :: recep1, recep2
    Integer                                   :: i, j, m, p

    Do i=i1,iN
      m=modulo(i-1,Nx)+1
      p=((i-1)/(Nx))+1
      b(i)=u0(i)+dt*f(case_choice,m*dx,p*dy,t)

      !Conditions de bord
      If ((m==1).or.(m==Nx)) Then
        b(i) = b(i) + D*dt*h(case_choice,m*dx,p*dy,t)/(dx**2)
      End If
      If ((p==1).or.(p==Ny)) Then
        b(i) = b(i) + D*dt*g(case_choice,m*dx,p*dy,t)/(dy**2)
      End If
    End Do

    If(alpha_rob.eq.0) Then
      If (Me.ne.0) Then
        Do i=i1,i1+Nx-1
          b(i)=recep1(i-i1+1)/(dx**2+dy**2)
        End Do
      End If
      If (Me.ne.(Np-1)) Then
        Do i=iN-Nx+1,iN
          b(i)=recep2(i-(iN-Nx))/(dx**2+dy**2)
        End Do
      End If
    Else
      If(Me.ne.0) Then
        Do i=i1,i1+Nx-1
          b(i)=b(i)+D*dt*(2._PR*beta_rob*recep1(i-i1+1+Nx)/(alpha_rob*dy)+(recep1(i-i1+1)-recep1(i-i1+1+2*Nx))/(dy)**2)
        End Do
      End If
      If(Me.ne.(Np-1)) Then
        Do i=iN-Nx+1,iN
          b(i)=b(i)+D*dt*(2._PR*beta_rob*recep2(i-iN+2*Nx)/(alpha_rob*dy)+(recep2(i-iN+3*Nx)-recep2(i-iN+Nx))/(dy)**2)
        End Do
      End If
    End If

  End Subroutine second_membre

  Function Prod_scalaire(x,y,i1,iN) Result(p)
    !Effectue un produit scalaire
    Integer, Intent(In)                     :: i1, iN
    Real(PR), Dimension(i1:iN), Intent(In)  :: x, y
    Real(PR)                                :: p
    Integer                                 :: i

    p=0
    Do i=i1,iN
       p=p+x(i)*y(i)
    End Do

  End Function Prod_scalaire


End Module matrice
