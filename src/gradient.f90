Module gradient
  Use param
  Use matrice
  Use MPI
  Use charge
  Implicit None

CONTAINS

  Subroutine grad(alpha,beta,gamma,b,x,i1,iN,dx,dy,Me,Np)
    !Renvoie la partie de solution du processeur par la méthode du gradient conjugué
    Integer, Intent(In)                       :: i1, iN, Me, Np
    Real(PR),dimension(i1:iN), Intent(In)     :: b
    Real(PR),dimension(i1:iN), Intent(InOut)  :: x
    Real(PR),Intent(In)                       :: alpha, beta, gamma, dx, dy
    Real(PR),dimension(i1:iN)                 :: r, p, z
    Real(PR)                                  :: bet, alph, gamm, psr, somme, somme2, psz
    Integer                                   :: k
    Integer                                   :: Statinfo
    Integer, dimension(MPI_STATUS_SIZE)       :: Status

    x=1
    r=b-MatVect(alpha,beta,gamma,x,i1,iN,dx,dy,Me,Np)
    p=r
    psr=Prod_scalaire(r,r,i1,iN)
    bet=sqrt(psr)
    k=0

    !Boucle de tolérance
    Do While ((bet>epsilon).and.(k<=kmax))
       z=Matvect(alpha,beta,gamma,p,i1,iN,dx,dy,Me,Np)
       psz=Prod_scalaire(z,p,i1,iN)

       !Réduction pour avoir alph
       alph=psr/psz
       x=x+alph*p
       r=r-alph*z
       psr=Prod_scalaire(r,r,i1,iN)

       !Réduction pour avoir gamm
       gamm=psr/(bet**2)
       p=r+gamm*p
       bet=sqrt(psr)
       k=k+1
    End Do

  End Subroutine grad

  Subroutine bicgstab(alpha,beta,gamma,b,x,i1,iN,dx,dy,Me,Np)

      Integer, Intent(In)                       :: i1, iN, Me, Np
      Real(PR),dimension(i1:iN), Intent(In)     :: b
      Real(PR),dimension(i1:iN), Intent(InOut)  :: x
      Real(PR),Intent(In)                       :: alpha, beta, gamma, dx, dy
      Real(PR),dimension(i1:iN)                 :: r, p, v, rs, s, t
      Real(PR)                                  :: bet, alph, gamm, psr, somme, somme2
      Real(PR)                                  :: rho, omega, norm_r, norm_b, rho_prev
      Integer                                   :: k
      Integer                                   :: Statinfo
      Integer, dimension(MPI_STATUS_SIZE)       :: Status

      k = 0
      x = 0._PR
      r = b-MatVect(alpha,beta,gamma,x,i1,iN,dx,dy,Me,Np)
      rs = r
      rho = 1._PR
      alph = 1._PR
      omega = 1._PR
      v = 0._PR
      p = 0._PR

      norm_r = dsqrt(dot_product(r,r))                !
      norm_b = dsqrt(dot_product(b,b))

      Do While((norm_r .GT. epsilon*norm_b) .And. (k < kmax))
          rho_prev = rho
          rho = dot_product(rs,r)
          bet = (rho/rho_prev)*(alph/omega)
          p = r + bet*(p-omega*v)
          v = MatVect(alpha,beta,gamma,p,i1,iN,dx,dy,Me,Np)
          alph = rho/dot_product(rs,v)
          s = r - alph*v
          t = MatVect(alpha,beta,gamma,s,i1,iN,dx,dy,Me,Np)
          omega = dot_product(t,s)/dot_product(t,t)
          x = x + alph*p + omega*s
          r = s - omega*t
          norm_r = dsqrt(dot_product(r,r))
          norm_b = dsqrt(dot_product(b,b))

          k = k+1
      End Do

  End Subroutine bicgstab

End Module gradient
