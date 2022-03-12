Module charge
  Use param
  Implicit None

CONTAINS

  Subroutine charge_(n,Np,me,i1,iN)
    !Répartit la charge totale entre les processeurs
    Integer, Intent(In)  :: n, Np, me
    Integer, Intent(Out) :: i1, iN
    Integer              :: r, q

    q = n/Np
    r = n-q*Np

    If (me<r) Then
       i1 = me*(q+1) + 1
       iN = (1+me)*(q+1)
    Else
       i1 = 1+r+me*q
       iN = i1+q-1
    End if

    If (me.ne.0) Then
      i1=i1-recouvrement
    End If
    If (me.ne.(Np-1)) Then
      iN=iN+recouvrement
    End If

    i1=(i1-1)*Nx+1
    iN=iN*Nx

  End subroutine charge_

End Module charge
