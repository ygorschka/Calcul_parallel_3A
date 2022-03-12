Module communications
  Use param
  Use MPI
  Implicit None

CONTAINS

  Subroutine comms(u0,i1,iN,Me,Np,recep1,recep2)
    Integer, Intent(In)                       :: i1, iN, Me, Np
    Real(PR), Dimension(i1:iN), Intent(In)    :: u0
    Real(PR), Dimension(:), Intent(InOut)     :: recep1, recep2
    Integer                                   :: Statinfo
    Integer, Dimension(MPI_STATUS_size)       :: Status

    If(alpha_rob.eq.0) Then

      If (Me/=(Np-1)) Then
        Call MPI_Send(u0(iN-2*recouvrement*Nx+1:iN-(2*recouvrement-1)*Nx),Nx,MPI_REAL8,Me+1,1,MPI_COMM_WORLD,Statinfo)
      End If

      !Réception
      If (Me/=0) Then
        Call MPI_Recv(recep1,Nx,MPI_REAL8,Me-1,1,MPI_COMM_WORLD,Status,Statinfo)
      End If

      !Envoi des Nx premiers éléments au processeur d'avant
      If (Me/=0) Then
        Call MPI_Send(u0(i1+(2*recouvrement-1)*Nx:i1+2*recouvrement*Nx-1),Nx,MPI_REAL8,Me-1,0,MPI_COMM_WORLD,Statinfo)
      End If

      !Réception
      If (Me/=(Np-1)) Then
        Call MPI_Recv(recep2,Nx,MPI_REAL8,Me+1,0,MPI_COMM_WORLD,Status,Statinfo)
      End If

    Else

      If (Me/=(Np-1)) Then
        !Call MPI_Send(u0(iN-3*Nx+1:iN),3*Nx,MPI_REAL8,Me+1,1,MPI_COMM_WORLD,Statinfo)
        Call MPI_Send(u0(iN-(2*recouvrement+1)*Nx+1:iN-(2*recouvrement-2)*Nx),3*Nx,MPI_REAL8,Me+1,1,MPI_COMM_WORLD,Statinfo)
      End If

      !Réception
      If (Me/=0) Then
        Call MPI_Recv(recep1,3*Nx,MPI_REAL8,Me-1,1,MPI_COMM_WORLD,Status,Statinfo)
      End If

      !Envoi des Nx premiers éléments au processeur d'avant
      If (Me/=0) Then
        !Call MPI_Send(u0(i1:i1+3*Nx-1),3*Nx,MPI_REAL8,Me-1,0,MPI_COMM_WORLD,Statinfo)
        Call MPI_Send(u0(i1+(2*recouvrement-2)*Nx:i1+(2*recouvrement+1)*Nx-1),3*Nx,MPI_REAL8,Me-1,0,MPI_COMM_WORLD,Statinfo)
      End If

      !Réception
      If (Me/=(Np-1)) Then
        Call MPI_Recv(recep2,3*Nx,MPI_REAL8,Me+1,0,MPI_COMM_WORLD,Status,Statinfo)
      End If

    End If

  End Subroutine comms

  Function error_schwarz(u,i1,iN,Me,Np,recep1,recep2)
    Integer, Intent(In)                    :: i1, iN, Me, Np
    Real(PR), Dimension(i1:iN), Intent(In) :: u
    Real(PR), Dimension(:), Intent(In)     :: recep1, recep2
    Integer                                :: Statinfo
    Real(PR)                               :: error_schwarz
    Real(PR)                               :: error_left, error_right, error_loc
    Real(PR), Dimension(:), Allocatable    :: recep1_temp, recep2_temp

    Allocate(recep1_temp(1:Nx),recep2_temp(1:Nx))
    If(alpha_rob.eq.0) Then
      recep1_temp=recep1
      recep2_temp=recep2
    Else
      recep1_temp=recep1(Nx+1:2*Nx)
      recep2_temp=recep2(Nx+1:2*Nx)
    End If

    If(Me==0) Then
      error_left=0._PR
      error_right=maxval(abs(u(iN-Nx+1:iN)-recep2_temp))
    ElseIf(Me==(Np-1)) Then
      error_left=maxval(abs(u(i1:i1+Nx-1)-recep1_temp))
      error_right=0._PR
    Else
      error_left=maxval(abs(u(i1:i1+Nx-1)-recep1_temp))
      error_right=maxval(abs(u(iN-Nx+1:iN)-recep2_temp))
    End If
    error_loc=max(error_left,error_right)
    Call MPI_Allreduce(error_loc,error_schwarz,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,Statinfo)
    Deallocate(recep1_temp,recep2_temp)

  End Function error_schwarz

End Module communications
