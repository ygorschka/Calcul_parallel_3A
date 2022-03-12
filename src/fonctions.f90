Module fonctions
  Use param
  Implicit None

CONTAINS

  Function f(case_choice,x,y,t)
    !Terme source
    Integer, Intent(In)  :: case_choice
    Real(PR), Intent(In) :: x, y, t
    Real(PR)             :: f

    !Sélection du cas
    Select case(case_choice)
    Case(1)
      f=2._PR*(y-y**2+x-x**2)
    Case(2)
      f=sin(x)+cos(y)
    Case(3)
      f=exp(-(x-(Lx/(2._PR)))**2)*exp(-(y-(Ly/(2._PR)))**2)*cos(Pi*t/(2._PR))
    Case Default
      Print *, "Choix de cas non valide"
    End Select
  End Function f

  Function g(case_choice,x,y,t)
    !Condition de bords en haut et en bas
    Integer, Intent(In)  :: case_choice
    Real(PR), Intent(In) :: x, y, t
    Real(PR)             :: g

    !Sélection du cas
    Select case(case_choice)
    Case(1)
      g=0._PR
    Case(2)
      g=sin(x)+cos(y)
    Case(3)
      g=0._PR
    Case Default
      Print *, "Choix de cas non valide"
    End Select
  End Function g

  Function h(case_choice,x,y,t)
    !Condition sur les bords latéraux
    Integer, Intent(In)  :: case_choice
    Real(PR), Intent(In) :: x, y, t
    Real(PR)             :: h

    !Sélection du cas
    Select case(case_choice)
    Case(1)
      h=0._PR
    Case(2)
      h=sin(x)+cos(y)
    Case(3)
      h=1._PR
    Case Default
      Print *, "Choix de cas non valide"
    End Select
  End Function h

End Module fonctions
