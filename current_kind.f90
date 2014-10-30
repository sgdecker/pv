module current_kind
  use kind_types, only: cp => dp
  implicit none
  save

  real(cp), parameter :: Zero = 0._cp, One = 1._cp, Quarter = .25_cp,  &
                         Four = 4._cp, Half = .5_cp, Eighth = .125_cp,  &
                         Two = 2._cp, Third = 1._cp / 3._cp

end module current_kind
