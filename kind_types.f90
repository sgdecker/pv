module kind_types
  implicit none
  save

  integer, parameter :: sp = selected_real_kind(5)
  integer, parameter :: dp = selected_real_kind(10)
end module kind_types
