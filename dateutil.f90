module dateutil
  implicit none
  save
  private

  public :: hr_diff
  public :: operator(<), operator(>), operator(==)

  type, public :: datestamp
     integer :: yr, mth, day, hr
  end type datestamp
  
  interface operator(<)
     module procedure datestamp_lt
  end interface
  
  interface operator(>)
     module procedure datestamp_gt
  end interface

  interface operator(==)
     module procedure datestamp_eq
  end interface
  
  
contains  ! ===================================================================
  
  
  ! hr_diff returns the number of hours elapsed from d1 to d2.
  ! If the second time is before the first, -1 is returned.
  integer function hr_diff(d2, d1)
    type(datestamp), intent(in) :: d2, d1
    
    integer, parameter :: MaxHours = 999
    
    type(datestamp) :: d
    integer :: i
    
    ! Make sure the second time doesn't come before the first time
    if (d2 < d1) then
       hr_diff = -1
    else if (d2 == d1) then  ! Times match
       hr_diff = 0
    else
       d = d1
       hr_diff = 0
       do i = 1, MaxHours  ! Progress forward in time hour by hour
          hr_diff = hr_diff + 1
          d%hr = d%hr + 1
          if (d%hr > 23) then  ! New day
             d%hr = 0
             d%day = d%day + 1
             if (d%day > days_in_month(d%mth, d%yr)) then  ! New month
                d%day = 1
                d%mth = d%mth + 1
                if (d%mth > 12) then  ! New year
                   d%mth = 1
                   d%yr = d%yr + 1
                end if
             end if
          end if
          
          ! Check to see if it is the same time
          if (d == d2) exit
       end do
    end if
  end function hr_diff
  
  ! ===========================================================================
  
  ! days_in_month returns the number of days month has given the year.
  integer function days_in_month(month, year)
    integer, intent(in) :: month, year
    
    select case (month)
    case (4,6,9,11)
       days_in_month = 30
    case (1,3,5,7,8,10,12)
       days_in_month = 31
    case (2)
       if (mod(year, 4) == 0 .and.  &
            (mod(year, 100) /= 0 .or. mod(year, 400) == 0)) then
          days_in_month = 29
       else
          days_in_month = 28
       end if
    case default
       stop "Invalid input to days_in_month!"
    end select
  end function days_in_month
  
  ! ===========================================================================
  
  ! datestamp_lt returns .true. if d1 occurs before d2, .false. otherwise.
  logical function datestamp_lt(d1, d2)
    type(datestamp), intent(in) :: d1, d2
    
    if (d1%yr > d2%yr .or. d1%mth > d2%mth .and. d1%yr == d2%yr  &
         .or. d1%day > d2%day .and. d1%mth == d2%mth  .and. d1%yr == d2%yr  &
         .or. d1%hr >= d2%hr .and. d1%day == d2%day .and. d1%mth == d2%mth  &
         .and. d1%yr == d2%yr) then
       datestamp_lt = .false.
    else
       datestamp_lt = .true.
    end if
  end function datestamp_lt
  
  ! ===========================================================================

  ! datestamp_gt returns .true. if d1 occurs after d2, .false. otherwise.
   logical function datestamp_gt(d1, d2)
    type(datestamp), intent(in) :: d1, d2
    
    if (d1%yr > d2%yr .or. d1%mth > d2%mth .and. d1%yr == d2%yr  &
         .or. d1%day > d2%day .and. d1%mth == d2%mth  .and. d1%yr == d2%yr  &
         .or. d1%hr > d2%hr .and. d1%day == d2%day .and. d1%mth == d2%mth  &
         .and. d1%yr == d2%yr) then
       datestamp_gt = .true.
    else
       datestamp_gt = .false.
    end if
  end function datestamp_gt

  ! ===========================================================================
  
  ! datestamp_eq returns .true. if d1 equals d2, .false. otherwise.
  logical function datestamp_eq(d1, d2)
    type(datestamp), intent(in) :: d1, d2
    
    if (d1%yr == d2%yr .and. d1%mth == d2%mth .and. d1%day == d2%day  &
         .and. d1%hr == d2%hr) then
       datestamp_eq = .true.
    else
       datestamp_eq = .false.
    end if
  end function datestamp_eq
end module dateutil
