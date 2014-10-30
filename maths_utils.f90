module maths_utils
  implicit none
  private
  public :: cds1, cds2
  
contains
  
  pure integer function cds1(nPoints, dim)
    integer, intent(in) :: nPoints, dim
    
    if (dim == 1) then
       cds1 = nPoints - 1
    else
       cds1 = nPoints
    end if
  end function cds1
  
  !............................................................................
  
  pure integer function cds2(nPoints, dim)
    integer, intent(in) :: nPoints, dim
    
    if (dim == 2) then
       cds2 = nPoints - 1
    else
       cds2 = nPoints
    end if
  end function cds2
end module maths_utils
