module maths
  use maths_utils, only: cds1, cds2
  implicit none
  private
  public :: mean, cen_diff_stag, cen_diff, uneven_deriv, calc_mfpsi

  real, parameter, public :: Half = 0.5, One = 1.
  real, parameter         :: Quarter = 0.25

contains

  real function mean(data)
    real, dimension(:,:), intent(in) :: data

    mean = sum(data) / size(data)
  end function mean

  ! ===========================================================================

  function cen_diff_stag(data, rdx, dim)
    real, dimension(:,:), intent(in) :: data
    real,                 intent(in) :: rdx
    integer,              intent(in) :: dim
    real, dimension(cds1(size(data,1),dim),cds2(size(data,2),dim)) :: cen_diff_stag

    integer :: nPoints
    
    nPoints = size(data,dim) - 1

    select case (dim)
    case (1)
       cen_diff_stag = rdx * (data(2:nPoints+1,:) - data(1:nPoints,:))
    case (2)
       cen_diff_stag = rdx * (data(:,2:nPoints+1) - data(:,1:nPoints))
    case default
       stop "Error in dim value in cen_diff_stag!"
    end select
  end function cen_diff_stag

  ! ===========================================================================

  function cen_diff(d, dxr, dim)
    real, dimension(:,:,:,:), intent(in) :: d
    real,                     intent(in) :: dxr
    integer,                  intent(in) :: dim
    real, dimension(size(d,1),size(d,2),size(d,3),size(d,4)) :: cen_diff

    integer :: n

    n = size(d,dim)

    select case (dim)
    case (1)
       cen_diff(2:n-1,:,:,:) = Half * dxr * (d(3:n,:,:,:) - d(1:n-2,:,:,:))
       cen_diff((/1,n/),:,:,:) = dxr * (d((/2,n/),:,:,:) - d((/1,n-1/),:,:,:))
    case (2)
       cen_diff(:,2:n-1,:,:) = Half * dxr * (d(:,3:n,:,:) - d(:,1:n-2,:,:))
       cen_diff(:,(/1,n/),:,:) = dxr * (d(:,(/2,n/),:,:) - d(:,(/1,n-1/),:,:))
    case (4)
       cen_diff(:,:,:,2:n-1) = Half * dxr * (d(:,:,:,3:n) - d(:,:,:,1:n-2))
       cen_diff(:,:,:,(/1,n/)) = dxr * (d(:,:,:,(/2,n/)) - d(:,:,:,(/1,n-1/)))
    case default
       stop "Bad dim in cen_diff!"
    end select
  end function cen_diff

  ! ===========================================================================

  function uneven_deriv(d, coordVal, dim)
    real, dimension(:,:,:,:), intent(in) :: d
    real, dimension(:),       intent(in) :: coordVal
    integer,                  intent(in) :: dim
    real, dimension(size(d,1),size(d,2),size(d,3),size(d,4)) :: uneven_deriv

    real, dimension(2:size(coordVal)) :: deltaCoord
    real, dimension(3)                :: c
    integer                           :: n, i
    
    n = size(d,dim)
    if (n /= size(coordVal)) stop "Arrays don't conform in uneven_deriv!"
    
    deltaCoord = coordVal(2:) - coordVal(:n-1)
    select case (dim)
    case (3)
       do i = 2, n - 1
          c = deriv_coeff(deltaCoord(i), deltaCoord(i+1))
          uneven_deriv(:,:,i,:) = c(1) * d(:,:,i-1,:) + c(2) * d(:,:,i,:) +  &
               c(3)*d(:,:,i+1,:)
       end do
       uneven_deriv(:,:,(/1,n/),:) = (d(:,:,(/2,n/),:) -  &
            d(:,:,(/1,n-1/),:)) / spread(spread(spread(coordVal((/2,n/)) -  &
            coordVal((/1,n-1/)),1,size(d,1)),2,size(d,2)),4,size(d,4))
    case default
       stop "This dim value not implemented in uneven_deriv!"
    end select
  end function uneven_deriv

  ! ===========================================================================

  function calc_mfpsi(mfu, mfv)
    real, dimension(:,:), intent(in)             :: mfu, mfv
    real, dimension(size(mfv,1)-1,size(mfu,2)-1) :: calc_mfpsi

    integer :: nx, ny!, i, j

    nx = size(mfv,1)
    ny = size(mfu,2)
    if (size(mfv,2) /= ny-1 .or. size(mfu,1) /= nx-1) stop "Mismatch in mfpsi!"
    
    calc_mfpsi = Quarter * (mfu(:,:ny-1)+mfu(:,2:) + mfv(:nx-1,:)+mfv(2:,:))
!!$    do j = 1, ny-1
!!$       do i = 1, nx-1
!!$          calc_mfpsi(i,j) = Quarter *  &
!!$               (mfu(i,j) + mfu(i,j+1) + mfv(i,j) + mfv(i+1,j))
!!$       end do
!!$    end do
  end function calc_mfpsi

  ! ===========================================================================

  ! A function to calculate the coefficients for taking a first derivative
  ! with variable grid spacing
  pure function deriv_coeff(d1, d2)
    real, intent(in)   :: d1, d2      ! Grid spacings adj. to derivative pt.
    real, dimension(3) :: deriv_coeff
    
    deriv_coeff(1) = -d2 / (d1 * (d1+d2))
    deriv_coeff(2) = (d2**2 - d1**2) / (d1 * d2 * (d1+d2))
    deriv_coeff(3) = d1 / (d2 * (d1+d2))
  end function deriv_coeff
  
end module maths
