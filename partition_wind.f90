module partition_wind
  use current_kind, only: cp, Zero, One, Quarter, Four, Half, Eighth, Two, Third
  use diagnostics,  only: RGas, Grav, Kappa
  implicit none
  private
  public :: partition, balance_lhs, balance_lhs_exp, balance_lhs_tlm,  &
            balance_lhs_adj, balance_rhs_simp, balance_rhs_simp_exp, balance_rhs_simp_tlm,  &
            balance_rhs_simp_adj, imbalance_simp_adj, imbalance_exp, pv_exp, pv_tlm, pv_adj,  &
            solve_phi_psi

contains
  
  ! Extents should be as follows:
  ! us/mu - (0:nx,ny) ; vs/mv - (nx,0:ny) ; div/mchi - (nx,ny)
  ! vort/mpsi - (nx-1,ny-1) ; chi - (0:nx+1,0:ny+1) ; psi - (0:nx,0:ny)
  subroutine partition(us, vs, div, vort, mchi, mpsi, mu, mv, phill, rdx, chi, psi)
    real(cp), dimension(0:,:), intent(in) :: us, mu
    real(cp), dimension(:,0:), intent(in) :: vs, mv
    real(cp), dimension(:,:),  intent(in) :: div, vort, mchi, mpsi
    real(cp),                  intent(in) :: phill, rdx
    real(cp), dimension(0:size(us,1),0:size(vs,2)), intent(out) :: chi
    real(cp), dimension(0:size(vs,1),0:size(us,2)), intent(out) :: psi

    real(cp),         parameter :: Eps = 1e-4_cp
    integer,          parameter :: MaxIter = 50000
    character(len=*), parameter :: Mismatch = "mismatch in partition!"

    real(cp), dimension(:,:,:),        allocatable :: ps, ch
    real(cp), dimension(size(div,1),size(div,2))   :: divForcing, divCheck
    real(cp), dimension(0:size(div,1),size(div,2)) :: uu
    real(cp), dimension(size(div,1),0:size(div,2)) :: vv
    real(cp), dimension(size(div,1)-1,size(div,2)) :: u
    real(cp), dimension(size(div,1),size(div,2)-1) :: v
    real(cp), dimension(size(vort,1),size(vort,2)) :: vortForcing, vortCheck
    real(cp) :: dx, dx2, err
    integer  :: nx, ny, prev, cur, i, j, k

    nx = size(vs,1)
    ny = size(us,2)
    dx = One / rdx
    dx2 = dx**2
    divForcing = dx2 * div / mchi**2
    vortForcing = dx2 * vort / mpsi**2
    
    ! Sanity checks
    if (any(shape(us) /= shape(mu))) stop "us mu " // Mismatch
    if (any(shape(vs) /= shape(mv))) stop "vs mv " // Mismatch
    if (any(shape(div) /= shape(mchi))) stop "div " // Mismatch
    if (any(shape(vort) /= shape(mpsi))) stop "vort " // Mismatch
    if (any(shape(mchi) /= (/ nx, ny /))) stop "mchi " // Mismatch
    if (any(shape(mpsi) /= (/ nx-1, ny-1 /))) stop "mpsi " // Mismatch
    if (any(shape(chi) /= (/ nx+2, ny+2 /))) stop "chi " // Mismatch
    if (any(shape(psi) /= (/ nx+1, ny+1 /))) stop "psi " // Mismatch

    ! Initialize
    allocate(ch(0:nx+1,0:ny+1,2), ps(0:nx,0:ny,2))
    prev = 1; cur = 2
    
    ! First guess
    ch = Zero ; ps = Zero
    
    do k = 1, MaxIter
       prev = 2 - mod(k,2)
       cur = 2 - mod(k+1,2)

       ! Compute boundary values of psi
!       ps(0,1,cur) = (One - Eighth) * ps(0,1,prev) + Eighth * (Two * phill - ps(1,1,prev))
!       print *, phill, ps(0,1,cur), ps(1,1,prev)
!       do j = 2, ny   ! Set ps(0,1) to zero
!!!!       do j = 1, ny
!          ps(0,j,cur) = ps(1,j,prev) +  &
!               (ch(1,j+1,prev) - ch(1,j,prev)) - dx * vs(1,j) / mv(1,j)
!       end do
       ps(0,2:,cur) = ps(1,2:,prev) +  &
            (ch(1,3:,prev) - ch(1,2:ny,prev)) - dx * vs(1,2:) / mv(1,2:)
!!$       do i = 1, nx
!!$          ps(i,ny,cur) = ps(i,ny-1,prev) +  &
!!$               (ch(i+1,ny,prev) - ch(i,ny,prev)) - dx * us(i,ny) / mu(i,ny)
!!$       end do
       ps(1:,ny,cur) = ps(1:,ny-1,prev) +  &
            (ch(2:,ny,prev) - ch(1:nx,ny,prev)) - dx * us(1:,ny) / mu(1:,ny)
!!$       do j = ny-1, 0, -1
!!$          ps(nx,j,cur) = ps(nx-1,j,prev) -  &
!!$               (ch(nx,j+1,prev) - ch(nx,j,prev)) + dx * vs(nx,j) / mv(nx,j)
!!$       end do
       ps(nx,:ny-1,cur) = ps(nx-1,:ny-1,prev) -  &
            (ch(nx,1:ny,prev)-ch(nx,:ny-1,prev)) + dx*vs(nx,:ny-1)/mv(nx,:ny-1)
!!$       do i = nx-1, 1, -1
!!$          ps(i,0,cur) = ps(i,1,prev) -  &
!!$               (ch(i+1,1,prev) - ch(i,1,prev)) + dx * us(i,1) / mu(i,1)
!!$       end do
       ps(1:nx-1,0,cur) = ps(1:nx-1,1,prev) -  &
            (ch(2:nx,1,prev)-ch(1:nx-1,1,prev)) + dx*us(1:nx-1,1)/mu(1:nx-1,1)
!!$       do j = 1, ny
!!$          ps(0,j,cur) = ps(0,j-1,cur) +  &
!!$               (ch(1,j,prev) - ch(0,j,prev)) - dx * us(0,j) / mu(0,j)
!!$       end do
!!$       do i = 1, nx-1
!!$          ps(i,ny,cur) = ps(i-1,ny,cur) -  &
!!$               (ch(i,ny+1,prev) - ch(i,ny,prev)) + dx * vs(i,ny) / mv(i,ny)
!!$       end do
!!$       do i = 1, nx
!!$          ps(i,0,cur) = ps(i-1,0,cur) -  &
!!$               (ch(i,1,prev) - ch(i,0,prev)) + dx * vs(i,0) / mv(i,0)
!!$       end do
!!$       do j = 1, ny-1
!!$          ps(nx,j,cur) = ps(nx,j-1,cur) +  &
!!$               (ch(nx+1,j,prev) - ch(nx,j,prev)) - dx * us(nx,j) / mu(nx,j)
!!$       end do

       ! Compute interior values of psi
       call update_sor(ps, vortForcing, prev, cur)
!!$       do j = 1, ny-2
!!$          do i = 1, nx-2
!!$             ps(i,j,cur) = ps(i,j,prev) + Quarter * Omega * (  &
!!$                  ps(i-1,j,cur) + ps(i+1,j,prev) + ps(i,j-1,cur) +  &
!!$                  ps(i,j+1,prev) - Four * ps(i,j,prev) -  &
!!$                  dx2 * vort(i,j) / mpsi(i,j)**2 )
!!$          end do
!!$          ! i = nx-1 special
!!$          ps(i,j,cur) = ps(i,j,prev) + Quarter * Omega * (  &
!!$               ps(i-1,j,cur) + ps(i+1,j,cur) + ps(i,j-1,cur) +  &
!!$               ps(i,j+1,prev) - Four * ps(i,j,prev) -  &
!!$               dx2 * vort(i,j) / mpsi(i,j)**2 )
!!$       end do
!!$       ! j = ny-1 special
!!$       do i = 1, nx-2
!!$          ps(i,j,cur) = ps(i,j,prev) + Quarter * Omega * (  &
!!$               ps(i-1,j,cur) + ps(i+1,j,prev) + ps(i,j-1,cur) +  &
!!$               ps(i,j+1,cur) - Four * ps(i,j,prev) -  &
!!$               dx2 * vort(i,j) / mpsi(i,j)**2 )
!!$       end do
!!$       ! (nx-1,ny-1) special
!!$       ps(i,j,cur) = ps(i,j,prev) + Quarter * Omega * (  &
!!$            ps(i-1,j,cur) + ps(i+1,j,cur) + ps(i,j-1,cur) +  &
!!$            ps(i,j+1,cur) - Four * ps(i,j,prev) -  &
!!$            dx2 * vort(i,j) / mpsi(i,j)**2 )

       ! Print out psi
!!$       print *, dx
!!$       print *, ps(4:5,10,cur)
!!$       print *, ch(5,10:11,cur)
!!$       print *, vs(5,10)

       ! Compute boundary values of chi
!       do j = 1, ny-1   ! Let ch(nx+1,ny) be fixed at zero
!!!!       do j = 1, ny  ! Don't fix ch
!          ch(nx+1,j,cur) = ch(nx,j,prev) +  &
!               (ps(nx,j,cur) - ps(nx,j-1,cur)) + dx * us(nx,j) / mu(nx,j)
!       end do
!sgd       ch(nx+1,1:ny-1,cur) = ch(nx,1:ny-1,prev) +  &
!sgd            (ps(nx,1:ny-1,cur)-ps(nx,:ny-2,cur)) + dx*us(nx,:ny-1)/mu(nx,:ny-1)
!!$       do i = 1, nx
!!$          ch(i,0,cur) = ch(i,1,prev) +  &
!!$               (ps(i,0,cur) - ps(i-1,0,cur)) - dx * vs(i,0) / mv(i,0)
!!$       end do
!sgd       ch(1:nx,0,cur) = ch(1:nx,1,prev) +  &
!sgd            (ps(1:,0,cur) - ps(:nx-1,0,cur)) - dx * vs(:,0) / mv(:,0)
!!$       do j = 1, ny
!!$          ch(0,j,cur) = ch(1,j,prev) -  &
!!$               (ps(0,j,cur) - ps(0,j-1,cur)) - dx * us(0,j) / mu(0,j)
!!$       end do
!sgd       ch(0,1:ny,cur) = ch(1,1:ny,prev) -  &
!sgd            (ps(0,1:,cur) - ps(0,:ny-1,cur)) - dx * us(0,:) / mu(0,:)
!!$       do i = 1, nx
!!$          ch(i,ny+1,cur) = ch(i,ny,prev) -  &
!!$               (ps(i,ny,cur) - ps(i-1,ny,cur)) + dx * vs(i,ny) / mv(i,ny)
!!$       end do
!sgd       ch(1:nx,ny+1,cur) = ch(1:nx,ny,prev) -  &
!sgd            (ps(1:,ny,cur) - ps(:nx-1,ny,cur)) + dx * vs(:,ny) / mv(:,ny)

       ! Compute interior values of chi
       call update_sor(ch, divForcing, prev, cur)

       ! Check current errors
       ! First u
!!$       do j = 1, ny
!!$          do i = 1, nx-1
!!$             u(i,j) = mu(i,j) * (ch(i+1,j,cur) - ch(i,j,cur) -  &
!!$                  (ps(i,j,cur) - ps(i,j-1,cur))) / dx
!!$          end do
!!$       end do
       u = mu(1:nx-1,:) * (ch(2:nx,1:ny,cur) - ch(1:nx-1,1:ny,cur) -  &
            (ps(1:nx-1,1:,cur) - ps(1:nx-1,:ny-1,cur))) / dx
       ! Then v
!!$       do j = 1, ny-1
!!$          do i = 1, nx
!!$             v(i,j) = mv(i,j) * (ps(i,j,cur) - ps(i-1,j,cur) +  &
!!$                  ch(i,j+1,cur) - ch(i,j,cur)) / dx
!!$          end do
!!$       end do
       v = mv(:,1:ny-1) * (ps(1:,1:ny-1,cur) - ps(:nx-1,1:ny-1,cur) +  &
            ch(1:nx,2:ny,cur) - ch(1:nx,1:ny-1,cur)) / dx
       err = sqrt(sum(abs(u-us(1:nx-1,:))**2) + sum(abs(v-vs(:,1:ny-1))**2))
!       if (mod(k,10) == 0) print *, "Wind:", k, err
       if(err < Eps) exit

       ! Check for errors in vorticity and divergence
!!$       do j = 1, ny-1
!!$          do i = 1, nx-1
!!$             vortCheck(i,j) = mpsi(i,j)**2 * (  &
!!$                  ps(i-1,j,cur) + ps(i+1,j,cur) + ps(i,j-1,cur) +  &
!!$                  ps(i,j+1,cur) - Four * ps(i,j,cur) ) / dx2
!!$          end do
!!$       end do
       vortCheck = mpsi**2 * (  &
            ps(:nx-2,1:ny-1,cur) + ps(2:,1:ny-1,cur) + ps(1:nx-1,:ny-2,cur) + &
            ps(1:nx-1,2:,cur) - Four * ps(1:nx-1,1:ny-1,cur) ) / dx2
!       if (mod(k,10) == 0) print *, "Vort:", k, maxval(abs(vortCheck-vort))
!!$       do j = 1, ny
!!$          do i = 1, nx
!!$             divCheck(i,j) = mchi(i,j)**2 * (  &
!!$                  ch(i-1,j,cur) + ch(i+1,j,cur) + ch(i,j-1,cur) +  &
!!$                  ch(i,j+1,cur) - Four * ch(i,j,cur) ) / dx2
!!$          end do
!!$       end do
       divCheck = mchi**2 * (  &
            ch(:nx-1,1:ny,cur) + ch(2:,1:ny,cur) + ch(1:nx,:ny-1,cur) +  &
            ch(1:nx,2:,cur) - Four * ch(1:nx,1:ny,cur) ) / dx2
!       if (mod(k,10) == 0) print *, "Div: ", k, maxval(abs(divCheck-div))
       if (mod(k,250) == 0) write (*, "(a,i6,a,3es12.5)")  &
            "Iteration", k, ". Errors in Wnd/Vor/Div: ", err,  &
            maxval(abs(vortCheck-vort)), maxval(abs(divCheck-div))

       ! Finally, calculate divergence from psi and vorticity from chi
!!$       do j = 1, ny
!!$          do i = 1, nx-1
!!$             u(i,j) = -mu(i,j) * (ps(i,j,cur) - ps(i,j-1,cur)) / dx
!!$          end do
!!$       end do
!!$       do j = 1, ny-1
!!$          do i = 1, nx
!!$             v(i,j) = mv(i,j) * (ps(i,j,cur) - ps(i-1,j,cur)) / dx
!!$          end do
!!$       end do
!!$
!!$       divCheck(1,1) = mchi(1,1)**2 * rdx * (  &
!!$            u(1,1)/mu(1,1) - us(0,1)/mu(0,1) +  &
!!$            v(1,1)/mv(1,1) - vs(1,0)/mv(1,0) )
!!$       do i = 2, nx-1
!!$          divCheck(i,1) = mchi(i,1)**2 * rdx * (  &
!!$               u(i,1)/mu(i,1) - u(i-1,1)/mu(i-1,1) +  &
!!$               v(i,1)/mv(i,1) - vs(i,0)/mv(i,0) )
!!$       end do
!!$       divCheck(i,1) = mchi(i,1)**2 * rdx * (  &
!!$            us(i,1)/mu(i,1) - u(i-1,1)/mu(i-1,1) +  &
!!$            v(i,1)/mv(i,1) - vs(i,0)/mv(i,0) )
!!$       do j = 2, ny-1
!!$          divCheck(1,j) = mchi(1,j)**2 * rdx * (  &
!!$               u(1,j)/mu(1,j) - us(0,j)/mu(0,j) +  &
!!$               v(1,j)/mv(1,j) - v(1,j-1)/mv(1,j-1) )
!!$          do i = 2, nx-1
!!$             divCheck(i,j) = mchi(i,j)**2 * rdx * (  &
!!$                  u(i,j)/mu(i,j) - u(i-1,j)/mu(i-1,j) +  &
!!$                  v(i,j)/mv(i,j) - v(i,j-1)/mv(i,j-1) )
!!$          end do
!!$          divCheck(i,j) = mchi(i,j)**2 * rdx * (  &
!!$               us(i,j)/mu(i,j) - u(i-1,j)/mu(i-1,j) +  &
!!$               v(i,j)/mv(i,j) - v(i,j-1)/mv(i,j-1) )          
!!$       end do
!!$       divCheck(1,j) = mchi(1,j)**2 * rdx * (  &
!!$            u(1,j)/mu(1,j) - us(0,j)/mu(0,j) +  &
!!$            vs(1,j)/mv(1,j) - v(1,j-1)/mv(1,j-1) )
!!$       do i = 2, nx-1
!!$          divCheck(i,j) = mchi(i,j)**2 * rdx * (  &
!!$               u(i,j)/mu(i,j) - u(i-1,j)/mu(i-1,j) +  &
!!$               vs(i,j)/mv(i,j) - v(i,j-1)/mv(i,j-1) )
!!$       end do
!!$       divCheck(i,j) = mchi(i,j)**2 * rdx * (  &
!!$            us(i,j)/mu(i,j) - u(i-1,j)/mu(i-1,j) +  &
!!$            vs(i,j)/mv(i,j) - v(i,j-1)/mv(i,j-1) )       

!!$       do j = 1, ny
!!$          do i = 0, nx
!!$             uu(i,j) = -mu(i,j) * (ps(i,j,cur) - ps(i,j-1,cur)) / dx
!!$          end do
!!$       end do
       uu = -mu * (ps(:,1:,cur) - ps(:,:ny-1,cur)) / dx
!!$       do j = 0, ny
!!$          do i = 1, nx
!!$             vv(i,j) = mv(i,j) * (ps(i,j,cur) - ps(i-1,j,cur)) / dx
!!$          end do
!!$       end do
       vv = mv * (ps(1:,:,cur) - ps(:nx-1,:,cur)) / dx
!!$       do j = 1, ny
!!$          do i = 1, nx
!!$             divCheck(i,j) = mchi(i,j)**2 * rdx * (  &
!!$                  uu(i,j)/mu(i,j) - uu(i-1,j)/mu(i-1,j) +  &
!!$                  vv(i,j)/mv(i,j) - vv(i,j-1)/mv(i,j-1) )
!!$          end do
!!$       end do
       divCheck = mchi**2 * rdx * (  &
            uu(1:,:)/mu(1:,:) - uu(:nx-1,:)/mu(:nx-1,:) +  &
            vv(:,1:)/mv(:,1:) - vv(:,:ny-1)/mv(:,:ny-1) )

!!$       do j = 1, ny
!!$          do i = 0, nx
!!$             uu(i,j) = mu(i,j) * (ch(i+1,j,cur) - ch(i,j,cur)) / dx
!!$          end do
!!$       end do
       uu = mu * (ch(1:,1:ny,cur) - ch(:nx,1:ny,cur)) / dx
!!$       do j = 0, ny
!!$          do i = 1, nx
!!$             vv(i,j) = mv(i,j) * (ch(i,j+1,cur) - ch(i,j,cur)) / dx
!!$          end do
!!$       end do
       vv = mv * (ch(1:nx,1:,cur) - ch(1:nx,:ny,cur)) / dx
!!$       do j = 1, ny-1
!!$          do i = 1, nx-1
!!$             vortCheck(i,j) = mpsi(i,j)**2 * rdx * &
!!$                  (vv(i+1,j)/mv(i+1,j) - vv(i,j)/mv(i,j) -  &
!!$                  (uu(i,j+1)/mu(i,j+1) - uu(i,j)/mu(i,j)))
!!$          end do
!!$       end do
       vortCheck = mpsi**2 * rdx *  &
            (vv(2:,1:ny-1)/mv(2:,1:ny-1)-vv(:nx-1,1:ny-1)/mv(:nx-1,1:ny-1) -  &
            (uu(1:nx-1,2:)/mu(1:nx-1,2:)-uu(1:nx-1,:ny-1)/mu(1:nx-1,:ny-1)))
            
       if (mod(k,250) == 0) write (*, "(a,2es12.5)")  &
            "          Div from Psi / Vor from Chi: ", maxval(abs(divCheck)), &
            maxval(abs(vortCheck))
!!$       print *, vortCheck(1,:)
!!$       print *, vortCheck(5,:)
!!$       print *, nx-1, vortCheck(nx-1,:)
!!$       print *
    end do
    if (k > MaxIter) stop "Too many iterations!"

!!$!    do j = 2, 0, -1
!!$!       print *, ps(0:2,j,cur)
!!$!    end do
!!$    do j = 2, 1, -1
!!$       print *, uu(0:2,j)
!!$    end do
!!$    do j = 2, 0, -1
!!$       print *, vv(1:2,j)
!!$    end do
!!$!    do j = 2, 1, -1
!!$!       print *, vort(1:2,j)
!!$!    end do
!!$!    print *
!!$    do j = 2, 1, -1
!!$       print *, mpsi(1:2,j)
!!$    end do
    

    chi = ch(:,:,cur)
    psi = ps(:,:,cur)
  end subroutine partition

  ! ===========================================================================
  
  subroutine update_sor(data, forcing, prev, cur)
    real(cp), dimension(0:,0:,:), intent(inout) :: data
    real(cp), dimension(:,:),     intent(in)    :: forcing
    integer,                      intent(in)    :: prev, cur

    real(cp), parameter :: Omega = 1.8_cp 

    integer :: nx, ny, j, i

    nx = size(forcing,1)
    ny = size(forcing,2)

    ! Sanity checks
    if (all(shape(data) /= (/ nx+2, ny+2, 2/))) stop "Mismatch in update_sor!"
    if (prev + cur /= 3 .or. cur /= 1 .and. cur /= 2) stop "prev/cur error!"

    do j = 1, ny-1
       do i = 1, nx-1
          data(i,j,cur) = data(i,j,prev) + Quarter * Omega * (  &
               data(i-1,j,cur) + data(i+1,j,prev) + data(i,j-1,cur) +  &
               data(i,j+1,prev) - Four * data(i,j,prev) - forcing(i,j) )
       end do
       ! i = nx special
       data(i,j,cur) = data(i,j,prev) + Quarter * Omega * (  &
            data(i-1,j,cur) + data(i+1,j,cur) + data(i,j-1,cur) +  &
            data(i,j+1,prev) - Four * data(i,j,prev) - forcing(i,j) )
    end do
    
    ! j = ny special
    do i = 1, nx-1
       data(i,j,cur) = data(i,j,prev) + Quarter * Omega * (  &
            data(i-1,j,cur) + data(i+1,j,prev) + data(i,j-1,cur) +  &
            data(i,j+1,cur) - Four * data(i,j,prev) - forcing(i,j) )
    end do
    ! (nx,ny) special
    data(i,j,cur) = data(i,j,prev) + Quarter * Omega * (  &
         data(i-1,j,cur) + data(i+1,j,cur) + data(i,j-1,cur) +  &
         data(i,j+1,cur) - Four * data(i,j,prev) - forcing(i,j) )
  end subroutine update_sor

  ! ===========================================================================

  subroutine balance_lhs(psi, f, mCap, m, mu, mv, rdx, lhs)
    real(cp), dimension(0:,0:,:), intent(in) :: psi          ! (0:nx,0:ny)
    real(cp), dimension(:,:),     intent(in) :: f, mCap, m   ! (1:nx,1:ny)
    real(cp), dimension(0:,:),    intent(in) :: mu
    real(cp), dimension(:,0:),    intent(in) :: mv           ! (1:nx,0:ny)
    real(cp),                     intent(in) :: rdx
    real, dimension(size(f,1),size(f,2),size(psi,3)), intent(out) :: lhs

    character(len=*), parameter :: Mismatch = "mismatch in balance_lhs!"
    
    real(cp), dimension(:,:), allocatable :: mdxp, mdyp, c, fmvdxp, fmudyp, t1,  &
                                             mdx1, mdx2, t2, mdy1, mdy2, t3
    integer :: nx, ny, nz, k, i, j

    nx = size(psi,1) - 1
    ny = size(psi,2) - 1
    nz = size(psi,3)

    ! Sanity checks
    if (any(shape(f)  /= (/ nx,   ny   /) )) stop "f "  // Mismatch
    if (any(shape(m)  /= (/ nx,   ny   /) )) stop "m "  // Mismatch
    if (any(shape(mu) /= (/ nx+1, ny   /) )) stop "mu " // Mismatch
    if (any(shape(mv) /= (/ nx,   ny+1 /) )) stop "mv " // Mismatch

    allocate(mdxp(nx,0:ny), mdyp(0:nx,ny), c(nx-1,ny-1))
    allocate(fmvdxp(nx-1,2:ny-1), fmudyp(2:nx-1,ny-1))
    allocate(t1(2:nx-1,2:ny-1), t2(2:nx-1,2:ny-1), t3(2:nx-1,2:ny-1))
    allocate(mdx1(nx-1,2:ny-1), mdx2(nx-1,2:ny-1))
    allocate(mdy1(2:nx-1,ny-1), mdy2(2:nx-1,ny-1))

    lhs = Zero
    do k = 1, nz
       mdxp = mv * xdiff(psi(:,:,k), rdx)
       mdyp = mu * ydiff(psi(:,:,k), rdx)
       fmvdxp = xavg(f(:,2:ny-1)) * xyavg(mdxp(:,1:ny-1)) / mu(1:nx-1,2:ny-1)
       fmudyp = yavg(f(2:nx-1,:)) * xyavg(mdyp(1:nx-1,:)) / mv(2:nx-1,1:ny-1)
       t1 = m(2:nx-1,2:ny-1)**2 * (xdiff(fmvdxp, rdx) + ydiff(fmudyp, rdx))
       
!!$       mdx1 = mu(1:nx-1,2:ny-1) * mdyp(1:nx-1,2:ny-1) * xdiff(xavg(mdyp(:,2:ny-1)) / mCap(:,2:ny-1), rdx)
       forall (i = 1:nx-1, j = 2:ny-1)
          mdx1(i,j) = Half * mu(i,j)**2 * rdx**3  &
               * ( mu(i,j) * (One/mCap(i+1,j) - One/mCap(i,j)) * psi(i,j,k)**2  &
               - Two * mu(i,j) * (One/mCap(i+1,j) - One/mCap(i,j)) * psi(i,j-1,k)  &
                                                                          * psi(i,j,k)  &
               + mu(i,j) * (One/mCap(i+1,j) - One/mCap(i,j)) * psi(i,j-1,k)**2  &
               + mu(i+1,j)/mCap(i+1,j) * psi(i,j,k) * psi(i+1,j,k)  &
               - mu(i+1,j)/mCap(i+1,j) * psi(i,j-1,k) * psi(i+1,j,k)  &
               - mu(i+1,j)/mCap(i+1,j) * psi(i,j,k) * psi(i+1,j-1,k)  &
               + mu(i+1,j)/mCap(i+1,j) * psi(i,j-1,k) * psi(i+1,j-1,k)  &
               - mu(i-1,j)/mCap(i,j) * psi(i,j,k) * psi(i-1,j,k)  &
               + mu(i-1,j)/mCap(i,j) * psi(i,j-1,k) * psi(i-1,j,k)  &
               + mu(i-1,j)/mCap(i,j) * psi(i,j,k) * psi(i-1,j-1,k)  &
               - mu(i-1,j)/mCap(i,j) * psi(i,j-1,k) * psi(i-1,j-1,k) )
       end forall
!!$       mdx2 = mu(1:nx-1,2:ny-1) * xyavg(mdxp(:,1:ny-1)) * ydiff(yavg(mdyp(1:nx-1,:))  &
!!$            / xyavg(mCap), rdx)
       forall (i = 1:nx-1, j = 1:ny-1)
          c(i,j) = One / (mCap(i,j) + mCap(i+1,j) + mCap(i,j+1) + mCap(i+1,j+1))
       end forall
       forall (i = 1:nx-1, j = 2:ny-1)
          mdx2(i,j) = Half * mu(i,j) * rdx**3  &
               * ( mv(i,j-1) * psi(i,j-1,k) - mv(i,j-1) * psi(i-1,j-1,k)  &
               + mv(i+1,j-1) * psi(i+1,j-1,k) - mv(i+1,j-1) * psi(i,j-1,k)  &
               + mv(i,j) * psi(i,j,k) - mv(i,j) * psi(i-1,j,k)  &
               + mv(i+1,j) * psi(i+1,j,k) - mv(i+1,j) * psi(i,j,k) )  &
               * ( c(i,j) * mu(i,j+1) * psi(i,j+1,k) - c(i,j) * mu(i,j+1) * psi(i,j,k)  &
               + c(i,j) * mu(i,j) * psi(i,j,k) - c(i,j) * mu(i,j) * psi(i,j-1,k)  &
               - c(i,j-1) * mu(i,j) * psi(i,j,k) + c(i,j-1) * mu(i,j) * psi(i,j-1,k)  &
               - c(i,j-1) * mu(i,j-1) * psi(i,j-1,k)  &
               + c(i,j-1) * mu(i,j-1) * psi(i,j-2,k) )
       end forall

!       t2 = m(2:nx-1,2:ny-1) * xdiff(mdx1 - mdx2, rdx)
       forall (i = 2:nx-1, j = 2:ny-1)
          t2(i,j) = Half * rdx**4 * m(i,j) * (  &
               mu(i,j) * (mu(i,j) * (psi(i,j,k) - psi(i,j-1,k)) * (  &
               mu(i,j) * (One/mCap(i+1,j)-One/mCap(i,j)) * (psi(i,j,k) - psi(i,j-1,k)) +&
               mu(i+1,j)/mCap(i+1,j) * (psi(i+1,j,k) - psi(i+1,j-1,k)) -  &
               mu(i-1,j)/mCap(i,j) * (psi(i-1,j,k) - psi(i-1,j-1,k))) - (  &
               mv(i,j-1) * (psi(i,j-1,k) - psi(i-1,j-1,k)) +  &
               mv(i+1,j-1) * (psi(i+1,j-1,k) - psi(i,j-1,k)) +  &
               mv(i,j) * (psi(i,j,k) - psi(i-1,j,k)) +  &
               mv(i+1,j) * (psi(i+1,j,k) - psi(i,j,k))) * (  &
               c(i,j) * (mu(i,j+1) * (psi(i,j+1,k) - psi(i,j,k)) +  &
               mu(i,j) * (psi(i,j,k) - psi(i,j-1,k))) -  &
               c(i,j-1) * (mu(i,j) * (psi(i,j,k) - psi(i,j-1,k)) +  &
               mu(i,j-1) * (psi(i,j-1,k) - psi(i,j-2,k))))) -  &
               mu(i-1,j) * (mu(i-1,j) * (psi(i-1,j,k) - psi(i-1,j-1,k)) * (  &
               mu(i-1,j)*(One/mCap(i,j)-One/mCap(i-1,j))*(psi(i-1,j,k)-psi(i-1,j-1,k)) +&
               mu(i,j)/mCap(i,j) * (psi(i,j,k) - psi(i,j-1,k)) -  &
               mu(i-2,j)/mCap(i-1,j) * (psi(i-2,j,k) - psi(i-2,j-1,k))) - (  &
               mv(i-1,j-1) * (psi(i-1,j-1,k) - psi(i-2,j-1,k)) +  &
               mv(i,j-1) * (psi(i,j-1,k) - psi(i-1,j-1,k)) +  &
               mv(i-1,j) * (psi(i-1,j,k) - psi(i-2,j,k)) +  &
               mv(i,j) * (psi(i,j,k) - psi(i-1,j,k))) * (  &
               c(i-1,j) * (mu(i-1,j+1) * (psi(i-1,j+1,k) - psi(i-1,j,k)) +  &
               mu(i-1,j) * (psi(i-1,j,k) - psi(i-1,j-1,k))) -  &
               c(i-1,j-1) * (mu(i-1,j) * (psi(i-1,j,k) - psi(i-1,j-1,k)) +  &
               mu(i-1,j-1) * (psi(i-1,j-1,k) - psi(i-1,j-2,k))))))
       end forall

       mdy1 = mv(2:nx-1,1:ny-1) * mdxp(2:nx-1,1:ny-1) * ydiff(yavg(mdxp(2:nx-1,:))  &
            / mCap(2:nx-1,:), rdx)
       mdy2 = mv(2:nx-1,1:ny-1) * xyavg(mdyp(1:nx-1,:)) * xdiff(xavg(mdxp(:,1:ny-1))  &
            / xyavg(mCap), rdx)
       t3 = m(2:nx-1,2:ny-1) * ydiff(mdy1 - mdy2, rdx)
       lhs(2:nx-1,2:ny-1,k) = t1 - t2 - t3
    end do
    deallocate(mdy2, mdy1, mdx2, mdx1, t3, t2, t1, fmudyp, fmvdxp, mdyp, mdxp)
  end subroutine balance_lhs

  ! ===========================================================================

  subroutine balance_lhs_exp(psi, f, mCap, m, mu, mv, rdx, lhs)
    real(cp), dimension(0:,0:,:), intent(in) :: psi          ! (0:nx,0:ny)
    real(cp), dimension(:,:),     intent(in) :: f, mCap, m   ! (1:nx,1:ny)
    real(cp), dimension(0:,:),    intent(in) :: mu
    real(cp), dimension(:,0:),    intent(in) :: mv           ! (1:nx,0:ny)
    real(cp),                     intent(in) :: rdx
    real, dimension(size(f,1),size(f,2),size(psi,3)), intent(out) :: lhs

    character(len=*), parameter :: Mismatch = "mismatch in balance_lhs!"
    
    real(cp), dimension(:,:,:), allocatable :: t1, t2, t3
    real(cp), dimension(:,:),   allocatable :: c    
    integer :: nx, ny, nz, i, j, k
    
    nx = size(psi,1) - 1
    ny = size(psi,2) - 1
    nz = size(psi,3)

    ! Sanity checks
    if (any(shape(f)  /= (/ nx,   ny   /) )) stop "f "  // Mismatch
    if (any(shape(m)  /= (/ nx,   ny   /) )) stop "m "  // Mismatch
    if (any(shape(mu) /= (/ nx+1, ny   /) )) stop "mu " // Mismatch
    if (any(shape(mv) /= (/ nx,   ny+1 /) )) stop "mv " // Mismatch

    allocate(t1(2:nx-1,2:ny-1,nz), t2(2:nx-1,2:ny-1,nz), t3(2:nx-1,2:ny-1,nz))
    allocate(c(nx-1,ny-1))
    
    forall (i = 1:nx-1, j = 1:ny-1)
       c(i,j) = One / (mCap(i,j) + mCap(i+1,j) + mCap(i,j+1) + mCap(i+1,j+1))
    end forall

    lhs = Zero
    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       t1(i,j,k) = Eighth * rdx**2 * m(i,j)**2 * (                           &
          (f(i,  j) + f(i+1,j)) * (                                          &
            mv(i,  j-1) * (psi(i,  j-1,k) - psi(i-1,j-1,k)) +                &
            mv(i+1,j-1) * (psi(i+1,j-1,k) - psi(i,  j-1,k)) +                &
            mv(i,  j  ) * (psi(i,  j,  k) - psi(i-1,j,  k)) +                &
            mv(i+1,j  ) * (psi(i+1,j,  k) - psi(i,  j,  k)) ) / mu(i,  j) -  &
          (f(i-1,j) + f(i,j  )) * (                                          &
            mv(i-1,j-1) * (psi(i-1,j-1,k) - psi(i-2,j-1,k)) +                &
            mv(i,  j-1) * (psi(i,  j-1,k) - psi(i-1,j-1,k)) +                &
            mv(i-1,j  ) * (psi(i-1,j,  k) - psi(i-2,j,  k)) +                &
            mv(i,  j  ) * (psi(i,  j,  k) - psi(i-1,j,  k)) ) / mu(i-1,j) +  &
          (f(i,  j) + f(i,j+1)) * (                                          &
            mu(i-1,j  ) * (psi(i-1,j,  k) - psi(i-1,j-1,k)) +                &
            mu(i,  j  ) * (psi(i,  j,  k) - psi(i,  j-1,k)) +                &
            mu(i-1,j+1) * (psi(i-1,j+1,k) - psi(i-1,j,  k)) +                &
            mu(i,  j+1) * (psi(i,  j+1,k) - psi(i,  j,  k)) ) / mv(i,  j) -  &
          (f(i,j-1) + f(i,j  )) * (                                          &
            mu(i-1,j-1) * (psi(i-1,j-1,k) - psi(i-1,j-2,k)) +                &
            mu(i,  j-1) * (psi(i,  j-1,k) - psi(i,  j-2,k)) +                &
            mu(i-1,j  ) * (psi(i-1,j,  k) - psi(i-1,j-1,k)) +                &
            mu(i,  j  ) * (psi(i,  j,  k) - psi(i,  j-1,k)) ) / mv(i,j-1) )

       t2(i,j,k) = Half * rdx**4 * m(i,j) * (  &
            mu(i,j) * (mu(i,j) * (psi(i,j,k) - psi(i,j-1,k))  &
            * (mu(i,j) * (One/mCap(i+1,j)-One/mCap(i,j)) * (psi(i,j,k) - psi(i,j-1,k))  &
            + mu(i+1,j)/mCap(i+1,j) * (psi(i+1,j,k) - psi(i+1,j-1,k))  &
            - mu(i-1,j)/mCap(i,j) * (psi(i-1,j,k) - psi(i-1,j-1,k)))  &
            - (mv(i,j-1) * (psi(i,j-1,k) - psi(i-1,j-1,k))  &
            + mv(i+1,j-1) * (psi(i+1,j-1,k) - psi(i,j-1,k))  &
            + mv(i,j) * (psi(i,j,k) - psi(i-1,j,k)) + mv(i+1,j) * (psi(i+1,j,k) - psi(i,j,k)))  &
            * (c(i,j) * (mu(i,j+1) * (psi(i,j+1,k) - psi(i,j,k)) + mu(i,j)  &
            * (psi(i,j,k) - psi(i,j-1,k))) - c(i,j-1) * (mu(i,j) * (psi(i,j,k) - psi(i,j-1,k))  &
            + mu(i,j-1) * (psi(i,j-1,k) - psi(i,j-2,k)))))  &
            - mu(i-1,j) * (mu(i-1,j) * (psi(i-1,j,k) - psi(i-1,j-1,k))  &
            * (mu(i-1,j) * (One/mCap(i,j)-One/mCap(i-1,j)) * (psi(i-1,j,k) - psi(i-1,j-1,k))  &
            + mu(i,j)/mCap(i,j) * (psi(i,j,k) - psi(i,j-1,k))  &
            - mu(i-2,j)/mCap(i-1,j) * (psi(i-2,j,k) - psi(i-2,j-1,k)))  &
            - (mv(i-1,j-1) * (psi(i-1,j-1,k) - psi(i-2,j-1,k))  &
            + mv(i,j-1) * (psi(i,j-1,k) - psi(i-1,j-1,k))  &
            + mv(i-1,j) * (psi(i-1,j,k) - psi(i-2,j,k)) + mv(i,j) * (psi(i,j,k) - psi(i-1,j,k)))  &
            * (c(i-1,j) * (mu(i-1,j+1) * (psi(i-1,j+1,k) - psi(i-1,j,k))  &
            + mu(i-1,j) * (psi(i-1,j,k) - psi(i-1,j-1,k)))  &
            - c(i-1,j-1) * (mu(i-1,j) * (psi(i-1,j,k) - psi(i-1,j-1,k))  &
            + mu(i-1,j-1) * (psi(i-1,j-1,k) - psi(i-1,j-2,k))))))

       t3(i,j,k) = Half * rdx**4 * m(i,j) * (  &
            mv(i,j) * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) * (One/mCap(i,j+1)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i,j+1) * (psi(i,j+1,k)-psi(i-1,j+1,k)))  &
            - One/mCap(i,j) * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))  &
            + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)))) - (mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))  &
            + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k))  &
            + mu(i,j+1) * (psi(i,j+1,k)-psi(i,j,k))) * (c(i,j)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i+1,j) * (psi(i+1,j,k)-psi(i,j,k)))  &
            - c(i-1,j) * (mv(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))  &
            + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)))))  &
            - mv(i,j-1) * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)) * (One/mCap(i,j)  &
            * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)))  &
            - One/mCap(i,j-1) * (mv(i,j-2) * (psi(i,j-2,k)-psi(i-1,j-2,k))  &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))))  &
            - (mu(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-1,j-2,k))  &
            + mu(i,j-1) * (psi(i,j-1,k)-psi(i,j-2,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)) &
            + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k))) * (c(i,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k)))  &
            - c(i-1,j-1) * (mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))  &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))))))
    end forall
    lhs(2:nx-1,2:ny-1,:) = t1 - t2 - t3
    deallocate(c, t3, t2, t1)
  end subroutine balance_lhs_exp
  
  ! ===========================================================================

  subroutine balance_lhs_tlm(dps, psi, f, mCap, m, mu, mv, rdx, dlhs)
    real(cp), dimension(0:,0:,:), intent(in) :: dps, psi     ! (0:nx,0:ny)
    real(cp), dimension(:,:),     intent(in) :: f, mCap, m   ! (1:nx,1:ny)
    real(cp), dimension(0:,:),    intent(in) :: mu
    real(cp), dimension(:,0:),    intent(in) :: mv           ! (1:nx,0:ny)
    real(cp),                     intent(in) :: rdx
    real, dimension(size(f,1),size(f,2),size(dps,3)), intent(out) :: dlhs

    character(len=*), parameter :: Mismatch = "mismatch in balance_lhs_tlm!"
    
    real(cp), dimension(:,:,:), allocatable :: dpsi, t1, t2, t3
    real(cp), dimension(:,:),   allocatable :: c
    integer :: nx, ny, nz, i, j, k
    
    nx = size(dps,1) - 1
    ny = size(dps,2) - 1
    nz = size(dps,3)

    ! Sanity checks
    if (any(shape(f)  /= (/ nx,   ny   /) )) stop "f "  // Mismatch
    if (any(shape(m)  /= (/ nx,   ny   /) )) stop "m "  // Mismatch
    if (any(shape(mu) /= (/ nx+1, ny   /) )) stop "mu " // Mismatch
    if (any(shape(mv) /= (/ nx,   ny+1 /) )) stop "mv " // Mismatch

    allocate(t1(2:nx-1,2:ny-1,nz), t2(2:nx-1,2:ny-1,nz), t3(2:nx-1,2:ny-1,nz))
    allocate(dpsi(0:nx,0:ny,nz))
    allocate(c(nx-1,ny-1))

    dpsi = dps
    dpsi((/0,nx/),:,:) = Zero
    dpsi(1:nx-1,(/0,ny/),:) = Zero
    dlhs = Zero  ! Really just setting the i/j boundaries.

    forall (i = 1:nx-1, j = 1:ny-1)
       c(i,j) = One / (mCap(i,j) + mCap(i+1,j) + mCap(i,j+1) + mCap(i+1,j+1))
    end forall
    
    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       t1(i,j,k) = m(i,j)**2 * Eighth * rdx**2 * (  &
            dpsi(i-1,j-2,k) * mu(i-1,j-1) / mv(i,j-1) * (f(i,j-1) + f(i,j)) + &
            dpsi(i,  j-2,k) * mu(i,  j-1) / mv(i,j-1) * (f(i,j-1) + f(i,j)) + &
            dpsi(i-2,j-1,k) * mv(i-1,j-1) / mu(i-1,j) * (f(i-1,j) + f(i,j)) - &
            dpsi(i-1,j-1,k) * ( mv(i,j-1) / mu(i,  j) * (f(i,j) + f(i+1,j)) + &
               (f(i-1,j) + f(i,j)) / mu(i-1,j) * (mv(i-1,j-1) - mv(i,j-1)) +  &
                                mu(i-1,j) / mv(i,  j) * (f(i,j) + f(i,j+1)) + &
               (f(i,j-1) + f(i,j)) / mv(i,j-1) * (mu(i-1,j-1) - mu(i-1,j))) + &
            dpsi(i,  j-1,k) * (-mv(i,j-1) / mu(i-1,j) * (f(i-1,j) + f(i,j)) + &
               (f(i,j) + f(i+1,j)) / mu(i,  j) * (mv(i,j-1) - mv(i+1,j-1)) -  &
                                mu(i,  j) / mv(i,  j) * (f(i,j) + f(i,j+1)) + &
               (f(i,j-1) + f(i,j)) / mv(i,j-1) * (mu(i,  j) - mu(i,  j-1))) + &
            dpsi(i+1,j-1,k) * mv(i+1,j-1) / mu(i,  j) * (f(i,j) + f(i+1,j)) + &
            dpsi(i-2,j,  k) * mv(i-1,  j) / mu(i-1,j) * (f(i-1,j) + f(i,j)) - &
            dpsi(i-1,j,  k) * ( mv(i,  j) / mu(i,  j) * (f(i,j) + f(i+1,j)) + &
               (f(i-1,j) + f(i,j)) / mu(i-1,j) * (mv(i-1,j) - mv(i,    j)) +  &
                                mu(i-1,j) / mv(i,j-1) * (f(i,j-1) + f(i,j)) + &
               (f(i,j) + f(i,j+1)) / mv(i,  j) * (mu(i-1,j+1) - mu(i-1,j))) + &
            dpsi(i,  j,  k) * (-mv(i,  j) / mu(i-1,j) * (f(i-1,j) + f(i,j)) + &
               (f(i,j) + f(i+1,j)) / mu(i,  j) * (mv(i,    j) - mv(i+1,j)) -  &
                                mu(i,  j) / mv(i,j-1) * (f(i,j-1) + f(i,j)) + &
               (f(i,j) + f(i,j+1)) / mv(i,  j) * (mu(i,    j) - mu(i,j+1))) + &
            dpsi(i+1,j,  k) * mv(i+1,  j) / mu(i,  j) * (f(i,j) + f(i+1,j)) + &
            dpsi(i-1,j+1,k) * mu(i-1,j+1) / mv(i,  j) * (f(i,j) + f(i,j+1)) + &
            dpsi(i,  j+1,k) * mu(i,  j+1) / mv(i,  j) * (f(i,j) + f(i,j+1)))
       
       t2(i,j,k) = m(i,j) * Half * rdx**4 * (  &
            dpsi(i-1,j-2,k) * c(i-1,j-1) * mu(i-1,j-1) * mu(i-1,j) * (  &
            (psi(i,j,k)-psi(i-1,j,k)) * mv(i,j) + (psi(i-1,j,k)-psi(i-2,j,k)) * mv(i-1,j) +  &
            (psi(i,j-1,k)-psi(i-1,j-1,k))*mv(i,j-1) + (psi(i-1,j-1,k)-psi(i-2,j-1,k))*mv(i-1,j-1)))
       t2(i,j,k) = t2(i,j,k) - m(i,j) * Half * rdx**4 * (  &
            dpsi(i,j-2,k) * c(i,j-1) * mu(i,j-1) * mu(i,j) * (  &
            (psi(i+1,j,k)-psi(i,j,k)) * mv(i+1,j) + (psi(i,j,k)-psi(i-1,j,k)) * mv(i,j) +  &
            (psi(i+1,j-1,k)-psi(i,j-1,k))*mv(i+1,j-1) + (psi(i,j-1,k)-psi(i-1,j-1,k))*mv(i,j-1)))
       t2(i,j,k) = t2(i,j,k) - m(i,j) * Half * rdx**4 * (  &
            dpsi(i-2,j-1,k) * mu(i-1,j) * ((c(i-1,j) * (  &
            (psi(i-1,j+1,k)-psi(i-1,j,k))*mu(i-1,j+1) + (psi(i-1,j,k)-psi(i-1,j-1,k))*mu(i-1,j))  &
            - c(i-1,j-1) * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j) +  &
            (psi(i-1,j-1,k)-psi(i-1,j-2,k))*mu(i-1,j-1))) * mv(i-1,j-1)  &
            + mu(i-2,j)*mu(i-1,j)/mCap(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))))
 
       t2(i,j,k) = t2(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i-1,j-1,k) * (mu(i,j) * ((  &
            c(i,j) * ((psi(i,j+1,k)-psi(i,j,k))*mu(i,j+1) + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) -  &
            c(i,j-1) * ((psi(i,j,k)-psi(i,j-1,k))*mu(i,j) + (psi(i,j-1,k)-psi(i,j-2,k))*mu(i,j-1) &
            )) * mv(i,j-1) + mu(i-1,j)*mu(i,j)/mCap(i,j) * (psi(i,j,k)-psi(i,j-1,k)))  &
            - mu(i-1,j) * ((c(i-1,j) * mu(i-1,j) + c(i-1,j-1) * (mu(i-1,j-1) - mu(i-1,j))) * (  &
            (psi(i,j,k)-psi(i-1,j,k)) * mv(i,j) + (psi(i-1,j,k)-psi(i-2,j,k)) * mv(i-1,j) +  &
            (psi(i,j-1,k)-psi(i-1,j-1,k))*mv(i,j-1) + (psi(i-1,j-1,k)-psi(i-2,j-1,k))*mv(i-1,j-1))&
            - (c(i-1,j) * ((psi(i-1,j+1,k)-psi(i-1,j,k)) * mu(i-1,j+1)  &
                                                    + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j))  &
            - c(i-1,j-1) * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)  &
                                               + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1)))  &
            * (mv(i-1,j-1) - mv(i,j-1)) - mu(i-1,j) * (mu(i,j)/mCap(i,j)  &
            * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (One/mCap(i,j) - One/mCap(i-1,j))  &
            * (psi(i-1,j,k)-psi(i-1,j-1,k))-mu(i-2,j)/mCap(i-1,j)*(psi(i-2,j,k)-psi(i-2,j-1,k)))  &
            - mu(i-1,j)**2 * (One/mCap(i,j) - One/mCap(i-1,j)) * (psi(i-1,j,k)-psi(i-1,j-1,k)))))

       t2(i,j,k) = t2(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i,j-1,k) * (mu(i,j) * ((c(i,j) * mu(i,j) + c(i,j-1) * (mu(i,j-1) - mu(i,j))) *  &
            ((psi(i+1,j,k)-psi(i,j,k)) * mv(i+1,j) + (psi(i,j,k)-psi(i-1,j,k)) * mv(i,j) +  &
            (psi(i+1,j-1,k)-psi(i,j-1,k))*mv(i+1,j-1) + (psi(i,j-1,k)-psi(i-1,j-1,k))*mv(i,j-1))  &
            - (c(i,j) * ((psi(i,j+1,k)-psi(i,j,k))*mu(i,j+1) + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) &
            - c(i,j-1)*((psi(i,j,k)-psi(i,j-1,k))*mu(i,j)+(psi(i,j-1,k)-psi(i,j-2,k))*mu(i,j-1))) &
            * (mv(i,j-1) - mv(i+1,j-1))  &
            - mu(i,j) * (mu(i+1,j)/mCap(i+1,j) * (psi(i+1,j,k)-psi(i+1,j-1,k))  &
            + mu(i,j) * (One/mCap(i+1,j) - One/mCap(i,j)) * (psi(i,j,k)-psi(i,j-1,k))  &
            - mu(i-1,j)/mCap(i,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) - mu(i,j)**2  &
            * (One/mCap(i+1,j) - One/mCap(i,j)) * (psi(i,j,k)-psi(i,j-1,k))) - mu(i-1,j) * ( -(  &
            c(i-1,j) * ((psi(i-1,j+1,k)-psi(i-1,j,k)) * mu(i-1,j+1)  &
            + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)) - c(i-1,j-1)  &
            * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)  &
            + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1))) * mv(i,j-1)  &
            - mu(i-1,j)*mu(i,j)/mCap(i,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)))))

       t2(i,j,k) = t2(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i+1,j-1,k) * mu(i,j) * ( -(c(i,j) * ((psi(i,j+1,k)-psi(i,j,k)) * mu(i,j+1)  &
            + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) - c(i,j-1) * ((psi(i,j,k)-psi(i,j-1,k))*mu(i,j)  &
            + (psi(i,j-1,k)-psi(i,j-2,k)) * mu(i,j-1))) * mv(i+1,j-1)  &
            - mu(i,j)*mu(i+1,j)/mCap(i+1,j) * (psi(i,j,k)-psi(i,j-1,k))))
       t2(i,j,k) = t2(i,j,k) - m(i,j) * Half * rdx**4 * (  &
            dpsi(i-2,j,k) * mu(i-1,j) * ((c(i-1,j) * ((psi(i-1,j+1,k)-psi(i-1,j,k)) * mu(i-1,j+1) &
            + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)) - c(i-1,j-1)  &
            * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)  &
            + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1))) * mv(i-1,j)  &
            - mu(i-2,j)*mu(i-1,j)/mCap(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))))

       t2(i,j,k) = t2(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i-1,j,k) * (mu(i,j) * ((c(i,j) * ((psi(i,j+1,k)-psi(i,j,k)) * mu(i,j+1)  &
            + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) - c(i,j-1) * ((psi(i,j,k)-psi(i,j-1,k))*mu(i,j)  &
            + (psi(i,j-1,k)-psi(i,j-2,k)) * mu(i,j-1))) * mv(i,j) - mu(i-1,j)*mu(i,j)/mCap(i,j)  &
            * (psi(i,j,k)-psi(i,j-1,k))) - mu(i-1,j) * (-(c(i-1,j) * (mu(i-1,j) - mu(i-1,j+1))  &
            - c(i-1,j-1) * mu(i-1,j)) * ((psi(i,j,k)-psi(i-1,j,k)) * mv(i,j)  &
            + (psi(i-1,j,k)-psi(i-2,j,k)) * mv(i-1,j) + (psi(i,j-1,k)-psi(i-1,j-1,k)) * mv(i,j-1) &
            + (psi(i-1,j-1,k)-psi(i-2,j-1,k)) * mv(i-1,j-1))  &
            - (c(i-1,j) * ((psi(i-1,j+1,k)-psi(i-1,j,k)) * mu(i-1,j+1)  &
            + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j))  &
            - c(i-1,j-1) * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)  &
            + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1))) * (mv(i-1,j) - mv(i,j))  &
            + mu(i-1,j) * (mu(i,j)/mCap(i,j) * (psi(i,j,k)-psi(i,j-1,k))  &
            + mu(i-1,j) * (One/mCap(i,j) - One/mCap(i-1,j)) * (psi(i-1,j,k)-psi(i-1,j-1,k))  &
            - mu(i-2,j)/mCap(i-1,j) * (psi(i-2,j,k)-psi(i-2,j-1,k)))  &
            + mu(i-1,j)**2 * (One/mCap(i,j) - One/mCap(i-1,j)) * (psi(i-1,j,k)-psi(i-1,j-1,k)))))

       t2(i,j,k) = t2(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i,j,k) * (mu(i,j) * (-(c(i,j) * (mu(i,j)-mu(i,j+1)) - c(i,j-1)*mu(i,j))  &
            * ((psi(i+1,j,k)-psi(i,j,k)) * mv(i+1,j) + (psi(i,j,k)-psi(i-1,j,k)) * mv(i,j)  &
            + (psi(i+1,j-1,k)-psi(i,j-1,k))*mv(i+1,j-1) + (psi(i,j-1,k)-psi(i-1,j-1,k))*mv(i,j-1))&
            - (c(i,j) * ((psi(i,j+1,k)-psi(i,j,k))*mu(i,j+1) + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) &
            - c(i,j-1) * ((psi(i,j,k)-psi(i,j-1,k)) * mu(i,j)  &
            + (psi(i,j-1,k)-psi(i,j-2,k)) * mu(i,j-1))) * (mv(i,j) - mv(i+1,j)) + mu(i,j)  &
            * (mu(i+1,j)/mCap(i+1,j) * (psi(i+1,j,k)-psi(i+1,j-1,k)) + mu(i,j)  &
            * (One/mCap(i+1,j) - One/mCap(i,j)) * (psi(i,j,k)-psi(i,j-1,k)) - mu(i-1,j)/mCap(i,j) &
            * (psi(i-1,j,k)-psi(i-1,j-1,k))) + mu(i,j)**2 * (One/mCap(i+1,j) - One/mCap(i,j))  &
            * (psi(i,j,k)-psi(i,j-1,k))) - mu(i-1,j) * (mu(i-1,j)*mu(i,j)/mCap(i,j)  &
            * (psi(i-1,j,k)-psi(i-1,j-1,k)) - (c(i,j) * ((psi(i-1,j+1,k)-psi(i-1,j,k))*mu(i-1,j+1)&
            + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)) - c(i-1,j-1)  &
            * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)  &
            + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1))) * mv(i,j))))

       t2(i,j,k) = t2(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i+1,j,k) * mu(i,j) * (mu(i,j)*mu(i+1,j)/mCap(i+1,j) * (psi(i,j,k)-psi(i,j-1,k))  &
            - (c(i,j) * ((psi(i,j+1,k)-psi(i,j,k))*mu(i,j+1) + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) &
            - c(i,j-1) * ((psi(i,j,k)-psi(i,j-1,k)) * mu(i,j)  &
            + (psi(i,j-1,k)-psi(i,j-2,k)) * mu(i,j-1))) * mv(i+1,j)))
       t2(i,j,k) = t2(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i-1,j+1,k) * c(i-1,j) * mu(i-1,j) * mu(i-1,j+1)  &
            * ((psi(i,j,k)-psi(i-1,j,k)) * mv(i,j) + (psi(i-1,j,k)-psi(i-2,j,k)) * mv(i-1,j)  &
            + (psi(i,j-1,k)-psi(i-1,j-1,k))*mv(i,j-1)+(psi(i-1,j-1,k)-psi(i-2,j-1,k))*mv(i-1,j-1)))
       t2(i,j,k) = t2(i,j,k) - m(i,j) * Half * rdx**4 * (  &
            dpsi(i,j+1,k) * c(i,j) * mu(i,j) * mu(i,j+1) * ((psi(i+1,j,k)-psi(i,j,k)) * mv(i+1,j) &
            + (psi(i,j,k)-psi(i-1,j,k)) * mv(i,j) + (psi(i+1,j-1,k)-psi(i,j-1,k)) * mv(i+1,j-1)  &
            + (psi(i,j-1,k)-psi(i-1,j-1,k)) * mv(i,j-1)))

       t3(i,j,k) = m(i,j) * Half * rdx**4 * (  &
            -dpsi(i-1,j-2,k) * mv(i,j-1) * (mu(i-1,j-1) * (c(i,j-1) * (mv(i+1,j-1)  &
            * (psi(i+1,j-1,k)-psi(i,j-1,k)) + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)))  &
            - c(i-1,j-1) * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))  &
            + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))))  &
            + mv(i,j-2)/mCap(i,j-1) * mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))))
       t3(i,j,k) = t3(i,j,k) - m(i,j) * Half * rdx**4 * (  &
            dpsi(i,j-2,k) * mv(i,j-1) * (mu(i,j-1) * (c(i,j-1) * (mv(i+1,j-1)  &
            * (psi(i+1,j-1,k)-psi(i,j-1,k)) + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)))  &
            - c(i-1,j-1) * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))  &
            + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))))  &
            - mv(i,j-2)/mCap(i,j-1) * mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))))
       t3(i,j,k) = t3(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i-2,j-1,k) * c(i-1,j-1) * mv(i-1,j-1) * mv(i,j-1) * (mu(i,j)  &
            * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)) + mu(i,j-1)  &
            * (psi(i,j-1,k)-psi(i,j-2,k)) + mu(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-1,j-2,k))))

       t3(i,j,k) = t3(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i-1,j-1,k) * (mv(i,j) * (mu(i-1,j)  &
            * (c(i,j) * (mv(i+1,j)*(psi(i+1,j,k)-psi(i,j,k)) + mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))) &
            - c(i-1,j)*(mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))+mv(i-1,j)*(psi(i-1,j,k)-psi(i-2,j,k)))) &
            + mv(i,j)*mv(i,j-1)/mCap(i,j) * (psi(i,j,k)-psi(i-1,j,k))) - mv(i,j-1) * (-mv(i,j-1)  &
            * (One/mCap(i,j) * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k))) - One/mCap(i,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i,j-2) * (psi(i,j-2,k)-psi(i-1,j-2,k))))  &
            - (mu(i-1,j-1) - mu(i-1,j)) * (c(i,j-1) * (mv(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k))&
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))) - c(i-1,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))))  &
            + mv(i,j-1)**2 * (psi(i,j-1,k)-psi(i-1,j-1,k)) * (One/mCap(i,j-1) - One/mCap(i,j))  &
            + (mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))  &
            + mu(i,j-1)*(psi(i,j-1,k)-psi(i,j-2,k)) + mu(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-1,j-2,k)))&
            * (c(i,j-1) * mv(i,j-1) + c(i-1,j-1) * (mv(i-1,j-1) - mv(i,j-1))))))

       t3(i,j,k) = t3(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i,j-1,k) * (mv(i,j) * (mu(i,j)  &
            * (c(i,j) * (mv(i+1,j)*(psi(i+1,j,k)-psi(i,j,k)) + mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))) &
            - c(i-1,j)*(mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))+mv(i-1,j)*(psi(i-1,j,k)-psi(i-2,j,k)))) &
            - mv(i,j)*mv(i,j-1)/mCap(i,j) * (psi(i,j,k)-psi(i-1,j,k))) - mv(i,j-1) * (mv(i,j-1)  &
            * (One/mCap(i,j) * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k))) - One/mCap(i,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i,j-2) * (psi(i,j-2,k)-psi(i-1,j-2,k))))  &
            - (mu(i,j-1) - mu(i,j)) * (c(i,j-1) * (mv(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k))  &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))) - c(i-1,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))))  &
            - (mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))  &
            + mu(i,j-1)*(psi(i,j-1,k)-psi(i,j-2,k)) + mu(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-1,j-2,k)))&
            * (c(i,j-1) * (mv(i,j-1) - mv(i+1,j-1)) - c(i-1,j-1) * mv(i,j-1))  &
            + m(i,j-1)**2 * (psi(i,j-1,k)-psi(i-1,j-1,k)) * (One/mCap(i,j) - One/mCap(i,j-1)))))

       t3(i,j,k) = t3(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i+1,j-1,k) * c(i,j-1) * mv(i,j-1) * mv(i+1,j-1) * (mu(i,j)  &
            * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)) + mu(i,j-1)  &
            * (psi(i,j-1,k)-psi(i,j-2,k)) + mu(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-1,j-2,k))))
       t3(i,j,k) = t3(i,j,k) - m(i,j) * Half * rdx**4 * (  &
            dpsi(i-2,j,k) * c(i-1,j) * mv(i-1,j) * mv(i,j)  &
            * (mu(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) + mu(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k))&
            + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))))

       t3(i,j,k) = t3(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i-1,j,k) * (mv(i,j) * (-(mu(i-1,j) - mu(i-1,j+1)) * (c(i,j) * (mv(i+1,j)  &
            * (psi(i+1,j,k)-psi(i,j,k)) + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k))) - c(i-1,j)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))))  &
            - mv(i,j) * (One/mCap(i,j+1) * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            + mv(i,j+1) * (psi(i,j+1,k)-psi(i-1,j+1,k))) - One/mCap(i,j)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))))  &
            + mv(i,j)**2 * (psi(i,j,k)-psi(i-1,j,k)) * (One/mCap(i,j) - One/mCap(i,j+1))  &
            - (mu(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) + mu(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k))&
            + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)))  &
            * (-c(i,j) * mv(i,j) - c(i-1,j) * (mv(i-1,j) - mv(i,j))))  &
            - mv(i,j-1) * (-mv(i,j-1)*mv(i,j)/mCap(i,j) * (psi(i,j-1,k)-psi(i-1,j-1,k))  &
            - mu(i-1,j) * (c(i,j-1) * (mv(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k))  &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))) - c(i-1,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k)))))))

       t3(i,j,k) = t3(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i,j,k) * (mv(i,j) * (-(mu(i,j) - mu(i,j+1)) * (c(i,j) * (mv(i+1,j)  &
            * (psi(i+1,j,k)-psi(i,j,k)) + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k))) - c(i-1,j)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))))  &
            - (mu(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) + mu(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k))&
            + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)))  &
            * (c(i,j) * (mv(i,j) - mv(i+1,j)) - c(i-1,j) * mv(i,j)) + mv(i,j) * (One/mCap(i,j+1)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i,j+1) * (psi(i,j+1,k)-psi(i-1,j+1,k)))  &
            - One/mCap(i,j) * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)))) - mv(i,j)**2 * (psi(i,j,k)-psi(i-1,j,k))&
            * (One/mCap(i,j) - One/mCap(i,j+1))) - mv(i,j-1) * (mv(i,j-1)*mv(i,j)/mCap(i,j)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) - mu(i,j) * (c(i,j-1) * (mv(i+1,j-1)  &
            * (psi(i+1,j-1,k)-psi(i,j-1,k)) + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)))  &
            - c(i-1,j-1) * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))  &
            + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k)))))))

       t3(i,j,k) = t3(i,j,k) - m(i,j) * Half * rdx**4 * (  &
            dpsi(i+1,j,k) * c(i,j) * mv(i,j) * mv(i+1,j)  &
            * (mu(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) + mu(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k))&
            + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))))
       t3(i,j,k) = t3(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i-1,j+1,k) * mv(i,j) * (-mu(i-1,j+1) * (c(i,j)  &
            * (mv(i+1,j) * (psi(i+1,j,k)-psi(i,j,k)) + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)))  &
            - c(i-1,j) * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            + mv(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))))  &
            - mv(i,j)*mv(i,j+1)/mCap(i,j+1) * (psi(i,j,k)-psi(i-1,j,k))))
       t3(i,j,k) = t3(i,j,k) + m(i,j) * Half * rdx**4 * (  &
            dpsi(i,j+1,k) * mv(i,j) * (mv(i,j)*mv(i,j+1)/mCap(i,j+1) * (psi(i,j,k)-psi(i-1,j,k))  &
            - mu(i,j+1) * (c(i,j) * (mv(i+1,j) * (psi(i+1,j,k)-psi(i,j,k))  &
            + mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))) - c(i-1,j) * (mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))  &
            + mv(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))))))
    end forall
    dlhs(2:nx-1,2:ny-1,:) = t1 - t2 - t3
    deallocate(c, dpsi, t3, t2, t1)
  end subroutine balance_lhs_tlm

  ! ===========================================================================

  subroutine balance_lhs_adj(t1, psi, f, mCap, m, mu, mv, rdx, ps)
    real(cp), dimension(:,:,:),   intent(in) :: t1         ! (1:nx,1:ny)
    real(cp), dimension(0:,0:,:), intent(in) :: psi        ! (0:nx,0:ny)
    real(cp), dimension(:,:),     intent(in) :: f, mCap, m ! (1:nx,1:ny)
    real(cp), dimension(0:,:),    intent(in) :: mu
    real(cp), dimension(:,0:),    intent(in) :: mv         ! (1:nx,0:ny)
    real(cp),                     intent(in) :: rdx
    real(cp), dimension(0:size(t1,1),0:size(t1,2),size(t1,3)), intent(out) :: ps

    character(len=*), parameter :: Mismatch = "mismatch in balance_lhs_adj!"
    
    real(cp), dimension(:,:), allocatable :: c, cc, dd
    integer :: nx, ny, nz, i, j, k
    
    nx = size(t1,1)
    ny = size(t1,2)
    nz = size(t1,3)

    ! Sanity checks
    if (any(shape(psi) /= (/ nx+1, ny+1, nz /) )) stop "psi "  // Mismatch   
    if (any(shape(f)  /= (/ nx,   ny   /) )) stop "f "  // Mismatch
    if (any(shape(mCap) /= (/ nx, ny   /) )) stop "mCap " // Mismatch    
    if (any(shape(m)  /= (/ nx,   ny   /) )) stop "m "  // Mismatch
    if (any(shape(mu) /= (/ nx+1, ny   /) )) stop "mu " // Mismatch
    if (any(shape(mv) /= (/ nx,   ny+1 /) )) stop "mv " // Mismatch

    allocate(c(nx-1,ny-1), cc(nx,ny), dd(nx,ny))

    forall (i = 1:nx-1, j = 1:ny-1)
       c(i,j) = One / (mCap(i,j) + mCap(i+1,j) + mCap(i,j+1) + mCap(i+1,j+1))
    end forall
    ps = Zero
    cc = m**2 * Eighth * rdx**2
    dd = m * Half * rdx**4

    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       ps(i-1,j-2,k) = ps(i-1,j-2,k) + cc(i,j) *                           &
            mu(i-1,j-1) / mv(i,j-1) * (f(i,j-1) + f(i,  j)) * t1(i,j,k)
       ps(i-1,j-2,k) = ps(i-1,j-2,k) - dd(i,j) *                           &
            c(i-1,j-1) * mu(i-1,j-1) * mu(i-1,j) * ((psi(i,j,k)-psi(i-1,j,k)) * mv(i,j)  &
            + (psi(i-1,j,k)-psi(i-2,j,k)) * mv(i-1,j) + (psi(i,j-1,k)-psi(i-1,j-1,k)) * mv(i,j-1) &
            + (psi(i-1,j-1,k)-psi(i-2,j-1,k)) * mv(i-1,j-1)) * t1(i,j,k)
       ps(i-1,j-2,k) = ps(i-1,j-2,k) + dd(i,j) *                           &
            mv(i,j-1) * (mu(i-1,j-1) * (c(i,j-1) * (mv(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k))  &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))) - c(i-1,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))))  &
            + mv(i,j-2)/mCap(i,j-1) * mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))) * t1(i,j,k)

       ps(i  ,j-2,k) = ps(i  ,j-2,k) + cc(i,j) *                           &
            mu(i,  j-1) / mv(i,j-1) * (f(i,j-1) + f(i,  j)) * t1(i,j,k)
       ps(i  ,j-2,k) = ps(i  ,j-2,k) + dd(i,j) *                           &
            c(i,j-1) * mu(i,j-1) * mu(i,j) * ((psi(i+1,j,k)-psi(i,j,k)) * mv(i+1,j)  &
            + (psi(i,j,k)-psi(i-1,j,k)) * mv(i,j) + (psi(i+1,j-1,k)-psi(i,j-1,k)) * mv(i+1,j-1)  &
            + (psi(i,j-1,k)-psi(i-1,j-1,k)) * mv(i,j-1)) * t1(i,j,k)
       ps(i  ,j-2,k) = ps(i  ,j-2,k) + dd(i,j) *                           &
            mv(i,j-1) * (mu(i,j-1) * (c(i,j-1) * (mv(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k))  &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))) - c(i-1,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))))  &
            - mv(i,j-2)/mCap(i,j-1) * mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))) * t1(i,j,k)

       ps(i-2,j-1,k) = ps(i-2,j-1,k) + cc(i,j) *                           &
            mv(i-1,j-1) / mu(i-1,j) * (f(i-1,j) + f(i,  j)) * t1(i,j,k)
       ps(i-2,j-1,k) = ps(i-2,j-1,k) + dd(i,j) *                           &
            mu(i-1,j) * ((c(i-1,j) * ((psi(i-1,j+1,k)-psi(i-1,j,k)) * mu(i-1,j+1)  &
            + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)) - c(i-1,j-1)      &
            * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)                   &
            + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1))) * mv(i-1,j-1)  &
            + mu(i-2,j)*mu(i-1,j)/mCap(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) * t1(i,j,k)
       ps(i-2,j-1,k) = ps(i-2,j-1,k) - dd(i,j) *                           &
            c(i-1,j-1) * mv(i-1,j-1) * mv(i,j-1) * (mu(i,j) * (psi(i,j,k)-psi(i,j-1,k))  &
            + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)) + mu(i,j-1) * (psi(i,j-1,k)-psi(i,j-2,k)) &
            + mu(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-1,j-2,k))) * t1(i,j,k)

       ps(i-1,j-1,k) = ps(i-1,j-1,k) - cc(i,j) * (                         &
            mv(i,  j-1) / mu(i,  j) * (f(i,j) + f(i+1,j)) +                &
            (f(i-1,j) + f(i,j)) / mu(i-1,j) * (mv(i-1,j-1) - mv(i,j-1)) +  &
            mu(i-1,  j) / mv(i,  j) * (f(i,j) + f(i,j+1)) +                &
            (f(i,j-1) + f(i,j)) / mv(i,j-1) * (mu(i-1,j-1) - mu(i-1,j)) )  &
            * t1(i,j,k)
       ps(i-1,j-1,k) = ps(i-1,j-1,k) - dd(i,j) *                           &
            (mu(i,j) * ((c(i,j) * ((psi(i,j+1,k)-psi(i,j,k)) * mu(i,j+1)   &
            + (psi(i,j,k)-psi(i,j-1,k)) * mu(i,j)) - c(i,j-1) * ((psi(i,j,k)-psi(i,j-1,k))*mu(i,j)&
            + (psi(i,j-1,k)-psi(i,j-2,k)) * mu(i,j-1))) * mv(i,j-1) + mu(i-1,j)*mu(i,j)/mCap(i,j) &
            * (psi(i,j,k)-psi(i,j-1,k))) - mu(i-1,j) * ((c(i-1,j) * mu(i-1,j) + c(i-1,j-1)  &
            * (mu(i-1,j-1) - mu(i-1,j))) * ((psi(i,j,k)-psi(i-1,j,k)) * mv(i,j)  &
            + (psi(i-1,j,k)-psi(i-2,j,k)) * mv(i-1,j) + (psi(i,j-1,k)-psi(i-1,j-1,k)) * mv(i,j-1) &
            + (psi(i-1,j-1,k)-psi(i-2,j-1,k)) * mv(i-1,j-1)) - (c(i-1,j)   &
            * ((psi(i-1,j+1,k)-psi(i-1,j,k)) * mu(i-1,j+1)                 &
                             + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j))  &
            - c(i-1,j-1) * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)      &
                        + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1)))  &
            * (mv(i-1,j-1) - mv(i,j-1)) - mu(i-1,j) * (mu(i,j)/mCap(i,j)   &
            * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (One/mCap(i,j) - One/mCap(i-1,j))  &
            * (psi(i-1,j,k)-psi(i-1,j-1,k))-mu(i-2,j)/mCap(i-1,j)*(psi(i-2,j,k)-psi(i-2,j-1,k)))  &
            - mu(i-1,j)**2 * (One/mCap(i,j) - One/mCap(i-1,j)) * (psi(i-1,j,k)-psi(i-1,j-1,k))))  &
            * t1(i,j,k)
       ps(i-1,j-1,k) = ps(i-1,j-1,k) - dd(i,j) *                           &
            (mv(i,j) * (mu(i-1,j) * (c(i,j) * (mv(i+1,j)*(psi(i+1,j,k)-psi(i,j,k))  &
            + mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))) - c(i-1,j)*(mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))  &
            + mv(i-1,j)*(psi(i-1,j,k)-psi(i-2,j,k)))) + mv(i,j)*mv(i,j-1)/mCap(i,j)  &
            * (psi(i,j,k)-psi(i-1,j,k))) - mv(i,j-1) * (-mv(i,j-1) * (One/mCap(i,j) * (mv(i,j)  &
            * (psi(i,j,k)-psi(i-1,j,k)) + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)))  &
            - One/mCap(i,j-1) * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))  &
            + mv(i,j-2) * (psi(i,j-2,k)-psi(i-1,j-2,k)))) - (mu(i-1,j-1) - mu(i-1,j)) * (c(i,j-1) &
            * (mv(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k))                 &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))) - c(i-1,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))))  &
            + mv(i,j-1)**2 * (psi(i,j-1,k)-psi(i-1,j-1,k)) * (One/mCap(i,j-1) - One/mCap(i,j))  &
            + (mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))  &
            + mu(i,j-1)*(psi(i,j-1,k)-psi(i,j-2,k)) + mu(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-1,j-2,k)))&
            * (c(i,j-1) * mv(i,j-1) + c(i-1,j-1) * (mv(i-1,j-1) - mv(i,j-1))))) * t1(i,j,k)

       ps(i  ,j-1,k) = ps(i  ,j-1,k) + cc(i,j) * (                         &
           -mv(i,  j-1) / mu(i-1,j) * (f(i-1,j) + f(i,j)) +                &
            (f(i,j) + f(i+1,j)) / mu(i,  j) * (mv(i,j-1) - mv(i+1,j-1)) -  &
            mu(i,    j) / mv(i,  j) * (f(i,j) + f(i,j+1)) +                &
            (f(i,j-1) + f(i,j)) / mv(i,j-1) * (mu(i,  j) - mu(i,  j-1)) )  &
            * t1(i,j,k)
       ps(i  ,j-1,k) = ps(i  ,j-1,k) - dd(i,j) *                           &
            (mu(i,j) * ((c(i,j) * mu(i,j) + c(i,j-1) * (mu(i,j-1) - mu(i,j)))  &
            * ((psi(i+1,j,k)-psi(i,j,k)) * mv(i+1,j) + (psi(i,j,k)-psi(i-1,j,k)) * mv(i,j)  &
            + (psi(i+1,j-1,k)-psi(i,j-1,k)) * mv(i+1,j-1)                  &
            + (psi(i,j-1,k)-psi(i-1,j-1,k))*mv(i,j-1))                     &
            - (c(i,j) * ((psi(i,j+1,k)-psi(i,j,k))*mu(i,j+1) + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) &
            - c(i,j-1)*((psi(i,j,k)-psi(i,j-1,k))*mu(i,j)+(psi(i,j-1,k)-psi(i,j-2,k))*mu(i,j-1))) &
            * (mv(i,j-1) - mv(i+1,j-1))                                    &
            - mu(i,j) * (mu(i+1,j)/mCap(i+1,j) * (psi(i+1,j,k)-psi(i+1,j-1,k))  &
            + mu(i,j) * (One/mCap(i+1,j) - One/mCap(i,j)) * (psi(i,j,k)-psi(i,j-1,k))  &
            - mu(i-1,j)/mCap(i,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) - mu(i,j)**2  &
            * (One/mCap(i+1,j) - One/mCap(i,j)) * (psi(i,j,k)-psi(i,j-1,k))) - mu(i-1,j) * ( -(  &
            c(i-1,j) * ((psi(i-1,j+1,k)-psi(i-1,j,k)) * mu(i-1,j+1)        &
            + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)) - c(i-1,j-1)      &
            * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)                   &
            + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1))) * mv(i,j-1)  &
            - mu(i-1,j)*mu(i,j)/mCap(i,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)))) * t1(i,j,k)
       ps(i  ,j-1,k) = ps(i  ,j-1,k) - dd(i,j) *                           &
            (mv(i,j) * (mu(i,j) * (c(i,j) * (mv(i+1,j)*(psi(i+1,j,k)-psi(i,j,k))  &
            + mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))) - c(i-1,j)*(mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))  &
            + mv(i-1,j)*(psi(i-1,j,k)-psi(i-2,j,k)))) - mv(i,j)*mv(i,j-1)/mCap(i,j)  &
            * (psi(i,j,k)-psi(i-1,j,k))) - mv(i,j-1) * (mv(i,j-1) * (One/mCap(i,j)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)))  &
            - One/mCap(i,j-1) * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))  &
            + mv(i,j-2) * (psi(i,j-2,k)-psi(i-1,j-2,k)))) - (mu(i,j-1) - mu(i,j)) * (c(i,j-1)  &
            * (mv(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k))                 &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))) - c(i-1,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))))  &
            - (mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))  &
            + mu(i,j-1)*(psi(i,j-1,k)-psi(i,j-2,k)) + mu(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-1,j-2,k)))&
            * (c(i,j-1) * (mv(i,j-1) - mv(i+1,j-1)) - c(i-1,j-1) * mv(i,j-1))  &
            + m(i,j-1)**2 * (psi(i,j-1,k)-psi(i-1,j-1,k)) * (One/mCap(i,j) - One/mCap(i,j-1))))  &
            * t1(i,j,k)

       ps(i+1,j-1,k) = ps(i+1,j-1,k) + cc(i,j) *                           &
            mv(i+1,j-1) / mu(i,  j) * (f(i,j) + f(i+1,j)) * t1(i,j,k)
       ps(i+1,j-1,k) = ps(i+1,j-1,k) - dd(i,j) *                           &
            mu(i,j) * ( -(c(i,j) * ((psi(i,j+1,k)-psi(i,j,k)) * mu(i,j+1)  &
            + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) - c(i,j-1) * ((psi(i,j,k)-psi(i,j-1,k))*mu(i,j)  &
            + (psi(i,j-1,k)-psi(i,j-2,k)) * mu(i,j-1))) * mv(i+1,j-1)      &
            - mu(i,j)*mu(i+1,j)/mCap(i+1,j) * (psi(i,j,k)-psi(i,j-1,k))) * t1(i,j,k)
       ps(i+1,j-1,k) = ps(i+1,j-1,k) - dd(i,j) *                           &
            c(i,j-1) * mv(i,j-1) * mv(i+1,j-1) * (mu(i,j) * (psi(i,j,k)-psi(i,j-1,k))  &
            + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)) + mu(i,j-1) * (psi(i,j-1,k)-psi(i,j-2,k)) &
            + mu(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-1,j-2,k))) * t1(i,j,k)

       ps(i-2,j  ,k) = ps(i-2,j  ,k) + cc(i,j) *                           &
            mv(i-1,j  ) / mu(i-1,j) * (f(i-1,j) + f(i,j)) * t1(i,j,k)
       ps(i-2,j  ,k) = ps(i-2,j  ,k) + dd(i,j) *                           &
            mu(i-1,j) * ((c(i-1,j) * ((psi(i-1,j+1,k)-psi(i-1,j,k)) * mu(i-1,j+1) &
            + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)) - c(i-1,j-1)      &
            * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)                   &
            + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1))) * mv(i-1,j)  &
            - mu(i-2,j)*mu(i-1,j)/mCap(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) * t1(i,j,k)
       ps(i-2,j  ,k) = ps(i-2,j  ,k) + dd(i,j) *                           &
            c(i-1,j) * mv(i-1,j) * mv(i,j) * (mu(i,j+1) * (psi(i,j+1,k)-psi(i,j,k))  &
            + mu(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k)) + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k))  &
            + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) * t1(i,j,k)

       ps(i-1,j  ,k) = ps(i-1,j  ,k) - cc(i,j) * (                         &
            mv(i  ,j  ) / mu(i,  j) * (f(i,j) + f(i+1,j)) +                &
            (f(i-1,j) + f(i,j)) / mu(i-1,j) * (mv(i-1,j  ) - mv(i,  j)) +  &
            mu(i-1,j  ) / mv(i,j-1) * (f(i,j-1) + f(i,j)) +                &
            (f(i,j) + f(i,j+1)) / mv(i,  j) * (mu(i-1,j+1) - mu(i-1,j)) )  &
            * t1(i,j,k)
       ps(i-1,j  ,k) = ps(i-1,j  ,k) - dd(i,j) *                           &
            (mu(i,j) * ((c(i,j) * ((psi(i,j+1,k)-psi(i,j,k)) * mu(i,j+1)   &
            + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) - c(i,j-1) * ((psi(i,j,k)-psi(i,j-1,k))*mu(i,j)  &
            + (psi(i,j-1,k)-psi(i,j-2,k)) * mu(i,j-1))) * mv(i,j) - mu(i-1,j)*mu(i,j)/mCap(i,j)  &
            * (psi(i,j,k)-psi(i,j-1,k))) - mu(i-1,j) * (-(c(i-1,j) * (mu(i-1,j) - mu(i-1,j+1))  &
            - c(i-1,j-1) * mu(i-1,j)) * ((psi(i,j,k)-psi(i-1,j,k)) * mv(i,j)  &
            + (psi(i-1,j,k)-psi(i-2,j,k)) * mv(i-1,j) + (psi(i,j-1,k)-psi(i-1,j-1,k)) * mv(i,j-1) &
            + (psi(i-1,j-1,k)-psi(i-2,j-1,k)) * mv(i-1,j-1))               &
            - (c(i-1,j) * ((psi(i-1,j+1,k)-psi(i-1,j,k)) * mu(i-1,j+1)     &
            + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j))                   &
            - c(i-1,j-1) * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)      &
            + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1))) * (mv(i-1,j) - mv(i,j))  &
            + mu(i-1,j) * (mu(i,j)/mCap(i,j) * (psi(i,j,k)-psi(i,j-1,k))   &
            + mu(i-1,j) * (One/mCap(i,j) - One/mCap(i-1,j)) * (psi(i-1,j,k)-psi(i-1,j-1,k))  &
            - mu(i-2,j)/mCap(i-1,j) * (psi(i-2,j,k)-psi(i-2,j-1,k)))       &
            + mu(i-1,j)**2 * (One/mCap(i,j) - One/mCap(i-1,j)) * (psi(i-1,j,k)-psi(i-1,j-1,k))))  &
            * t1(i,j,k)
       ps(i-1,j  ,k) = ps(i-1,j  ,k) - dd(i,j) *                           &
            (mv(i,j) * (-(mu(i-1,j) - mu(i-1,j+1)) * (c(i,j) * (mv(i+1,j)  &
            * (psi(i+1,j,k)-psi(i,j,k)) + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k))) - c(i-1,j)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))))  &
            - mv(i,j) * (One/mCap(i,j+1) * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            + mv(i,j+1) * (psi(i,j+1,k)-psi(i-1,j+1,k))) - One/mCap(i,j)   &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))))  &
            + mv(i,j)**2 * (psi(i,j,k)-psi(i-1,j,k)) * (One/mCap(i,j) - One/mCap(i,j+1))  &
            - (mu(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) + mu(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k))&
            + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)))  &
            * (-c(i,j) * mv(i,j) - c(i-1,j) * (mv(i-1,j) - mv(i,j))))      &
            - mv(i,j-1) * (-mv(i,j-1)*mv(i,j)/mCap(i,j) * (psi(i,j-1,k)-psi(i-1,j-1,k))  &
            - mu(i-1,j) * (c(i,j-1) * (mv(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k))  &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))) - c(i-1,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))))))  &
            * t1(i,j,k)

       ps(i  ,j  ,k) = ps(i  ,j  ,k) + cc(i,j) * (                         &
           -mv(i  ,j  ) / mu(i-1,j) * (f(i-1,j) + f(i,j)) +                &
            (f(i,j) + f(i+1,j)) / mu(i,  j) * (mv(i,    j) - mv(i+1,j)) -  &
            mu(i  ,j  ) / mv(i,j-1) * (f(i,j-1) + f(i,j)) +                &
            (f(i,j) + f(i,j+1)) / mv(i,  j) * (mu(i,    j) - mu(i,j+1)) )  &
            * t1(i,j,k)
       ps(i  ,j  ,k) = ps(i  ,j  ,k) - dd(i,j) *                           &
            (mu(i,j) * (-(c(i,j) * (mu(i,j)-mu(i,j+1)) - c(i,j-1)*mu(i,j))  &
            * ((psi(i+1,j,k)-psi(i,j,k)) * mv(i+1,j) + (psi(i,j,k)-psi(i-1,j,k)) * mv(i,j)  &
            + (psi(i+1,j-1,k)-psi(i,j-1,k))*mv(i+1,j-1) + (psi(i,j-1,k)-psi(i-1,j-1,k))*mv(i,j-1))&
            - (c(i,j) * ((psi(i,j+1,k)-psi(i,j,k))*mu(i,j+1) + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) &
            - c(i,j-1) * ((psi(i,j,k)-psi(i,j-1,k)) * mu(i,j)              &
            + (psi(i,j-1,k)-psi(i,j-2,k)) * mu(i,j-1))) * (mv(i,j) - mv(i+1,j)) + mu(i,j)  &
            * (mu(i+1,j)/mCap(i+1,j) * (psi(i+1,j,k)-psi(i+1,j-1,k)) + mu(i,j)  &
            * (One/mCap(i+1,j) - One/mCap(i,j)) * (psi(i,j,k)-psi(i,j-1,k)) - mu(i-1,j)/mCap(i,j) &
            * (psi(i-1,j,k)-psi(i-1,j-1,k))) + mu(i,j)**2 * (One/mCap(i+1,j) - One/mCap(i,j))  &
            * (psi(i,j,k)-psi(i,j-1,k))) - mu(i-1,j) * (mu(i-1,j)*mu(i,j)/mCap(i,j)  &
            * (psi(i-1,j,k)-psi(i-1,j-1,k)) - (c(i,j) * ((psi(i-1,j+1,k)-psi(i-1,j,k))*mu(i-1,j+1)&
            + (psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)) - c(i-1,j-1)      &
            * ((psi(i-1,j,k)-psi(i-1,j-1,k)) * mu(i-1,j)                   &
            + (psi(i-1,j-1,k)-psi(i-1,j-2,k)) * mu(i-1,j-1))) * mv(i,j))) * t1(i,j,k)
       ps(i  ,j  ,k) = ps(i  ,j  ,k) - dd(i,j) *                           &
            (mv(i,j) * (-(mu(i,j) - mu(i,j+1)) * (c(i,j) * (mv(i+1,j) * (psi(i+1,j,k)-psi(i,j,k)) &
            + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k))) - c(i-1,j)              &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))))  &
            - (mu(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) + mu(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k))&
            + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)))  &
            * (c(i,j) * (mv(i,j) - mv(i+1,j)) - c(i-1,j) * mv(i,j)) + mv(i,j) * (One/mCap(i,j+1)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i,j+1) * (psi(i,j+1,k)-psi(i-1,j+1,k)))  &
            - One/mCap(i,j) * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k))         &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)))) - mv(i,j)**2 * (psi(i,j,k)-psi(i-1,j,k))&
            * (One/mCap(i,j) - One/mCap(i,j+1))) - mv(i,j-1) * (mv(i,j-1)*mv(i,j)/mCap(i,j)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) - mu(i,j) * (c(i,j-1) * (mv(i+1,j-1)  &
            * (psi(i+1,j-1,k)-psi(i,j-1,k)) + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)))  &
            - c(i-1,j-1) * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))      &
            + mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k)))))) * t1(i,j,k)

       ps(i+1,j  ,k) = ps(i+1,j  ,k) + cc(i,j) *                           &
            mv(i+1,j  ) / mu(i,  j) * (f(i,j) + f(i+1,j)) * t1(i,j,k)
       ps(i+1,j  ,k) = ps(i+1,j  ,k) - dd(i,j) *                           &
            mu(i,j) * (mu(i,j)*mu(i+1,j)/mCap(i+1,j) * (psi(i,j,k)-psi(i,j-1,k))  &
            - (c(i,j) * ((psi(i,j+1,k)-psi(i,j,k))*mu(i,j+1) + (psi(i,j,k)-psi(i,j-1,k))*mu(i,j)) &
            - c(i,j-1) * ((psi(i,j,k)-psi(i,j-1,k)) * mu(i,j)              &
            + (psi(i,j-1,k)-psi(i,j-2,k)) * mu(i,j-1))) * mv(i+1,j)) * t1(i,j,k)
       ps(i+1,j  ,k) = ps(i+1,j  ,k) + dd(i,j) *                           &
            c(i,j) * mv(i,j) * mv(i+1,j) * (mu(i,j+1) * (psi(i,j+1,k)-psi(i,j,k))  &
            + mu(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k)) + mu(i,j)        &
            * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) * t1(i,j,k)

       ps(i-1,j+1,k) = ps(i-1,j+1,k) + cc(i,j) *                           &
            mu(i-1,j+1) / mv(i,  j) * (f(i,j) + f(i,j+1)) * t1(i,j,k)
       ps(i-1,j+1,k) = ps(i-1,j+1,k) - dd(i,j) *                           &       
            c(i-1,j) * mu(i-1,j) * mu(i-1,j+1) * ((psi(i,j,k)-psi(i-1,j,k)) * mv(i,j)  &
            + (psi(i-1,j,k)-psi(i-2,j,k)) * mv(i-1,j) + (psi(i,j-1,k)-psi(i-1,j-1,k)) * mv(i,j-1) &
            + (psi(i-1,j-1,k)-psi(i-2,j-1,k)) * mv(i-1,j-1)) * t1(i,j,k)
       ps(i-1,j+1,k) = ps(i-1,j+1,k) - dd(i,j) *                           &
            mv(i,j) * (-mu(i-1,j+1) * (c(i,j) * (mv(i+1,j) * (psi(i+1,j,k)-psi(i,j,k))  &
            + mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))) - c(i-1,j) * (mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))  &
            + mv(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))))                    &
            - mv(i,j)*mv(i,j+1)/mCap(i,j+1) * (psi(i,j,k)-psi(i-1,j,k))) * t1(i,j,k)

       ps(i  ,j+1,k) = ps(i  ,j+1,k) + cc(i,j) *                           &
            mu(i,  j+1) / mv(i,  j) * (f(i,j) + f(i,j+1)) * t1(i,j,k)
       ps(i  ,j+1,k) = ps(i  ,j+1,k) + dd(i,j) *                           &      
            c(i,j) * mu(i,j) * mu(i,j+1) * ((psi(i+1,j,k)-psi(i,j,k)) * mv(i+1,j) &
            + (psi(i,j,k)-psi(i-1,j,k)) * mv(i,j) + (psi(i+1,j-1,k)-psi(i,j-1,k)) * mv(i+1,j-1)  &
            + (psi(i,j-1,k)-psi(i-1,j-1,k)) * mv(i,j-1)) * t1(i,j,k)
       ps(i  ,j+1,k) = ps(i  ,j+1,k) - dd(i,j) *                           &
            mv(i,j) * (mv(i,j)*mv(i,j+1)/mCap(i,j+1) * (psi(i,j,k)-psi(i-1,j,k))  &
            - mu(i,j+1) * (c(i,j) * (mv(i+1,j) * (psi(i+1,j,k)-psi(i,j,k))  &
            + mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))) - c(i-1,j) * (mv(i,j)*(psi(i,j,k)-psi(i-1,j,k))  &
            + mv(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))))) * t1(i,j,k)
    end forall
    ps((/0,nx/),:,:) = Zero
    ps(1:nx-1,(/0,ny/),:) = Zero
    deallocate(dd, cc, c)
  end subroutine balance_lhs_adj

  ! ===========================================================================

  subroutine balance_rhs_simp(phi, mCap, m, etaw, etam, rdx, rhs)
    real(cp), dimension(:,:,0:), intent(in) :: phi
    real(cp), dimension(:,:),    intent(in) :: mCap, m
    real(cp), dimension(0:),     intent(in) :: etaw
    real(cp), dimension(:),      intent(in) :: etam
    real(cp),                    intent(in) :: rdx
    real, dimension(size(phi,1),size(phi,2),size(etam)), intent(out) :: rhs

    character(len=*), parameter :: Mismatch = "mismatch in balance_rhs_simp!"

    real(cp), dimension(:,:), allocatable :: dpde, e1, e2
    integer :: nx, ny, nz, k

    nx = size(phi,1)
    ny = size(phi,2)
    nz = size(phi,3) - 1

    ! Sanity checks
    if (any(shape(mCap) /= (/ nx, ny /) )) stop "mCap " // Mismatch
    if (any(shape(m) /= (/ nx, ny /) )) stop "m " // Mismatch
    if (size(etaw) /= nz+1) stop "etaw " // Mismatch
    if (size(etam) /= nz) stop "etam " // Mismatch

    allocate(dpde(nx,ny), e1(1:nx-1,2:ny-1), e2(2:nx-1,1:ny-1))

    rhs = Zero
    do k = 1, nz
       dpde = (phi(:,:,k) - phi(:,:,k-1)) / (etaw(k) - etaw(k-1))
       
       e1 = Half*(xdiff(phi(:,2:ny-1,k-1),rdx)+xdiff(phi(:,2:ny-1,k),rdx)) *  &
            xavg(mCap(:,2:ny-1)) - etam(k) * xdiff(mCap(:,2:ny-1), rdx) *  &
            xavg(dpde(:,2:ny-1))
       e2 = Half*(ydiff(phi(2:nx-1,:,k-1),rdx)+ydiff(phi(2:nx-1,:,k),rdx)) *  &
            yavg(mCap(2:nx-1,:)) - etam(k) * ydiff(mCap(2:nx-1,:), rdx) *  &
            yavg(dpde(2:nx-1,:))
       rhs(2:nx-1,2:ny-1,k) = m(2:nx-1,2:ny-1)**2 *  &
            (xdiff(e1, rdx) + ydiff(e2, rdx))
    end do
    deallocate(e2, e1, dpde)
  end subroutine balance_rhs_simp

  ! ===========================================================================

  subroutine balance_rhs_simp_exp(phi, mCap, m, etaw, etam, rdx, rhs)
    real(cp), dimension(:,:,0:), intent(in) :: phi
    real(cp), dimension(:,:),    intent(in) :: mCap, m
    real(cp), dimension(0:),     intent(in) :: etaw
    real(cp), dimension(:),      intent(in) :: etam
    real(cp),                    intent(in) :: rdx
    real, dimension(size(phi,1),size(phi,2),size(etam)), intent(out) :: rhs

    character(len=*), parameter :: Mismatch = "mismatch in balance_rhs_simp_exp!"
    
    integer :: nx, ny, nz, i, j, k
    
    nx = size(phi,1)
    ny = size(phi,2)
    nz = size(phi,3) - 1

    ! Sanity checks
    if (any(shape(mCap) /= (/ nx, ny /) )) stop "mCap " // Mismatch
    if (any(shape(m) /= (/ nx, ny /) )) stop "m " // Mismatch
    if (size(etaw) /= nz+1) stop "etaw " // Mismatch
    if (size(etam) /= nz) stop "etam " // Mismatch

    rhs = Zero
    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       rhs(i,j,k) = m(i,j)**2 * Half * rdx**2 * (Half  &
            * ((phi(i+1,j,k-1)-phi(i,j,k-1) + phi(i+1,j,k)-phi(i,j,k)) * (mCap(i,j) + mCap(i+1,j))&
            - (phi(i,j,k-1)-phi(i-1,j,k-1) + phi(i,j,k)-phi(i-1,j,k)) * (mCap(i-1,j) + mCap(i,j)) &
            + (phi(i,j+1,k-1)-phi(i,j,k-1) + phi(i,j+1,k)-phi(i,j,k)) * (mCap(i,j) + mCap(i,j+1)) &
            - (phi(i,j,k-1)-phi(i,j-1,k-1) + phi(i,j,k)-phi(i,j-1,k)) * (mCap(i,j-1) + mCap(i,j)))&
            + etam(k) / (etaw(k) - etaw(k-1)) * (-(mCap(i+1,j) - mCap(i,j))  &
            * (phi(i,j,k)-phi(i,j,k-1) + phi(i+1,j,k)-phi(i+1,j,k-1)) + (mCap(i,j) - mCap(i-1,j)) &
            * (phi(i-1,j,k)-phi(i-1,j,k-1) + phi(i,j,k)-phi(i,j,k-1)) - (mCap(i,j+1) - mCap(i,j)) &
            * (phi(i,j,k)-phi(i,j,k-1) + phi(i,j+1,k)-phi(i,j+1,k-1)) + (mCap(i,j) - mCap(i,j-1)) &
            * (phi(i,j-1,k)-phi(i,j-1,k-1) + phi(i,j,k)-phi(i,j,k-1))))
    end forall
  end subroutine balance_rhs_simp_exp

  ! ===========================================================================

  subroutine balance_rhs_simp_tlm(dph, mCap, m, etaw, etam, rdx, drhs)
    real(cp), dimension(:,:,0:), intent(in) :: dph
    real(cp), dimension(:,:),    intent(in) :: mCap, m
    real(cp), dimension(0:),     intent(in) :: etaw
    real(cp), dimension(:),      intent(in) :: etam
    real(cp),                    intent(in) :: rdx
    real, dimension(size(dph,1),size(dph,2),size(etam)), intent(out) :: drhs
    
    character(len=*), parameter :: Mismatch = "mismatch in balance_rhs_simp_tlm!"
    
    real(cp), dimension(:,:,:), allocatable :: dphi
    integer :: nx, ny, nz, i, j, k

    nx = size(dph,1)
    ny = size(dph,2)
    nz = size(dph,3) - 1
    
    ! Sanity checks
    if (any(shape(mCap) /= (/ nx, ny /) )) stop "mCap " // Mismatch
    if (any(shape(m) /= (/ nx, ny /) )) stop "m " // Mismatch
    if (size(etaw) /= nz+1) stop "etaw " // Mismatch
    if (size(etam) /= nz) stop "etam " // Mismatch

    allocate(dphi(nx,ny,0:nz))
    
    dphi = dph
    dphi((/1,nx/),:,:) = Zero
    dphi(2:nx-1,(/1,ny/),:) = Zero
    dphi(2:nx-1,2:ny-1,(/0,nz/)) = Zero
    drhs = Zero  ! Really just setting the boundaries

    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       drhs(i,j,k) = m(i,j)**2 * Half * rdx**2 * (  &
            dphi(i,j-1,k-1) * (etam(k) / (etaw(k) - etaw(k-1)) * (mCap(i,j-1) - mCap(i,j))  &
                                                        + Half * (mCap(i,j) + mCap(i,j-1)))  &
            + dphi(i-1,j,k-1) * (etam(k) / (etaw(k) - etaw(k-1)) * (mCap(i-1,j) - mCap(i,j))  &
                                                          + Half * (mCap(i,j) + mCap(i-1,j)))  &
            + dphi(i,j,k-1) * (etam(k) / (etaw(k) - etaw(k-1))  &
                    * (mCap(i,j-1) + mCap(i-1,j) - Four * mCap(i,j) + mCap(i+1,j) + mCap(i,j+1))  &
             - Half * (mCap(i,j-1) + mCap(i-1,j) + Four * mCap(i,j) + mCap(i+1,j) + mCap(i,j+1))) &
            + dphi(i+1,j,k-1) * (etam(k) / (etaw(k) - etaw(k-1)) * (mCap(i+1,j) - mCap(i,j))  &
                                                          + Half * (mCap(i,j) + mCap(i+1,j)))  &
            + dphi(i,j+1,k-1) * (etam(k) / (etaw(k) - etaw(k-1)) * (mCap(i,j+1) - mCap(i,j))  &
                                                          + Half * (mCap(i,j+1) + mCap(i,j)))  &
            + dphi(i,j-1,k) * (etam(k) / (etaw(k) - etaw(k-1)) * (mCap(i,j) - mCap(i,j-1))  &
                                                          + Half * (mCap(i,j) + mCap(i,j-1)))  &
            + dphi(i-1,j,k) * (etam(k) / (etaw(k) - etaw(k-1)) * (mCap(i,j) - mCap(i-1,j))  &
                                                          + Half * (mCap(i,j) + mCap(i-1,j)))  &
            + dphi(i,j,k) * (etam(k) / (etaw(k) - etaw(k-1))  &
                    * (Four * mCap(i,j) - mCap(i,j-1) - mCap(i-1,j) - mCap(i+1,j) - mCap(i,j+1))  &
             - Half * (mCap(i,j-1) + mCap(i-1,j) + Four * mCap(i,j) + mCap(i+1,j) + mCap(i,j+1))) &
            + dphi(i+1,j,k) * (etam(k) / (etaw(k) - etaw(k-1)) * (mCap(i,j) - mCap(i+1,j))  &
                                                        + Half * (mCap(i,j) + mCap(i+1,j)))  &
            + dphi(i,j+1,k) * (etam(k) / (etaw(k) - etaw(k-1)) * (mCap(i,j) - mCap(i,j+1))  &
                                                        + Half * (mCap(i,j) + mCap(i,j+1))))
    end forall
    deallocate(dphi)
  end subroutine balance_rhs_simp_tlm

  ! ===========================================================================

  subroutine balance_rhs_simp_adj(adjrhs, mCap, m, etaw, etam, rdx, ph)
    real(cp), dimension(:,:,:),  intent(in) :: adjrhs
    real(cp), dimension(:,:),    intent(in) :: mCap, m
    real(cp), dimension(0:),     intent(in) :: etaw
    real(cp), dimension(:),      intent(in) :: etam
    real(cp),                    intent(in) :: rdx
    real(cp), dimension(size(adjrhs,1),size(adjrhs,2),0:size(adjrhs,3)), intent(out) :: ph

    character(len=*), parameter :: Mismatch = "mismatch in balance_rhs_simp_adj!"
    
    real(cp), dimension(:,:), allocatable :: c
    integer :: nx, ny, nz, i, j, k

    nx = size(adjrhs,1)
    ny = size(adjrhs,2)
    nz = size(adjrhs,3)
    
    ! Sanity checks
    if (any(shape(mCap) /= (/ nx, ny /) )) stop "mCap " // Mismatch
    if (any(shape(m) /= (/ nx, ny /) )) stop "m " // Mismatch
    if (size(etaw) /= nz+1) stop "etaw " // Mismatch
    if (size(etam) /= nz) stop "etam " // Mismatch

    allocate(c(nx,ny))

    ph = Zero
    c = m**2 * Half * rdx**2
    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       ph(i,j-1,k-1) = ph(i,j-1,k-1) + c(i,j) * (etam(k) / (etaw(k) - etaw(k-1))  &
            * (mCap(i,j-1) - mCap(i,j)) + Half * (mCap(i,j) + mCap(i,j-1))) * adjrhs(i,j,k)
       ph(i-1,j,k-1) = ph(i-1,j,k-1) + c(i,j) * (etam(k) / (etaw(k) - etaw(k-1))  &
            * (mCap(i-1,j) - mCap(i,j)) + Half * (mCap(i,j) + mCap(i-1,j))) * adjrhs(i,j,k)
       ph(i,j,k-1) = ph(i,j,k-1) + c(i,j) * (etam(k) / (etaw(k) - etaw(k-1))  &
            * (mCap(i,j-1) + mCap(i-1,j) - Four * mCap(i,j) + mCap(i+1,j) + mCap(i,j+1))  &
            - Half * (mCap(i,j-1) + mCap(i-1,j) + Four * mCap(i,j) + mCap(i+1,j) + mCap(i,j+1)))  &
            * adjrhs(i,j,k)
       ph(i+1,j,k-1) = ph(i+1,j,k-1) + c(i,j) * (etam(k) / (etaw(k) - etaw(k-1))  &
            * (mCap(i+1,j) - mCap(i,j)) + Half * (mCap(i,j) + mCap(i+1,j))) * adjrhs(i,j,k)
       ph(i,j+1,k-1) = ph(i,j+1,k-1) + c(i,j) * (etam(k) / (etaw(k) - etaw(k-1))  &
            * (mCap(i,j+1) - mCap(i,j)) + Half * (mCap(i,j+1) + mCap(i,j))) * adjrhs(i,j,k)
       ph(i,j-1,k) = ph(i,j-1,k) + c(i,j) * (etam(k) / (etaw(k) - etaw(k-1))  &
            * (mCap(i,j) - mCap(i,j-1)) + Half * (mCap(i,j) + mCap(i,j-1))) * adjrhs(i,j,k)
       ph(i-1,j,k) = ph(i-1,j,k) + c(i,j) * (etam(k) / (etaw(k) - etaw(k-1))  &
            * (mCap(i,j) - mCap(i-1,j)) + Half * (mCap(i,j) + mCap(i-1,j))) * adjrhs(i,j,k)
       ph(i,j,k) = ph(i,j,k) + c(i,j) * (etam(k) / (etaw(k) - etaw(k-1))  &
            * (Four * mCap(i,j) - mCap(i,j-1) - mCap(i-1,j) - mCap(i+1,j) - mCap(i,j+1))  &
            - Half * (mCap(i,j-1) + mCap(i-1,j) + Four * mCap(i,j) + mCap(i+1,j) + mCap(i,j+1)))  &
            * adjrhs(i,j,k)
       ph(i+1,j,k) = ph(i+1,j,k) + c(i,j) * (etam(k) / (etaw(k) - etaw(k-1))  &
            * (mCap(i,j) - mCap(i+1,j)) + Half * (mCap(i,j) + mCap(i+1,j))) * adjrhs(i,j,k)
       ph(i,j+1,k) = ph(i,j+1,k) + c(i,j) * (etam(k) / (etaw(k) - etaw(k-1))  &
            * (mCap(i,j) - mCap(i,j+1)) + Half * (mCap(i,j) + mCap(i,j+1))) * adjrhs(i,j,k)
    end forall
    ph((/1,nx/),:,:) = Zero
    ph(2:nx-1,(/1,ny/),:) = Zero
    ph(2:nx-1,2:ny-1,(/0,nz/)) = Zero
    deallocate(c)
  end subroutine balance_rhs_simp_adj

  ! ===========================================================================

  subroutine imbalance_simp_adj(adjimb, psi, f, mCap, m, mu, mv, etaw, etam, rdx, adjpsi, adjphi)
    real(cp), dimension(:,:,:),   intent(in) :: adjimb
    real(cp), dimension(0:,0:,:), intent(in) :: psi
    real(cp), dimension(:,:),     intent(in) :: f, mCap, m
    real(cp), dimension(0:,:),    intent(in) :: mu
    real(cp), dimension(:,0:),    intent(in) :: mv
    real(cp), dimension(0:),      intent(in) :: etaw
    real(cp), dimension(:),       intent(in) :: etam
    real(cp),                     intent(in) :: rdx
    real(cp), dimension(0:size(adjimb,1),0:size(adjimb,2),size(adjimb,3)), intent(out) :: adjpsi
    real(cp), dimension(size(adjimb,1),size(adjimb,2),0:size(adjimb,3)),   intent(out) :: adjphi

    character(len=*), parameter :: Mismatch = "mismatch in balance_simp_adj!"

    integer :: nx, ny, nz

    nx = size(adjimb,1)
    ny = size(adjimb,2)
    nz = size(adjimb,3)

    ! Sanity checks
    if (any(shape(psi) /= (/ nx+1, ny+1, nz /) )) stop "psi "  // Mismatch   
    if (any(shape(f)  /= (/ nx,   ny   /) )) stop "f "  // Mismatch
    if (any(shape(mCap) /= (/ nx, ny   /) )) stop "mCap " // Mismatch    
    if (any(shape(m)  /= (/ nx,   ny   /) )) stop "m "  // Mismatch
    if (any(shape(mu) /= (/ nx+1, ny   /) )) stop "mu " // Mismatch
    if (any(shape(mv) /= (/ nx,   ny+1 /) )) stop "mv " // Mismatch
    if (size(etaw) /= nz+1) stop "etaw " // Mismatch
    if (size(etam) /= nz)   stop "etam " // Mismatch

    call balance_lhs_adj(adjimb, psi, f, mCap, m, mu, mv, rdx, adjpsi)
    call balance_rhs_simp_adj(adjimb, mCap, m, etaw, etam, rdx, adjphi)
    adjphi = -adjphi
  end subroutine imbalance_simp_adj

  ! ===========================================================================

  subroutine imbalance_exp(psi, phi, f, mCap, m, mu, mv, etaw, etam, rdx, imb)
    real(cp), dimension(0:,0:,:), intent(in) :: psi          ! (0:nx,0:ny)
    real(cp), dimension(:,:,0:),  intent(in) :: phi
    real(cp), dimension(:,:),     intent(in) :: f, mCap, m   ! (1:nx,1:ny)
    real(cp), dimension(0:,:),    intent(in) :: mu
    real(cp), dimension(:,0:),    intent(in) :: mv           ! (1:nx,0:ny)
    real(cp), dimension(0:),      intent(in) :: etaw
    real(cp), dimension(:),       intent(in) :: etam
    real(cp),                     intent(in) :: rdx
    real(cp), dimension(size(f,1),size(f,2),size(psi,3)), intent(out) :: imb
    
    character(len=*), parameter :: Mismatch = "mismatch in imbalance_exp!"
    
    real(cp), dimension(:,:,:), allocatable :: t1, t2, t3
    real(cp), dimension(:,:),   allocatable :: c
    integer :: nx, ny, nz, i, j, k
    
    nx = size(psi,1) - 1
    ny = size(psi,2) - 1
    nz = size(psi,3)
    
    ! Sanity checks
    if (any(shape(f)  /= (/ nx,   ny   /) )) stop "f "  // Mismatch
    if (any(shape(mCap) /= (/ nx, ny /) )) stop "mCap " // Mismatch
    if (any(shape(m)  /= (/ nx,   ny   /) )) stop "m "  // Mismatch
    if (any(shape(mu) /= (/ nx+1, ny   /) )) stop "mu " // Mismatch
    if (any(shape(mv) /= (/ nx,   ny+1 /) )) stop "mv " // Mismatch
    if (any(shape(mCap) /= (/ nx, ny /) )) stop "mCap " // Mismatch
    if (size(etaw) /= nz+1) stop "etaw " // Mismatch
    if (size(etam) /= nz) stop "etam " // Mismatch
    
    allocate(t1(2:nx-1,2:ny-1,nz), t2(2:nx-1,2:ny-1,nz), t3(2:nx-1,2:ny-1,nz))
    allocate(c(nx-1,ny-1))
    
    forall (i = 1:nx-1, j = 1:ny-1)
       c(i,j) = One / (mCap(i,j) + mCap(i+1,j) + mCap(i,j+1) + mCap(i+1,j+1))
    end forall

    imb = Zero
    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       t1(i,j,k) = Eighth * rdx**2 * m(i,j)**2 * (                           &
          (f(i,  j) + f(i+1,j)) * (                                          &
            mv(i,  j-1) * (psi(i,  j-1,k) - psi(i-1,j-1,k)) +                &
            mv(i+1,j-1) * (psi(i+1,j-1,k) - psi(i,  j-1,k)) +                &
            mv(i,  j  ) * (psi(i,  j,  k) - psi(i-1,j,  k)) +                &
            mv(i+1,j  ) * (psi(i+1,j,  k) - psi(i,  j,  k)) ) / mu(i,  j) -  &
          (f(i-1,j) + f(i,j  )) * (                                          &
            mv(i-1,j-1) * (psi(i-1,j-1,k) - psi(i-2,j-1,k)) +                &
            mv(i,  j-1) * (psi(i,  j-1,k) - psi(i-1,j-1,k)) +                &
            mv(i-1,j  ) * (psi(i-1,j,  k) - psi(i-2,j,  k)) +                &
            mv(i,  j  ) * (psi(i,  j,  k) - psi(i-1,j,  k)) ) / mu(i-1,j) +  &
          (f(i,  j) + f(i,j+1)) * (                                          &
            mu(i-1,j  ) * (psi(i-1,j,  k) - psi(i-1,j-1,k)) +                &
            mu(i,  j  ) * (psi(i,  j,  k) - psi(i,  j-1,k)) +                &
            mu(i-1,j+1) * (psi(i-1,j+1,k) - psi(i-1,j,  k)) +                &
            mu(i,  j+1) * (psi(i,  j+1,k) - psi(i,  j,  k)) ) / mv(i,  j) -  &
          (f(i,j-1) + f(i,j  )) * (                                          &
            mu(i-1,j-1) * (psi(i-1,j-1,k) - psi(i-1,j-2,k)) +                &
            mu(i,  j-1) * (psi(i,  j-1,k) - psi(i,  j-2,k)) +                &
            mu(i-1,j  ) * (psi(i-1,j,  k) - psi(i-1,j-1,k)) +                &
            mu(i,  j  ) * (psi(i,  j,  k) - psi(i,  j-1,k)) ) / mv(i,j-1) )

       t2(i,j,k) = Half * rdx**4 * m(i,j) * (  &
            mu(i,j) * (mu(i,j) * (psi(i,j,k) - psi(i,j-1,k))  &
            * (mu(i,j) * (One/mCap(i+1,j)-One/mCap(i,j)) * (psi(i,j,k) - psi(i,j-1,k))  &
            + mu(i+1,j)/mCap(i+1,j) * (psi(i+1,j,k) - psi(i+1,j-1,k))  &
            - mu(i-1,j)/mCap(i,j) * (psi(i-1,j,k) - psi(i-1,j-1,k)))  &
            - (mv(i,j-1) * (psi(i,j-1,k) - psi(i-1,j-1,k))  &
            + mv(i+1,j-1) * (psi(i+1,j-1,k) - psi(i,j-1,k))  &
            + mv(i,j) * (psi(i,j,k) - psi(i-1,j,k)) + mv(i+1,j) * (psi(i+1,j,k) - psi(i,j,k)))  &
            * (c(i,j) * (mu(i,j+1) * (psi(i,j+1,k) - psi(i,j,k)) + mu(i,j)  &
            * (psi(i,j,k) - psi(i,j-1,k))) - c(i,j-1) * (mu(i,j) * (psi(i,j,k) - psi(i,j-1,k))  &
            + mu(i,j-1) * (psi(i,j-1,k) - psi(i,j-2,k)))))  &
            - mu(i-1,j) * (mu(i-1,j) * (psi(i-1,j,k) - psi(i-1,j-1,k))  &
            * (mu(i-1,j) * (One/mCap(i,j)-One/mCap(i-1,j)) * (psi(i-1,j,k) - psi(i-1,j-1,k))  &
            + mu(i,j)/mCap(i,j) * (psi(i,j,k) - psi(i,j-1,k))  &
            - mu(i-2,j)/mCap(i-1,j) * (psi(i-2,j,k) - psi(i-2,j-1,k)))  &
            - (mv(i-1,j-1) * (psi(i-1,j-1,k) - psi(i-2,j-1,k))  &
            + mv(i,j-1) * (psi(i,j-1,k) - psi(i-1,j-1,k))  &
            + mv(i-1,j) * (psi(i-1,j,k) - psi(i-2,j,k)) + mv(i,j) * (psi(i,j,k) - psi(i-1,j,k)))  &
            * (c(i-1,j) * (mu(i-1,j+1) * (psi(i-1,j+1,k) - psi(i-1,j,k))  &
            + mu(i-1,j) * (psi(i-1,j,k) - psi(i-1,j-1,k)))  &
            - c(i-1,j-1) * (mu(i-1,j) * (psi(i-1,j,k) - psi(i-1,j-1,k))  &
            + mu(i-1,j-1) * (psi(i-1,j-1,k) - psi(i-1,j-2,k))))))

       t3(i,j,k) = Half * rdx**4 * m(i,j) * (  &
            mv(i,j) * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) * (One/mCap(i,j+1)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i,j+1) * (psi(i,j+1,k)-psi(i-1,j+1,k)))  &
            - One/mCap(i,j) * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))  &
            + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)))) - (mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))  &
            + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k)) + mu(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k))  &
            + mu(i,j+1) * (psi(i,j+1,k)-psi(i,j,k))) * (c(i,j)  &
            * (mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)) + mv(i+1,j) * (psi(i+1,j,k)-psi(i,j,k)))  &
            - c(i-1,j) * (mv(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))  &
            + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)))))  &
            - mv(i,j-1) * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)) * (One/mCap(i,j)  &
            * (mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i,j) * (psi(i,j,k)-psi(i-1,j,k)))  &
            - One/mCap(i,j-1) * (mv(i,j-2) * (psi(i,j-2,k)-psi(i-1,j-2,k))  &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))))  &
            - (mu(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-1,j-2,k))  &
            + mu(i,j-1) * (psi(i,j-1,k)-psi(i,j-2,k)) + mu(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k)) &
            + mu(i,j) * (psi(i,j,k)-psi(i,j-1,k))) * (c(i,j-1) * (mv(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + mv(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k)))  &
            - c(i-1,j-1) * (mv(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-2,j-1,k))  &
            + mv(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))))))
    end forall
    imb(2:nx-1,2:ny-1,:) = t1 - t2 - t3
    deallocate(c, t3, t2, t1)

    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       imb(i,j,k) = imb(i,j,k) - m(i,j)**2 * Half * rdx**2 * (Half  &
            * ((phi(i+1,j,k-1)-phi(i,j,k-1) + phi(i+1,j,k)-phi(i,j,k)) * (mCap(i,j) + mCap(i+1,j))&
            - (phi(i,j,k-1)-phi(i-1,j,k-1) + phi(i,j,k)-phi(i-1,j,k)) * (mCap(i-1,j) + mCap(i,j)) &
            + (phi(i,j+1,k-1)-phi(i,j,k-1) + phi(i,j+1,k)-phi(i,j,k)) * (mCap(i,j) + mCap(i,j+1)) &
            - (phi(i,j,k-1)-phi(i,j-1,k-1) + phi(i,j,k)-phi(i,j-1,k)) * (mCap(i,j-1) + mCap(i,j)))&
            + etam(k) / (etaw(k) - etaw(k-1)) * (-(mCap(i+1,j) - mCap(i,j))  &
            * (phi(i,j,k)-phi(i,j,k-1) + phi(i+1,j,k)-phi(i+1,j,k-1)) + (mCap(i,j) - mCap(i-1,j)) &
            * (phi(i-1,j,k)-phi(i-1,j,k-1) + phi(i,j,k)-phi(i,j,k-1)) - (mCap(i,j+1) - mCap(i,j)) &
            * (phi(i,j,k)-phi(i,j,k-1) + phi(i,j+1,k)-phi(i,j+1,k-1)) + (mCap(i,j) - mCap(i,j-1)) &
            * (phi(i,j-1,k)-phi(i,j-1,k-1) + phi(i,j,k)-phi(i,j,k-1))))
    end forall
  end subroutine imbalance_exp

  ! ===========================================================================

  subroutine pv_exp(psi, phi, f, mCap, m, mu, mv, etaz, etam, pt, rdx, q)
    real(cp), dimension(0:,0:,:), intent(in) :: psi          ! (0:nx,0:ny)
    real(cp), dimension(:,:,0:),  intent(in) :: phi
    real(cp), dimension(:,:),     intent(in) :: f, mCap, m   ! (1:nx,1:ny)
    real(cp), dimension(0:,:),    intent(in) :: mu
    real(cp), dimension(:,0:),    intent(in) :: mv           ! (1:nx,0:ny)
    real(cp), dimension(0:),      intent(in) :: etaz
    real(cp), dimension(:),       intent(in) :: etam
    real(cp),                     intent(in) :: pt, rdx
    real, dimension(size(f,1),size(f,2),size(psi,3)), intent(out) :: q

    real(cp), parameter :: p0 = 100000._cp
    character(len=*), parameter :: Mismatch = "mismatch in pv_exp!"
    
    real(cp), dimension(:,:,:), allocatable :: e
    real(cp), dimension(:,:),   allocatable :: mpsi, c, d
    integer :: nx, ny, nz, i, j, k
    
    nx = size(psi,1) - 1
    ny = size(psi,2) - 1
    nz = size(psi,3)

    ! Sanity checks
    if (any(shape(phi) /= (/ nx,  ny, nz+1 /) )) stop "phi " // Mismatch   
    if (any(shape(f)  /= (/ nx,   ny   /) )) stop "f "  // Mismatch
    if (any(shape(mCap) /= (/ nx, ny   /) )) stop "mCap " // Mismatch    
    if (any(shape(m)  /= (/ nx,   ny   /) )) stop "m "  // Mismatch
    if (any(shape(mu) /= (/ nx+1, ny   /) )) stop "mu " // Mismatch
    if (any(shape(mv) /= (/ nx,   ny+1 /) )) stop "mv " // Mismatch
    if (size(etaz) /= nz+1) stop "etaz " // Mismatch
    if (size(etam) /= nz) stop "etam " // Mismatch

    allocate(e(2:nx-1,2:ny-1,nz))
    allocate(mpsi(nx-1,ny-1), c(nx-1,ny), d(nx,ny-1))

    forall (i = 1:nx-1, j = 1:ny-1)
       mpsi(i,j) = Quarter * (mu(i,j) + mu(i,j+1) + mv(i,j) + mv(i+1,j))
    end forall
    forall (i = 1:nx-1, j = 1:ny)
       c(i,j) = mu(i,j) / (mCap(i,j) + mCap(i+1,j))
    end forall
    forall (i = 1:nx, j = 1:ny-1)
       d(i,j) = mv(i,j) / (mCap(i,j) + mCap(i,j+1))
    end forall
    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       e(i,j,k) = (p0 / (pt + etam(k) * m(i,j)))**Kappa  &
            / (RGas * log((pt + etaz(k-1) * m(i,j)) / (pt + etaz(k) * m(i,j))))
    end forall

    forall (i = 2:nx-1, j = 2:ny-1, k = 2:nz-1)
       q(i,j,k) = Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1)))  &
            * ((phi(i,j,k-1)-phi(i,j,k-2)) * e(i,j,k-1) - (phi(i,j,k+1)-phi(i,j,k)) * e(i,j,k+1)) &
            * (f(i,j) + Half * rdx**2 * (mpsi(i-1,j-1) * (d(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))&
            - d(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-2,j-1,k)) + c(i-1,j)*(psi(i-1,j,k)-psi(i-1,j-1,k)) &
            - c(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-1,j-2,k))) + mpsi(i,j-1)  &
            * (d(i+1,j-1)*(psi(i+1,j-1,k)-psi(i,j-1,k)) - d(i,j-1)*(psi(i,j-1,k)-psi(i-1,j-1,k))  &
            + c(i,j) * (psi(i,j,k)-psi(i,j-1,k)) - c(i,j-1) * (psi(i,j-1,k)-psi(i,j-2,k)))  &
            + mpsi(i-1,j) * (d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            - d(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k)) + c(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k)) &
            - c(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) + mpsi(i,j)  &
            * (d(i+1,j) * (psi(i+1,j,k)-psi(i,j,k)) - d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            + c(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) - c(i,j) * (psi(i,j,k)-psi(i,j-1,k)))  &
            - Quarter / (phi(i,j,k)-phi(i,j,k-1))  &
            * ((mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k))  &
            + mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)))  &
            * (d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1)) &
            + d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1)))  &
            + (mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k))  &
            + mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k)))  &
            * (c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1)) &
            + c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))))))
    end forall
    q((/1,nx/),2:ny-1,2:nz-1) = q((/2,nx-1/),2:ny-1,2:nz-1)
    q(:,(/1,ny/),2:nz-1) = q(:,(/2,ny-1/),2:nz-1)
    q(:,:,(/1,nz/)) = q(:,:,(/2,nz-1/))
  end subroutine pv_exp

  ! ===========================================================================

  subroutine pv_tlm(dps, dph, psi, phi, f, mCap, m, mu, mv, etaz, etam, pt, rdx, dpv)
    real(cp), dimension(0:,0:,:), intent(in) :: dps, psi          ! (0:nx,0:ny)
    real(cp), dimension(:,:,0:),  intent(in) :: dph, phi
    real(cp), dimension(:,:),     intent(in) :: f, mCap, m   ! (1:nx,1:ny)
    real(cp), dimension(0:,:),    intent(in) :: mu
    real(cp), dimension(:,0:),    intent(in) :: mv           ! (1:nx,0:ny)
    real(cp), dimension(0:),      intent(in) :: etaz
    real(cp), dimension(:),       intent(in) :: etam
    real(cp),                     intent(in) :: pt, rdx
    real, dimension(size(f,1),size(f,2),size(psi,3)), intent(out) :: dpv

    real(cp), parameter :: p0 = 100000._cp
    character(len=*), parameter :: Mismatch = "mismatch in pv_tlm!"
    
    real(cp), dimension(:,:,:), allocatable :: dpsi, dphi, e
    real(cp), dimension(:,:),   allocatable :: mpsi, c, d
    integer :: nx, ny, nz, i, j, k
    
    nx = size(psi,1) - 1
    ny = size(psi,2) - 1
    nz = size(psi,3)

    ! Sanity checks
    if (any(shape(dps) /= (/ nx+1,  ny+1, nz /) )) stop "dps " // Mismatch     
    if (any(shape(psi) /= (/ nx+1,  ny+1, nz /) )) stop "psi " // Mismatch    
    if (any(shape(dph) /= (/ nx,  ny, nz+1 /) )) stop "dph " // Mismatch     
    if (any(shape(phi) /= (/ nx,  ny, nz+1 /) )) stop "phi " // Mismatch   
    if (any(shape(f)  /= (/ nx,   ny   /) )) stop "f "  // Mismatch
    if (any(shape(mCap) /= (/ nx, ny   /) )) stop "mCap " // Mismatch    
    if (any(shape(m)  /= (/ nx,   ny   /) )) stop "m "  // Mismatch
    if (any(shape(mu) /= (/ nx+1, ny   /) )) stop "mu " // Mismatch
    if (any(shape(mv) /= (/ nx,   ny+1 /) )) stop "mv " // Mismatch
    if (size(etaz) /= nz+1) stop "etaz " // Mismatch
    if (size(etam) /= nz) stop "etam " // Mismatch

    allocate(dpsi(0:nx,0:ny,nz), dphi(nx,ny,0:nz), e(2:nx-1,2:ny-1,nz))
    allocate(mpsi(nx-1,ny-1), c(nx-1,ny), d(nx,ny-1))

    dpsi = dps
    dpsi((/0,nx/),:,:) = Zero
    dpsi(1:nx-1,(/0,ny/),:) = Zero
    dphi = dph
    dphi((/1,nx/),:,:) = Zero
    dphi(2:nx-1,(/1,ny/),:) = Zero
    dphi(2:nx-1,2:ny-1,(/0,nz/)) = Zero
    dpv = Zero   ! Really just setting the boundaries.

    forall (i = 1:nx-1, j = 1:ny-1)
       mpsi(i,j) = Quarter * (mu(i,j) + mu(i,j+1) + mv(i,j) + mv(i+1,j))
    end forall
    forall (i = 1:nx-1, j = 1:ny)
       c(i,j) = mu(i,j) / (mCap(i,j) + mCap(i+1,j))
    end forall
    forall (i = 1:nx, j = 1:ny-1)
       d(i,j) = mv(i,j) / (mCap(i,j) + mCap(i,j+1))
    end forall
    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       e(i,j,k) = (p0 / (pt + etam(k) * m(i,j)))**Kappa  &
            / (RGas * log((pt + etaz(k-1) * m(i,j)) / (pt + etaz(k) * m(i,j))))
    end forall

    forall (i = 2:nx-1, j = 2:ny-1, k = 2:nz-1)
       dpv(i,j,k) = Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dphi(i,j,k-2) * e(i,j,k-1) * (f(i,j) + Half * rdx**2 * (mpsi(i-1,j-1)  &
            * (d(i,j-1)*(psi(i,j-1,k)-psi(i-1,j-1,k)) - d(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-2,j-1,k))&
            + c(i-1,j)*(psi(i-1,j,k)-psi(i-1,j-1,k)) - c(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-1,j-2,k)))&
            + mpsi(i,j-1) * (d(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k)) - d(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + c(i,j) * (psi(i,j,k)-psi(i,j-1,k)) - c(i,j-1)  &
            * (psi(i,j-1,k)-psi(i,j-2,k))) + mpsi(i-1,j) * (d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            - d(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k)) + c(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k)) &
            - c(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) + mpsi(i,j)  &
            * (d(i+1,j) * (psi(i+1,j,k)-psi(i,j,k)) - d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            + c(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) - c(i,j) * (psi(i,j,k)-psi(i,j-1,k)))  &
            - Quarter / (phi(i,j,k)-phi(i,j,k-1))  &
            * ((mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k))  &
            + mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)))  &
            * (d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1)) &
            + d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1)))  &
            + (mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k))  &
            + mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k)))  &
            * (c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1)) &
            + c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1)))))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dphi(i,j-1,k-1) * Eighth * rdx**2 * mv(i,j-1) / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (c(i,j) * (-psi(i,j,k-1) + psi(i,j-1,k-1) + psi(i,j,k+1) - psi(i,j-1,k+1))  &
            + c(i-1,j) * (-psi(i-1,j,k-1) + psi(i-1,j-1,k-1) + psi(i-1,j,k+1) - psi(i-1,j-1,k+1))))
       
       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dphi(i-1,j,k-1) * Eighth * rdx**2 * mu(i-1,j) / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j) * (-psi(i,j,k-1) + psi(i-1,j,k-1) + psi(i,j,k+1) - psi(i-1,j,k+1))  &
            + d(i,j-1) * (-psi(i,j-1,k-1) + psi(i-1,j-1,k-1) + psi(i,j-1,k+1) - psi(i-1,j-1,k+1))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dphi(i,j,k-1) * (e(i,j,k-1) * (f(i,j) + Half * rdx**2 * (mpsi(i-1,j-1)  &
            * (d(i,j-1)*(psi(i,j-1,k)-psi(i-1,j-1,k)) - d(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-2,j-1,k))&
            + c(i-1,j)*(psi(i-1,j,k)-psi(i-1,j-1,k)) - c(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-1,j-2,k)))&
            + mpsi(i,j-1) * (d(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k)) - d(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + c(i,j) * (psi(i,j,k)-psi(i,j-1,k)) - c(i,j-1)  &
            * (psi(i,j-1,k)-psi(i,j-2,k))) + mpsi(i-1,j) * (d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            - d(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k)) + c(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k)) &
            - c(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) + mpsi(i,j)  &
            * (d(i+1,j) * (psi(i+1,j,k)-psi(i,j,k)) - d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            + c(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) - c(i,j) * (psi(i,j,k)-psi(i,j-1,k)))  &
            - Quarter / (phi(i,j,k)-phi(i,j,k-1))  &
            * ((mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k))  &
            + mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)))  &
            * (d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1)) &
            + d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1)))  &
            + (mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k))  &
            + mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k)))  &
            * (c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1)) &
            + c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))))))  &
            + Half * rdx**2 * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (-Quarter / (phi(i,j,k)-phi(i,j,k-1))**2  &
            * ((mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            * (d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1))  &
            + d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1))) &
            + (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))  &
            * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))  &
            + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1))))&
            - Quarter / (phi(i,j,k)-phi(i,j,k-1)) * ((mu(i-1,j) - mu(i,j))  &
            * (d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1))  &
            + d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1))) &
            + (mv(i,j-1) - mv(i,j)) * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1)  &
            + psi(i,j-1,k-1)) + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1)  &
            + psi(i-1,j-1,k-1)))))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dphi(i+1,j,k-1) * Eighth * rdx**2 * mu(i,j) / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j) * (-psi(i,j,k-1) + psi(i-1,j,k-1) + psi(i,j,k+1) - psi(i-1,j,k+1))  &
            + d(i,j-1) * (-psi(i,j-1,k-1) + psi(i-1,j-1,k-1) + psi(i,j-1,k+1) - psi(i-1,j-1,k+1))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dphi(i,j+1,k-1) * Eighth * rdx**2 * mv(i,j) / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))  &
            + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dphi(i,j-1,k) * Eighth * rdx**2 * mv(i,j-1) / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (c(i,j) * (-psi(i,j,k-1) + psi(i,j-1,k-1) + psi(i,j,k+1) - psi(i,j-1,k+1))  &
            + c(i-1,j) * (-psi(i-1,j,k-1) + psi(i-1,j-1,k-1) + psi(i-1,j,k+1) - psi(i-1,j-1,k+1))))
 
       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dphi(i-1,j,k) * Eighth * rdx**2 * mu(i-1,j) / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j) * (-psi(i,j,k-1) + psi(i-1,j,k-1) + psi(i,j,k+1) - psi(i-1,j,k+1))  &
            + d(i,j-1) * (-psi(i,j-1,k-1) + psi(i-1,j-1,k-1) + psi(i,j-1,k+1) - psi(i-1,j-1,k+1))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dphi(i,j,k) * (e(i,j,k+1) * (f(i,j) + Half * rdx**2 * (mpsi(i-1,j-1)  &
            * (d(i,j-1)*(psi(i,j-1,k)-psi(i-1,j-1,k)) - d(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-2,j-1,k))&
            + c(i-1,j)*(psi(i-1,j,k)-psi(i-1,j-1,k)) - c(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-1,j-2,k)))&
            + mpsi(i,j-1) * (d(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k)) - d(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + c(i,j) * (psi(i,j,k)-psi(i,j-1,k)) - c(i,j-1)  &
            * (psi(i,j-1,k)-psi(i,j-2,k))) + mpsi(i-1,j) * (d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            - d(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k)) + c(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k)) &
            - c(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) + mpsi(i,j)  &
            * (d(i+1,j) * (psi(i+1,j,k)-psi(i,j,k)) - d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            + c(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) - c(i,j) * (psi(i,j,k)-psi(i,j-1,k)))  &
            - Quarter / (phi(i,j,k)-phi(i,j,k-1))  &
            * ((mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k))  &
            + mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)))  &
            * (d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1)) &
            + d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1)))  &
            + (mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k))  &
            + mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k)))  &
            * (c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1)) &
            + c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))))))  &
            + Half * rdx**2 * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (Quarter / (phi(i,j,k)-phi(i,j,k-1))**2  &
            * ((mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            * (d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1))  &
            + d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1))) &
            + (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))  &
            * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))  &
            + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1))))&
            - Quarter / (phi(i,j,k)-phi(i,j,k-1)) * ((mu(i-1,j) - mu(i,j))  &
            * (d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1))  &
            + d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1))) &
            + (mv(i,j-1) - mv(i,j)) * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1)  &
            + psi(i,j-1,k-1)) + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1)  &
            + psi(i-1,j-1,k-1)))))))
       
       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dphi(i+1,j,k) * Eighth * rdx**2 * mu(i,j) / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j) * (-psi(i,j,k-1) + psi(i-1,j,k-1) + psi(i,j,k+1) - psi(i-1,j,k+1))  &
            + d(i,j-1) * (-psi(i,j-1,k-1) + psi(i-1,j-1,k-1) + psi(i,j-1,k+1) - psi(i-1,j-1,k+1))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dphi(i,j+1,k) * Eighth * rdx**2 * mv(i,j) / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))  &
            + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1))))
       
       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dphi(i,j,k+1) * e(i,j,k+1) * (f(i,j) + Half * rdx**2 * (mpsi(i-1,j-1)  &
            * (d(i,j-1)*(psi(i,j-1,k)-psi(i-1,j-1,k)) - d(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-2,j-1,k))&
            + c(i-1,j)*(psi(i-1,j,k)-psi(i-1,j-1,k)) - c(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-1,j-2,k)))&
            + mpsi(i,j-1) * (d(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k)) - d(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + c(i,j) * (psi(i,j,k)-psi(i,j-1,k)) - c(i,j-1)  &
            * (psi(i,j-1,k)-psi(i,j-2,k))) + mpsi(i-1,j) * (d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            - d(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k)) + c(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k)) &
            - c(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) + mpsi(i,j)  &
            * (d(i+1,j) * (psi(i+1,j,k)-psi(i,j,k)) - d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            + c(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) - c(i,j) * (psi(i,j,k)-psi(i,j-1,k)))  &
            - Quarter / (phi(i,j,k)-phi(i,j,k-1))  &
            * ((mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k))  &
            + mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)))  &
            * (d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1)) &
            + d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1)))  &
            + (mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k))  &
            + mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k)))  &
            * (c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1)) &
            + c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1)))))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dpsi(i-1,j-1,k-1) * Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j-1) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            + c(i-1,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dpsi(i,j-1,k-1) * Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (c(i,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))  &
            - d(i,j-1) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dpsi(i-1,j,k-1) * Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            - c(i-1,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dpsi(i,j,k-1) * Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (-d(i,j) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            - c(i,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i-1,j-2,k) * Half * rdx**2 * mpsi(i-1,j-1) * c(i-1,j-1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i,j-2,k) * Half * rdx**2 * mpsi(i,j-1) * c(i,j-1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i-2,j-1,k) * Half * rdx**2 * mpsi(i-1,j-1) * d(i-1,j-1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i-1,j-1,k) * Half * rdx**2 * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2))  &
            - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) * (mpsi(i,j-1) * d(i,j-1) + mpsi(i-1,j-1)  &
            * (-d(i,j-1) - d(i-1,j-1) - c(i-1,j) - c(i-1,j-1)) + mpsi(i-1,j) * c(i-1,j)))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i,j-1,k) * Half * rdx**2 * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2))  &
            - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) * (mpsi(i,j-1) * (-d(i+1,j-1) - d(i,j-1)  &
            - c(i,j) - c(i,j-1)) + mpsi(i-1,j-1) * d(i,j-1) + mpsi(i,j) * c(i,j)))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i+1,j-1,k) * Half * rdx**2 * mpsi(i,j-1) * d(i+1,j-1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i-2,j,k) * Half * rdx**2 * mpsi(i-1,j) * d(i-1,j)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i-1,j,k) * Half * rdx**2 * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2))  &
            - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) * (mpsi(i,j) * d(i,j) + mpsi(i-1,j)  &
            * (-d(i,j) - d(i-1,j) - c(i-1,j+1) - c(i-1,j)) + mpsi(i-1,j-1) * c(i-1,j)))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i,j,k) * Half * rdx**2 * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2))  &
            - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) * (mpsi(i,j) * (-d(i+1,j) - d(i,j)  &
            - c(i,j+1) - c(i,j)) + mpsi(i-1,j) * d(i,j) + mpsi(i,j-1) * c(i,j)))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i+1,j,k) * Half * rdx**2 * mpsi(i,j) * d(i+1,j)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))
    
       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i-1,j+1,k) * Half * rdx**2 * mpsi(i-1,j) * c(i-1,j+1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))
  
       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            dpsi(i,j+1,k) * Half * rdx**2 * mpsi(i,j) * c(i,j+1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dpsi(i-1,j-1,k+1) * Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (-d(i,j-1) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)) &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            - c(i-1,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dpsi(i,j-1,k+1) * Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j-1) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)) &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            - c(i,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dpsi(i-1,j,k+1) * Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (c(i-1,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))  &
            - d(i,j) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))))

       dpv(i,j,k) = dpv(i,j,k) + Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1))) * (  &
            -dpsi(i,j,k+1) * Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            + c(i,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))
    end forall
    deallocate(d, c, mpsi, e, dphi, dpsi)
  end subroutine pv_tlm

  ! ===========================================================================

  subroutine pv_adj(adjpv, psi, phi, f, mCap, m, mu, mv, etaz, etam, pt, rdx, ps, ph)
    real(cp), dimension(:,:,:),   intent(in) :: adjpv
    real(cp), dimension(0:,0:,:), intent(in) :: psi          ! (0:nx,0:ny)
    real(cp), dimension(:,:,0:),  intent(in) :: phi
    real(cp), dimension(:,:),     intent(in) :: f, mCap, m   ! (1:nx,1:ny)
    real(cp), dimension(0:,:),    intent(in) :: mu
    real(cp), dimension(:,0:),    intent(in) :: mv           ! (1:nx,0:ny)
    real(cp), dimension(0:),      intent(in) :: etaz
    real(cp), dimension(:),       intent(in) :: etam
    real(cp),                     intent(in) :: pt, rdx
    real(cp), dimension(0:size(phi,1),0:size(phi,2),size(psi,3)), intent(out) :: ps
    real(cp), dimension(size(phi,1),size(phi,2),0:size(psi,3)),   intent(out) :: ph    

    real(cp), parameter :: p0 = 100000._cp
    character(len=*), parameter :: Mismatch = "mismatch in pv_adj!"
    
    real(cp), dimension(:,:,:), allocatable :: b, e
    real(cp), dimension(:,:),   allocatable :: mpsi, c, d
    integer :: nx, ny, nz, i, j, k
    
    nx = size(psi,1) - 1
    ny = size(psi,2) - 1
    nz = size(psi,3)

    ! Sanity checks
    if (any(shape(adjpv) /= (/ nx, ny, nz /) )) stop "adjpv " // Mismatch     
    if (any(shape(psi) /= (/ nx+1,  ny+1, nz /) )) stop "psi " // Mismatch    
    if (any(shape(phi) /= (/ nx,  ny, nz+1 /) )) stop "phi " // Mismatch   
    if (any(shape(f)  /= (/ nx,   ny   /) )) stop "f "  // Mismatch
    if (any(shape(mCap) /= (/ nx, ny   /) )) stop "mCap " // Mismatch    
    if (any(shape(m)  /= (/ nx,   ny   /) )) stop "m "  // Mismatch
    if (any(shape(mu) /= (/ nx+1, ny   /) )) stop "mu " // Mismatch
    if (any(shape(mv) /= (/ nx,   ny+1 /) )) stop "mv " // Mismatch
    if (size(etaz) /= nz+1) stop "etaz " // Mismatch
    if (size(etam) /= nz) stop "etam " // Mismatch

    allocate(b(2:nx-1,2:ny-1,2:nz-1), e(2:nx-1,2:ny-1,nz))
    allocate(mpsi(nx-1,ny-1), c(nx-1,ny), d(nx,ny-1))

    ps = 0.
    ph = 0.
    forall (i = 1:nx-1, j = 1:ny-1)
       mpsi(i,j) = Quarter * (mu(i,j) + mu(i,j+1) + mv(i,j) + mv(i+1,j))
    end forall
    forall (i = 2:nx-1, j = 2:ny-1, k = 2:nz-1)
       b(i,j,k) = Half * Grav / (m(i,j) * (etaz(k)-etaz(k-1)))
    end forall
    forall (i = 1:nx-1, j = 1:ny)
       c(i,j) = mu(i,j) / (mCap(i,j) + mCap(i+1,j))
    end forall
    forall (i = 1:nx, j = 1:ny-1)
       d(i,j) = mv(i,j) / (mCap(i,j) + mCap(i,j+1))
    end forall
    forall (i = 2:nx-1, j = 2:ny-1, k = 1:nz)
       e(i,j,k) = (p0 / (pt + etam(k) * m(i,j)))**Kappa  &
            / (RGas * log((pt + etaz(k-1) * m(i,j)) / (pt + etaz(k) * m(i,j))))
    end forall

    forall (i = 2:nx-1, j = 2:ny-1, k = 2:nz-1)
       ph(i,j,k-2) = ph(i,j,k-2) - b(i,j,k) * (e(i,j,k-1) * (f(i,j) + Half * rdx**2  &
            * (mpsi(i-1,j-1) * (d(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k))  &
            - d(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-2,j-1,k)) + c(i-1,j)*(psi(i-1,j,k)-psi(i-1,j-1,k)) &
            - c(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-1,j-2,k))) + mpsi(i,j-1)  &
            * (d(i+1,j-1)*(psi(i+1,j-1,k)-psi(i,j-1,k)) - d(i,j-1)*(psi(i,j-1,k)-psi(i-1,j-1,k))  &
            + c(i,j) * (psi(i,j,k)-psi(i,j-1,k)) - c(i,j-1) * (psi(i,j-1,k)-psi(i,j-2,k)))  &
            + mpsi(i-1,j)*(d(i,j)*(psi(i,j,k)-psi(i-1,j,k))-d(i-1,j)*(psi(i-1,j,k)-psi(i-2,j,k))  &
            + c(i-1,j+1)*(psi(i-1,j+1,k)-psi(i-1,j,k)) - c(i-1,j)*(psi(i-1,j,k)-psi(i-1,j-1,k)))  &
            + mpsi(i,j) * (d(i+1,j)*(psi(i+1,j,k)-psi(i,j,k)) - d(i,j)*(psi(i,j,k)-psi(i-1,j,k))  &
            + c(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) - c(i,j) * (psi(i,j,k)-psi(i,j-1,k)))  &
            - Quarter / (phi(i,j,k)-phi(i,j,k-1))  &
            * ((mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k))  &
            + mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)))  &
            * (d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1)) &
            + d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1)))  &
            + (mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k))  &
            + mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k)))  &
            * (c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1)) &
            + c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1)))))))  &
            * adjpv(i,j,k)

       ph(i,j-1,k-1) = ph(i,j-1,k-1) + b(i,j,k) * (Eighth * rdx**2 * mv(i,j-1)  &
            / (phi(i,j,k)-phi(i,j,k-1)) * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (c(i,j) * (-psi(i,j,k-1) + psi(i,j-1,k-1)  &
            + psi(i,j,k+1) - psi(i,j-1,k+1)) + c(i-1,j) * (-psi(i-1,j,k-1) + psi(i-1,j-1,k-1)  &
            + psi(i-1,j,k+1) - psi(i-1,j-1,k+1)))) * adjpv(i,j,k)
       ph(i-1,j,k-1) = ph(i-1,j,k-1) + b(i,j,k) * (Eighth * rdx**2 * mu(i-1,j)  &
            / (phi(i,j,k)-phi(i,j,k-1)) * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (d(i,j) * (-psi(i,j,k-1) + psi(i-1,j,k-1)  &
            + psi(i,j,k+1) - psi(i-1,j,k+1)) + d(i,j-1) * (-psi(i,j-1,k-1) + psi(i-1,j-1,k-1)  &
            + psi(i,j-1,k+1) - psi(i-1,j-1,k+1)))) * adjpv(i,j,k)

       ph(i,j,k-1) = ph(i,j,k-1) + b(i,j,k) * (e(i,j,k-1) * (f(i,j) + Half * rdx**2  &
            * (mpsi(i-1,j-1) * (d(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)) - d(i-1,j-1)  &
            * (psi(i-1,j-1,k)-psi(i-2,j-1,k)) + c(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))  &
            - c(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-1,j-2,k))) + mpsi(i,j-1) * (d(i+1,j-1)  &
            * (psi(i+1,j-1,k)-psi(i,j-1,k)) - d(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)) + c(i,j)  &
            * (psi(i,j,k)-psi(i,j-1,k)) - c(i,j-1) * (psi(i,j-1,k)-psi(i,j-2,k))) + mpsi(i-1,j)  &
            * (d(i,j) * (psi(i,j,k)-psi(i-1,j,k)) - d(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))  &
            + c(i-1,j+1)*(psi(i-1,j+1,k)-psi(i-1,j,k)) - c(i-1,j)*(psi(i-1,j,k)-psi(i-1,j-1,k))) &
            + mpsi(i,j) * (d(i+1,j)*(psi(i+1,j,k)-psi(i,j,k)) - d(i,j)*(psi(i,j,k)-psi(i-1,j,k))  &
            + c(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) - c(i,j) * (psi(i,j,k)-psi(i,j-1,k)))  &
            - Quarter / (phi(i,j,k)-phi(i,j,k-1))  &
            * ((mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k))  &
            + mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)))  &
            * (d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1)) &
            + d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1)))  &
            + (mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k))  &
            + mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k)))  &
            * (c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1)) &
            + c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))))))  &
            + Half * rdx**2 * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (-Quarter / (phi(i,j,k)-phi(i,j,k-1))**2  &
            * ((mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            * (d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1))  &
            + d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1))) &
            + (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))  &
            * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))  &
            + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1))))&
            - Quarter / (phi(i,j,k)-phi(i,j,k-1)) * ((mu(i-1,j) - mu(i,j))  &
            * (d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1))  &
            + d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1))) &
            + (mv(i,j-1) - mv(i,j)) * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1)  &
            + psi(i,j-1,k-1)) + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1)  &
            + psi(i-1,j-1,k-1)))))) * adjpv(i,j,k)

       ph(i+1,j,k-1) = ph(i+1,j,k-1) - b(i,j,k) * (Eighth * rdx**2 * mu(i,j)  &
            / (phi(i,j,k)-phi(i,j,k-1)) * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (d(i,j) * (-psi(i,j,k-1) + psi(i-1,j,k-1)  &
            + psi(i,j,k+1) - psi(i-1,j,k+1)) + d(i,j-1) * (-psi(i,j-1,k-1) + psi(i-1,j-1,k-1)  &
            + psi(i,j-1,k+1) - psi(i-1,j-1,k+1)))) * adjpv(i,j,k)
       ph(i,j+1,k-1) = ph(i,j+1,k-1) - b(i,j,k) * (Eighth * rdx**2 * mv(i,j)  &
            / (phi(i,j,k)-phi(i,j,k-1)) * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1)&
            + psi(i,j-1,k-1)) + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1)  &
            + psi(i-1,j-1,k-1)))) * adjpv(i,j,k)
       ph(i,j-1,k) = ph(i,j-1,k) + b(i,j,k) * (Eighth * rdx**2 * mv(i,j-1)  &
            / (phi(i,j,k)-phi(i,j,k-1)) * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (c(i,j) * (-psi(i,j,k-1) + psi(i,j-1,k-1)  &
            + psi(i,j,k+1) - psi(i,j-1,k+1)) + c(i-1,j) * (-psi(i-1,j,k-1) + psi(i-1,j-1,k-1)  &
            + psi(i-1,j,k+1) - psi(i-1,j-1,k+1)))) * adjpv(i,j,k)
       ph(i-1,j,k) = ph(i-1,j,k) + b(i,j,k) * (Eighth * rdx**2 * mu(i-1,j)  &
            / (phi(i,j,k)-phi(i,j,k-1)) * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (d(i,j) * (-psi(i,j,k-1) + psi(i-1,j,k-1)  &
            + psi(i,j,k+1) - psi(i-1,j,k+1)) + d(i,j-1) * (-psi(i,j-1,k-1) + psi(i-1,j-1,k-1)  &
            + psi(i,j-1,k+1) - psi(i-1,j-1,k+1)))) * adjpv(i,j,k)

       ph(i,j,k) = ph(i,j,k) + b(i,j,k) * (e(i,j,k+1) * (f(i,j) + Half * rdx**2 * (mpsi(i-1,j-1)  &
            * (d(i,j-1)*(psi(i,j-1,k)-psi(i-1,j-1,k)) - d(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-2,j-1,k))&
            + c(i-1,j)*(psi(i-1,j,k)-psi(i-1,j-1,k)) - c(i-1,j-1)*(psi(i-1,j-1,k)-psi(i-1,j-2,k)))&
            + mpsi(i,j-1) * (d(i+1,j-1) * (psi(i+1,j-1,k)-psi(i,j-1,k)) - d(i,j-1)  &
            * (psi(i,j-1,k)-psi(i-1,j-1,k)) + c(i,j) * (psi(i,j,k)-psi(i,j-1,k)) - c(i,j-1)  &
            * (psi(i,j-1,k)-psi(i,j-2,k))) + mpsi(i-1,j) * (d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            - d(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k)) + c(i-1,j+1) * (psi(i-1,j+1,k)-psi(i-1,j,k)) &
            - c(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))) + mpsi(i,j)  &
            * (d(i+1,j) * (psi(i+1,j,k)-psi(i,j,k)) - d(i,j) * (psi(i,j,k)-psi(i-1,j,k))  &
            + c(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) - c(i,j) * (psi(i,j,k)-psi(i,j-1,k)))  &
            - Quarter / (phi(i,j,k)-phi(i,j,k-1))  &
            * ((mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k))  &
            + mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)))  &
            * (d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1)) &
            + d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1)))  &
            + (mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k))  &
            + mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k)))  &
            * (c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1)) &
            + c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))))))  &
            + Half * rdx**2 * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (Quarter / (phi(i,j,k)-phi(i,j,k-1))**2  &
            * ((mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            * (d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1))  &
            + d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1))) &
            + (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))  &
            * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1))  &
            + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1))))&
            - Quarter / (phi(i,j,k)-phi(i,j,k-1)) * ((mu(i-1,j) - mu(i,j))  &
            * (d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1))  &
            + d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1))) &
            + (mv(i,j-1) - mv(i,j)) * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1)  &
            + psi(i,j-1,k-1)) + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1)  &
            + psi(i-1,j-1,k-1)))))) * adjpv(i,j,k)

       ph(i+1,j,k) = ph(i+1,j,k) - b(i,j,k) * (Eighth * rdx**2 * mu(i,j)  &
            / (phi(i,j,k)-phi(i,j,k-1)) * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (d(i,j) * (-psi(i,j,k-1) + psi(i-1,j,k-1)  &
            + psi(i,j,k+1) - psi(i-1,j,k+1)) + d(i,j-1) * (-psi(i,j-1,k-1) + psi(i-1,j-1,k-1)  &
            + psi(i,j-1,k+1) - psi(i-1,j-1,k+1)))) * adjpv(i,j,k)
       ph(i,j+1,k) = ph(i,j+1,k) - b(i,j,k) * (Eighth * rdx**2 * mv(i,j)  &
            / (phi(i,j,k)-phi(i,j,k-1)) * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1)  &
            * (phi(i,j,k+1)-phi(i,j,k))) * (c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1)  &
            - psi(i,j,k-1) + psi(i,j-1,k-1)) + c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1)  &
            - psi(i-1,j,k-1) + psi(i-1,j-1,k-1)))) * adjpv(i,j,k)

       ph(i,j,k+1) = ph(i,j,k+1) - b(i,j,k) * (e(i,j,k+1) * (f(i,j) + Half * rdx**2  &
            * (mpsi(i-1,j-1) * (d(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)) - d(i-1,j-1)  &
            * (psi(i-1,j-1,k)-psi(i-2,j-1,k)) + c(i-1,j) * (psi(i-1,j,k)-psi(i-1,j-1,k))  &
            - c(i-1,j-1) * (psi(i-1,j-1,k)-psi(i-1,j-2,k))) + mpsi(i,j-1) * (d(i+1,j-1)  &
            * (psi(i+1,j-1,k)-psi(i,j-1,k)) - d(i,j-1) * (psi(i,j-1,k)-psi(i-1,j-1,k)) + c(i,j)  &
            * (psi(i,j,k)-psi(i,j-1,k)) - c(i,j-1) * (psi(i,j-1,k)-psi(i,j-2,k))) + mpsi(i-1,j)  &
            * (d(i,j) * (psi(i,j,k)-psi(i-1,j,k)) - d(i-1,j) * (psi(i-1,j,k)-psi(i-2,j,k))  &
            + c(i-1,j+1)*(psi(i-1,j+1,k)-psi(i-1,j,k)) - c(i-1,j)*(psi(i-1,j,k)-psi(i-1,j-1,k)))  &
            + mpsi(i,j) * (d(i+1,j)*(psi(i+1,j,k)-psi(i,j,k)) - d(i,j)*(psi(i,j,k)-psi(i-1,j,k))  &
            + c(i,j+1) * (psi(i,j+1,k)-psi(i,j,k)) - c(i,j) * (psi(i,j,k)-psi(i,j-1,k)))  &
            - Quarter / (phi(i,j,k)-phi(i,j,k-1))  &
            * ((mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k))  &
            + mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)))  &
            * (d(i,j-1) * (psi(i,j-1,k+1) - psi(i-1,j-1,k+1) - psi(i,j-1,k-1) + psi(i-1,j-1,k-1)) &
            + d(i,j) * (psi(i,j,k+1) - psi(i-1,j,k+1) - psi(i,j,k-1) + psi(i-1,j,k-1)))  &
            + (mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k))  &
            + mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k)))  &
            * (c(i-1,j) * (psi(i-1,j,k+1) - psi(i-1,j-1,k+1) - psi(i-1,j,k-1) + psi(i-1,j-1,k-1)) &
            + c(i,j) * (psi(i,j,k+1) - psi(i,j-1,k+1) - psi(i,j,k-1) + psi(i,j-1,k-1)))))))  &
            * adjpv(i,j,k)

       ps(i-1,j-1,k-1) = ps(i-1,j-1,k-1) - b(i,j,k) * (Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))&
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j-1) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            + c(i-1,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))  &
            * adjpv(i,j,k)
       
       ps(i,j-1,k-1) = ps(i,j-1,k-1) - b(i,j,k) * (Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (c(i,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))  &
            - d(i,j-1) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))))  &
            * adjpv(i,j,k)

       ps(i-1,j,k-1) = ps(i-1,j,k-1) - b(i,j,k) * (Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            - c(i-1,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))  &
            * adjpv(i,j,k)

       ps(i,j,k-1) = ps(i,j,k-1) - b(i,j,k) * (Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (-d(i,j) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            - c(i,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))  &
            * adjpv(i,j,k)

       ps(i-1,j-2,k) = ps(i-1,j-2,k) + b(i,j,k) * (Half * rdx**2 * mpsi(i-1,j-1) * c(i-1,j-1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))&
            * adjpv(i,j,k)
       ps(i,j-2,k) = ps(i,j-2,k) + b(i,j,k) * (Half * rdx**2 * mpsi(i,j-1) * c(i,j-1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))&
            * adjpv(i,j,k)
       ps(i-2,j-1,k) = ps(i-2,j-1,k) + b(i,j,k) * (Half * rdx**2 * mpsi(i-1,j-1) * d(i-1,j-1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))&
            * adjpv(i,j,k)

       ps(i-1,j-1,k) = ps(i-1,j-1,k) + b(i,j,k) * (Half * rdx**2 * (e(i,j,k-1)  &
            * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) * (mpsi(i,j-1)&
            * d(i,j-1) + mpsi(i-1,j-1) * (-d(i,j-1) - d(i-1,j-1) - c(i-1,j) - c(i-1,j-1))  &
            + mpsi(i-1,j) * c(i-1,j))) * adjpv(i,j,k)
       ps(i,j-1,k) = ps(i,j-1,k) + b(i,j,k) * (Half * rdx**2 * (e(i,j,k-1)  &
            * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) * (mpsi(i,j-1)&
            * (-d(i+1,j-1) - d(i,j-1) - c(i,j) - c(i,j-1)) + mpsi(i-1,j-1) * d(i,j-1) + mpsi(i,j) &
            * c(i,j))) * adjpv(i,j,k)

       ps(i+1,j-1,k) = ps(i+1,j-1,k) + b(i,j,k) * (Half * rdx**2 * mpsi(i,j-1) * d(i+1,j-1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))&
            * adjpv(i,j,k)
       ps(i-2,j,k) = ps(i-2,j,k) + b(i,j,k) * (Half * rdx**2 * mpsi(i-1,j) * d(i-1,j)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))&
            * adjpv(i,j,k)

       ps(i-1,j,k) = ps(i-1,j,k) + b(i,j,k) * (Half * rdx**2 * (e(i,j,k-1)  &
            * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) * (mpsi(i,j)  &
            * d(i,j) + mpsi(i-1,j) * (-d(i,j) - d(i-1,j) - c(i-1,j+1) - c(i-1,j)) + mpsi(i-1,j-1) &
            * c(i-1,j))) * adjpv(i,j,k)
       ps(i,j,k) = ps(i,j,k) + b(i,j,k) * (Half * rdx**2 * (e(i,j,k-1)  &
            * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) * (mpsi(i,j)  &
            * (-d(i+1,j) - d(i,j) - c(i,j+1) - c(i,j)) + mpsi(i-1,j) * d(i,j) + mpsi(i,j-1)  &
            * c(i,j))) * adjpv(i,j,k)

       ps(i+1,j,k) = ps(i+1,j,k) + b(i,j,k) * (Half * rdx**2 * mpsi(i,j) * d(i+1,j) * (e(i,j,k-1) &
            * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k)))) * adjpv(i,j,k)
       ps(i-1,j+1,k) = ps(i-1,j+1,k) + b(i,j,k) * (Half * rdx**2 * mpsi(i-1,j) * c(i-1,j+1)  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))))&
            * adjpv(i,j,k)
       ps(i,j+1,k) = ps(i,j+1,k) + b(i,j,k) * (Half * rdx**2 * mpsi(i,j) * c(i,j+1) * (e(i,j,k-1) &
            * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k)))) * adjpv(i,j,k)

       ps(i-1,j-1,k+1) = ps(i-1,j-1,k+1) - b(i,j,k) * (Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))&
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (-d(i,j-1) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)) &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            - c(i-1,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))  &
            * adjpv(i,j,k)
       
       ps(i,j-1,k+1) = ps(i,j-1,k+1) - b(i,j,k) * (Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j-1) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k)) &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            - c(i,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))  &
            * adjpv(i,j,k)

       ps(i-1,j,k+1) = ps(i-1,j,k+1) - b(i,j,k) * (Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (c(i-1,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))  &
            - d(i,j) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))))  &
            * adjpv(i,j,k)
    
       ps(i,j,k+1) = ps(i,j,k+1) - b(i,j,k) * (Eighth * rdx**2 / (phi(i,j,k)-phi(i,j,k-1))  &
            * (e(i,j,k-1) * (phi(i,j,k-1)-phi(i,j,k-2)) - e(i,j,k+1) * (phi(i,j,k+1)-phi(i,j,k))) &
            * (d(i,j) * (mu(i,j) * (phi(i+1,j,k-1) - phi(i,j,k-1) + phi(i+1,j,k) - phi(i,j,k))  &
            + mu(i-1,j) * (phi(i,j,k-1) - phi(i-1,j,k-1) + phi(i,j,k) - phi(i-1,j,k)))  &
            + c(i,j) * (mv(i,j) * (phi(i,j+1,k-1) - phi(i,j,k-1) + phi(i,j+1,k) - phi(i,j,k))  &
            + mv(i,j-1) * (phi(i,j,k-1) - phi(i,j-1,k-1) + phi(i,j,k) - phi(i,j-1,k)))))  &
            * adjpv(i,j,k)       
    end forall
    ph((/1,nx/),:,:) = Zero
    ph(2:nx-1,(/1,ny/),:) = Zero
    ph(2:nx-1,2:ny-1,(/0,nz/)) = Zero
    ps((/0,nx/),:,:) = Zero
    ps(1:nx-1,(/0,ny/),:) = Zero
    deallocate(d, c, mpsi, e, b)
  end subroutine pv_adj

  ! ===========================================================================

  subroutine solve_phi_psi(pvGiven, psi, phi, f, mCap, mu, mfm, mfu, mfv, etaz, etam, pt, rdx)
    real(cp), dimension(:,:,:),   intent(in) :: pvGiven
    real(cp), dimension(0:,0:,:), intent(inout) :: psi          ! (0:nx,0:ny)
    real(cp), dimension(:,:,0:),  intent(inout) :: phi
    real(cp), dimension(:,:),     intent(in) :: f, mCap, mu, mfm   ! (1:nx,1:ny)
    real(cp), dimension(0:,:),    intent(in) :: mfu
    real(cp), dimension(:,0:),    intent(in) :: mfv           ! (1:nx,0:ny)
    real(cp), dimension(0:),      intent(in) :: etaz
    real(cp), dimension(:),       intent(in) :: etam
    real(cp),                     intent(in) :: pt, rdx

    real(cp),         parameter :: e = 150000., MaxStep = 2187.
    character(len=*), parameter :: Mismatch = "mismatch in solve_phi_psi!"

    real(cp), dimension(:,:,:), allocatable :: imb, pv
    real(cp), dimension(:,:,:), allocatable :: gradphi1, gradphi2, gradphi, pphi, phiTrial
    real(cp), dimension(:,:,:), allocatable :: gradpsi1, gradpsi2, gradpsi, ppsi, psiTrial
    real, dimension(:,:,:), allocatable :: pv4
    real(cp) :: costa, costb, cost, costlb, costrb, costleft, costright, oldcost
    real(cp) :: alpha, alphalb, alpharb, alphaleft, alpharight
    integer :: nx, ny, nz, i, j, k
    
    nx = size(psi,1) - 1
    ny = size(psi,2) - 1
    nz = size(psi,3)

    ! Sanity checks
    if (any(shape(pvGiven) /= (/ nx, ny, nz /) )) stop "pvGiven " // Mismatch     
    if (any(shape(psi) /= (/ nx+1,  ny+1, nz /) )) stop "psi " // Mismatch    
    if (any(shape(phi) /= (/ nx,  ny, nz+1 /) )) stop "phi " // Mismatch   
    if (any(shape(f)  /= (/ nx,   ny   /) )) stop "f "  // Mismatch
    if (any(shape(mCap) /= (/ nx, ny   /) )) stop "mCap " // Mismatch    
    if (any(shape(mu) /= (/ nx, ny /) )) stop "mu "  // Mismatch    
    if (any(shape(mfm)  /= (/ nx,   ny   /) )) stop "mfm "  // Mismatch
    if (any(shape(mfu) /= (/ nx+1, ny   /) )) stop "mfu " // Mismatch
    if (any(shape(mfv) /= (/ nx,   ny+1 /) )) stop "mfv " // Mismatch
    if (size(etaz) /= nz+1) stop "etaz " // Mismatch
    if (size(etam) /= nz) stop "etam " // Mismatch

    allocate(imb(nx,ny,nz), pv(nx,ny,nz))
    allocate(gradphi1(nx,ny,0:nz), gradphi2(nx,ny,0:nz), gradphi(nx,ny,0:nz), pphi(nx,ny,0:nz),  &
         phiTrial(nx,ny,0:nz))
    allocate(gradpsi1(0:nx,0:ny,nz), gradpsi2(0:nx,0:ny,nz), gradpsi(0:nx,0:ny,nz),  &
         ppsi(0:nx,0:ny,nz), psiTrial(0:nx,0:ny,nz))
    allocate(pv4(nx,ny,nz))

    i = 1
    do
       print *, "Iteration: ", i
       ! Calculate cost function
       call imbalance_exp(psi, phi, f, mCap, mfm, mfu, mfv, etaz, etam, rdx, imb)
       costb = Half * e * sum(imb**2)

       call pv_exp(psi, phi, f, mCap, mu, mfu, mfv, etaz, etam, pt, rdx, pv4)
       pv = pv4
       costa = Half * sum( (pv - pvGiven)**2 )

       cost = costa + costb
       write(*,*) " Cost A & B:", costa, costb
       write(*,*) " Total cost:", cost

       write(55,*), costa, costb, cost

       ! Calculate gradient of cost function
       call pv_adj(pv-pvGiven, psi, phi, f, mCap, mu, mfu, mfv, etaz, etam, pt, rdx, gradpsi1,  &
            gradphi1)
       call imbalance_simp_adj(imb, psi, f, mCap, mfm, mfu, mfv, etaz, etam, rdx, gradpsi2,  &
            gradphi2)
       gradphi = gradphi1 + e * gradphi2
       gradpsi = gradpsi1 + e * gradpsi2

!       pphi = gradphi / sqrt(sum(gradphi**2)+sum(gradpsi**2))
!       ppsi = gradpsi / sqrt(sum(gradphi**2)+sum(gradpsi**2))
       pphi = gradphi / sqrt(sum(gradphi**2))
       ppsi = gradpsi / sqrt(sum(gradpsi**2))

       ! Perform line search
       alphalb = Zero
       costlb = cost

       alpharb = Maxstep
       
       phiTrial = phi - alpha * pphi
       psiTrial = psi - alpha * ppsi
       ! Calculate cost function
       call imbalance_exp(psiTrial, phiTrial, f, mCap, mfm, mfu, mfv, etaz, etam, rdx, imb)
       costb = Half * e * sum(imb**2)
       call pv_exp(psiTrial, phiTrial, f, mCap, mu, mfu, mfv, etaz, etam, pt, rdx, pv4)
       pv = pv4
       costa = Half * sum( (pv - pvGiven)**2 )
       costrb = costa + costb

       oldcost = Half * (costlb + costrb)
       do j = 1, 100
          alphaleft = alphalb + Third * (alpharb - alphalb)
          alpharight = alphalb + Two * Third * (alpharb - alphalb)
!          print *, alphalb, alphaleft, alpharight, alpharb
          
          ! Try left
          phiTrial = phi - alphaleft * pphi
          psiTrial = psi - alphaleft * ppsi
          call imbalance_exp(psiTrial, phiTrial, f, mCap, mfm, mfu, mfv, etaz, etam, rdx, imb)
          costb = Half * e * sum(imb**2)
          call pv_exp(psiTrial, phiTrial, f, mCap, mu, mfu, mfv, etaz, etam, pt, rdx, pv4)
          pv = pv4
          costa = Half * sum( (pv - pvGiven)**2 )
          costleft = costa + costb

          ! Try right
          phiTrial = phi - alpharight * pphi
          psiTrial = psi - alpharight * ppsi
          call imbalance_exp(psiTrial, phiTrial, f, mCap, mfm, mfu, mfv, etaz, etam, rdx, imb)
          costb = Half * e * sum(imb**2)
          call pv_exp(psiTrial, phiTrial, f, mCap, mu, mfu, mfv, etaz, etam, pt, rdx, pv4)
          pv = pv4
          costa = Half * sum( (pv - pvGiven)**2 )
          costright = costa + costb

          ! Which is best?
!          print *, "left:", alphaleft, costleft
!          print *, "right:", alpharight, costright
          if (costleft < costright) then
             alpharb = alpharight
!             print *, "Bounds:", alphalb, alpharb
!             print *, j, Half*(alphalb+alpharb), costleft
          else
             alphalb = alphaleft
!             print *, "Bounds:", alphalb, alpharb
!             print *, j, Half*(alphalb+alpharb), costright
          end if
          alpha = Half * (alphalb + alpharb)
          if (abs(alpharb-alphalb)/(Half*(alphalb+alpharb)) < 1e-5) then
!          if (abs(costright - costleft)/(Half*(costright+costleft)) < 1e-6) then
             print *, j, alpha
             exit
          end if
       end do

       if (Half * (costleft + costright) > oldcost) exit
       oldcost = Half * (costleft + costright)

       ! Update phi and psi
       phi = phi - alpharight * pphi
       psi = psi - alpharight * ppsi
       i = i + 1
       if (i > 10000) exit
    end do
    open(34, file="psphiinvert.dat", status="new", form="unformatted", action="write")
    write (34) nx, ny, nz
    write (34) psi
    write (34) phi
    close(34)
  end subroutine solve_phi_psi

  ! ===========================================================================

  function xdiff(data, rdx) result(res)
    real(cp), dimension(:,:), intent(in) :: data
    real(cp),                 intent(in) :: rdx
    real(cp), dimension(size(data,1)-1,size(data,2)) :: res
    
    integer :: nx
    
    nx = size(data,1)-1
    
    res = rdx * (data(2:,:) - data(:nx,:))
  end function xdiff

  ! ===========================================================================

  function ydiff(data, rdx) result(res)
    real(cp), dimension(:,:), intent(in) :: data
    real(cp),                 intent(in) :: rdx
    real(cp), dimension(size(data,1),size(data,2)-1) :: res
    
    integer :: ny
    
    ny = size(data,2)-1
    
    res = rdx * (data(:,2:) - data(:,:ny))
  end function ydiff

  ! ===========================================================================

  function xyavg(data) result(res)
    real(cp), dimension(:,:), intent(in) :: data
    real(cp), dimension(size(data,1)-1,size(data,2)-1) :: res
 
    integer :: nx, ny

    nx = size(data,1)-1
    ny = size(data,2)-1

    res = Quarter * (data(:nx,:ny) + data(2:,:ny) + data(:nx,2:) + data(2:,2:))
  end function xyavg

  ! ===========================================================================

  function xavg(data) result(res)
    real(cp), dimension(:,:), intent(in) :: data
    real(cp), dimension(size(data,1)-1,size(data,2)) :: res
 
    integer :: nx

    nx = size(data,1)-1

    res = Half * (data(:nx,:) + data(2:,:))
  end function xavg

  ! ===========================================================================

  function yavg(data) result(res)
    real(cp), dimension(:,:), intent(in) :: data
    real(cp), dimension(size(data,1),size(data,2)-1) :: res
 
    integer :: ny

    ny = size(data,2)-1
    
    res = Half * (data(:,:ny) + data(:,2:))
  end function yavg
end module partition_wind
