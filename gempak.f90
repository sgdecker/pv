module gempak
  use wrf2gem_parameters, only: NLUnit, NLFName
  implicit none
  save
  private
  public :: init_gem, create_gemfile, write_gempak, close_and_exit_gem

  integer, parameter      :: navsz = 256, ianlsz = 128, ihdrsz = 2

  real, dimension(navsz)  :: rnvblk, rnvblk2
  real, dimension(ianlsz) :: anlblk = 0
  integer, dimension(2)           :: level, ighdr = 0
  character(len=20), dimension(2) :: gdattm
  logical                         :: pack = .true., rewrit = .true.


contains  ! ===================================================================


  ! init_gem initializes the GEMPAK interface.
  subroutine init_gem()
    integer :: status

    call in_bdta(status)
    call handle_gerr(status, "ID_BDTA")
  end subroutine init_gem
  
  ! ===========================================================================

  ! create_gemfile creates a new GEMPAK file, using the following information:
  ! o Projection type (proj)
  ! o Grid extent in x and y dimensions (nx, ny)
  ! o Latitude/longitude in lower-left and upper-right grid corners (lat, lon)
  ! o Additional angles necessary for the map projection (ang)
  ! If successful, a file handle (gid) is assigned and provided as output.
  ! If the GEMPAK file already exists, the program will either use it or abort
  ! depending on the user's preference.  However, if gsFlag = .false., the grid
  ! navigation of the current GEMPAK file does not match what's needed;
  ! otherwise there's no problem.
  subroutine create_gemfile(proj, nx, ny, lat, lon, ang, gid, gsFlag)
    integer,            intent(in)  :: proj, nx, ny
    real, dimension(2), intent(in)  :: lat, lon
    real, dimension(3), intent(in)  :: ang
    integer,            intent(out) :: gid
    logical,            intent(out) :: gsFlag

    integer                 :: maxGrd = 5000, status
    character(len=80)       :: gFName = "wrf.gem"
    character(len=3)        :: projName
    logical                 :: overwrite = .false., nlExist, gemExist, angFlg
    
    namelist /GEMFILE/ gFName, maxGrd, pack, overwrite

    ! Check for namelist existence, then read it.
    inquire (file=NLFName, exist=nlExist)    
    if (nlExist) then
       open (unit=NLUnit, file=NLFName, iostat=status)
       if (status /= 0) stop "There is a problem opening the namelist file."
       read (unit=NLUnit, nml=GEMFILE, iostat=status)
       if (status /= 0) stop "There is a problem reading the namelist file."
       close (unit=NLUnit)
    end if

    ! If the GEMPAK file already exists, follow user's wishes.
    inquire (file=gFName, exist=gemExist)
    if (gemExist) then
       if (overwrite) then
          print *, "GEMPAK file already exists! wrf2gem will attempt to add to&
               & and/or overwrite it."
       else
          stop "GEMPAK file already exists!"
       end if
    end if

    ! Set up map projection
    select case(proj)
       case (0)
          projName = "CED"
          angFlg = .false.
       case (1)
          if (ang(1) < 0 .and. ang(3) < 0) then
             projName = "SCC"
          else
             projName = "LCC"
          end if
          angFlg = .true.
       case (3)
          projName = "MER"
          angFlg = .false.
       case default
          print *, "This projection not yet implemented."
          stop
    end select
       
    ! Call appropriate gemlib routines to create or add to the GEMPAK file.
    call gr_mnav(projName, nx, ny, lat(1), lon(1), lat(2), lon(2), ang(1),  &
         ang(2), ang(3), angFlg, rnvblk, status)
    call handle_gerr(status, "GR_MNAV")
    if (gemExist) then
       call gd_opnr(gFName, gid, navsz, rnvblk2, ianlsz, anlblk, ihdrsz,  &
            maxGrd, status)
       call handle_gerr(status, "GD_OPNR")
       ! Make sure grid navigation blocks match.
       call gr_cnav(rnvblk, rnvblk2, navsz, gsFlag, status)
       call handle_gerr(status, "GR_CNAV")
    else
       call gd_cref(gFName, navsz, rnvblk, ianlsz, anlblk, ihdrsz, maxGrd,  &
            gid, status)
       call handle_gerr(status, "GD_CREF")
       gsFlag = .true.
    end if
  end subroutine create_gemfile

  ! ===========================================================================

  ! write_gempak writes a 2D array of data (grid) to the GEMPAK file associated
  ! with id. The extents of grid are specified by nx and ny. The time at which
  ! the array is valid is specified in GEMPAK form in gdat1. It is assumed that
  ! the grid does not have two times associated with it. GEMPAK level
  ! information is provided by lev1 and lev2, and the GEMPAK parameter name is
  ! given by parm.
  subroutine write_gempak(id, grid, nx, ny, gdat1, lev1, lev2, cord, parm)
    integer,                intent(in) :: id, nx, ny, lev1, lev2, cord
    real, dimension(nx,ny), intent(in) :: grid
    character(len=20),      intent(in) :: gdat1
    character(len=12),      intent(in) :: parm

    integer :: status

    ! Initialize GEMPAK level and date/time information
    level(1) = lev1 ; level(2) = lev2
    gdattm(1) = gdat1 ; gdattm(2) = " "

    ! Call appropriate gemlib routines
    if (pack) then
       call gd_wpgd(id, grid, nx, ny, ighdr, gdattm, level, cord, parm,  &
            rewrit, 1, 16, status)
       call handle_gerr(status, "GD_WPGD")
    else
       call gd_wdat(id, grid, nx, ny, ighdr, gdattm, level, cord, parm,  &
            rewrit, status)
       call handle_gerr(status, "GD_WDAT")
    end if
  end subroutine write_gempak

  ! ===========================================================================

  ! close_and_exit_gem closes the GEMPAK file associated with igdfln as well as
  ! the gemlib interface.
  subroutine close_and_exit_gem(igdfln)
    integer, intent(in) :: igdfln

    integer :: status
    
    call gd_clos(igdfln, status)
    call handle_gerr(status, "GD_CLOS")
    call ip_exit(status)
    call handle_gerr(status, "IP_EXIT")
  end subroutine close_and_exit_gem

  ! ===========================================================================

  ! handle_gerr examines the result code from a gemlib routine. If necessary,
  ! it then prints an error message and aborts the program.
  subroutine handle_gerr(status, routine)
    integer,          intent(in) :: status
    character(len=*), intent(in) :: routine
    
    if (status /= 0) then
       write (*,"(a,a,i3)") routine, ": status = ", status
       stop
    end if
  end subroutine handle_gerr
end module gempak
