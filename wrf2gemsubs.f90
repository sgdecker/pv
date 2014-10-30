module wrf2gem_subs
  use netcdf, only: nf90_open, nf90_inq_dimid, nf90_inquire_dimension,  &
                    nf90_get_att, nf90_inq_varid, nf90_get_var, nf90_close,  &
                    nf90_inquire_variable
  use netcdf, only: NF90_NoWrite, NF90_NoErr, NF90_Global, NF90_StrError,  &
                    NF90_Max_Var_Dims
  use dateutil, only: hr_diff
  use dateutil, only: operator(<), operator(>)
  use dateutil, only: datestamp
  use wrf2gem_parameters, only: NLUnit, NLFName
  implicit none
  save
  private

  public :: get_out_fields, get_num_files, open_wrf_file, get_grid_info,  &
            update_grid_info, get_times, close_wrf_file, handle_err,  &
            mem_error, read_values, read_values_const, freq_divisible
  public :: NLUnit, NLFName

  integer, parameter :: MaxOuts = 200

  integer,           dimension(MaxOuts), target :: outFields
  character(len=80), dimension(MaxOuts), target :: fName

  interface read_values
     module procedure rv0d
     module procedure rv1d
     module procedure rv2d
     module procedure rv3d
     module procedure rv4d
  end interface

  interface read_values_const
     module procedure rv0dc
     module procedure rv1dc
     module procedure rv2dc
     module procedure rv3dc
  end interface

contains  ! ===================================================================

  
  ! get_out_fields reads the namelist to determine what variables should be
  ! output to the GEMPAK file.  It returns the number of variables to output
  ! (num) and a pointer to an array of variable IDs (ptr).  If there are no
  ! variables to output, the status of ptr does not change.  The subroutine
  ! also determines how pressure interpolation should be handled, giving
  ! values for the bottom and top pressure levels (pb, pt) and the pressure
  ! interval (dp).
  subroutine get_out_fields(num, ptr, pb, pt, dp)
    integer, intent(out) :: num
    integer, dimension(:), pointer :: ptr
    integer, intent(out) :: pb, pt, dp

    integer :: status
    logical :: nlExist

    namelist /NUMOUTS/ num
    namelist /VARS/ outFields, pb, pt, dp

    ! The default if there is no namelist information is to output nothing.
    num = 0

    ! Defaults for pressure level interpolation information
    pb = 1000; pt = 200; dp = 50

    ! Check for namelist existence, then read it.
    inquire (file=NLFName, exist=nlExist)
    if (nlExist) then
       open (unit=NLUnit, file=NLFName, iostat=status)
       if (status /= 0) stop "There is a problem opening the namelist file."
       read (unit=NLUnit, nml=NUMOUTS, iostat=status)
       if (status /= 0) stop "There is a problem reading NUMOUTS from the &
            &namelist file."
       
       if (num > MaxOuts) stop "There are too many outputs!"

       read (unit=NLUnit, nml=VARS, iostat=status)
       if (status /= 0) stop "There is a problem reading VARS from the &
            &namelist file."
       close (unit=NLUnit)

       if (num > 0) ptr => outFields(1:num)
    end if
  end subroutine get_out_fields

  ! ===========================================================================
  
  ! get_num_files reads the namelist to determine what WRF files should be
  ! read for conversion.  It returns the number of files to convert (numFiles)
  ! and a pointer to an array of filenames (point).  If there are no files to
  ! convert, the status of point does not change.
  subroutine get_num_files(numFiles, point)
    integer, intent(out) :: numFiles
    character(len=80), dimension(:), pointer :: point
    
    integer :: status
    logical :: nlExist

    namelist /OUTFILES/ numFiles
    namelist /WRFFILE/ fName

    ! The default if there is no namelist information is to output nothing.
    numFiles = 0

    ! Check for namelist existence, then read it.
    inquire (file=NLFName, exist=nlExist)
    if (nlExist) then
       open (unit=NLUnit, file=NLFName, iostat=status)
       if (status /= 0) stop "There is a problem opening the namelist file."
       read (unit=NLUnit, nml=OUTFILES, iostat=status)
       if (status /= 0) stop "There is a problem reading OUTFILES from the &
            &namelist file."

       if (numFiles > MaxOuts) stop "There are too many files to process!"

       read (unit=NLUnit, nml=WRFFILE, iostat=status)
       if (status /= 0) stop "There is a problem reading WRFFILE from the &
            &namelist file."
       close (unit=NLUnit)

       if (numFiles > 0) point => fName(1:numFiles)
    end if
  end subroutine get_num_files

  ! ===========================================================================

  ! open_wrf_file opens the WRF history file (fil) and returns its handle (id).
  subroutine open_wrf_file(fil, id)
    character(len=80), intent(in)  :: fil
    integer,           intent(out) :: id
    
    integer :: status
    
    status = nf90_open(fil, NF90_NoWrite, id)
    call handle_err(status)    
  end subroutine open_wrf_file

  ! ===========================================================================

  ! get_grid_info gathers the following info from the WRF history file (ncid):
  ! o Model grid extent in x, y, and z directions (nx, ny, nz)
  ! o Number of times in the file (nt)
  ! o Map projection type (proj)
  ! o Latitude/longitude in lower-left and upper-right grid corners (lat, lon)
  ! o Projection angle info (ang)
  subroutine get_grid_info(ncid, nx, ny, nz, nt, proj, lat, lon, ang)
    integer,            intent(in)  :: ncid
    integer,            intent(out) :: nx, ny, nz, nt, proj
    real, dimension(2), intent(out) :: lat, lon
    real, dimension(3), intent(out) :: ang
    
    character(len=11), dimension(4), parameter :: DimName = (/ "west_east  ", &
                                 "south_north", "bottom_top ", "Time       " /)
    integer, dimension(4) :: temp
    integer               :: i, varid, status

    do i = 1, 4
       status = nf90_inq_dimid(ncid, trim(DimName(i)), varid)
       call handle_err(status)  
       status = nf90_inquire_dimension(ncid, varid, len=temp(i))
       call handle_err(status)
    end do
    nx = temp(1)
    ny = temp(2)
    nz = temp(3)
    nt = temp(4)
        
    status = nf90_get_att(ncid, NF90_Global, "MAP_PROJ", proj)
    call handle_err(status) 
    
    status = nf90_inq_varid(ncid, "XLONG", varid)
    call handle_err(status)
    status = nf90_get_var(ncid, varid, lon(1), start = (/ 1, 1, 1 /))
    call handle_err(status)
    status = nf90_get_var(ncid, varid, lon(2), start = (/ nx, ny, 1 /))
    call handle_err(status)
    
    status = nf90_inq_varid(ncid, "XLAT", varid)
    call handle_err(status)
    status = nf90_get_var(ncid, varid, lat(1), start = (/ 1, 1, 1 /))
    call handle_err(status)
    status = nf90_get_var(ncid, varid, lat(2), start = (/ nx, ny, 1 /))
    call handle_err(status)
    
    select case (proj)
    case (0) ! Ideal x/y
       ang = 0
       lon = (/ -100. - nx/2., -100. + nx/2. /)
       lat = (/ 45. - ny/2., 45. + ny/2. /)
    case (1) ! LCC
       status = nf90_get_att(ncid, NF90_Global, "TRUELAT1", ang(1))
       call handle_err(status) 
       status = nf90_get_att(ncid, NF90_Global, "CEN_LON", ang(2))
       call handle_err(status) 
       status = nf90_get_att(ncid, NF90_Global, "TRUELAT2", ang(3))
       call handle_err(status)
    case (3) ! Mercator
       ang = 0
    case default
       write (*,"(a,i2)") "proj = ", proj
       print *, "This map projection has not yet been implemented."
       stop
    end select
  end subroutine get_grid_info

  ! ===========================================================================
  
  ! update_grid_info compares grid size and prejection settings from the WRF
  ! history file (ncid) with the following previously saved settings:
  ! o Model grid extent in x, y, and z directions (nx, ny, nz)
  ! o Map projection type (proj)
  ! o Latitude/longitude in lower-left and upper-right grid corners (lat, lon)
  ! o Projection angle info (ang)
  ! If there is not an exact match, ok is set to .false.; otherwise, ok is
  ! .true. The number of times in the WRF file (nt) is also provided as output.
  subroutine update_grid_info(ncid, nx, ny, nz, proj, lat, lon, ang, nt, ok)
    integer,            intent(in)  :: ncid, nx, ny, nz, proj
    real, dimension(2), intent(in)  :: lat, lon
    real, dimension(3), intent(in)  :: ang
    integer,            intent(out) :: nt
    logical,            intent(out) :: ok
    
    character(len=11), dimension(4), parameter :: DimName = (/ "west_east  ", &
                                 "south_north", "bottom_top ", "Time       " /)
    real,    dimension(3) :: angl
    real,    dimension(2) :: lonl, latl
    integer, dimension(4) :: temp
    integer               :: i, varid, status, nxl, nyl, nzl, projl

    do i = 1, 4
       status = nf90_inq_dimid(ncid, trim(DimName(i)), varid)
       call handle_err(status)  
       status = nf90_inquire_dimension(ncid, varid, len=temp(i))
       call handle_err(status)
    end do
    nxl = temp(1)
    nyl = temp(2)
    nzl = temp(3)
    nt = temp(4)
        
    status = nf90_get_att(ncid, NF90_Global, "MAP_PROJ", projl)
    call handle_err(status) 
    
    status = nf90_inq_varid(ncid, "XLONG", varid)
    call handle_err(status)
    status = nf90_get_var(ncid, varid, lonl(1), start = (/ 1, 1, 1 /))
    call handle_err(status)
    status = nf90_get_var(ncid, varid, lonl(2), start = (/ nxl, nyl, 1 /))
    call handle_err(status)
    
    status = nf90_inq_varid(ncid, "XLAT", varid)
    call handle_err(status)
    status = nf90_get_var(ncid, varid, latl(1), start = (/ 1, 1, 1 /))
    call handle_err(status)
    status = nf90_get_var(ncid, varid, latl(2), start = (/ nxl, nyl, 1 /))
    call handle_err(status)
    
    select case (projl)
    case (0)  ! Ideal x/y
       angl = 0
       lonl = (/ -100. - nx/2., -100. + nx/2. /)
       latl = (/ 45. - ny/2., 45. + ny/2. /)
    case (1) ! LCC
       status = nf90_get_att(ncid, NF90_Global, "TRUELAT1", angl(1))
       call handle_err(status) 
       status = nf90_get_att(ncid, NF90_Global, "CEN_LON", angl(2))
       call handle_err(status) 
       status = nf90_get_att(ncid, NF90_Global, "TRUELAT2", angl(3))
       call handle_err(status)
    case (3) ! Mercator
       angl = 0
    case default
       write (*,"(a,i2)") "proj = ", proj
       print *, "This map projection has not yet been implemented."
       stop
    end select

    ! Check for exact match
    ok = (nxl == nx .and. nyl == ny .and. nzl == nz .and.  &
         projl == proj .and. all(latl == lat) .and. all(lonl == lon) .and.  &
         all(angl == ang))
  end subroutine update_grid_info

  ! ===========================================================================

  ! get_times uses information from the WRF history file (id) and the namelist
  ! to create an appropriate array of output times in GEMPAK form (times).  In
  ! addition, the number of outputs per 12 hours (timesPer12hr) is provided for
  ! output.  Also provided as input is the list of history files (outFiles),
  ! and the index of the current history file being processed (outInd).  This
  ! subroutine assumes the history file contains either output every n hours,
  ! where n is some positive integer, or output for one time.
  subroutine get_times(outInd, outFiles, id, times, timesPer12hr)
    integer,                         intent(in)    :: outInd
    character(len=80), dimension(:), intent(in)    :: outFiles
    integer,                         intent(inout) :: id
    character(len=*),  dimension(:), intent(out)   :: times
    integer,                         intent(out)   :: timesPer12hr

    integer, parameter :: y1 = 1, y2 = 4, m1 = 6, m2 = 7, d1 = 9, d2 = 10,  &
                          h1 = 12, h2 = 13

    character(len=19), dimension(size(times,1)) :: wrfTimes
    type(datestamp) :: fileInit, runInit, adj, temp
    integer :: numTimes, varid, status, adjID, interval, ihr1, i
    character(len=19) :: adjTimes, startDate
    character(len=3)  :: fhr
    character(len=2)  :: chr1, chr2

    ! Initialization
    numTimes = size(times,1)
    
    ! Read WRF history file times
    ! The format for these times is yyyy-mm-dd_hh:mm:ss
    status = nf90_inq_varid(id, "Times", varid)
    call handle_err(status)
    status = nf90_get_var(id, varid, wrfTimes)
    call handle_err(status)

    ! Determine initial time from the history file
    status = nf90_get_att(id, nf90_global, "START_DATE", startDate)
    call handle_err(status)    
    
    fileInit%yr  = char_to_int(wrfTimes(1)(y1:y2))
    runInit%yr   = char_to_int(  startDate(y1:y2))
    fileInit%mth = char_to_int(wrfTimes(1)(m1:m2))
    runInit%mth  = char_to_int(  startDate(m1:m2))
    fileInit%day = char_to_int(wrfTimes(1)(d1:d2))
    runInit%day  = char_to_int(  startDate(d1:d2))
    fileInit%hr  = char_to_int(wrfTimes(1)(h1:h2))
    runInit%hr   = char_to_int(  startDate(h1:h2))

    ! Calculate interval between outputs in hours, and number of outputs per 12
    ! hours.
    if (numTimes > 1) then
       adj%yr  = char_to_int(wrfTimes(2)(y1:y2))
       adj%mth = char_to_int(wrfTimes(2)(m1:m2))
       adj%day = char_to_int(wrfTimes(2)(d1:d2))
       adj%hr  = char_to_int(wrfTimes(2)(h1:h2))
    else if (size(outFiles) > 1) then
       call close_wrf_file(id)  ! My system can't handle two open netCDF files!
       if (outInd > 1) then
          call open_wrf_file(outFiles(outInd-1), adjID)
       else
          call open_wrf_file(outFiles(2), adjID)
       end if
       status = nf90_inq_varid(adjID, "Times", varid)
       call handle_err(status)
       status = nf90_get_var(id, varid, adjTimes)
       call handle_err(status)
       call close_wrf_file(adjID)
       call open_wrf_file(outFiles(outInd), id)

       adj%yr  = char_to_int(adjTimes(y1:y2))
       adj%mth = char_to_int(adjTimes(m1:m2))
       adj%day = char_to_int(adjTimes(d1:d2))
       adj%hr  = char_to_int(adjTimes(h1:h2))
    else
       adj = fileInit
    end if

    if (adj > fileInit) then
       interval = hr_diff(adj, fileInit)
    else if (adj < fileInit) then
       interval = hr_diff(fileInit, adj)
    else
       interval = 12  ! Set so that timesPer12hr = 1
    end if
    timesPer12hr = 12 / interval   

    ! Calculate forecast hour to which the first time in history file
    ! corresponds
    ihr1 = hr_diff(fileInit, runInit)
    if (ihr1 < 0) then
       print *, "WARNING: Can't handle forecast hour calculation properly!"
       ihr1 = 0
    end if

    ! Form GEMPAK times for all of the times in the history file
    ! GEMPAK format is yymmdd/hhmmFfff
    do i = 1, numTimes
       times(i)(1:4) = int_to_char2(runInit%yr) // int_to_char2(runInit%mth)
       times(i)(5:7) = int_to_char2(runInit%day) // "/"
       times(i)(8:12) = int_to_char2(runInit%hr) // wrfTimes(1)(15:16) // "F"

       write (unit=fhr, fmt="(i3)") ihr1
       select case (ihr1)
       case (0:9)
          fhr(1:2) = "00"
       case (10:99)
          fhr(1:1) = "0"
       end select

       times(i)(13:15) = fhr
       ihr1 = ihr1 + interval
    end do
  end subroutine get_times

  ! ===========================================================================

  ! close_wrf_file closes the WRF history file specified by id.
  subroutine close_wrf_file(id)
    integer, intent(in) :: id
    
    integer :: status
    
    status = nf90_close(id)
    call handle_err(status)
  end subroutine close_wrf_file
    
  ! ===========================================================================

  ! handle_err examines the result code from a netCDF function (status). If
  ! necessary, it then prints a netCDF error message and aborts the program.
  subroutine handle_err(status)
    integer, intent(in) :: status
    
    if (status /= NF90_NoErr) then
       print *, trim(NF90_StrError(status))
       print *, "Stopped"
       stop
    end if
  end subroutine handle_err

  ! ===========================================================================
  
  ! mem_error examines the result code from a memory allocation/deallocation
  ! (status).  id determines whether allocation or deallocation was attempted.
  ! string specifies the program unit that attempted the memory operation.  If
  ! necessary, an error message is printed and the program is aborted.
  subroutine mem_error(status, id, string)
    integer, intent(in) :: status, id
    character(len=*), intent(in) :: string

    if (status /= 0) then
       select case (id)
       case (1)
          print *, "Memory allocation error in ", string
       case (2)
          print *, "Memory deallocation error in ", string
       case default
          print *, "Memory error in ", string
       end select
       
       stop
    end if
  end subroutine mem_error

  ! ===========================================================================

  ! rv0d supplies a scalar (value) containing the variable specified by varName
  ! from the netCDF file associated with id.
  subroutine rv0d(id, varName, value)
    integer,          intent(in)  :: id
    character(len=*), intent(in)  :: varName
    real,             intent(out) :: value
    
    integer, dimension(NF90_Max_Var_Dims) :: dimIDs
    integer :: varid, status, nDims, ncSize
    
    ! Get array extents from the netCDF file
    status = nf90_inq_varid(id, varName, varid)
    call handle_err(status)
    status = nf90_inquire_variable(id, varid, ndims = nDims, dimids = dimIDs)
    call handle_err(status)
    status = nf90_inquire_dimension(id, dimIDs(1), len = ncSize)
    call handle_err(status)
    
    ! Make sure array extents match!
    if (ncSize /= 1 .or. nDims /= 0) then
       print *, "Array dimensions are incorrect in read_values."
       stop
    end if
    
    ! Read the variable from the netCDF file
    status = nf90_get_var(id, varid, value)
    call handle_err(status)
  end subroutine rv0d

  ! ...........................................................................

  ! rv1d supplies a 1D array (values) containing the variable specified by
  ! varName from the netCDF file associated with id.
  subroutine rv1d(id, varName, values)
    integer,            intent(in)  :: id
    character(len=*),   intent(in)  :: varName
    real, dimension(:), intent(out) :: values
    
    integer, parameter :: NumDims = 1
    integer, dimension(NF90_Max_Var_Dims) :: dimIDs
    integer, dimension(NumDims)           :: sizes, ncSizes
    integer :: i, varid, status, nDims
    
    ! Determine array extents that wrf2gem expects
    sizes = shape(values)
    
    ! Get array extents from the netCDF file
    status = nf90_inq_varid(id, varName, varid)
    call handle_err(status)
    status = nf90_inquire_variable(id, varid, ndims = nDims, dimids = dimIDs)
    call handle_err(status)
    do i = 1, NumDims
       status = nf90_inquire_dimension(id, dimIDs(i), len = ncSizes(i))
       call handle_err(status)
    end do
    
    ! Make sure array extents match!
    if (any(ncSizes /= sizes) .or. nDims /= NumDims) stop  &
         "Array dimensions are incorrect in read_values."
    
    ! Read the variable from the netCDF file
    status = nf90_get_var(id, varid, values)
    call handle_err(status)
  end subroutine rv1d

  ! ...........................................................................

  ! rv2d is identical to rv1d, but for a 2D array.
  subroutine rv2d(id, varName, values)
    integer,              intent(in)  :: id
    character(len=*),     intent(in)  :: varName
    real, dimension(:,:), intent(out) :: values
    
    integer, parameter :: NumDims = 2
    integer, dimension(NF90_Max_Var_Dims) :: dimIDs
    integer, dimension(NumDims)           :: sizes, ncSizes
    integer :: i, varid, status, nDims
    
    ! Determine array extents that wrf2gem expects
    sizes = shape(values)
    
    ! Get array extents from the netCDF file
    status = nf90_inq_varid(id, varName, varid)
    call handle_err(status)
    status = nf90_inquire_variable(id, varid, ndims = nDims, dimids = dimIDs)
    call handle_err(status)
    do i = 1, NumDims
       status = nf90_inquire_dimension(id, dimIDs(i), len = ncSizes(i))
       call handle_err(status)
    end do
    
    ! Make sure array extents match!
    if (any(ncSizes /= sizes) .or. nDims /= NumDims) stop  &
         "Array dimensions are incorrect in read_values."

    ! Read the variable from the netCDF file
    status = nf90_get_var(id, varid, values)
    call handle_err(status)
  end subroutine rv2d

  ! ...........................................................................
  
  ! rv3d is identical to rv2d, but for a 3D array.
  subroutine rv3d(id, varName, values)
    integer,                intent(in)  :: id
    character(len=*),       intent(in)  :: varName
    real, dimension(:,:,:), intent(out) :: values

    integer, parameter :: NumDims = 3
    integer, dimension(NF90_Max_Var_Dims) :: dimIDs
    integer, dimension(NumDims)           :: sizes, ncSizes
    integer :: i, varid, status, nDims

    ! Determine array extents that wrf2gem expects
    sizes = shape(values)
    
    ! Get array extents from the netCDF file
    status = nf90_inq_varid(id, varName, varid)
    call handle_err(status)
    status = nf90_inquire_variable(id, varid, ndims = nDims, dimids = dimIDs)
    call handle_err(status)
    do i = 1, NumDims
       status = nf90_inquire_dimension(id, dimIDs(i), len = ncSizes(i))
       call handle_err(status)
    end do

    ! Make sure array extents match!
    if (any(ncSizes /= sizes) .or. nDims /= NumDims) stop  &
         "Array dimensions are incorrect in read_values."

    ! Read the variable from the netCDF file
    status = nf90_get_var(id, varid, values)
    call handle_err(status)
  end subroutine rv3d

  ! ...........................................................................

  ! rv4d is identical to rv2d, but for a 4D array.
  subroutine rv4d(id, varName, values)
    integer,                  intent(in)  :: id
    character(len=*),         intent(in)  :: varName
    real, dimension(:,:,:,:), intent(out) :: values
    
    integer, parameter :: NumDims = 4
    integer, dimension(NF90_Max_Var_Dims) :: dimIDs   
    integer, dimension(NumDims)           :: sizes, ncsizes
    integer :: i, varid, status, nDims

    ! Determine array extents that wrf2gem expects
    sizes = shape(values)
    
    ! Get array extents from the netCDF file
    status = nf90_inq_varid(id, varName, varid)
    call handle_err(status)
    status = nf90_inquire_variable(id, varid, ndims = nDims, dimids = dimIDs)
    call handle_err(status)
    do i = 1, NumDims
       status = nf90_inquire_dimension(id, dimIDs(i), len = ncSizes(i))
       call handle_err(status)
    end do
    
    ! Make sure array extents match!
    if (any(ncSizes /= sizes) .or. nDims /= NumDims) stop  &
         "Array dimensions are incorrect in read_values."

    ! Read the variable from the netCDF file
    status = nf90_get_var(id, varid, values)
    call handle_err(status)
  end subroutine rv4d

  ! ===========================================================================

  ! Same as rv0dc, but accounting for the fact that the value doesn't change
  ! with time.
  subroutine rv0dc(id, varName, value)
    integer,          intent(in)  :: id
    character(len=*), intent(in)  :: varName
    real,             intent(out) :: value
    
    integer, parameter :: NumDims = 2
    real,    dimension(:,:), allocatable  :: temp
    integer, dimension(NF90_Max_Var_Dims) :: dimIDs
    integer, dimension(NumDims)           :: sizes, ncSizes
    integer :: status, varid, nt, i, nDims
    
    ! Get number of times in data file
    status = nf90_inq_dimid(id, "Time", varid)
    call handle_err(status)  
    status = nf90_inquire_dimension(id, varid, len=nt)
    call handle_err(status)

    ! Determine array extents that wrf2gem expects
    allocate(temp(1,nt), stat=status)
    call mem_error(status, 1, "rv0dc")
    sizes = shape(temp)
    
    ! Get array extents from the netCDF file
    status = nf90_inq_varid(id, varName, varid)
    call handle_err(status)
    status = nf90_inquire_variable(id, varid, ndims = nDims, dimids = dimIDs)
    call handle_err(status)
    do i = 1, NumDims
       status = nf90_inquire_dimension(id, dimIDs(i), len = ncSizes(i))
       call handle_err(status)
    end do

    ! Make sure array extents match!
    if (any(ncSizes /= sizes) .or. nDims /= NumDims) stop  &
         "Array dimensions are incorrect in read_values."
    
    ! Read the variable from the netCDF file
    status = nf90_get_var(id, varid, temp)
    call handle_err(status)
    value = temp(1,1)
  end subroutine rv0dc

  ! ...........................................................................

  ! Same as rv0dc, but for a 1D array.
  subroutine rv1dc(id, varName, values)
    integer,            intent(in)  :: id
    character(len=*),   intent(in)  :: varName
    real, dimension(:), intent(out) :: values
    
    integer, parameter :: NumDims = 2
    real,    dimension(:,:), allocatable  :: temp
    integer, dimension(NF90_Max_Var_Dims) :: dimIDs
    integer, dimension(NumDims)           :: sizes, ncSizes
    integer :: status, varid, nt, i, nDims
    
    ! Get number of times in data file
    status = nf90_inq_dimid(id, "Time", varid)
    call handle_err(status)  
    status = nf90_inquire_dimension(id, varid, len=nt)
    call handle_err(status)

    ! Determine array extents that wrf2gem expects
    allocate(temp(size(values),nt), stat=status)
    call mem_error(status, 1, "rv1dc")
    sizes = shape(temp)
    
    ! Get array extents from the netCDF file
    status = nf90_inq_varid(id, varName, varid)
    call handle_err(status)
    status = nf90_inquire_variable(id, varid, ndims = nDims, dimids = dimIDs)
    call handle_err(status)
    do i = 1, NumDims
       status = nf90_inquire_dimension(id, dimIDs(i), len = ncSizes(i))
       call handle_err(status)
    end do

    ! Make sure array extents match!
    if (any(ncSizes /= sizes) .or. nDims /= NumDims) stop  &
         "Array dimensions are incorrect in read_values."

    ! Read the variable from the netCDF file
    status = nf90_get_var(id, varid, temp)
    call handle_err(status)
    values = temp(:,1)
  end subroutine rv1dc

  ! ...........................................................................

  ! Same as rv0dc, but for a 2D array.
  subroutine rv2dc(id, varName, values)
    integer,              intent(in)  :: id
    character(len=*),     intent(in)  :: varName
    real, dimension(:,:), intent(out) :: values

    integer, parameter :: NumDims = 3
    real,    dimension(:,:,:), allocatable :: temp
    integer, dimension(NF90_Max_Var_Dims)  :: dimIDs
    integer, dimension(NumDims)            :: sizes, ncSizes
    integer :: status, varid, nt, i, nDims
    
    ! Get number of times in data file
    status = nf90_inq_dimid(id, "Time", varid)
    call handle_err(status)  
    status = nf90_inquire_dimension(id, varid, len=nt)
    call handle_err(status)

    ! Determine array extents that wrf2gem expects
    allocate(temp(size(values,1),size(values,2),nt), stat=status)
    call mem_error(status, 1, "rv2dc")
    sizes = shape(temp)
    
    ! Get array extents from the netCDF file
    status = nf90_inq_varid(id, varName, varid)
    call handle_err(status)
    status = nf90_inquire_variable(id, varid, ndims = nDims, dimids = dimIDs)
    call handle_err(status)
    do i = 1, NumDims
       status = nf90_inquire_dimension(id, dimIDs(i), len = ncSizes(i))
       call handle_err(status)
    end do

    ! Make sure array extents match!
    if (any(ncSizes /= sizes) .or. nDims /= NumDims) stop  &
         "Array dimensions are incorrect in read_values."

    ! Read the variable from the netCDF file
    status = nf90_get_var(id, varid, temp)
    call handle_err(status)
    values = temp(:,:,1)
  end subroutine rv2dc

  ! ...........................................................................

  ! Same as rv0dc, but for a 3D array.
  subroutine rv3dc(id, varName, values)
    integer,                intent(in)  :: id
    character(len=*),       intent(in)  :: varName
    real, dimension(:,:,:), intent(out) :: values

    integer, parameter :: NumDims = 4
    real,    dimension(:,:,:,:), allocatable :: temp
    integer, dimension(NF90_Max_Var_Dims)    :: dimIDs
    integer, dimension(NumDims)              :: sizes, ncSizes
    integer :: status, varid, nt, i, nDims
    
    ! Get number of times in data file
    status = nf90_inq_dimid(id, "Time", varid)
    call handle_err(status)  
    status = nf90_inquire_dimension(id, varid, len=nt)
    call handle_err(status)

    ! Determine array extents that wrf2gem expects
    allocate(temp(size(values,1),size(values,2),size(values,3),nt),stat=status)
    call mem_error(status, 1, "rv3dc")
    sizes = shape(temp)
    
    ! Get array extents from the netCDF file
    status = nf90_inq_varid(id, varName, varid)
    call handle_err(status)
    status = nf90_inquire_variable(id, varid, ndims = nDims, dimids = dimIDs)
    call handle_err(status)
    do i = 1, NumDims
       status = nf90_inquire_dimension(id, dimIDs(i), len = ncSizes(i))
       call handle_err(status)
    end do

    ! Make sure array extents match!
    if (any(ncSizes /= sizes) .or. nDims /= NumDims) stop  &
         "Array dimensions are incorrect in read_values."

    ! Read the variable from the netCDF file
    status = nf90_get_var(id, varid, temp)
    call handle_err(status)
    values = temp(:,:,:,1)
  end subroutine rv3dc

  ! ===========================================================================
  
  ! Given that there are p outputs per q hours, can the outputs be subdivided
  ! into blocks of r hours? This function answers that question.
  ! Yes = .true., and No = .false.
  ! Negative values for p, q, or r will also give .false. as output.
  logical function freq_divisible(p, q, r)
    integer, intent(in) :: p, q, r
    
    integer :: gcf, pReduced, qReduced, outputTimestep
    
    ! Make sure p, q, and r are all positive
    if (p <= 0 .or. q <= 0 .or. r <= 0) then
       freq_divisible = .false.
       return
    end if
    
    ! Reduce the fraction p/q to lowest terms.
    gcf = greatest_common_factor(p, q)
    pReduced = p / gcf
    qReduced = q / gcf
    
    ! First test: If r is divisible by q, then certainly there is a subdivision
    ! into r; in fact, there are pr blocks.
    ! Second test: If p doesn't divide into q, then the output timestep is not
    ! an integer, which is not allowed if r is not divisible by q.
    ! Third test: If r is a multiple of the output timestep, we've found the 
    ! subdivision we're looking for.
    
    if (mod(r,qReduced) == 0) then
       freq_divisible = .true.
    else if (mod(qReduced,pReduced) /= 0) then
       freq_divisible = .false.
    else
       outputTimestep = qReduced / pReduced
       if (mod(r,outputTimestep) == 0) then
          freq_divisible = .true.
       else
          freq_divisible = .false.
       end if
    end if
  end function freq_divisible

  ! ===========================================================================

  ! greatest_common_factor uses the Euclidean algorithm to compute the GCF of
  ! integers a and b. It returns zero if both a and b are zero. Otherwise, it
  ! returns a positive integer.
  recursive function greatest_common_factor(a, b) result (gcf)
    integer, intent(in) :: a, b
    integer             :: gcf
    
    if (b == 0) then
       gcf = a
    else
       gcf = greatest_common_factor(b, mod(a,b))
    end if
  end function greatest_common_factor
  
  ! ===========================================================================

  ! This function converts the input character string (ch) to an integer.
  integer function char_to_int(ch)
    character(len=*), intent(in) :: ch

    character(len=*), parameter :: digs = "0123456789 -", frmt = "(i9)"
    integer :: chLen, i, ind
    character(len=len(ch)) :: chl

    chLen = len(ch)
    chl = ch

    if (chLen < 1) then
       char_to_int = 0
       return
    end if

    ! Replace bad characters with spaces
    do i = 1, chLen
       if (scan(chl(i:i),digs) == 0) chl(i:i) = " "
    end do

    ! Make sure only spaces precede a minus sign
    ind = scan(chl,"-")
    chl(1:ind-1) = " "

    ! Make sure no additional minus signs occur
    do i = ind+1, chLen
       if (chl(i:i) == "-") chl(i:i) = " "
    end do

    read (chl, frmt) char_to_int
  end function char_to_int

  ! ===========================================================================
  
  ! This function converts a positive integer (num) into a character of 
  ! length 2. If num < 0, it's absolute value is used. If num has more than
  ! 2 digits, it is truncated such that only the ones and tens digits convert.
  ! A leading zero is added if necessary.
  character(len=2) function int_to_char2(num)
    integer, intent(in) :: num

    integer :: numl

    numl = abs(num)
    if (numl > 99) numl = mod(numl,100)
    
    write(unit=int_to_char2, fmt="(i2)") numl
    if (numl < 10) int_to_char2(1:1) = "0"
  end function int_to_char2
end module wrf2gem_subs
