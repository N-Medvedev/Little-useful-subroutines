Program Convolve
! Compilation:
! for DEBUG:
! ifx.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /fpp /Qtrapuv /dbglibs  Convolution_with_Gaussian.f90 -o Gauss_Convolve.exe /link /stack:9999999999 
! for RELEASE:
! ifx.exe /F9999999999 /O3 /fpp  Convolution_with_Gaussian.f90 -o Gauss_Convolve.exe /link /stack:9999999999 
!<===========================================

integer i, j, FN, FN2, N, M, Reason, N_skip, j_run
real(8), dimension(:), allocatable :: Read_data
real(8), dimension(:,:), allocatable :: Spectr
real(8), dimension(:,:), allocatable :: Conv_Spectr
real(8) temp, temp2(3), Gaus, t_g, ave_val, dt
character(100) :: File_name, File_name2, First_line
character(500) :: Read_line
!character(100), dimension(:), allocatable :: First_line
logical file_exists

t_g = 1.0d0    ! default gauss width for convolution

call get_add_data(t_g, File_name)	! below


inquire(file=trim(adjustl(File_name)),exist=file_exists) ! check if input file is there
exists:if (file_exists) then
   ! Input file:
   open (newunit=FN, file=trim(adjustl(File_name)), readonly)
   call Count_lines_in_file(FN, N)
   if (N > 2) then 
      N_skip = 2
	  call Count_columns_in_file(FN, M, skip_lines=2)
   elseif(N > 1) then 
	  N_skip = 1
	  call Count_columns_in_file(FN, M, skip_lines=1)
   else
	  N_skip = 0
      call Count_columns_in_file(FN, M)
   endif
   allocate(Spectr(N,M))
   allocate(Conv_Spectr(N,M))
   allocate(Read_data(M), source = 0.0d0)
   !allocate(First_line(M))
   print*, 'Size:', size(Spectr,1), size(Spectr,2)
   Spectr = 0.0d0
   Conv_Spectr = 0.0d0
   
   ! Output file:
   File_name2 = 'Convolved_'//trim(adjustl(File_name))
   open (newunit=FN2, file=trim(adjustl(File_name2))) ! output file
   
   !read(FN,'(a)',IOSTAT=Reason) First_line
   !read(FN,*,IOSTAT=Reason)
   do i = 1, N
	  read(FN,'(a)',IOSTAT=Reason) Read_line
      IF (Reason .GT. 0)  THEN ! ... something wrong ...
        write(*,'(a,i4,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
        goto 1010
      ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
        write(*,'(a,i4,a)') 'Problem reading input file in line ', i, ', unexpected END of file'
        goto 1010
      ELSE   ! normal reading
      END IF
	  
	  if (trim(adjustl(Read_line(1:1))) == '#') then 
		write(FN2,'(a)')  trim(adjustl(Read_line))		! transfer comment line to output file
	  else
		read(Read_line,*,IOSTAT=Reason) Read_data
		Spectr(i,:) = Read_data(:)
	  endif
	  
      !print*, Spectr(i,1:2)
   enddo
   close(FN)
   
   Conv_Spectr(:,1) = Spectr(:,1)
   dt = Conv_Spectr(2+N_skip,1) - Conv_Spectr(1+N_skip,1)		! time step, assuming constant

   ! Convolution:
   do i = 1+N_skip,N ! all time steps
      do j_run = i-floor(3.0d0*t_g/dt), i+ceiling(3.0d0*t_g/dt)  ! convolve
        if (j_run < 1+N_skip) then
			j = N_skip + abs(j_run)	! count back from the bottom
		elseif (j_run > N) then
			j = j_run - N	! count back from the top
		else
			j = j_run	! count normal
		endif
		call Gaussian(mu=Spectr(i,1), sigma=t_g, x=Spectr(j,1), Gaus=Gaus)  ! Gaussian weight
        Conv_Spectr(i,2:) = Conv_Spectr(i,2:) + Spectr(j,2:)*Gaus*dt   ! Convolution of all arrays
        !write(*,'(f,$)') Spectr(j,1), Spectr(j,2), Gaus, Conv_Spectr(i,2)
      enddo
	  Read_data(:) = Conv_Spectr(i,:)
	  !print*, i, Conv_Spectr(i,1), Conv_Spectr(i,2)
      WRITE(FN2,'(e,$)') Read_data(:)
      write(FN2,'(a)') ' '
      !pause 'TEST PAUSE'
   enddo
   close(FN2)
else
   print*, 'File ', trim(adjustl(File_name)), ' not found.'
endif exists


1010 continue
print*, 'I am done convoling'
!pause 'PROGRAM FINISHED'
!STOP

Contains



! Reads additional data from the command line passed along:
subroutine get_add_data(Sigma, File_name)
   real(8), intent(out) :: Sigma
   character(*), intent(inout) :: File_name
   !---------------------------------------
   integer, dimension(:), allocatable :: stop_markers ! how many different options are passed?
   integer :: i, count_dash
   character(500) :: string, char1
   logical :: read_well  ! did we read ok?
   
   read_well = .true. ! to start with
   
   string = '' ! to start with empty line
   ! Read all arguments passed to the code:
   do i = 1, iargc() ! for as many arguments as we passed
      call getarg(i, char1)
      string = trim(adjustl(string))//' '//trim(adjustl(char1)) ! collect them all into one line
   enddo
   
   read(string, *, IOSTAT=Reason) sigma, File_name
   call read_file_here(Reason, 1, read_well) ! below
   
   if (.not.read_well) then
      write(*,'(a)') '***************************************************************************'
      write(*,'(a)') 'Could not interpret the passed string: '//trim(adjustl(string))
	  write(*,'(a)') 'The call format must be as follows: Gauss_Convolve.exe number filename'
	  write(*,'(a)') 'Parameters "number" or "filename" were wrongly specified'
	  write(*,'(a)') '***************************************************************************'
	  write(*,'(a)') 'Add the parameters maually.'
	  write(*,'(a)') 'Enter the convolution width parameter:'
	  read(*,*) sigma
	  write(*,'(a)') 'Enter the name of the file to convolve'
	  read(*,*) File_name
	  write(*,'(a)') '***************************************************************************'
   endif
end subroutine get_add_data



subroutine read_file_here(Reason, i, read_well)
   integer, intent(in) :: Reason    ! file number where to read from
   integer, intent(in) :: i      ! line number
   logical, intent(inout) :: read_well  ! did we read ok?
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       write(*,'(a,i3,a)') 'Problem reading input line ', i, ', wrong type of variable'
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       write(*,'(a,i3,a)') 'Problem reading input line ', i, ', unexpected END'
       read_well = .false.
   ELSE   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   END IF
end subroutine read_file_here
 

subroutine Gaussian(mu, sigma, x, Gaus) ! at the time x according to Gaussian shape
   real(8), intent(in) :: mu, sigma, x
   real(8), intent(out) :: Gaus
   real(8), parameter :: g_Pi = 3.1415926535897932384626433832795d0
   Gaus = 1.0d0/(sqrt(2.0d0*g_Pi)*sigma)*dexp(-0.5d0*((x-mu)/sigma)**2)
end subroutine Gaussian


subroutine Count_columns_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of columns in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    real(8) temp
    integer i, Reason
    if (present(skip_lines)) then
       do i=1,skip_lines
          read(File_num,*, end=601) 
       enddo
       601 continue
    endif
    i = 0
    do
        read(File_num, '(e)', advance='no', IOSTAT=Reason) temp
        if (Reason .NE. 0) exit
        i = i + 1
    enddo
    602 continue
    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    N = i
end subroutine Count_columns_in_file


subroutine Count_lines_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    integer i
    if (present(skip_lines)) then
       do i=1,skip_lines
          read(File_num,*, end=604) 
       enddo
       604 continue
    endif
    i = 0
    do
        read(File_num,*, end=603)
        i = i + 1
    enddo
    603 continue
    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    N = i
end subroutine Count_lines_in_file

End Program Convolve