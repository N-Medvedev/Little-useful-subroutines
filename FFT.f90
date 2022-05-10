PROGRAM FFT
! Compilation:
! for DEBUG:
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /Qvec-report1 /fpp /Qtrapuv /dbglibs FFT.f90 -o FFT.exe /link /stack:9999999999 
! for RELEASE:
! ifort.exe /F9999999999 /O3 /Qipo /Qvec-report1 /fpp /Qopenmp /heap-arrays FFT.f90 -o FFT.exe /link /stack:9999999999 
!<===========================================


! Initiate modules
use Universal_constants

! For OpenMP external library
USE OMP_LIB

implicit none

! Variables used:
real(8), dimension(:,:), allocatable :: Array	! data from the file provided by user
complex(8), dimension(:), allocatable :: FFT_vec	! Vector after FFT
integer :: col_num	! number of column in the file containing the data to by FFT
complex(8), dimension(:), allocatable :: FFT_vector	! discrete Fourier-transformed data
complex(8), dimension(:), allocatable :: IFFT_vector	! inverse Fourier-transformed data
complex(8), dimension(:), allocatable :: FFT_Gauss	! Fourier transrofm of a Guassian function (used for deconvolution)
complex(8), dimension(:), allocatable :: G_matrix	! Fourier transrofm of a Guassian function (used for deconvolution)
integer :: FN, Reason, i
character(200) :: File_Name, temp
logical :: file_exists, file_opened, do_inverse, do_deconvolution
real(8) :: mu, sigma
real(8) :: start  ! to measure time spent
!---------------------------------------
! Get the current time to use later to count the wall-time elapsed:
start = omp_get_wtime() ! intrinsic subroutine of the library OMP_LIB
   
File_name = ''	! just to start
do_inverse = .false.	! by default, do direct DFT, not inverse
do_deconvolution = .false.	! by default, do not do deconvolution
mu = 0.0d0
sigma = 1.0d0

! Read the passed parameters: file name and the column number:
call get_add_data(File_name, col_num, do_inverse, do_deconvolution, mu, sigma)

call  print_elapsed_time(start, 'Reading options took:') ! below

if ( LEN(trim(adjustl(File_Name))) < 1 ) then 
   print*, 'File name is undefined'
   goto 2012
endif

! Read the data from the file into the Array:
call read_from_file(File_Name, Array)

call  print_elapsed_time(start, 'Reading data file took:') ! below

if (col_num > size(Array,2)) then
   write(temp,'(i6)') size(Array,2)
   print*, 'The file '//trim(adjustl(File_Name))//' has '//trim(adjustl(temp))//' columns only.'
   write(temp,'(i6)') col_num
   print*, 'Cannot use specified column number '//trim(adjustl(temp))//'. Program terminates.'
   goto 2012
endif

! if (.not.do_inverse) then ! direct DFT:
   ! Make fast Fourier transform with the data from the chosen column in the given file:
   if (size(Array(:,col_num)) < 500) then ! direct
      call make_DFT(Array(:,col_num),FFT_vector)
   else ! fast
      call make_FFT(Array(:,col_num),FFT_vector)
   endif
   File_name = 'OUTPUT_FFT.dat'
   
call  print_elapsed_time(start, 'Fourier transform took:') ! below

! Set the gaussian Fourier transform to be used for deconvolution
if (do_deconvolution) then ! do deconvolution
   if (.not.allocated(FFT_Gauss)) allocate(FFT_Gauss(size(FFT_vector)))
   if (.not.allocated(G_matrix)) allocate(G_matrix(size(FFT_vector)))
   FFT_Gauss = (0.0d0,0.0d0)
   do i = 1, size(FFT_Gauss)
      FFT_Gauss(i) = Gaussian(mu, sigma, dble(i)) ! below
!       write(*,'(i,es,es,es,es)') i, FFT_Gauss(i), FFT_of_Gaussian(mu, sigma, dble(i))
   enddo
   call make_FFT(DBLE(FFT_Gauss),FFT_Gauss)
   
!     call Wiener_deconvolution(G_matrix, FFT_Gauss, FFT_vector)	! below
!     FFT_vector(:) = FFT_vector(:)*G_matrix(:)	! in the Fourier domain, deconvolution is just a division

    where (FFT_Gauss(:) /= (0.0d0,0.0d0))
       FFT_vector(:) = FFT_vector(:)/FFT_Gauss(:)	! in the Fourier domain, deconvolution is just a division
    endwhere

   call  print_elapsed_time(start, 'Deconvolution took:') ! below
endif

! else ! Inverse DFT:
!    if (.not.allocated(FFT_vector)) allocate(FFT_vector(size(Array(:,col_num))))
!    FFT_vector(:) = DCMPLX(Array(:,col_num),0.0d0)
   ! Test if Fourier was done correctly by making inverse DFT:
   if (size(Array(:,col_num)) < 500) then ! direct
      call make_IDFT (FFT_vector, IFFT_vector)
   else ! fast
      call inverse_FFT (FFT_vector, IFFT_vector)
   endif
!    File_name = 'OUTPUT_IFFT.dat'
! endif

call  print_elapsed_time(start, 'Inverse Fourier  took:') ! below

open (newunit=FN, file=trim(adjustl(File_name)))

do i = 1, size(Array,1) 
     write(FN,'(i,es,es,es,es,es)') i, Array(i,col_num), FFT_vector(i), IFFT_vector(i)
 enddo

! if (.not.do_inverse) then ! direct DFT:
!    do i = 1, size(Array,1) 
!       write(FN,'(i,es,es,es)') i, Array(i,col_num), FFT_vector(i)
!    enddo
! else ! inverse DFT:
!    do i = 1, size(Array,1) 
!       write(FN,'(i,es,es,es)') i, Array(i,col_num), IFFT_vector(i)
!    enddo
! endif

call  print_elapsed_time(start, 'Writing output took:') ! below

close(FN)

2012 continue


!cccccccccccccccccccccc
 contains ! all subroutines:

! The Wiener deconvolution: 
subroutine Wiener_deconvolution(G, H, S, N)
   complex(8), dimension(:), intent(inout) :: G	! matrix
   complex(8), dimension(:), intent(in) :: H, S
   complex(8), dimension(:), intent(in), optional :: N 
   !-----------------------------
   if (present(N)) then
      where (H(:) /= (0.0d0,0.0d0) .and. S(:) /= (0.0d0,0.0d0) .and. N(:) /= (0.0d0,0.0d0))
         G(:) = Conjg(H(:))*S(:)/(ABS(H(:))*ABS(H(:))*S(:) + N(:))
      endwhere
   else
      where (H(:) /= (0.0d0,0.0d0) .and. S(:) /= (0.0d0,0.0d0))
         G(:) = Conjg(H(:))*S(:)/(ABS(H(:))*ABS(H(:))*S(:))
      endwhere
   endif
end subroutine Wiener_deconvolution


! FFT at the time x according to Gaussian shape:
pure function FFT_of_Gaussian(mu, sigma, ksi) result(FFT_Gaus) ! THIS DOES NOT WORK CORRECTLY, FOR SOME REASON...
   real(8), intent(in) :: mu, sigma, ksi
   complex(8) :: FFT_Gaus	! NOTE: it is only real function of mu=0, otherwise it is complex!
   FFT_Gaus = 1.0d0/sqrt(2.0d0*g_Pi)*exp(g_I*mu*ksi)*dexp(-ksi*ksi*sigma*sigma/2.0d0)
   !FFT_Gaus = exp(g_I*mu*ksi)*dexp(-(ksi)*(ksi)*sigma*sigma/2.0d0)
end function FFT_of_Gaussian


! at the time x according to Gaussian shape:
pure function Gaussian(mu, sigma, x) result(Gaus)
   real(8), intent(in) :: mu, sigma, x
   real(8) :: Gaus
!    Gaus = 1.0d0/(sqrt(2.0d0*g_Pi)*sigma)*dexp(-(x-mu)*(x-mu)/(2.0d0*sigma*sigma))
   Gaus = 1.0d0/sigma*dexp(-(x-mu)*(x-mu)/(2.0d0*sigma*sigma))
end function Gaussian

 

subroutine inverse_FFT (Vector, IFFT_vector) ! Fast Fourier Transform
   complex(8), dimension(:), intent(in) :: Vector	! input
   complex(8), dimension(:), allocatable, intent(out) :: IFFT_vector	! output: inverse discrete Fourier-transformed Vector
   !-----------------------------
   integer :: N, half_N, i, j, k
   complex(8), dimension(:), allocatable :: odd_vec, even_vec
   complex(8), dimension(:,:), allocatable :: exp_vec
   complex(8) :: exp_prefac,  part_1, part_2
   real(8) :: two_pi_N, fac
   N = size(Vector)
   if (.not.is_even(N)) N = N - 1 ! make it odd by skipping the last point, assuming this point is unimportant
   half_N = N/2
      
   if (.not.allocated(IFFT_vector)) allocate(IFFT_vector(N))	! the same dimension
   ! Sub-arrays used in FFT:
   allocate(odd_vec(half_N))
   allocate(even_vec(half_N))
   allocate(exp_vec(half_N,half_N))
   
   two_pi_N = 2.0d0*g_Pi/dble(half_N) ! factor enterring all exponents
   ! Set the odd and even parts of vector separately:
   do i = 0, half_N-1
      even_vec(i+1) =  Vector(2*i+1)	! X0
      odd_vec(i+1) = Vector(2*(i+1))	! X1
      fac = two_pi_N*dble(i)
      FORALL(j=0:half_N-1)  exp_vec(i+1,j+1) = exp(g_I*fac*dble(j))
   enddo
   
   ! Construct re FFT:
   IFFT_vector = (0.0d0,0.0d0)
   do i = 0, half_N-1
      part_1 = (0.0d0,0.0d0)
      part_2 = (0.0d0,0.0d0)
      do j = 0, half_N-1
          !j = k - (half_N-1)
          part_1 = part_1 + even_vec(j+1)*exp_vec(i+1,j+1)
          part_2 = part_2 + odd_vec(j+1)*exp_vec(i+1,j+1)
      enddo
      fac = 2.0d0*g_Pi/dble(N)*dble(i)
      exp_prefac =  exp(g_I*fac)
      ! Since Fourier transform is periodic function in time, we can shift all the outcome by half of a period:
      IFFT_vector(i+1) = part_1 + exp_prefac*part_2
      IFFT_vector(i+1+half_N) = part_1 - exp_prefac*part_2
   enddo
   IFFT_vector = IFFT_vector/dble(N)
end subroutine inverse_FFT 
 
 
subroutine make_IDFT(Vector, FFT_vector) ! straightforward Discrete Fourier Transform
   complex(8), dimension(:), intent(in) :: Vector	! input
   complex(8), dimension(:), allocatable, intent(out) :: FFT_vector	! output: discrete Fourier-transformed Vector
   !-----------------------------
   integer :: N, i, j
   real(8) :: two_pi_N, fac
   N = size(Vector)
!    if (.not.is_even(N)) N = N - 1 ! make it odd by skipping the last point, assuming this point is unimportant

   if (.not.allocated(FFT_vector)) allocate(FFT_vector(N))	! the same dimension
   FFT_vector = 0.0d0
   two_pi_N = 2.0d0*g_Pi/dble(N) ! factor enterring all exponents
   ! Set the odd and even parts of vector separately:
   do i = 0, N-1
      fac = two_pi_N*dble(i)
      do j = 0, N-1
         FFT_vector(i+1) = FFT_vector(i+1) + Vector(j+1) * exp(g_I*fac*dble(j))
      enddo
   enddo
   FFT_vector = FFT_vector/dble(N)
end subroutine make_IDFT
 
 
 
subroutine make_DFT(Vector, FFT_vector) ! straightforward Discrete Fourier Transform
   real(8), dimension(:), intent(in) :: Vector	! input
   complex(8), dimension(:), allocatable, intent(out) :: FFT_vector	! output: discrete Fourier-transformed Vector
   !-----------------------------
   integer :: N, i, j
   real(8) :: two_pi_N, fac
   N = size(Vector)
!    if (.not.is_even(N)) N = N - 1 ! make it odd by skipping the last point, assuming this point is unimportant

   if (.not.allocated(FFT_vector)) allocate(FFT_vector(N))	! the same dimension
   FFT_vector = 0.0d0
   two_pi_N = 2.0d0*g_Pi/dble(N) ! factor enterring all exponents
   ! Set the odd and even parts of vector separately:
   do i = 0, N-1
      fac = two_pi_N*dble(i)
      do j = 0, N-1
         FFT_vector(i+1) = FFT_vector(i+1) + Vector(j+1) * exp(-g_I*fac*dble(j))
      enddo
   enddo
end subroutine make_DFT

 
 
subroutine make_FFT(Vector, FFT_vector) ! Fast Fourier Transform
   real(8), dimension(:), intent(in) :: Vector	! input
   complex(8), dimension(:), allocatable, intent(out) :: FFT_vector	! output: discrete Fourier-transformed Vector
   !-----------------------------
   integer :: N, half_N, i, j, k
   real(8), dimension(:), allocatable :: odd_vec, even_vec
   complex(8), dimension(:,:), allocatable :: exp_vec
   complex(8) :: exp_prefac,  part_1, part_2
   real(8) :: two_pi_N, fac
   N = size(Vector)
   if (.not.is_even(N)) N = N - 1 ! make it odd by skipping the last point, assuming this point is unimportant
   half_N = N/2
      
   if (.not.allocated(FFT_vector)) allocate(FFT_vector(N))	! the same dimension
   ! Sub-arrays used in FFT:
   allocate(odd_vec(half_N))
   allocate(even_vec(half_N))
   allocate(exp_vec(half_N,half_N))
   
   two_pi_N = 2.0d0*g_Pi/dble(half_N) ! factor enterring all exponents
   ! Set the odd and even parts of vector separately:
   do i = 0, half_N-1
      even_vec(i+1) =  Vector(2*i+1)	! X0
      odd_vec(i+1) = Vector(2*(i+1))	! X1
!       even_vec(i+1) =  Vector(2*i+1)	! X0
!       odd_vec(i+1) = Vector(2*(i))	! X1
      fac = two_pi_N*dble(i)
      FORALL(j=0:half_N-1)  exp_vec(i+1,j+1) = exp(-g_I*fac*dble(j))
   enddo
   
   ! Construct re FFT:
   FFT_vector = (0.0d0,0.0d0)
   do i = 0, half_N-1
      part_1 = (0.0d0,0.0d0)
      part_2 = (0.0d0,0.0d0)
      do j = 0, half_N-1
          !j = k - (half_N-1)
          part_1 = part_1 + even_vec(j+1)*exp_vec(i+1,j+1)
          part_2 = part_2 + odd_vec(j+1)*exp_vec(i+1,j+1)
      enddo
      fac = 2.0d0*g_Pi/dble(N)*dble(i)
      exp_prefac =  exp(-g_I*fac)
      ! Since Fourier transform is periodic function in time, we can shift all the outcome by half of a period:
      FFT_vector(i+1) = part_1 + exp_prefac*part_2
      FFT_vector(i+1+half_N) = part_1 - exp_prefac*part_2
   enddo
end subroutine make_FFT


! This function determines whether given integer N is odd (returns .false.) or even (returns .true.)
pure function is_even(N) result(evn)
   logical :: evn
   integer, intent(in) :: N
   if (MOD(N,2) == 0) then
      evn = .true.
   else
      evn = .false.
   endif
end function is_even


subroutine read_from_file(File_name, Array)
   character(*), intent(in) :: File_Name 
   real(8), dimension(:,:), allocatable :: Array
   !--------------------------
   logical :: file_exists, file_opened
   integer :: FN, N, M, i, Reason
   !--------------------------
   inquire(file=trim(adjustl(File_name)),exist=file_exists) ! check if input file is there
   exsts:if (file_exists) then
      ! Input file:
      open (newunit=FN, file=trim(adjustl(File_name)), readonly)
      call Count_columns_in_file(FN, M)
      call Count_lines_in_file(FN, N)
      allocate(Array(N,M))
!       print*, 'Size:', size(Array,1), size(Array,2)
      Array = 0.0d0
      do i = 1, N
         read(FN,*,IOSTAT=Reason) Array(i,:)
         IF (Reason .GT. 0)  THEN ! ... something wrong ...
           write(*,'(a,i4,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
           goto 1010
         ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
           write(*,'(a,i4,a)') 'Problem reading input file in line ', i, ', unexpected END of file'
           goto 1010
         ELSE
             ! normal reading
         END IF
         !print*, Array(i,:)
      enddo
1010  close(FN)
   else exsts
      print*, 'File '//trim(adjustl(File_name))//' not found. Cannot proceed with FFT.'
   endif exsts
end subroutine read_from_file

 
 

! Reads additional data from the command line passed along:
subroutine get_add_data(File_name, col_num, do_inverse, do_deconvolution, mu, sigma)
   character(*), intent(inout) :: File_Name
   integer, intent(inout) :: col_num
   logical, intent(inout) :: do_inverse, do_deconvolution
   real(8) , intent(inout):: mu, sigma
   !---------------------------------------
   integer, dimension(:), allocatable :: stop_markers ! how many different optionas are passed?
   integer :: i, N
   character(500) :: string, char1
   
   string = '' ! to start with empty line
   ! Read all arguments passed to the code:
   do i = 1, iargc() ! for as many arguments as we passed
      call getarg(i, char1)
      string = trim(adjustl(string))//' '//trim(adjustl(char1)) ! collect them all into one line
   enddo
   
   ! Separate them into independent variables, find numbers of character which start new command:
   call parse_add_param_line(string, stop_markers) ! see below
   
   N = size(stop_markers)
   if ( (N < 1) .and. (LEN(trim(adjustl(string))) > 0) ) then 
      write(*,'(a)') '*************************************************************'
      print*, 'Could not interpret the passed string: ', trim(adjustl(string))
      print*, 'Continue with defaul options'
      write(*,'(a)') '*************************************************************'
   endif
   ! Interpret them, for each command - read it, interpret and execute:
   do i = 1, N
      if (i < N) then
         call react_to_command_passed(trim(adjustl(string(stop_markers(i):stop_markers(i+1)-1))), File_name, col_num, do_inverse, do_deconvolution, mu, sigma)
      else
         call react_to_command_passed(trim(adjustl(string(stop_markers(i):LEN(trim(adjustl(string)))+1))),  File_name, col_num, do_inverse, do_deconvolution, mu, sigma)
      endif
   enddo
end subroutine get_add_data


pure subroutine parse_add_param_line(string, stop_markers)
   character(*), intent(in) :: string ! the line we read from the 
   integer, dimension(:), allocatable, intent(inout) :: stop_markers ! how many different optionas are passed?
   integer :: leng, i, count_dash
   character(LEN(trim(adjustl(string)))) :: string_cur
   character(1) :: separator
   separator = '/'
   string_cur = trim(adjustl(string))
   leng = LEN(trim(adjustl(string_cur)))
   count_dash = 0
   do i = 1,leng
      if (trim(adjustl(string_cur(i:i))) == separator) then
         count_dash = count_dash + 1
      endif
   enddo
   if (allocated(stop_markers)) deallocate(stop_markers)
   allocate(stop_markers(count_dash))
   count_dash = 0
   do i = 1,leng
      if (trim(adjustl(string_cur(i:i))) == separator) then
         count_dash = count_dash + 1  ! next marker
         stop_markers(count_dash) = i ! save where the separation is
      endif
   enddo
end subroutine parse_add_param_line



! interpret command passed by the user, and act on it accordingly:
subroutine react_to_command_passed(string, File_name, col_num, do_inverse, do_deconvolution, mu, sigma)
   character(*), intent(in) :: string	! line with the command to analyze
   character(*), intent(inout) :: File_Name
   integer, intent(inout) :: col_num
   logical, intent(inout) :: do_inverse, do_deconvolution
   real(8), intent(inout) :: mu, sigma
   !-----------
   real(8) :: temp_real
   character(LEN(trim(adjustl(string)))) :: string_cur
   character(1000) :: printline
   character(100) :: text_part, char_temp
   integer :: leng, i, count_i, Reason, temp_int
   logical :: yes

   leng = LEN(trim(adjustl(string)))
   string_cur = trim(adjustl(string))
   
   Reason = 0
   printline = ''
   text_part = ''
   count_i = 1
   LIN:do i = 1, leng
      count_i = i
      if (trim(adjustl(string_cur(i:i))) == ' ') exit LIN
   enddo LIN
   text_part(1:count_i) = trim(adjustl(string_cur(2:count_i))) ! that's the text part of the command to be interpreted
   
   select case (trim(adjustl(text_part)))
   case ('FILE', 'file', 'File', 'File_name', 'FILE_NAME', 'file_name', 'DATA', 'Data', 'data')
      read(string_cur(count_i:),*,iostat=Reason) File_Name	! read in the file name
   case ('COL', 'Col', 'col', 'Column', 'COLUMN', 'column')
      read(string_cur(count_i:),*,iostat=Reason) col_num		! read in the number of column
   case ('Direct', 'DIRECT', 'direct', 'DFT', 'Dft', 'dft', 'FAST', 'Fast', 'fast', 'FFT', 'Fft', 'fft')
      do_inverse = .false. ! not inverse by direct Fourier
   case ('INVERSE', 'Inverse', 'inverse', 'IDFT', 'Idft', 'idft',  'IFFT', 'Ifft', 'ifft')
      do_inverse = .true. ! inverse Fourier
   case ('DECONVOLVE', 'Deconvolve', 'deconvolve', 'DECONVOLUTION', 'Deconvolution', 'deconvolution',  'Deconv', 'DECONV', 'deconv')
      do_deconvolution = .true. ! deconvolution with a gaussian function
      read(string_cur(count_i:),*,iostat=Reason) mu, sigma	! read in the gaussian mu and sigma for the function to deconvolve from
      sigma = sigma/2.3548d0	! convert FWHM into Gaussian-sigma parameter
   case ('info', 'INFO', 'Info')
      printline = trim(adjustl(printline))//' The code is written in 2017 by' !//NEW_LINE('c')
      printline = trim(adjustl(printline))//' Dr. Nikita Medvedev'//NEW_LINE('A')
      printline = trim(adjustl(printline))//' as a part of the research at'//NEW_LINE('A')
      printline = trim(adjustl(printline))//' Institute of Physics and Institute of Plasma Physics,'//NEW_LINE('A') 
      printline = trim(adjustl(printline))//' Academy of Science of Czech Republic,'//NEW_LINE('A') 
      printline = trim(adjustl(printline))//' Na Slovance 1999/2, 18221 Prague 8, Czech Republic'//NEW_LINE('A')
      printline = trim(adjustl(printline))//' Should you have any questions, contact the author: nikita.medvedev@fzu.cz'//NEW_LINE('A')
      printline = trim(adjustl(printline))//' Or by private email: n.a.medvedev@gmail.com'//NEW_LINE('A')
      write(*,'(a)') '*************************************************************'
      write(*,'(a)') trim(adjustl(printline))
      write(*,'(a)') '*************************************************************'
   case default
      write(*,'(a)') '*************************************************************'
      print*, 'Could not interpret the following command: ', trim(adjustl(text_part))
      print*, 'Continue calculations with the default parameters'
      write(*,'(a)') '*************************************************************'
   end select
end subroutine react_to_command_passed


subroutine Count_columns_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of columns in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    real(8) temp
    character(1000) temp_ch
    integer i, Reason
    integer :: temp_i
    if (present(skip_lines)) then ! in case you want to skip some comment lines and count column only in a line with data
       do i=1,skip_lines
          read(File_num,*, end=601) 
       enddo
       601 continue
    endif

    read(File_num,'(a)', IOSTAT=Reason) temp_ch ! count columns in this line
    N = number_of_columns(trim(adjustl(temp_ch))) ! see below

    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
end subroutine Count_columns_in_file


pure function number_of_columns(line)
   integer :: number_of_columns
   character(*), intent(in) :: line
   integer i, n, toks
   logical :: same_space
   same_space = .false.
   i = 0
   n = len(line)
   number_of_columns = 0
   do while(i < n) ! scan through all the line
      i = i + 1
      selectcase (line(i:I))
      case (' ', '	') ! space or tab can be a separator between the columns
         if (.not.same_space) number_of_columns = number_of_columns + 1
         same_space = .true. ! in case columns are separated by more than one space or tab
      case default ! column data themselves, not a space inbetween
         same_space = .false.
      endselect
   enddo
   number_of_columns = number_of_columns + 1	! number of columns is by 1 more than number of spaces inbetween
end function number_of_columns 




subroutine Count_lines_in_file(File_num, N, skip_lines)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer, intent(in), optional :: skip_lines ! if you want to start not from the first line
    integer i
    if (present(skip_lines)) then ! in case you want to skip some comment lines and count only lines with data
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


subroutine print_elapsed_time(start, text)
   real(8), intent(inout) :: start ! starting time
   character(*), intent(in), optional :: text	! text to print out
   character(200) :: string
   real(8) :: finish
   finish = omp_get_wtime ( ) ! OMP provided subroutine
   if (present(text)) then
      string = trim(adjustl(text))
   else
      string = 'Elapsed time: '
   endif
   write(*,'(a,es12.2,a)')  trim(adjustl(string))//'	', finish - start, ' [sec]'
   start = finish
end subroutine print_elapsed_time


END PROGRAM FFT