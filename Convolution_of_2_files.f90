program Convolution
! Compilation:
! for DEBUG:
! ifort.exe /debug:all /check:all /check:bounds /fp:precise /fpe-all:0 /Qopenmp /Qftz- /Qfp-stack-check /Od /Zi /traceback /gen-interfaces /warn:all /warn:nounused /Qvec-report1 /fpp /Qtrapuv /dbglibs Convolution_of_2_files.f90 -o Convolve.exe /link /stack:9999999999 
! for RELEASE:
! ifort.exe /F9999999999 /O3 /Qipo /Qvec-report1 /fpp /Qopenmp /heap-arrays Convolution_of_2_files.f90 -o Convolve.exe /link /stack:9999999999 
!<===========================================

implicit none

real(8), dimension(:,:), allocatable :: F_to_conv, F_conv_with
real(8), dimension(:,:), allocatable :: F_out
real(8) :: grid_size, grid_dx, grid_start, F1, F2, f_diff, F_0, Norm_F1, Norm_F2, Norm_out
integer :: N_f1, N_f2, N_out, i, j, k
integer :: FN, File_1_col_x, File_1_col_y, File_2_col_x, File_2_col_y
character(200) :: File_name_1, File_name_2, Output_file
logical :: read_well, file_exists, file_opened
!<===========================================

FN = 9999 ! output file
Output_file = 'OUTPUT_Convolved.dat'

!<===========================================
! get the files and parameters:
call get_add_data(File_name_1, File_name_2, File_1_col_x, File_1_col_y, File_2_col_x, File_2_col_y, read_well)	! below

!<===========================================
! If user did not provide necessary parameters, inform the user on the format of the program call:
if (.not.read_well) then
   write(*,'(a)') '***************************************************************************'
   write(*,'(a)') 'Call the program with the following inputs:'
   write(*,'(a)') 'Convolve.exe  /File_#1 col_1 col_2  /File_#2 col_1 col_2'
   write(*,'(a)') 'File_#1 is the function to be convolved'
   write(*,'(a)') 'File_#2 is the function to be convolved with (will be normalized to 1)'
   write(*,'(a)') 'col_1 is the column with the x axis (grid)'
   write(*,'(a)') 'col_2 is the column with the y axis (data)'
   write(*,'(a)') 'Entry separator " / " must be used before each file name'
   write(*,'(a)') '***************************************************************************'
   goto 2016
else ! continue the program
   ! First file:
   inquire(file=trim(adjustl(File_name_1)),exist=file_exists)    ! check if input file exists
   if (file_exists) then
      write(*,'(a,a)') 'File with function to be convolved: ', trim(adjustl(File_name_1))
      call read_file_function(File_name_1, File_1_col_x, File_1_col_y, F_to_conv, read_well) ! below
      if (.not.read_well) goto 2016
   else 
      write(*,'(a)') '***************************************************************************'
      write(*,'(a,a,a)') 'File ', trim(adjustl(File_name_1)), ' not found'
      write(*,'(a)') '***************************************************************************'
      goto 2016
   endif
   
   ! Second file:
   inquire(file=trim(adjustl(File_name_2)),exist=file_exists)    ! check if input file exists
   if (file_exists) then
      write(*,'(a,a)') 'File with function to be convolved with (normalized): ', trim(adjustl(File_name_2))
      call read_file_function(File_name_2, File_2_col_x, File_2_col_y, F_conv_with, read_well) ! below
      if (.not.read_well) goto 2016
   else
      write(*,'(a)') '***************************************************************************'
      write(*,'(a,a,a)') 'File ', trim(adjustl(File_name_2)), ' not found'
      write(*,'(a)') '***************************************************************************'
      goto 2016
   endif
endif

!<===========================================
! Grid sizes of input functions:
N_f1 = size(F_to_conv,2)
N_f2 = size(F_conv_with,2)
! Find grid parameters for the output function:
grid_start = min( F_to_conv(1,1), F_conv_with(1,1) )
grid_size = max((F_to_conv(1,N_f1)-F_to_conv(1,1)), (F_conv_with(1,N_f2)-F_conv_with(1,1)))
grid_dx = min( (F_to_conv(1,2)-F_to_conv(1,1)), (F_conv_with(1,2)-F_conv_with(1,1)) )
N_out = CEILING(grid_size/grid_dx)
allocate(F_out(2,N_out))
F_out = 0.0d0
! Set grid for the output function:
do i = 1, N_out
   F_out(1,i) = grid_start + dble(i-1)*grid_dx
enddo

!<===========================================
! Normalization of the input function:
F_0 = 0.0d0
Norm_F2 = 0.0d0
do i = 1, N_f2
   if (i > 1) then
      Norm_F2 = Norm_F2 + (F_conv_with(2,i) + F_conv_with(2,i-1))/2.0d0*(F_conv_with(1,i) - F_conv_with(1,i-1))
   else
      Norm_F2 = Norm_F2 + (F_conv_with(2,i) + 0.0d0)/2.0d0*(F_conv_with(1,i) - 0.0d0)
   endif
enddo
F_conv_with(2,:) = F_conv_with(2,:)/Norm_F2	! normalize to 1

!<===========================================
! Convolve the two functions:
! Convolution (t) = int ( F1(tau) * F2(t - tau) d tau)
! https://en.wikipedia.org/wiki/Convolution#Definition
do i = 1, N_out
   F_0 = 0.0d0 ! integrant on the last point
   do j = 1, N_out
      ! Interpolate first function onto the grid point:
      if (F_out(1,j) <= F_to_conv(1,1)) then
         F1 = F_to_conv(2,1)
      elseif (F_out(1,j) >= F_to_conv(1,N_f1)) then
         F1 = F_to_conv(2,N_f1)
      else
         call Find_in_monotonous_1D_array(F_to_conv(1,:), F_out(1,j), k)	! below
         call Linear_interpolation(F_to_conv(2,k-1), F_to_conv(2,k), F_to_conv(1,k-1), F_to_conv(1,k), F_out(1,j), F1)	! below
      endif
      
      f_diff = F_out(1,i) - F_out(1,j)
      if ( f_diff  <= F_conv_with(1,1) ) then
         F2 = F_conv_with(2,1)
      elseif  ( f_diff  >= F_conv_with(1,N_f2) ) then
         F2 = F_conv_with(2,N_f2)
      else
         call Find_in_monotonous_1D_array(F_conv_with(1,:), f_diff, k)	! below
         call Linear_interpolation(F_conv_with(2,k-1), F_conv_with(2,k), F_conv_with(1,k-1), F_conv_with(1,k), f_diff, F2)	! below
      endif

!       if (F1 < 0.0d0) print*, 'F1', F_to_conv(2,k-1), F_to_conv(2,k), F_to_conv(1,k-1), F_to_conv(1,k), F_out(1,j), F1
!       if (F2 < 0.0d0) print*, 'F2', F_conv_with(2,k-1), F_conv_with(2,k), F_conv_with(1,k-1), F_conv_with(1,k), f_diff, F2
      
      ! Integrate:
      F_out(2,i) = F_out(2,i) + (F1*F2 + F_0)/2.0d0*grid_dx
      
      F_0 = F1*F2	! save the last point
   enddo ! j
enddo ! i
!<===========================================

! Check normalization after convolution:
F_0 = 0.0d0
Norm_F1 = 0.0d0
do i = 1, N_f1
   if (i > 1) then
      Norm_F1 = Norm_F1 + (F_to_conv(2,i) + F_to_conv(2,i-1))/2.0d0*(F_to_conv(1,i) - F_to_conv(1,i-1))
   else
      Norm_F1 = Norm_F1 + (F_to_conv(2,i) + 0.0d0)/2.0d0*(F_to_conv(1,i) - 0.0d0)
   endif
enddo
write(*,'(a,es)') 'Normalization of the original function: ', Norm_F1

F_0 = 0.0d0
Norm_out = 0.0d0
do i = 1, N_out
   if (i > 1) then
      Norm_out = Norm_out + (F_out(2,i) + F_out(2,i-1))/2.0d0*(F_out(1,i) - F_out(1,i-1))
   else
      Norm_out = Norm_out + (F_out(2,i) + 0.0d0)/2.0d0*(F_out(1,i) - 0.0d0)
   endif
enddo
write(*,'(a,es)') 'Normalization of the convolved function: ', Norm_out
write(*,'(a,es,a)') 'Relative difference between them: ', ABS(Norm_F1-Norm_out)/Norm_F1*100.0d0, ' %'

!<===========================================

! OUTPUT file:
open(unit = FN, FILE = trim(adjustl(Output_file)))   ! if yes, open it and read
! write output into the file:
do i = 1, N_out
   write(FN,'(f,es)') F_out(:,i)
!    write(*,'(f,es)') F_out(:,i)
enddo

write(*,'(a)') 'Output file was created: '//trim(adjustl(Output_file))

2016 inquire(unit=FN,opened=file_opened)    ! check if this file is opened
if (file_opened) close(FN)             ! and if it is, close it


STOP

!<===========================================
 contains
 

subroutine Linear_interpolation(y1, y2, x1, x2, x, y)
   real(8) y1, y2, x1, x2, x, y
   y = y1 + (y2-y1)/(x2-x1)*(x-x1)
end subroutine Linear_interpolation

 
subroutine read_file_function(File_name, col_x, col_y, funct, read_well)
   character(*), intent(in) :: File_name
   integer, intent(in) :: col_x, col_y
   real(8), dimension(:,:), allocatable, intent(inout) :: funct
   logical, intent(inout) :: read_well
   !-------------------
   real(8), dimension(:), allocatable :: temp_to_read
   integer :: N_cols, N_lines, i, Reason
   integer :: FN ! file number
   FN = 100
   open(unit = FN, FILE = trim(adjustl(File_name)), status = 'old', readonly)   ! if yes, open it and read
   ! Get the file format:
   call Count_lines_in_file(FN, N_lines) ! below
   call Count_columns_in_file(FN, N_cols) ! below
   ! Allocate the array with function:
   allocate(funct(2,N_lines))
   ! Allocate temporary array to read data to
   allocate(temp_to_read(N_cols))
   ! Read file line by line:
   do i = 1, N_lines
      read(FN,*, IOSTAT=Reason) temp_to_read(:)
      call read_file_here(Reason, i, read_well) ! below
      if (.not.read_well) then
         write(*,'(a)') 'Problem in the file '//trim(adjustl(File_name))
         goto 2013
      endif
      ! Now read into array the x and y data:
      funct(1,i) = temp_to_read(col_x)
      funct(2,i) = temp_to_read(col_y)
!       print*, 'TEST ', funct(:,i)
   enddo
2013   inquire(unit=FN,opened=file_opened)    ! check if this file is opened
   if (file_opened) close(FN)             ! and if it is, close it
end subroutine 
 

! Reads additional data from the command line passed along:
subroutine get_add_data(File_name_1, File_name_2, File_1_col_x, File_1_col_y, File_2_col_x, File_2_col_y, read_well)
   integer, intent(inout) :: File_1_col_x, File_1_col_y, File_2_col_x, File_2_col_y
   character(*), intent(inout) :: File_name_1, File_name_2
   logical, intent(inout) :: read_well  ! did we read ok?
   !---------------------------------------
   integer, dimension(:), allocatable :: stop_markers ! how many different options are passed?
   integer :: i, count_dash
   character(500) :: string, char1
   read_well = .true. ! to start with
   
   string = '' ! to start with empty line
   ! Read all arguments passed to the code:
   do i = 1, iargc() ! for as many arguments as we passed
      call getarg(i, char1)
      string = trim(adjustl(string))//' '//trim(adjustl(char1)) ! collect them all into one line
   enddo
   
   ! Separate them into independent variables, find numbers of character which start new command:
   call parse_add_param_line(string, count_dash, stop_markers) ! see below
   
   if (count_dash == 2) then
      ! Interpret them, for each command - read it, interpret and execute:
      ! File #1:
      call react_to_command_passed(trim(adjustl(string(stop_markers(1):stop_markers(2)))), File_name_1, File_1_col_x, File_1_col_y, read_well) ! below
      if (.not.read_well) goto 2017
      ! File #2:
      call react_to_command_passed(trim(adjustl(string(stop_markers(2):LEN(trim(adjustl(string)))+1))), File_name_2, File_2_col_x, File_2_col_y, read_well) ! below
      if (.not.read_well) goto 2017
   else ! The input format is wrong, nothing else to do here
      read_well = .false. ! didn't read the input correctly
   endif
   
2017   if (.not.read_well) then
      write(*,'(a)') '***************************************************************************'
      write(*,'(a)') 'Could not interpret the passed string: '//trim(adjustl(string))
   endif
end subroutine get_add_data
 

subroutine react_to_command_passed(string, File_name, col_x, col_y, read_well)
   character(*), intent(in) :: string
   character(*), intent(out) :: File_name
   integer, intent(out) :: col_x, col_y
   logical, intent(inout) :: read_well  ! did we read ok?
   character(1) :: temp_ch
   integer :: spaces_count, temp_i, i, Reason
   spaces_count = 0
   File_name = ''
   col_x = 0
   col_y = 0
   do i = 1,LEN(TRIM(adjustl(string)))
      read(string(i:i),'(a)', IOSTAT=Reason) temp_ch
      call read_file_here(Reason, i, read_well) ! below
      if (.not.read_well) goto 2012
      if (temp_ch == ' ') then
         spaces_count = spaces_count + 1
      else
         selectcase(spaces_count)
         case (1) ! col_x
            read(temp_ch,*,IOSTAT=Reason) temp_i
            call read_file_here(Reason, i, read_well) ! below
            if (.not.read_well) goto 2012
            col_x = temp_i + col_x*10
         case (2) ! col y
            read(temp_ch,*,IOSTAT=Reason) temp_i
            call read_file_here(Reason, i, read_well) ! below
            if (.not.read_well) goto 2012
            col_y = temp_i + col_y*10
         case default ! file name
            if (temp_ch(1:1) /= '/') File_name(i-1:i-1) = temp_ch
         endselect
      endif
   enddo

2012   if (spaces_count < 2) then
      read_well = .false.
   endif
end subroutine react_to_command_passed


subroutine read_file_here(Reason, i, read_well)
   integer, intent(in) :: Reason    ! file number where to read from
   integer, intent(in) :: i      ! line number
   logical, intent(inout) :: read_well  ! did we read ok?
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', unexpected END of file'
       read_well = .false.
   ELSE   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   END IF
end subroutine read_file_here

 
pure subroutine parse_add_param_line(string, count_dash, stop_markers)
   character(*), intent(in) :: string ! the line we read from the 
   integer, dimension(:), allocatable, intent(inout) :: stop_markers ! how many different optionas are passed?
   integer, intent(inout) :: count_dash	! how many separators
   integer :: leng, i
   character(LEN(trim(adjustl(string)))) :: string_cur
   character(1) :: separator
   
   separator = '/'	! this is the mark how to separate new entry
   
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
   logical :: same_space, there_is_number
   same_space = .false.
   there_is_number = .false.
   i = 0
   n = len(line)
   number_of_columns = 0
   do while(i < n) ! scan through all the line
      i = i + 1
      selectcase (line(i:I))
      case (' ', '	') ! space or tab can be a separator between the columns
         if (.not.same_space) number_of_columns = number_of_columns + 1
         same_space = .true. ! in case columns are separated by more than one space or tab
         there_is_number = .false.	! so far it is only empty space, there is no number there
      case default ! column data themselves, not a space in-between
         same_space = .false.
         there_is_number = .true.	! there is a number to be read, so number of columns must be increased by one, see below
      endselect
   enddo
   if (there_is_number) number_of_columns = number_of_columns + 1	! number of columns is by 1 more than number of spaces in-between
end function number_of_columns 



subroutine Find_in_monotonous_1D_array(Array, Value0, Number)
   REAL(8), dimension(:), INTENT(in) :: Array ! in which we are looking for the Value
   REAL(8), INTENT(in) :: Value0   ! to be found in the array as near as possible
   integer, INTENT(out) :: Number ! number of the element which we are looking for 
   integer i, N, i_cur, i_1, i_2, coun
   real(8) temp_val, val_1, val_2
   
   N = size(Array)
   i_1 = 1
   val_1 = Array(i_1)
   i_2 = N
   val_2 = Array(i_2)
   i_cur = FLOOR((i_1+i_2)/2.0)
   temp_val = Array(i_cur)
   if (isnan(Value0)) then
        print*, 'The subroutine Find_in_monotonous_1D_array'
        print*, 'cannot proceed, the value of Value0 is', Value0
        write(*, '(f,f,f,f)') Value0, Array(i_cur), Array(i_1), Array(i_2)
        pause 'STOPPED WORKING...'
   else
       if (Value0 .LT. Array(1)) then ! it's the first value, no need to search
           i_cur = 0
       else if (Value0 .GE. Array(N)) then ! it's the last value, no need to search
           i_cur = N-1
       else
           coun = 0
           do ! until the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if ((Value0 .GE. Array(i_cur)) .AND. (Value0 .LE. Array(i_cur+1))) exit ! when the Value is in between Array(i_cur) and Array(i_cur+1) => we found i_cur
                if (temp_val .LE. Value0) then
                   i_1 = i_cur
                   val_1 = Array(i_1)
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(i_cur)
                else
                   i_2 = i_cur
                   val_2 = temp_val
                   i_cur = FLOOR((i_1+i_2)/2.0)
                   temp_val = Array(i_cur)
                endif
                coun = coun + 1
                if (coun .GT. 1e3) then
                    print*, 'PROBLEM WITH CONVERGANCE IN'
                    print*, 'Find_in_monotonous_1D_array', coun
                    write(*, '(f,f,f,f)') Value0, Array(i_cur), Array(i_1), Array(i_2)
                    pause 'STOPPED WORKING...'
                endif
           enddo
       endif
   endif    ! isnan
   Number = i_cur+1
end subroutine Find_in_monotonous_1D_array




end program Convolution