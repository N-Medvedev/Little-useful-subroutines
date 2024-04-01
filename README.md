# Little useful subroutines
Little utilities, not belonging to a particular project but useful in certain cases

 ## Disclaimer

Although we endeavour to ensure that all the codes and results delivered are correct, no warranty is given as to their accuracy. We assume no responsibility for possible errors or omissions. We shall not be liable for any damage arising from the use of these codes or their parts or any results produced with them, or from any action or decision taken as a result of using these codes or any related material.

These codes are distributed as is for non-commercial peaceful purposes only, such as research and education. It is explicitly prohibited to use the codes, their parts, their results or any related material for military-related and other than peaceful purposes. By using these codes or their materials, you agree with these terms and conditions.

## Universal_constants.f90
Contains universal constants that can be used in a code as global variables


## Convolution_of_2_files.f90
The code convolves two arbitrary data arrays read from user-defined files. To execute, call the compiled file (e.g. Convolve.exe) in the following format: 
* Convolve.exe  /File_#1 col_1 col_2  /File_#2 col_1 col_2
with the following key-words:
* File_#1 is the function to be convolved
* File_#2 is the function to be convolved with (will be normalized to 1)
* col_1 is the column with the x axis (grid)
* col_2 is the column with the y axis (data)

Note that Entry separator " / " must be used before each file name


## FFT.f90
The code performes a Fourier (or inverse Fourier) transform on the data from a user-specified file. To execute, call the compiled code (e.g. FFT.exe) with the following key-words:
* FILE -- to specify the filename with the data
* COL -- to specify the column in the file to be used for Fourier transform
* DIRECT -- to specify that it should be Fourier (not inverse Fourier) transform
* INVERSE  -- to specify that it should be inverse Fourier transform
* DECONVOLVE -- to specify if deconvolution of the data with a Gaussian function should be performed (occasionally useful feature)



## Convolution_with_Gaussian.f90
The code convolves an arbitrary data read from user-defined file  with a Gaussian function of a given width. To execute, call the compiled file (e.g. Gaussian_Convolve.exe) in the following format:
* Convolve.exe Sigma FileName
with the following key-words:
* Sigma is the Gaussian width (real)
* Filename is the file wit the data, assuming the first column is the time; all the other columns will be convolved
