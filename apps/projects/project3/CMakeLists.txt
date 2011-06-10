
# Make sure the compiler can find include files from our library
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib)

# Add executable called "helloDemo" that is built from the source files
ADD_EXECUTABLE(testparams testparams.cxx)

# Add executable called "helloDemo" that is built from the source files
ADD_EXECUTABLE(testtransform transform.cxx testtransform.cxx)

# Add project 3 executable
ADD_EXECUTABLE(project3 transform.cxx project3.cxx)

# Link the executable to the required library.
TARGET_LINK_LIBRARIES(project3
image 
jpeg 
tiff 
fitsio 
util 
matrix
param 
lapack
lapack1
lapack2
blas
F77
m
)

# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(testparams
image 
jpeg 
tiff 
fitsio 
util 
matrix
param 
lapack
lapack1
lapack2
blas
F77
m
)

# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(testtransform
image 
jpeg 
tiff 
fitsio 
util 
matrix
param 
lapack
lapack1
lapack2
blas
F77
m
)
