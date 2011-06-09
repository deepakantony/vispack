SET(CMAKE_CXX_FLAGS -g)

# Add executable called "helloDemo" that is built from the source files
ADD_EXECUTABLE(tomofit tomofit.cxx projection.cxx)
# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(tomofit
vol
image 
jpeg 
tiff 
fitsio 
util 
matrix
param 
lapack
blas
F77
m
)


ADD_EXECUTABLE(testtomo testtomo.cxx projection.cxx)
TARGET_LINK_LIBRARIES(testtomo
vol
image 
jpeg 
tiff 
fitsio 
util 
matrix
param 
lapack
blas
F77
m
)




