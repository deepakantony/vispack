# Make sure the compiler can find include files from our library
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib)

# Add executable called "helloDemo" that is built from the source files
ADD_EXECUTABLE(project4 fourn.cxx transform.cxx mosaic.cxx project4.cxx)

# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(project4
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

