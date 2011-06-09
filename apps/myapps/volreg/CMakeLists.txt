SET(CMAKE_CXX_FLAGS -g)

# Make sure the compiler can find include files from our library
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib)

# Add executable called "helloDemo" that is built from the source files
ADD_EXECUTABLE(transvol transvol.cxx transform.cxx)
# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(transvol
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


ADD_EXECUTABLE(cropvol cropvol.cxx)
# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(cropvol
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

ADD_EXECUTABLE(colormap colormap.cxx)
# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(colormap
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

ADD_EXECUTABLE(convertslices convertslices.cxx)
# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(convertslices
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

ADD_EXECUTABLE(resample resample.cxx)
# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(resample
scan
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



