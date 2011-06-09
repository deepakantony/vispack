SUBDIRS(backup)

# Make sure the compiler can find include files from our library
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

# Make sure the linker can find the library once it is built.
LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib)

# Add executables
ADD_EXECUTABLE(MIregister MIregister.cxx registrator.h)
ADD_EXECUTABLE(transformImages transformImages.cxx)
ADD_EXECUTABLE(viewer viewer.cxx)
ADD_EXECUTABLE(grapher grapher.cxx)
ADD_EXECUTABLE(graphCorrelation graphCorrelation.cxx)
ADD_EXECUTABLE(graphExperimental graphExperimental.cxx)
ADD_EXECUTABLE(graphMI graphMI.cxx)
ADD_EXECUTABLE(makeRawFloat makeRawFloat.cxx)
ADD_EXECUTABLE(converter converter.cxx)


# Link the executables to the library.
TARGET_LINK_LIBRARIES(MIregister
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
TARGET_LINK_LIBRARIES(transformImages
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
TARGET_LINK_LIBRARIES(viewer
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
TARGET_LINK_LIBRARIES(grapher
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
TARGET_LINK_LIBRARIES(graphCorrelation
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
TARGET_LINK_LIBRARIES(graphExperimental
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
TARGET_LINK_LIBRARIES(graphMI
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
TARGET_LINK_LIBRARIES(makeRawFloat
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
TARGET_LINK_LIBRARIES(converter
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