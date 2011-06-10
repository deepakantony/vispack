SET(CMAKE_CXX_FLAGS -g)
# Make sure the compiler can find include files from the library.
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib)

# Link the executable to the libraries.
ADD_EXECUTABLE(edgedetector edgedetector.cxx util.cxx edgemain.cxx)

TARGET_LINK_LIBRARIES(edgedetector matrix param image lapack lapack1 lapack2 jpeg tiff fitsio util)

ADD_EXECUTABLE(htwithradius edgedetector.cxx util.cxx houghtransform.cxx htwithradius.cxx)

TARGET_LINK_LIBRARIES(htwithradius matrix param image lapack lapack1 lapack2 jpeg tiff fitsio util)

ADD_EXECUTABLE(htnoradius edgedetector.cxx util.cxx houghtransform.cxx htmanyradius.cxx)

TARGET_LINK_LIBRARIES(htnoradius matrix param image lapack lapack1 lapack2 jpeg tiff fitsio util)

ADD_EXECUTABLE(anisodiff util.cxx anisodiff.cxx)

TARGET_LINK_LIBRARIES(anisodiff matrix param image lapack lapack1 lapack2 jpeg tiff fitsio util)
