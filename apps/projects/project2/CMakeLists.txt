SET(CMAKE_CXX_FLAGS -g)
# Make sure the compiler can find include files from the library.
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib)

# Link the executable to the libraries.
ADD_EXECUTABLE(project2 project2.cxx)
TARGET_LINK_LIBRARIES(project2 image jpeg tiff fitsio util)


