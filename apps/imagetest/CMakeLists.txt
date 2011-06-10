SET(CMAKE_CXX_FLAGS -g)
# Make sure the compiler can find include files from the library.
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib)

# Link the executable to the libraries.
ADD_EXECUTABLE(imagetest imagetest.cxx)
TARGET_LINK_LIBRARIES(imagetest image jpeg tiff fitsio util)

ADD_EXECUTABLE(read_text read_text.cxx)
TARGET_LINK_LIBRARIES(read_text image jpeg tiff fitsio util)


