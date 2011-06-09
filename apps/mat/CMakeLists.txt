
# Make sure the compiler can find include files from our Hellolibrary.
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib)

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx".  
ADD_EXECUTABLE(mattest mattest.cxx)

# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(mattest util matrix lapack1 lapack2 blas F77 m)
