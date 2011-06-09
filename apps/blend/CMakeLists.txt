
# Make sure the compiler can find include files from our Hellolibrary.
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib)

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx".  
ADD_EXECUTABLE(blend blend.C morph.C)

# Link the executable to the Hello library.
TARGET_LINK_LIBRARIES(blend matrix lapack1 lapack2 blas m F77 image jpeg
tiff fitsio util m)
