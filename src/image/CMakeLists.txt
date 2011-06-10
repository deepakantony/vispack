#Then for each subdirectory listed in the SUBDIRS command,
# CMakeLists.txt files are created. In the ./Hello directory, the
# following CMakeLists.txt file is created:

# Create a library called "Hello" which includes the source file "hello.cxx".
# Any number of sources could be listed here.
ADD_LIBRARY(image image.cxx imagefile.cxx rgba.cxx imageRGBA.cxx
im.cxx indexlist.cxx vispacktiffio.cxx
)

INCLUDE_DIRECTORIES(${VISPACK_PACKAGES_DIR}/tiff ${VISPACK_PACKAGES_DIR}/fits ${VISPACK_PACKAGES_DIR}/jpeg)

