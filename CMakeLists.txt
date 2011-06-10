#The first, top-level directory contains the following CMakeLists.txt file.
# The name of our project is "VISPACK".  CMakeLists files in this project can
# refer to the root source directory of the project as ${VISPACK_SOURCE_DIR} and
# to the root binary directory of the project as ${VISPACK_BINARY_DIR}.
PROJECT(VISPACK)

# Recurse into the "Hello" and "Demo" subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.
#SUBDIRS(src packages apps doc)
SUBDIRS(src packages apps)

#FIND_PACKAGE(TEEM REQUIRED)

#SET (NRRDIO_LIB /home/sci/whitaker/code/NrrdIO-1.9.0-linux)
#SET (NRRDIO_INCLUDE /home/sci/whitaker/code/NrrdIO-1.9.0-src)

SET (VISPACK_PACKAGES_DIR "${VISPACK_SOURCE_DIR}/packages")
SET(LIBRARY_OUTPUT_PATH ${VISPACK_BINARY_DIR}/lib)
SET(VISPACK_INCLUDE_DIRECTORIES ${VISPACK_SOURCE_DIR}/src ${VISPACK_SOURCE_DIR}/src/util ${VISPACK_SOURCE_DIR}/src/image ${VISPACK_SOURCE_DIR}/src/vol ${VISPACK_SOURCE_DIR}/src/mat ${VISPACK_SOURCE_DIR}/src/voxmodel ${VISPACK_SOURCE_DIR}/src/scan ${VISPACK_SOURCE_DIR}/src/param ${VISPACK_SOURCE_DIR}/src/diffuse
)

#Finally, in the ./Demo directory, the third and final CMakeLists.txt
#file is created: 
# Make sure the compiler can find include files from our Hello library.
#INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES} ${TEEM_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

#Need this for the sgi
#SET(CMAKE_CXX_FLAGS -LANG:std)
SET(CMAKE_CXX_FLAGS -g)
#SET(CMAKE_CXX_FLAGS -O2)





