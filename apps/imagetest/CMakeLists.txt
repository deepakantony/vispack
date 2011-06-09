
#SET(CMAKE_CXX_FLAGS -g)

# Make sure the compiler can find include files from our Hellolibrary.
INCLUDE_DIRECTORIES(. ${VISPACK_INCLUDE_DIRECTORIES})

# Make sure the linker can find the Hello library once it is built.
LINK_DIRECTORIES(${VISPACK_BINARY_DIR}/lib)

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx".  
#ADD_EXECUTABLE(imagetest imagetest.cxx)
#ADD_EXECUTABLE(paramtest paramtest.cxx)

ADD_EXECUTABLE(conor conor.cxx)
TARGET_LINK_LIBRARIES(conor image jpeg tiff fitsio util matrix lapack1 lapack2 blas F77 m)

ADD_EXECUTABLE(k_nn k_nn.cxx)
TARGET_LINK_LIBRARIES(k_nn image jpeg tiff fitsio util matrix lapack1 lapack2 blas F77 m)

ADD_EXECUTABLE(k_nn_test k_nn_test.cxx)
TARGET_LINK_LIBRARIES(k_nn_test image jpeg tiff fitsio util matrix lapack1 lapack2 blas F77 m)

ADD_EXECUTABLE(totalvariation totalvariation.cxx)
TARGET_LINK_LIBRARIES(totalvariation image jpeg tiff fitsio util m)

ADD_EXECUTABLE(imagedifference imagedifference.cxx)
TARGET_LINK_LIBRARIES(imagedifference image jpeg tiff fitsio util m)

ADD_EXECUTABLE(divflow divflow.cxx)
TARGET_LINK_LIBRARIES(divflow image jpeg tiff fitsio util m)

ADD_EXECUTABLE(modflow modflow.cxx mtxlib.cxx)
TARGET_LINK_LIBRARIES(modflow image jpeg tiff fitsio util m)

ADD_EXECUTABLE(modflow_old modflow_old.cxx)
TARGET_LINK_LIBRARIES(modflow_old image jpeg tiff fitsio util m)

ADD_EXECUTABLE(noise_tmp noise_tmp.cxx)
TARGET_LINK_LIBRARIES(noise_tmp image jpeg tiff fitsio util m)

ADD_EXECUTABLE(read_text read_text.cxx)
TARGET_LINK_LIBRARIES(read_text image jpeg tiff fitsio util m)

# Link the executable to the Hello library.
#TARGET_LINK_LIBRARIES(imagetest image jpeg tiff fitsio util m)
#TARGET_LINK_LIBRARIES(paramtest param)

ADD_EXECUTABLE(medianmode medianmode.cxx)
TARGET_LINK_LIBRARIES(medianmode param util matrix lapack1 lapack2 blas F77 m)

ADD_EXECUTABLE(process_tori process_tori.cxx)
TARGET_LINK_LIBRARIES(process_tori param util matrix lapack1 lapack2 blas F77 m)

ADD_EXECUTABLE(setpoint setpoint.cxx)
TARGET_LINK_LIBRARIES(setpoint image jpeg tiff fitsio util m)

ADD_EXECUTABLE(removeyellow removeyellow.cxx)
TARGET_LINK_LIBRARIES(removeyellow image jpeg tiff fitsio util m)

ADD_EXECUTABLE(overlaycolor overlaycolor.cxx)
TARGET_LINK_LIBRARIES(overlaycolor image jpeg tiff fitsio util m)

 ADD_EXECUTABLE(fourier fourier.cxx fourn.cxx)
 TARGET_LINK_LIBRARIES(fourier
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

