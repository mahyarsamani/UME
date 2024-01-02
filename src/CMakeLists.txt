# Ume/src/CMakeLists.txt

add_subdirectory(Ume)

add_executable(ume_serial ume_serial.cc)
target_include_directories(ume_serial PRIVATE ${ROI_BASE_INC} ${ROI_TYPE_INC})
target_link_directories(ume_serial PRIVATE ${ROI_BASE_LIB_DIR} ${ROI_TYPE_LIB_DIR})
target_link_libraries(ume_serial Ume ${ROI_LIB_OBJS})
set_target_properties(ume_serial PROPERTIES
  COMPILE_DEFINITIONS "${COMMON_COMPILE_DEFINITIONS};${ROI_CXX_FLAGS}"
  COMPILE_OPTIONS "${WARNING_FLAGS}"
  )

add_executable(txt2bin txt2bin.cc)
target_link_libraries(txt2bin Ume)
target_compile_options(txt2bin PRIVATE ${WARNING_FLAGS}
  ${COMMON_COMPILE_OPTIONS})


#############################################################################
# Things that depend on MPI
#############################################################################
if(USE_MPI)
add_executable(ume_mpi ume_mpi.cc)
target_include_directories(ume_mpi PRIVATE ${ROI_BASE_INC} ${ROI_TYPE_INC})
target_link_directories(ume_mpi PRIVATE ${ROI_BASE_LIB_DIR} ${ROI_TYPE_LIB_DIR})
target_link_libraries(ume_mpi Ume ${ROI_LIB_OBJS})
set_target_properties(ume_mpi PROPERTIES
  COMPILE_DEFINITIONS "${COMMON_COMPILE_DEFINITIONS};${ROI_CXX_FLAGS}"
  COMPILE_OPTIONS "${WARNING_FLAGS}"
  )
endif()


#############################################################################
# Things that depend on Kokkos
#############################################################################

if(USE_KOKKOS)
  add_executable(kokkos_hello kokkos_hello.cc)
  target_link_libraries(kokkos_hello
    Ume
    Kokkos::kokkos
    )
endif()


#############################################################################
# Things that depend on OpenACC
#############################################################################

if(USE_OPENACC)
  add_executable(openacc_hello openacc_hello.cc)
  target_link_libraries(openacc_hello
    Ume
    OpenACC::OpenACC_CXX
    )
endif()
