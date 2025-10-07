# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  if( atlas_HAVE_ECTRANS )
      
    add_fctest( TARGET atlas_fctest_trans
      LINKER_LANGUAGE Fortran
      SOURCES         fctest_trans.F90
      LIBS            atlas_f
      CONDITION       eckit_HAVE_MPI AND ( transi_HAVE_MPI OR ectrans_HAVE_MPI )
      MPI             4
      ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
    )
  endif()

  if( atlas_HAVE_ECTRANS )
    add_fctest( TARGET atlas_fctest_trans_invtrans_grad
      LINKER_LANGUAGE Fortran
      SOURCES         fctest_trans_invtrans_grad.F90
      LIBS            atlas_f
      ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
    )
  endif()

  add_fctest( TARGET atlas_fctest_trans_unstructured
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_trans_unstructured.F90
    LIBS            atlas_f
  )

endif()



ecbuild_add_test( TARGET atlas_test_trans_localcache
  SOURCES   test_trans_localcache.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_TRACE_REPORT=1
  CONDITION atlas_HAVE_ATLAS_TRANS
)


#--

# ...XX...
ecbuild_add_test( TARGET atlas_test_trans_mpi_setup_01
  SOURCES      test_trans_mpi_setup_01.cc
  LIBS         atlas transi
  ENVIRONMENT  ${ATLAS_TEST_ENVIRONMENT} ATLAS_FPE=FE_OVERFLOW,FE_DIVBYZERO
  CONDITION    atlas_HAVE_ECTRANS AND eckit_HAVE_MPI AND ( transi_HAVE_MPI OR ectrans_HAVE_MPI )
  MPI          4
)

ecbuild_add_test( TARGET atlas_test_trans_mpi_setup_02
  SOURCES      test_trans_mpi_setup_02.cc
  LIBS         atlas transi
  ENVIRONMENT  ${ATLAS_TEST_ENVIRONMENT} ATLAS_FPE=FE_OVERFLOW,FE_DIVBYZERO
  CONDITION    atlas_HAVE_ECTRANS AND eckit_HAVE_MPI AND ( transi_HAVE_MPI OR ectrans_HAVE_MPI )
  MPI          4
)
