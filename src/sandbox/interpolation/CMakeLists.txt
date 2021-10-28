# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

list (APPEND Data_yaml
  interp.yaml
)
foreach(FILENAME ${Data_yaml})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

ecbuild_add_executable(
    TARGET  atlas-parallel-interpolation-fieldset
    SOURCES atlas-parallel-interpolation-fieldsets.cc
            InterpRedistrWrapper.h
            InterpRedistrWrapper.cc
            InterpRedistrWrapper2.h
            InterpRedistrWrapper2.cc
            InterpRedistrWrapper3.h
            InterpRedistrWrapper3.cc
    LIBS    atlas
    NOINSTALL
)

ecbuild_add_executable(
    TARGET  atlas-parallel-interpolation
    SOURCES atlas-parallel-interpolation.cc
            PartitionedMesh.h
            PartitionedMesh.cc
    LIBS    atlas
    NOINSTALL
)

ecbuild_add_executable(
    TARGET  atlas-parallel-structured-interpolation
    SOURCES atlas-parallel-structured-interpolation.cc
    LIBS    atlas
    NOINSTALL
)

