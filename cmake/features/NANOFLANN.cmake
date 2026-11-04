### nanoflann k-d tree backend
#
# The nanoflann backend (atlas::util::detail::KDTree_nanoflann, using the bundled header-only library
# src/atlas/util/detail/nanoflann.hpp) is always available, and can be selected explicitly by constructing
# atlas::util::KDTree with `new atlas::util::detail::KDTree_nanoflann<Payload,Point>(...)`.
# This feature makes it the default backend of atlas::util::KDTree, instead of eckit::KDTreeMemory.

ecbuild_add_option( FEATURE NANOFLANN
                    DEFAULT OFF
                    DESCRIPTION "Use nanoflann instead of eckit as default k-d tree backend for atlas::util::KDTree" )
