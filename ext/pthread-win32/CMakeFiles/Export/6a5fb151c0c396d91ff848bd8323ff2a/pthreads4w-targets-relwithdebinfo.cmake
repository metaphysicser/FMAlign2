#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "pthreads4w::pthreadVCE3" for configuration "RelWithDebInfo"
set_property(TARGET pthreads4w::pthreadVCE3 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(pthreads4w::pthreadVCE3 PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/pthreadVCE3.lib"
  )

list(APPEND _cmake_import_check_targets pthreads4w::pthreadVCE3 )
list(APPEND _cmake_import_check_files_for_pthreads4w::pthreadVCE3 "${_IMPORT_PREFIX}/lib/pthreadVCE3.lib" )

# Import target "pthreads4w::pthreadVSE3" for configuration "RelWithDebInfo"
set_property(TARGET pthreads4w::pthreadVSE3 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(pthreads4w::pthreadVSE3 PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/pthreadVSE3.lib"
  )

list(APPEND _cmake_import_check_targets pthreads4w::pthreadVSE3 )
list(APPEND _cmake_import_check_files_for_pthreads4w::pthreadVSE3 "${_IMPORT_PREFIX}/lib/pthreadVSE3.lib" )

# Import target "pthreads4w::pthreadVC3" for configuration "RelWithDebInfo"
set_property(TARGET pthreads4w::pthreadVC3 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(pthreads4w::pthreadVC3 PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "C"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/pthreadVC3.lib"
  )

list(APPEND _cmake_import_check_targets pthreads4w::pthreadVC3 )
list(APPEND _cmake_import_check_files_for_pthreads4w::pthreadVC3 "${_IMPORT_PREFIX}/lib/pthreadVC3.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
