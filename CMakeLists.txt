cmake_minimum_required(VERSION 3.5)

# Set extension name here
set(TARGET_NAME geosquare)

if(NOT EXISTS ${CMAKE_BINARY_DIR}/deps)
  execute_process(
    # Generate project
    COMMAND
      ${CMAKE_COMMAND} -G ${CMAKE_GENERATOR}
      -DDUCKDB_ENABLE_DEPRECATED_API=1
      -DWASM_LOADABLE_EXTENSIONS=1
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DOSX_BUILD_ARCH=${OSX_BUILD_ARCH}
      -DSPATIAL_USE_NETWORK=${SPATIAL_USE_NETWORK}
      -DOPENSSL_ROOT_DIR=$ENV{OPENSSL_ROOT_DIR}
      -DVCPKG_MANIFEST_DIR='${VCPKG_MANIFEST_DIR}'
      -DVCPKG_INSTALLED_DIR='${CMAKE_BINARY_DIR}/vcpkg_installed'
      -DCMAKE_TOOLCHAIN_FILE='${CMAKE_TOOLCHAIN_FILE}'
      -DVCPKG_TARGET_TRIPLET='${VCPKG_TARGET_TRIPLET}'
      -S ${CMAKE_CURRENT_SOURCE_DIR}/deps -B ${CMAKE_BINARY_DIR}/deps
    RESULT_VARIABLE DEPENDENCIES_GENERATE_RESULT)
  if(NOT DEPENDENCIES_GENERATE_RESULT EQUAL 0)
    message(FATAL_ERROR "Could not generate dependencies project")
  endif()

  execute_process(
    # Build project
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/deps --config
            ${CMAKE_BUILD_TYPE} RESULT_VARIABLE DEPENDENCIES_BUILD_RESULT)
  if(NOT DEPENDENCIES_BUILD_RESULT EQUAL 0)
    message(FATAL_ERROR "Could not build dependencies project")
  endif()
endif()

# Add dependencies
list(APPEND CMAKE_FIND_ROOT_PATH "${CMAKE_BINARY_DIR}/deps/local/")
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}/deps/local/")

# DuckDB's extension distribution supports vcpkg. As such, dependencies can be added in ./vcpkg.json and then
# used in cmake with find_package. Feel free to remove or replace with other dependencies.
# Note that it should also be removed from vcpkg.json to prevent needlessly installing it..
#find_package(OpenSSL REQUIRED)
find_package(GEOS REQUIRED)

set(EXTENSION_NAME ${TARGET_NAME}_extension)
set(LOADABLE_EXTENSION_NAME ${TARGET_NAME}_loadable_extension)

project(${TARGET_NAME})
include_directories(src/include)

set(EXTENSION_SOURCES src/geosquare_extension.cpp)

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})
build_loadable_extension(${TARGET_NAME} " " ${EXTENSION_SOURCES})

# Link OpenSSL in both the static library as the loadable extension
#target_link_libraries(${EXTENSION_NAME} OpenSSL::SSL OpenSSL::Crypto)
#target_link_libraries(${LOADABLE_EXTENSION_NAME} OpenSSL::SSL OpenSSL::Crypto)
target_link_libraries(${EXTENSION_NAME} GEOS::geos_c)
target_link_libraries(${LOADABLE_EXTENSION_NAME} GEOS::geos_c)

install(
  TARGETS ${EXTENSION_NAME}
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")
