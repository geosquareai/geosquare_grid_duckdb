cmake_minimum_required(VERSION 3.5)

# Set extension name here
set(TARGET_NAME geosquare)

# Check if building for WASM
if(WASM_LOADABLE_EXTENSIONS)
  set(BUILDING_FOR_WASM TRUE)
  message(STATUS "Building for WASM environment")
endif()

# Check for ARM architecture
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm|aarch64)")
  set(BUILDING_FOR_ARM TRUE)
  message(STATUS "Building for ARM architecture: ${CMAKE_SYSTEM_PROCESSOR}")
  
  # ARM-specific compiler flags
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mcpu=native")
  endif()
endif()

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

# Handle GEOS dependency differently based on build target
if(BUILDING_FOR_WASM OR BUILDING_FOR_ARM)
  # For WASM and ARM builds, we need special handling
  set(GEOS_LIBRARIES "${CMAKE_BINARY_DIR}/deps/local/lib/libgeos_c.a")
  set(GEOS_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/deps/local/include")
  
  if(NOT EXISTS ${GEOS_LIBRARIES})
    # Try alternative locations for ARM
    if(BUILDING_FOR_ARM)
      set(GEOS_LIBRARIES "${CMAKE_BINARY_DIR}/deps/local/lib64/libgeos_c.a")
      if(NOT EXISTS ${GEOS_LIBRARIES})
        message(WARNING "GEOS library for ARM not found. Trying dynamic linking.")
        set(GEOS_LIBRARIES "geos_c")
      endif()
    else()
      message(WARNING "GEOS library for WASM not found at ${GEOS_LIBRARIES}. Linking may fail.")
    endif()
  endif()
  
  message(STATUS "Using specialized GEOS: ${GEOS_LIBRARIES}")
else()
  # Standard build - use find_package
  find_package(GEOS REQUIRED)
endif()

set(EXTENSION_NAME ${TARGET_NAME}_extension)
set(LOADABLE_EXTENSION_NAME ${TARGET_NAME}_loadable_extension)

project(${TARGET_NAME})
include_directories(src/include)

# Add GEOS include directory if we're in WASM or ARM mode
if(BUILDING_FOR_WASM OR BUILDING_FOR_ARM)
  include_directories(${GEOS_INCLUDE_DIRS})
endif()

set(EXTENSION_SOURCES src/geosquare_extension.cpp)

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})
build_loadable_extension(${TARGET_NAME} " " ${EXTENSION_SOURCES})

# Link GEOS libraries based on build target
if(BUILDING_FOR_WASM OR BUILDING_FOR_ARM)
  # For WASM and ARM, link directly to the static libraries or use system libs
  target_link_libraries(${EXTENSION_NAME} ${GEOS_LIBRARIES})
  target_link_libraries(${LOADABLE_EXTENSION_NAME} ${GEOS_LIBRARIES})
  
  # Add math library for ARM builds
  if(BUILDING_FOR_ARM)
    target_link_libraries(${EXTENSION_NAME} m)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} m)
  endif()
else()
  # Standard build uses the imported targets
  target_link_libraries(${EXTENSION_NAME} GEOS::geos_c)
  target_link_libraries(${LOADABLE_EXTENSION_NAME} GEOS::geos_c)
endif()

# Set C++ standard
set_property(TARGET ${EXTENSION_NAME} PROPERTY CXX_STANDARD 17)
set_property(TARGET ${LOADABLE_EXTENSION_NAME} PROPERTY CXX_STANDARD 17)

install(
  TARGETS ${EXTENSION_NAME}
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")
