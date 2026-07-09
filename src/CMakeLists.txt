add_library(minisketch_field_sources INTERFACE)
target_sources(minisketch_field_sources
  INTERFACE
    fields/generic_1byte.cpp
    fields/generic_2bytes.cpp
    fields/generic_3bytes.cpp
    fields/generic_4bytes.cpp
    fields/generic_5bytes.cpp
    fields/generic_6bytes.cpp
    fields/generic_7bytes.cpp
    fields/generic_8bytes.cpp
)
if(HAVE_CLMUL)
  set(clmul_sources
    fields/clmul_1byte.cpp
    fields/clmul_2bytes.cpp
    fields/clmul_3bytes.cpp
    fields/clmul_4bytes.cpp
    fields/clmul_5bytes.cpp
    fields/clmul_6bytes.cpp
    fields/clmul_7bytes.cpp
    fields/clmul_8bytes.cpp
  )
  set_property(SOURCE minisketch.cpp
    PROPERTY COMPILE_DEFINITIONS HAVE_CLMUL
  )
  set_property(SOURCE ${clmul_sources}
    PROPERTY COMPILE_OPTIONS ${CLMUL_CXXFLAGS}
  )
  target_sources(minisketch_field_sources INTERFACE ${clmul_sources})
endif()

add_compile_definitions($<$<BOOL:${HAVE_CLZ}>:HAVE_CLZ>)

add_library(minisketch minisketch.cpp)
set_property(TARGET minisketch PROPERTY PUBLIC_HEADER
  ${PROJECT_SOURCE_DIR}/include/minisketch.h
)
target_include_directories(minisketch
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)
target_link_libraries(minisketch PRIVATE minisketch_field_sources)

add_library(minisketch_verify EXCLUDE_FROM_ALL minisketch.cpp)
target_compile_definitions(minisketch_verify
  PUBLIC
    MINISKETCH_VERIFY
)
target_link_libraries(minisketch_verify PRIVATE minisketch_field_sources)
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  target_compile_options(minisketch_verify
    PRIVATE
      /wd4702
  )
endif()

if(MINISKETCH_BUILD_TESTS)
  add_executable(test-noverify test.cpp)
  target_link_libraries(test-noverify PRIVATE minisketch)
  add_test(NAME ${PROJECT_NAME}_test_noverify COMMAND test-noverify)

  add_executable(test-verify test.cpp)
  target_link_libraries(test-verify PRIVATE minisketch_verify)
  add_test(NAME ${PROJECT_NAME}_test_verify COMMAND test-verify)

  # Unit tests for internal algorithm stages (Boost.Test). These compile the
  # field sources directly (not linked against libminisketch) so they can
  # exercise the sketch_impl.h templates; they are built with MINISKETCH_VERIFY
  # so internal CHECK_SAFE assertions fire.
  find_package(Boost CONFIG COMPONENTS unit_test_framework)
  if(TARGET Boost::unit_test_framework)
    add_executable(unit-tests unit_tests.cpp unit_tests_clmul.cpp)
    target_compile_definitions(unit-tests PRIVATE MINISKETCH_VERIFY BOOST_TEST_DYN_LINK)
    set_target_properties(unit-tests PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
    target_link_libraries(unit-tests PRIVATE Boost::unit_test_framework)
    if(HAVE_CLMUL)
      set_property(SOURCE unit_tests_clmul.cpp PROPERTY COMPILE_OPTIONS ${CLMUL_CXXFLAGS})
      set_property(SOURCE unit_tests_clmul.cpp APPEND PROPERTY COMPILE_DEFINITIONS HAVE_CLMUL)
    endif()
    add_test(NAME ${PROJECT_NAME}_unit_tests COMMAND unit-tests)
  else()
    message(WARNING "Boost.Test (unit_test_framework) not found; the unit-tests target will not be built.")
  endif()

  # Python test suites: pyminisketch's own unit tests, and a differential
  # test comparing the C library against the pure-Python reimplementation
  # through a shared library loaded with ctypes.
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND AND NOT WIN32)
    add_library(minisketch-testlib SHARED minisketch.cpp)
    set_target_properties(minisketch-testlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(minisketch-testlib PRIVATE minisketch_field_sources)

    add_test(NAME ${PROJECT_NAME}_pyminisketch
      COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/pyminisketch.py
    )
    add_test(NAME ${PROJECT_NAME}_differential
      COMMAND Python3::Interpreter ${PROJECT_SOURCE_DIR}/tests/differential.py
        --lib $<TARGET_FILE:minisketch-testlib>
    )
  endif()
endif()

if(MINISKETCH_BUILD_BENCHMARK)
  add_executable(bench bench.cpp)
  target_link_libraries(bench PRIVATE minisketch)
endif()

if(MINISKETCH_BUILD_FUZZ)
  # One fuzz binary containing all targets (select with the FUZZ env var),
  # linked against the verify library so CHECK_SAFE assertions fire.
  add_executable(fuzz
    fuzz/fuzz.cpp
    fuzz/decode.cpp
    fuzz/roundtrip.cpp
    fuzz/poly_ops.cpp
    fuzz/reconcile.cpp
  )
  target_compile_definitions(fuzz PRIVATE MINISKETCH_VERIFY)
  set_target_properties(fuzz PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
  target_link_libraries(fuzz PRIVATE minisketch_verify)
  target_link_options(fuzz PRIVATE -fsanitize=fuzzer)

  # Smoke tests: a short bounded run of each target so plain `ctest` exercises
  # them (real fuzzing sessions run the binary directly with a corpus dir).
  foreach(fuzz_target decode roundtrip poly_ops reconcile)
    add_test(NAME ${PROJECT_NAME}_fuzz_${fuzz_target}_smoke
      COMMAND fuzz -runs=2000 -max_len=512 -seed=1
    )
    set_tests_properties(${PROJECT_NAME}_fuzz_${fuzz_target}_smoke
      PROPERTIES ENVIRONMENT FUZZ=${fuzz_target}
    )
  endforeach()
endif()

if(MINISKETCH_INSTALL)
  include(GNUInstallDirs)
  install(TARGETS minisketch)
endif()
