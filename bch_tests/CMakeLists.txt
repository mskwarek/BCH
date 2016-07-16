add_executable(runBasicTests
        runBasicTests.cpp test_cases.h)

target_link_libraries(runBasicTests BCH_codes)
target_link_libraries(runBasicTests gtest gtest_main)
