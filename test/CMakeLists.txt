# Test suite provided by google test (GTest)
INCLUDE(FetchContent)
FETCHCONTENT_DECLARE(
        googletest
        URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
SET(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FETCHCONTENT_MAKEAVAILABLE(googletest)
ENABLE_TESTING()

FILE(GLOB_RECURSE TEST_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
ADD_EXECUTABLE(
        BLS_TEST
        ${TEST_FILES}
)

TARGET_LINK_LIBRARIES(
        BLS_TEST
        GTest::gtest_main
        BLS12_381
)

INCLUDE(GoogleTest)
GTEST_DISCOVER_TESTS(BLS_TEST)