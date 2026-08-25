cmake_minimum_required(VERSION 3.16)

project(simDER LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(LAPACK REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

add_library(simder_core
    world.cpp
    elasticRod.cpp
    elasticStretchingForce.cpp
    elasticBendingForce.cpp
    elasticTwistingForce.cpp
    externalGravityForce.cpp
    inertialForce.cpp
    dampingForce.cpp
    externalMagneticForce.cpp
    timeStepper.cpp
    setInput.cpp
    externalContactForce.cpp
)

target_include_directories(simder_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(simder_core PUBLIC Eigen3::Eigen ${LAPACK_LIBRARIES})

add_executable(simDER main.cpp)
target_link_libraries(simDER PRIVATE simder_core OpenGL::GL OpenGL::GLU GLUT::GLUT)

include(CTest)

if(BUILD_TESTING)
	add_executable(static_equilibrium_test tests/static_equilibrium_test.cpp)
	target_link_libraries(static_equilibrium_test PRIVATE simder_core)

    add_test(
        NAME legacy_headless_smoke
        COMMAND simDER ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke_options.txt
    )
	add_test(NAME zero_load_static_equilibrium COMMAND static_equilibrium_test)

    set_tests_properties(
        legacy_headless_smoke
        zero_load_static_equilibrium
        PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()
