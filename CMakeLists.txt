cmake_minimum_required(VERSION 3.16)

project(simDER LANGUAGES CXX)

option(SIMDER_BUILD_PYTHON "Build the optional pybind11 planning module" OFF)

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
	tipMagneticForce.cpp
    timeStepper.cpp
    setInput.cpp
    externalContactForce.cpp
	geometry/triangleGeometry.cpp
	geometry/triangleMeshSurface.cpp
	geometry/planarSlabDomain.cpp
	geometry/sphericalShellDomain.cpp
	geometry/sphericalObstacleDomain.cpp
	geometry/doubleSphericalObstacleDomain.cpp
	contact/contactDetector.cpp
	contact/barrierPotential.cpp
	contact/planarBarrierContactForce.cpp
	contact/contactKktSystem.cpp
	contact/contactConstraintAnalysis.cpp
	insertion/proximalGuideInsertionModel.cpp
	planning/mechanicsSession.cpp
)

target_include_directories(simder_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(simder_core PUBLIC Eigen3::Eigen ${LAPACK_LIBRARIES})
set_target_properties(simder_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(simDER main.cpp)
target_link_libraries(simDER PRIVATE simder_core OpenGL::GL OpenGL::GLU GLUT::GLUT)

include(CTest)

if(BUILD_TESTING)
	add_executable(static_equilibrium_test tests/static_equilibrium_test.cpp)
	target_link_libraries(static_equilibrium_test PRIVATE simder_core)
	add_executable(rod_state_test tests/rod_state_test.cpp)
	target_link_libraries(rod_state_test PRIVATE simder_core)
	add_executable(static_solver_test tests/static_solver_test.cpp)
	target_link_libraries(static_solver_test PRIVATE simder_core)
	add_executable(elastic_jacobian_test tests/elastic_jacobian_test.cpp)
	target_link_libraries(elastic_jacobian_test PRIVATE simder_core)
	add_executable(relaxation_comparison_test tests/relaxation_comparison_test.cpp)
	target_link_libraries(relaxation_comparison_test PRIVATE simder_core)
	add_executable(tip_magnetic_test tests/tip_magnetic_test.cpp)
	target_link_libraries(tip_magnetic_test PRIVATE simder_core)
	add_executable(field_sensitivity_test tests/field_sensitivity_test.cpp)
	target_link_libraries(field_sensitivity_test PRIVATE simder_core)
	add_executable(triangle_geometry_test tests/triangle_geometry_test.cpp)
	target_link_libraries(triangle_geometry_test PRIVATE simder_core)
	add_executable(planar_slab_test tests/planar_slab_test.cpp)
	target_link_libraries(planar_slab_test PRIVATE simder_core)
	add_executable(spherical_shell_test tests/spherical_shell_test.cpp)
	target_link_libraries(spherical_shell_test PRIVATE simder_core)
	add_executable(spherical_obstacle_test tests/spherical_obstacle_test.cpp)
	target_link_libraries(spherical_obstacle_test PRIVATE simder_core)
	add_executable(double_spherical_obstacle_test
		tests/double_spherical_obstacle_test.cpp)
	target_link_libraries(double_spherical_obstacle_test PRIVATE simder_core)
	add_executable(contact_detector_test tests/contact_detector_test.cpp)
	target_link_libraries(contact_detector_test PRIVATE simder_core)
	add_executable(barrier_contact_test tests/barrier_contact_test.cpp)
	target_link_libraries(barrier_contact_test PRIVATE simder_core)
	add_executable(contact_equilibrium_test tests/contact_equilibrium_test.cpp)
	target_link_libraries(contact_equilibrium_test PRIVATE simder_core)
	add_executable(contact_kkt_test tests/contact_kkt_test.cpp)
	target_link_libraries(contact_kkt_test PRIVATE simder_core)
	add_executable(contact_kkt_equilibrium_test tests/contact_kkt_equilibrium_test.cpp)
	target_link_libraries(contact_kkt_equilibrium_test PRIVATE simder_core)
	add_executable(contact_constraint_analysis_test tests/contact_constraint_analysis_test.cpp)
	target_link_libraries(contact_constraint_analysis_test PRIVATE simder_core)
	add_executable(contact_field_sensitivity_test tests/contact_field_sensitivity_test.cpp)
	target_link_libraries(contact_field_sensitivity_test PRIVATE simder_core)
	add_executable(spherical_contact_kkt_test tests/spherical_contact_kkt_test.cpp)
	target_link_libraries(spherical_contact_kkt_test PRIVATE simder_core)
	add_executable(field_continuation_test tests/field_continuation_test.cpp)
	target_link_libraries(field_continuation_test PRIVATE simder_core)
	add_executable(insertion_model_test tests/insertion_model_test.cpp)
	target_link_libraries(insertion_model_test PRIVATE simder_core)
	add_executable(actuation_continuation_test tests/actuation_continuation_test.cpp)
	target_link_libraries(actuation_continuation_test PRIVATE simder_core)
	add_executable(mechanics_session_test tests/mechanics_session_test.cpp)
	target_link_libraries(mechanics_session_test PRIVATE simder_core)

    add_test(
        NAME legacy_headless_smoke
        COMMAND simDER ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke_options.txt
    )
	add_test(NAME zero_load_static_equilibrium COMMAND static_equilibrium_test)
	add_test(
		NAME rod_state_round_trip
		COMMAND rod_state_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke_options.txt
	)
	add_test(
		NAME contact_free_static_newton
		COMMAND static_solver_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/static_failure_options.txt
	)
	add_test(NAME elastic_jacobian_finite_difference COMMAND elastic_jacobian_test)
	add_test(
		NAME static_dynamic_relaxation_comparison
		COMMAND relaxation_comparison_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/relaxation_options.txt
	)
	add_test(
		NAME axial_tip_magnetic_validation
		COMMAND tip_magnetic_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/tip_magnetic_options.txt
	)
	add_test(
		NAME equilibrium_field_sensitivity
		COMMAND field_sensitivity_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/tip_magnetic_options.txt
	)
	add_test(
		NAME triangle_surface_projection
		COMMAND triangle_geometry_test
			${CMAKE_CURRENT_SOURCE_DIR}/mesh_1/nodesInput.txt
			${CMAKE_CURRENT_SOURCE_DIR}/mesh_1/triangleInput.txt
	)
	add_test(NAME planar_slab_clearance COMMAND planar_slab_test)
	add_test(NAME spherical_shell_clearance COMMAND spherical_shell_test)
	add_test(NAME spherical_obstacle_clearance COMMAND spherical_obstacle_test)
	add_test(NAME double_spherical_obstacle_clearance
		COMMAND double_spherical_obstacle_test)
	add_test(NAME contact_candidate_detection COMMAND contact_detector_test)
	add_test(NAME planar_barrier_contact_validation COMMAND barrier_contact_test)
	add_test(NAME loaded_planar_contact_equilibrium COMMAND contact_equilibrium_test)
	add_test(NAME planar_contact_kkt_linearization COMMAND contact_kkt_test)
	add_test(NAME nonlinear_planar_contact_kkt COMMAND contact_kkt_equilibrium_test)
	add_test(NAME contact_constraint_analysis COMMAND contact_constraint_analysis_test)
	add_test(NAME contact_field_sensitivity COMMAND contact_field_sensitivity_test)
	add_test(NAME spherical_contact_kkt_mechanics COMMAND spherical_contact_kkt_test)
	add_test(NAME adaptive_field_continuation COMMAND field_continuation_test)
	add_test(NAME proximal_guide_insertion_derivatives COMMAND insertion_model_test)
	add_test(NAME combined_actuation_continuation COMMAND actuation_continuation_test)
	add_test(NAME mechanics_session_facade COMMAND mechanics_session_test)

    set_tests_properties(
        legacy_headless_smoke
        zero_load_static_equilibrium
		rod_state_round_trip
		contact_free_static_newton
		elastic_jacobian_finite_difference
		static_dynamic_relaxation_comparison
		axial_tip_magnetic_validation
		equilibrium_field_sensitivity
		triangle_surface_projection
		planar_slab_clearance
		spherical_shell_clearance
		spherical_obstacle_clearance
		double_spherical_obstacle_clearance
		contact_candidate_detection
		planar_barrier_contact_validation
		loaded_planar_contact_equilibrium
		planar_contact_kkt_linearization
		nonlinear_planar_contact_kkt
		contact_constraint_analysis
		contact_field_sensitivity
		spherical_contact_kkt_mechanics
		adaptive_field_continuation
		proximal_guide_insertion_derivatives
		combined_actuation_continuation
		mechanics_session_facade
        PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

if(SIMDER_BUILD_PYTHON)
	find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
	find_package(pybind11 CONFIG QUIET)
	if(NOT pybind11_FOUND)
		execute_process(
			COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
			OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE
			RESULT_VARIABLE PYBIND11_CMAKE_RESULT
		)
		if(PYBIND11_CMAKE_RESULT EQUAL 0)
			find_package(pybind11 CONFIG REQUIRED PATHS ${PYBIND11_CMAKE_DIR})
		else()
			message(FATAL_ERROR
				"SIMDER_BUILD_PYTHON=ON requires pybind11 for ${Python3_EXECUTABLE}")
		endif()
	endif()

	pybind11_add_module(simder MODULE python/bindings.cpp)
	target_link_libraries(simder PRIVATE simder_core)

	if(BUILD_TESTING)
		add_test(
			NAME python_mechanics_session_smoke
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_bindings_smoke.py
		)
		add_test(
			NAME python_baseline_rrt
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_rrt_test.py
		)
		add_test(
			NAME python_terminal_connector
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_terminal_connector_test.py
		)
		add_test(
			NAME python_task_a_experiment
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_task_a_experiment_test.py
		)
		add_test(
			NAME python_task_a_comparison
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_task_a_comparison_test.py
		)
		add_test(
			NAME python_planar_contact_task
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_planar_contact_task_test.py
		)
		add_test(
			NAME python_spherical_shell_task
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_spherical_shell_task_test.py
		)
		add_test(
			NAME python_dead_end_task
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_dead_end_task_test.py
		)
		add_test(
			NAME python_dead_end_sweep
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_dead_end_sweep_test.py
		)
		add_test(
			NAME python_contact_branching_mechanics
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_contact_branching_mechanics_test.py
		)
		add_test(
			NAME python_contact_branching_budget
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/tests/python_contact_branching_budget_test.py
		)
		set_tests_properties(
			python_mechanics_session_smoke
			python_baseline_rrt
			python_terminal_connector
			python_task_a_experiment
			python_task_a_comparison
			python_planar_contact_task
			python_spherical_shell_task
			python_dead_end_task
			python_dead_end_sweep
			python_contact_branching_mechanics
			python_contact_branching_budget
			PROPERTIES
				WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
				ENVIRONMENT
					"PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:$<TARGET_FILE_DIR:simder>"
		)
	endif()
endif()
