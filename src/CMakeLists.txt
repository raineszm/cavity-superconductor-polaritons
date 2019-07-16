if(CMAKE_BUILD_TYPE STREQUAL "Release")
pybind11_add_module(_bardasis_schrieffer
    bindings.cpp
    bindsystem.cpp
    bindstate.cpp
    bindcoupling.cpp
    bindbs.cpp
    bindpolariton.cpp
    bindcavity.cpp)

  target_link_libraries(_bardasis_schrieffer PUBLIC Eigen3::Eigen  GSL::gsl)
  target_compile_features(_bardasis_schrieffer PUBLIC cxx_std_17)
  target_include_directories(_bardasis_schrieffer PUBLIC ${CMAKE_SOURCE_DIR}/include)
  install(TARGETS _bardasis_schrieffer DESTINATION ${CMAKE_SOURCE_DIR}/python/lib/bardasis_schrieffer)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_executable(debug debug.cpp)
  target_link_libraries(debug PUBLIC Eigen3::Eigen  GSL::gsl)
  target_compile_features(debug PUBLIC cxx_std_17)
  target_include_directories(debug PUBLIC ${CMAKE_SOURCE_DIR}/include)
endif()
