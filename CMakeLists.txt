cmake_minimum_required(VERSION 3.14)

# Project name and version
project(SimpleMPI VERSION 1.0
	DESCRIPTION "Extended Application of Huckel's Theory"
		  LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Add the executable
add_executable(extended_huckel src/main.cpp)

# Find Armadillo and json
find_package(Armadillo REQUIRED)
find_package(nlohmann_json 3.9.1 REQUIRED)
#target_link_libraries(extended_huckel PUBLIC Armadillo::Armadillo)
target_link_libraries(extended_huckel 
PUBLIC 
${ARMADILLO_LIBRARIES}
nlohmann_json::nlohmann_json
)