cmake_minimum_required(VERSION 3.15)

project(LAB-MC658
        VERSION     1.0
        LANGUAGES   CXX)

set(CMAKE_CXX_STANDARD 11)

# Set a default build type if none was specified
set(default_build_type "Release")

# Headers directories
include_directories("include")
include_directories("lemon/lemon-1.3.1/include")
include_directories("brkga/src")

# Utils used in all targets
add_library(MY_LIB src/mygraphlib.cpp src/geompack.cpp src/myutils.cpp
        src/mycolor.cpp src/deprecated.cpp src/pickup_delivery_utils.cpp)

# Lemon lib directory
link_directories(${CMAKE_SOURCE_DIR}/lemon/lemon-1.3.1/lib)

# Compiler flags
set(CMAKE_CXX_FLAGS "-m64 -Ofast -Wall -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 -Wunused-local-typedefs -lemon")

# Target executable and linking
add_executable(lab1_pickup_delivery src/lab1_pickup_delivery.cpp)
target_link_libraries(lab1_pickup_delivery m)
target_link_libraries(lab1_pickup_delivery pthread)
target_link_libraries(lab1_pickup_delivery MY_LIB)

add_executable(lab2_pickup_delivery src/lab2_pickup_delivery.cpp src/PickupDeliveryDecoder.cpp)
target_link_libraries(lab2_pickup_delivery m)
target_link_libraries(lab2_pickup_delivery pthread)
target_link_libraries(lab2_pickup_delivery MY_LIB)
