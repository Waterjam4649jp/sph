cmake_minimum_required(VERSION 3.10)
project(sph)

set(CMAKE_CXX_STANDARD 17)

find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)
add_executable(
    main
    src/main.cpp
    src/sph.cpp
    src/sph.h)
target_link_libraries(main sfml-graphics sfml-window sfml-system)

