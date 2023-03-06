cmake_minimum_required(VERSION 3.22)
project(FCPMS)

set(CMAKE_CXX_STANDARD 14)

#find_package(Qt5Core REQUIRED)
#find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)


add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/qcustomplot)
add_executable(FCPMS
        gui.cpp
        gui.ui
        main.cpp GenerationSignal.cpp)
add_dependencies(FCPMS qcustomplot)



target_compile_definitions(FCPMS PRIVATE QCUSTOMPLOT_USE_LIBRARY)
target_link_libraries(FCPMS PRIVATE Qt5::Widgets qcustomplot)
