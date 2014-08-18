CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

PROJECT(mbmesh)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wno-shadow -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wundef -Wno-unused-parameter -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option")
ENDIF()

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

ADD_LIBRARY(mbmesh mbmesh.cpp)

ADD_EXECUTABLE(mbmesh2text mbmesh2text.cpp)
TARGET_LINK_LIBRARIES(mbmesh2text mbmesh)

ADD_EXECUTABLE(mbanim2text mbanim2text.cpp)
TARGET_LINK_LIBRARIES(mbanim2text mbmesh)

FIND_PACKAGE(PhysFS REQUIRED)
INCLUDE_DIRECTORIES(${PHYSFS_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(mbmesh ${PHYSFS_LIBRARY})

