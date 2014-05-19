
cmake_minimum_required(VERSION 2.8)

project(libvio)

FILE(GLOB SOURCE_FILES src/codec/h264/**/*.cpp)

add_executable(
  libvio
  ${SOURCE_FILES}
)

set (CMAKE_INSTALL_PREFIX ..)
set (CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++")
add_definitions(-std=c++11 -Wno-deprecated-declarations)

include_directories(src/codec/h264/core src/codec/h264/decoder src/codec/h264/framebuf src/codec/h264/parser)

target_link_libraries(libvio)

# add the intstall targets

install(TARGETS libvio DESTINATION bin)
