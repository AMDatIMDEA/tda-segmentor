cmake_minimum_required(VERSION 3.3 FATAL_ERROR)


project(poretda)



#list(APPEND CMAKE_PREFIX_PATH "/Users/aditya/Programs/ParaView-v5.10.1/build/lib/cmake/paraview-5.10/vtk")
#list(APPEND CMAKE_PREFIX_PATH "/Users/aditya/Programs/ParaView-v5.10.1/build/lib/cmake/paraview-5.10/vtk/vtkm")
#list(APPEND CMAKE_PREFIX_PATH "/Users/aditya/Programs/ParaView-v5.10.1/build")
#list(APPEND CMAKE_PREFIX_PATH "/Users/aditya/Programs/ttk/ttk_install/lib/cmake/ttkVTK")







find_package(VTK COMPONENTS 
  vtkCommonColor
  vtkCommonCore
  vtkCommonDataModel
  vtkDomainsChemistry
  vtkDomainsChemistryOpenGL2
  vtkIOChemistry
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  QUIET
)
find_package(TTKVTK REQUIRED)

if (NOT VTK_FOUND)
  message("Skipping ReadCML: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(poretda MACOSX_BUNDLE main.cpp )
  target_link_libraries(ReadCML PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(poretda MACOSX_BUNDLE src/main.cpp src/instructions.cpp src/parameters.cpp src/logger.cpp)
  target_link_libraries(poretda PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS poretda
    MODULES ${VTK_LIBRARIES}
    )
endif ()

target_link_libraries(poretda
  PUBLIC
    ttk::vtk::ttkAll
    )
