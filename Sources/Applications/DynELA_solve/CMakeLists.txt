FILE(GLOB SRCS *.h *.C)

# Projet global : Application DynELA solve
project(DynELA_solve)

# définitions générales

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} 
	${SOURCE_LIBRARY_DIRECTORY}/basicTools 
	${SOURCE_LIBRARY_DIRECTORY}/lapack 
	${SOURCE_LIBRARY_DIRECTORY}/linearAlgebra 
	${SOURCE_LIBRARY_DIRECTORY}/finiteElement
	${SOURCE_LIBRARY_DIRECTORY}/interpretor
  ${SOURCE_LIBRARY_DIRECTORY}/lsdyna)
  
LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH})

#ADD_DEFINITIONS(-Wall -O2)

# Création d'un exécutable
# ADD_EXECUTABLE(DynELA_solve ${SRCS})
# TARGET_LINK_LIBRARIES(DynELA_solve dnlFE dnlLA dnlBT lapack )
# #INSTALL(TARGETS DynELA_solve DESTINATION ${DynELA_base}/bin)
link_directories("${MMG_DIR}/lib")
include_directories("${MMG_DIR}/include")

ADD_EXECUTABLE(test element_2d.cpp)
TARGET_LINK_LIBRARIES(test dnlFE mmg dnlLA dnlBT lapack lsdyna)


ADD_EXECUTABLE(test_contact element_2d_contact.cpp)
TARGET_LINK_LIBRARIES(test_contact dnlFE mmg dnlLA dnlBT lapack lsdyna)


#SINGLE 1 time remesh
#ADD_EXECUTABLE(test_contact_1time_rmsh test_contact_1time_rmsh.cpp)
#TARGET_LINK_LIBRARIES(test_contact_1time_rmsh dnlFE mmg dnlLA dnlBT lapack lsdyna)

ADD_EXECUTABLE(test_contact_solver_amr test_contact_solver_amr.cpp)
TARGET_LINK_LIBRARIES(test_contact_solver_amr dnlFE mmg dnlLA dnlBT lapack lsdyna)

ADD_EXECUTABLE(test_tri test_tri.cpp)
TARGET_LINK_LIBRARIES(test_tri dnlFE mmg dnlLA dnlBT lapack lsdyna)


