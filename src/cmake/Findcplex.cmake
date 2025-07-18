SET(CPLEX_ROOT_DIR "" CACHE PATH "CPLEX root directory")
if(NOT CPLEX_ROOT_DIR)
  set(CPLEX_ROOT_DIR $ENV{CPLEX_ROOT_DIR})
endif()

find_path(CPLEX_INCLUDE_DIR
  ilcplex/cplex.h
  PATHS "C:/ILOG/CPLEX/include"
  PATHS "/opt/ilog/cplex/include"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/cplex/include"
  PATHS "~/ILOG/CPLEX_Studio1210/cplex/include"
  PATHS "/Applications/CPLEX_Studio1210/cplex/include"
  PATHS "/Applications/CPLEX_Studio221/cplex/include"
  PATHS "/Applications/CPLEX_Studio2211/cplex/include"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio1210/cplex/include"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio221/cplex/include"
  HINTS ${CPLEX_ROOT_DIR}/include
  HINTS ${CPLEX_ROOT_DIR}/cplex/include
)

find_path(CONCERT_INCLUDE_DIR
  ilconcert/ilomodel.h
  PATHS "C:/ILOG/CONCERT/include"
  PATHS "/opt/ilog/concert/include"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/concert/include"
  PATHS "~/ILOG/CPLEX_Studio1210/concert/include"
  PATHS "/Applications/CPLEX_Studio1210/concert/include"
  PATHS "/Applications/CPLEX_Studio221/concert/include"
  PATHS "/Applications/CPLEX_Studio2211/concert/include"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio1210/concert/include"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio221/concert/include"
  HINTS ${CPLEX_ROOT_DIR}/include
  HINTS ${CPLEX_ROOT_DIR}/concert/include
)

find_path(CPOPTIMIZER_INCLUDE_DIR
  ilcp/cp.h
  PATHS "C:/ILOG/CPOPTIMIZER/include"
  PATHS "/opt/ilog/cpoptimizer/include"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/cpoptimizer/include"
  PATHS "~/ILOG/CPLEX_Studio1210/cpoptimizer/include"
  PATHS "/Applications/CPLEX_Studio1210/cpoptimizer/include"
  PATHS "/Applications/CPLEX_Studio221/cpoptimizer/include"
  PATHS "/Applications/CPLEX_Studio2211/cpoptimizer/include"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio1210/cpoptimizer/include"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio221/cpoptimizer/include"
  HINTS ${CPLEX_ROOT_DIR}/include
  HINTS ${CPLEX_ROOT_DIR}/cpoptimizer/include
)

FIND_LIBRARY(CPLEX_LIBRARY
  cplex
  PATHS "C:/ILOG/CPLEX/lib/msvc7/stat_mda"
  PATHS "/opt/ilog/cplex/bin"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/cplex/lib/x86-64_sles10_4.1/static_pic/"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/cplex/bin"
  PATHS "~/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/"
  PATHS "/Applications/CPLEX_Studio1210/cplex/bin"
  PATHS "/Applications/CPLEX_Studio221/cplex/bin"
  PATHS "/Applications/CPLEX_Studio2211/cplex/bin"
  PATHS "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/"
  PATHS "/Applications/CPLEX_Studio221/cplex/bin/x86-64_osx/"
  PATHS "/Applications/CPLEX_Studio2211/cplex/bin/x86-64_osx/"
  PATHS "/Applications/CPLEX_Studio2211/cplex/bin/arm64_osx/"
  PATHS "/Applications/CPLEX_Studio1210/cplex/lib/x86-64_osx/static_pic/"
  PATHS "/Applications/CPLEX_Studio221/cplex/lib/x86-64_osx/static_pic/"
  PATHS "/Applications/CPLEX_Studio2211/cplex/lib/x86-64_osx/static_pic/"
  PATHS "/Applications/CPLEX_Studio2211/cplex/lib/arm64_osx/static_pic/"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio221/cplex/lib/x86-64_linux/static_pic/"

  HINTS ${CPLEX_ROOT_DIR}/bin
  HINTS ${CPLEX_ROOT_DIR}/cplex/bin
  HINTS ${CPLEX_ROOT_DIR}/lib
  HINTS ${CPLEX_ROOT_DIR}/cplex/lib
  HINTS ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_sles10_4.1/static_pic
  HINTS ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux/static_pic
)

FIND_LIBRARY(ILOCPLEX_LIBRARY
  ilocplex
  PATHS "C:/ILOG/CPLEX/lib/msvc7/stat_mda"
  PATHS "/opt/ilog/cplex/bin"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/cplex/lib/x86-64_sles10_4.1/static_pic/"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/cplex/bin"
  PATHS "~/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/"
  PATHS "/Applications/CPLEX_Studio1210/cplex/bin"
  PATHS "/Applications/CPLEX_Studio221/cplex/bin"
  PATHS "/Applications/CPLEX_Studio1210/cplex/lib/x86-64_osx/static_pic/"
  PATHS "/Applications/CPLEX_Studio2211/cplex/lib/x86-64_osx/static_pic/"
  PATHS "/Applications/CPLEX_Studio2211/cplex/lib/arm64_osx/static_pic/"
  PATHS "/Applications/CPLEX_Studio221/cplex/lib/x86-64_osx/static_pic/"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio221/cplex/lib/x86-64_linux/static_pic/"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic/"
  HINTS ${CPLEX_ROOT_DIR}/bin
  HINTS ${CPLEX_ROOT_DIR}/cplex/bin
  HINTS ${CPLEX_ROOT_DIR}/lib
  HINTS ${CPLEX_ROOT_DIR}/cplex/lib
  HINTS ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_sles10_4.1/static_pic
  HINTS ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux/static_pic
)

FIND_LIBRARY(CONCERT_LIBRARY
  concert
  PATHS "C:/ILOG/CONCERT/lib/msvc7/stat_mda"
  PATHS "/opt/ilog/concert/bin"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/concert/lib/x86-64_sles10_4.1/static_pic/"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/concert/bin"
  PATHS "~/ILOG/CPLEX_Studio1210/concert/lib/x86-64_linux/static_pic/"
  PATHS "/Applications/CPLEX_Studio1210/concert/bin"
  PATHS "/Applications/CPLEX_Studio221/concert/bin"
  PATHS "/Applications/CPLEX_Studio2211/concert/bin"
  PATHS "/Applications/CPLEX_Studio1210/concert/lib/x86-64_osx/static_pic/"
  PATHS "/Applications/CPLEX_Studio221/concert/lib/x86-64_osx/static_pic/"
  PATHS "/Applications/CPLEX_Studio2211/concert/lib/x86-64_osx/static_pic/"
  PATHS "/Applications/CPLEX_Studio2211/concert/lib/arm64_osx/static_pic/"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio1210/concert/lib/x86-64_linux/static_pic/"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio221/concert/lib/x86-64_linux/static_pic/"
  HINTS ${CPLEX_ROOT_DIR}/bin
  HINTS ${CPLEX_ROOT_DIR}/concert/bin
  HINTS ${CPLEX_ROOT_DIR}/lib
  HINTS ${CPLEX_ROOT_DIR}/concert/lib
  HINTS ${CPLEX_ROOT_DIR}/concert/lib/x86-64_sles10_4.1/static_pic
  HINTS ${CPLEX_ROOT_DIR}/concert/lib/x86-64_linux/static_pic
)

FIND_LIBRARY(CPOPTIMIZER_LIBRARY
  cp
  PATHS "C:/ILOG/CPOPTIMIZER/lib/msvc7/stat_mda"
  PATHS "/opt/ilog/cpoptimizer/bin"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/cpoptimizer/lib/x86-64_sles10_4.1/static_pic/"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/cpoptimizer/bin"
  PATHS "~/ILOG/CPLEX_Studio1210/cpoptimizer/lib/x86-64_linux/static_pic/"
  PATHS "/Applications/CPLEX_Studio1210/cpoptimizer/bin"
  PATHS "/Applications/CPLEX_Studio221/cpoptimizer/bin"
  PATHS "/Applications/CPLEX_Studio2211/cpoptimizer/bin"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio1210/cpoptimizer/bin"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio221/cpoptimizer/bin"
  HINTS ${CPLEX_ROOT_DIR}/bin
  HINTS ${CPLEX_ROOT_DIR}/cpoptimizer/bin
  HINTS ${CPLEX_ROOT_DIR}/lib
  HINTS ${CPLEX_ROOT_DIR}/cpoptimizer/lib
  HINTS ${CPLEX_ROOT_DIR}/cpoptimizer/lib/x86-64_sles10_4.1/static_pic
  HINTS ${CPLEX_ROOT_DIR}/cpoptimizer/lib/x86-64_linux/static_pic
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(cplex DEFAULT_MSG CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

find_path(CPLEX_BIN_DIR
  cplex.dll
  PATHS "C:/ILOG/CPLEX/bin/x86_win32"
  PATHS "~/ILOG/CPLEX_Studio_AcademicResearch122/cplex/bin"
  PATHS "~/ILOG/CPLEX_Studio1210/cplex/bin"
  PATHS "/Applications/CPLEX_Studio1210/cplex/bin"
  PATHS "/Applications/CPLEX_Studio221/cplex/bin"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin"
  PATHS "/opt/ibm/ILOG/CPLEX_Studio221/cplex/bin"
  HINTS ${CPLEX_ROOT_DIR}/bin
  HINTS ${CPLEX_ROOT_DIR}/cplex/bin
)

IF(CPLEX_FOUND)
  SET(CPLEX_INCLUDE_DIRS "${CPLEX_INCLUDE_DIR};${CONCERT_INCLUDE_DIR};${CPOPTIMIZER_INCLUDE_DIR}")
  SET(CPLEX_LIBRARIES "${ILOCPLEX_LIBRARY};${CPLEX_LIBRARY};${CONCERT_LIBRARY};pthread;z")
  IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    SET(CPLEX_LIBRARIES "${CPLEX_LIBRARIES};m;pthread")
  ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
ENDIF(CPLEX_FOUND)

MARK_AS_ADVANCED(CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_BIN_DIR)

IF(CPLEX_FOUND)
  SET(LEMON_HAVE_LP TRUE)
  SET(LEMON_HAVE_MIP TRUE)
  SET(LEMON_HAVE_CPLEX TRUE)
ENDIF(CPLEX_FOUND)
