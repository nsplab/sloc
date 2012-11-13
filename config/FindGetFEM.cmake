# FindGetFEM.cmake
#
if(NOT (GETFEM_CFLAGS AND GETFEM_LIBS))

    # XXX: how to check cache for GETFEM_PREFIX?

    if(NOT DEFINED GETFEM_PREFIX)
        set(GETFEM_PREFIX "/usr/local")
    endif()

    set(GETFEM_CONFIG "${GETFEM_PREFIX}/bin/getfem-config")
    #message(STATUS "GETFEM_CONFIG = ${GETFEM_CONFIG}")

    execute_process(COMMAND "${GETFEM_CONFIG}" --cflags OUTPUT_VARIABLE GETFEM_CFLAGS)
    string(STRIP ${GETFEM_CFLAGS} GETFEM_CFLAGS)

    execute_process(COMMAND "${GETFEM_CONFIG}" --libs OUTPUT_VARIABLE GETFEM_LIBS)
    string(STRIP ${GETFEM_LIBS} GETFEM_LIBS)

    if(GETFEM_CFLAGS AND GETFEM_LIBS)
        set(GETFEM_FOUND true)
    endif()

    message(STATUS "Found these GetFEM settings:")
    message(STATUS "  GETFEM_CFLAGS = ${GETFEM_CFLAGS}")
    message(STATUS "  GETFEM_LIBS   = ${GETFEM_LIBS}")

endif()

