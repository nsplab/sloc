# FindDealII.cmake
#
if(NOT DEALII_CFLAGS)

    if(NOT DEFINED DEALII_PREFIX)
        set(DEALII_PREFIX "/usr/local")
    endif()

    set(DEALII_INCLUDE_DIR "${DEALII_PREFIX}/include")

    set(DEALII_CFLAGS
        "-I${DEALII_PREFIX}/include"
        "-I${DEALII_PREFIX}/build/include"
        "-I${DEALII_PREFIX}/bundled/tbb30_104oss/include"
    )

    set(DEALII_LDFLAGS
        "-L${DEALII_PREFIX}/build/source"
        "-lpthread"
        "-lz"
    )
    if(APPLE)
        set(DEALII_LDFLAGS ${DEALII_LDFLAGS}
            ${DEALII_PREFIX}/build/source/libdeal_II.dylib
        )
    else()
        set(DEALII_LDFLAGS ${DEALII_LDFLAGS}
            ${DEALII_PREFIX}/build/source/libdeal_II.so
        )
    endif()

    #set(DEALII_LIBS "-lz -lboost_thread-mt -lboost_serialization-mt -lboost_system-mt")

    if(DEALII_CFLAGS)
        set(DEALII_FOUND true)
    endif()

    message(STATUS "Found these deal.II settings:")
    message(STATUS "  DEALII_PREFIX = ${DEALII_PREFIX}")
    message(STATUS "  DEALII_CFLAGS = ${DEALII_CFLAGS}")
    message(STATUS "  DEALII_LDFLAGS = ${DEALII_LDFLAGS}")

endif()
