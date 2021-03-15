TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        cell.cpp \
        cell_block.cpp \
        main.cpp \
        md.cpp \
        neighbor.cpp \
        output.cpp \
        particle.cpp \
        particle_ops.cpp \
        slicer.cpp \
        vector3.cpp

HEADERS += \
    cell.h \
    cell_block.h \
    cell_slice_2d.h \
    extent.h \
    index.h \
    md.h \
    mesh3d.h \
    neighbor.h \
    output.h \
    particle.h \
    particle_ops.h \
    slicer.h \
    vector3.h

unix {
    QMAKE_CXX = mpicxx
    QMAKE_LINK = mpicxx
    QMAKE_CC = mpicc
    INCLUDEPATH += "/usr/lib/x86_64-linux-gnu/openmpi/include"
}

