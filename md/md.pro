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
        vector3.cpp

HEADERS += \
    cell.h \
    cell_block.h \
    extent.h \
    index.h \
    md.h \
    mesh3d.h \
    neighbor.h \
    output.h \
    particle.h \
    particle_ops.h \
    vector3.h
