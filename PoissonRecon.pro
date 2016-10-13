INCLUDEPATH += ../../../../vcglib
INCLUDEPATH -= .
CONFIG += console stl
TEMPLATE = app
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
CONFIG += console warn_off
SOURCES +=  Src/PoissonRecon.cpp \
            Src/MarchingCubes.cpp \
            Src/PlyFile.cpp \
            Src/CmdLineParser.cpp \
            Src/Factor.cpp \
            Src/Geometry.cpp

QMAKE_CXXFLAGS+=-fopenmp -Wsign-compare
macx:QMAKE_CXXFLAGS-= -fopenmp

TARGET=PoissonRecon

