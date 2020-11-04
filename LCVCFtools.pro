TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++14

unix: LIBS += -lboost_iostreams
unix|win32: LIBS += -lz

HEADERS += \
    src/lcvcftools.h

SOURCES += \
    src/lcvcftools.cpp

CONFIG(debug, debug|release) {
    DESTDIR = bin/debug
}
CONFIG(release, debug|release) {
    DESTDIR = bin
}
OBJECTS_DIR = $$DESTDIR
