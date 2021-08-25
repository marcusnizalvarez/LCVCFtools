TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -std=c++11
QMAKE_LFLAGS += -static
LIBS += -lboost_iostreams
LIBS += -lz
HEADERS += \
    src/lcvcftools.h
SOURCES += \
    src/lcvcftools.cpp \
    src/main.cpp
CONFIG(debug, debug|release) {
    OBJECTS_DIR = build/debug
}
CONFIG(release, debug|release) {
    OBJECTS_DIR = build/release
}
