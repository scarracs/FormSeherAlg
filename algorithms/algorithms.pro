#-------------------------------------------------
#
# Project created by QtCreator 2014-05-13T04:07:07
#
#-------------------------------------------------

include(../global.pri)
QT       -= core gui

TARGET = algorithms01
TEMPLATE = lib
VERSION = 0.1
CONFIG += staticlib
DESTDIR = $${ALG_BIN_DIR}

DEFINES += ALGORITHMS_LIBRARY

SOURCES += $${ALG_SRC_DIR}/algorithm.cpp \
    $${ALG_SRC_DIR}/line.cpp \
    $${ALG_SRC_DIR}/edl/edl.cpp \
    $${ALG_SRC_DIR}/mathutil.cpp

HEADERS += $${ALG_INCL_DIR}/algorithm.h \
    $${ALG_INCL_DIR}/line.h \
    $${ALG_INCL_DIR}/edl/edl.h \
    $${ALG_INCL_DIR}/mathutil.h

INCLUDEPATH += $${ALG_INCL_DIR}


#-------------------------------------------------
#               win32 specifics
#-------------------------------------------------
win32{
INCLUDEPATH += $(OPENCV_DIR_INCLUDE)
LIBS += -L$(OPENCV_DIR_LIB)
LIBS += -llibopencv_core248 -llibopencv_imgproc248 -llibopencv_highgui248
}
#-------------------------------------------------
#               Linux/Unix specifics
#-------------------------------------------------
linux{
LIBS += -lopencv_core -lopencv_imgproc
}
