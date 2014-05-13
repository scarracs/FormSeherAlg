include(../global.pri)
QT       += core gui
TEMPLATE = app

CONFIG += console

CONFIG -= qt

SOURCES += main.cpp

#-------------------------------------------------
#               win32 specifics
#-------------------------------------------------
win32{
INCLUDEPATH += $(OPENCV_DIR_INCLUDE)
INCLUDEPATH += $${ALG_INCL_DIR}

LIBS += -L$(OPENCV_DIR_LIB)
LIBS += -llibopencv_core248 -llibopencv_imgproc248 -llibopencv_highgui248
LIBS += -L$${ALG_BIN_DIR} -lalgorithms01
}
#-------------------------------------------------
#               Linux/Unix specifics
#-------------------------------------------------
linux{
LIBS += -lopencv_core -lopencv_imgproc
}
