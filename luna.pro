#-------------------------------------------------
#
# Project created by QtCreator 2016-09-19T22:09:39
#
#-------------------------------------------------

#QT       -= gui
win32:LIBS += -L$$PWD//WpdPack/Lib/ -lPacket
win32:LIBS += -L$$PWD//WpdPack/Lib/ -lwpcap
TARGET = luna
TEMPLATE = lib

DEFINES += LUNA_LIBRARY

SOURCES += luna.cpp

HEADERS += luna.h\
        luna_global.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
