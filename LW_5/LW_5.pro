#-------------------------------------------------
#
# Project created by QtCreator 2018-10-24T18:12:19
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = LW_5
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++17

QMAKE_CXXFLAGS += \
               -std=c++17 \
               -fexceptions \
               -Wall \
               -Werror \
               -Wextra \
               -pedantic \
               -Wfatal-errors \
               -pedantic-errors \
               -Winit-self \
               -Wnon-virtual-dtor \
               -Winline \
               -Wmissing-include-dirs \
               -Wredundant-decls \
               -Wfloat-equal \
               -Wmain \
               -Wunreachable-code \
               -Wshadow \
               -Wcast-align \
               -Wswitch-enum

win32-g++ {
    INCLUDEPATH += C:\Qt\boost_1_68_0\mingw530_32\include\boost-1_68
}

SOURCES += \
        main.cpp \
        mainwindow.cpp

HEADERS += \
        mainwindow.h \
    parabolic_pde.hpp

FORMS += \
        mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
