#-------------------------------------------------
#
# Project created by QtCreator 2017-12-11T16:12:07
#
#-------------------------------------------------

# Uncomment the line below to enable debug STL (with more checks on iterator invalidation etc.)
# NOTE 1: Enabling debug STL mode will make the performance WORSE. So don't enable it when running performance tests!
# NOTE 2: If you uncomment or recomment the line, remember to recompile EVERYTHING by selecting
# "Rebuild all" from the Build menu
#QMAKE_CXXFLAGS += -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC

QT       += core gui

CONFIG += c++17 warn_on

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = prg1
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
    datastructures.cc \
    mainwindow.cc \
    mainprogram.cc

HEADERS += \
    datastructures.hh \
    mainwindow.hh \
    mainprogram.hh

FORMS += \
    mainwindow.ui


# If you uncomment the lines below and recompile EVERYTHING (by selecting "Rebuild all" from the Build menu),
# you'll get a non-graphical command line version of the program (just like if you had compiled the program
# directly using the g++ command.
#FORMS -= mainprogram.ui
#CONFIG -= core gui qt
#CONFIG += console

DISTFILES += \
    example-all-in.txt \
    example-all-out.txt \
    example-areas.txt \
    example-compulsory-in.txt \
    example-compulsory-out.txt \
    example-places.txt \
    helv-areas.txt \
    helv-places.txt \
    hervanta-south-areas.txt \
    hervanta-south-places.txt \
    kintulammi-areas.txt \
    kintulammi-places.txt \
    kintulammi-test-all-in.txt \
    kintulammi-test-all-out.txt \
    kintulammi-test-compulsory-in.txt \
    kintulammi-test-compulsory-out.txt \
    perftest-access.txt \
    perftest-all.txt \
    perftest-all_subareas.txt \
    perftest-change.txt \
    perftest-common_area.txt \
    perftest-compulsory.txt \
    perftest-find_places.txt \
    perftest-places_closest_to.txt \
    perftest-remove.txt \
    perftest-sorting.txt \
    perftest-subareas.txt \
    simpletest-all-in.txt \
    simpletest-all-out.txt \
    simpletest-compulsory-in.txt \
    simpletest-compulsory-out.txt
