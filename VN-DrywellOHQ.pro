
TEMPLATE = app

QT -= gui
QT += core

CONFIG += console
CONFIG -= app_bundle

CONFIG += c++14

DEFINES += GSL

CONFIG += Behzad
DEFINES += Behzad


#CONFIG += PowerEdge
#DEFINES += PowerEdge


#CONFIG += Arash
#DEFINES += Arash

Behzad {
    OHQPATH = /home/behzad/Projects/OpenHydroQual/aquifolium
    VTKBUILDPATH = /home/behzad/Projects/VTK-9.3.1/VTK-build
    VTKHEADERPATH = /home/behzad/Projects/VTK-9.3.1
    VTK_V = -9.3
}

PowerEdge {
    OHQPATH = ../OpenHydroQual/aquifolium
    VTKBUILDPATH = ../VTK-build
    VTKHEADERPATH = ../VTK
    VTK_V = -9.0
}

Arash {
    OHQPATH = /home/arash/Projects/QAquifolium/aquifolium
    VTKBUILDPATH = /home/arash/Projects/VTK/VTK-build
    VTKHEADERPATH = /home/arash/Projects/VTK
    VTK_V = -9.0
}

DEFINES += use_VTK ARMA_USE_SUPERLU
CONFIG += use_VTK

message("DEFINES: $$DEFINES")
message("VTK: $$VTKHEADERPATH")

INCLUDEPATH += $$OHQPATH/include
INCLUDEPATH += $$OHQPATH/src
INCLUDEPATH += $$OHQPATH/include/GA
INCLUDEPATH += $$OHQPATH/include/MCMC
INCLUDEPATH += ../jsoncpp/include/



if==macx:CONFIG += staticlib
macx: DEFINES +=mac_version
linux: DEFINES +=ubuntu_version
win32: DEFINES +=windows_version

DEFINES += Terminal_version

TARGET = DrywellOHQ

win32:QMAKE_CXXFLAGS += /MP

macx: {
    QMAKE_CXXFLAGS += -Xpreprocessor -fopenmp -lomp -Iusr/local/lib/
}

macx: {
    QMAKE_LFLAGS += -lomp
}

macx: {
    LIBS += -L /usr/local/lib /usr/local/lib/libomp.dylib
}

macx: {
    INCLUDEPATH += /usr/local/include/
}


CONFIG(debug, debug|release) {
    message(Building in debug mode)
    !macx: QMAKE_CXXFLAGS *= "-Xpreprocessor -fopenmp"
    !macx: QMAKE_LFLAGS +=  -fopenmp
    !macx: LIBS += -lgomp -lpthread
    LIBS += -lpthread
    DEFINES += NO_OPENMP DEBUG

} else {
    message(Building in release mode)
    !macx: QMAKE_CXXFLAGS *= "-Xpreprocessor -fopenmp"
    !macx: QMAKE_LFLAGS +=  -fopenmp
    # QMAKE_CFLAGS+=-pg
    # QMAKE_CXXFLAGS+=-pg
    # QMAKE_LFLAGS+=-pg
    # macx: DEFINES += NO_OPENMP
    ! macx: LIBS += -lgomp -lpthread
    macx: LIBS += -lpthread
}




SOURCES += \
        $$OHQPATH/src/Block.cpp \
        $$OHQPATH/src/Command.cpp \
        $$OHQPATH/src/Condition.cpp \
        $$OHQPATH/src/ErrorHandler.cpp \
        $$OHQPATH/src/Expression.cpp \
        $$OHQPATH/src/Link.cpp \
        $$OHQPATH/src/Matrix.cpp \
        $$OHQPATH/src/Matrix_arma.cpp \
        $$OHQPATH/src/MetaModel.cpp \
        $$OHQPATH/src/NormalDist.cpp \
        $$OHQPATH/src/Object.cpp \
        $$OHQPATH/src/Objective_Function.cpp \
        $$OHQPATH/src/Objective_Function_Set.cpp \
        $$OHQPATH/src/Parameter.cpp \
        $$OHQPATH/src/Parameter_Set.cpp \
        $$OHQPATH/src/Precipitation.cpp \
        $$OHQPATH/src/Quan.cpp \
        $$OHQPATH/src/QuanSet.cpp \
        $$OHQPATH/src/QuickSort.cpp \
        $$OHQPATH/src/Rule.cpp \
        $$OHQPATH/src/RxnParameter.cpp \
        $$OHQPATH/src/Script.cpp \
        $$OHQPATH/src/Source.cpp \
        $$OHQPATH/src/System.cpp \
        $$OHQPATH/src/Utilities.cpp \
        $$OHQPATH/src/Vector.cpp \
        $$OHQPATH/src/Vector_arma.cpp \
        $$OHQPATH/src/constituent.cpp \
        $$OHQPATH/src/observation.cpp \
        $$OHQPATH/src/precalculatedfunction.cpp \
        $$OHQPATH/src/reaction.cpp \
        $$OHQPATH/src/restorepoint.cpp \
        $$OHQPATH/src/solutionlogger.cpp \
        $$OHQPATH/src/GA/Binary.cpp \
        $$OHQPATH/src/GA/Individual.cpp \
        $$OHQPATH/src/GA/DistributionNUnif.cpp \
        $$OHQPATH/src/GA/Distribution.cpp \
        ../jsoncpp/src/lib_json/json_reader.cpp \
        ../jsoncpp/src/lib_json/json_value.cpp \
        ../jsoncpp/src/lib_json/json_writer.cpp \
        main.cpp \
        modelcreator.cpp \
        resultgrid.cpp

HEADERS += \
    $$OHQPATH/include/Objective_Function.h \
    $$OHQPATH/include/Objective_Function_Set.h \
    $$OHQPATH/include/Precipitation.h \
    $$OHQPATH/include/RxnParameter.h \
    $$OHQPATH/include/constituent.h \
    $$OHQPATH/include/observation.h \
    $$OHQPATH/include/precalculatedfunction.h \
    $$OHQPATH/include/solutionlogger.h \
    $$OHQPATH/include/GA/GA.h \
    $$OHQPATH/include/MCMC/MCMC.h \
    $$OHQPATH/include/MCMC/MCMC.hpp \
    $$OHQPATH/include/Utilities.h \
    $$OHQPATH/include/restorepoint.h \
    $$OHQPATH/include/safevector.h \
    $$OHQPATH/include/safevector.hpp \
    $$OHQPATH/include/Block.h \
    $$OHQPATH/include/BTC.h \
    $$OHQPATH/include/BTCSet.h \
    $$OHQPATH/include/Expression.h \
    $$OHQPATH/include/Link.h \
    $$OHQPATH/include/Matrix.h \
    $$OHQPATH/include/Matrix_arma.h \
    $$OHQPATH/include/MetaModel.h \
    $$OHQPATH/include/NormalDist.h \
    $$OHQPATH/include/Object.h \
    $$OHQPATH/include/Quan.h \
    $$OHQPATH/include/QuanSet.h \
    $$OHQPATH/include/QuickSort.h \
    $$OHQPATH/include/StringOP.h \
    $$OHQPATH/include/System.h \
    $$OHQPATH/include/Vector.h \
    $$OHQPATH/include/Vector_arma.h \
    ../jsoncpp/include/json/allocator.h \
    ../jsoncpp/include/json/assertions.h \
    ../jsoncpp/include/json/autolink.h \
    ../jsoncpp/include/json/config.h \
    ../jsoncpp/include/json/features.h \
    ../jsoncpp/include/json/forwards.h \
    ../jsoncpp/include/json/json.h \
    ../jsoncpp/include/json/reader.h \
    ../jsoncpp/include/json/value.h \
    ../jsoncpp/include/json/version.h \
    ../jsoncpp/include/json/writer.h \
    ../jsoncpp/src/lib_json/json_tool.h \
    ../jsoncpp/src/lib_json/version.h.in \
    $$OHQPATH/include/Parameter.h \
    $$OHQPATH/include/Parameter_Set.h \
    $$OHQPATH/include/Command.h \
    $$OHQPATH/include/Script.h \
    $$OHQPATH/include/GA/Binary.h \
    $$OHQPATH/include/GA/Distribution.h \
    $$OHQPATH/include/GA/DistributionNUnif.h \
    $$OHQPATH/include/GA/Individual.h \
    $$OHQPATH/include/Objective_Function.h \
    $$OHQPATH/include/Objective_Function_Set.h \
    $$OHQPATH/include/GA/GA.hpp \
    $$OHQPATH/src/BTC.hpp \
    $$OHQPATH/src/BTCSet.hpp \
    $$OHQPATH/include/reaction.h \
    modelcreator.h \
    resultgrid.h


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

# LAPACK — Linear Algebra PACKage lib and include locations

win32 {

    LAPACK_INCLUDE = $$PWD/include
    #64 bits build
    contains(QMAKE_TARGET.arch, x86_64) {
        #debug
        CONFIG(debug, debug|release) {
            LAPACK_LIB_DIR = $$PWD/libs/lapack-blas_lib_win64/debug
            LIBS +=  -L$${LAPACK_LIB_DIR} -llapack_win64_MTd \
                    -lblas_win64_MTd
        }
        #release
        CONFIG(release, debug|release) {
            LAPACK_LIB_DIR = $$PWD/libs/lapack-blas_lib_win64/release
            LIBS +=  -L$${LAPACK_LIB_DIR} -llapack_win64_MT \
                    -lblas_win64_MT
        }
    }

    INCLUDEPATH += $${LAPACK_INCLUDE}

    DEFINES += ARMA_USE_LAPACK ARMA_USE_BLAS

}

linux {
    #sudo apt-get install libblas-dev liblapack-dev
     DEFINES += ARMA_USE_LAPACK ARMA_USE_BLAS
     LIBS += -larmadillo -llapack -lblas -lgsl
}

macx {

    LIBS += /opt/homebrew/Cellar/armadillo/11.4.2/lib/libarmadillo.dylib
    INCLUDEPATH += $$PWD/../../../../opt/homebrew/Cellar/armadillo/11.4.2/include
    DEPENDPATH += $$PWD/../../../../opt/homebrew/Cellar/armadillo/11.4.2/include

}


use_VTK {
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkChartsCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonColor$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonComputationalGeometry$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonDataModel$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonExecutionModel$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonMath$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonMisc$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonSystem$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonTransforms$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkDICOMParser$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkexpat$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersAMR$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersExtraction$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersFlowPaths$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersGeneral$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersGeneric$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersGeometry$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersHybrid$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersHyperTree$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersImaging$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersModeling$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersParallel$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersParallelImaging$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersPoints$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersProgrammable$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersSelection$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersSMP$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersSources$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersStatistics$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersTexture$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersTopology$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersVerdict$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkfreetype$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkGeovisCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkgl2ps$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkglew$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkhdf5$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkhdf5_hl$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingColor$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingFourier$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingGeneral$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingHybrid$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingMath$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingMorphological$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingSources$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingStatistics$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingStencil$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInfovisCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInfovisLayout$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInteractionImage$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInteractionStyle$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInteractionWidgets$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOAMR$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOEnSight$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOExodus$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOGeometry$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOImage$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOImport$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOInfovis$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOLegacy$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOLSDyna$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOMINC$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOMovie$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIONetCDF$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOParallel$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOParallelXML$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOPLY$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOSQL$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOTecplotTable$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOVideo$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOXML$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOXMLParser$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkjpeg$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkjsoncpp$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtklibharu$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtklibxml2$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtklz4$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkmetaio$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkParallelCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkpng$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingAnnotation$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingContext2D$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingFreeType$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingGL2PSOpenGL2$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingImage$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingLabel$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingLOD$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingOpenGL2$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingVolume$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingVolumeOpenGL2$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtksqlite$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtksys$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtktiff$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkverdict$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkViewsContext2D$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkViewsCore$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkViewsInfovis$$VTK_V
    LIBS += -L$$VTKBUILDPATH/lib/ -lvtkzlib$$VTK_V
    LIBS += -L"/usr/local/lib/ -lsuperlu.so"

    #VTK Include files
    INCLUDEPATH +=$${VTKHEADERPATH}
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/Core
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/Color
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/Transforms
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/Transforms
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/Color
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/DataModel
    INCLUDEPATH +=$${VTKBUILDPATH}/Utilities/KWIML
    INCLUDEPATH +=$${VTKHEADERPATH}/Utilities/KWIML
    INCLUDEPATH +=$${VTKHEADERPATH}/Rendering/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Rendering/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Core
    INCLUDEPATH +=$${VTKHEADERPATH}/Charts/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Charts/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Filters/General
    INCLUDEPATH +=$${VTKBUILDPATH}/Rendering/Context2D
    INCLUDEPATH +=$${VTKHEADERPATH}/Rendering/Context2D
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/DataModel
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/Math
    INCLUDEPATH +=$${VTKHEADERPATH}/Views/Context2D
    INCLUDEPATH +=$${VTKBUILDPATH}/Views/Context2D
    INCLUDEPATH +=$${VTKBUILDPATH}/Views/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Interaction/Widgets
    INCLUDEPATH +=$${VTKHEADERPATH}/Views/Core
    INCLUDEPATH +=$${VTKHEADERPATH}/Interaction/Style
    INCLUDEPATH +=$${VTKBUILDPATH}/Interaction/Style
    INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Modeling
    INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Modeling
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/ExecutionModel
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/ExecutionModel
    INCLUDEPATH +=$${VTKHEADERPATH}/Interaction/Widgets/
    INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Core/
    INCLUDEPATH +=$${VTKHEADERPATH}/Common/Misc/
    INCLUDEPATH +=$${VTKBUILDPATH}/Common/Misc
    INCLUDEPATH +=$${VTKHEADERPATH}/IO/XML/
    INCLUDEPATH +=$${VTKBUILDPATH}/IO/XML
    INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Sources
    INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Sources
    INCLUDEPATH +=$${VTKHEADERPATH}/Filters/General
    INCLUDEPATH +=$${VTKHEADERPATH}/IO/Image
    INCLUDEPATH +=$${VTKBUILDPATH}/IO/Image
    INCLUDEPATH +=$${VTKHEADERPATH}/Imaging/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Imaging/Core
    INCLUDEPATH +=$${VTKBUILDPATH}/Utilities/KWSys

}

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    $$OHQPATH/src/BTC.hpp \
    $$OHQPATH/src/BTCSet.hpp



