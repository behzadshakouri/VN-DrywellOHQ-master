#pragma once
//#include "stdafx.h"
#include "vtkAutoInit.h" 
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);


#ifdef use_VTK
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#endif

// For compatibility with new VTK generic data arrays
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif
