#include "resultgrid.h"
#include <QString>
#include <QStringList>
#include <string>
#include <System.h>
#include "vtk.h"

using namespace std;

ResultGrid::ResultGrid():TimeSeriesSet<double>()
{

}

ResultGrid::ResultGrid(const ResultGrid& rhs):TimeSeriesSet<double>(rhs)
{
    Positions = rhs.Positions;
}

ResultGrid::~ResultGrid()
{

}

ResultGrid& ResultGrid::operator=(const ResultGrid &rhs)
{
    TimeSeriesSet::operator=(rhs);
    Positions = rhs.Positions;
    return *this;
}

ResultGrid::ResultGrid(const TimeSeriesSet<double> &cts, const vector<string> &components, const string &quantity)
{
    for (unsigned int i=0; i<components.size(); i++)
    {
        for (unsigned int j=0; j<cts.size(); j++)
            if (cts.getSeriesName(j)==components[i] + "_" + quantity)
                append(cts[j],components[i]);
    }
}

TimeSeries<double> ResultGrid::Sum()
{
    TimeSeries<double> out;
    out = operator[](0);
    for (unsigned int i = 1; i<size(); i++)
    {
        out %= operator[](i);
    }
    return out;
}

TimeSeries<double> ResultGrid::SumIntegrate()
{
    return Sum().integrate();
}

ResultGrid::ResultGrid(const TimeSeriesSet<double> &cts, const string &quantity, System *system)
{
    for (int i=0; i<cts.size(); i++)
    {
        string block_name = QString::fromStdString(cts.getSeriesName(i)).split("_")[0].toStdString();
        string quan = QString::fromStdString(cts.getSeriesName(i)).split("_")[1].toStdString();
        if (quan==quantity)
        {
            point pt;
            if (system->block(block_name))
            {
                pt.x = system->block(block_name)->GetVal("act_X");
                pt.y = system->block(block_name)->GetVal("act_Y");
                Positions.push_back(pt);
                append(cts[i],block_name);
            }

        }
    }
}

ResultGrid::ResultGrid(const string &quantity, System *system)
{
    for (int i=0; i<system->BlockCount(); i++)
    {
        string block_name = system->block(i)->GetName();

        if (system->block(i)->HasQuantity(quantity))
        {
            point pt;
            if (system->block(block_name))
            {
                TimeSeries<double> value;
                value.append(0,system->block(block_name)->GetVal(quantity,Expression::timing::past));
                pt.x = system->block(block_name)->GetVal("act_X");
                pt.y = system->block(block_name)->GetVal("act_Y");
                Positions.push_back(pt);
                append(value,block_name);
            }

        }
    }
}

#ifdef use_VTK
void ResultGrid::WriteToVTP(const std::string &quanname, const std::string &filename, int i, const double &scale) const
{
    vtkSmartPointer<vtkPoints> points_3 =
            vtkSmartPointer<vtkPoints>::New();

        double xx, yy, zz;
        vtkSmartPointer<vtkFloatArray> values =
            vtkSmartPointer<vtkFloatArray>::New();

        values->SetNumberOfComponents(1);

        values->SetName(quanname.c_str());
        //cout<<quanname<<endl;



        for (unsigned int j = 0; j < size(); j++)
        {
            //cout<<"Positions "<<j<<endl;
            yy = Positions[j].y;
            xx = Positions[j].x;
            zz = operator[](j).getValue(i)*scale;
            //cout<<"Positions "<<j<<":"<<xx<<":"<<yy<<":"<<zz<<" done!"<<endl;
            float tt = float(operator[](j).getValue(i));
            float t[1] = { tt };
            //cout<<"0.1"<<endl;
            points_3->InsertNextPoint(xx, yy, zz);
            //cout<<"0.2"<<endl;
            //cout<<tt<<":"<<t<<endl;
            values->InsertNextTuple(t);
            //cout<<"1"<<endl;
        }


        // Add the grid points to a polydata object
        vtkSmartPointer<vtkPolyData> inputPolyData =
            vtkSmartPointer<vtkPolyData>::New();
        inputPolyData->SetPoints(points_3);

        // Triangulate the grid points
        vtkSmartPointer<vtkDelaunay2D> delaunay =
            vtkSmartPointer<vtkDelaunay2D>::New();
    #if VTK_MAJOR_VERSION <= 5
        delaunay->SetInput(inputPolyData);
    #else
        delaunay->SetInputData(inputPolyData);
    #endif
        delaunay->Update();
        vtkPolyData* outputPolyData = delaunay->GetOutput();

        double bounds[6];
        outputPolyData->GetBounds(bounds);

                outputPolyData->GetPointData()->SetScalars(values);


        //Append the two meshes
        vtkSmartPointer<vtkAppendPolyData> appendFilter =
            vtkSmartPointer<vtkAppendPolyData>::New();
    #if VTK_MAJOR_VERSION <= 5
        appendFilter->AddInputConnection(input1->GetProducerPort());
        appendFilter->AddInputConnection(input2->GetProducerPort());
    #else
        appendFilter->AddInputData(outputPolyData);
    #endif
        appendFilter->Update();
        vtkSmartPointer<vtkPolyDataMapper> mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
        mapper->SetInputConnection(polydata->GetProducerPort());
    #else
        mapper->SetInputConnection(appendFilter->GetOutputPort());
        //mapper->SetInputData(polydata_1);
    #endif

        vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetPointSize(5);

        vtkSmartPointer<vtkXMLPolyDataWriter> writer =
            vtkSmartPointer<vtkXMLPolyDataWriter>::New();

        writer->SetFileName(filename.c_str());
        writer->SetInputData(mapper->GetInput());

        writer->SetDataModeToAscii();
        writer->Write();


}

void ResultGrid::WriteToVTP(const std::string &quanname, const std::string &filename, const double &scale) const
{
    for (unsigned k=0; k<operator[](0).size(); k++)
    {
        cout<<"snapshot: "<<k<<endl;
        std::string name = aquiutils::split(filename,'.')[0]+"_"+aquiutils::numbertostring(k+1,4)+".vtp";
        cout<<name<<endl;
        WriteToVTP(quanname, name, k, scale);
    }
}
#endif
