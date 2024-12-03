#include "resultgrid.h"
#include <QString>
#include <QStringList>
#include <string>
#include <System.h>
#include "vtk.h"

using namespace std;

ResultGrid::ResultGrid():CTimeSeriesSet<double>()
{

}

ResultGrid::ResultGrid(const ResultGrid& rhs):CTimeSeriesSet<double>(rhs)
{
    Positions = rhs.Positions;
}
ResultGrid::~ResultGrid()
{

}

ResultGrid& ResultGrid::operator=(const ResultGrid &rhs)
{
    CTimeSeriesSet::operator=(rhs);
    Positions = rhs.Positions;
    return *this;
}

ResultGrid::ResultGrid(const CTimeSeriesSet<double> &cts, const vector<string> &components, const string &quantity)
{
    for (unsigned int i=0; i<components.size(); i++)
    {
        if (cts.names[i]==components[i] + "_" + quantity)
            append(cts.BTC[i],components[i]);
    }
}

CTimeSeries<double> ResultGrid::Sum()
{
    CTimeSeries<double> out;
    out = BTC[0];
    for (unsigned int i = 1; i<BTC.size(); i++)
    {
        out %= BTC[i];
    }
    return out;
}

CTimeSeries<double> ResultGrid::SumIntegrate()
{
    return Sum().integrate();
}

ResultGrid::ResultGrid(const CTimeSeriesSet<double> &cts, const string &quantity, System *system)
{
    for (int i=0; i<cts.nvars; i++)
    {
        string block_name = QString::fromStdString(cts.names[i]).split("_")[0].toStdString();
        string quan = QString::fromStdString(cts.names[i]).split("_")[1].toStdString();
        if (quan==quantity)
        {
            point pt;
            if (system->block(block_name))
            {
                pt.x = system->block(block_name)->GetVal("act_X");
                pt.y = system->block(block_name)->GetVal("act_Y");
                Positions.push_back(pt);
                append(cts.BTC[i],block_name);
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



        for (unsigned int j = 0; j < nvars; j++)
        {
            //cout<<"Positions "<<j<<endl;
            yy = Positions[j].y;
            xx = Positions[j].x;
            zz = BTC[j].GetC(i)*scale;
            //cout<<"Positions "<<j<<":"<<xx<<":"<<yy<<":"<<zz<<" done!"<<endl;
            float tt = float(BTC[j].GetC(i));
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
    for (unsigned k=0; k<BTC[0].n; k++)
    {
        cout<<"snapshot: "<<k<<endl;
        std::string name = aquiutils::split(filename,'.')[0]+"_"+aquiutils::numbertostring(k+1,4)+".vtp";
        cout<<name<<endl;
        WriteToVTP(quanname, name, k, scale);
    }
}
#endif
