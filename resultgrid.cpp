#include "resultgrid.h"
#include <QString>
#include <QStringList>
#include <string>
#include <System.h>

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
    TimeSeriesSet<double>::operator=(rhs);
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
void ResultGrid::WriteToVTPSnapShot(const std::string &quanname, const std::string &filename, int i, const double &scale) const
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

void ResultGrid::WriteToVTP(const std::string &quanname, const std::string &filename, const double &scale, int start_id) const
{
    for (unsigned k=0; k<operator[](0).size(); k++)
    {
        cout<<"snapshot: "<<k<<endl;
        std::string name = aquiutils::split(filename,'.')[0]+"_"+aquiutils::numbertostring(k+1+start_id,4)+".vtp";
        cout<<name<<endl;
        WriteToVTPSnapShot(quanname, name, k, scale);
    }
}

vtkSmartPointer<vtkPolyData> ResultGrid::MakeCylinder(
    double radius, double height, double zCenter
    ) {
    vtkNew<vtkCylinderSource> cyl;
    cyl->SetRadius(radius);
    cyl->SetHeight(height);
    cyl->SetResolution(64);
    cyl->CappingOn();

    // Rotate X→Z and move to zCenter
    vtkNew<vtkTransform> tf;
    tf->RotateZ(90);
    tf->Translate(0.0, 0.0, zCenter);

    vtkNew<vtkTransformPolyDataFilter> tfFilter;
    tfFilter->SetTransform(tf);
    tfFilter->SetInputConnection(cyl->GetOutputPort());
    tfFilter->Update();

    vtkNew<vtkTriangleFilter> triang;
    triang->SetInputConnection(tfFilter->GetOutputPort());
    triang->Update();

    vtkNew<vtkCleanPolyData> cleaner;
    cleaner->SetInputConnection(triang->GetOutputPort());
    cleaner->Update();

    return cleaner->GetOutput();
}

vtkSmartPointer<vtkPolyData> ResultGrid::MakeTube(
    double radius, double height, double zCenter,
    std::map<std::string, float> values) {

    // Validate inputs
    if (radius <= 0 || height <= 0) {
        std::cerr << "Invalid tube parameters: radius=" << radius << ", height=" << height << std::endl;
        return nullptr;
    }

    std::cout << "Creating solid cylinder with radius=" << radius
              << ", height=" << height << ", zCenter=" << zCenter << std::endl;

    // Create solid cylinder
    vtkNew<vtkCylinderSource> cylinder;
    cylinder->SetRadius(radius);
    cylinder->SetHeight(height);
    cylinder->SetResolution(32);
    cylinder->CappingOff(); // Open ends for tube/pipe
    cylinder->SetCenter(0, 0, zCenter);

    try {
        cylinder->Update();
    } catch (...) {
        std::cerr << "Failed to create solid cylinder" << std::endl;
        return nullptr;
    }

    vtkPolyData* rawResult = cylinder->GetOutput();
    if (!rawResult || rawResult->GetNumberOfPoints() == 0) {
        std::cerr << "Cylinder source produced no geometry" << std::endl;
        return nullptr;
    }

    // Clean the mesh
    vtkNew<vtkCleanPolyData> cleaner;
    cleaner->SetInputData(rawResult);
    cleaner->Update();

    vtkSmartPointer<vtkPolyData> result = cleaner->GetOutput();

    if (!result || result->GetNumberOfPoints() == 0) {
        std::cerr << "Cleaning resulted in empty mesh" << std::endl;
        return nullptr;
    }

    // Add scalar arrays with validation
    vtkIdType npts = result->GetNumberOfPoints();
    for (const auto& [name, val] : values) {
        // Validate the value
        if (!std::isfinite(val)) {
            std::cerr << "Warning: Non-finite value for " << name << ": " << val << std::endl;
            continue; // Skip this array
        }

        vtkNew<vtkFloatArray> array;
        array->SetName(name.c_str());
        array->SetNumberOfTuples(npts);
        array->FillComponent(0, val);
        result->GetPointData()->AddArray(array);
        std::cout << "Added array '" << name << "' with value " << val << std::endl;
    }

    std::cout << "Created solid cylinder with " << result->GetNumberOfPoints()
              << " points and " << result->GetNumberOfCells() << " cells" << std::endl;

    return result;
}

vtkSmartPointer<vtkPolyData> ResultGrid::MakeTube(
    double outerRadius, double innerRadius, double height, double zCenter,
    std::map<std::string, float> values) {

    // Validate inputs
    if (outerRadius <= 0 || height <= 0) {
        std::cerr << "Invalid tube parameters: outerRadius=" << outerRadius
                  << ", height=" << height << std::endl;
        return nullptr;
    }

    if (innerRadius < 0) {
        std::cerr << "Inner radius cannot be negative: " << innerRadius << std::endl;
        return nullptr;
    }

    if (innerRadius >= outerRadius) {
        std::cerr << "Inner radius (" << innerRadius
                  << ") must be less than outer radius (" << outerRadius << ")" << std::endl;
        return nullptr;
    }

    // Case 1: Solid cylinder (innerRadius == 0 or very small)
    if (innerRadius < 1e-6) {
        std::cout << "Creating solid cylinder with radius=" << outerRadius
                  << ", height=" << height << std::endl;

        vtkNew<vtkCylinderSource> cylinder;
        cylinder->SetRadius(outerRadius);
        cylinder->SetHeight(height);
        cylinder->SetResolution(32);
        cylinder->CappingOff(); // Open ends
        cylinder->SetCenter(0, 0, zCenter);

        try {
            cylinder->Update();
        } catch (...) {
            std::cerr << "Failed to create solid cylinder" << std::endl;
            return nullptr;
        }

        vtkSmartPointer<vtkPolyData> result = cylinder->GetOutput();

        if (!result || result->GetNumberOfPoints() == 0) {
            std::cerr << "Cylinder source produced no geometry" << std::endl;
            return nullptr;
        }

        // Add scalar arrays
        vtkIdType npts = result->GetNumberOfPoints();
        for (const auto& [name, val] : values) {
            if (!std::isfinite(val)) {
                std::cerr << "Warning: Non-finite value for " << name << ": " << val << std::endl;
                continue;
            }

            vtkNew<vtkFloatArray> array;
            array->SetName(name.c_str());
            array->SetNumberOfTuples(npts);
            array->FillComponent(0, val);
            result->GetPointData()->AddArray(array);
        }

        std::cout << "Created solid cylinder with " << result->GetNumberOfPoints()
                  << " points and " << result->GetNumberOfCells() << " cells" << std::endl;

        return result;
    }

    // Case 2: Hollow tube - CREATE MANUALLY (no boolean operations)
    std::cout << "Creating hollow tube manually with outerRadius=" << outerRadius
              << ", innerRadius=" << innerRadius << ", height=" << height << std::endl;

    const int numSides = 32;
    const double halfHeight = height / 2.0;
    const double bottomZ = zCenter - halfHeight;
    const double topZ = zCenter + halfHeight;

    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> polys;

    // Create points for all rings
    // Bottom outer ring (indices 0 to numSides-1)
    for (int i = 0; i < numSides; i++) {
        double angle = 2.0 * M_PI * i / numSides;
        double x = outerRadius * cos(angle);
        double y = outerRadius * sin(angle);
        points->InsertNextPoint(x, y, bottomZ);
    }

    // Bottom inner ring (indices numSides to 2*numSides-1)
    for (int i = 0; i < numSides; i++) {
        double angle = 2.0 * M_PI * i / numSides;
        double x = innerRadius * cos(angle);
        double y = innerRadius * sin(angle);
        points->InsertNextPoint(x, y, bottomZ);
    }

    // Top outer ring (indices 2*numSides to 3*numSides-1)
    for (int i = 0; i < numSides; i++) {
        double angle = 2.0 * M_PI * i / numSides;
        double x = outerRadius * cos(angle);
        double y = outerRadius * sin(angle);
        points->InsertNextPoint(x, y, topZ);
    }

    // Top inner ring (indices 3*numSides to 4*numSides-1)
    for (int i = 0; i < numSides; i++) {
        double angle = 2.0 * M_PI * i / numSides;
        double x = innerRadius * cos(angle);
        double y = innerRadius * sin(angle);
        points->InsertNextPoint(x, y, topZ);
    }

    // Create faces
    for (int i = 0; i < numSides; i++) {
        int next = (i + 1) % numSides;

        // Outer surface (quads)
        vtkNew<vtkQuad> outerQuad;
        outerQuad->GetPointIds()->SetId(0, i);                      // bottom outer
        outerQuad->GetPointIds()->SetId(1, next);                   // bottom outer next
        outerQuad->GetPointIds()->SetId(2, 2 * numSides + next);    // top outer next
        outerQuad->GetPointIds()->SetId(3, 2 * numSides + i);       // top outer
        polys->InsertNextCell(outerQuad);

        // Inner surface (quads, reversed winding for inward normals)
        vtkNew<vtkQuad> innerQuad;
        innerQuad->GetPointIds()->SetId(0, numSides + next);        // bottom inner next
        innerQuad->GetPointIds()->SetId(1, numSides + i);           // bottom inner
        innerQuad->GetPointIds()->SetId(2, 3 * numSides + i);       // top inner
        innerQuad->GetPointIds()->SetId(3, 3 * numSides + next);    // top inner next
        polys->InsertNextCell(innerQuad);
    }

    // Note: We're NOT adding end caps to keep the tube open (like a pipe)
    // If you want closed ends, you would add triangular faces here

    vtkNew<vtkPolyData> result;
    result->SetPoints(points);
    result->SetPolys(polys);

    // Generate normals for proper shading
    vtkNew<vtkPolyDataNormals> normalGenerator;
    normalGenerator->SetInputData(result);
    normalGenerator->ConsistencyOn();
    normalGenerator->SetAutoOrientNormals(true);
    normalGenerator->Update();

    vtkSmartPointer<vtkPolyData> finalResult = normalGenerator->GetOutput();

    if (!finalResult || finalResult->GetNumberOfPoints() == 0) {
        std::cerr << "Manual tube creation failed" << std::endl;
        return nullptr;
    }

    // Add scalar arrays
    vtkIdType npts = finalResult->GetNumberOfPoints();
    for (const auto& [name, val] : values) {
        if (!std::isfinite(val)) {
            std::cerr << "Warning: Non-finite value for " << name << ": " << val << std::endl;
            continue;
        }

        vtkNew<vtkFloatArray> array;
        array->SetName(name.c_str());
        array->SetNumberOfTuples(npts);
        array->FillComponent(0, val);
        finalResult->GetPointData()->AddArray(array);
    }

    std::cout << "Created hollow tube manually with " << finalResult->GetNumberOfPoints()
              << " points and " << finalResult->GetNumberOfCells() << " cells" << std::endl;

    return finalResult;
}


vtkSmartPointer<vtkPolyData> ResultGrid::MakeHollowTube(
    double outerRadius, double innerRadius, double height, double zCenter,
    std::map<std::string, float> values) {

    // Create a line along Z-axis
    vtkNew<vtkLineSource> line;
    line->SetPoint1(0, 0, zCenter - height/2);
    line->SetPoint2(0, 0, zCenter + height/2);
    line->SetResolution(20);
    line->Update();

    // Create OUTER tube
    vtkNew<vtkTubeFilter> outerTube;
    outerTube->SetInputConnection(line->GetOutputPort());
    outerTube->SetRadius(outerRadius);
    outerTube->SetNumberOfSides(64); // Higher resolution for better boolean
    outerTube->CappingOff();
    outerTube->Update();

    // Create INNER tube - make it slightly longer to ensure intersection
    vtkNew<vtkLineSource> innerLine;
    innerLine->SetPoint1(0, 0, zCenter - height/2 - 0.1); // Extend slightly
    innerLine->SetPoint2(0, 0, zCenter + height/2 + 0.1); // Extend slightly
    innerLine->SetResolution(20);
    innerLine->Update();

    vtkNew<vtkTubeFilter> innerTube;
    innerTube->SetInputConnection(innerLine->GetOutputPort());
    innerTube->SetRadius(innerRadius);
    innerTube->SetNumberOfSides(64); // Same resolution
    innerTube->CappingOff();
    innerTube->Update();

    // Clean the meshes before boolean operation
    vtkNew<vtkCleanPolyData> cleanOuter;
    cleanOuter->SetInputData(outerTube->GetOutput());
    cleanOuter->Update();

    vtkNew<vtkCleanPolyData> cleanInner;
    cleanInner->SetInputData(innerTube->GetOutput());
    cleanInner->Update();

    // Boolean difference with better tolerance
    vtkNew<vtkBooleanOperationPolyDataFilter> booleanFilter;
    booleanFilter->SetOperationToDifference();
    booleanFilter->SetInputData(0, cleanOuter->GetOutput());
    booleanFilter->SetInputData(1, cleanInner->GetOutput());
    booleanFilter->SetTolerance(1e-3); // Much larger tolerance
    booleanFilter->Update();

    vtkSmartPointer<vtkPolyData> result = booleanFilter->GetOutput();

    // Check if boolean operation worked
    if (!result || result->GetNumberOfPoints() == 0) {
        std::cerr << "Boolean operation failed, returning solid tube" << std::endl;
        result = outerTube->GetOutput(); // Fallback to solid tube
    }

    // Add scalar arrays
    vtkIdType npts = result->GetNumberOfPoints();
    for (const auto& [name, val] : values) {
        vtkNew<vtkFloatArray> array;
        array->SetName(name.c_str());
        array->SetNumberOfTuples(npts);
        array->FillComponent(0, val);
        result->GetPointData()->AddArray(array);
    }

    return result;
}

vtkSmartPointer<vtkPolyData> ResultGrid::MakeHollowTubeManual(
    double outerRadius, double innerRadius, double height, double zCenter,
    std::map<std::string, float> values) {

    const int numSides = 32;
    const double halfHeight = height / 2.0;

    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> polys;

    // Create points: bottom outer, bottom inner, top outer, top inner
    // Bottom outer ring
    for (int i = 0; i < numSides; i++) {
        double angle = 2.0 * M_PI * i / numSides;
        double x = outerRadius * cos(angle);
        double y = outerRadius * sin(angle);
        points->InsertNextPoint(x, y, zCenter - halfHeight);
    }

    // Bottom inner ring
    for (int i = 0; i < numSides; i++) {
        double angle = 2.0 * M_PI * i / numSides;
        double x = innerRadius * cos(angle);
        double y = innerRadius * sin(angle);
        points->InsertNextPoint(x, y, zCenter - halfHeight);
    }

    // Top outer ring
    for (int i = 0; i < numSides; i++) {
        double angle = 2.0 * M_PI * i / numSides;
        double x = outerRadius * cos(angle);
        double y = outerRadius * sin(angle);
        points->InsertNextPoint(x, y, zCenter + halfHeight);
    }

    // Top inner ring
    for (int i = 0; i < numSides; i++) {
        double angle = 2.0 * M_PI * i / numSides;
        double x = innerRadius * cos(angle);
        double y = innerRadius * sin(angle);
        points->InsertNextPoint(x, y, zCenter + halfHeight);
    }

    // Create faces
    // Outer surface quads
    for (int i = 0; i < numSides; i++) {
        int next = (i + 1) % numSides;

        vtkNew<vtkQuad> quad;
        quad->GetPointIds()->SetId(0, i);                    // bottom outer
        quad->GetPointIds()->SetId(1, next);                 // bottom outer next
        quad->GetPointIds()->SetId(2, 2 * numSides + next);  // top outer next
        quad->GetPointIds()->SetId(3, 2 * numSides + i);     // top outer
        polys->InsertNextCell(quad);
    }

    // Inner surface quads (reversed winding for inward-facing normals)
    for (int i = 0; i < numSides; i++) {
        int next = (i + 1) % numSides;

        vtkNew<vtkQuad> quad;
        quad->GetPointIds()->SetId(0, numSides + next);      // bottom inner next
        quad->GetPointIds()->SetId(1, numSides + i);         // bottom inner
        quad->GetPointIds()->SetId(2, 3 * numSides + i);     // top inner
        quad->GetPointIds()->SetId(3, 3 * numSides + next);  // top inner next
        polys->InsertNextCell(quad);
    }

    vtkNew<vtkPolyData> result;
    result->SetPoints(points);
    result->SetPolys(polys);

    // Add scalar arrays
    vtkIdType npts = result->GetNumberOfPoints();
    for (const auto& [name, val] : values) {
        vtkNew<vtkFloatArray> array;
        array->SetName(name.c_str());
        array->SetNumberOfTuples(npts);
        array->FillComponent(0, val);
        result->GetPointData()->AddArray(array);
    }

    return result;
}

void ResultGrid::Make3DVTK(vector<string> quantity, double dr, System *system, string filename) {
    if (!system) {
        std::cerr << "System is null!" << std::endl;
        return;
    }

    vtkNew<vtkMultiBlockDataSet> mb;
    int validBlockCount = 0;

    for (int i = 0; i < system->BlockCount(); i++) {
        string block_name = system->block(i)->GetName();

        if (!system->block(block_name) || system->block(block_name)->GetType()!="Soil") {
            std::cout << "Skipping non-Soil block: " << block_name << std::endl;
            continue;
        }

        std::cout << "Processing Block: " << block_name << std::endl;

        try {
            double ptx = system->block(block_name)->GetVal("act_X");
            double pty = system->block(block_name)->GetVal("act_Y");
            double depth = system->block(block_name)->GetVal("depth");

            std::cout << "  ptx=" << ptx << ", pty=" << pty << ", depth=" << depth << ", dr=" << dr << std::endl;

            // Validate parameters
            if (!std::isfinite(ptx) || !std::isfinite(pty) || !std::isfinite(depth) || depth <= 0) {
                std::cerr << "  Invalid geometry parameters, skipping block" << std::endl;
                continue;
            }

            // Calculate radii
            double outerRadius = std::abs(ptx) + dr/2.0;
            double innerRadius = std::abs(ptx) - dr/2.0;

            // Ensure inner radius is not negative - if it would be, make it zero for solid cylinder
            if (innerRadius < 0) {
                std::cout << "  Inner radius would be negative (" << innerRadius
                          << "), setting to 0 for solid cylinder" << std::endl;
                innerRadius = 0;
            }

            // Ensure outer radius is positive
            if (outerRadius <= 0) {
                std::cerr << "  Invalid outer radius (" << outerRadius << "), skipping block" << std::endl;
                continue;
            }

            std::cout << "  outerRadius=" << outerRadius << ", innerRadius=" << innerRadius << std::endl;

            // Collect quantity values
            map<string, float> quantityvaluemap;
            for (int j = 0; j < quantity.size(); j++) {
                if (system->block(i)->HasQuantity(quantity[j])) {
                    float val = system->block(block_name)->GetVal(quantity[j], Expression::timing::past);
                    if (std::isfinite(val)) {
                        quantityvaluemap[quantity[j]] = val;
                        std::cout << "    " << quantity[j] << " = " << val << std::endl;
                    } else {
                        std::cerr << "    Warning: Non-finite value for " << quantity[j] << std::endl;
                    }
                }
            }

            // Create tube (solid cylinder if innerRadius=0, hollow tube otherwise)
            auto tube = MakeTube(outerRadius, innerRadius, depth, pty, quantityvaluemap);

            if (tube && tube->GetNumberOfPoints() > 0) {
                mb->SetBlock(validBlockCount, tube);
                mb->GetMetaData(validBlockCount)->Set(vtkCompositeDataSet::NAME(), block_name.c_str());
                validBlockCount++;
                std::cout << "  Successfully added block " << validBlockCount << std::endl;
            } else {
                std::cerr << "  Failed to create tube for block: " << block_name << std::endl;
            }

        } catch (const std::exception& e) {
            std::cerr << "  Exception processing block " << block_name << ": " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "  Unknown exception processing block " << block_name << std::endl;
        }
    }

    if (validBlockCount == 0) {
        std::cerr << "No valid blocks to write!" << std::endl;
        return;
    }

    std::cout << "Writing " << validBlockCount << " blocks to " << filename << std::endl;

    try {
        vtkNew<vtkXMLMultiBlockDataWriter> writer;
        writer->SetFileName(filename.c_str());
        writer->SetInputData(mb);
        writer->SetDataModeToAscii(); // Use ASCII mode for better debugging
        writer->Write();

        std::cout << "Successfully saved " << validBlockCount << " blocks to " << filename << std::endl;
    } catch (...) {
        std::cerr << "Failed to write VTK file: " << filename << std::endl;
    }
}
#endif
