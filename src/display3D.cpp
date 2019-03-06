#include <DGtal/base/Common.h>
#include <DGtal/io/Display3D.h>
#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/viewers/Viewer3D.h>

#define LOG(X) std::cout << X << std::endl

using namespace DGtal;

// TODO: Bounding box
// TODO: Gaussian voxelization
// TODO: Raytracing / Raymarching
// TODO: Equations plane & lines
// TODO: Intersections
// NOTE: Vector vertices / faces & display
int main(int argc, char **argv)
{
    int voxelizationXRes;
    int voxelizationYRes;
    int voxelizationZRes;

    QApplication application(argc, argv);
    Viewer3D<> viewer;
    viewer.show();
    std::string inputFile = "../Mesh/bunny.off";

    if (argc > 1)
    {
        inputFile = argv[1];
    }

    if (argc == 5)
    {
        voxelizationXRes = std::atoi(argv[2]);
        voxelizationYRes = std::atoi(argv[3]);
        voxelizationZRes = std::atoi(argv[4]);
    }

    // Since the input points are not necessary integers we use the PointD3D from Display3D.

    Mesh<Viewer3D<>::RealPoint> mesh;

    Mesh<Viewer3D<>::RealPoint>::Iterator itr;

    mesh << inputFile;
    trace.info() << "Importing done..." << std::endl;

    LOG("Begin: " << mesh.vertexBegin().base());
    LOG("End: " << mesh.vertexEnd().base());

    LOG(mesh.getVertex(0)[0]);
    mesh.getFace(0);

    /* for (itr = mesh.vertexBegin(); itr < mesh.vertexEnd(); itr++)
    {
        trace.info() << *itr << std::endl;
    } */

    // viewer.setLineColor(DGtal::Color(150, 0, 0, 254));
    viewer << mesh;
    viewer << Viewer3D<>::updateDisplay;
    viewer >> "/home/raitoning/bunny.off";

    return application.exec();
}