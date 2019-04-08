// Axes:
// X: Horizontal
// Y: Vertical (Y up)
// Z: Forward/Backward (-Z forward)

#include <DGtal/base/Common.h>
#include <DGtal/io/Display3D.h>
#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/viewers/Viewer3D.h>
#include "Algorithms.h"

#define LOG(X) std::cout << X << std::endl

using namespace DGtal;

// TODO: Gaussian voxelization
// TODO: Equations plane & lines
// TODO: Test -> Direction

int main(int argc, char **argv)
{
    // File name of the mesh to import.
    std::string inputFile;

    // Define if the voxelisation will be Gaussian or resolution based.
    bool gaussian = true;

    // Resolution on the X axis.
    int horizontalResolution;

    // Resolution on the Y axis.
    int verticalResolution;

    // Resolution on the Z axis.
    int forwardResolution;

    vector<Z3i::RealPoint> intersectionPoints;

    // Checking the arguments.
    // If there isn't a file name as argument, stop the execution.
    if (argc != 2 && argc != 5)
    {
        trace.info() << "Usage: display3D file [hres] [vres] [zres]" << std::endl
                     << "Now exiting..." << std::endl;

        return 0;
    }
    else
    {
        // Get the file name to import the mesh.
        inputFile = argv[1];
        trace.info() << "Input file: " << inputFile << std::endl;
    }

    if (argc == 5)
    {
        horizontalResolution = atoi(argv[2]);
        verticalResolution = atoi(argv[3]);
        forwardResolution = atoi(argv[4]);
        gaussian = false;

        LOG("X, Y, Z:" << horizontalResolution << " " << verticalResolution << " " << forwardResolution);
    }

    // qT Application hosting the viewer.
    QApplication application(argc, argv);

    // 3D Viewer.
    Viewer3D<> viewer;
    viewer.show();

    // Since the input points are not necessary integers we use the PointD3D from Display3D.
    Mesh<Viewer3D<>::RealPoint> mesh;

    // Importing the file
    trace.info() << "Importing..." << std::endl;

    mesh << inputFile;

    trace.info() << "Importing done..." << std::endl;
    trace.info() << "Number of vertices: " << mesh.nbVertex() << std::endl;
    trace.info() << "Number of faces: " << mesh.nbFaces() << std::endl;

    // Getting the bounding box.
    std::pair<Viewer3D<>::RealPoint, Viewer3D<>::RealPoint>
        boundingBox = mesh.getBoundingBox();

    double xSize = abs(boundingBox.second[0] - boundingBox.first[0]);
    double ySize = abs(boundingBox.second[1] - boundingBox.first[1]);
    double zSize = abs(boundingBox.second[2] - boundingBox.first[2]);

    trace.info() << "Mesh size: " << xSize << "; " << ySize << "; " << zSize << std::endl;

    double scaleFactor = 1.;

    scaleFactor = 1. / (std::min(xSize, std::min(ySize, zSize)));

    trace.info() << "Scale factor: " << scaleFactor << std::endl;

    // Change the scale of the mesh if it's too small.
    mesh.changeScale(scaleFactor);

    boundingBox = mesh.getBoundingBox();

    xSize = abs(boundingBox.second[0] - boundingBox.first[0]);
    ySize = abs(boundingBox.second[1] - boundingBox.first[1]);
    zSize = abs(boundingBox.second[2] - boundingBox.first[2]);

    trace.info() << "Mesh size: " << xSize << "; " << ySize << "; " << zSize << std::endl;

    DisplayBoundingBox(viewer, boundingBox.first, boundingBox.second);

    trace.info() << "Bounding box: " << std::endl
                 << boundingBox.first << std::endl
                 << boundingBox.second << std::endl;

    Z3i::RealPoint intersection;

    double xStep = ((boundingBox.second[0] - boundingBox.first[0]) / double(horizontalResolution - 1));
    double yStep = ((boundingBox.second[1] - boundingBox.first[1]) / double(verticalResolution - 1));
    double zStep = ((boundingBox.second[2] - boundingBox.first[2]) / double(forwardResolution - 1));

    if (!gaussian)
    {
        // Raytracing from under.
        Z3i::RealPoint rayDirection = Z3i::RealPoint(0, 1, 0);

        LOG("Min bbox:" << boundingBox.first[0] << " " << boundingBox.first[1] << " " << boundingBox.first[2]);
        LOG("Max bbox:" << boundingBox.second[0] << " " << boundingBox.second[1] << " " << boundingBox.second[2]);

        for (double x = boundingBox.first[0]; x <= boundingBox.second[0]; x += xStep)
        {
            for (double z = boundingBox.first[2]; z <= boundingBox.second[2]; z += zStep)
            {
                Z3i::RealPoint rayOrigin(x, boundingBox.first[1] - 1, z);

                intersectionPoints.push_back(rayOrigin);
                intersectionPoints.push_back(rayOrigin + rayDirection * 3);
                // Test if the test ray can intersect anything.
                for (int i = 0; i < mesh.nbFaces(); i++)
                {
                    // If a face is intersected, set it's color to red.
                    if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                    {
                        mesh.setFaceColor(i, Color(255, 0, 0));
                        trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
                    }
                }
            }
        }

        // Raytracing from in front of.
        rayDirection = Z3i::RealPoint(0, 0, -1);

        for (double x = boundingBox.first[0]; x <= boundingBox.second[0]; x += xStep)
        {
            for (double y = boundingBox.first[1]; y <= boundingBox.second[1]; y += yStep)
            {
                Z3i::RealPoint rayOrigin(x, y, boundingBox.first[2] + 1);

                intersectionPoints.push_back(rayOrigin);
                intersectionPoints.push_back(rayOrigin + rayDirection * 3);
                // Test if the test ray can intersect anything.
                for (int i = 0; i < mesh.nbFaces(); i++)
                {
                    // If a face is intersected, set it's color to red.
                    if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                    {
                        // intersectionPoints.push_back(rayOrigin);
                        // intersectionPoints.push_back(rayOrigin + rayDirection);
                        mesh.setFaceColor(i, Color(255, 0, 0));
                        trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
                    }
                }
            }
        }

        // Raytracing from the left.
        rayDirection = Z3i::RealPoint(1, 0, 0);

        for (double y = boundingBox.first[1]; y <= boundingBox.second[1]; y += yStep)
        {
            for (double z = boundingBox.first[2]; z <= boundingBox.second[2]; z += zStep)
            {
                Z3i::RealPoint rayOrigin(boundingBox.first[0] - 1, y, z);

                intersectionPoints.push_back(rayOrigin);
                intersectionPoints.push_back(rayOrigin + rayDirection * 3);
                // Test if the test ray can intersect anything.
                for (int i = 0; i < mesh.nbFaces(); i++)
                {
                    // If a face is intersected, set it's color to red.
                    if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                    {
                        mesh.setFaceColor(i, Color(255, 0, 0));
                        trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
                    }
                }
            }
        }
    }
    else
    {
        // Gaussian voxelization
        // TODO: Do the actual voxelization after unit tests
        trace.warning() << "Gaussian voxelization not implemented yet. Doing nothing." << std::endl;
    }

    // Push the mesh into the viewer.
    viewer << mesh;

    intersectionPoints.push_back(boundingBox.first);
    intersectionPoints.push_back(boundingBox.second);

    for (int i = 0; i < intersectionPoints.size(); i += 2)
    {
        viewer.addLine(intersectionPoints[i], intersectionPoints[i + 1], 0.03);
    }

    viewer << Viewer3D<>::updateDisplay;

    // Return the qT application.
    return application.exec();
}
