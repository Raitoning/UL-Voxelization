// Axes:
// X: Horizontal
// Y: Vertical (Y up)
// Z: Forward/Backward (-Z forward)

#include <string>
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

    // Map containing all the argments.
    std::map<std::string, std::string> argments;
    std::map<std::string, std::string>::iterator argsIterator;

    // Define if the mesh must be normalized.
    bool normalize = false;

    // Resolution on the X axis.
    int horizontalResolution;

    // Resolution on the Y axis.
    int verticalResolution;

    // Resolution on the Z axis.
    int forwardResolution;

    // Scale factor used in normalization
    double requiredScale = 1.;

    vector<Z3i::RealPoint> intersectionPoints;

    // Checking the arguments.
    // If there isn't a file name as argument, stop the execution.
    if (argc < 2)
    {
        DGtal::trace.info() << "Usage: display3D file [hres] [vres] [zres] [--normalize:size]" << std::endl
                            << "Now exiting..." << std::endl;

        return 0;
    }

    // If there are at least 2 argments (0 = program name, 1 = file name)
    // Create a map to store each argments for easy use afterward.
    for (int i = 2; i < argc; i++)
    {
        std::stringstream stream(argv[i]);
        std::string value;
        std::vector<std::string> strings;

        // If it's an argments with a value (eg. --normalize=5)
        // Split it in two strings: "--normalize" and "5" and put it in the map

        // TODO: Find a way to split either on '=' or ':'
        while (std::getline(stream, value, '='))
        {
            strings.push_back(value);
        }

        // if the split strings size is 1, it's not an argment with a value
        // Simply put it as <argment, argment> in the map.
        if (strings.size() == 1)
        {
            argments[argv[i]] = argv[i];
        }
        // Else, put the name of the argment and the value as <name, value>
        // In the map (eg. <"--normalize", "5">).
        else
        {
            argments[strings[0]] = strings[1];
        }
    }

    // Get the file name to import the mesh.
    inputFile = argv[1];
    DGtal::trace.info() << "Input file: " << inputFile << std::endl;

    argsIterator = argments.find("--normalize");

    if (argsIterator != argments.end())
    {
        normalize = true;

        std::string value = argments["--normalize"];

        requiredScale = atof(value.c_str());
    }

    argsIterator = argments.find("--resolution");

    if (argsIterator != argments.end())
    {
        std::stringstream stream(argsIterator->second);
        std::string value;
        std::vector<std::string> strings;

        while (std::getline(stream, value, ' '))
        {
            strings.push_back(value);
        }

        horizontalResolution = atoi(strings[0].c_str());
        verticalResolution = atoi(strings[1].c_str());
        forwardResolution = atoi(strings[2].c_str());

        DGtal::trace.info() << "X, Y, Z:" << horizontalResolution << " " << verticalResolution << " " << forwardResolution << std::endl;

        if (horizontalResolution < 1 || verticalResolution < 1 || forwardResolution < 1)
        {
            DGtal::trace.error() << "A resolution must be at least of two." << std::endl;
            DGtal::trace.info() << "Now exiting..." << std::endl;

            return 0;
        }
        gaussian = false;
    }

    // qT Application hosting the viewer.
    QApplication application(argc, argv);

    // 3D Viewer.
    Viewer3D<> viewer;
    viewer.show();

    // Since the input points are not necessary integers we use the PointD3D from Display3D.
    Mesh<Viewer3D<>::RealPoint> mesh;

    // Importing the file
    DGtal::trace.info() << "Importing..." << std::endl;

    mesh << inputFile;

    DGtal::trace.info() << "Importing done..." << std::endl;
    DGtal::trace.info() << "Number of vertices: " << mesh.nbVertex() << std::endl;
    DGtal::trace.info() << "Number of faces: " << mesh.nbFaces() << std::endl;

    // Getting the bounding box.
    std::pair<Viewer3D<>::RealPoint, Viewer3D<>::RealPoint>
        boundingBox = mesh.getBoundingBox();

    double xSize = abs(boundingBox.second[0] - boundingBox.first[0]);
    double ySize = abs(boundingBox.second[1] - boundingBox.first[1]);
    double zSize = abs(boundingBox.second[2] - boundingBox.first[2]);

    DGtal::trace.info() << "Mesh size: " << xSize << "; " << ySize << "; " << zSize << std::endl;

    if (normalize)
    {

        double scaleFactor = requiredScale / (std::min(xSize, std::min(ySize, zSize)));
        DGtal::trace.info() << "Scale factor: " << scaleFactor << std::endl;

        // Change the scale of the mesh if it's too small.
        mesh.changeScale(scaleFactor);
    }

    boundingBox = mesh.getBoundingBox();

    xSize = abs(boundingBox.second[0] - boundingBox.first[0]);
    ySize = abs(boundingBox.second[1] - boundingBox.first[1]);
    zSize = abs(boundingBox.second[2] - boundingBox.first[2]);

    DGtal::trace.info() << "Mesh size: " << xSize << "; " << ySize << "; " << zSize << std::endl;

    DisplayBoundingBox(viewer, boundingBox.first, boundingBox.second);

    DGtal::trace.info() << "Bounding box: " << std::endl
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
                        DGtal::trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
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
                        DGtal::trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
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
                        DGtal::trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
                    }
                }
            }
        }
    }
    else
    {
        // Gaussian voxelization
        // TODO: Do the actual voxelization after unit tests
        DGtal::trace.warning() << "Gaussian voxelization not implemented yet. Doing nothing." << std::endl;
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
