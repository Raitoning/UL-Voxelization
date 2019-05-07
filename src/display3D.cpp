// Axes:
// X: Horizontal
// Y: Vertical (Y up)
// Z: Forward/Backward (-Z forward)

#include <string>
#include <DGtal/base/Common.h>
#include <DGtal/io/Display3D.h>
#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/viewers/Viewer3D.h>

#include <omp.h>
#include <chrono>

#include "Algorithms.h"

#define LOG(X) std::cout << X << std::endl
#define epsilon 1e-10

using namespace DGtal;

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

    // RealPoint for the non orthogonal plane origin
    Z3i::RealPoint planeOrigin;

    // RealPoint for the non orthogonal plane direction
    Z3i::RealPoint planeDirection;

    // Delta for the non orthogonal plane
    double planeDelta = 1.0;

    vector<Z3i::RealPoint> intersectionPoints;

    std::vector<std::vector<Z3i::RealPoint>> pointInterieurs;

    // Checking the arguments.
    // If there isn't a file name as argument, stop the execution.
    if (argc < 2)
    {
        DGtal::trace.info() << "Usage: display3D file" << std::endl
                            << "Options:" << std::endl
                            << "[--resolution=\"x y z\"]\t\tDefine the resolution used to voxelize" << std::endl
                            << "[--normalize=size]\t\tDefine the minimal size of the mesh" << std::endl
                            << "[--threaded=# of threads]\tDefine the number of threads to use. Default is the number of threads of the CPU." << std::endl
                            << "[--origin=\"x y z\"]\t\tDefine the origin of a non orthogonal plane." << std::endl
                            << "[--direction=\"x y z\"]\t\tDefine the direction of a non orthogonal plane." << std::endl
                            << "[--delta=value]\t\t\tDefine the delta of a non orthogonal plane." << std::endl
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

    // Handling the normalization argment.
    argsIterator = argments.find("--normalize");

    if (argsIterator != argments.end())
    {
        normalize = true;

        std::string value = argments["--normalize"];

        requiredScale = atof(value.c_str());
    }

    // Handling the resolution argment.
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

    // Handling the threading argment.
    argsIterator = argments.find("--threaded");

    if (argsIterator != argments.end())
    {
        std::string value = argments["--threaded"];
        int nbThreads = atoi(value.c_str());

        if (nbThreads < 1)
        {
            DGtal::trace.error() << "At least one thread must be used." << std::endl;
            DGtal::trace.info() << "Now exiting..." << std::endl;

            return 0;
        }

        omp_set_num_threads(nbThreads);
    }

    // Handling the plane origin argment.
    argsIterator = argments.find("--origin");

    if (argsIterator != argments.end())
    {
        std::stringstream stream(argsIterator->second);
        std::string value;
        std::vector<std::string> strings;

        while (std::getline(stream, value, ' '))
        {
            strings.push_back(value);
        }

        double x = atof(strings[0].c_str());
        double y = atof(strings[1].c_str());
        double z = atof(strings[2].c_str());

        DGtal::trace.info() << "X, Y, Z:" << x << " " << y << " " << z << std::endl;

        planeOrigin = Z3i::RealPoint(x, y, z);

        DGtal::trace.info() << "Plane origin: " << planeOrigin << std::endl;
    }

    // Handling the plane direction argment.
    argsIterator = argments.find("--direction");

    if (argsIterator != argments.end())
    {
        std::stringstream stream(argsIterator->second);
        std::string value;
        std::vector<std::string> strings;

        while (std::getline(stream, value, ' '))
        {
            strings.push_back(value);
        }

        double x = atof(strings[0].c_str());
        double y = atof(strings[1].c_str());
        double z = atof(strings[2].c_str());

        DGtal::trace.info() << "X, Y, Z:" << x << " " << y << " " << z << std::endl;

        planeDirection = Z3i::RealPoint(x, y, z);

        DGtal::trace.info() << "Plane direction: " << planeDirection << std::endl;
    }

    // Handling the delta argment.
    argsIterator = argments.find("--delta");

    if (argsIterator != argments.end())
    {
        std::string value = argments["--delta"];
        planeDelta = atof(value.c_str());

        DGtal::trace.info() << "Plane delta: " << planeDelta << std::endl;
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

    auto startTime = std::chrono::high_resolution_clock::now();

    if (!gaussian)
    {
        // Raytracing from under.
        Z3i::RealPoint rayDirection = Z3i::RealPoint(0, 1, 0);

        LOG("Min bbox:" << boundingBox.first[0] << " " << boundingBox.first[1] << " " << boundingBox.first[2]);
        LOG("Max bbox:" << boundingBox.second[0] << " " << boundingBox.second[1] << " " << boundingBox.second[2]);

        for (double x = boundingBox.first[0]; x <= boundingBox.second[0] + epsilon; x += xStep)
        {
            for (double z = boundingBox.first[2]; z <= boundingBox.second[2] + epsilon; z += zStep)
            {
                Z3i::RealPoint rayOrigin(x, boundingBox.first[1] - 1, z);

                intersectionPoints.push_back(rayOrigin);
                intersectionPoints.push_back(rayOrigin + rayDirection * 3);

                // Test if the test ray can intersect anything.
#pragma omp parallel for ordered
                for (uint i = 0; i < mesh.nbFaces(); i++)
                {
                    // If a face is intersected, set it's color to red.
                    if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                    {
                        mesh.setFaceColor(i, Color(255, 0, 0));
                        DGtal::trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")"
                                            << "(Thread #" << omp_get_thread_num() << ")" << std::endl;
                    }
                }
            }
        }

        // Raytracing from in front of.
        rayDirection = Z3i::RealPoint(0, 0, -1);

        for (double x = boundingBox.first[0]; x <= boundingBox.second[0] + epsilon; x += xStep)
        {
            for (double y = boundingBox.first[1]; y <= boundingBox.second[1] + epsilon; y += yStep)
            {
                Z3i::RealPoint rayOrigin(x, y, boundingBox.first[2] + 1);

                intersectionPoints.push_back(rayOrigin);
                intersectionPoints.push_back(rayOrigin + rayDirection * 3);
                // Test if the test ray can intersect anything.
#pragma omp parallel for ordered
                for (uint i = 0; i < mesh.nbFaces(); i++)
                {
                    // If a face is intersected, set it's color to red.
                    if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                    {
                        // intersectionPoints.push_back(rayOrigin);
                        // intersectionPoints.push_back(rayOrigin + rayDirection);
                        mesh.setFaceColor(i, Color(255, 0, 0));
                        DGtal::trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")"
                                            << "(Thread #" << omp_get_thread_num() << ")" << std::endl;
                    }
                }
            }
        }

        // Raytracing from the left.
        rayDirection = Z3i::RealPoint(1, 0, 0);

        for (double y = boundingBox.first[1]; y <= boundingBox.second[1] + epsilon; y += yStep)
        {
            for (double z = boundingBox.first[2]; z <= boundingBox.second[2] + epsilon; z += zStep)
            {
                Z3i::RealPoint rayOrigin(boundingBox.first[0] - 1, y, z);

                intersectionPoints.push_back(rayOrigin);
                intersectionPoints.push_back(rayOrigin + rayDirection * 3);
                // Test if the test ray can intersect anything.
#pragma omp parallel for ordered
                for (uint i = 0; i < mesh.nbFaces(); i++)
                {
                    // If a face is intersected, set it's color to red.
                    if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                    {
                        mesh.setFaceColor(i, Color(255, 0, 0));
                        DGtal::trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")"
                                            << "(Thread #" << omp_get_thread_num() << ")" << std::endl;
                    }
                }
            }
        }
    }
    else
    {
        // Gaussian voxelization
        // Raytracing from under.
        Z3i::RealPoint rayDirection = Z3i::RealPoint(0, 1, 0);

        Z3i::RealPoint stepInterieur = createStep(rayDirection, 1, 1, 1);
        std::vector<Z3i::RealPoint> intersectionsVecteur;

        for (int x = boundingBox.first[0] - 1; x <= boundingBox.second[0] + 1; x++)
        {

            for (int z = boundingBox.first[2] - 1; z <= boundingBox.second[2] + 1; z++)
            {
                Z3i::RealPoint rayOrigin(x, boundingBox.first[1] - 1, z);

                intersectionPoints.push_back(rayOrigin);
                intersectionPoints.push_back(rayOrigin + rayDirection * (boundingBox.second[1] - boundingBox.first[1] + 2));

                // Test if the test ray can intersect anything.
#pragma omp parallel for ordered
                for (uint i = 0; i < mesh.nbFaces(); i++)
                {
                    // If a face is intersected, set it's color to red.
                    if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                    {
                        mesh.setFaceColor(i, Color(255, 0, 0));
                        DGtal::trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")"
                                            << "(Thread #" << omp_get_thread_num() << ")" << std::endl;

                        //on recupère les points d'intersection d'un vecteur
#pragma omp ordered
                        {
                            intersectionsVecteur.push_back(intersection);
                        }
                    }
                }

                //calcul des point a l'interieur
                pointInterieurs.push_back(pointInterieur(rayOrigin, rayDirection, intersectionsVecteur, stepInterieur));
                intersectionsVecteur.clear();
            }
        }

        // Raytracing from in front of.
        rayDirection = Z3i::RealPoint(0, 0, 1);

        stepInterieur = createStep(rayDirection, 1, 1, 1);

        for (int x = boundingBox.first[0] - 1; x <= boundingBox.second[0] + 1; x++)
        {
            for (int y = boundingBox.first[1] - 1; y <= boundingBox.second[1] + 1; y++)
            {
                Z3i::RealPoint rayOrigin(x, y, boundingBox.first[2] - 1);

                intersectionPoints.push_back(rayOrigin);
                intersectionPoints.push_back(rayOrigin + rayDirection * (boundingBox.second[2] - boundingBox.first[2] + 2));

                // Test if the test ray can intersect anything.
#pragma omp parallel for ordered
                for (uint i = 0; i < mesh.nbFaces(); i++)
                {
                    // If a face is intersected, set it's color to red.
                    if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                    {
                        mesh.setFaceColor(i, Color(255, 0, 0));
                        DGtal::trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")"
                                            << "(Thread #" << omp_get_thread_num() << ")" << std::endl;

                        //on recupère les points d'intersection d'un vecteur
#pragma omp ordered
                        {
                            intersectionsVecteur.push_back(intersection);
                        }
                    }
                }

                //calcul des point a l'interieur
                pointInterieurs.push_back(pointInterieur(rayOrigin, rayDirection, intersectionsVecteur, stepInterieur));
                intersectionsVecteur.clear();
            }
        }

        // Raytracing from the left.
        rayDirection = Z3i::RealPoint(1, 0, 0);

        stepInterieur = createStep(rayDirection, 1, 1, 1);

        for (int y = boundingBox.first[1] - 1; y <= boundingBox.second[1] + 1; y++)
        {

            for (int z = int(boundingBox.first[2] - 1); z <= int(boundingBox.second[2] + 1); z++)
            {
                Z3i::RealPoint rayOrigin(boundingBox.first[0] - 1, y, z);

                intersectionPoints.push_back(rayOrigin);
                intersectionPoints.push_back(rayOrigin + rayDirection * (boundingBox.second[0] - boundingBox.first[0] + 2));

                // Test if the test ray can intersect anything.

#pragma omp parallel for ordered
                for (uint i = 0; i < mesh.nbFaces(); i++)
                {
                    // If a face is intersected, set it's color to red.
                    if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                    {
                        mesh.setFaceColor(i, Color(255, 0, 0));
                        DGtal::trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")"
                                            << "(Thread #" << omp_get_thread_num() << ")" << std::endl;

                        //on recupère les points d'intersection d'un vecteur
#pragma omp ordered
                        {
                            intersectionsVecteur.push_back(intersection);
                        }
                    }
                }

                //calcul des point a l'interieur
                pointInterieurs.push_back(pointInterieur(rayOrigin, rayDirection, intersectionsVecteur, stepInterieur));
                intersectionsVecteur.clear();
            }
        }
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    //intersectionPoints.push_back(boundingBox.first);
    //intersectionPoints.push_back(boundingBox.second);
    // Viewer3D<>::RealPoint d = boundingBox.second - boundingBox.first;
    // vector<Viewer3D<>::RealPoint> testingVectors = originPoints(boundingBox.first, d, boundingBox.first, boundingBox.second);

    // for (uint k = 0; k < testingVectors.size(); k++)
    // {
    //     for (uint i = 0; i < mesh.nbFaces(); i++)
    //     {
    //         //If a face is intersected, set it's color to red.
    //         if (RayIntersectsTriangle(testingVectors[k], d, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
    //         {
    //             mesh.setFaceColor(i, Color(255, 0, 0));
    //             trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
    //         }
    //     }
    // }

    // Push the mesh into the viewer.
    viewer << mesh;

    for (uint i = 0; i < intersectionPoints.size(); i += 2)
    {
        viewer.addLine(intersectionPoints[i], intersectionPoints[i + 1], 0.03);
    }

    //gestion des points interieur/voxels
    int count = 0;
    int seuil = 0;

    int tabX = (int)abs(boundingBox.first[0] - boundingBox.second[0]) + 1;
    int tabY = (int)abs(boundingBox.first[1] - boundingBox.second[1]) + 1;
    int tabZ = (int)abs(boundingBox.first[2] - boundingBox.second[2]) + 1;

    int voxels[tabX][tabY][tabZ] = {};
    //print des points a l'interieur
    for (int g = 0; g < pointInterieurs.size(); g++)
    {
        for (int h = 0; h < pointInterieurs[g].size(); h++)
        {
            viewer.addCube(pointInterieurs[g][h]);
            //on stock au coordonné g,h,X les points (et leurs nombre d'occurences)
            voxels[(int)round(pointInterieurs[g][h][0] - boundingBox.first[0])][(int)round(pointInterieurs[g][h][1] - boundingBox.first[1])][(int)round(pointInterieurs[g][h][2] - boundingBox.first[2])] += 1;
            count++;
        }
    }

    //conservationSurface(voxels,tabX,tabY,tabZ,seuil);
    //todo retrouver offset
    for (int i = 0; i < tabX; i++)
    {
        for (int j = 0; j < tabY; j++)
        {
            for (int k = 0; k < tabZ; k++)
            {
                if (voxels[i][j][k] > seuil)
                    viewer.addCube(Z3i::RealPoint(i, j, k));
            }
        }
    }

    std::cout
        << "Computation time: " << duration.count() << " milliseconds with " << omp_get_max_threads() << " threads." << std::endl;

    std::ofstream file("result.vox");

    if (file.is_open())
    {
        for (uint i = 0; i < pointInterieurs.size(); i++)
        {
            for (uint j = 0; j < pointInterieurs[i].size(); j++)
            {
                file << "" << pointInterieurs[i][j] << std::endl;
            }
        }

        file.close();

        std::cout << "Succesfully exported voxels to result.vox" << std::endl;
    }
    else
    {
        std::cout << "Failed to export voxels." << std::endl;
    }

    viewer << Viewer3D<>::updateDisplay;

    // Return the qT application.
    return application.exec();
}
