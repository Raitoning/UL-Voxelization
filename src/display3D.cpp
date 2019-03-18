// Axes:
// X: Horizontal
// Y: Vertical (Y up)
// Z: Forward/Backward (-Z forward)

#include <DGtal/base/Common.h>
#include <DGtal/io/Display3D.h>
#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/viewers/Viewer3D.h>

using namespace DGtal;

// TODO: Gaussian voxelization
// TODO: Raytracing / Raymarching
// TODO: Equations plane & lines
// TODO: Intersections
// TODO: Vector vertices / faces & display
// TODO: Test -> Direction
// TODO: Resolution
// TODO: Unit tests

// Direction Vector from 2 points (A->B)

//Viewer3D<>::RealPoint directionVector(Viewer3D<>::RealPoint a,
//                                      Viewer3D<>::RealPoint b) {
//  Viewer3D<>::RealPoint s = b-a;
//  return s;
//}

int DisplayBoundingBox(Viewer3D<> &view, Viewer3D<>::RealPoint min,
                                         Viewer3D<>::RealPoint max) {
    Viewer3D<>::RealPoint A = max;
    A[0] = min[0]; 
    Viewer3D<>::RealPoint B = max;
    
    Viewer3D<>::RealPoint C =max;
    C[2] = min[2];
    Viewer3D<>::RealPoint D =min;
    D[1] = max[1];
    Viewer3D<>::RealPoint E =min;
    
    Viewer3D<>::RealPoint F =min;
    F[0] = max[0];
    Viewer3D<>::RealPoint G =max;
    G[1] = min[1];
    Viewer3D<>::RealPoint H =min;
    H[2] = max[2];
    
    view.addLine(A,B,0.05);
    view.addLine(B,C,0.05);
    view.addLine(C,D,0.05);
    view.addLine(D,A,0.05);
    view.addLine(H,E,0.05);
    view.addLine(E,F,0.05);
    view.addLine(F,G,0.05);
    view.addLine(G,H,0.05);
    view.addLine(A,H,0.05);
    view.addLine(D,E,0.05);
    view.addLine(C,F,0.05);
    view.addLine(B,G,0.05);
    
    return 0;
    
}

// Möller-Trumborne algorithm
bool RayIntersectsTriangle(Viewer3D<>::RealPoint rayOrigin,
                           Viewer3D<>::RealPoint rayVector,
                           Viewer3D<>::RealPoint a,
                           Viewer3D<>::RealPoint b,
                           Viewer3D<>::RealPoint c,
                           Viewer3D<>::RealPoint &outIntersectionPoint)
{

    // Normal
    Z3i::RealPoint ab = b - a;
    Z3i::RealPoint bc = c - b;
    Z3i::RealPoint ac = c - a;
    Z3i::RealPoint ca = a - c;

    Z3i::RealPoint normal = ab.crossProduct(ac);

    float dot = normal.dot(rayVector);

    // Rayon parallèle
    if (dot < FLT_EPSILON && dot > -FLT_EPSILON)
    {
        return false;
    }

    // NOTE: Ajouter test triangle derrière l'origine du rayon
    float distance = normal.dot(a);

    float t = (normal.dot(rayOrigin) + distance / normal.dot(rayVector));

    Z3i::RealPoint intersection = rayOrigin + rayVector * t;

    Z3i::RealPoint planeNormal;

    Z3i::RealPoint pa = intersection - a;
    planeNormal = ab.crossProduct(pa);
    if (normal.dot(planeNormal) < 0)
    {
        return false;
    }

    Z3i::RealPoint pb = intersection - b;
    planeNormal = bc.crossProduct(pb);
    if (normal.dot(planeNormal) < 0)
    {
        return false;
    }

    Z3i::RealPoint pc = intersection - c;
    planeNormal = ca.crossProduct(pc);
    if (normal.dot(planeNormal) < 0)
    {
        return false;
    }

    outIntersectionPoint = rayOrigin + rayVector * t;

    return true;
}

int main(int argc, char **argv)
{
    // File name of the mesh to import.
    std::string inputFile;

    // Checking the arguments.
    // If there isn't a file name as argument, stop the execution.
    if (argc != 2)
    {
        trace.info() << "Usage: display3D file" << std::endl
                     << "Now exiting..." << std::endl;

        return 0;
    }
    else
    {
        // Get the file name to import the mesh.
        inputFile = argv[1];
        trace.info() << "Input file: " << inputFile << std::endl;
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

    // Change the scale of the mesh if it's too small.
    mesh.changeScale(100);

    // Getting the bounding box.
    std::pair<Viewer3D<>::RealPoint, Viewer3D<>::RealPoint> boundingBox = mesh.getBoundingBox();
    Z3i::Point minBoundingBox = boundingBox.first;
    Z3i::Point maxBoundingBox = boundingBox.second;
  
  DisplayBoundingBox(viewer, boundingBox.first, boundingBox.second);

    // Setting the bounding box to the nearest integer values.
    minBoundingBox[0] = int(minBoundingBox[0]) - 1;
    minBoundingBox[1] = int(minBoundingBox[1]) - 1;

    maxBoundingBox[0] = int(maxBoundingBox[0]) + 1;
    maxBoundingBox[1] = int(maxBoundingBox[1]) + 1;

    // Drawing a domain based on the bouding box.
    // Z3i::Domain domain(minBoundingBox, maxBoundingBox);
    // viewer << domain;

    trace.info() << "Bounding box: " << std::endl
                 << boundingBox.first << std::endl
                 << boundingBox.second << std::endl;

    // Creating a single test ray from below,
    // At the center of the bounding box,
    // And facing upwards.
    Z3i::RealPoint rayOrigin((minBoundingBox[0] + maxBoundingBox[0]) / 2, minBoundingBox[1] - 1, (minBoundingBox[2] + maxBoundingBox[2]) / 2);
    Z3i::RealPoint rayDirection(0, 1, 0);
    Z3i::RealPoint intersection;

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

    // Push the mesh into the viewer.
    viewer << mesh;
    viewer << Viewer3D<>::updateDisplay;

    // Return the qT application.
    return application.exec();
}
