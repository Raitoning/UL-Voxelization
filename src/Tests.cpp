// Axes:
// X: Horizontal
// Y: Vertical (Y up)
// Z: Forward/Backward (-Z forward)

#include <DGtal/base/Common.h>
#include <DGtal/io/Display3D.h>
#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/viewers/Viewer3D.h>
#include "Algorithms.h"

using namespace DGtal;

bool RunTest(Mesh<Viewer3D<>::RealPoint> &mesh, Z3i::RealPoint &rayOrigin, Z3i::RealPoint &rayDirection, std::vector<Z3i::RealPoint> &intersectionPoints)
{
    intersectionPoints.push_back(rayOrigin);
    intersectionPoints.push_back(rayOrigin + rayDirection * 3);

    bool hasIntersections = false;

    std::pair<Viewer3D<>::RealPoint, Viewer3D<>::RealPoint> boundingBox = mesh.getBoundingBox();

    Z3i::RealPoint intersection;

    // Test if the test ray can intersect anything.
    for (uint i = 0; i < mesh.nbFaces(); i++)
    {
        //If a face is intersected, set it's color to red.
        if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
        {
            hasIntersections = true;
            // mesh.setFaceColor(i, Color(255, 0, 0));
            trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
        }
    }

    return hasIntersections;
}

int main(int argc, char **argv)
{
    // File name of the mesh to import.
    std::string inputFile;

    std::vector<Z3i::RealPoint> intersectionPoints;

    inputFile = "../Mesh/Triangle.off";

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

    DisplayBoundingBox(viewer, boundingBox.first, boundingBox.second);

    trace.info() << "Bounding box: " << std::endl
                 << boundingBox.first << std::endl
                 << boundingBox.second << std::endl;

    trace.info() << std::endl
                 << "Now running unit tests..." << std::endl
                 << std::endl;
    // NOTE: Unit tests
    Z3i::RealPoint rayOrigin;
    Z3i::RealPoint rayDirection;
    Z3i::RealPoint intersection;

    trace.info() << "Tests with a single triangle..." << std::endl
                 << std::endl;

    // NOTE: ray orthogonal to a face
    rayOrigin = Z3i::RealPoint(0, 0, 1);
    rayDirection = Z3i::RealPoint(0, 0, -1);
    if (!RunTest(mesh, rayOrigin, rayDirection, intersectionPoints))
    {
        trace.error() << "No intersection when a ray is orthogonal to a face" << std::endl;
    }
    else
    {
        trace.warning() << "Ray orthogonal to a face: Success" << std::endl;
    }

    // NOTE: ray orthogonal to an edge
    rayOrigin = Z3i::RealPoint(0, -0.5, 0);
    rayDirection = Z3i::RealPoint(0, 0, -1);
    if (!RunTest(mesh, rayOrigin, rayDirection, intersectionPoints))
    {
        trace.error() << "No intersection when a ray is orthogonal to an edge" << std::endl;
    }
    else
    {
        trace.warning() << "Ray orthogonal to an edge: Success" << std::endl;
    }

    // NOTE: ray tangent to an edge
    rayOrigin = Z3i::RealPoint(-1, -0.5, 0);
    rayDirection = Z3i::RealPoint(1, 0, 0);
    if (RunTest(mesh, rayOrigin, rayDirection, intersectionPoints))
    {
        trace.error() << "Untersection when a ray is tangent to an edge" << std::endl;
    }
    else
    {
        trace.warning() << "Ray tangent to an edge: Success" << std::endl;
    }

    // NOTE: ray tangent to a vertex
    rayOrigin = Z3i::RealPoint(-1, 0.5, 0);
    rayDirection = Z3i::RealPoint(1, 0, 0);
    if (RunTest(mesh, rayOrigin, rayDirection, intersectionPoints))
    {
        trace.error() << "Intersection when a ray is tangent to a vertex" << std::endl;
    }
    else
    {
        trace.warning() << "Ray tangent to a vertex: Success" << std::endl;
    }

    trace.info() << std::endl
                 << "Tests with two triangles..." << std::endl
                 << std::endl;

    inputFile = "../Mesh/2Triangles.off";

    mesh = Mesh<Z3i::RealPoint>();
    mesh << inputFile;

    // NOTE: ray orthogonal to the edge of 2 triangles
    rayOrigin = Z3i::RealPoint(0, 0, 1);
    rayDirection = Z3i::RealPoint(0, 0, -1);
    if (!RunTest(mesh, rayOrigin, rayDirection, intersectionPoints))
    {
        trace.error() << "No intersection when a ray is orthogonal to the edge of 2 triangles" << std::endl;
    }
    else
    {
        trace.warning() << "Ray orthogonal to the edge of 2 triangles: Success" << std::endl;
    }

    // NOTE: ray tangent to the edge of 2 triangles
    rayOrigin = Z3i::RealPoint(0, -1, 0);
    rayDirection = Z3i::RealPoint(0, 1, 0);
    if (RunTest(mesh, rayOrigin, rayDirection, intersectionPoints))
    {
        trace.error() << "Intersection when a ray is tangent to the edge of 2 triangles" << std::endl;
    }
    else
    {
        trace.warning() << "Ray tangent to the edge of 2 triangles: Success" << std::endl;
    }

    // NOTE: ray tangent to a vertex of 2 triangles
    rayOrigin = Z3i::RealPoint(-1, 0.5, 0);
    rayDirection = Z3i::RealPoint(1, 0, 0);
    if (RunTest(mesh, rayOrigin, rayDirection, intersectionPoints))
    {
        trace.error() << "Intersection when a ray is tangent to a vertex of 2 triangles" << std::endl;
    }
    else
    {
        trace.warning() << "Ray tangent to a vertex of 2 triangles: Success" << std::endl;
    }

    // NOTE: ray orthogonal to a vertex of 2 triangles
    rayOrigin = Z3i::RealPoint(0, 0.5, 1);
    rayDirection = Z3i::RealPoint(0, 0, -1);
    if (!RunTest(mesh, rayOrigin, rayDirection, intersectionPoints))
    {
        trace.error() << "No intersection when a ray is orthogonal to a vertex of 2 triangles" << std::endl;
    }
    else
    {
        trace.warning() << "Ray orthogonal to a vertex of 2 triangles: Success" << std::endl;
    }

    // NOTE: ray not hitting anything
    rayOrigin = Z3i::RealPoint(0.5, 0.5, 0.5);
    rayDirection = Z3i::RealPoint(0, 0, -1);

    if (RunTest(mesh, rayOrigin, rayDirection, intersectionPoints))
    {
        trace.error() << "Intersection when a ray doesn't hit anything" << std::endl;
    }
    else
    {
        trace.warning() << "Ray not hitting anything: Success" << std::endl;
    }

    // Push the mesh into the viewer.
    viewer << mesh;

    intersectionPoints.push_back(boundingBox.first);
    intersectionPoints.push_back(boundingBox.second);

    for (uint i = 0; i < intersectionPoints.size(); i += 2)
    {
        viewer.addLine(intersectionPoints[i], intersectionPoints[i + 1], 0.03);
    }

    viewer << Viewer3D<>::updateDisplay;

    // Return the qT application.
    return application.exec();
}
