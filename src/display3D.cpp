// Axes:
// X: Horizontal
// Y: Vertical (Y up)
// Z: Forward/Backward (-Z forward)

#include <DGtal/base/Common.h>
#include <DGtal/io/Display3D.h>
#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/viewers/Viewer3D.h>

#define LOG(X) std::cout << X << std::endl

using namespace DGtal;

// TODO: Gaussian voxelization
// TODO: Equations plane & lines

Viewer3D<>::RealPoint planDirection(Viewer3D<>::RealPoint a, Viewer3D<>::RealPoint b,
                                    Viewer3D<>::RealPoint c)
{
    Viewer3D<>::RealPoint u = b - a;
    Viewer3D<>::RealPoint v = c - a;

    // A x, B y, C z du plan.
    float mA = u[1] * v[2] - u[2] * v[1];
    float mB = u[2] * v[0] - u[0] * v[2];
    float mC = u[0] * v[1] - u[1] * v[0];

    float mD = -(mA * a[0] + mB * a[1] + mC * a[2]);

    Viewer3D<>::RealPoint d = a;

    d[0] = mA;
    d[1] = mB;
    d[2] = mC;

    return d;
}

int DisplayBoundingBox(Viewer3D<> &view, Viewer3D<>::RealPoint min,
                       Viewer3D<>::RealPoint max)
{
    Viewer3D<>::RealPoint A = max;
    A[0] = min[0];
    Viewer3D<>::RealPoint B = max;

    Viewer3D<>::RealPoint C = max;
    C[2] = min[2];
    Viewer3D<>::RealPoint D = min;
    D[1] = max[1];
    Viewer3D<>::RealPoint E = min;

    Viewer3D<>::RealPoint F = min;
    F[0] = max[0];
    Viewer3D<>::RealPoint G = max;
    G[1] = min[1];
    Viewer3D<>::RealPoint H = min;
    H[2] = max[2];

    view.addLine(A, B, 0.05);
    view.addLine(B, C, 0.05);
    view.addLine(C, D, 0.05);
    view.addLine(D, A, 0.05);
    view.addLine(H, E, 0.05);
    view.addLine(E, F, 0.05);
    view.addLine(F, G, 0.05);
    view.addLine(G, H, 0.05);
    view.addLine(A, H, 0.05);
    view.addLine(D, E, 0.05);
    view.addLine(C, F, 0.05);
    view.addLine(B, G, 0.05);

    return 0;
}

// Badouel's algorithm
// https://graphics.stanford.edu/courses/cs348b-98/gg/intersect.html
bool RayIntersectsTriangle(Viewer3D<>::RealPoint rayOrigin,
                           Viewer3D<>::RealPoint rayDirection,
                           Viewer3D<>::RealPoint a,
                           Viewer3D<>::RealPoint b,
                           Viewer3D<>::RealPoint c,
                           Viewer3D<>::RealPoint &outIntersectionPoint)
{
    Viewer3D<>::RealPoint ab = b - a;
    Viewer3D<>::RealPoint ac = c - a;

    Viewer3D<>::RealPoint cb = b - c;
    Viewer3D<>::RealPoint ca = a - c;

    Viewer3D<>::RealPoint normal = ab.crossProduct(ac);

    float d = -a.dot(normal);

    float angle = normal.dot(rayDirection);

    // Test parallèle
    // Si angle entre -Epsilon et Epsilon
    if (angle < FLT_EPSILON && angle > -FLT_EPSILON)
    {
        return false;
    }

    float t = -((normal.dot(rayOrigin) + d) / normal.dot(rayDirection));

    // Calcul point intersection
    Viewer3D<>::RealPoint point = rayOrigin + (rayDirection * t);

    Viewer3D<>::RealPoint bc = c - b;
    Viewer3D<>::RealPoint cp = point - c;

    float det = cb[1] * ca[0] + bc[0] * ca[1];
    float factor_alpha = cb[1] * cp[0] + bc[0] * cp[1];
    float factor_beta = ac[1] * cp[0] + ca[0] * cp[1];
    float alpha = factor_alpha / det;
    float beta = factor_beta / det;
    float gamma = 1.0f - alpha - beta;

    if (alpha >= 0 && alpha <= 1)
    {
        if (beta >= 0 && beta <= 1)
        {
            if (gamma >= 0 && gamma <= 1)
            {
                outIntersectionPoint = point;
                return true;
            }
        }
    }

    return false;
}

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

    double scaleFactor = 1;

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

    float xStep = ((boundingBox.second[0] - boundingBox.first[0]) / float(horizontalResolution - 1));
    float yStep = ((boundingBox.second[1] - boundingBox.first[1]) / float(verticalResolution - 1));
    float zStep = ((boundingBox.second[2] - boundingBox.first[2]) / float(forwardResolution - 1));

    if (!gaussian)
    {
        // Raytracing from under.
        Z3i::RealPoint rayDirection = Z3i::RealPoint(0, 1, 0);

        LOG("Min bbox:" << boundingBox.first[0] << " " << boundingBox.first[1] << " " << boundingBox.first[2]);
        LOG("Max bbox:" << boundingBox.second[0] << " " << boundingBox.second[1] << " " << boundingBox.second[2]);

        for (float x = boundingBox.first[0]; x <= boundingBox.second[0]; x += xStep)
        {
            for (float z = boundingBox.first[2]; z <= boundingBox.second[2]; z += zStep)
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

        for (float x = boundingBox.first[0]; x <= boundingBox.second[0]; x += xStep)
        {
            for (float y = boundingBox.first[1]; y <= boundingBox.second[1]; y += yStep)
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

        for (float y = boundingBox.first[1]; y <= boundingBox.second[1]; y += yStep)
        {
            for (float z = boundingBox.first[2]; z <= boundingBox.second[2]; z += zStep)
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

        // NOTE: Unit tests
        Z3i::RealPoint rayOrigin;
        Z3i::RealPoint rayDirection;

        // NOTE: ray orthogonal to a face
        // rayOrigin = Z3i::RealPoint(0, 0, 1);
        // rayDirection = Z3i::RealPoint(0, 0, -1);

        // NOTE: ray orthogonal to an edge
        // rayOrigin = Z3i::RealPoint(0, -0.5, 0);
        // rayDirection = Z3i::RealPoint(0, 0, -1);

        // NOTE: ray tangent to an edge
        // rayOrigin = Z3i::RealPoint(-1, -0.5, 0);
        // rayDirection = Z3i::RealPoint(1, 0, 0);

        // NOTE: ray tangent to a vertex
        // rayOrigin = Z3i::RealPoint(-1, 0.5, 0);
        // rayDirection = Z3i::RealPoint(1, 0, 0);

        // NOTE: ray orthogonal to the edge of 2 triangles
        // rayOrigin = Z3i::RealPoint(0, 0, 1);
        // rayDirection = Z3i::RealPoint(0, 0, -1);

        // NOTE: ray tangent to the edge of 2 triangles
        // rayOrigin = Z3i::RealPoint(0, -1, 0);
        // rayDirection = Z3i::RealPoint(0, 1, 0);

        // NOTE: ray tangent to the edge of 2 triangles
        // rayOrigin = Z3i::RealPoint(-1, 0.5, 0);
        // rayDirection = Z3i::RealPoint(1, 0, 0);

        // NOTE: ray orthogonal to the edge of 2 triangles
        // rayOrigin = Z3i::RealPoint(0, 0.5, 1);
        // rayDirection = Z3i::RealPoint(0, 0, -1);

        // NOTE: ray not hitting anything
        rayOrigin = Z3i::RealPoint(0.5, 0.5, 0.5);
        rayDirection = Z3i::RealPoint(0, 0, -1);

        intersectionPoints.push_back(rayOrigin);
        intersectionPoints.push_back(rayOrigin + rayDirection * 3);

        // Test if the test ray can intersect anything.
        for (int i = 0; i < mesh.nbFaces(); i++)
        {
            //If a face is intersected, set it's color to red.
            if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
            {
                mesh.setFaceColor(i, Color(255, 0, 0));
                trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
            }
        }
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
