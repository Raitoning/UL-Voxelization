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

// Möller-Trumborne algorithm
// bool RayIntersectsTriangle(Viewer3D<>::RealPoint rayOrigin,
//                            Viewer3D<>::RealPoint rayVector,
//                            Viewer3D<>::RealPoint a,
//                            Viewer3D<>::RealPoint b,
//                            Viewer3D<>::RealPoint c,
//                            Viewer3D<>::RealPoint &outIntersectionPoint)
// {

//     // Normal
//     Z3i::RealPoint ab = b - a;
//     Z3i::RealPoint bc = c - b;
//     Z3i::RealPoint ac = c - a;
//     Z3i::RealPoint ca = a - c;

//     Z3i::RealPoint normal = ab.crossProduct(ac);

//     float dot = normal.dot(rayVector);

//     // Rayon parallèle
//     if (dot < FLT_EPSILON && dot > -FLT_EPSILON)
//     {
//         return false;
//     }

//     // NOTE: Ajouter test triangle derrière l'origine du rayon
//     float distance = normal.dot(a);

//     float t = (normal.dot(rayOrigin) + distance / normal.dot(rayVector));

//     Z3i::RealPoint intersection = rayOrigin + rayVector * t;

//     Z3i::RealPoint planeNormal;

//     Z3i::RealPoint pa = intersection - a;
//     planeNormal = ab.crossProduct(pa);
//     if (normal.dot(planeNormal) < 0)
//     {
//         return false;
//     }

//     Z3i::RealPoint pb = intersection - b;
//     planeNormal = bc.crossProduct(pb);
//     if (normal.dot(planeNormal) < 0)
//     {
//         return false;
//     }

//     Z3i::RealPoint pc = intersection - c;
//     planeNormal = ca.crossProduct(pc);
//     if (normal.dot(planeNormal) < 0)
//     {
//         return false;
//     }

//     outIntersectionPoint = rayOrigin + rayVector * t;

//     return true;
// }

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

    // TODO: calculer distance
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

    // Change the scale of the mesh if it's too small.
    // mesh.changeScale(100);

    // Getting the bounding box.
    std::pair<Viewer3D<>::RealPoint, Viewer3D<>::RealPoint> boundingBox = mesh.getBoundingBox();
    DisplayBoundingBox(viewer, boundingBox.first, boundingBox.second);

    trace.info() << "Bounding box: " << std::endl
                 << boundingBox.first << std::endl
                 << boundingBox.second << std::endl;

    // Creating a single test ray from below,
    // At the center of the bounding box,
    // And facing upwards.
    Z3i::RealPoint rayOrigin((boundingBox.second[0] + boundingBox.first[0]) / 2, boundingBox.first[1] - 1, (boundingBox.first[2] + boundingBox.first[2]) / 2);
    Z3i::RealPoint rayDirection(0, 1, 0);
    Z3i::RealPoint intersection;

    float xStep = ((boundingBox.second[0] - boundingBox.first[0]) / float(horizontalResolution));
    float yStep = ((boundingBox.second[1] - boundingBox.first[1]) / float(verticalResolution));
    float zStep = ((boundingBox.second[2] - boundingBox.first[2]) / float(forwardResolution));

    if (!gaussian)
    {
        // Raytracing from under.
        rayDirection = Z3i::RealPoint(0, 1, 0);

        LOG("Min bbox:" << boundingBox.first[0] << " " << boundingBox.first[1] << " " << boundingBox.first[2]);
        LOG("Max bbox:" << boundingBox.second[0] << " " << boundingBox.second[1] << " " << boundingBox.second[2]);

        for (float x = boundingBox.first[0] - xStep; x < boundingBox.second[0]; x += xStep)
        {
            for (float z = boundingBox.first[2] - zStep; z < boundingBox.second[2]; z += zStep)
            {
                LOG("x;z: " << x << " " << z);
                Z3i::RealPoint rayOrigin(x, boundingBox.first[1] - 1, z);

                // intersectionPoints.push_back(rayOrigin);
                // intersectionPoints.push_back(rayOrigin + rayDirection * 2);
                // Test if the test ray can intersect anything.
                // for (int i = 0; i < mesh.nbFaces(); i++)
                // {
                //     // If a face is intersected, set it's color to red.
                //     if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                //     {
                //         mesh.setFaceColor(i, Color(255, 0, 0));
                //         trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
                //     }
                // }
            }
        }

        // Raytracing from in front of.
        rayDirection = Z3i::RealPoint(0, 0, -1);

        for (float x = boundingBox.first[0]; x < boundingBox.second[0]; x += xStep)
        {
            for (float y = boundingBox.first[1]; y < boundingBox.second[1]; y += yStep)
            {
                Z3i::RealPoint rayOrigin(x, y, boundingBox.first[2] + 1);

                // intersectionPoints.push_back(rayOrigin);
                // intersectionPoints.push_back(rayOrigin + rayDirection * 10);
                // Test if the test ray can intersect anything.
                // for (int i = 0; i < mesh.nbFaces(); i++)
                // {
                //     // If a face is intersected, set it's color to red.
                //     if (RayIntersectsTriangle(rayOrigin, rayDirection, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
                //     {
                //         intersectionPoints.push_back(rayOrigin);
                //         intersectionPoints.push_back(rayOrigin + rayDirection * 10);
                //         mesh.setFaceColor(i, Color(255, 0, 0));
                //         trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
                //     }
                // }
            }
        }

        // Raytracing from the left.
        rayDirection = Z3i::RealPoint(1, 0, 0);

        for (float y = boundingBox.first[1]; y < boundingBox.second[1]; y += yStep)
        {
            for (float z = boundingBox.first[2]; z < boundingBox.second[2]; z += zStep)
            {
                Z3i::RealPoint rayOrigin(boundingBox.first[0] - 1, y, z);

                intersectionPoints.push_back(rayOrigin);
                intersectionPoints.push_back(rayOrigin + rayDirection);
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
