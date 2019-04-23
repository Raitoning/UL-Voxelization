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
// TODO: Test -> Direction
bool AbetweenBandC( Viewer3D<>::RealPoint A,Viewer3D<>::RealPoint B, 
                                            Viewer3D<>::RealPoint C){

    if (A[0]<std::min(B[0],C[0])||A[0]>std::max(B[0],C[0])){
        return false;
    }
    if (A[1]<std::min(B[1],C[1])||A[0]>std::max(B[1],C[1])){
        return false;
    }
    if (A[2]<std::min(B[2],C[2])||A[2]>std::max(B[2],C[2])){
        return false;
    }

    return true;
}

float distance_to_origin(Viewer3D<>::RealPoint p){

    return sqrtf(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

Viewer3D<>::RealPoint normalisePoint(Viewer3D<>::RealPoint p){

    if(distance_to_origin(p)!=0){
        return p * (1/distance_to_origin(p));
    }
    return p;
}

bool AinsideBoundingBox( Viewer3D<>::RealPoint A,Viewer3D<>::RealPoint min, 
                                            Viewer3D<>::RealPoint max){
    if (A[0]<min[0]||A[0]>max[0]){
        return false;
    }
    if (A[1]<min[1]||A[1]>max[1]){
        return false;
    }
    if (A[2]<min[2]||A[2]>max[2]){
        return false;
    }
    return true;
} 

Viewer3D<>::RealPoint planDirection(Viewer3D<>::RealPoint a, Viewer3D<>::RealPoint b,
                       Viewer3D<>::RealPoint c, Viewer3D<> &view, Viewer3D<>::RealPoint min,
                       Viewer3D<>::RealPoint max){
    Viewer3D<>::RealPoint u = b - a;
    Viewer3D<>::RealPoint v = c - a;
  
  
    // A x, B y, C z du plan.
    double mA = u[1]*v[2] - u[2]*v[1];
    double mB = u[2]*v[0] - u[0]*v[2];
    double mC = u[0]*v[1] - u[1]*v[0];
  
    double mD = -(mA*a[0] + mB*a[1] + mC*a[2]);
  
    Viewer3D<>::RealPoint d = u.crossProduct(v);

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

    vector<Viewer3D<>::RealPoint> bBoxArcs;

    vector<Viewer3D<>::RealPoint> planCut;

    bBoxArcs.push_back(B-A);
    bBoxArcs.push_back(A);
    bBoxArcs.push_back(C-B);
    bBoxArcs.push_back(B);
    bBoxArcs.push_back(D-C);
    bBoxArcs.push_back(C);
    bBoxArcs.push_back(D-A);
    bBoxArcs.push_back(A);
    bBoxArcs.push_back(E-H);
    bBoxArcs.push_back(H);
    bBoxArcs.push_back(F-E);
    bBoxArcs.push_back(E);
    bBoxArcs.push_back(G-F);
    bBoxArcs.push_back(F);
    bBoxArcs.push_back(H-G);
    bBoxArcs.push_back(G);
    bBoxArcs.push_back(H-A);
    bBoxArcs.push_back(A);
    bBoxArcs.push_back(E-D);
    bBoxArcs.push_back(D);
    bBoxArcs.push_back(F-C);
    bBoxArcs.push_back(C);
    bBoxArcs.push_back(G-B);
    bBoxArcs.push_back(B);

    for (int i = 0; i < bBoxArcs.size(); i += 2)
    {
        double t = -((d.dot(bBoxArcs[i+1]) + mD) / d.dot(bBoxArcs[i]));
        Viewer3D<>::RealPoint point = bBoxArcs[i+1] + (bBoxArcs[i] * t);

        if (AbetweenBandC(point,bBoxArcs[i+1],bBoxArcs[i+1]+bBoxArcs[i])){
            if (AinsideBoundingBox(point, min, max) ){
            planCut.push_back(point);
            }
        }

    }

    LOG("plancut size:" << planCut.size());
    for (int i = 0; i < planCut.size(); i++){
        if (i == planCut.size()-1) {
            view.addLine(planCut[i], planCut[0], 0.05);
        }
        else {
            view.addLine(planCut[i], planCut[i+1], 0.05);
        }
    }

    return d;
}

vector<Viewer3D<>::RealPoint> orthogonalDirections(Viewer3D<>::RealPoint direction){
    vector<Viewer3D<>::RealPoint> result;

    if (abs(direction[0])>abs(direction[1])&&abs(direction[0])>abs(direction[2])){
        Viewer3D<>::RealPoint a(direction[0], direction[2], direction[1]);
        Viewer3D<>::RealPoint b = a.crossProduct(direction);
        result.push_back(normalisePoint(a));
        result.push_back(normalisePoint(b));
    }
    else if (abs(direction[1])>abs(direction[0])&&abs(direction[1])>abs(direction[2])){
        Viewer3D<>::RealPoint a(direction[2], direction[1], direction[0]);
        Viewer3D<>::RealPoint b = a.crossProduct(direction);
        result.push_back(normalisePoint(a));
        result.push_back(normalisePoint(b));
    }
    else if (abs(direction[2])>abs(direction[0])&&abs(direction[2])>abs(direction[1])){

        Viewer3D<>::RealPoint a(direction[2], direction[1], direction[0]);
        Viewer3D<>::RealPoint b = a.crossProduct(direction);
        result.push_back(normalisePoint(a));
        result.push_back(normalisePoint(b));
        
    } else {//Most unlikely case, but still a posiblity.
        Viewer3D<>::RealPoint a(direction[0], direction[2], direction[1]);
        Viewer3D<>::RealPoint b = a.crossProduct(direction);
        a = direction.crossProduct(b);
        result.push_back(normalisePoint(a));
        result.push_back(normalisePoint(b));
    }

    return result;
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

bool intersectsBoundingBox(Viewer3D<>::RealPoint rayOrigin, Viewer3D<>::RealPoint rayDirection,
                            Viewer3D<>::RealPoint min, Viewer3D<>::RealPoint max)
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

    vector<Viewer3D<>::RealPoint> bBoxVertexes;

    bBoxVertexes.push_back(A);bBoxVertexes.push_back(B);bBoxVertexes.push_back(D);
    bBoxVertexes.push_back(C);bBoxVertexes.push_back(B);bBoxVertexes.push_back(D);
    
    bBoxVertexes.push_back(A);bBoxVertexes.push_back(B);bBoxVertexes.push_back(H);
    bBoxVertexes.push_back(G);bBoxVertexes.push_back(B);bBoxVertexes.push_back(H);
    
    bBoxVertexes.push_back(A);bBoxVertexes.push_back(E);bBoxVertexes.push_back(H);
    bBoxVertexes.push_back(A);bBoxVertexes.push_back(E);bBoxVertexes.push_back(D);

    bBoxVertexes.push_back(C);bBoxVertexes.push_back(F);bBoxVertexes.push_back(D);
    bBoxVertexes.push_back(E);bBoxVertexes.push_back(F);bBoxVertexes.push_back(D);

    bBoxVertexes.push_back(G);bBoxVertexes.push_back(F);bBoxVertexes.push_back(H);
    bBoxVertexes.push_back(E);bBoxVertexes.push_back(F);bBoxVertexes.push_back(H);

    bBoxVertexes.push_back(C);bBoxVertexes.push_back(G);bBoxVertexes.push_back(B);
    bBoxVertexes.push_back(C);bBoxVertexes.push_back(G);bBoxVertexes.push_back(F);

    int i =0;
    bool match = false;
    
     while ((!match) && (i < bBoxVertexes.size())){
        if (RayIntersectsTriangle(rayOrigin, rayDirection, bBoxVertexes[i], bBoxVertexes[i+1], bBoxVertexes[i+2],H))
        {
            match = true;
        }
        i+=3;
    }

    return match;
}

void originPointsRecursive(Viewer3D<>::RealPoint o, Viewer3D<>::RealPoint dir, Viewer3D<>::RealPoint a, Viewer3D<>::RealPoint b,
                            vector<Viewer3D<>::RealPoint> &folder, Viewer3D<>::RealPoint min, Viewer3D<>::RealPoint max){

    if(std::find(folder.begin(), folder.end(), o) == folder.end()) {
        if(intersectsBoundingBox(o,dir,min,max)){
            folder.push_back(o);
        LOG("Folder"<<folder.size());
            originPointsRecursive(o+a,dir,a,b,folder,min,max);
            originPointsRecursive(o-a,dir,a,b,folder,min,max);
            originPointsRecursive(o+b,dir,a,b,folder,min,max);
            originPointsRecursive(o-b,dir,a,b,folder,min,max);
        }   
    }

}

vector<Viewer3D<>::RealPoint> originPoints(Viewer3D<>::RealPoint origin, Viewer3D<>::RealPoint normale,
                                                Viewer3D<>::RealPoint min, Viewer3D<>::RealPoint max){
    vector<Viewer3D<>::RealPoint> result;

    vector<Viewer3D<>::RealPoint> d1d2=orthogonalDirections(normale);

    originPointsRecursive(origin,normalisePoint(normale),d1d2[0],d1d2[1],result,min,max);

    return result;    
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

    //DisplayBoundingBox(viewer, boundingBox.first, boundingBox.second);

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

        //planDirection(mesh.getVertex(mesh.getFace(0)[0]), mesh.getVertex(mesh.getFace(0)[1]), mesh.getVertex(mesh.getFace(0)[2]),viewer,boundingBox.first, boundingBox.second);

      
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

    

    //intersectionPoints.push_back(boundingBox.first);
    //intersectionPoints.push_back(boundingBox.second);
    Viewer3D<>::RealPoint d = boundingBox.second-boundingBox.first;
    vector<Viewer3D<>::RealPoint> testingVectors = originPoints(boundingBox.first, d, boundingBox.first, boundingBox.second);

    LOG("Size"<<testingVectors.size());

    for (int k=0;k < testingVectors.size();k++){
        for (int i = 0; i < mesh.nbFaces(); i++)
        {
            //If a face is intersected, set it's color to red.
            if (RayIntersectsTriangle(testingVectors[k], d, mesh.getVertex(mesh.getFace(i)[0]), mesh.getVertex(mesh.getFace(i)[1]), mesh.getVertex(mesh.getFace(i)[2]), intersection))
            {
                mesh.setFaceColor(i, Color(255, 0, 0));
                trace.info() << "Intersection at: (" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")" << std::endl;
            }
        }
    }

    // Push the mesh into the viewer.
    viewer << mesh;


    for (int i = 0; i < testingVectors.size(); i += 2)
    {
        viewer.addLine(testingVectors[i], testingVectors[i]+(d*2), 0.03+i);
    }

    viewer << Viewer3D<>::updateDisplay;

    // Return the qT application.
    return application.exec();
}
