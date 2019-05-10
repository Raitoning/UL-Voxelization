#pragma once

#include <DGtal/io/viewers/Viewer3D.h>

using namespace DGtal;

// FIXME: Do not create homemade structures
struct indexation
{
    int index;
    double value;
    bool operator<(const indexation &a) const
    {
        return value < a.value;
    }
};

// struct stockage
// {
//     Viewer3D<>::RealPoint point;
//     int qte;
// };

void DisplayBoundingBox(Viewer3D<> &view, Viewer3D<>::RealPoint min,
                        Viewer3D<>::RealPoint max);

// Badouel's algorithm
// https://graphics.stanford.edu/courses/cs348b-98/gg/intersect.html
bool RayIntersectsTriangle(Viewer3D<>::RealPoint rayOrigin,
                           Viewer3D<>::RealPoint rayDirection,
                           Viewer3D<>::RealPoint a,
                           Viewer3D<>::RealPoint b,
                           Viewer3D<>::RealPoint c,
                           Viewer3D<>::RealPoint &outIntersectionPoint);

bool realPointEquals(Viewer3D<>::RealPoint pointA, Viewer3D<>::RealPoint pointB);

Viewer3D<>::RealPoint createStep(Viewer3D<>::RealPoint dir, double ratioX, double ratioY, double ratioZ);

vector<Viewer3D<>::RealPoint> pointInterieur(Viewer3D<>::RealPoint origin, Viewer3D<>::RealPoint dir, vector<Viewer3D<>::RealPoint> intersects, Viewer3D<>::RealPoint step, std::pair<Viewer3D<>::RealPoint, Viewer3D<>::RealPoint> boundingBox);

std::vector<Viewer3D<>::RealPoint> retirerDouble(std::vector<Viewer3D<>::RealPoint> valeurs);

void conservationSurface(int ***voxels, int a, int b, int c, int seuil);

// FIXME: Use std::maps instead of homemade functions.
// bool realPointEquals(Viewer3D<>::RealPoint pointA, Viewer3D<>::RealPoint pointB);

// void addResult(vector<stockage> result, vector<Viewer3D<>::RealPoint> listePoint);

// vector<Viewer3D<>::RealPoint> computeVote(vector<stockage> resultats, int seuil);

// Algorithms used to create non orthogonal planes.

bool AbetweenBandC(Viewer3D<>::RealPoint A, Viewer3D<>::RealPoint B,
                   Viewer3D<>::RealPoint C);

bool AinsideBoundingBox(Viewer3D<>::RealPoint A, Viewer3D<>::RealPoint min,
                        Viewer3D<>::RealPoint max);

Viewer3D<>::RealPoint planDirection(Viewer3D<>::RealPoint a, Viewer3D<>::RealPoint b,
                                    Viewer3D<>::RealPoint c, Viewer3D<> &view, Viewer3D<>::RealPoint min,
                                    Viewer3D<>::RealPoint max);

float distance_to_origin(Viewer3D<>::RealPoint p);

Viewer3D<>::RealPoint normalisePoint(Viewer3D<>::RealPoint p);

vector<Viewer3D<>::RealPoint> orthogonalDirections(Viewer3D<>::RealPoint direction);

vector<Viewer3D<>::RealPoint> intersectsBoundingBoxCore(Viewer3D<>::RealPoint rayOrigin, Viewer3D<>::RealPoint rayDirection,
                                                        Viewer3D<>::RealPoint min, Viewer3D<>::RealPoint max);

bool intersectsBoundingBox(Viewer3D<>::RealPoint rayOrigin, Viewer3D<>::RealPoint rayDirection,
                           Viewer3D<>::RealPoint min, Viewer3D<>::RealPoint max);

vector<Viewer3D<>::RealPoint> intersectsBoundingBoxReturnsPoint(Viewer3D<>::RealPoint rayOrigin, Viewer3D<>::RealPoint rayDirection,
                                                                Viewer3D<>::RealPoint min, Viewer3D<>::RealPoint max);

void originPointsRecursive(Viewer3D<>::RealPoint o, Viewer3D<>::RealPoint dir, Viewer3D<>::RealPoint a, Viewer3D<>::RealPoint b,
                           vector<Viewer3D<>::RealPoint> &folder, Viewer3D<>::RealPoint min, Viewer3D<>::RealPoint max);

vector<Viewer3D<>::RealPoint> originPoints(Viewer3D<>::RealPoint origin, Viewer3D<>::RealPoint normale,
                                           Viewer3D<>::RealPoint min, Viewer3D<>::RealPoint max, double delta = 1.0);
