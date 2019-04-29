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

bool AbetweenBandC(Viewer3D<>::RealPoint A, Viewer3D<>::RealPoint B,
                   Viewer3D<>::RealPoint C);

bool AinsideBoundingBox(Viewer3D<>::RealPoint A, Viewer3D<>::RealPoint min,
                        Viewer3D<>::RealPoint max);

Viewer3D<>::RealPoint planDirection(Viewer3D<>::RealPoint a, Viewer3D<>::RealPoint b,
                                    Viewer3D<>::RealPoint c, Viewer3D<> &view, Viewer3D<>::RealPoint min,
                                    Viewer3D<>::RealPoint max);

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

// FIXME: Use std::maps instead of homemade functions.

// bool realPointEquals(Viewer3D<>::RealPoint pointA, Viewer3D<>::RealPoint pointB);

Viewer3D<>::RealPoint createStep(Viewer3D<>::RealPoint dir, double ratioX, double ratioY, double ratioZ);

vector<Viewer3D<>::RealPoint> pointInterieur(Viewer3D<>::RealPoint origin, Viewer3D<>::RealPoint dir, vector<Viewer3D<>::RealPoint> intersects, Viewer3D<>::RealPoint step);

// void addResult(vector<stockage> result, vector<Viewer3D<>::RealPoint> listePoint);

// vector<Viewer3D<>::RealPoint> computeVote(vector<stockage> resultats, int seuil);
