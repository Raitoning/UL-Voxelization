#pragma once

#include <DGtal/io/viewers/Viewer3D.h>

bool AbetweenBandC(DGtal::Viewer3D<>::RealPoint A, DGtal::Viewer3D<>::RealPoint B,
                   DGtal::Viewer3D<>::RealPoint C);

bool AinsideBoundingBox(DGtal::Viewer3D<>::RealPoint A, DGtal::Viewer3D<>::RealPoint min,
                        DGtal::Viewer3D<>::RealPoint max);

DGtal::Viewer3D<>::RealPoint planDirection(DGtal::Viewer3D<>::RealPoint a, DGtal::Viewer3D<>::RealPoint b,
                                           DGtal::Viewer3D<>::RealPoint c, DGtal::Viewer3D<> &view, DGtal::Viewer3D<>::RealPoint min,
                                           DGtal::Viewer3D<>::RealPoint max);

void DisplayBoundingBox(DGtal::Viewer3D<> &view, DGtal::Viewer3D<>::RealPoint min,
                        DGtal::Viewer3D<>::RealPoint max);

// Badouel's algorithm
// https://graphics.stanford.edu/courses/cs348b-98/gg/intersect.html
bool RayIntersectsTriangle(DGtal::Viewer3D<>::RealPoint rayOrigin,
                           DGtal::Viewer3D<>::RealPoint rayDirection,
                           DGtal::Viewer3D<>::RealPoint a,
                           DGtal::Viewer3D<>::RealPoint b,
                           DGtal::Viewer3D<>::RealPoint c,
                           DGtal::Viewer3D<>::RealPoint &outIntersectionPoint);

bool realPointEquals(Viewer3D<>::RealPoint pointA, Viewer3D<>::RealPoint pointB);

Viewer3D<>::RealPoint createStep(Viewer3D<>::RealPoint dir, double ratioX, double ratioY, double ratioZ);

vector<Viewer3D<>::RealPoint> pointInterieur(Viewer3D<>::RealPoint origin, Viewer3D<>::RealPoint dir, vector<Viewer3D<>::RealPoint> intersects, Viewer3D<>::RealPoint step);

void addResult(vector<stockage> result, vector<Viewer3D<>::RealPoint> listePoint);

struct indexation{
    int index;
    double value;
    bool operator<(const indexation& a) const{
        return value < a.value;
    }
};

struct stockage{
    Viewer3D<>::RealPoint point;
    int qte;
};