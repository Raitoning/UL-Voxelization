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
