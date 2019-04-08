#include "Algorithms.h"

using namespace DGtal;

bool AbetweenBandC(Viewer3D<>::RealPoint A, Viewer3D<>::RealPoint B,
                   Viewer3D<>::RealPoint C)
{

    if (A[0] < std::min(B[0], C[0]) || A[0] > std::max(B[0], C[0]))
    {
        return false;
    }
    if (A[1] < std::min(B[1], C[1]) || A[0] > std::max(B[1], C[1]))
    {
        return false;
    }
    if (A[2] < std::min(B[2], C[2]) || A[2] > std::max(B[2], C[2]))
    {
        return false;
    }

    return true;
}

bool AinsideBoundingBox(Viewer3D<>::RealPoint A, Viewer3D<>::RealPoint min,
                        Viewer3D<>::RealPoint max)
{
    if (A[0] < min[0] || A[0] > max[0])
    {
        return false;
    }
    if (A[1] < min[1] || A[1] > max[1])
    {
        return false;
    }
    if (A[2] < min[2] || A[2] > max[2])
    {
        return false;
    }
    return true;
}

Viewer3D<>::RealPoint planDirection(Viewer3D<>::RealPoint a, Viewer3D<>::RealPoint b,
                                    Viewer3D<>::RealPoint c, Viewer3D<> &view, Viewer3D<>::RealPoint min,
                                    Viewer3D<>::RealPoint max)
{
    Viewer3D<>::RealPoint u = b - a;
    Viewer3D<>::RealPoint v = c - a;

    // A x, B y, C z du plan.
    double mA = u[1] * v[2] - u[2] * v[1];
    double mB = u[2] * v[0] - u[0] * v[2];
    double mC = u[0] * v[1] - u[1] * v[0];

    double mD = -(mA * a[0] + mB * a[1] + mC * a[2]);

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

    bBoxArcs.push_back(B - A);
    bBoxArcs.push_back(A);
    bBoxArcs.push_back(C - B);
    bBoxArcs.push_back(B);
    bBoxArcs.push_back(D - C);
    bBoxArcs.push_back(C);
    bBoxArcs.push_back(D - A);
    bBoxArcs.push_back(A);
    bBoxArcs.push_back(E - H);
    bBoxArcs.push_back(H);
    bBoxArcs.push_back(F - E);
    bBoxArcs.push_back(E);
    bBoxArcs.push_back(G - F);
    bBoxArcs.push_back(F);
    bBoxArcs.push_back(H - G);
    bBoxArcs.push_back(G);
    bBoxArcs.push_back(H - A);
    bBoxArcs.push_back(A);
    bBoxArcs.push_back(E - D);
    bBoxArcs.push_back(D);
    bBoxArcs.push_back(F - C);
    bBoxArcs.push_back(C);
    bBoxArcs.push_back(G - B);
    bBoxArcs.push_back(B);

    for (int i = 0; i < bBoxArcs.size(); i += 2)
    {
        double t = -((d.dot(bBoxArcs[i + 1]) + mD) / d.dot(bBoxArcs[i]));
        Viewer3D<>::RealPoint point = bBoxArcs[i + 1] + (bBoxArcs[i] * t);

        if (AbetweenBandC(point, bBoxArcs[i + 1], bBoxArcs[i + 1] + bBoxArcs[i]))
        {
            if (AinsideBoundingBox(point, min, max))
            {
                planCut.push_back(point);
            }
        }
    }

    trace.info() << "plancut size:" << planCut.size();
    for (int i = 0; i < planCut.size(); i++)
    {
        if (i == planCut.size() - 1)
        {
            view.addLine(planCut[i], planCut[0], 0.05);
        }
        else
        {
            view.addLine(planCut[i], planCut[i + 1], 0.05);
        }
    }

    return d;
}

void DisplayBoundingBox(Viewer3D<> &view, Viewer3D<>::RealPoint min,
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
