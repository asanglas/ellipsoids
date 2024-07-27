#pragma once
#include "stdio.h"
#include "math.h"
#include "ellipsoids.h"

// ----------------------------
//
// References
//
// ----------------------------

// Mainly from this one (first algorithm): https://www.researchgate.net/publication/4352232_A_Note_on_Approximate_Minimum_Volume_Enclosing_Ellipsoid_of_Ellipsoids

// Borrow things from here https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=7ba3adaeeb8f411a75904a9395ec61d73404f94e

// ----------------------------
//
// Linear Algebra functions
//
// ----------------------------

Vec2 mat2x2vec2_mul(Mat2 mat, Vec2 vec)
{
    Vec2 result;
    result.x = mat.a11 * vec.x + mat.a12 * vec.y;
    result.y = mat.a21 * vec.x + mat.a22 * vec.y;
    return result;
}

float vec2vec2_dot(Vec2 vec1, Vec2 vec2)
{
    return vec1.x * vec2.x + vec1.y * vec2.y;
}

Mat2 vec2vec2_outer(Vec2 vec1, Vec2 vec2)
{
    Mat2 outer;
    outer.a11 = vec1.x * vec2.x;
    outer.a12 = vec1.x * vec2.y;
    outer.a21 = vec1.y * vec2.x;
    outer.a22 = vec1.y * vec2.y;
    return outer;
}

Mat2 mat2_inverse(Mat2 m)
{
    // Use an analytical formula. https://www.dr-lex.be/random/matrix-inv.html
    float det = m.a11 * m.a22 - m.a12 * m.a21;
    Mat2 inv;
    inv.a11 = m.a22 / det;
    inv.a12 = -m.a12 / det;
    inv.a21 = -m.a21 / det;
    inv.a22 = m.a11 / det;
    return inv;
}

Mat3 mat3_inverse(Mat3 m)
{
    // Use an analytical formula. https://www.dr-lex.be/random/matrix-inv.html
    // We are dealing only with symmetric matrices so it could be optimized.
    // For 4x4 maybe this is useful https://math.stackexchange.com/questions/1033611/inverse-4x4-matrix

    float det = m.a11 * (m.a33 * m.a22 - m.a32 * m.a23) - m.a21 * (m.a33 * m.a12 - m.a32 * m.a13) + m.a32 * (m.a23 * m.a12 - m.a22 * m.a12);

    Mat3 inv;
    inv.a11 = (m.a33 * m.a22 - m.a32 * m.a23) / det;
    inv.a12 = -(m.a33 * m.a12)
}

void volume_approximation(Vec2* initial_points, Vec2* initial_set, u32 n)

{
}

Ellipsoid2d compute_bounding_ellipsoid_iterative_2d(Ellipsoid2d* ellipsoids, u32 count)
{
    Ellipsoid2d bounding_ellipsoid;

    return bounding_ellipsoid;
}
