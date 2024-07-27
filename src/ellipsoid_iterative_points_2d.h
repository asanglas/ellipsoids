#pragma once
#include "stdio.h"
#include "math.h"
#include "ellipsoids.h"

// ----------------------------
//
// References
//
// ----------------------------

// An adaptation of https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=7ba3adaeeb8f411a75904a9395ec61d73404f94e for points only.
// Possible optimization: https://people.orie.cornell.edu/miketodd/TYKhach.pdf (modification of the algorithm)
// Another optimization: https://www.es.mdh.se/pdf_publications/5680.pdf
// Another (didn't check): https://journal.austms.org.au/ojs/index.php/ANZIAMJ/article/view/17956/2361
// TODO next: The algorithm for ellipses (as intented in https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=7ba3adaeeb8f411a75904a9395ec61d73404f94e)

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

Vec2 vec2_normalize(Vec2 vec)
{
    Vec2 n;
    float norm = sqrt(vec.x * vec.x + vec.y * vec.y);
    n.x = vec.x / norm;
    n.y = vec.y / norm;

    return n;
}

Mat2 mat2x2_sum(Mat2 mat1, Mat2 mat2)
{
    Mat2 sum;
    sum.a11 = mat1.a11 + mat2.a11;
    sum.a12 = mat1.a12 + mat2.a12;
    sum.a21 = mat1.a21 + mat2.a21;
    sum.a22 = mat1.a22 + mat2.a22;
    return sum;
}

Mat2 mat2x2_muls(Mat2 mat, float s)
{
    Mat2 m;
    m.a11 = mat.a11 * s;
    m.a12 = mat.a12 * s;
    m.a21 = mat.a21 * s;
    m.a22 = mat.a22 * s;
    return m;
}

Vec2 mat2x2_eigenvalues(Mat2 mat)
{
    // assume symmetric
    float a = mat.a11, b = mat.a22, c = mat.a12;
    float delta = sqrtf(4 * pow(c, 2) + pow(a - b, 2));
    float lambda1 = (a + b - delta) / 2;
    float lambda2 = (a + b + delta) / 2;
    return (Vec2){.x = lambda1, .y = lambda2};
}

Mat2 mat2x2_eigenvectors(Mat2 mat)
{
    // assume symmetric
    float a = mat.a11, b = mat.a22, c = mat.a12;
    float delta = sqrtf(4 * pow(c, 2) + pow(a - b, 2));
    float lambda1 = (a + b - delta) / 2;
    float lambda2 = (a + b + delta) / 2;

    float v1_x = (lambda2 - b) / c;
    float v1_y = 1;
    float v2_x = (lambda1 - b) / c;
    float v2_y = 1;

    return (Mat2){.a11 = v1_x, .a12 = v1_y, .a21 = v2_x, .a22 = v2_y};
}

// ----------------------------
//
// Convert Q, c to center, halfdims and angle
//
// ----------------------------

Ellipsoid2d ellispoid_from_QcForm(Mat2 Q, Vec2 c)
{

    Ellipsoid2d el;
    Vec2 eigenvalues = mat2x2_eigenvalues(Q);
    Mat2 eigenvectors = mat2x2_eigenvectors(Q);

    float angle = atan2(eigenvectors.a22, eigenvectors.a21);
    Vec2 half_dims = {.x = 1 / sqrt(eigenvalues.x), .y = 1 / sqrt(eigenvalues.y)};

    el.center = c;
    el.half_dims = half_dims;
    el.angle = angle * 180 / PI;
    return el;
}

// ----------------------------
//
// The algorithm
//
// ----------------------------

void find_supporting_points(Vec2* points, u32 n_points, Vec2 direction, Vec2* min_point, Vec2* max_point)
{
    float min_projection = INFINITY, max_projection = -INFINITY;
    u32 min_idx = -1, max_idx = -1;

    for (u32 i = 0; i < n_points; i++) {
        float projection = vec2vec2_dot(points[i], direction);
        if (projection < min_projection) {
            min_projection = projection;
            min_idx = i;
        }
        if (projection > max_projection) {
            max_projection = projection;
            max_idx = i;
        }
    }

    *min_point = points[min_idx];
    *max_point = points[max_idx];
}

void initial_volume_approximation(Vec2* aprox, Vec2* points, u32 n_points)
{
    // Generate 4 points
    Vec2 direction1 = vec2_normalize((Vec2){.x = ((float)rand() / (float)RAND_MAX) * 2 - 1,
                                            .y = ((float)rand() / (float)RAND_MAX) * 2 - 1});
    // generate the first 2
    Vec2 min_point1, max_point1;
    find_supporting_points(points, n_points, direction1, &min_point1, &max_point1);
    aprox[0] = min_point1;
    aprox[1] = max_point1;

    // generate the last 2
    Vec2 direction2 = {.x = -direction1.y, .y = direction1.x};
    Vec2 min_point2, max_point2;
    find_supporting_points(points, n_points, direction2, &min_point2, &max_point2);
    aprox[2] = min_point2;
    aprox[3] = max_point2;
}

Ellipsoid2d compute_bounding_ellipsoid_iterative_2d(Ellipsoid2d* ellipsoids, u32 count, float epsilon)
{

    if (count == 0) {
        Ellipsoid2d be;
        return be;
    }

    // 1. Compute the vertices of the bounding boxes for each children ellipsoid
    Vec2 points[count * 4];
    for (u32 i = 0; i < count; i += 1) {
        //  The xmin/xmax and ymin/ymax for a non rotated and centered ellipsoid
        float xmin = -ellipsoids[i].half_dims.x, xmax = +ellipsoids[i].half_dims.x, ymin = -ellipsoids[i].half_dims.y, ymax = +ellipsoids[i].half_dims.y;
        float angle = ellipsoids[i].angle * PI / 180;

        Vec2 v[4] = {(Vec2){xmin, ymin}, (Vec2){xmax, ymin}, (Vec2){xmax, ymax}, (Vec2){xmin, ymax}};
        for (u32 j = 0; j < 4; j += 1) {
            points[4 * i + j].x = cosf(angle) * v[j].x - sinf(angle) * v[j].y + ellipsoids[i].center.x;
            points[4 * i + j].y = sinf(angle) * v[j].x + cosf(angle) * v[j].y + ellipsoids[i].center.y;
        }
    };

    // Algorithm 3.1 of the paper

    // Compute the X0 set
    Vec2 X[1000]; // TODO: actually only is needed 2*d = 4; But we will keep adding... do it dynamic?

    initial_volume_approximation(X, points, ARRAY_SIZE(points));
    // for (u32 i = 0; i < 4; i += 1) {
    //     printf("point %u: {%f, %f}\n", i, X[i].x, X[i].y);
    // }

    // Algorithm 4.1 of the paper (adapted to points)

    // --- Initial ellipse

    // 1. Compute u0
    float u[1000]; // same as X0
    u[0] = 0.25f, u[1] = 0.25f, u[2] = 0.25f, u[3] = 0.25f;

    // 2. Compute w0
    Vec2 w0 = {.x = 0.f, .y = 0.f};
    for (u32 i = 0; i < 4; i += 1) {
        w0.x += X[i].x * u[i];
        w0.y += X[i].y * u[i];
    }
    // printf("w0: {%f, %f}\n", w0.x, w0.y);

    // 3. Compute M0
    Mat2 M0 = {.a11 = 0.f, .a12 = 0.f, .a21 = 0.f, .a22 = 0.f};
    for (u32 i = 0; i < 4; i++) {
        Vec2 v1 = {.x = X[i].x - w0.x, .y = X[i].y - w0.y};
        Mat2 m = mat2x2_muls(vec2vec2_outer(v1, v1), u[i]);
        M0 = mat2x2_sum(M0, m);
    }
    M0 = mat2x2_muls(mat2x2_muls(M0, 2), 0.001f); // trick for the inverse to not get small values.
    M0 = mat2_inverse(M0);                        // The real M should be multiplied by 0.001;
    // printf("M: { {%.8f, %.8f}, {%.8f, %.8f}}\n", M0.a11, M0.a12, M0.a21, M0.a22);

    // 4. Compute x^2n+1. Here we adapt it to compute the farthest point from the center of the ellipse
    float max_distance = -INFINITY, distance = -INFINITY;
    u32 max_idx;
    for (u32 i = 0; i < 4 * count; i += 1) {
        Vec2 v = {.x = points[i].x - w0.x, .y = points[i].y - w0.y};
        distance = vec2vec2_dot(v, mat2x2vec2_mul(M0, v));
        if (distance > max_distance) {
            max_distance = distance;
            max_idx = i;
        }
    }
    Vec2 x2n1 = points[max_idx];
    // printf("x2n1: {%f, %f}\n", x2n1.x, x2n1.y);

    // 5. Add this farthest point to X0
    X[4] = x2n1;
    // TODO: check uniques in X!

    // for (u32 i = 0; i < ARRAY_SIZE(points); i++) {
    //     printf(" X[%u]: (%f, %f)", i, X[i].x, X[i].y);
    // }

    // 6. Compute eps0, to check the error made. This has a problem. And it is the floating point precision.... how to solve it?
    Vec2 v = {.x = x2n1.x - w0.x, .y = x2n1.y - w0.y};
    float epsk = -1 + vec2vec2_dot(v, mat2x2vec2_mul(M0, v)) * 1e-3;
    // printf("eps0: %f\n", epsk); // we get 0.479007 for the test. In Mathematica we get 0.478976

    u32 k = 0;
    Mat2 M = {.a11 = 0.f, .a12 = 0.f, .a21 = 0.f, .a22 = 0.f};
    Vec2 w = {.x = 0.f, .y = 0.f};
    while (1) {
        if (epsk < epsilon)
            break;
        // compute beta
        float beta = epsk / (3 * (1 + epsk));
        // printf("beta: %f\n", beta);
        // update k
        k += 1;
        // update u
        for (u32 i = 0; i < 4 + k; i += 1) {
            u[i] = (1 - beta) * u[i];
        }
        u[4 + k - 1] = beta;

        // for (u32 i = 0; i < 4 + k; i += 1) {
        //     printf(" u%u: %f", i, u[i]);
        // }
        // printf("\n");

        // compute wk
        w.x = 0.f, w.y = 0.f;
        for (u32 i = 0; i < 4 + k; i += 1) {
            w.x += X[i].x * u[i];
            w.y += X[i].y * u[i];
        }
        // printf("w%u: {%f, %f}\n", k, w.x, w.y);

        // comute Mk
        M.a11 = 0.f, M.a12 = 0.f, M.a21 = 0.f, M.a22 = 0.f;
        for (u32 i = 0; i < 4 + k; i++) {
            Vec2 v1 = {.x = X[i].x - w.x, .y = X[i].y - w.y};
            Mat2 m = mat2x2_muls(vec2vec2_outer(v1, v1), u[i]);
            M = mat2x2_sum(M, m);
        }
        M = mat2x2_muls(mat2x2_muls(M, 2), 0.001f); // trick for the inverse to not get small values.
        M = mat2_inverse(M);                        // The real M should be multiplied by 0.001;
        // printf("M%u: { {%.8f, %.8f}, {%.8f, %.8f}}\n", k, M.a11, M.a12, M.a21, M.a22);

        // compute x^2n+1+k. Here we adapt it to compute the farthest point from the center of the ellipse
        float max_distance = -INFINITY, distance = -INFINITY;
        u32 max_idx;
        for (u32 i = 0; i < 4 * count; i += 1) {
            Vec2 v = {.x = points[i].x - w.x, .y = points[i].y - w.y};
            distance = vec2vec2_dot(v, mat2x2vec2_mul(M, v));
            if (distance > max_distance) {
                max_distance = distance;
                max_idx = i;
            }
        }
        Vec2 x2n1k = points[max_idx];
        // printf("x2n1%u: {%f, %f}\n", k, x2n1k.x, x2n1k.y);

        // add the x2n1k to the X set
        X[4 + k] = x2n1k;
        for (u32 i = 0; i < ARRAY_SIZE(points); i++) {
            // printf(" X[%u]: (%f, %f)", i, X[i].x, X[i].y);
        }

        // update the epsk
        Vec2 v = {.x = x2n1k.x - w.x, .y = x2n1k.y - w.y};
        epsk = -1 + vec2vec2_dot(v, mat2x2vec2_mul(M, v)) * 1e-3;
        // printf("eps%u: %f\n", k, epsk); // floating point error

        // printf("---------------------\n");
        // if (k == ARRAY_SIZE(points) + 6)
        //     break;
    }

    M = mat2x2_muls(M, 0.001);
    M = mat2x2_muls(M, 1.f / (1 + epsk));
    // printf("Number of iterations: %u\n", k);
    // printf("epsk: %f\n", epsk);
    // printf("M%u: { {%.8f, %.8f}, {%.8f, %.8f}}\n", k, M.a11, M.a12, M.a21, M.a22);
    // printf("w%u: {%f, %f}\n", k, w.x, w.y);

    Ellipsoid2d bounding_ellipsoid = ellispoid_from_QcForm(M, w);

    // printf("-----------------\n");
    // printf("\tCenter: (%f, %f)\n", bounding_ellipsoid.center.x, bounding_ellipsoid.center.y);
    // printf("\tHalf dims: (%f, %f)\n", bounding_ellipsoid.half_dims.x, bounding_ellipsoid.half_dims.y);
    // printf("\tAngle (rad): %f\n", bounding_ellipsoid.angle / 180 * PI);

    return bounding_ellipsoid;
}
