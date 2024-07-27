#pragma once
#include "stdio.h"
#include "math.h"
#include "ellipsoids.h"

Ellipsoid2d compute_bounding_ellipsoid_obb_2d(Ellipsoid2d* ellipsoids, u32 count)
{
    Ellipsoid2d bounding_ellipsoid;

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

    // 2. Compute covariance matrix of the vertices
    float sum_x = 0, sum_y = 0, sum2_x = 0, sum2_y = 0, sum2_xy = 0;
    for (u32 i = 0; i < (4 * count); i++) {
        sum_x += points[i].x;
        sum_y += points[i].y;
        sum2_x += points[i].x * points[i].x;
        sum2_y += points[i].y * points[i].y;
        sum2_xy += points[i].x * points[i].y;
    }
    float cov_xx = sum2_x / (4 * count) - (sum_x / (4 * count)) * (sum_x / (4 * count));
    float cov_yy = sum2_y / (4 * count) - (sum_y / (4 * count)) * (sum_y / (4 * count));
    float cov_xy = sum2_xy / (4 * count) - (sum_x / (4 * count)) * (sum_y / (4 * count));

    // 3. Calculate eigenvectors (we follow the notation of )
    float a = cov_xx, b = cov_yy, c = cov_xy;
    float delta = sqrtf(4 * pow(c, 2) + pow(a - b, 2));
    float lambda1 = (a + b - delta) / 2;
    float lambda2 = (a + b + delta) / 2;

    float v1_x = (lambda2 - b) / c;
    float v1_y = 1;
    float v2_x = (lambda1 - b) / c;
    float v2_y = 1;

    // 4. Project the points to the base given by the eigenvectors. Simultaneously, we store the minimum/maximum projected x and y. These will be the corners of the bounding box
    // The projection matrix P is given by the eigenvectors in columns.
    float p11 = v1_x, p12 = v1_y, p21 = v2_x, p22 = v2_y;
    float xmin_proj = INFINITY, ymin_proj = INFINITY, xmax_proj = -INFINITY, ymax_proj = -INFINITY;
    for (u32 i = 0; i < count * 4; i += 1) {
        // initial points
        float x0 = points[i].x, y0 = points[i].y;

        // Project. For a point x, the projection is given by P^T x (T = transposed)
        float xproj = p11 * x0 + p12 * y0;
        float yproj = p21 * x0 + p22 * y0;

        // Check if that point is the min/max
        xmin_proj = MIN(xproj, xmin_proj);
        ymin_proj = MIN(yproj, ymin_proj);
        xmax_proj = MAX(xproj, xmax_proj);
        ymax_proj = MAX(yproj, ymax_proj);
    }

    // 5. Project back the 4 vertices of the box. The inverse projection is given by (P^-1)x
    Vec2 vp[4] = {(Vec2){xmin_proj, ymin_proj}, (Vec2){xmax_proj, ymin_proj}, (Vec2){xmax_proj, ymax_proj}, (Vec2){xmin_proj, ymax_proj}};
    // The inverse projection matrix is given by
    float detp = p11 * p22 - p12 * p21;
    float ip11 = p22 / detp, ip12 = -p12 / detp, ip21 = -p21 / detp, ip22 = p11 / detp;
    // Project
    Vec2 vbox[4];
    for (u32 i = 0; i < 4; i += 1) {
        vbox[i].x = ip11 * vp[i].x + ip12 * vp[i].y;
        vbox[i].y = ip21 * vp[i].x + ip22 * vp[i].y;
    }

    // 6. The halfdims of the projected back box will be the distance between vbox[1]-vbox[0] and vbox[3]-vbox[0], over 2
    Vec2 bbox_half_dims;
    bbox_half_dims.x = sqrtf(pow(vbox[1].x - vbox[0].x, 2) + pow(vbox[1].y - vbox[0].y, 2)) / 2;
    bbox_half_dims.y = sqrtf(pow(vbox[3].x - vbox[0].x, 2) + pow(vbox[3].y - vbox[0].y, 2)) / 2;

    // 7. The bounding ellipsoid
    bounding_ellipsoid.half_dims.x = bbox_half_dims.x * sqrtf(2); // the half dim * sqrt(2)
    bounding_ellipsoid.half_dims.y = bbox_half_dims.y * sqrtf(2);
    bounding_ellipsoid.angle = atan2(1, v1_x) * 180 / PI;                              // the angle of the new v1 axis with respect the (1, 0) axis
    bounding_ellipsoid.center.x = (vbox[0].x + vbox[1].x + vbox[2].x + vbox[3].x) / 4; // the mean of the box vertices
    bounding_ellipsoid.center.y = (vbox[0].y + vbox[1].y + vbox[2].y + vbox[3].y) / 4;

    return bounding_ellipsoid;
}
