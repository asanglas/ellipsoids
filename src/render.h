#pragma once
#include "stdio.h"
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"

#include "ellipsoids.h"

// ----------------------------------
//
//  Raylib
//
// ----------------------------------

u32 screen_width = 1000;
u32 screen_height = 1000;
const char* title = "Ellipsoids";

void init_raylib()
{
    InitWindow(screen_width, screen_height, title);
    SetTargetFPS(60);
}

void render_frame()
{
}

// ----------------------------------
//
// Drawing functions
//
// ----------------------------------

void draw_ellipsoid_2d(Ellipsoid2d e, Color color, u32 jump)
{
    if (jump <= 0)
        return;

    int segments = 360;
    int step = 5;
    float dt = 2 * PI / segments;

    // Set line width
    rlEnableTexture(0); // Ensure no texture is enabled
    rlSetLineWidth(2.0f);
    rlBegin(RL_LINES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    for (u32 i = 0; i <= segments; i += jump * step) {

        float t1 = i * dt;
        float t2 = (i + step) * dt;
        float angle = e.angle * PI / 180;

        float x1 = e.center.x + e.half_dims.x * cosf(t1) * cosf(angle) - e.half_dims.y * sinf(t1) * sinf(angle);
        float y1 = e.center.y + e.half_dims.y * cosf(angle) * sinf(t1) + e.half_dims.x * cosf(t1) * sinf(angle);

        float x2 = e.center.x + e.half_dims.x * cosf(t2) * cosf(angle) - e.half_dims.y * sinf(t2) * sinf(angle);
        float y2 = e.center.y + e.half_dims.y * cosf(angle) * sinf(t2) + e.half_dims.x * cosf(t2) * sinf(angle);

        // to screen (the orgin is at the middle)
        rlVertex2f(x1 + screen_width / 2, -y1 + screen_height / 2);
        rlVertex2f(x2 + screen_width / 2, -y2 + screen_height / 2);
    }
    rlEnd();
}

void draw_rect_2d(Box2d b, Color color, u32 jump)
{
    if (jump <= 0)
        return;

    int segments = 100;
    int step = 5;
    float dt = 1.0 / segments;

    // Set line width
    rlEnableTexture(0); // Ensure no texture is enabled
    rlSetLineWidth(2.0f);
    rlBegin(RL_LINES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    Vec2 v[4];
    v[0] = (Vec2){.x = -b.half_dims.x, .y = -b.half_dims.y};
    v[1] = (Vec2){.x = b.half_dims.x, .y = -b.half_dims.y};
    v[2] = (Vec2){.x = b.half_dims.x, .y = b.half_dims.y};
    v[3] = (Vec2){.x = -b.half_dims.x, .y = b.half_dims.y};
    float angle = b.angle * PI / 180;

    for (u32 i = 0; i < 4; i += 1) {
        v[i] = (Vec2){
            .x = cosf(angle) * v[i].x - sinf(angle) * v[i].y + b.center.x,
            .y = sinf(angle) * v[i].x + cosf(angle) * v[i].y + b.center.y,
        };
    }

    for (u32 i = 0; i < 4; i += 1) {
        for (u32 j = 0; j < segments; j += jump * step) {
            float t1 = j * dt;
            float t2 = (j + step) * dt;
            float x1 = v[i].x + t1 * (v[(i + 1) % 4].x - v[i].x);
            float y1 = v[i].y + t1 * (v[(i + 1) % 4].y - v[i].y);
            float x2 = v[i].x + t2 * (v[(i + 1) % 4].x - v[i].x);
            float y2 = v[i].y + t2 * (v[(i + 1) % 4].y - v[i].y);
            rlVertex2f(x1 + screen_width / 2, -y1 + screen_height / 2);
            rlVertex2f(x2 + screen_width / 2, -y2 + screen_height / 2);
        }
    }
    rlEnd();
}
