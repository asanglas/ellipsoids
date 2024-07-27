#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ellipsoids.h"
#include "ellipsoid_obb_2d.h"

#define DEBUG 0

#ifdef RENDER
#include "render.h"
u32 RENDER_DEBUG = 1;
u32 RENDER_INFO = 1;
#endif

u32 NUM_ELLIPSOIDS = 2;

void generate_random_ellipses(Ellipsoid2d* ellipses, u32 count, u32 seed)
{
    srand(seed);
    for (u32 i = 0; i < count; i += 1) {
        ellipses[i].center = (Vec2){.x = GENERATE_RANDF(-100, 250), .y = GENERATE_RANDF(-250, 250)};
        ellipses[i].half_dims = (Vec2){.x = GENERATE_RANDF(10, 50), .y = GENERATE_RANDF(10, 50)};
        ellipses[i].angle = GENERATE_RANDF(0, 360);
#if DEBUG
        printf("-----------------\n");
        printf("Ellipsoid %u\n", i);
        printf("\tCenter: (%f, %f)\n", ellipses[i].center.x, ellipses[i].center.y);
        printf("\tHalf dims: (%f, %f)\n", ellipses[i].half_dims.x, ellipses[i].half_dims.y);
        printf("\tAngle (rad): %f\n", ellipses[i].angle / 180 * PI);
#endif
    }
}

int main()
{

    u32 seed = 1721908430;

#ifdef RENDER // Draw

    init_raylib();
    while (!WindowShouldClose()) {

        // Keyboard input
        if (IsKeyPressed(KEY_Q))
            RENDER_DEBUG = !RENDER_DEBUG;
        if (IsKeyPressed(KEY_I))
            RENDER_INFO = !RENDER_INFO;
        if (IsKeyPressed(KEY_R)) {
            seed = (u32)rand();
        }
        if (IsKeyPressed(KEY_M)) {
            NUM_ELLIPSOIDS += 1;
        }
        if (IsKeyPressed(KEY_L)) {
            NUM_ELLIPSOIDS = NUM_ELLIPSOIDS > 0 ? NUM_ELLIPSOIDS - 1 : 0;
        }

        Ellipsoid2d ellipses[NUM_ELLIPSOIDS];
        generate_random_ellipses(ellipses, NUM_ELLIPSOIDS, seed);
        // Compute bounding ellipsoid each frame
        float t1 = GetTime();
        u32 count = ARRAY_SIZE(ellipses);
        Ellipsoid2d bounding_ellipsoid = compute_bounding_ellipsoid_obb_2d(ellipses, count);
        float et1 = GetTime() - t1;

        BeginDrawing();
        ClearBackground(RAYWHITE);

        // Draw ellipsoids
        for (u32 i = 0; i < ARRAY_SIZE(ellipses); i++) {
            draw_ellipsoid_2d(ellipses[i], BLACK, 1);
        }
        // Draw bounding ellipsoid
        draw_ellipsoid_2d(bounding_ellipsoid, (Color){200, 122, 255, 255}, 1);

        // Compute the area
        float area = PI * bounding_ellipsoid.half_dims.x * bounding_ellipsoid.half_dims.y;

        // Draw Debug stuff (boxes, etc)
        if (RENDER_DEBUG) {
            for (u32 i = 0; i < ARRAY_SIZE(ellipses); i++) {
                Box2d box;
                box.center = ellipses[i].center;
                box.half_dims = ellipses[i].half_dims;
                box.angle = ellipses[i].angle;
                draw_rect_2d(box, (Color){200, 122, 255, 100}, 1);
            }

            Box2d box;
            box.center = bounding_ellipsoid.center;
            box.half_dims.x = bounding_ellipsoid.half_dims.x / sqrtf(2);
            box.half_dims.y = bounding_ellipsoid.half_dims.y / sqrtf(2);
            box.angle = bounding_ellipsoid.angle;
            draw_rect_2d(box, (Color){200, 122, 255, 100}, 1);
        }

        // Draw UI
        DrawText(TextFormat("DEBUG (Q): "), 10, 10, 20, BLACK);
        DrawText(TextFormat(RENDER_DEBUG ? "ON" : "OFF"), 125, 10, 20, RENDER_DEBUG ? GREEN : RED);
        DrawText(TextFormat("INFO  (I): "), 10, 30, 20, BLACK);
        DrawText(TextFormat(RENDER_INFO ? "ON" : "OFF"), 125, 30, 20, RENDER_INFO ? GREEN : RED);
        DrawText(TextFormat("RANDOMIZE (R): %u", seed), 10, 50, 20, BLACK);
        DrawText(TextFormat("# Ellipsoids (M/L): %u", count), 10, 70, 20, BLACK);

        if (RENDER_INFO) {
            DrawText(TextFormat("--- INFO --- "), 10, 110, 20, BLUE);
            for (u32 method = 0; method < 1; method += 1) {
                u32 offset = 130 + method * 130;
                DrawText(TextFormat("-- Method 1: OBB --"), 10, offset + 10, 20, (Color){200, 122, 255, 255});
                DrawText(TextFormat("C: (%.02f, %.02f)", bounding_ellipsoid.center.x, bounding_ellipsoid.center.y), 10, offset + 30, 20, BLACK);
                DrawText(TextFormat("HD: (%.02f, %.02f)", bounding_ellipsoid.half_dims.x, bounding_ellipsoid.half_dims.y), 10, offset + 50, 20, BLACK);
                DrawText(TextFormat("A: %.02f", bounding_ellipsoid.angle), 10, offset + 70, 20, BLACK);
                DrawText(TextFormat("Area: %.02f", area), 10, offset + 90, 20, BLACK);
                DrawText(TextFormat("ET: %.06f ms", et1 * 1e3), 10, offset + 110, 20, (Color){200, 122, 255, 255});
            }
        }

        EndDrawing();
    }
    CloseWindow();
#endif

    // TESTS

    Vec2 points[20];

    return 0;
}