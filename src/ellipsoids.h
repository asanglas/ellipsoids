#pragma once
#include "stdint.h"

// ----------------------------------
//
//  Structs
//
// ----------------------------------

typedef struct Vec2 Vec2;
struct Vec2 {
    float x, y;
};

typedef struct Mat2 Mat2;
struct Mat2 {
    float a11, a12;
    float a21, a22;
};

typedef struct Mat3 Mat3;
struct Mat3 {
    float a11, a12, a13;
    float a21, a22, a23;
    float a31, a32, a33;
};

typedef struct Ellipsoid2d Ellipsoid2d;
struct Ellipsoid2d {
    Vec2 center;
    Vec2 half_dims;
    float angle;
};

typedef struct Box2d Box2d;
struct Box2d {
    Vec2 center;
    Vec2 half_dims;
    float angle;
};

// ----------------------------------
//
//  Types
//
// ----------------------------------

typedef uint32_t u32;

// ----------------------------------
//
//  Constants
//
// ----------------------------------

#define PI 3.14159265358979323846f

// ----------------------------------
//
//  Macros
//
// ----------------------------------
#ifndef ARRAY_SIZE
#define ARRAY_SIZE(_arr) ((u32)(sizeof(_arr) / sizeof(_arr[0])))
#endif // ARRAY_SIZE

#ifndef MIN
#define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#endif

#ifndef MAX
#define MAX(_a, _b) (((_b) < (_a)) ? (_a) : (_b))
#endif

#ifndef GENERATE_RANDF
#define GENERATE_RANDF(_min, _max) (_min + (float)rand() / (float)(RAND_MAX) * (_max - _min))
#endif