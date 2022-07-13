#include "v3d.h"

#include <math.h>

v3d normalize (v3d a)
    {
    double length = module(a); 
    
    return length == 0? a : div (a, length);
    }

v3d FromPolar (double a, double thetta, double fi)
    {
    return (v3d){a * sin (thetta) * cos (fi), 
                 a * sin (thetta) * sin (fi),
                 a * cos (thetta)};
    }

v3d add (v3d a, v3d b)
    {
    return (v3d){a.x + b.x, a.y + b.y, a.z + b.z};
    }

v3d sub (v3d a, v3d b)
    {
    return (v3d){a.x - b.x, a.y - b.y, a.z - b.z};
    }

v3d mul (v3d a, double b)
    {
    return (v3d){b * a.x, b * a.y, b * a.z};
    }

v3d div (v3d a, double b)
    {
    return (v3d){a.x / b, a.y / b, a.z / b};
    }

v3d vec_mul (v3d a, v3d b)
    {
    return (v3d){a.y * b.z - a.z * b.y, 
              - (a.x * b.z - a.z * b.x),
                 a.x * b.y - a.y * b.x};
    }

double scl_mul (v3d a, v3d b)
    {
    return a.x * b.x + a.y * b.y + a.z * b.z; 
    }

double module (v3d a)
    {
    return sqrt (a.x * a.x + a.y * a.y + a.z * a.z);
    }

