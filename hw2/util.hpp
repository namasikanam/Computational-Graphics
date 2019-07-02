#pragma once
#include <cmath>
#include <cstdio>

struct Vec
{
    double x, y, z;
    Vec(double _x = 0, double _y = 0, double _z = 0) { x = _x, y = _y, z = _z; }
    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec &unit() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
    Vec operator%(Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
    double norm() { return sqrt(x * x + y * y * z * z); }
};

void outVec(Vec v)
{
    printf("(%.6f, %.6f, %.6f)", v.x, v.y, v.z);
}

struct Ray
{
    Vec origin, direction;
    Ray(Vec _origin, Vec _direction) : origin(_origin), direction(_direction) {}
};

enum Reflection_type
{
    DIFFUSE,
    SPECULAR,
    REFRACTIVE,
    SOURCE
};

const double eps = 1e-4;
const double inf = 1e20;

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

inline double BRDF(Vec in_d, Vec out_d) { return fabs(in_d.unit().dot(out_d.unit())); }

inline double distance(Vec a, Vec b) { return (a - b).norm(); }

inline int factorial(int n)
{
    int ans = 1;
    for (int i = n; i; --i)
        ans = ans * i;
    return ans;
}
inline int binomial(int n, int m)
{
    return factorial(n) / factorial(m) / factorial(n - m);
}