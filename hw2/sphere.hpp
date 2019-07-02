#pragma once
#include "util.hpp"

struct Sphere
{
    double radius;
    Vec position, emission, color;
    // 我不确定这个emission是否是必要的。
    // 之后我可能会添加光源，将物体本身定义为不发光的。
    Reflection_type reflection;

    Sphere(double _radius, Vec _position, Vec _emission, Vec _color, Reflection_type _reflection) : radius(_radius), position(_position), emission(_emission), color(_color), reflection(_reflection) {}
    double intersect(const Ray &ray) const
    {
        Vec op = position - ray.origin;
        double b = op.dot(ray.direction), det = b * b - op.dot(op) + radius * radius;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        double t;
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};