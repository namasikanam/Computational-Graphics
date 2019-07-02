#pragma once

#include "sphere.hpp"
#include <vector>
#include "bezier.hpp"

Sphere spheres[] = {
    // Scene: radius, position, emission, color, material
    // x in [1, 99], y in [0, 81.6], z in [0, 170]
    Sphere(1, Vec(50, 50, 180), Vec(), Vec(1, 1, 1), SOURCE),                    //光源
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFFUSE),   //Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFFUSE), //Rght
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFFUSE),         //Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFFUSE),               //Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFFUSE),         //Botm
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFFUSE), //Top
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1), SPECULAR),              //Mirr
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1), REFRACTIVE),            //Glas
    Sphere(16.5, Vec(20, 60, 20), Vec(), Vec(1, 1, 1), DIFFUSE)                  //再加个球看看效果
};

Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).unit());
Vec source(50, 50, 180); // emm在内部丢了一个点光源

inline bool scene_intersect(const Ray &ray, double &t, int &id)
{
    double d;
    t = inf;
    for (int i = sizeof(spheres) / sizeof(Sphere); i--;)
        if ((d = spheres[i].intersect(ray)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < inf;
}

inline bool scene_photon_intersect(const Ray &ray, double &t, int &id)
{
    double d;
    t = inf;
    for (int i = sizeof(spheres) / sizeof(Sphere); --i;)
        if ((d = spheres[i].intersect(ray)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < inf;
}

inline double sphere_and_bezier_distance(Ray r)
{
    double t;
    int id;
    scene_intersect(r, t, id);
    double s;
    Vec x, normal;
    bezier_intersect(r, s, x, normal);
    return std::min(t, s);
}

inline bool sphere_and_bezier_intersect(Ray r, Vec &x, Vec &n, Vec &color, Reflection_type &refl, bool photon_flag)
{
    // x = r.origin + r.direction * t;
    // n = (x - obj.position).unit();

    // printf("> sphere_and_bezier_intersect\n");

    double sphere_t, bezier_t;
    int id;
    bool intersect_with_sphere = photon_flag ? scene_photon_intersect(r, sphere_t, id) : scene_intersect(r, sphere_t, id);
    bool intersect_with_bezier = bezier_intersect(r, bezier_t, x, n);
    // bool intersect_with_bezier = false;

    // printf("< sphere_and_bezier_intersect\n");

    if (!intersect_with_sphere && !intersect_with_bezier)
        return 0;
    else if (intersect_with_sphere && (!intersect_with_bezier || sphere_t < bezier_t))
    // else if (true)
    {
        // if (photon_flag && intersect_with_bezier)
        // {
        //     printf("=================\n");
        //     printf("origin = "), outVec(r.origin), puts("");
        //     printf("direction = "), outVec(r.direction), puts("");
        //     printf("sphere_t = %.6f\n", sphere_t);
        //     printf("bezier_t = %.6f\n", bezier_t);
        // }

        Sphere &obj = spheres[id];
        x = r.origin + r.direction * sphere_t;
        n = (x - obj.position).unit();
        color = obj.color;
        refl = obj.reflection;
        return 1;
    }
    else
    {
        color = Vec(1, 1, 1);
        refl = DIFFUSE;

        // {
        //     static int tot = 0;
        //     if (photon_flag)
        //         printf("Intersect with bezier!\n");
        // }

        return 1;
    }
}