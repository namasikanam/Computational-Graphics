#pragma once
#include "scene.hpp"
#include "sphere.hpp"
#include "util.hpp"
#include <cmath>
#include <vector>
#include <cstdio>

struct HitPoint
{
    Vec pos;
    Vec dir;
    Vec coe;   // 贡献系数。
    int pixel; // 贡献给谁。
    double r;
    double n; // 虽然paper里用的是integer，但我觉得是double比较合适，更加精确。
    Vec tau;
    bool source;

    int m = 0;
    Vec tauM = Vec(0, 0, 0);
};

void outHitPoint(HitPoint p)
{
    printf("{\n");
    printf("\tpos = "), outVec(p.pos), puts("");
    printf("\tdir = "), outVec(p.dir), puts("");
    printf("\tcoe = "), outVec(p.coe), puts("");
    printf("\tpixel = %d, r = %.6f, n = %.6f\n", p.pixel, p.r, p.n);
    printf("\ttau = "), outVec(p.tau), puts("");
    puts("}");
}

void rayTracing(std::vector<HitPoint> &hitPoints, Vec coe, int pixel, const Ray &r, int depth, unsigned short *Xi)
{
    if (coe.norm() < eps)
        return;

    Vec x, n, f;
    Reflection_type refl;
    if (!sphere_and_bezier_intersect(r, x, n, f, refl, 0))
        return; // if miss, return black

    if (refl == SOURCE)
    {
        HitPoint hitPoint;
        hitPoint.source = true, hitPoint.coe = coe, hitPoint.pixel = pixel;
        hitPoints.push_back(hitPoint);
    }

    bool into = n.dot(r.direction) < 0;
    Vec nl = into ? n : n * -1;           // 与ray夹角为钝角的那个法向量。
    double p = fmax(f.x, fmax(f.y, f.z)); // 直接以最大值作为RR的概率，其实蛮让人困惑的。
    if (++depth > 5)
        if (erand48(Xi) < p)
            f = f * (1 / p);
        else
            return; // 我不太确定这里是放一个hit point比较好还是直接return比较好。
                    // 但我觉得，在理想情况下，ray tracing和photon mapping应该是完全对称的。
    if (refl == DIFFUSE)
    {
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).unit(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).unit(); // 我还并未理解Sampling Unit Hemisphere，但就暂且相信它的正确性吧！
        hitPoints.push_back(HitPoint({x, r.direction, coe, pixel}));

        // outHitPoint(hitPoints.back());
        // {
        //     static int tot = 0;
        //     if (++tot % 10000 == 0)
        //         printf("the %dth hitPoint!\n", tot);
        // }

        rayTracing(hitPoints, f.mult(coe), pixel, Ray(x, d), depth, Xi);
    }
    else if (refl == SPECULAR)
    {
        rayTracing(hitPoints, f.mult(coe), pixel, Ray(x, r.direction - n * 2 * n.dot(r.direction)), depth, Xi);
    }
    else
    {
        Ray reflection_ray(x, r.direction - n * 2 * n.dot(r.direction));
        double nc = 1, nt = 1.5; // nt是折射率，其实是一个可以放到scene里的参数。
        double nnt = into ? nc / nt : nt / nc, coss = -r.direction.dot(nl), ddn = r.direction.dot(nl), cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - coss * coss)) < 0)
            rayTracing(hitPoints, f.mult(coe), pixel, reflection_ray, depth, Xi);
        Vec tdir = (r.direction * nnt - nl * (ddn * nnt + sqrt(cos2t))).unit();
        double R0 = (nt - nc) * (nt - nc) / (nt + nc) / (nt + nc);
        double c = 1 + (into ? ddn : -tdir.dot(n)); // 这个c也蛮奇怪的，我并未看到解释。
        double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re;
        p = .25 + .5 * Re;
        double Rp = Re / p, Tp = Tr / (1 - p);
        if (depth > 2)
            if (erand48(Xi) < p)
                rayTracing(hitPoints, f.mult(coe) * Rp, pixel, reflection_ray, depth, Xi);
            else
                rayTracing(hitPoints, f.mult(coe) * Tp, pixel, Ray(x, tdir), depth, Xi);
        else
        {
            rayTracing(hitPoints, f.mult(coe) * Re, pixel, reflection_ray, depth, Xi);
            rayTracing(hitPoints, f.mult(coe) * Tr, pixel, Ray(x, tdir), depth, Xi);
        }
    }
}

std::vector<HitPoint> hitPointMapping(int w, int h, int samps, double *depth_of_field) // 如果是PPM的话，其实根本撑不了很大的samps。
{
    std::vector<HitPoint> hitPoints;
    Vec cx = Vec(.5135 * w / h), cy = (cx % cam.direction).unit() * .5135, r;
    for (int y = 0; y < h; y++)
    {
        // printf("======= y=%d ==========\n", y);
        for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y}; x < w; x++)
        { // Loop cols
            int i = (h - y - 1) * w + x;
            for (int sy = 0; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec())
                { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.direction;
                        rayTracing(hitPoints, Vec(1, 1, 1) * (.25 / samps), i, Ray(cam.origin + d * 140, d.unit()), 0, Xi);

                        depth_of_field[i] += sphere_and_bezier_distance(Ray(cam.origin + d * 140, d.unit())) / 4;
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    // 这里的clamp真的蛮奇怪，radiance出来可能在[0, 1]之外么？
                    // 把这里的clamp删掉了，之后别忘了clamp！
                }
        }
    }
    return hitPoints;
}