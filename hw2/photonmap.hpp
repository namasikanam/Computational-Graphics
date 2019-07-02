#pragma once

#include "util.hpp"
#include "scene.hpp"
#include <vector>
#include <cstdlib>
#include <algorithm>

struct Photon
{
    Vec pos, dir;
    Vec color; // photon power of every color
};

void outPhoton(Photon p)
{
    printf("{\n");
    printf("\tpos = "), outVec(p.pos), puts("");
    printf("\tdir = "), outVec(p.dir), puts("");
    printf("\tcolor = "), outVec(p.color), puts("");
    puts("}");
}

void photonTracing(std::vector<Photon> &photons, Vec color, const Ray &r, int depth, unsigned short *Xi)
{
    if (color.norm() < eps)
        return;

    Vec x, n, f;
    Reflection_type refl;
    if (!sphere_and_bezier_intersect(r, x, n, f, refl, 1)) // 总之是要求出一个交点。
        return;

    bool into = n.dot(r.direction) < 0;
    Vec nl = into ? n : n * -1; // 与ray夹角为钝角的那个法向量。
    double p = fmax(f.x, fmax(f.y, f.z));
    if (++depth > 10)
        if (erand48(Xi) < p)
            f = f * (1 / p);
        else
            return;
    if (refl == DIFFUSE)
    {
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).unit(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).unit();
        photons.push_back(Photon({x, r.direction, color}));          // 类比于Yuzheng Gu的做法，只有在漫反射表面才有光子停留。
        photonTracing(photons, f.mult(color), Ray(x, d), depth, Xi); // 我看Yuzheng Gu是这么写的，不过我自己并不是太理解。
    }
    else if (refl == SPECULAR)
    {
        photonTracing(photons, f.mult(color), Ray(x, r.direction - n * 2 * n.dot(r.direction)), depth, Xi);
    }
    else
    {
        Ray reflection_ray(x, r.direction - n * 2 * n.dot(r.direction));
        double nc = 1, nt = 1.5; // nt是折射率，其实是一个可以放到scene里的参数。
        double nnt = into ? nc / nt : nt / nc, coss = -r.direction.dot(nl), ddn = r.direction.dot(nl), cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - coss * coss)) < 0)
        {
            photonTracing(photons, f.mult(color), reflection_ray, depth, Xi);
            return;
        }
        Vec tdir = (r.direction * nnt - nl * (ddn * nnt + sqrt(cos2t))).unit();
        double R0 = (nt - nc) * (nt - nc) / (nt + nc) / (nt + nc);
        double c = 1 + (into ? ddn : -tdir.dot(n)); // 这个c也蛮奇怪的，我并未看到解释。
        double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re;
        p = .25 + .5 * Re;
        double Rp = Re / p, Tp = Tr / (1 - p);
        if (depth > 2)
            if (erand48(Xi) < p)
                photonTracing(photons, f.mult(color) * Rp, reflection_ray, depth, Xi);
            else
                photonTracing(photons, f.mult(color) * Tp, Ray(x, tdir), depth, Xi);
        else
        {
            photonTracing(photons, f.mult(color) * Re, reflection_ray, depth, Xi);
            photonTracing(photons, f.mult(color) * Tr, Ray(x, tdir), depth, Xi);
        }
    }
}

struct PhotonMap
{
    std::vector<std::vector<std::vector<std::vector<Photon>>>> photonList; // size * size * size * arbitrary
    int size;
    std::vector<int> x, y, z;
};

PhotonMap photonMapping(int number_of_photons, Vec source, unsigned short *Xi) // 从一个点光源发射光子。
{
    std::vector<Photon> photons;
    while (number_of_photons--)
    {
        Vec d(erand48(Xi), erand48(Xi), erand48(Xi)); // 这里随机的可能不是很正确，我还没太想明白如何正确地随机
        // 这里直接做如此粗暴的随机很可能会影响效率，是一个可以改进的地方。

        // {
        //     static int tot = 0;
        //     if (tot++ % 1000 == 0)
        //         printf("%dth photon tracing\n", tot);
        // }

        photonTracing(photons, Vec(1, 1, 1), Ray(source, d.unit()), 0, Xi);
    }
    PhotonMap photonMap;
    int V = photons.size();
    int block_size = pow(V, 2.0 / 3);
    photonMap.size = (V - 1) / block_size + 1;

    printf("V = %d\n", V);
    // for (int i = 3; i--;)
    // {
    //     int x = rand() % V;
    //     outPhoton(photons[x]);
    // }

    // 不妨假设每一个光子的各维坐标都是互不相同的。
    std::vector<int> photon_id;
    std::vector<int> photon_i, photon_j, photon_k;
    for (int i = V; i--;)
        photon_id.push_back(i);
    photon_i.resize(V), photon_j.resize(V), photon_k.resize(V);

    std::sort(photon_id.begin(), photon_id.end(), [&photons](int u, int v) { return photons[u].pos.x < photons[v].pos.x; });
    for (int i = 0; i < V; ++i)
    {
        if (i % block_size == 0)
            photonMap.x.push_back(photons[photon_id[i]].pos.x);
        photon_i[photon_id[i]] = photonMap.x.size() - 1;
    }
    std::sort(photon_id.begin(), photon_id.end(), [&photons](int u, int v) { return photons[u].pos.y < photons[v].pos.y; });
    for (int i = 0; i < V; ++i)
    {
        if (i % block_size == 0)
            photonMap.y.push_back(photons[photon_id[i]].pos.y);
        photon_j[photon_id[i]] = photonMap.y.size() - 1;
    }
    std::sort(photon_id.begin(), photon_id.end(), [&photons](int u, int v) { return photons[u].pos.z < photons[v].pos.z; });
    for (int i = 0; i < V; ++i)
    {
        if (i % block_size == 0)
            photonMap.z.push_back(photons[photon_id[i]].pos.z);
        photon_k[photon_id[i]] = photonMap.z.size() - 1;
    }

    photonMap.photonList.resize(photonMap.size);
    for (int i = photonMap.size; i--;)
    {
        photonMap.photonList[i].resize(photonMap.size);
        for (int j = photonMap.size; j--;)
            photonMap.photonList[i][j].resize(photonMap.size);
    }
    for (int i = 0; i < photons.size(); ++i)
        photonMap.photonList[photon_i[i]][photon_j[i]][photon_k[i]].push_back(photons[i]);

    return photonMap;
}