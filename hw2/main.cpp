#include "raytrace.hpp"
#include "photonmap.hpp"
#include <cstdio>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <ctime>
#include <functional>

const int number_of_photons = 1e6;
const int number_of_photon_mappings = 1000;
const double alpha = 0.7; // 看到paper里提到一次实验中用了这个参数……
const double init_r = 1;  // 感觉这个大小应该还可以？
const double beta = 0.3;

int main(int argc, char *argv[])
{
    // freopen("debug.out", "w", stdout);

    int w = 1024, h = 768, samps = 1;
    Vec *c = new Vec[w * h];
    double *depth_of_field = new double[w * h];

    std::vector<HitPoint> hitPoints;
    hitPoints = hitPointMapping(w, h, samps, depth_of_field);

    printf("number of hitPoints = %d\n", hitPoints.size());
    // for (int i = 3; i--;)
    // {
    //     int x = rand() % hitPoints.size();
    //     outHitPoint(hitPoints[x]);
    // }

    // 从HitPoint找PhotonMapping的初始半径要多大才合适？不知道，过会儿看看Yuzheng Gu怎么做的。
    for (auto &hitPoint : hitPoints)
    {
        hitPoint.r = init_r;
        hitPoint.tau = Vec(0, 0, 0);
    }
    unsigned short Xi[3] = {0, 0, 1};
    for (int _number_of_photon_mappings = 1; _number_of_photon_mappings <= number_of_photon_mappings; ++_number_of_photon_mappings)
    {
        PhotonMap photonMap = photonMapping(number_of_photons, source, Xi);

        printf("size of photonMap = %d\n", photonMap.size);

        for (auto &hitPoint : hitPoints)
            if (!hitPoint.source)
            {
                // 求M, tauM
                int xl = std::max((int)(std::upper_bound(photonMap.x.begin(), photonMap.x.end(), hitPoint.pos.x - hitPoint.r) - photonMap.x.begin() - 1), 0);
                int xr = std::lower_bound(photonMap.x.begin(), photonMap.x.end(), hitPoint.pos.x + hitPoint.r) - photonMap.x.begin();
                int yl = std::max((int)(std::upper_bound(photonMap.y.begin(), photonMap.y.end(), hitPoint.pos.y - hitPoint.r) - photonMap.y.begin() - 1), 0);
                int yr = std::lower_bound(photonMap.y.begin(), photonMap.y.end(), hitPoint.pos.y + hitPoint.r) - photonMap.y.begin();
                int zl = std::max((int)(std::upper_bound(photonMap.z.begin(), photonMap.z.end(), hitPoint.pos.z - hitPoint.r) - photonMap.z.begin() - 1), 0);
                int zr = std::lower_bound(photonMap.z.begin(), photonMap.z.end(), hitPoint.pos.z + hitPoint.r) - photonMap.z.begin();

                // outHitPoint(hitPoint);
                // printf("x in [%d, %d), y in [%d, %d), z in [%d, %d)\n", xl, xr, yl, yr, zl, zr);

                for (int i = xl; i < xr; ++i)
                    for (int j = yl; j < yr; ++j)
                        for (int k = zl; k < zr; ++k)
                        {
                            for (Photon photon : photonMap.photonList.at(i).at(j).at(k))
                            {
                                if (distance(hitPoint.pos, photon.pos) < hitPoint.r)
                                {
                                    ++hitPoint.m;
                                    hitPoint.tauM = hitPoint.tauM + photon.color * BRDF(hitPoint.dir, photon.dir);

                                    // {
                                    //     static int tot = 0;
                                    //     printf("Photon:\n");
                                    //     outPhoton(photon);
                                    // }
                                }
                            }
                        }

                // 更新R, N, tau

                // printf("M = %d\n", M);
                // printf("tauM = "), outVec(tauM), puts("");

                if (hitPoint.m >= 10)
                {
                    hitPoint.r *= sqrt((hitPoint.n + alpha * hitPoint.m) / (hitPoint.n + hitPoint.m));
                    hitPoint.tau = (hitPoint.tau + hitPoint.tauM) * ((hitPoint.n + alpha * hitPoint.m) / (hitPoint.n + hitPoint.m));
                    hitPoint.n = hitPoint.n + alpha * hitPoint.m;

                    hitPoint.m = 0;
                    hitPoint.tauM = Vec(0, 0, 0);
                }
            }

        printf("The %dth photon mapping.\n", _number_of_photon_mappings);
        printf("Time = %.6f\n", (double)clock() / CLOCKS_PER_SEC);

        memset(c, 0, sizeof(Vec) * w * h);

        for (HitPoint hitPoint : hitPoints)
            if (hitPoint.source)
                c[hitPoint.pixel] = c[hitPoint.pixel] + hitPoint.coe.mult(Vec(1, 1, 1));
            else
                c[hitPoint.pixel] = c[hitPoint.pixel] + hitPoint.coe.mult(hitPoint.tau) * (1 / M_PI / hitPoint.r / hitPoint.r);

        // 处理景深
        // {
        //     double *td = new double[w * h];
        //     memcpy(td, depth_of_field, sizeof(double) * w * h);
        //     std::sort(td, td + w * h, std::less<double>());
        //     Vec *new_c = new Vec[w * h];
        //     for (int i = w * h; i--;)
        //         new_c[i] = c[i];

        //     for (int k = 1, _k = 3; k < _k; ++k)
        //         for (int x = k; x < w - k; ++x)
        //             for (int y = k; y < h - k; ++y)
        //             {
        //                 int i = (h - y - 1) * w + x;
        //                 if (depth_of_field[i] < td[w * h * k / _k] && (k == _k - 1 || depth_of_field[i] > td[w * h * (k + 1) / 10]))
        //                 {
        //                     new_c[i] = 0;
        //                     for (int tx = x - k; tx <= x + k; ++tx)
        //                         for (int ty = y - k; ty <= y + k; ++ty)
        //                         {
        //                             int j = (h - ty - 1) * w + tx;
        //                             new_c[i] = new_c[i] + c[j];
        //                         }
        //                     new_c[i] = new_c[i] * (1.0 / (2 * k + 1) / (2 * k + 1));
        //                 }
        //             }
        //     for (int i = w * h; i--;)
        //         c[i] = new_c[i];
        // }

        // 调一下亮度
        double xmax = 0, ymax = 0, zmax = 0;
        for (int i = w * h; i--;)
        {
            xmax = std::max(xmax, c[i].x);
            ymax = std::max(ymax, c[i].y);
            zmax = std::max(zmax, c[i].z);
        }
        xmax *= beta, ymax *= beta, zmax *= beta;
        for (int i = w * h; i--;)
        {
            c[i].x /= xmax;
            c[i].y /= ymax;
            c[i].z /= zmax;
        }

        char debug_filename[100];
        sprintf(debug_filename, "debug%d.out", _number_of_photon_mappings);
        FILE *fdebug = fopen(debug_filename, "w");
        for (int i = 0; i < w * h; ++i)
            fprintf(fdebug, "%.10f %.10f %.10f\n", c[i].x, c[i].y, c[i].z);
        fclose(fdebug);

        char image_filename[100];
        sprintf(image_filename, "image%d.ppm", _number_of_photon_mappings);
        FILE *f = fopen(image_filename, "w"); // Write image to PPM file.
        fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
        for (int i = 0; i < w * h; i++)
            fprintf(f, "%d %d %d ", toInt(clamp(c[i].x)), toInt(clamp(c[i].y)), toInt(clamp(c[i].z)));
        fclose(f);
    }
}