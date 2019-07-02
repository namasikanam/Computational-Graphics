#pragma once

#include <vector>
#include "util.hpp"
#include <algorithm>

struct Bezier
{
    std::vector<std::pair<double, double>> p;
    std::vector<double> x, y;
    double x0 = 40, z0 = 90;
    int n = 3;

    Bezier()
    {
        p.push_back(std::make_pair(5.808, 0));
        p.push_back(std::make_pair(15.895, 12.141));
        p.push_back(std::make_pair(0.741, 21.708));
        p.push_back(std::make_pair(6.27, 28.706));

        x.resize(n + 1), y.resize(n + 1);
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n - i; ++j)
            {
                int coe = binomial(n, i) * binomial(n - i, j) * (j & 1 ? -1 : 1);
                x[i + j] += coe * p[i].first;
                y[i + j] += coe * p[i].second;
            }
    }
    double getX(double t)
    {
        double ans = 0;
        for (int i = n; i >= 0; --i)
            ans = ans * t + x[i];
        return ans;
    }
    double getdX(double t)
    {
        double ans = 0;
        for (int i = n; i; --i)
            ans = ans * t + i * x[i];
        return ans;
    }
    double getY(double t)
    {
        double ans = 0;
        for (int i = n; i >= 0; --i)
            ans = ans * t + y[i];
        return ans;
    }
    double getdY(double t)
    {
        double ans = 0;
        for (int i = n; i; --i)
            ans = ans * t + i * y[i];
        return ans;
    }
};

Bezier bezier;

bool bezier_intersect(Ray r, double &s, Vec &x, Vec &normal)
{
    if (fabs(r.direction.y) < eps) // 战略放弃掉与y轴垂直的光线了……
        return 0;

    // 首先算出多项式
    const int n = 4, deg = 2 * n + 1;
    double coe[deg] = {0}; // 多项式的系数
    double a = r.direction.x / r.direction.y, c = r.direction.z / r.direction.y;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            coe[i + j] += (a * a + c * c) * bezier.y[i] * bezier.y[j] - bezier.x[i] * bezier.x[j];
    for (int i = 0; i <= n; ++i)
        coe[i] += 2 * (a * (r.origin.x - bezier.x0 - a * r.origin.y) + c * (r.origin.z - bezier.z0 - c * r.origin.y)) * bezier.y[i];
    coe[0] += pow(r.origin.x - bezier.x0 - a * r.origin.y, 2) + pow(r.origin.z - bezier.z0 - c * r.origin.y, 2);

    // 牛顿迭代
    // 凡是算不出来的情况都统统放弃吧……
    auto f = [coe](double x) {
        double ans = 0;
        for (int i = deg; i--;)
            ans = ans * x + coe[i];
        return ans;
    };
    auto df = [coe](double x) {
        double ans = 0;
        for (int i = deg; --i;)
            ans = ans * x + i * coe[i];
        return ans;
    };
    auto newtonMethod = [f, df, coe](double t0, double &t) -> bool {
        if (fabs(f(t0)) < eps)
        {
            t = t0;
            return 1;
        }
        if (fabs(df(t0)) < eps)
            return 0;

        int restDepth = 100;
        double t1 = t0 - f(t0) / df(t0);
        while (--restDepth && fabs(f(t1)) > eps && fabs(t0 - t1) > eps)
        {
            t0 = t1;
            t1 = t1 - f(t1) / df(t1);
        }
        if (restDepth == 0)
            return 0;
        else
        {
            t = t1;
            return t >= 0 && t <= 1;
        }
    };
    double t, t0, t1;
    double y, y0, y1;
    if (!newtonMethod(0.5, t0))
        return 0;
    y0 = bezier.getY(t0);

    if (f(0) * f(1) > 0)
    { // 有两个解，需要在里面选择一个入射的。
        auto binSolve = [f, coe](double l, double r) {
            double m;
            while (r - l > eps)
            {
                m = (l + r) / 2;
                if (f(m) * f(r) < 0)
                    l = m;
                else
                    r = m;
            }
            return m;
        };
        if (f(t0 + eps) * f(1) < 0)
            t1 = binSolve(t0 + eps, 1);
        else
            t1 = binSolve(0, t0 - eps);
        y1 = bezier.getY(t1);
        if (y0 > y1)
        {
            std::swap(y0, y1);
            std::swap(t0, t1);
        }

        // printf("y0 = %.6f, y1 = %.6f\n", y0, y1);

        if (r.direction.y > 0)
            if (r.origin.y < y0)
                t = t0, y = y0;
            else if (r.origin.y < y1)
                t = t1, y = y1;
            else
                return 0;
        else if (r.origin.y > t1)
            t = t1, y = y1;
        else if (r.origin.y > t0)
            t = t0, y = y0;
        else
            return 0;
    }
    else
        t = t0, y = y0;

    s = (y - r.origin.y) / r.direction.y;
    if (s < 0)
        return 0;

    x = r.origin + r.direction * s;
    normal = Vec(bezier.getdY(t) * (x.x - bezier.x0), -bezier.getX(t) * bezier.getdX(t), bezier.getdY(t) * (x.z - bezier.z0)).unit();

    // printf("=================\n");
    // printf("origin = "), outVec(r.origin), puts("");
    // printf("direction = "), outVec(r.direction), puts("");
    // printf("intersection = "), outVec(x), puts("");
    // printf("normal = "), outVec(normal), puts("");

    // printf("Find a intersection!: s = %.6f\n", s);

    return 1;
}