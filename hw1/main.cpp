#include <iostream>
#include <opencv2/opencv.hpp>
#include <bits/stdc++.h>
using namespace cv;
using namespace std;

const int width = 300;
const int height = 300;
const int n = 17;

class Color
{
  public:
    int b, g, r;
    bool operator==(const Color &o) const
    {
        return b == o.b && g == o.g && r == o.r;
    }
};

vector<pair<int, int>> getLine(int x0, int y0, int x1, int y1)
{
    // printf("getLine: (%d, %d) - (%d, %d)\n", x0, y0, x1, y1);

    bool swapped = 0;
    if (abs(x1 - x0) < abs(y1 - y0))
    {
        swap(x0, y0), swap(x1, y1);
        swapped = 1;
    }
    if (x0 > x1)
        swap(x0, x1), swap(y0, y1);

    // printf("-> (%d, %d) - (%d, %d)\n", x0, y0, x1, y1);

    vector<pair<int, int>> ans;

    int x = x0,
        y = y0, dx = x1 - x0, dy = y1 - y0, e = -dx;
    if (y1 > y0)
        for (int i = 0; i <= dx; ++i)
        {
            ans.push_back(make_pair(x, y));
            ++x, e += 2 * dy;
            if (e >= 0)
                ++y, e -= 2 * dx;
        }
    else
        for (int i = 0; i <= dx; ++i)
        {
            ans.push_back(make_pair(x, y));
            ++x, e += 2 * dy;
            if (e <= 0)
                --y, e += 2 * dx;
        }

    if (swapped)
        for (pair<int, int> &pp : ans)
            swap(pp.first, pp.second);

    return ans;
}

void draw(Mat &canvas, vector<pair<int, int>> points, Color c)
{
    for (pair<int, int> pp : points)
    {
        int x = pp.first, y = pp.second;
        canvas.at<Vec3b>(x, y)[0] = c.b, canvas.at<Vec3b>(x, y)[1] = c.g, canvas.at<Vec3b>(x, y)[2] = c.r;
    }

    // printf("First point: (%d, %d)\n", points[0].first, points[0].second);
    // printf("Last point: (%d, %d)\n", points.back().first, points.back().second);
}

void areafill(Mat &canvas, int x, int y, Color oldColor, Color newColor)
{
    stack<pair<int, int>> s;
    s.push(make_pair(x, y));

    while (!s.empty())
    {
        x = s.top().first, y = s.top().second;
        s.pop();

        if (0 <= x && x < height && 0 <= y && y < width && canvas.at<Vec3b>(x, y)[0] == oldColor.b && canvas.at<Vec3b>(x, y)[1] == oldColor.g && canvas.at<Vec3b>(x, y)[2] == oldColor.r)
        {
            canvas.at<Vec3b>(x, y)[0] = newColor.b, canvas.at<Vec3b>(x, y)[1] = newColor.g, canvas.at<Vec3b>(x, y)[2] = newColor.r;

            // static int dx[] = {-1, -1, -1, 0, 0, 0, 1, 1, 1}, dy[] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
            // for (int i = 0; i < 8; ++i)

            static int dx[] = {-1, 0, 1, 0}, dy[] = {0, -1, 0, 1};
            for (int i = 0; i < 4; ++i)
                s.push(make_pair(x + dx[i], y + dy[i]));
        }
    }
}

int main()
{
    Mat canvas(height, width, CV_8UC3, cv::Scalar(0, 0, 0));

    // Draw a n-polygon
    {
        double phi = 0, dphi = 2 * acos(-1) / n;
        int d = 100;
        for (int i = 0; i < n; ++i)
        {
            draw(canvas, getLine(150 + d * cos(phi), 150 + d * sin(phi), 150 + d * cos(phi + dphi), 150 + d * sin(phi + dphi)), Color({255, 255, 255}));
            phi += dphi;
        }
    }

    imwrite("p1.png", canvas);

    areafill(canvas, 150, 150, Color({0, 0, 0}), Color({255, 255, 255}));

    imwrite("p2.png", canvas);
}