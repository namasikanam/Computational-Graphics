#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <set>
#include <cstdlib>
#include <cstring>
#include <cassert>
using namespace std;
typedef long long ll;

const double t = 0.1;

const double eps = 1e-5;

struct Point
{
    double x, y, z;
    double Q[4][4];
    Point() : x(0), y(0), z(0) { memset(Q, 0, sizeof(Q)); }
    Point(char *s)
    {
        sscanf(s, "%lf%lf%lf", &x, &y, &z);
        memset(Q, 0, sizeof(Q));
    }
    void addQ(double _Q[4][4])
    {
        for (int i = 4; i--;)
            for (int j = 4; j--;)
                Q[i][j] += _Q[i][j];
    }
    double calcError()
    {
        double ans = 0;
        double v[4] = {x, y, z, 1};
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                ans += Q[i][j] * v[i] * v[j];
        return ans;
    }
};
vector<Point> points;
vector<bool> point_deleted;
vector<int> point_final_id;
inline double distance(int u, int v) { return sqrt(pow(points[u].x - points[v].x, 2) + pow(points[u].y - points[v].y, 2) + pow(points[u].z - points[v].z, 2)); }

struct Vector
{
    double x, y, z;
    Vector() : x(0), y(0), z(0) {}
    Vector(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    Vector(int i) : x(points[i].x), y(points[i].y), z(points[i].z) {}
    Vector(Point p) : x(p.x), y(p.y), z(p.z) {}
    Vector operator+(Vector o) const { return Vector(x + o.x, y + o.y, z + o.z); }
    Vector operator-(Vector o) const { return Vector(x - o.x, y - o.y, z - o.z); }
    Vector operator%(Vector o) const { return Vector(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x); }
    double operator*(Vector o) const { return x * o.x + y * o.y + z * o.z; }
    Vector operator*(double o) const { return Vector(x * o, y * o, z * o); }
    void copy_to(Point &p) { p.x = x, p.y = y, p.z = z; }
};

struct Face
{
    int point[3];
    Face(char *s)
    {
        sscanf(s, "%d%d%d", point, point + 1, point + 2);
        for (int i = 3; i--;)
            --point[i];
    }
};
vector<Face> faces;
vector<bool> face_deleted;
vector<Vector> normals;

set<pair<int, int>> existing_edge;
vector<bool> is_boundary;
vector<bool> in_union;

vector<vector<int>> neighbored_faces;

vector<pair<int, int>> valid_pairs; // (u, v)有保证u < v. 且是并查集的根节点。
vector<bool> pair_deleted;
vector<vector<int>> neighbored_pairs;

vector<int> heap_id;
struct HeapNode
{
    int pair_id;
    double cost;
    bool inversion;
    Point newnode;
    HeapNode() {}
    HeapNode(int _pair_id) : pair_id(_pair_id), newnode(Point())
    {
        // printf("New Heap Node with %d(pair_id)\n", pair_id);

        int u = valid_pairs[pair_id].first, v = valid_pairs[pair_id].second;
        newnode.addQ(points[u].Q), newnode.addQ(points[v].Q);

        // 试图高斯消元
        double mat[4][5];
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 4; ++j)
                mat[i][j] = newnode.Q[i][j];
            mat[i][4] = 0;
        }
        memset(mat[3], 0, sizeof(double) * 3);
        fill(mat[3] + 3, mat[3] + 5, 1.0);
        if ([&mat]() -> bool {
                for (int i = 0; i < 3; ++i)
                {
                    if (!([&mat, i]() -> bool { // 列主元是否非零
                            int _j = i;
                            for (int j = i + 1; j < 3; ++j)
                                if (fabs(mat[j][i]) > fabs(mat[_j][i]))
                                    _j = j;
                            if (fabs(mat[_j][i] > eps))
                            {
                                for (int k = i; k < 5; ++k)
                                    swap(mat[i][k], mat[_j][k]);
                                return 1;
                            }
                            else
                                return 0;
                        }()))
                        return 0;
                    for (int j = i + 1; j < 3; ++j)
                    {
                        double tmp = mat[j][i] / mat[i][i];
                        for (int k = j; k < 5; ++k)
                            mat[j][k] -= mat[i][k] * tmp;
                    }
                }
                for (int i = 4; i--;)
                {
                    mat[i][4] /= mat[i][i];
                    for (int j = i; j--;)
                        mat[j][4] -= mat[j][i] * mat[i][4];
                }
                return 1; // 最终，答案存在mat[i][4]中。
            }())
        {
            newnode.x = mat[0][4], newnode.y = mat[1][4], newnode.z = mat[2][4];
        }
        else
        {
            Vector v1 = u, v2 = v;
            double a = 0, b = 0, c = 0;
            for (int i = 3; i--;)
                for (int j = 3; j--;)
                    a += mat[i][j];
            for (int i = 3; i--;)
                b += mat[i][3] * 2;
            c = mat[3][3];
            double x;
            if (a <= eps || [a, b, c, &x]() -> bool {
                    double delta = b * b - 4 * a * c;
                    if (delta < 0)
                        return 1;
                    double x1 = (-b + sqrt(delta)) / (2 * a), x2 = (-b - sqrt(delta)) / (2 * a);
                    if (x1 >= 0 && x1 <= 1)
                    {
                        x = x1;
                        return 0;
                    }
                    if (x2 >= 0 && x2 <= 1)
                    {
                        x = x2;
                        return 0;
                    }
                    return 1;
                }())
            {
                v1.copy_to(newnode);
                double cost1 = newnode.calcError();
                v2.copy_to(newnode);
                double cost2 = newnode.calcError();
                if (cost1 < cost2)
                    v1.copy_to(newnode);
            }
            else
            {
                (v1 + (v2 - v1) * x).copy_to(newnode);
            }
        }

        cost = newnode.calcError() - points[u].calcError() - points[v].calcError(); // 没有按paper上的来，算是改进了一下吧（不知道这个改进有没有用）。

        // 下面要来计算inversion：是否出现了面反转
        inversion = false;
        auto check_inverse = [this, u, v](int face_id) {
            Vector point[3];
            for (int i = 0; i < 3; ++i)
                if (faces[face_id].point[i] == u || faces[face_id].point[i] == v)
                    point[i] = Vector(newnode);
                else
                    point[i] = Vector(faces[face_id].point[i]);
            return (point[1] - point[0]) % (point[2] - point[0]) * normals[face_id] < -eps;
        };
        for (int face_id : neighbored_faces[u])
            if (!face_deleted[face_id] && check_inverse(face_id))
            {
                inversion = true;
                break;
            }
        if (!inversion)
            for (int face_id : neighbored_faces[v])
                if (!face_deleted[face_id] && check_inverse(face_id))
                {
                    inversion = false;
                    break;
                }

        // printf("Finishing building new heap node.\n");
    }
    bool operator<(const HeapNode &right_element) const
    {
        return inversion != right_element.inversion ? inversion < right_element.inversion : cost < right_element.cost;
    }
};
vector<HeapNode> heap;
inline void upHeap(int node)
{
    for (int next = node >> 1; next && heap[node] < heap[next]; node >>= 1, next >>= 1)
    {
        swap(heap_id[heap[node].pair_id], heap_id[heap[next].pair_id]);
        swap(heap[node], heap[next]);
    }
}
inline void downHeap(int node)
{
    for (int next = node << 1; next < heap.size(); node = next, next <<= 1)
    {
        if (next + 1 < heap.size() && heap[next + 1] < heap[next])
            ++next;
        if (heap[node] < heap[next])
            break;
        swap(heap_id[heap[node].pair_id], heap_id[heap[next].pair_id]);
        swap(heap[node], heap[next]);
    }
}
inline void initHeap()
{
    heap.push_back(HeapNode()); // 1-based heap
    for (int i = 0; i < valid_pairs.size(); ++i)
    {
        heap_id.push_back(i + 1);
        heap.push_back(HeapNode(i));
    }
    for (int i = heap.size(); --i;)
        upHeap(i);
}
inline void delHeap(int pair_id)
{
    pair_deleted[pair_id] = 1;
    int heap_pos = heap_id[pair_id];
    heap[heap_pos] = heap.back();
    heap_id[heap[heap_pos].pair_id] = heap_pos;
    heap.pop_back();
    upHeap(heap_pos), downHeap(heap_pos);
}
inline void updateHeap(int pair_id)
{
    // printf("Start updating heap (%d)\n", pair_id);

    int heap_pos = heap_id[pair_id];

    // printf("heap_pos = %d\n", heap_pos);

    heap[heap_pos] = HeapNode(pair_id);
    upHeap(heap_pos), downHeap(heap_pos);

    // printf("End updating heap %d(pair_id) -> %d(heap_id)\n", pair_id, heap_id[pair_id]);
}

int main(int argv, char **argc)
{
    bool flag = 0;
    freopen("meshsimplifier.out", "w", stdout);
    // 读
    FILE *fin, *fout;
    double simp_ratio;
    if (argv != 4)
    {
        fin = fopen("Arma.obj", "r");
        fout = fopen("Arma_out.obj", "w");
        simp_ratio = 0.1;
    }
    else
    {
        fin = fopen(argc[1], "r");
        fout = fopen(argc[2], "w");
        simp_ratio = atof(argc[3]);
    }

    // 把.opt读进来
    for (char s[100]; fgets(s, 100, fin) != nullptr;)
        if (s[0] == 'v')
            points.push_back(Point(s + 1));
        else if (s[0] == 'f')
            faces.push_back(Face(s + 1));

    int V = points.size();
    point_deleted.resize(V);
    is_boundary.resize(V);
    neighbored_pairs.resize(V);
    in_union.resize(V);

    face_deleted.resize(faces.size());

    for (int i = 0; i < V; ++i)
    {
        printf("points[%d] = (%.6f, %.6f, %.6f)\n", i, points[i].x, points[i].y, points[i].z);
    }
    puts("=============Finish Reading================");

    // 求出boundary node
    vector<int> point_counter(V);
    for (auto face : faces)
        for (int i = 3; i--;)
        {
            pair<int, int> pp = make_pair(face.point[i], face.point[(i + 1) % 3]);
            if (pp.first > pp.second)
                swap(pp.first, pp.second);

            if (existing_edge.find(pp) == existing_edge.end())
            {
                ++point_counter[pp.first], ++point_counter[pp.second];
                existing_edge.insert(pp);
            }
            else
                --point_counter[pp.first], --point_counter[pp.second];
        }
    for (int i = V; i--;)
        is_boundary[i] = point_counter[i];

    puts("========Finding all boundary nodes===============");

    // 求neighbored_faces
    neighbored_faces.resize(V);
    for (int i = faces.size(); i--;)
        for (int j = 3; j--;)
            if (!is_boundary[faces[i].point[j]])
                neighbored_faces[faces[i].point[j]].push_back(i);
    for (int i = points.size(); i--;)
        sort(neighbored_faces[i].begin(), neighbored_faces[i].end());

    // 求initial error quadrics和normals
    for (auto face : faces)
    {
        Vector normal = (Vector(face.point[1]) - Vector(face.point[0])) % (Vector(face.point[2]) - Vector(face.point[0]));
        normals.push_back(normal);

        double P[4] = {normal.x, normal.y, normal.z, normal * Vector(face.point[0]) * -1};
        double Q[4][4] = {0};
        for (int i = 4; i--;)
            for (int j = 4; j--;)
                Q[i][j] = P[i] * P[j];

        for (int i = 3; i--;)
            if (!is_boundary[face.point[i]])
                points[face.point[i]].addQ(Q);
    }

    puts("=======Initialzing error quadrics and normals========");

    // 求出valid pairs
    vector<int> sorted_point_y, unsorted_point_y, sorted_point_x;
    for (int i = V; i--;)
        if (!is_boundary[i])
            sorted_point_x.push_back(i);
    auto cmp = [](int u, int v) {
        return points[u].x != points[v].x ? points[u].x < points[v].x : u < v;
    };
    sort(sorted_point_x.begin(), sorted_point_x.end(), cmp);
    auto addPair = [](int u, int v) {
        neighbored_pairs[u].push_back(valid_pairs.size()), neighbored_pairs[v].push_back(valid_pairs.size());
        valid_pairs.push_back(u < v ? make_pair(u, v) : make_pair(v, u));
    };
    for (int _V = sorted_point_x.size(), B = sqrt(_V), i = 0, j = 0; i < _V; ++i)
    {
        // printf("====Calculate: sorted_point_x[%d]=%d====\n", i, sorted_point_x[i]);

        // 先把需要的点弄进unsorted_point_y里
        while (j < _V && points[sorted_point_x[j]].x <= points[sorted_point_x[i]].x + t)
            unsorted_point_y.push_back(sorted_point_x[j++]);
        // 如果unsorted_point_y够多了，就合并进sorted_point_y，记得把sorted_point_x[i]左边的点删掉。
        if (unsorted_point_y.size() > B)
        {
            vector<int> tmp_vector;
            sort(unsorted_point_y.begin(), unsorted_point_y.end(), cmp);
            for (auto it = sorted_point_y.begin(), jt = unsorted_point_y.begin(); it != sorted_point_y.end() || jt != unsorted_point_y.end();)
            {
                int tmp;
                if (it == sorted_point_y.end())
                    tmp = *jt++;
                else if (jt == unsorted_point_y.end())
                    tmp = *it++;
                else if (cmp(*it, *jt))
                    tmp = *it++;
                else
                    tmp = *jt++;
                if (cmp(sorted_point_x[i], tmp))
                    tmp_vector.push_back(tmp);
            }
            sorted_point_y = tmp_vector;
            unsorted_point_y.clear();
        }

        // puts("===After merge===");

        // 找与sorted_point_x[i]对应的有效对，记得别考虑左边的点。
        int p = sorted_point_x[i];
        for (auto it = lower_bound(sorted_point_y.begin(), sorted_point_y.end(), p, [](int u, int v) { return points[u].y < points[v].y; }); it != sorted_point_y.end() && points[*it].y < points[p].y + t; ++it)
            if (cmp(p, *it) && distance(*it, p) < t)
                addPair(*it, p);
        for (auto it = unsorted_point_y.begin(); it != unsorted_point_y.end(); ++it)
            if (cmp(p, *it) && distance(*it, p) < t)
                addPair(*it, p);
    }

    existing_edge.clear();
    for (auto face : faces)
        for (int i = 3; i--;)
        {
            pair<int, int> pp = make_pair(face.point[i], face.point[(i + 1) % 3]);
            if (!is_boundary[pp.first] && !is_boundary[pp.second])
            {
                if (pp.first > pp.second)
                    swap(pp.first, pp.second);

                if (existing_edge.find(pp) == existing_edge.end())
                    existing_edge.insert(pp);
                else
                    addPair(pp.first, pp.second);
            }
        }

    for (int i = 0; i < valid_pairs.size(); ++i)
        printf("valid_pairs[%d] = (%d, %d)\n", i, valid_pairs[i].first, valid_pairs[i].second);
    puts("===All valid pairs are got===");

    pair_deleted.resize(valid_pairs.size());

    puts("=======Finding all valid pairs=========");

    // 数据结构初始化
    initHeap();

    puts("==========Finish Heap Initialization===============");

    // 贪心
    auto getAno = [](int pair_id, int u) { return valid_pairs[pair_id].first ^ valid_pairs[pair_id].second ^ u; };
    for (int _V = V * simp_ratio; _V < V; ++_V)
    {
        // 从堆中取出cost最小的点对
        if (heap.size() == 1 || heap[1].inversion)
            break;

        // 更新
        int u = valid_pairs[heap[1].pair_id].first, v = valid_pairs[heap[1].pair_id].second;
        if (u > v)
            swap(u, v);
        // assert(!pair_deleted[heap[1].pair_id]);

        printf("------Aggregate %d %d--------\n", u, v);

        points[u] = heap[1].newnode;

        point_deleted[v] = 1;
        { // 更新面片
            vector<int> new_neighbored_faces;
            for (int i = 0, j = 0; i < neighbored_faces[u].size() || j < neighbored_faces[v].size();)
            {
                int face_id;
                if (i == neighbored_faces[u].size())
                    face_id = neighbored_faces[v][j++];
                else if (j == neighbored_faces[v].size())
                    face_id = neighbored_faces[u][i++];
                else if (neighbored_faces[u][i] < neighbored_faces[v][j])
                    face_id = neighbored_faces[u][i++];
                else if (neighbored_faces[u][i] > neighbored_faces[v][j])
                    face_id = neighbored_faces[v][j++];
                else
                {
                    face_id = neighbored_faces[u][i];
                    ++i, ++j;
                }

                int cnt = 0;
                for (int i = 3; i--;)
                {
                    if (faces[face_id].point[i] == v)
                        faces[face_id].point[i] = u;
                    cnt += faces[face_id].point[i] == u;
                }
                if (!face_deleted[face_id])
                    if (cnt == 1)
                        new_neighbored_faces.push_back(face_id);
                    else
                        face_deleted[face_id] = 1;
            }
            neighbored_faces[u] = new_neighbored_faces;
        }

        puts("----Neighbored Faces are updated----");

        { // 更新有效对及其堆
            delHeap(heap[1].pair_id);
            vector<int> new_neighbored_pairs;
            for (int pair_id : neighbored_pairs[v])
            {
                if (valid_pairs[pair_id].first == v)
                    valid_pairs[pair_id].first = u;
                if (valid_pairs[pair_id].second == v)
                    valid_pairs[pair_id].second = u;
                if (!pair_deleted[pair_id])
                {
                    new_neighbored_pairs.push_back(pair_id);
                    in_union[getAno(pair_id, u)] = 1;
                    updateHeap(pair_id);
                }
            }
            for (int pair_id : neighbored_pairs[u])
                if (!pair_deleted[pair_id])
                    if (in_union[getAno(pair_id, u)])
                        delHeap(pair_id);
                    else
                    {
                        new_neighbored_pairs.push_back(pair_id);
                        updateHeap(pair_id);
                    }
            for (int pair_id : neighbored_pairs[v])
                if (!pair_deleted[pair_id])
                    in_union[getAno(pair_id, u)] = 0;
            neighbored_pairs[u] = new_neighbored_pairs;
        }

        puts("----Neighbored Pairs are updated-----");
        // printf("valid_pairs[5679] = (%d, %d)\n", valid_pairs[5679].first, valid_pairs[5679].second);
        // printf("valid_pairs[5682] = (%d, %d)\n", valid_pairs[5682].first, valid_pairs[5682].second);

        // if (flag)
        //     exit(0);
    }

    puts("=========Finish Mesh Simplification===========");

    // 写答案
    {
        int cur_id = 0;
        point_final_id.resize(V);

        double minNorm = 1e9;

        for (int i = 0; i < points.size(); ++i)
            if (!point_deleted[i])
            {
                point_final_id[i] = ++cur_id;
                fprintf(fout, "v %.6f %.6f %.6f\n", points[i].x, points[i].y, points[i].z);

                minNorm = min(minNorm, points[i].x * points[i].x + points[i].y * points[i].y + points[i].z * points[i].z);
            }

        printf("minNorm = %.6f\n", minNorm);

        for (int i = 0; i < faces.size(); ++i)
            if (!face_deleted[i])
            {
                fprintf(fout, "f");
                for (int j = 0; j < 3; ++j)
                {
                    fprintf(fout, " %d", point_final_id[faces[i].point[j]]);

                    if (point_final_id[faces[i].point[j]] == 0)
                    {
                        printf("face_id = %d, point = %d\n", i, faces[i].point[j]);
                        exit(0);
                    }
                }
                fprintf(fout, "\n");

                // printf("f");
                // for (int j = 0; j < 3; ++j)
                //     printf(" %d", point_final_id[faces[i].point[j]]);
                // printf("\n");
            }
    }
    puts("All is output.");
}