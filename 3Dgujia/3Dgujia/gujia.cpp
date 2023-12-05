//#include <opencv2/opencv.hpp>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <chrono>
#include<map>
#include<algorithm>
#pragma warning(disable:4996)
using namespace std;
//using namespace cv;
using namespace chrono;
struct PM {//ƽ��
    double a, b, c;
    double x, y, z;
    PM() {}
    PM(double a, double b, double c, double x, double y, double z) :a(a), b(b), c(c), x(x), y(y), z(z) {}
}up, down;
struct point {
    double x, y, z;
    point() {}
    point(double x, double y, double z) :x(x), y(y), z(z) {}
    point operator+(const point& a)const& {
        return point(x + a.x, y + a.y, z + a.z);
    }
    point operator-(const point& a)const& {
        return point(x - a.x, y - a.y, z - a.z);
    }
    friend point operator*(const double& a, const point& b) {
        return point(a * b.x, a * b.y, a * b.z);
    }
};
struct triangle {
    point a, b, c;
    triangle() {}
    triangle(point a, point b, point c) :a(a), b(b), c(c) {}

};
vector<triangle>triangles;
map<double, int>mp;
int flag;
double xx, yy, zz;
float ma = -1e9, mi = 1e9;

void read_stl(string path, vector<double>& triangles)
{
    ifstream stl(path, fstream::binary);
    stl.seekg(80);
    char _head[4];
    stl.read(_head, 4);
    stl.seekg(84);
    uint32_t triangle_size;
    memmove(&triangle_size, _head, 4);
    //cout << triangle_size << '\n';
    vector<char> points(triangle_size * 50);
    stl.read((char*)&points[0], triangle_size * 50);
    triangles.resize(9 * triangle_size);
    float point3f_3[12];
    freopen("stl.txt", "w", stdout);
    for (int i = 0; i < triangle_size; i++)
    {
        memmove(point3f_3, &(points[i * 50 + 12]), 36);
        for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < 3; k++) {
                triangles[i * 9 + j * 3 + k] = point3f_3[j * 3 + k];
            }
            ma = max(ma, point3f_3[j * 3 + 2]);
            mi = min(mi, point3f_3[j * 3 + 2]);
        }
        for (size_t j = 0; j < 3; ++j) {
            cout << triangles[i * 9 + j * 3] << ' ' << triangles[i * 9 + j * 3 + 1] << ' ' << triangles[i * 9 + j * 3 + 2] << '\n';
        }
    }
    //freopen("CON", "w", stdout);
}

double findjiaodian(PM a, double x, double y, double z, double vx, double vy, double vz) {//��ӷ�������ʽ��ʾֱ��
    double vp1, vp2, vp3, n1, n2, n3, v1, v2, v3, m1, m2, m3, t, vpt;
    vp1 = a.a;vp2 = a.b;vp3 = a.c;n1 = a.x;n2 = a.y;n3 = a.z;v1 = vx;v2 = vy;v3 = vz;m1 = x;m2 = y;m3 = z;
    //��ֱ֪��L����m��m1��m2��m3�����ҷ�������ΪVL��v1��v2��v3����ƽ��P����n��n1��n2��n3�����ҷ��߷�������ΪVP��vp1��vp2��vp3��
    vpt = v1 * vp1 + v2 * vp2 + v3 * vp3;
    //�����ж�ֱ���Ƿ���ƽ��ƽ��
    if (vpt == 0)
    {
        return -1;//ֻ����0-1֮������߶��ϣ�����-1��ʾû�н���
    }
    else
    {
        t = ((n1 - m1) * vp1 + (n2 - m2) * vp2 + (n3 - m3) * vp3) / vpt;
        return t;
    }
}

point recursive_bezier(const std::vector<point>& control_points, float t)
{
    // TODO: Implement de Casteljau's algorithm

    if (control_points.size() == 1) return control_points[0];

    std::vector<point> a;
    for (int i = 0; i + 1 < control_points.size(); i++) {
        auto p = control_points[i] + t * (control_points[i + 1] - control_points[i]);
        a.push_back(p);
    }

    return recursive_bezier(a, t);
}

void bezier(const std::vector<point>& control_points)
{
    // TODO: Iterate through all t = 0 to t = 1 with small steps, and call de Casteljau's 
    // recursive Bezier algorithm.
    freopen("out.txt", "w", stdout);
    double delta = 0.001;
    for (double t = 0; t <= 1; t += delta) {
        auto point = recursive_bezier(control_points, t);
        cout << point.x << ' ' << point.y << ' ' << point.z << '\n';
        int w = 1;
        for (int i = -w + 1; i <= w; i++) {
            for (int j = -w + 1; j <= w; j++) {
                for (int k = -w + 1; k <= w; k++) {
                    double x = point.x + i, y = point.y + j, z = point.z + k;
                    //cout << x << ' ' << y << ' ' << z << '\n';
                    /*double dist = sqrt(pow(point.x - x, 2) + pow(point.y - y, 2) + pow(point.z - z, 2));
                    window.at<cv::Vec3b>(z, y, x)[1] = std::min(window.at<cv::Vec3b>(z, y, x)[1] + 255 * std::max(2 - exp(dist), 0.0), 255.0);*/
                }
                // auto k = abs(((int)(point.x+1-i))-point.x) * abs(((int)(point.y+1-j))-point.y);
                // window.at<cv::Vec3b>(y, x)[1] = std::min(window.at<cv::Vec3b>(y, x)[1] + 255 * k, 255.0f);
            }
        }
    }
}

void judgejiao(PM p, point a, point b) {
    double vx, vy, vz;
    vx = b.x - a.x;
    vy = b.y - a.y;
    vz = b.z - a.z;
    double f, f1;
    f = a.x;
    f = f * 1000 + a.y;
    f = f * 1000 + a.z;
    f = f * 1000 + vx + a.x;
    f = f * 1000 + vy + a.y;
    f = f * 1000 + vz + a.z;
    f1 = vx + a.x;
    f1 = f1 * 1000 + vy + a.y;
    f1 = f1 * 1000 + vz + a.z;
    f1 = f1 * 1000 + a.x;
    f1 = f1 * 1000 + a.y;
    f1 = f1 * 1000 + a.z;
    double t;
    if (mp.count(f) == 0 && mp.count(f1) == 0) {
        t = findjiaodian(p, a.x, a.y, a.z, vx, vy, vz);
        mp[f] = 1;
        mp[f1] = 1;
    }
    else {
        t = -1;
    }
    if (t >= 0 && t <= 1) {
        flag++;
        xx += a.x + t * vx;
        yy += a.y + t * vy;
        zz += a.z + t * vz;
    }
}

double dotproduct(point a, point b) {//���
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double moudle(point a) {//������ģ
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

bool pointcmp(PM p, point a, point b) {//�жϵ�a�͵�b��ƽ��p�������µĴ�С��ϵ����ƽ��p�����������µĴ�С��ϵ��
    double d;
    d = dotproduct(point(a.x - b.x, a.y - b.y, a.z - b.z), point(p.a, p.b, p.c));
    d = d / moudle(point(p.a, p.b, p.c));
    return d < 0;
}

bool cmp1(triangle a, triangle b) {//�ж����������ζ���downƽ����Ե�λ���Ⱥ��ϵ
    point ma, mb;
    ma = a.a;
    mb = b.a;
    if (pointcmp(down, a.b, ma))ma = a.b;
    if (pointcmp(down, a.c, ma))ma = a.c;
    if (pointcmp(down, b.b, mb))mb = b.b;
    if (pointcmp(down, b.c, mb))mb = b.c;
    return pointcmp(down, ma, mb);
}

bool cmp2(triangle a, triangle b) {//�ж����������ζ���upƽ����Ե�λ���Ⱥ��ϵ
    point ma, mb;
    ma = a.a;
    mb = b.a;
    if (pointcmp(up, a.b, ma))ma = a.b;
    if (pointcmp(up, a.c, ma))ma = a.c;
    if (pointcmp(up, b.b, mb))mb = b.b;
    if (pointcmp(up, b.c, mb))mb = b.c;
    return pointcmp(up, ma, mb);
}

bool cmp3(triangle a, PM p) {//�ж���������ƽ����ı�
    point m = a.a;
    if (pointcmp(p, a.b, m))m = a.b;
    if (pointcmp(p, a.c, m))m = a.c;
    //return pointcmp(p, m, point(p.x, p.y, p.z));
    point pp = point(p.x, p.y, p.z);
    double d;
    d = dotproduct(point(m.x - pp.x, m.y - pp.y, m.z - pp.z), point(p.a, p.b, p.c));
    d = d / moudle(point(p.a, p.b, p.c));
    return d <= 0;
}

int main(int argc, char const* argv[])
{
    
    //////////////////////////////���뿪ʼ////////////////////////////
    cout << fixed << setprecision(6);
    vector<double> triangles_all;
    read_stl("new.stl", triangles_all);
    for (int i = 0;i < triangles_all.size();i += 9) {
        point a, b, c;
        a = point(triangles_all[i], triangles_all[i + 1], triangles_all[i + 2]);
        b = point(triangles_all[i + 3], triangles_all[i + 4], triangles_all[i + 5]);
        c = point(triangles_all[i + 6], triangles_all[i + 7], triangles_all[i + 8]);
        triangle p = triangle(a, b, c);
        triangles.push_back(p);
    }
    //up��z���ĵ�,down��z��С�ĵ�aaa
    up.z = mi; down.z = ma;
    for (int i = 0; i < triangles.size(); ++i) {
        if (up.z < triangles[i].a.z) {
            up = PM(0, 0, 0, triangles[i].a.x, triangles[i].a.y, triangles[i].a.z);
        }
        if (up.z < triangles[i].b.z) {
            up = PM(0, 0, 0, triangles[i].b.x, triangles[i].b.y, triangles[i].b.z);
        }
        if (up.z < triangles[i].c.z) {
            up = PM(0, 0, 0, triangles[i].c.x, triangles[i].c.y, triangles[i].c.z);
        }
        if (down.z > triangles[i].a.z) {
            down = PM(0, 0, 1, triangles[i].a.x, triangles[i].a.y, triangles[i].a.z);
        }
        if (down.z > triangles[i].b.z) {
            down = PM(0, 0, 1, triangles[i].b.x, triangles[i].b.y, triangles[i].b.z);
        }
        if (down.z > triangles[i].c.z) {
            down = PM(0, 0, 1, triangles[i].c.x, triangles[i].c.y, triangles[i].c.z);
        }
    }//��up��down�㷨ʽ�еĵ�

    //����up�ķ�����
    double X = 0, Y = 466 - 458, Z = 348 - 361;
    double XX = 1, YY = 0, ZZ = 0;
    up.a = Y * ZZ - Z * YY;
    up.b = Z * XX - X * ZZ;
    up.c = X * YY - Y * XX;
    if (up.c > 0) {
        up.a = -up.a;
        up.b = -up.b;
        up.c = -up.c;
    }//up��Ӧ�����ƶ���������z��������0��ת

    vector<point>points, tmppoints;//��¼���ĵ�

    double wucha = 1e9;//����wucha����Ϊ��벿��û�б�Ҫ

    double len = 1;//ƽ��ʱ�ĳ���

    long long howmany = 90;//����ȡ���ٲ�

    //down������ƽ��
    //freopen("CON", "w", stdout);
    sort(triangles.begin(), triangles.end(), cmp1);
    double fp = sqrt(down.a * down.a + down.b * down.b + down.c * down.c);
    int ccnt = 0;
    double prex = -1000000007, prey = -1000000007, prez = -1000000007;
    flag = -1;
    xx = 0, yy = 0, zz = 0;
    while (1) {
        //cout << flag << endl;
        ccnt++;
        if (points.size() > howmany)break;
        //cout << flag << '\n';
        if (flag == 0)break;
        xx /= flag;yy /= flag;zz /= flag;
        if (prex == -1000000007 && prey == -1000000007 && prez == -1000000007) {
            prex = xx;
            prey = yy;
            prez = zz;
        }
        else {
            /*cout << xx << ' ' << yy << ' ' << zz << '\n';
            cout << prex << ' ' << prey << ' ' << prez << '\n';
            cout << (xx - prex) * (xx - prex) + (yy - prey) * (yy - prey) + (zz - prez) * (zz - prez) << '\n';*/
            if ((xx - prex) * (xx - prex) + (yy - prey) * (yy - prey) + (zz - prez) * (zz - prez) > wucha) {
                break;
            }
            else {
                prex = xx;
                prey = yy;
                prez = zz;
            }
        }
        mp.clear();
        if (flag != -1) {
            if (ccnt > 15)//��Ϊ����ƽ�治��ȫƽ�У�����ɾ��һЩ�����ٵ������
                points.push_back(point(xx, yy, zz));
        }
        flag = 0;
        xx = 0, yy = 0, zz = 0;
        for (int i = 0;i < triangles.size(); ++i) {
            if (!cmp3(triangles[i], down)) break;
            judgejiao(down, triangles[i].a, triangles[i].b);
            judgejiao(down, triangles[i].a, triangles[i].c);
            judgejiao(down, triangles[i].b, triangles[i].c);
        }
        down.x += down.a / fp * len;
        down.y += down.b / fp * len;
        down.z += down.c / fp * len;//�����������ƶ�len
    }

    //up����ƽ��
    //freopen("CON", "w", stdout);
    sort(triangles.begin(), triangles.end(), cmp2);
    fp = sqrt(up.a * up.a + up.b * up.b + up.c * up.c);
    ccnt = 0;
    prex = -1000000007, prey = -1000000007, prez = -1000000007;
    flag = -1;
    xx = 0, yy = 0, zz = 0;
    while (1) {
        //cout << flag << endl;
        ccnt++;
        if (tmppoints.size() > howmany)break;
        if (flag == 0)break;
        xx /= flag;yy /= flag;zz /= flag;
        if (prex == -1000000007 && prey == -1000000007 && prez == -1000000007) {
            prex = xx;
            prey = yy;
            prez = zz;
        }
        else {
            if ((xx - prex) * (xx - prex) + (yy - prey) * (yy - prey) + (zz - prez) * (zz - prez) > wucha) {
                break;
            }
        }
        mp.clear();
        if (flag != -1) {
            if (ccnt > 10)//��Ϊ����ƽ�治��ȫƽ�У�����ɾ��һЩ�����ٵ������
                tmppoints.push_back(point(xx, yy, zz));
        }
        flag = 0;
        xx = 0, yy = 0, zz = 0;
        for (int i = 0;i < triangles.size(); ++i) {
            if (!cmp3(triangles[i], up)) break;
            judgejiao(up, triangles[i].a, triangles[i].b);
            judgejiao(up, triangles[i].a, triangles[i].c);
            judgejiao(up, triangles[i].b, triangles[i].c);
        }
        up.x += up.a / fp * len;
        up.y += up.b / fp * len;
        up.z += up.c / fp * len;//�����������ƶ�len
    }

    reverse(tmppoints.begin(), tmppoints.end());
    std::vector<point>bp;
    double px=0;
    for (int i = 0; i < points.size(); i++) {
        px += points[i].x;
    }
    for (int i = 0; i < tmppoints.size(); i++) {
        px += tmppoints[i].x;
    }
    px /= points.size() + tmppoints.size();
    freopen("cut.txt", "w", stdout);
    for (int i = 0;i < points.size();i++) {
        cout << px << ' ' << points[i].y << ' ' << points[i].z << '\n';
        bp.emplace_back(point(px, points[i].y, points[i].z));
    }
    for (int i = 0;i < tmppoints.size();i++) {
        cout << px << ' ' << tmppoints[i].y << ' ' << tmppoints[i].z << '\n';
        bp.emplace_back(point(px, tmppoints[i].y, tmppoints[i].z));
        //points.emplace_back(i);
    }

    bezier(bp);

    return 0;
}