#include <fstream>
#include <iostream>
#include <stdio.h>
#include <chrono>
#include <vector>
#include <map>
#include <cmath>
#include <iomanip>
#pragma warning(disable:4996)
using namespace std;
using namespace chrono;
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
    point operator*(const double& a)const& {
        return point(x * a, y * a, z * a);
    }
    bool operator==(const point& a)const& {
        return fabs(x - a.x) <= 1e-9 && fabs(y - a.y) <= 1e-9 && fabs(z - a.z) <= 1e-9;
    }
};
struct PM {//平面
    point o, d;
    PM() {}
    PM(point o, point d) :o(o), d(d) {}
    PM(double a, double b, double c, double x, double y, double z) { d = point(a, b, c); o = point(x, y, z); }
}up, down;

struct line {
    point o, d;//o直线上一点，d方向向量
    line() {}
    line(point o, point d) :o(o), d(d) {}
};
struct triangle {
    point a, b, c;
    triangle() {}
    triangle(point a, point b, point c) :a(a), b(b), c(c) {}

};
vector<triangle>triangles;//三角形
vector<point>points;//骨架点
vector<point>np;//选出的法面切片的骨架点
vector<PM>fa;//法面
vector<point>qiedian[1005];
//map<double, int>mp;
int ad[1000005][3];//记录每个三角形的相邻三角形编号
int sad[1000005];//每个三角形相邻三角形的个数
map<int, point>mpp;
vector<int>hao[1005];
map<int, int>vis;
map<int, int>ss;
void show(triangle a) {
    cout << a.a.x << " " << a.a.y << " " << a.a.z << endl << a.b.x << " " << a.b.y << " " << a.b.z << endl << a.c.x << " " << a.c.y << " " << a.c.z << endl;
}

double dotproduct(point a, point b) {//点乘
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double findjiaodian(PM a, double x, double y, double z, double vx, double vy, double vz) {//点加法向量方式表示直线
    double vp1, vp2, vp3, n1, n2, n3, v1, v2, v3, m1, m2, m3, t, vpt;
    vp1 = a.d.x; vp2 = a.d.y; vp3 = a.d.z; n1 = a.o.x; n2 = a.o.y; n3 = a.o.z; v1 = vx; v2 = vy; v3 = vz; m1 = x; m2 = y; m3 = z;
    //已知直线L过点m（m1，m2，m3），且方向向量为VL（v1，v2，v3），平面P过点n（n1，n2，n3），且法线方向向量为VP（vp1，vp2，vp3）
    vpt = v1 * vp1 + v2 * vp2 + v3 * vp3;
    //首先判断直线是否与平面平行
    if (vpt == 0)
    {
        return -1;//只有在0-1之间才是线段上，返回-1表示没有交点
    }
    else
    {
        t = ((n1 - m1) * vp1 + (n2 - m2) * vp2 + (n3 - m3) * vp3) / vpt;
        return t;
    }
}

point judgejiao(int i, int j, PM p, point a, point b) {
    point oa = a - p.o, ob = b - p.o;
    int flag = 0;
    if (fabs(dotproduct(oa, p.d)) <= 1e-9 && fabs(dotproduct(ob, p.d)) <= 1e-9)flag = 1;
    if (flag)return point(1000000007, 1000000007, 1000000007);
    double vx, vy, vz;
    vx = b.x - a.x;
    vy = b.y - a.y;
    vz = b.z - a.z;
    /*double f, f1;
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
    f1 = f1 * 1000 + a.z;*/
    double t;
    //if (mp.count(f) == 0 && mp.count(f1) == 0) {
    t = findjiaodian(p, a.x, a.y, a.z, vx, vy, vz);
        //mp[f] = 1;
        //mp[f1] = 1;
    //}
    //else {
        //t = -1;
    //}
    if (t >= 0 && t <= 1) {
        //qiedian[i].push_back(point(a.x + t * vx, a.y + t * vy, a.z + t * vz));
        return point(a.x + t * vx, a.y + t * vy, a.z + t * vz);
        //hao[i].push_back(j);
    }
    else {
        return point(-1000000007, -1000000007, -1000000007);
    }
}

double moudle(point a) {//向量的模
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

bool pointcmp(PM p, point a, point b) {//判断点a和点b在平面p的意义下的大小关系（即平面p法向量意义下的大小关系）
    double d;
    d = dotproduct(point(a.x - b.x, a.y - b.y, a.z - b.z), p.d);
    d = d / moudle(p.d);
    return d < 0;
}

int jiaoshe(int x, int y) {
    PM p = fa[y];
    int num = 0;
    for (int i = 0; i < qiedian[x].size(); i++) {
        if (pointcmp(p, qiedian[x][i], p.o))num++;
    }
    return num;
}

point chacheng(point a, point b) {
    return point(a.y * b.z - b.y * a.z, -(a.x * b.z - b.x * a.z), a.x * b.y - b.x * a.y);
}

double dis(point a, point b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

void printline(point o, point d, double begin, double length, double step) {
    if (length < 0) {
        begin += length;
        length = -length;
    }
    for (double len = begin; len <= begin + length; len += step) {
        cout << o.x + d.x / moudle(d) * len << ' ' << o.y + d.y / moudle(d) * len << ' ' << o.z + d.z / moudle(d) * len << '\n';
    }
}

void arccut(int a, int b) {
    point X = point(1, 1, 1);
    cerr << a << " " << b << endl;
    point o, d;
    d = chacheng(fa[a].d, fa[b].d);
    double aa = fa[a].d.x, ab = fa[a].d.y, ac = fa[a].d.z;
    double ad = -(aa * fa[a].o.x + ab * fa[a].o.y + ac * fa[a].o.z);
    double ba = fa[b].d.x, bb = fa[b].d.y, bc = fa[b].d.z;
    double bd = -(ba * fa[b].o.x + bb * fa[b].o.y + bc * fa[b].o.z);
    if (fabs(d.x) >= 1e-9) {
        o.x = X.x;
        ad += aa; bd += ba;
        double dc = bc * ab - ac * bb;
        double dd = bc * ad - ac * bd;
        o.y = -dd / dc;
        if (fabs(ac) >= 1e-9)o.z = -(ab * o.y + ad) / ac;
        else if (fabs(bc) >= 1e-9)o.z = -(bb * o.y + bd) / bc;
        else o.z = X.z;
    }
    else if (fabs(d.y) >= 1e-9) {
        o.y = X.y;
        ad += ab; bd += bb;
        double da = bc * aa - ac * ba;
        double dd = bc * ad - ac * bd;
        o.x = -dd / da;
        if (fabs(ac) >= 1e-9)o.z = -(aa * o.x + ad) / ac;
        else if (fabs(bc) >= 1e-9)o.z = -(ba * o.x + bd) / bc;
        else o.z = X.z;
    }
    else {
        o.z = X.z;
        ad += ac; bd += bc;
        double da = bb * aa - ab * ba;
        double dd = bb * ad - ab * bd;
        o.x = -dd / da;
        if (fabs(ab) >= 1e-9)o.y = -(aa * o.x + ad) / ab;
        else if (fabs(bb) >= 1e-9)o.y = -(ba * o.x + bd) / bb;
        else o.y = X.y;
    }
    //printline(o, d, 0, /*200*/0, 0.1);
    //计算A，B骨骼点的在交线的投影
    /*cout << fa[a].o.x << ' ' << fa[a].o.y << ' ' << fa[a].o.z << '\n';
    cout << fa[b].o.x << ' ' << fa[b].o.y << ' ' << fa[b].o.z << '\n';*/
    //printline(fa[a].o, fa[a].d, -5, 10, 0.1);
    //printline(fa[b].o, fa[b].d, -5, 10, 0.1);
    point w1 = d;
    point w2 = fa[a].o - o;
    double len = dotproduct(w1, w2) / moudle(w1) / moudle(w1);
    //printline(o, d, 0, len, 1);
    //printline(o, w2, 0, 10, 1);
    point ao = o + d * len;
    w2 = fa[b].o - o;
    len = dotproduct(w1, w2) / moudle(w1) / moudle(w1);
    point bo = o + d * len;


    //cout << ao.x << ' ' << ao.y << ' ' << ao.z << '\n';
    //cout << bo.x << ' ' << bo.y << ' ' << bo.z << '\n';

    point oa = fa[a].o - ao, ob = fa[b].o - bo;
    //printline(fa[a].o, oa, -100, 200, 0.01);
    //printline(fa[b].o, ob, -100, 200, 0.01);
    double angle = acos(dotproduct(oa, ob) / moudle(oa) / moudle(ob));
    point c = chacheng(oa, ob);
    if (dotproduct(c, d) < 0)
        d = d * -1;
    d = d * (1 / moudle(d));
    double tt = angle / (b - a);
    //cerr << angle << '\n';
    int ff = b;

    for (double i = angle; i > 0; i -= tt) {
        point v = d * (1 - cos(-i)) * dotproduct(ob, d) + ob * cos(-i) + chacheng(d, ob) * sin(-i) + bo;
        PM p = PM(v, chacheng(d, v - o));
        //if (a == 0) 
        {
            //cerr << d.x << ' ' << d.y << ' ' << d.z << '\n';
            /*for (double len = 0.1; len <= 1; len += 0.1) {
                cout << v.x + p.d.x / moudle(p.d) * len << ' ' << v.y + p.d.y / moudle(p.d) * len << ' ' << v.z + p.d.z / moudle(p.d) * len << '\n';
            }*/
            //cout << '\n';
        }
        //mp.clear();
        qiedian[ff].clear();
        hao[ff].clear();
        for (int j = 0; j < triangles.size(); j++) {
            point p1 = judgejiao(ff, j, p, triangles[j].a, triangles[j].b);
            point p2 = judgejiao(ff, j, p, triangles[j].a, triangles[j].c);
            point p3 = judgejiao(ff, j, p, triangles[j].b, triangles[j].c);
            if (p1 == point(1000000007, 1000000007, 1000000007)) {
                qiedian[ff].push_back(triangles[j].a);
                qiedian[ff].push_back(triangles[j].b);
                hao[ff].push_back(j);
                hao[ff].push_back(j);
            }
            else if (p2 == point(1000000007, 1000000007, 1000000007)) {
                qiedian[ff].push_back(triangles[j].a);
                qiedian[ff].push_back(triangles[j].c);
                hao[ff].push_back(j);
                hao[ff].push_back(j);
            }
            else if (p3 == point(1000000007, 1000000007, 1000000007)) {
                qiedian[ff].push_back(triangles[j].b);
                qiedian[ff].push_back(triangles[j].c);
                hao[ff].push_back(j);
                hao[ff].push_back(j);
            }
            else {
                if (!(p1 == point(-1000000007, -1000000007, -1000000007))) {
                    qiedian[ff].push_back(p1);
                    hao[ff].push_back(j);
                }
                if (!(p2 == point(-1000000007, -1000000007, -1000000007))) {
                    qiedian[ff].push_back(p2);
                    hao[ff].push_back(j);
                }
                if (!(p3 == point(-1000000007, -1000000007, -1000000007))) {
                    qiedian[ff].push_back(p3);
                    hao[ff].push_back(j);
                }
            }
        }
        ff--;
    }
    //fclose(stdout);
}

void lu(int i, int j) {
    cout << "    " << i << " " << j << endl;
    qiedian[i].push_back(mpp[j]);
    vis[j] = 1;
    for (int p = 0; p < sad[j]; p++) {
        cout << j << " " << p << endl;
        if (mpp.count(ad[j][p]) && vis[ad[j][p]] == 0) {
            lu(i, ad[j][p]);
        }
    }
}

void findlu(int i) {
    mpp.clear();
    vis.clear();
    ss.clear();
    for (int j = 0; j < qiedian[i].size(); j++) {
        mpp[hao[i][j]] = qiedian[i][j];
    }
    for (int j = 0; j < qiedian[i].size(); j++) {
        for (int k = 0; k < sad[hao[i][j]]; k++) {
            if (mpp.count(ad[hao[i][j]][k]))ss[j]++;
        }
    }
    int cnt = 0;
    for (int j = 0; j < qiedian[i].size(); j++) {
        if (ss[j] == 0)cnt++;
    }
    cerr << i << " " << cnt << endl;
    qiedian[i].clear();
    for (int j = 0; j < hao[i].size(); j++) {
        if (ss[j] == 1 && vis[hao[i][j]] == 0) {
            //lu(i, hao[i][j]);
            break;
        }
    }
}

int main() {
    cout << fixed << setprecision(6);
    freopen("../../3Dgujia/3Dgujia/stl.txt", "r", stdin);
    //freopen("CON", "w", stdout);
    double ax, ay, az, bx, by, bz, cx, cy, cz;
    while (cin >> ax >> ay >> az >> bx >> by >> bz >> cx >> cy >> cz) {
        triangles.push_back(triangle(point(ax, ay, az), point(bx, by, bz), point(cx, cy, cz)));
    }
    /*for (int i = 0;i < triangles.size();i++) {
        show(triangles[i]);
    }*/
    cin.clear();
    freopen("../../3Dgujia/3Dgujia/out.txt", "r", stdin);
    while (cin >> ax >> ay >> az) {
        points.push_back(point(ax, ay, az));
    }
    /*for (int i = 0;i < points.size();i++) {
        cout << points[i].x << " " << points[i].y << " " << points[i].z << endl;
    }*/
    //cout << points.size();
    cerr << triangles.size() << endl;
    for (int i = 0; i < triangles.size() ; i++) {
        for (int j = 0; j < triangles.size(); j++) {
            if (i == j)continue;
            point a1, a2, b1, b2, c1, c2;
            a1 = triangles[i].a;
            b1 = triangles[i].b;
            c1 = triangles[i].c;
            a2 = triangles[j].a;
            b2 = triangles[j].b;
            c2 = triangles[j].c;
            int flag = 0;
            if (a1 == a2 && b1 == b2)flag = 1;
            if (a1 == a2 && b1 == c2)flag = 1;
            if (a1 == b2 && b1 == a2)flag = 1;
            if (a1 == b2 && b1 == c2)flag = 1;
            if (a1 == c2 && b1 == a2)flag = 1;
            if (a1 == c2 && b1 == b2)flag = 1;

            if (a1 == a2 && c1 == b2)flag = 1;
            if (a1 == a2 && c1 == c2)flag = 1;
            if (a1 == b2 && c1 == a2)flag = 1;
            if (a1 == b2 && c1 == c2)flag = 1;
            if (a1 == c2 && c1 == a2)flag = 1;
            if (a1 == c2 && c1 == b2)flag = 1;

            if (b1 == a2 && c1 == b2)flag = 1;
            if (b1 == a2 && c1 == c2)flag = 1;
            if (b1 == b2 && c1 == a2)flag = 1;
            if (b1 == b2 && c1 == c2)flag = 1;
            if (b1 == c2 && c1 == a2)flag = 1;
            if (b1 == c2 && c1 == b2)flag = 1;

            if (flag) {
                for (int k = 0; k < sad[i]; k++) {
                    if (ad[i][k] == j)flag = 0;
                }
            }
            if (flag) {
                ad[i][sad[i]] = j;
                ad[j][sad[j]] = i;
                //if(i==2198)
                    //cout << sad[i] << " " << ad[i][sad[i]] << endl;
                sad[i]++;
                sad[j]++;
                //if (j == 2198) 
                {
                    //cout << i << " " << j << " " <<sad[j]<<" "<<ad[j][sad[j] - 1] << endl;
                }
            }
            /*if (ad[2198][0] != 0) {
                cout << "!!!!!! " << i << " " << j << endl;
                break;
            }*/
        }
    }
    //cout << "!" << ad[2198][0] << endl;
    //cout << "!" << sad[2198] << endl;
    for (int i = 0; i < triangles.size(); i++) {
        if (sad[i] > 3) {
            cout << i << endl;
            show(triangles[i]);
            for (int j = 0; j < sad[i]; j++) {
                cout << ad[i][j] << endl;
                show(triangles[ad[i][j]]);
            }
            break;
        }
    }
    //456 3788
    double d = 4;
    np.push_back(points[0]);
    for (int i = 1; i < points.size(); i++) {
        if (dis(np.back(), points[i]) >= d) {
            np.push_back(points[i]);
        }
    }
    //cout << np.size() << endl;
    for (int i = 0; i < np.size(); i++) {
        if (i == 0) {
            fa.push_back(PM(np[i + 1].x - np[i].x, np[i + 1].y - np[i].y, np[i + 1].z - np[i].z, np[i].x, np[i].y, np[i].z));
        }
        else if (i == np.size() - 1) {
            fa.push_back(PM(np[i].x - np[i - 1].x, np[i].y - np[i - 1].y, np[i].z - np[i - 1].z, np[i].x, np[i].y, np[i].z));
        }
        else {
            fa.push_back(PM(np[i + 1].x - np[i - 1].x, np[i + 1].y - np[i - 1].y, np[i + 1].z - np[i - 1].z, np[i].x, np[i].y, np[i].z));
        }
    }

    int n2 = np.size();
    //cout << n2 << endl;
    for (int i = 0; i < np.size(); i++) {
        //mp.clear();
        for (int j = 0; j < triangles.size(); j++) {
            point p1 = judgejiao(i, j, fa[i], triangles[j].a, triangles[j].b);
            point p2 = judgejiao(i, j, fa[i], triangles[j].a, triangles[j].c);
            point p3 = judgejiao(i, j, fa[i], triangles[j].b, triangles[j].c);
            if (p1 == point(1000000007, 1000000007, 1000000007)) {
                qiedian[i].push_back(triangles[j].a);
                qiedian[i].push_back(triangles[j].b);
                hao[i].push_back(j);
                hao[i].push_back(j);
            }
            else if (p2 == point(1000000007, 1000000007, 1000000007)) {
                qiedian[i].push_back(triangles[j].a);
                qiedian[i].push_back(triangles[j].c);
                hao[i].push_back(j);
                hao[i].push_back(j);
            }
            else if (p3 == point(1000000007, 1000000007, 1000000007)) {
                qiedian[i].push_back(triangles[j].b);
                qiedian[i].push_back(triangles[j].c);
                hao[i].push_back(j);
                hao[i].push_back(j);
            }
            else {
                if (!(p1 == point(-1000000007, -1000000007, -1000000007))) {
                    qiedian[i].push_back(p1);
                    hao[i].push_back(j);
                }
                if (!(p2 == point(-1000000007, -1000000007, -1000000007))) {
                    qiedian[i].push_back(p2);
                    hao[i].push_back(j);
                }
                if (!(p3 == point(-1000000007, -1000000007, -1000000007))) {
                    qiedian[i].push_back(p3);
                    hao[i].push_back(j);
                }
            }
        }
    }
    //freopen("../../3Dgujia/3Dgujia/cp.txt", "w", stdout);
    /*for (int i = 0; i < np.size(); i++) {
        for (int j = 0; j < qiedian[i].size(); j++) {
            if (i <= 20&&i>=18)
            {
                cout << qiedian[i][j].x << " " << qiedian[i][j].y << " " << qiedian[i][j].z << endl;
            }
        }
    }*/
    for (int i = 1; i < np.size() - 1; i++) {
        //mp.clear();
        int sum = qiedian[i].size();
        int flag1 = 1, flag2 = 1;
        if (jiaoshe(i, i - 1) != 0)flag1 = 0;
        if (jiaoshe(i, n2 - 1) != sum)flag2 = 0;
        if (flag1 && flag2) {//前后均无干涉
            cerr << i << endl;
            //for (int j = 0; j < qiedian[i].size(); j++) {
                //cout << qiedian[i][j].x << " " << qiedian[i][j].y << " " << qiedian[i][j].z << endl;
            //}
        }
        else if (flag1 && !flag2) {//前无干涉，后有干涉
            arccut(i - 1, n2 - 1);
            i = n2 - 1;
        }
        else if (!flag1 && flag2) {//前有干涉，后无干涉
            int a, b;
            a = i - 1;
            while (i < n2) {
                i++;
                int f1 = 1, f2 = 1;
                if (jiaoshe(i, a) != 0)f1 = 0;
                //cout << i << " " << n2 << endl;
                if (jiaoshe(i, n2 - 1) != qiedian[i].size())f2 = 0;
                if (f1 && f2) {
                    b = i;
                    break;
                }
                else if (f1 && !f2) {
                    b = n2 - 1;
                    i = n2 - 1;
                    break;
                }
            }
            if (i >= n2) {
                b = n2 - 1;
            }
            arccut(a, b);
        }
    }
    //fclose(stdout);
    /*for (int i = 0; i < np.size(); i++) {
        cerr << hao[i].size() << " " << qiedian[i].size() << endl;
    }*/
    for (int i = 0; i < np.size(); i++) {
        //findlu(i);
    }//求路径
    freopen("../../3Dgujia/3Dgujia/qiepian.txt", "w", stdout);
    for (int i = 0; i < np.size(); i++) {
        for (int j = 0; j < qiedian[i].size(); j++) {
            //if (i <= 0)
            {
                cout << qiedian[i][j].x << " " << qiedian[i][j].y << " " << qiedian[i][j].z << endl;
            }
        }
    }
    return 0;
}