#include<iostream>
#include<fstream>
#include<vector>
#include<array>
using namespace std;
const double PI = acos(-1);
struct Point {
	array<double, 3>a;
	double& x = a[0];
	double& y = a[1];
	double& z = a[2];
	Point() { a.fill(NAN); }
	Point(const double& x, const double& y, const double& z) { 
		this->a.at(0) = x;
		this->a.at(1) = y;
		this->a.at(2) = z;
	}
	friend double dot(const Point& a, const Point& b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}
	friend double crossyz(const Point& a, const Point& b) {
		return a.y * b.z - a.z * b.y;
	}
	Point operator-(const Point& p)const& {
		return Point(x - p.x, y - p.y, z - p.z);
	}
	Point operator*(const double& d)const& {
		return Point(x * d, y * d, z * d);
	}
	Point rotate(const double& a) const& {
		return Point(x, y * cos(a) + z * sin(a), -y * sin(a) + z * cos(a));
	}
	friend ofstream& operator<<(ofstream& out,const Point&p) {
		out << "X" << p.x << " Y" << p.y;
		return out;
	}
	void operator=(const Point& p)&{ 
		a = p.a;
	}
};
struct PM {
	Point o, d;
	PM(){}
	PM(Point o, Point d) :o(o), d(d) {}
	double disToPoint(const Point& p) {
		return dot((p - o), d) / sqrt(dot(d, d));
	}
};
ofstream gcode;
size_t layer_count;
vector<vector<Point>>layers;
vector<PM>pms;
Point move1;
double extruder, g10 = 1;
double layer_height, layer_width = 1;
double R;
void readRoute() {
	ifstream input;
	input.open("../../3Dgujia/3Dgujia/qiepian.txt", ios::in);
	input >> layer_count;
	for (size_t i = 0; i < layer_count; ++i) {
		size_t n; input >> n;
		vector<Point>layer;
		for (size_t j = 0; j < n; ++i) {
			Point in; input >> in.x >> in.y >> in.z;
			layer.emplace_back(in);
		}
		layers.emplace_back(layer);
		PM p;
		input >> p.o.x >> p.o.y >> p.o.z;
		input >> p.d.x >> p.d.y >> p.d.z;
		pms.emplace_back(p);
	}
	input.close();
}

inline double extrudedLength(const double& distance, const double& height) {
	return distance * height * layer_width;
}

void printGcode(const vector<Point>& l, PM p) {
	gcode << "G0 ";
	gcode << l.front() << " Z" << l.front().z << '\n';
	extruder += g10;
	gcode << "G1 E" << extruder << '\n';
	for (size_t i = 1; i < l.size(); ++i) {
		gcode << "G1 ";
		gcode << (l[i] - move1).rotate(R) << ' ';
		extruder += extrudedLength(sqrt(dot(l[i] - l[i - 1], l[i] - l[i - 1])), p.disToPoint(l[i]));
		gcode << "E" << extruder << '\n';
	}
	extruder -= g10;
	gcode << "G1 E" << extruder << '\n';
}

int main() {
	gcode.open("D.gcode", ios::out | ios::binary);
	gcode << ";Layer height:0.5" << '\n';
	gcode << "G0 X100 Y304\nG92 Y0 Z85.5 E0\n";
	gcode << ";LAYER_COUNT:" << layer_count << '\n';
	pms.emplace_back(PM(Point(1, 1, 0), Point(0, 0, 1)));
	readRoute();
	move1 = layers[0][0] - Point(0, 0, -85.5);
	for (size_t i = 1; i <= layer_count; ++i) {
		gcode << ";LAYER:" << i << '\n';
		double r = asin(crossyz(pms[i - 1].d, pms[i].d));
		gcode << "G1A " << r * 180 / PI / 0.72;
		R += r;
	}
	gcode << ";TIME_ELAPSED:???\n";
}