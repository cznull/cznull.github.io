template <int n>
struct doublen
{
	double x[n];
};

template <int n>
double dot(doublen<n> x1, doublen<n> x2);

template <int n>
double norm(doublen<n> x1);

template <int n>
doublen<n> operator-(doublen<n> x1, doublen<n> x2);

template <int n>
doublen<n> operator+(doublen<n> x1, doublen<n> x2);

struct double2 {
	double x, y;
};

double dot(double2 x1, double2 x2);

double norm(double2 x1);

double vol(double2 x1, double2 x2);

double2 operator-(double2 x1, double2 x2);

double2 operator+(double2 x1, double2 x2);

double2 operator*(double s, double2 x);

double2 operator*(double2 x, double s);

struct double3 {
	double x, y, z;
};

double dot(double3 x1, double3 x2);

double norm(double3 x1);

double vol(double3 x1, double3 x2, double3 x3);

double3 cross(double3 x1, double3 x2);

double3 operator-(double3 x1, double3 x2);

double3 operator+(double3 x1, double3 x2);

double3 operator*(double3 x, double s);

struct float3 {
	float x, y, z;
	float3(float a, float b, float c) :x(a),y(b),z(c){
	}

	float3(double3 vertex) {
		x = vertex.x;
		y = vertex.y;
		z = vertex.z;
	}
};

float dot(float3 x1, float3 x2);

float norm(float3 x1);

float vol(float3 x1, float3 x2, float3 x3);

float3 cross(float3 x1, float3 x2);

float3 operator-(float3 x1, float3 x2);

float3 operator+(float3 x1, float3 x2);