#ifndef VECTOR_H
#define VECTOR_H

struct vector2D{
	double x;
	double y;
};

struct vector3D{
	double x;
	double y;
	double z;
};

struct tensor22{
	double xx;
	double yx;
	double xy;
	double yy;
};

struct tensor33{
	double xx;
	double yx;
	double zx;
	double xy;
	double yy;
	double zy;
	double xz;
	double yz;
	double zz;
};


struct tensor3S{
	double xx;
	double yy;
	double zz;
	double xy;
	double yz;
	double zx;
};
#endif
