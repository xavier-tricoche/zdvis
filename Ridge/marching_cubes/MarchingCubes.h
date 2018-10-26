#include "Vec3f.hpp"

struct TriMeshFace
{
    TriMeshFace() {}
    TriMeshFace(uint I0, uint I1, uint I2)
    {
        I[0] = I0;
        I[1] = I1;
        I[2] = I2;
    }

    uint I[3];
};

struct GRIDCELL {
   Vec3f p[8];		//position of each corner of the grid in world space
   float val[8];	//value of the function at this grid corner
};

//given a grid cell, returns the set of triangles that approximates the region where val == 0.
int Polygonise(GRIDCELL &Grid, TriMeshFace *Triangles, int &NewVertexCount, Vec3f *Vertices);