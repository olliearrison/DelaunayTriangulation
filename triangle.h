/**
 * Delaunay Triangulation
 * Elisenda Lee-Palou(eleepalo), Ollie Arrison(darrison)
 */

#ifndef __WIREOPT_H__
#define __WIREOPT_H__

#include <cstdint>
//!#include <omp.h>
#include <vector>
#include <set>
#include <climits>
#include <unordered_map>

#define MAX_PTS_PER_WIRE 5
#define COST_REPORT_DEPTH 10

/** README(student):

*/
// struct validate_wire_t {
//   uint8_t num_pts;
//   struct {
//     uint16_t x;
//     uint16_t y;
//   } p[MAX_PTS_PER_WIRE];
//   validate_wire_t &cleanup(void);
//   void print_wire(void) const;
// };


struct Point {
  float x, y;
};

inline Point operator+(Point a, Point b){
    return Point{a.x + b.x, a.y + b.y};
}

inline Point operator-(Point a, Point b){
    return Point{a.x - b.x, a.y - b.y};
}

//for incremental algo
struct Triangle {
  int x, y, z;
  std::vector<int> E; //encroach set

  //DONT USE FOR PARALLEL VERSION
  //neighbors
  //index of the neighbor triangle in the mesh vector M
  //-1 if no neighbor
  int nbr_xy; //neighbor across edge (x,y)
  int nbr_yz;
  int nbr_zx;

  bool active;
};

struct Face {
 int a, b;
 
 //* page 14: https://www.cs.cmu.edu/~guyb/paralg/papers/BlellochGuShunSun.pdf
 //* "The hash-map is indexed on the d corners of a face in some canonical order."
 //* face init with min than max points
 Face(int x, int y) : a(std::min(x,y)), b(std::max(x,y)) {}
 
 bool operator==(const Face& other) const {
 return a == other.a && b == other.b;
 }
};

struct FaceHash {
 size_t operator()(const Face& f) const {
 return std::hash<long long>()(((long long)f.a << 32) | (unsigned int)f.b);
 }
};

struct Mesh {
 int n;
 std::vector<Triangle> triangles;
 std::unordered_map<Face, std::pair<int,int>, FaceHash> face_to_tri;
};

struct BoundaryResult {
 int ttInd;
 Face oldFace;
 Face side1; //* (a,v)
 Face side2; //* (b,v)
 int tInd;
 int t0Ind;
};

#include <iostream>

inline void printMesh(const Mesh& M, const std::vector<Point>& V) {
//  std::cout << "=== Mesh (" << M.triangles.size() << " total, ";
 int active = 0;
 for (const auto& t : M.triangles) if (t.active) active++;
//  std::cout << active << " active) ===\n";

 for (int i = 0; i < (int)M.triangles.size(); i++) {
 const Triangle& t = M.triangles[i];
 if (!t.active){
 // std::cout << "INNACTIVE T" << i << ": ("
 // << t.x << "[" << (int)V[t.x].x << "," << (int)V[t.x].y << "], "
 // << t.y << "[" << (int)V[t.y].x << "," << (int)V[t.y].y << "], "
 // << t.z << "[" << (int)V[t.z].x << "," << (int)V[t.z].y << "])"
 // << " E={";

 std::cout << " N/A: T" << i << ": ("
 << t.x << ", "
 << t.y << ", "
 << t.z << ") "
 << " E={";
 for (int j = 0; j < (int)t.E.size(); j++) {
 std::cout << t.E[j];
 if (j + 1 < (int)t.E.size()) std::cout << ",";
 }
 std::cout << "}\n";
 } else {
 // std::cout << " T" << i << ": ("
 // << t.x << "[" << (int)V[t.x].x << "," << (int)V[t.x].y << "], "
 // << t.y << "[" << (int)V[t.y].x << "," << (int)V[t.y].y << "], "
 // << t.z << "[" << (int)V[t.z].x << "," << (int)V[t.z].y << "])"
 // << " E={";

 std::cout << " T" << i << ": ("
 << t.x << ", "
 << t.y << ", "
 << t.z << ") "
 << " E={";
 for (int j = 0; j < (int)t.E.size(); j++) {
 std::cout << t.E[j];
 if (j + 1 < (int)t.E.size()) std::cout << ",";
 }
 std::cout << "}\n";
 }
 }

 std::cout << " face_to_tri (" << M.face_to_tri.size() << " faces):\n";
 for (const auto& [f, slots] : M.face_to_tri) {
 std::cout << " (" << f.a << "," << f.b << ") -> ["
 << slots.first << ", " << slots.second << "]\n";
 }

}

struct Task {
 int t0; // outside neighbor (-1 if boundary)
 int t; // inside triangle
 Face f; // boundary face
 int v; // point to insert

 Task(int _t0, int _t, Face _f, int _v)
 : t0(_t0), t(_t), f(_f), v(_v) {}

 Task() : t0(-1), t(-1), f(0,0), v(-1) {}

 bool operator<(const Task& other) const {
 if (v != other.v) return v < other.v;
 if (f.a != other.f.a) return f.a < other.f.a;
 if (f.b != other.f.b) return f.b < other.f.b;
 if (t != other.t) return t < other.t;
 return t0 < other.t0;
 }
};

//for divide and conquer algo
struct TriangleDC {
  int v[3]; //vertices
  TriangleDC* adj[3]; //adj[i] is opposite v[i] 

  TriangleDC(int a, int b, int c) { //a,b,c are indices for their points into V
        v[0] = a;
        v[1] = b;
        v[2] = c;

        adj[0] = nullptr;
        adj[1] = nullptr;
        adj[2] = nullptr;
    }
};
struct TriangulationDC {
  // std::vector<TriangleDC*> triangles;

  //graph[a] contains all point ind connected to a by an active edge
  std::vector<std::vector<int>> graph;

};


// Definition of the wire checker
// struct wr_checker {
//   std::vector<Wire> wires;
//   std::vector<std::vector<int>> occupancies;
//   const int nwires;
//   const int dim_x;
//   const int dim_y;
//   wr_checker(std::vector<Wire> &wires,
//              std::vector<std::vector<int>> &occupancies)
//       : wires(wires), occupancies(occupancies), nwires(wires.size()),
//         dim_x(occupancies[0].size()), dim_y(occupancies.size()) {}
//   void validate() const;
// };

const char *get_option_string(const char *option_name,
                              const char *default_value);
int get_option_int(const char *option_name, int default_value);
float get_option_float(const char *option_name, float default_value);

#endif
