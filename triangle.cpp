/**
 * Parallel Delaunay Triangulation
 * Name 1(eleepalo), Name 2(darrison)
 */

#include "triangle.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <set>
#include <utility>

//!#include <omp.h>
#include <unistd.h>

void print_stats(const std::vector<std::vector<int>> &occupancy)
{
  int max_occupancy = 0;
  long long total_cost = 0;

  for (const auto &row : occupancy)
  {
    for (const int count : row)
    {
      max_occupancy = std::max(max_occupancy, count);
      total_cost += count * count;
    }
  }

  std::cout << "Max occupancy: " << max_occupancy << '\n';
  std::cout << "Total cost: " << total_cost << '\n';
}


//? Source for orientation math: https://www.cs.cmu.edu/~quake/robust.html
//* returns positive if CCW, negative if CW, 0 otherwise (colinear)
float orientation(const Point& a, const Point& b, const Point& c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

//? Source for inCircle math: https://www.cs.cmu.edu/~quake/robust.html
bool inCircle(int v, Triangle t, const std::vector<Point>& V) {
  assert(v >= 0 && v < V.size());
  assert(t.x >= 0 && t.x < V.size());
  assert(t.y >= 0 && t.y < V.size());
  assert(t.z >= 0 && t.z < V.size());

  if (v == t.x || v == t.y || v == t.z) {
    return false;
  }

  const Point& p = V[v];
  const Point& a = V[t.x] - p;
  const Point& b = V[t.y] - p;
  const Point& c = V[t.z] - p;

  float det = ((a.x * a.x + a.y * a.y) * (b.x * c.y - b.y * c.x)
             - (b.x * b.x + b.y * b.y) * (a.x * c.y - a.y * c.x)
             + (c.x * c.x + c.y * c.y) * (a.x * b.y - a.y * b.x));

  float orient = orientation(a, b, c);

  if (orient > 0){
    return det > 0;
  } else {
    return det < 0;
  }
}

void E(Triangle t, const std::vector<Point>& V) {
  t.E.clear();
  for (int i = 0; i < V.size(); i++){
    if (i == t.x || i == t.y || i == t.z) continue;

    if (inCircle(i, t, V)){
      t.E.push_back(i);
    }
  }
}

//! In future might want to ensure all points in input are represented in output!
bool isDelaunay(const Mesh& M, const std::vector<Point>& V){
  for (const Triangle& t : M.triangles){
    if (!t.active) continue;

    const Point& a = V[t.x];
    const Point& b = V[t.y];
    const Point& c = V[t.z];

    //* ensure triangle isn't degenerate
    if (orientation(a,b,c) == 0.0f){
      printf("Degenerate triangle! %d %d %d\n", t.x, t.y, t.z);
      return false;
    }

    //* check circumcircle property
    for (int i = 0; i < V.size(); i++){
      if (i == t.x || i == t.y || i == t.z) continue;

      if (inCircle(i, t, V)){
        //* 
        printf("Point %d inside of circumcircle of triangle %d %d %d!\n", i, t.x, t.y, t.z);
        return false;
      }
    }
  }
  return true;
}

//detaches triangle t on face f
void clearNeighbor(Triangle& t, Face f) {
  //if f is edge (x,y)
  if ((t.x == f.a && t.y == f.b) || (t.x == f.b && t.y == f.a)) {
    t.nbr_xy = -1;
  }
  else if ((t.y == f.a && t.z == f.b) || (t.y == f.b && t.z == f.a)) {
    t.nbr_yz = -1;
  }
  else if ((t.z == f.a && t.x == f.b) || (t.z == f.b && t.x == f.a)) {
    t.nbr_zx = -1;
  }
}

//set t's neighbor as triangle M.triangles[neighborInd] on face f
void setNeighbor(Triangle& t, Face f, int neighborInd) {
  //if f is edge (x,y)
  if ((t.x == f.a && t.y == f.b) || (t.x == f.b && t.y == f.a)) {
    t.nbr_xy = neighborInd;
  }
  else if ((t.y == f.a && t.z == f.b) || (t.y == f.b && t.z == f.a)) {
    t.nbr_yz = neighborInd;
  }
  else if ((t.z == f.a && t.x == f.b) || (t.z == f.b && t.x == f.a)) {
    t.nbr_zx = neighborInd;
  }
}

//returns true if x is in vec
bool contains(const std::vector<int>& vec, int x) {
  return std::find(vec.begin(), vec.end(), x) != vec.end();
}

//returns true if edge (a,b) and (c,d) are the same edge
bool sameEdge(int a, int b, int c, int d) {
  return (a == c && b == d) || (a == d && b == c);
}

//returns 0 if t's edge (x,y) is same edges as (a,b)
//returns 1 if t's edge (y,z) is same edges as (a,b)
//returns 2 if t's edge (z,x) is same edges as (a,b)
//returns -1 if none of t's edges is same edges as (a,b)
int getEdgeInd(const Triangle& t, int a, int b) {
  //edge (x,y)
  if ((t.x == a && t.y == b) || (t.x == b && t.y == a)) {
    return 0;
  }
  //edge (y,z)
  if ((t.y == a && t.z == b) || (t.y == b && t.z == a)) {
    return 1;
  }
  //edge (z,x) 
  if ((t.z == a && t.x == b) || (t.z == b && t.x == a)) {
    return 2;
  }

  return -1; //none of t's edges match edge (a,b)
}


void setNeighborByEdgeInd(Triangle& t, int edgeInd, int neighborInd) {
  if (edgeInd == 0) {
    t.nbr_xy = neighborInd;
  }
  else if (edgeInd == 1) {
    t.nbr_yz = neighborInd;
  }
  else if (edgeInd == 2) {
    t.nbr_zx = neighborInd;
  }

  //nothing if edgeInd == -1
}

bool sharedEdgeInfo(const Triangle& t1, const Triangle& t2, int& t1EdgeInd, int& t2EdgeInd) {
  //check against t1 edge (x,y)
  t2EdgeInd = getEdgeInd(t2,t1.x,t1.y);
  if (t2EdgeInd != -1) {
    t1EdgeInd = 0; //indicates t2 is t1's xy neighbor
    return true;
  }

  //check against t1 edge (y,z)
  t2EdgeInd = getEdgeInd(t2,t1.y,t1.z);
  if (t2EdgeInd != -1) {
    t1EdgeInd = 1; //indicates t2 is t1's yz neighbor
    return true;
  }

  //check against t1 edge (z,x)
  t2EdgeInd = getEdgeInd(t2,t1.z,t1.x);
  if (t2EdgeInd != -1) {
    t1EdgeInd = 2; //indicates t2 is t1's zx neighbor
    return true;
  }

  return false; //if none of t1 and t2's edges match
}

void connectIfNeighbors(int t1Ind, int t2Ind, Mesh& M) {
  int e1, e2;

  if (sharedEdgeInfo(M.triangles[t1Ind], M.triangles[t2Ind], e1, e2)) {
    setNeighborByEdgeInd(M.triangles[t1Ind], e1, t2Ind);
    setNeighborByEdgeInd(M.triangles[t2Ind], e2, t1Ind);
  }
}


int replaceBoundary(int t0Ind, Face f, int tInd, int v, Mesh& M, std::vector<Point>& V) {

  Triangle& t = M.triangles[tInd];

  //make new triangle tt, from face f and point v
  Triangle tt;
  tt.x = f.a;
  tt.y = f.b;
  tt.z = v;

  tt.nbr_xy = -1;
  tt.nbr_yz = -1;
  tt.nbr_zx = -1;
  tt.active = true;
  
  //only check the points in E(t0) and E(t) to make E(tt)
  for (int i = 0; i < t.E.size(); i++) {
    int p = t.E[i];
    if (inCircle(p,tt,V)) {
      tt.E.push_back(p);
    }
  }

  if (t0Ind != -1) { //check there is a t0
    Triangle& t0 = M.triangles[t0Ind];
    for (int i = 0; i < t0.E.size(); i++) {
      int p = t0.E[i];

      //avoid duplicates
      bool duplicate = false;
      for (int j = 0; j < tt.E.size(); j++) {
        if (tt.E[j] == p) {
          duplicate = true;
          break;
        }
      }

      if (!duplicate && inCircle(p, tt,V)) {
        tt.E.push_back(p);
      }
    }
  }

  //add tt to M
  M.triangles.push_back(tt);
  int ttInd = M.triangles.size() -1;

  //detach t from ONLY face f here in M
  //this means: before t0 is adjacent to t across f
  //after t0 is adjacent to tt across f

  clearNeighbor(M.triangles[tInd], f);

  //reconnect t0 to tt across face f
  if (t0Ind != -1) {
    setNeighbor(M.triangles[t0Ind], f, ttInd);
    setNeighbor(M.triangles[ttInd],f,t0Ind);
  }

  //will do the actual detaching, neigbor stuff later in main (with all the new triangles)

  return ttInd;
}

//find cavity R: all active triangles whose E set has v's index
//find boundary faces of R: for each triangle in R, check its 3 edges, an edge on
//boundary if no neighbor across that edge or its not in R
//for each such face f, record the inside triangle t and the outside neighbor t0
//for each face: create new tri tt and get E(tt). add to the mesh

//after all these new tris added, reconnect neighbors: each outside triangle t0 should point across
//face f to corresponding new triangle, and new triangle shoud point back to t0
//connect new triangles to each other along the edges that include v, since adjacent boundary faces
//produce adjacent new triangles
//then mark all triangles in R active = false

//DELAUNAY END**


int main(int argc, char *argv[])
{
  const auto init_start = std::chrono::steady_clock::now();

  std::string input_filename;
  int num_threads = 0;
  double SA_prob = 0.1;
  int SA_iters = 5;
  char parallel_mode = '\0';
  int batch_size = 1;

  int opt;
  while ((opt = getopt(argc, argv, "f:n:p:i:m:b:")) != -1)
  {
    switch (opt)
    {
    case 'f':
      input_filename = optarg;
      break;
    case 'n':
      num_threads = atoi(optarg);
      break;
    case 'p':
      SA_prob = atof(optarg);
      break;
    case 'i':
      SA_iters = atoi(optarg);
      break;
    case 'm':
      parallel_mode = *optarg;
      break;
    case 'b':
      batch_size = atoi(optarg);
      break;
    default:
      std::cerr << "Usage: " << argv[0]
                << " -f input_filename -n num_threads [-p SA_prob] [-i "
                   "SA_iters] -m parallel_mode -b batch_size\n";
      exit(EXIT_FAILURE);
    }
  }

  // Check if required options are provided
  if (empty(input_filename) || num_threads <= 0 || SA_iters <= 0 ||
      (parallel_mode != 'A' && parallel_mode != 'W') || batch_size <= 0)
  {
    std::cerr << "Usage: " << argv[0]
              << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] "
                 "-m parallel_mode -b batch_size\n";
    exit(EXIT_FAILURE);
  }

  std::cout << "Number of threads: " << num_threads << '\n';
  std::cout << "Simulated annealing probability parameter: " << SA_prob << '\n';
  std::cout << "Simulated annealing iterations: " << SA_iters << '\n';
  std::cout << "Input file: " << input_filename << '\n';
  std::cout << "Parallel mode: " << parallel_mode << '\n';
  std::cout << "Batch size: " << batch_size << '\n';

  std::ifstream fin(input_filename);

  if (!fin)
  {
    std::cerr << "Unable to open file: " << input_filename << ".\n";
    exit(EXIT_FAILURE);
  }


  //ADDED FOR DELAUNAY START
  int n; //number of points
  fin >> n;

  std::vector<Point> V(n);

  for (int i = 0; i < n; i++) {
    fin >> V[i].x >> V[i].y;
  }

  std::cout << "Read " << n << " points:\n";
  for (int i = 0; i < n; i++)
  {
    std::cout << i << ": (" << V[i].x << ", " << V[i].y << ")\n";
  }

  const double init_time =
    std::chrono::duration_cast<std::chrono::duration<double>>(
      std::chrono::steady_clock::now() - init_start).count();
      std::cout << "Initialization time (sec): " << std::fixed
      << std::setprecision(10) << init_time << '\n';

  const auto compute_start = std::chrono::steady_clock::now();

  //* Timing Code Start


  Mesh M;

  Triangle tb;
  //!! JUST FOR DEBUGGING
  tb.x = 0;
  tb.y = 1;
  tb.z = 2;
  //!!make tb bounding triangle - convex hull

  tb.nbr_xy = -1;
  tb.nbr_yz = -1;
  tb.nbr_zx = -1;
  tb.active = true;

  //tb's encroach set is all of V
  for (int i = 0; i < n; i++) {
    tb.E.push_back(i);
  }

  //M = {tb}
  M.triangles.push_back(tb);

  //iterate through all points: V[i]
  for (int i = 0; i < n; i++) { //index of corresponding point into V
    std::vector<int> R;

    //build cavity R
    for (int j = 0; j < M.triangles.size(); j++) {
      if (!M.triangles[j].active) continue;

      if (contains(M.triangles[j].E,i)) { //if point i is in E(j), then add j to R 
        R.push_back(j);
      }
    }

    std::vector<int> newTriangles;

    for (int t = 0; t < R.size(); t++) { //for each triangle in R, check which are faces - if is, do replacing
      int tInd = R[t];

      //edge is a face is no other triangle uses it (-1)
      //or if only triangle outside of R shares it
      if (M.triangles[tInd].nbr_xy == -1 || !contains(R,M.triangles[tInd].nbr_xy)) { 
        Face f;
        f.a = M.triangles[tInd].x;
        f.b = M.triangles[tInd].y;
        int t0Ind = M.triangles[tInd].nbr_xy;

        newTriangles.push_back(replaceBoundary(t0Ind, f, tInd, i, M, V));
      }
      if (M.triangles[tInd].nbr_yz == -1 || !contains(R,M.triangles[tInd].nbr_yz)) {
        Face f;
        f.a = M.triangles[tInd].y;
        f.b = M.triangles[tInd].z;
        int t0Ind = M.triangles[tInd].nbr_yz;

        newTriangles.push_back(replaceBoundary(t0Ind, f, tInd, i, M, V));
      }
      if (M.triangles[tInd].nbr_zx == -1 || !contains(R,M.triangles[tInd].nbr_zx)) {
        Face f;
        f.a = M.triangles[tInd].z;
        f.b = M.triangles[tInd].x;
        int t0Ind = M.triangles[tInd].nbr_zx;

        newTriangles.push_back(replaceBoundary(t0Ind, f, tInd, i, M, V));
      }
    }

    //set neighbors for all the added triangles
    //they should only be with each other (the outside ones were taken care of in replaceBoundary already)
    for (int a = 0; a < newTriangles.size(); a++) {
      for (int b = a+1; b < newTriangles.size(); b++) {
        connectIfNeighbors(newTriangles[a], newTriangles[b], M);
      }
    }

     //deactivate all the triangles in R
     
     for (int t : R) {
      M.triangles[t].active = false;
      //is marking neighbors -1 needed? / wont hurt anything right?
      M.triangles[R[t]].nbr_xy = -1;
      M.triangles[R[t]].nbr_yz = -1;
      M.triangles[R[t]].nbr_zx = -1;
     }


    //need to do detaching for all triangles in R -- or should i still be doing this in replaceBoundary
    //plus update neighbors for all triangles...bc of new triangles
  }

  //return M

  //ADDED FOR DELAUNAY END


  /* Initialize any additional data structures needed in the algorithm */

  // Student code end





  //* Timing Code End

  const double compute_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - compute_start)
          .count();
  std::cout << "Computation time (sec): " << compute_time << "\n";

  bool valid = isDelaunay(M, V);

  if (valid){
    std::cout << "Valid Delaunay Triangulation :)\n";
  } else {
    std::cout << "Invalid Delaunay Triangulation :( \n";
  }


  std::ofstream out("outputs/dummy.txt");
  if (!out) {
      std::cerr << "Failed to open outputs/dummy.txt\n";
      exit(EXIT_FAILURE);
  }

  for (auto t : M.triangles){
    out << t.x << " " << t.y << " " << t.z << "\n";
  }
  
  out.close();

}

