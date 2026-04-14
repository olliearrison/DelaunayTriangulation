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

#include <omp.h>
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

/* This function write the output into 2 files
(1) It write occupancy grids into a file
(2) It convert wires from Wire to validate_wire_t by to_validate_format
(2) It write wires into another file
*/
void write_output(
    const std::vector<Wire> &wires, const int num_wires,
    const std::vector<std::vector<int>> &occupancy, const int dim_x,
    const int dim_y,
    std::string wires_output_file_path = "outputs/wire_output.txt",
    std::string occupancy_output_file_path = "outputs/occ_output.txt")
{

  std::ofstream out_occupancy(occupancy_output_file_path, std::fstream::out);
  if (!out_occupancy)
  {
    std::cerr << "Unable to open file: " << occupancy_output_file_path << '\n';
    exit(EXIT_FAILURE);
  }
  out_occupancy << dim_x << ' ' << dim_y << '\n';

  for (const auto &row : occupancy)
  {
    for (size_t i = 0; i < row.size(); ++i)
      out_occupancy << row[i] << (i == row.size() - 1 ? "" : " ");
    out_occupancy << '\n';
  }
  out_occupancy.close();

  std::ofstream out_wires(wires_output_file_path, std::fstream::out);
  if (!out_wires)
  {
    std::cerr << "Unable to open file: " << wires_output_file_path << '\n';
    exit(EXIT_FAILURE);
  }

  out_wires << dim_x << ' ' << dim_y << '\n';
  out_wires << num_wires << '\n';

  for (const auto &wire : wires)
  {
    // NOTICE: we convert to keypoint representation here, using
    // to_validate_format which need to be defined in the bottom of this file
    validate_wire_t keypoints = wire.to_validate_format();
    for (int i = 0; i < keypoints.num_pts; ++i)
    {
      out_wires << keypoints.p[i].x << ' ' << keypoints.p[i].y;
      if (i < keypoints.num_pts - 1)
        out_wires << ' ';
    }
    out_wires << '\n';
  }

  out_wires.close();
}

//DELAUNAY TRIANGULATION START**



void E(triangle t, const std::vector<Point>& V) {
  t.E.clear();
  for (int i = 0; i < V.size(); i++){
    if (i == t.x || i == t.y || i == t.z) continue;

    if (inCircle(V[i], t, V)){
      t.E.push_back(i);
    }
  }
}

//? Source for orientation math: https://www.cs.cmu.edu/~quake/robust.html
//* returns positive if CCW, negative if CW, 0 otherwise (colinear)
float orientation(const Point& a, const Point& b, const Point& c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

//? Source for inCircle math: https://www.cs.cmu.edu/~quake/robust.html
bool inCircle(int v, Triangle t, const std::vector<Point>& V) {
  const Point& p = V[v];
  const Point& a = V[t.x] - p;
  const Point& b = V[t.y] - p;
  const Point& c = V[t.z] -p;

  float det = ((a.x * a.x + a.y * a.y) + (b.x * c.y - b.y * c.x) -
              (b.x * b.x + b.y * b.y) + (a.x * c.y - a.y * c.x) +
              (c.x * c.x + c.y * c.y) + (a.x * b.y - a.y * b.x));

  float orient = orientation(a, b, c);

  if (orient > 0){
    return det > 0;
  } else {
    return det < 0;
  }
  }


//? Source for inCircle math: https://www.cs.cmu.edu/~quake/robust.html
bool inCircle(int v, triangle t, const std::vector<Point>& V) {
  const Point& p = V[v];
  //* put triangle so p is at origin
  const Point a = V[t.x] - p;
  const Point b = V[t.y] - p;
  const Point c = V[t.z] - p;

  float det = ((a.x * a.x + a.y * a.y) + (b.x * c.y - b.y * c.x) -
              (b.x * b.x + b.y * b.y) + (a.x * c.y - a.y * c.x) +
              (c.x * c.x + c.y * c.y) + (a.x * b.y - a.y * b.x));

  float orient = orientation(a, b, c);

  if (orient > 0){
    return det > 0;
  } else {
    return det < 0;
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
    if (orientation(a) == 0.0f){
      printf("Degenerate triangle! %d %d %d\n", t.x, t.y, t.z);
      return false;
    }

    //* check circumcircle property
    for (int i = 0; i < V.size(); i++){
      if (i == t.x || i == t.y || i == t.z) continue;

      if (inCircle(V[i], t, V)){
        //* 
        printf("Point %d inside of circumcircle of triangle %d %d %d!\n", i, t.x, t.y, t.z);
        return false;
      }
    }
  }
  return true;
}




void setNeighbor(Triangle&t, Face f, int neighborInd) {
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

//in algo, need to make sure to fill in encroach sets for t and t0

void replaceBoundary(int t0Ind, Face f, int tInd, int v, Mesh& M, std::vector<Point>& V) {

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
    if (inCircle(V[p],tt)) {
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

      if (!duplicate && inCircle(V[p], tt)) {
        tt.E.push_back(p);
      }
    }
  }


  //add tt to M
  M.triangles.push_back(tt);
  int ttInd = M.triangles.size() -1;

  //!!how to detach t from face f in M 
  //before: t0 is adjacent to t across f
  //after: t0 is adjacent to tt across f

  //detach t -- more complicated bc need to take care of neighbors...
  //**possibly do this at the end in main: detach all triangles in R...need to know about the
  //other new triangles
  t.active = false; //mark here just to be safe
  //will do the actual detaching, neigbor stuff later in main

  //reconnect t0 to tt across face f
  if (t0Ind != -1) {
    setNeighbor(M.triangles[t0Ind], f, ttInd);
    setNeighbor(M.triangles[ttInd],f,t0Ind);
  }

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

  std::vector<Points> V(n);

  for (int i = 0; i < n; i++) {
    fin >> V[i].x >> V[i].y;
  }

  std::cout << "Read " << n << " points:\n";
  for (int i = 0; i < n; i++)
  {
    std::cout << i << ": (" << V[i].x << ", " << V[i].y << ")\n";
  }


  Mesh M;

  Triangle tb;
  //make tb bounding triangle - convex hull

  //tb's encroach set is all of V
  for (int i = 0; i < n; i++) {
    tb.E.pushback(i);
  }

  //M = {tb}
  M.triangles.push_back(tb);

  //iterate through all points: V[i]
  for (int i = 0; i < n; i++) { //index of corresponding point into V
    //initialize R...
    std::vector<Triangle> R;
    for (Triangle t : M.triangles) {
      //if V[i] in E(t), then add t to R
      for (int j : t.E) { //j is index of corresponding point into V
        if (i == t.E[j]) {
          //add t to R
          R.pushback(t);
        }
      }
    }

    for (Triangle t : R) {
      //edge is a face is no other triangle uses it (-1)
      //or if only triangle outside of R shares it
      if (t.nbr_xy == -1 || !inR(t.nbr_xy,R)) { 
        Face f = (t.x,t.y);
        t0 = t.nbr_xy;
        //make sure references match up...
        replaceBoundary(t0, f, t, i, M, V);
      }
      if (t.nbr_yz == -1 || !inR(t.nbr_yz,R)) {
        Face f = (t.y, t.z);
        t0 = t.nbr_yz;
        replaceBoundary(t0, f, t, i, M, V);
      }
      if (t.nbr_zx == -1 || !inR(t.nbr_zx,R)) {
        Face f = (t.z, t.x);
        t0 = t.nbr_zx;
        replaceBoundary(t0, f, t, i, M, V);
      }
    }

    //need to do detaching for all triangles in R
    //plus update neighbors for all triangles...bc of new triangles
  }

  //return M

  //ADDED FOR DELAUNAY END


  /* Initialize any additional data structures needed in the algorithm */

  // Student code end
  const double init_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - init_start)
          .count();
  std::cout << "Initialization time (sec): " << std::fixed
            << std::setprecision(10) << init_time << '\n';

  const auto compute_start = std::chrono::steady_clock::now();

  /* TODO (student code start): Implement the wire routing algorithm here and
    feel free to structure the algorithm into different functions.
    Don't use global variables.
    Use OpenMP to parallelize the algorithm.
  */





  // Student code end
  // DON'T CHANGE THE FOLLOWING CODE
  const double compute_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - compute_start)
          .count();
  std::cout << "Computation time (sec): " << compute_time << '\n';

  /* wire to run check on wires and occupancy */
  wr_checker checker(wires, occupancy);
  checker.validate();

  /* Write wires and occupancy matrix to files */
  print_stats(occupancy);
  write_output(wires, num_wires, occupancy, dim_x, dim_y);
}

/* TODO (student): implement to_validate_format to convert Wire to
  validate_wire_t keypoint representation in order to run checker and
  write output
*/
validate_wire_t Wire::to_validate_format(void) const
{
  validate_wire_t w;

  // num of keypoints: start + bends + end
  w.num_pts = 2 + numBends;

  // start
  w.p[0].x = start_x;
  w.p[0].y = start_y;

  // bends
  for (int i = 0; i < numBends; i++)
  {
    w.p[i + 1].x = bendsX[i];
    w.p[i + 1].y = bendsY[i];
  }

  // end
  w.p[w.num_pts - 1].x = end_x;
  w.p[w.num_pts - 1].y = end_y;

  return w;
}
