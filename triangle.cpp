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

void E(Triangle &t, const std::vector<Point>& V) {
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

  std::vector<int> found(M.n, 0);

  //* check all triangles are valid delaunay triangles
  for (const Triangle& t : M.triangles){
    if (!t.active) continue;
    if (t.x >= M.n || t.y >= M.n || t.z >= M.n) continue;
    //* ignore any triangles containing points that are part of the initial 
    //* super triangle

    //* ensure all points are found
    found[t.x] ++;
    found[t.y] ++;
    found[t.z] ++;

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

  for (auto t: found){
    if (t == 0){
      return false;
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



//START DIVIDE AND CONQUER 

//orientation test: 2d cross product (b - a) x (c - a)
//
//if < 0: points (a, b, c) are in clockwise order
//if > 0: points (a, b, c) are in counterclockwise order
//if = 0: points are collinear
//
//for our triangleDC implementation, enforce triangle points to be stored in clockwise order - so swap points b, c when this outputs a positive value -- dont need to do anymore since not storing triangles...
double orient2d(const std::vector<Point>& V, const int a, const int b, const int c) {
    return (V[b].x - V[a].x) * (V[c].y - V[a].y)
         - (V[b].y - V[a].y) * (V[c].x - V[a].x);
}

void addGraphEdge(TriangulationDC& T, int a, int b) {
  T.graph[a].insert(b);
  T.graph[b].insert(a);
}

void deleteEdge(TriangulationDC& T, int a, int b) {
  T.graph[a].erase(b);
  T.graph[b].erase(a);
}

//returns true if rCand is inside the circumcircle(l,r,lCand)
//returns false if rCand is outside the circumcircle(l,r,lCand)
bool inCircleDC(const std::vector<Point>& V, int l, int r, int lCand, int rCand) {
  const Point& p = V[rCand]; //point you're checking if its in the circle

  //points making up the circumcircle's triangle
  const Point& a = V[l] - p;
  const Point& b = V[r] - p;
  const Point& c = V[lCand] - p;

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

//initialize subsets (of 2 or 3 points) as edge or "triangle" - automatic triangulation
TriangulationDC triangulateInit(const std::vector<Point>& V, int lo, int hi) {
  TriangulationDC T;
  T.graph.resize(V.size());
  int n = hi - lo; //number of points in subset

  if (n == 2) {
    int p0 = lo;
    int p1 = lo+1;

    addGraphEdge(T,p0,p1);

    return T;
  } else if (n == 3) {
    int p0 = lo;
    int p1 = lo+1;
    int p2 = lo+2;

    //check for clockwise condition
    double orient = orient2d(V, p0, p1, p2);

    if (orient == 0) { 
      //collinear, so degenerate
      addGraphEdge(T, p0, p1);
      addGraphEdge(T, p1, p2);
      return T;
    }

    //if orient < 0 (clockwise) or orient > 0 (counterclockwise)
    //graph edges are undirected, so orientation does not matter
    //dont need to swap points to fit clockwise ordering like you would if storing triangles
    addGraphEdge(T, p0, p1);
    addGraphEdge(T, p1, p2);
    addGraphEdge(T, p2, p0);
    return T;
  } else {
    printf("Error: subset <= 3 but not =2 or =3 ...split incorrectly");
    return T; //but will be incorrect...
  }
}

//assume p is collinear with (a,b)
//returns true if p is between a and b, not just somewhere on the same line
bool between(const std::vector<Point>& V, int a, int b, int p) {
  return std::min(V[a].x, V[b].x) <= V[p].x && V[p].x <= std::max(V[a].x, V[b].x) && std::min(V[a].y, V[b].y) <= V[p].y && V[p].y <= std::max(V[a].y, V[b].y);
}


//finds base LR edge
//i.e. the lower common tangent (edge from point l in L to point r in R s.t. all points in both triangulations L and R are on or above the line l,r)
std::pair<int, int> findBaseLR(const TriangulationDC& L,const TriangulationDC& R,const std::vector<Point>& V) {
  std::vector<int> LV;
  std::vector<int> RV;

  for (int i = 0; i < L.graph.size(); i++) {
    if (!L.graph[i].empty()) {
      LV.push_back(i);
    }
  }

  for (int i = 0; i < R.graph.size(); i++) {
    if (!R.graph[i].empty()) {
      RV.push_back(i);
    }
  }

  //try every possible edge from a left to a right vertex
  for (int l : LV) {
    for (int r : RV) {
      bool isTangent = true;

      //make sure no left point is below the line l,r
      for (int p : LV) {
        if (p == l) continue;

        double o = orient2d(V, l, r, p);

        
        if (o < 0) { //means p is below/to the right of l,r
          isTangent = false; //therefore l,r cant be lower tangent
          break;
        }

        //reject tangent that skips a collinear hull point
        if (o == 0 && between(V, l, r, p)) { //p collinear (o = 0) and is between l and r
          //this means l,r skips over a point that is between them
          isTangent = false; //so reject; want to choose an edge that connects to the closer point instead of skipping it
          break;
        }
      }

      if (!isTangent) continue;

      //check same conditions for all right points
      for (int p : RV) {
        if (p == r) continue;

        double o = orient2d(V, l, r, p);

        //if point below l,r then it is not the lower tangent
        if (o < 0) {
          isTangent = false;
          break;
        }

        //reject tangent that skips a collinear point
        if (o == 0 && between(V, l, r, p)) {
          isTangent = false;
          break;
        }
      }

      //all points are on/above l,r and no collinear point skipped
      if (isTangent) {
        return {l, r}; //base LR-edge
      }
    }
  }
  //no base found...should not happen
  return {-1, -1};
}

//returns index for left candidate point for merge
int findLeftCandidate(TriangulationDC& T, int l, int r, const std::vector<Point>& V) {
  std::vector<std::pair<double, int>> candidates;

  //look only at points directly connected to l
  for (int potCand : T.graph[l]) {
    if (potCand == r) continue; //skip current base LR point

    //vector for base edge l,r
    //"base LR vector"
    double baseX = V[r].x - V[l].x;
    double baseY = V[r].y - V[l].y;

    //vector from l to potCand 
    //"candidate vector"
    double candX = V[potCand].x - V[l].x;
    double candY = V[potCand].y - V[l].y;

    //cross product bw base LR vector and candidate vector
    //gives direction from LR base to candidate
    double cross = baseX * candY - baseY * candX;
    double dot = baseX * candX + baseY * candY; //dot product: helps get angle bw vectors

    //left candidate: counterclockwise from base l,r

    //this is the <180 degree angle check
    //if cross >0: counterclockwise angle (from base LR to l,potCand) bw 0 and 180 (exclusive)
    //if cross=0: angle is 0 or 180 (collinear)
    //if cross <0: angle is >180 (clockwise)
    if (cross <= 0) continue;
    //only want to consider candidates that are to the left of base LR (i.e. counterclockwise)

    double angle = atan2(cross, dot); //since cross>0 in this case, this is counterclockwise angle from l,r to l,potCand
    candidates.push_back({angle, potCand}); //store angle to sort for candidate order
  }

  //sorts by angle - want smallest angle to be first potential candidate
  std::sort(candidates.begin(), candidates.end()); 

  for (int i = 0; i < candidates.size(); i++) {
    int potCand = candidates[i].second;

    //no next potential candidate means potCand is final (bc dont have a next cand to check in circumcircle)
    if (i+1 >= candidates.size()) {
      return potCand;
    } 

    int nextPotCand = candidates[i+1].second;

    //if nextPotCand is inside circumcircle(l, r, potCand), then edge (l, potCand) violates Delaunay
    if (inCircleDC(V, l, r, potCand, nextPotCand)) {
      deleteEdge(T, l, potCand); //delete LL edge
      continue; //next potential candidate will become new potential candidate
    }
    return potCand; //otherwise, potCand passes both <180 and nextPotCand not in circumcircle, so potCand is final candidate
  }
  return -1; //no valid candidates
}

//returns index for right candidate point for merge
int findRightCandidate(TriangulationDC& T, int l, int r,const std::vector<Point>& V) {
  std::vector<std::pair<double, int>> candidates;

  //only look at points directly connected to r
  for (int potCand : T.graph[r]) {
    if (potCand == l) continue;

    double baseX = V[l].x - V[r].x;
    double baseY = V[l].y - V[r].y;

    double candX = V[potCand].x - V[r].x;
    double candY = V[potCand].y - V[r].y;

    //cross product: direction from (r,l) to (r,potCand)
    double cross = baseX * candY - baseY * candX;
    double dot = baseX * candX + baseY * candY;

    //right candidate: clockwise from r,l

    //this is the <180 degree angle check
    //cross < 0: clockwise angle bw 0 and 180 (exclusive)
    //cross =0: angle 0 or 180 (collinear)
    //cross >0: angle > 180 
    if (cross >= 0) continue;

    
    double angle = atan2(-cross, dot); //convert clockwise angle to positive for sorting
    candidates.push_back({angle, potCand});
  }

  //sort by angle
  std::sort(candidates.begin(), candidates.end());

  for (int i = 0; i < candidates.size(); i++) {
    int potCand = candidates[i].second;

    //no next potential candidate means potCand is final (bc dont have a next cand to check in circumcircle)
    if (i+1 >= candidates.size()) {
      return potCand;
    }

    int nextPotCand = candidates[i+1].second;

    //if nextPotCand is inside circumcircle(l, r, potCand), then edge (r, potCand) violates Delaunay
    if (inCircleDC(V, l, r, potCand, nextPotCand)) {
      deleteEdge(T, r, potCand); //delete RR edge
      continue; //nextPotCand becomes new potCand
    }

    //potCand passed the circumcircle test i.e. nextPotCand not in the circumcircle
    return potCand;
  }
  return -1; //no valid candidates
}

//add triangle l,r,cand to edge graph
void addTriangle(TriangulationDC& T, int l, int r, int cand,const std::vector<Point>& V) {
  addGraphEdge(T, l, r);
  addGraphEdge(T, r, cand);
  addGraphEdge(T, cand, l);
}

TriangulationDC mergeTriangulations(TriangulationDC L, TriangulationDC R, const std::vector<Point>& V) {
  TriangulationDC T;
  T.graph.resize(V.size());

  //combine edge graphs
  for (int i = 0; i <V.size(); i++) {
    T.graph[i].insert(L.graph[i].begin(), L.graph[i].end());
    T.graph[i].insert(R.graph[i].begin(), R.graph[i].end());
  }

  //find points for base LR edge
  std::pair<int, int> baseLR = findBaseLR(L, R, V);

  int l = baseLR.first; //left endpoint of base LR edge
  int r = baseLR.second; //right endpoint of base LR edge

  while (true) {
    //find left and right candidicates for the baseLR (l,r)
    int lCand = findLeftCandidate(T, l, r, V);
    int rCand = findRightCandidate(T, l, r, V);

    //print debugging statement
    // std::cout << "base=(" << l << "," << r << ") " << "lCand=" << lCand << " " << "rCand=" << rCand << "\n";

    if (lCand == -1 && rCand == -1) {
      break; //no candidates found, merge step is done
    } else if (lCand != -1 && rCand == -1) { //only left candidate selected
      addTriangle(T,l,r,lCand, V); //add baseLR along with the cand
      l = lCand; //new base LR edge is lcand with right point from prev base LR
    } else if (lCand == -1 && rCand != -1) { //only right candidate selected
      addTriangle(T,l,r,rCand, V); //add baseLR along with the cand
      r = rCand; //new base LR edge is rCand with left point from prev base LR
    } else { //both candidates picked
      if (!inCircleDC(V,l,r,lCand,rCand)) { //if right candidate is outside the circumcircle of baseLR (l,r) with left cand, then pick left cand to make new edge
        addTriangle(T,l,r,lCand,V);
        l = lCand; //new base LR edge is lCand with right point from prev base LR
      } else { //i.e. if left candidate is outside the circumcircle of baseLR (l,r) with right cand, then pick right cand to make new edge
        addTriangle(T,l,r,rCand,V);
        r = rCand; //new base LR edge is rCand with left point from prev base LR
      }
    }
  }
  return T;
}

//divide and conquer recursive algo
TriangulationDC divideAndConquer(const std::vector<Point>& V, int lo, int hi) {
  //number of points in this subset
  int n = hi - lo; //range is [lo, hi) (lo inclusive, hi exclusive)

  if (n <= 3) { //i.e. n = 2 (make line/edge) or n = 3 (make triangle)
    return triangulateInit(V, lo, hi);
  }

  //else: split again
  int mid = lo + n / 2;

  TriangulationDC L = divideAndConquer(V, lo, mid);
  TriangulationDC R = divideAndConquer(V, mid, hi);

  return mergeTriangulations(L, R, V);
}

//divide and conquer parallel recursive algo
TriangulationDC divideAndConquerPar(const std::vector<Point>& V, int lo, int hi) {
  //number of points in this subset
  int n = hi - lo; //range is [lo, hi) (lo inclusive, hi exclusive)

  if (n <= 3) { //i.e. n = 2 (make line/edge) or n = 3 (make triangle)
    return triangulateInit(V, lo, hi);
  }

  //else: split again
  int mid = lo + n / 2;

  TriangulationDC L;
  TriangulationDC R;

  #pragma omp task shared(L) if (n > 1024)
  {
    L = divideAndConquerPar(V, lo, mid);
  }

  #pragma omp task shared(R) if (n > 1024)
  {
    R = divideAndConquerPar(V, mid, hi);
  }


  #pragma omp taskwait

  return mergeTriangulations(L, R, V);
}


//**following functins are for divide and conquer validation:**

//returns angle of vector from point a to b
//need to sort point's neighbors in by angle, to traverse faces in planar graph
double getAngle(const std::vector<Point>& V, int a, int b) {
  return atan2(V[b].y - V[a].y, V[b].x - V[a].x);
}

//since edge graph only stores connectivity, need to reconstruct faces (regions bound by edges) from it
//outputs a vector of face; each face is list of vertex indices in traversal order
//faces with size three are triangles in the triangulation
std::vector<std::vector<int>> findFaces(const TriangulationDC& T,const std::vector<Point>& V) {
  int n = T.graph.size();

  //store sorted adjacency by angle
  std::vector<std::vector<int>> nbrs(n);

  for (int u = 0; u < n; u++) {
    //copy adjacency list from graph
    for (int v : T.graph[u]) {
      nbrs[u].push_back(v);
    }

    //sort neighbors in counterclockwise order around u
    std::sort(nbrs[u].begin(), nbrs[u].end(), [&](int a, int b) 
      { return getAngle(V, u, a) < getAngle(V, u, b);});
  }

  std::set<std::pair<int, int>> visited; //track (directed) edges already used in face walk
  std::vector<std::vector<int>> faces; //list of faces (face: list of vertex ind)

  //start a face walk from each directed edge (start, next)
  for (int start = 0; start < n; start++) {
    for (int next : nbrs[start]) {
      if (visited.count({start, next})) continue; //skip if already used edge

      std::vector<int> face;

      int u = start;
      int v = next;

      //walk around face
      while (true) {
        visited.insert({u, v});
        face.push_back(u); //add current vertex to face

        auto& nbrsV = nbrs[v];

        //position of prev vertex,u, in v's neighbor list
        auto prevPos = std::find(nbrsV.begin(), nbrsV.end(), u);
        if (prevPos == nbrsV.end()) break;

        int i = prevPos - nbrsV.begin();

        //choose previous neighbor in counterclockwise order
        //basically turns left and keeps face on left side of u,v
        int w = nbrsV[(i - 1 + nbrsV.size()) % nbrsV.size()];

        //move forward along face
        u = v;
        v = w;

        if (u == start && v == next) { //returned to starting edge, so face complete
          break;
        }
      }

      //only keep valid faces
      if (face.size() >= 3) {
        faces.push_back(face);
      }
    }
  }
  return faces;
}

//convertes faces from face walk to actual triangles
//only keeps the triangular faces
std::vector<TriangleDC> getTriangles(const TriangulationDC& T,const std::vector<Point>& V) {
  std::vector<std::vector<int>> faces = findFaces(T, V); //all faces, including outer one
  std::vector<TriangleDC> ts;

  for (auto& f : faces) {
    //only interior triangle faces correspond to actual triangles in the triangulation
    if (f.size() != 3) continue; //if face not size 3, prob outer boundary face

    int a = f[0];
    int b = f[1];
    int c = f[2];

    double o = orient2d(V, a, b, c);

    //face walk got faces in both orientations: one is bounded interior faces, other is outer face (opposite orient)
    //with the previous neighbor rule in findFaces,interior faces should be counterclockwise (so o >0)
    if (o <= 0) continue; //not counterclockwise interior face then

    //store triangles clockwise
    ts.emplace_back(a, c, b); //(a,b,c) counterclockwise (by o above)
  }
  return ts;
}

//check if T satisfies Delaunay:
//a. for each triangle, no point is inside the circumcircle of any triangle
//b. every point is in at least one triangle
//returns true if a and b hold
bool isDelaunayDC(const TriangulationDC& T, const std::vector<Point>& V) {
  std::vector<TriangleDC> ts = getTriangles(T, V); //reconstruct triangles from graph w/ facewalk
  std::vector<int> found(V.size(), 0); //counts # times each point found as part of a triangle

  if (ts.empty()) { //some problem w triangulation...
      std::cout << "No triangle faces found\n";
      return false;
  }

  //check delaunay condition for each triangle
  for (const TriangleDC& t : ts) {
    int a = t.v[0];
    int b = t.v[1];
    int c = t.v[2];

    //count each vertex i.e. is vertex of a triangle
    found[a]++;
    found[b]++;
    found[c]++;

    //check triangle not degenerate
    if (orientation(V[a],V[b],V[c]) == 0.0f) {
      printf("Degenerate triangle! %d %d %d\n", a, b, c);
      return false;
    }

    //check if point p is inside circumcircle of t
    for (int p = 0; p < V.size(); p++) {
      if (p == a || p == b || p == c) continue;

      if (inCircleDC(V, a, b, c, p)) {
        printf("Point %d inside of circumcircle of triangle %d %d %d!\n", p, a, b,c);
        return false;
      }
    }
  }

  //make sure each point is a vertex of at least one triangle
  for (int p = 0; p < V.size(); p++) {
    if (found[p] == 0) {
      printf("Point %d is not used in a triangle\n", p);
      return false;
    }
  }
  return true;
}

//END DIVIDE AND CONQUER

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

  //check if required options are provided
  //* I for incremental
  //* D for divide-and-conquer
  //* S for sequential (incremental)
  if ((input_filename.empty()) || num_threads <= 0 || SA_iters <= 0 ||
      (parallel_mode != 'I' && parallel_mode != 'D' && parallel_mode != 'S' && parallel_mode != 'P') || batch_size <= 0)
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

  Mesh M;
  M.n = n;

  TriangulationDC res;

  std::vector<Point> V(n);

  for (int i = 0; i < n; i++) {
    fin >> V[i].x >> V[i].y;
  }

  // std::cout << "Read " << n << " points:\n";
  // for (int i = 0; i < n; i++)
  // {
  //   std::cout << i << ": (" << V[i].x << ", " << V[i].y << ")\n";
  // }

  const double init_time =
    std::chrono::duration_cast<std::chrono::duration<double>>(
      std::chrono::steady_clock::now() - init_start).count();
      std::cout << "Initialization time (sec): " << std::fixed
      << std::setprecision(10) << init_time << '\n';

  //* Timing Code Start
  const auto compute_start = std::chrono::steady_clock::now();

  //* sequential
  if (parallel_mode == 'S'){
    Triangle tb;

    //* get bounding box
    float minx = V[0].x, maxx = V[0].x;
    float miny = V[0].y, maxy = V[0].y;
    
    for (int i = 0; i < n; i++){
      minx = std::min(V[i].x, minx);
      maxx = std::max(V[i].x, maxx);
      miny = std::min(V[i].y, miny);
      maxy = std::max(V[i].y, maxy);
    }
    float d = std::max(maxx-minx, maxy-miny);
    if (d == 0) d = 1.0;

    //* super triangle calculation: https://www.gorillasun.de/blog/bowyer-watson-algorithm-for-delaunay-triangulation/#the-super-triangle
    V.push_back({minx-10*d, miny-10*d});
    V.push_back({maxx+10*d, miny-10*d});
    V.push_back({(minx+maxx)/2, maxy+10*d});

    tb.x = V.size()-3;
    tb.y = V.size()-2;  tb.z = V.size()-1;

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
        M.triangles[t].nbr_xy = -1;
        M.triangles[t].nbr_yz = -1;
        M.triangles[t].nbr_zx = -1;
      }

    }
  //* parallel incremental
  } else if (parallel_mode == 'I'){
    printf("Starting incremental\n");
  } else if (parallel_mode == 'D'){
    //sequential divide and conquer
    printf("Starting seq divide-and-conquer \n");
    //algorithm based off of Guibas and Stolfi; used this as a reference: http://www.geom.uiuc.edu/~samuelp/del_project.html#algorithms

    //order points by x coordinates (use y coords to tie break)
    std::sort(V.begin(), V.end(), [](const Point& a, const Point& b) {
      if (a.x == b.x) {
        return (a.y < b.y);
      }
      else {
        return (a.x < b.x);
      } });

    // std::cout << "Read " << n << " points:\n";
    // for (int i = 0; i < n; i++) {
    //   std::cout << i << ": (" << V[i].x << ", " << V[i].y << ")\n";
    // }

    //important to do this because need to run the python visualization with the sorted points
    std::string name = input_filename.substr(input_filename.find_last_of("/\\") + 1);
    std::ofstream sortedPts("inputs/sorted_" + name);
    sortedPts << V.size() << "\n";
    for (const Point& p : V) {
      sortedPts << p.x << " " << p.y << "\n";
    }
    sortedPts.close();

    res = divideAndConquer(V, 0, V.size()); 

  } else if (parallel_mode == 'P') {
    //parallel divide and conquer
    printf("Starting par divide-and-conquer \n");

    //order points by x coordinates (use y coords to tie break)
    std::sort(V.begin(), V.end(), [](const Point& a, const Point& b) {
      if (a.x == b.x) {
        return (a.y < b.y);
      }
      else {
        return (a.x < b.x);
      } });

    // std::cout << "Read " << n << " points:\n";
    // for (int i = 0; i < n; i++) {
    //   std::cout << i << ": (" << V[i].x << ", " << V[i].y << ")\n";
    // }

    //important to do this because need to run the python visualization with the sorted points
    std::string name = input_filename.substr(input_filename.find_last_of("/\\") + 1);
    std::ofstream sortedPts("inputs/sorted_" + name);
    sortedPts << V.size() << "\n";
    for (const Point& p : V) {
      sortedPts << p.x << " " << p.y << "\n";
    }
    sortedPts.close();

    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
      #pragma omp single
      {
        res = divideAndConquerPar(V, 0, V.size()); 
      }
    }

  } else {
    assert(false);
  }

  //* Timing Code End

  const double compute_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - compute_start)
          .count();
  std::cout << "Computation time (sec): " << compute_time << "\n";

  bool valid = false;

  if (parallel_mode == 'S') {
    //FOR INCREMENTAL
    valid = isDelaunay(M, V);

  } else if (parallel_mode == 'D' || parallel_mode == 'P') {
    //FOR DIVIDE AND CONQUER
    valid = isDelaunayDC(res,V);
  }

  if (valid){
    std::cout << "Valid Delaunay Triangulation :)\n";
  } else {
    std::cout << "Invalid Delaunay Triangulation :( \n";
  }

  std::string name = input_filename.substr(input_filename.find_last_of("/\\") + 1);
  std::ofstream out("outputs/" + name);
  if (!out) {
      std::cerr << "Failed to open outputs" + name + "\n";
      exit(EXIT_FAILURE);
  }

  if (parallel_mode == 'S') {
    //FOR INCREMENTAL
    for (auto &t : M.triangles) {
      if (!t.active) continue;
      if (t.x >= n || t.y >= n || t.z >= n) continue;
      //* ignore any triangles containing points that are part of the initial 
      //* super triangle
      out << t.x << " " << t.y << " " << t.z << "\n";
    }
  } else if (parallel_mode == 'D' || parallel_mode == 'P') {
    //FOR DIVIDE AND CONQUER
    auto tris = getTriangles(res, V);

    for (const TriangleDC& t : tris) {
      out << t.v[0] << " " << t.v[1] << " " << t.v[2] << "\n";
    }
  }

  out.close();

}