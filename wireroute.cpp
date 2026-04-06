/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(eleepalo), Name 2(darrison)
 */

#include "wireroute.h"

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

int wireCost(const Wire &wire, const std::vector<std::vector<int>> &occupMatrix)
{
  int cost = 0;
  int x = wire.start_x;
  int y = wire.start_y;
  int endX = wire.end_x;
  int endY = wire.end_y;
  int numBends = wire.numBends;

  // add initial cell
  cost += occupMatrix[y][x] * occupMatrix[y][x];

  for (int b = 0; b <= numBends; b++)
  {
    int targetX, targetY;

    if (b < numBends)
    {
      targetX = wire.bendsX[b];
      targetY = wire.bendsY[b];
    }
    // need to account for segment between bend and endpoint
    else
    {
      targetX = endX;
      targetY = endY;
    }

    // works bc at any time, either y changing or x changing NOT both
    int moveX = 0;
    int moveY = 0;

    if (x < targetX)
      moveX = 1;
    if (x > targetX)
      moveX = -1;
    if (y < targetY)
      moveY = 1;
    if (y > targetY)
      moveY = -1;

    while (x != targetX || y != targetY)
    {
      x += moveX;
      y += moveY;
      cost += occupMatrix[y][x] * occupMatrix[y][x];
    }
  }
  return cost;
}

bool inWire(int x, int y, const Wire &wire)
{

  int prev_x = wire.start_x;
  int prev_y = wire.start_y;

  //* check each segment
  for (int i = 0; i <= wire.numBends; i++)
  {

    int next_x, next_y;

    //* next point
    if (i < wire.numBends)
    {
      next_x = wire.bendsX[i];
      next_y = wire.bendsY[i];
    }
    else
    {
      next_x = wire.end_x;
      next_y = wire.end_y;
    }

    //* check vertical segment
    if (prev_x == next_x)
    {
      if (x == prev_x &&
          y >= std::min(prev_y, next_y) &&
          y <= std::max(prev_y, next_y))
      {
        return true;
      }
    }

    //* check horizontal segment
    if (prev_y == next_y)
    {
      if (y == prev_y &&
          x >= std::min(prev_x, next_x) &&
          x <= std::max(prev_x, next_x))
      {
        return true;
      }
    }

    prev_x = next_x;
    prev_y = next_y;
  }

  return false;
}

void addCells(std::set<std::pair<int, int>> &originalCells, const Wire &wire)
{
  int x = wire.start_x;
  int y = wire.start_y;
  int endX = wire.end_x;
  int endY = wire.end_y;
  int numBends = wire.numBends;

  for (int b = 0; b <= numBends; b++)
  {
    int targetX, targetY;

    if (b < numBends)
    {
      targetX = wire.bendsX[b];
      targetY = wire.bendsY[b];
    }
    // need to account for segment between bend and endpoint
    else
    {
      targetX = endX;
      targetY = endY;
    }

    // works bc at any time, either y changing or x changing NOT both
    int moveX = 0;
    int moveY = 0;

    if (x < targetX)
      moveX = 1;
    if (x > targetX)
      moveX = -1;
    if (y < targetY)
      moveY = 1;
    if (y > targetY)
      moveY = -1;

    while (x != targetX || y != targetY)
    {
      x += moveX;
      y += moveY;
      originalCells.insert({x, y});
    }
  }
}

//* original cells must already be updated using add cells
int wireCostWithout(std::set<std::pair<int, int>> &originalCells, const Wire &wire, const std::vector<std::vector<int>> &occupMatrix)
{

  int cost = 0;
  int x = wire.start_x;
  int y = wire.start_y;
  int endX = wire.end_x;
  int endY = wire.end_y;
  int numBends = wire.numBends;

  for (int b = 0; b <= numBends; b++)
  {
    int targetX, targetY;

    if (b < numBends)
    {
      targetX = wire.bendsX[b];
      targetY = wire.bendsY[b];
    }
    // need to account for segment between bend and endpoint
    else
    {
      targetX = endX;
      targetY = endY;
    }

    // works bc at any time, either y changing or x changing NOT both
    int moveX = 0;
    int moveY = 0;

    if (x < targetX)
      moveX = 1;
    if (x > targetX)
      moveX = -1;
    if (y < targetY)
      moveY = 1;
    if (y > targetY)
      moveY = -1;

    while (x != targetX || y != targetY)
    {
      x += moveX;
      y += moveY;
      int occ = occupMatrix[y][x];

      // if (inWire(x,y,noWire))
      if (originalCells.count({x, y}))
        occ -= 1;

      cost += occ * occ;
    }
  }
  return cost;
}

int totalOccupCost(const std::vector<std::vector<int>> &occupMatrix)
{
  int cost = 0;
  for (int y = 0; y < occupMatrix.size(); y++)
  {
    for (int x = 0; x < occupMatrix[y].size(); x++)
    {
      cost += occupMatrix[y][x] * occupMatrix[y][x];
    }
  }
  return cost;
}

static int wAbs(int a)
{
  if (a < 0)
  {
    return -a;
  }
  else
  {
    return a;
  }
}

static void updateSeg(std::vector<std::vector<int>> &occupMatrix, int xS, int yS, int xE, int yE, int delta, bool skip)
{
  int x = xS;
  int y = yS;

  // skip false: first seg of wire, so include starting cell (xS,yS)
  // skip true: seg begins at bend point, was already added by prev seg
  if (!skip)
  {
    occupMatrix[y][x] += delta;
  }

  while (x < xE)
  {
    x++;
    occupMatrix[y][x] += delta;
  }
  while (x > xE)
  {
    x--;
    occupMatrix[y][x] += delta;
  }
  while (y < yE)
  {
    y++;
    occupMatrix[y][x] += delta;
  }
  while (y > yE)
  {
    y--;
    occupMatrix[y][x] += delta;
  }
}

static void updateWire(std::vector<std::vector<int>> &occupMatrix, const Wire &w, int delta)
{
  int xS = w.start_x;
  int yS = w.start_y;

  for (int i = 0; i < w.numBends; i++)
  {
    int xE = w.bendsX[i];
    int yE = w.bendsY[i];

    updateSeg(occupMatrix, xS, yS, xE, yE, delta, i != 0);
    xS = xE;
    yS = yE;
  }

  updateSeg(occupMatrix, xS, yS, w.end_x, w.end_y, delta, w.numBends != 0);
}

// build kth legal route
static void buildKRoute(const Wire &w, int k, Wire &newRoute)
{
  newRoute = w;
  newRoute.numBends = 0;

  int sX = w.start_x;
  int sY = w.start_y;
  int eX = w.end_x;
  int eY = w.end_y;

  int dx = wAbs(eX - sX);
  int dy = wAbs(eY - sY);

  // collinear, only one legal route
  if (dx == 0 || dy == 0)
  {
    newRoute.numBends = 0;
    return;
  }

  int xMin = sX;
  int xMax = sX;
  int yMin = sY;
  int yMax = sY;
  if (eX < xMin)
  {
    xMin = eX;
  }
  if (eX > xMax)
  {
    xMax = eX;
  }
  if (eY < yMin)
  {
    yMin = eY;
  }
  if (eY > yMax)
  {
    yMax = eY;
  }

  // first dx+dy routes: 1 bend, 2 bends
  // k=0: 1 bend horiz first (eX,sY)
  // k=1: 1 bend vert first (sX,eY)
  if (k == 0)
  {
    newRoute.numBends = 1;
    newRoute.bendsX[0] = eX;
    newRoute.bendsY[0] = sY;
    return;
  }
  if (k == 1)
  {
    newRoute.numBends = 1;
    newRoute.bendsX[0] = sX;
    newRoute.bendsY[0] = eY;
    return;
  }

  k -= 2;

  // next (dx-1): 2 bend horis vert horix
  if (k < dx - 1)
  {
    int xMid = xMin + 1 + k;
    newRoute.numBends = 2;
    newRoute.bendsX[0] = xMid;
    newRoute.bendsY[0] = sY;
    newRoute.bendsX[1] = xMid;
    newRoute.bendsY[1] = eY;
    return;
  }

  k -= (dx - 1);

  // next (dy-1): 2 bend vert horiz vert
  if (k < (dy - 1))
  {
    int yMid = yMin + 1 + k;
    newRoute.numBends = 2;
    newRoute.bendsX[0] = sX;
    newRoute.bendsY[0] = yMid;
    newRoute.bendsX[1] = eX;
    newRoute.bendsY[1] = yMid;
    return;
  }

  k -= (dy - 1);

  // 3 bend routes (2(dx-1)(dy-1))
  // split into 2 options: 1. horiz dir first 2. vert dir first
  // each dir has (dx-1)*(dy-1) routes (block gives size of one dir routes)
  int block = (dx - 1) * (dy - 1);
  int dir = 0;
  if (k >= block)
  { // in second dir (vert first), so sub block to reindex k
    dir = 1;
    k -= block;
  }

  int xi = k / (dy - 1); // 0, ..., dx-2
  int yi = k % (dy - 1); // 0, ..., dy-2
  int xMid = xMin + 1 + xi;
  int yMid = yMin + 1 + yi;

  newRoute.numBends = 3;

  if (dir == 0)
  { // horiz first (xMid, sY), (xMid, yMid), (eX, yMid)
    newRoute.bendsX[0] = xMid;
    newRoute.bendsY[0] = sY;

    newRoute.bendsX[1] = xMid;
    newRoute.bendsY[1] = yMid;

    newRoute.bendsX[2] = eX;
    newRoute.bendsY[2] = yMid;
  }
  else
  { // vert first (sX,yMid), (xMid, yMid), (xMid, eY)
    newRoute.bendsX[0] = sX;
    newRoute.bendsY[0] = yMid;

    newRoute.bendsX[1] = xMid;
    newRoute.bendsY[1] = yMid;

    newRoute.bendsX[2] = xMid;
    newRoute.bendsY[2] = eY;
  }
}

//* thread-safe update
void updateOccupancyTS(const Wire &wire, int delta)
{
}

int main(int argc, char *argv[])
{
  const auto init_start = std::chrono::steady_clock::now();

  srand(12345); // seed for randomness

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

  int dim_x, dim_y;
  int num_wires;

  /* Read the grid dimension and wire information from file */
  fin >> dim_x >> dim_y >> num_wires;

  std::vector<Wire> wires(num_wires);
  std::vector occupancy(dim_y, std::vector<int>(dim_x));

  //*
  // std::vector locks(dim_y, std::vector<omp_lock_t>(dim_x));

  std::cout << "Question Spec: dim_x=" << dim_x << ", dim_y=" << dim_y
            << ", number of wires=" << num_wires << '\n';

  // TODO (student code start): Read the wire information from file,
  // you may need to change this if you define the wire structure differently.
  for (auto &wire : wires)
  {
    fin >> wire.start_x >> wire.start_y >> wire.end_x >> wire.end_y;
    if (wire.start_x == wire.end_x || wire.start_y == wire.end_y)
    {
      wire.numBends = 0;
    }
    else
    {
      wire.numBends = 1;
      wire.bendsX[0] = wire.end_x;
      wire.bendsY[0] = wire.start_y;
    }
  }
  // above is just initializing init position for wire...want this to be quick and cheap, so just do simple horiz 1 bend

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

  // initialize wires
  // Within wires
  if (parallel_mode == 'W')
  {
    // within wires
    omp_set_num_threads(num_threads);

    for (int i = 0; i < num_wires; i++)
    {
      // build occupancy from init routes
      updateWire(occupancy, wires[i], +1);
    }

    float r = static_cast<float>(rand()) / RAND_MAX;

    for (int iter = 0; iter < SA_iters; iter++)
    {
      for (int i = 0; i < num_wires; i++)
      {
        Wire w = wires[i];

        int dx = wAbs(w.end_x - w.start_x);
        int dy = wAbs(w.end_y - w.start_y);

        // if collinear, do nothing
        if (dx == 0 || dy == 0)
          continue;

        // remove current route from occupancy
        updateWire(occupancy, w, -1);

        int numRoutes = dx + dy + (2 * (dx - 1) * (dy - 1));
        Wire bestWire = wires[i];

        // recompute curr cost after factoring out
        int bestCost = wireCost(wires[i], occupancy);

        // simulated annealing
        bool chooseRand = 0;
        // with prob P, choose rand route
        if (r < SA_prob)
        {
          chooseRand = 1;
        }

        if (chooseRand)
        {
          int k = rand() % numRoutes;

          Wire newRoute;
          buildKRoute(w, k, newRoute);
          wires[i] = newRoute;
        }
        else
        { // explore all routes in parallel, pick min
#pragma omp parallel
          {
            Wire localWire = bestWire;
            int localCost = bestCost;

#pragma omp for schedule(static)
            for (int k = 0; k < numRoutes; k++)
            {
              Wire newRoute;
              buildKRoute(w, k, newRoute);
              int c = wireCost(newRoute, occupancy);
              if (c < localCost)
              {
                localCost = c;
                localWire = newRoute;
              }
            }

#pragma omp critical
            {
              if (localCost < bestCost)
              {
                bestCost = localCost;
                bestWire = localWire;
              }
            }
          }

          wires[i] = bestWire;
        }
        // add (poss changed) route back into matrix
        updateWire(occupancy, wires[i], +1);
      }
    }
  }
  else
  {
    auto t1 = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_wires; i++)
    {
      // build occupancy from init routes
      updateWire(occupancy, wires[i], +1);
    }

    // across wires
    for (int iter = 0; iter < SA_iters; iter++)
    {
      std::vector<Wire> new_routes(num_wires);

      #pragma omp parallel
      {
        //* while (there is still work to do in this timestep)
        #pragma omp for schedule(dynamic, 1)
        for (int b = 0; b < (1 + (num_wires - 1) / batch_size); b++)
        {
          int batch_start = b * batch_size;
          int batch_end = std::min(batch_start + batch_size, num_wires);
          int local_batch_size = batch_end - batch_start;

          //* grab a batch of B wires to be routed;
          for (int i = 0; i < local_batch_size; i++)
          {
            
            //* determine the new route for wire i (within the batch);
            //* Note: remember the route, but do not update the occupancy matrix yet

            std::set<std::pair<int, int>> originalCells;

            Wire original = wires[batch_start + i];
            new_routes[batch_start + i] = original;
            addCells(originalCells, original);

            Wire bestRoute = original;

            int dx = abs(original.end_x - original.start_x);
            int dy = abs(original.end_y - original.start_y);

            if (dx != 0 && dy != 0)
            {
              int numRoutes = dx + dy + 2 * (dx - 1) * (dy - 1);
              int bestCost = wireCostWithout(originalCells, original, occupancy);

              for (int k = 0; k < numRoutes; k++)
              {
                Wire candidate;
                buildKRoute(original, k, candidate);

                int cost = wireCostWithout(originalCells, candidate, occupancy);

                if (cost < bestCost)
                {
                  bestCost = cost;
                  bestRoute = candidate;
                }
              }
            }

            new_routes[batch_start + i] = bestRoute;
          }
        }
      }

      auto t2 = std::chrono::high_resolution_clock::now();

      //* Note: the loop below does not start until after the loop above finishes
      for (int i = 0; i < num_wires; i++)
      {

        //* update the occupancy matrix for wire i (within the batch);
        //* remove old
          updateWire(occupancy, wires[i], -1);

          //* add new
          updateWire(occupancy, new_routes[i], +1);

          wires[i] = new_routes[i];
      }
      auto t3 = std::chrono::high_resolution_clock::now();
      auto routing_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
auto update_time  = std::chrono::duration_cast<std::chrono::duration<double>>(t3 - t2).count();

printf("%f %f\n", routing_time, update_time);
    }
  }

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
