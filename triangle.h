/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

#ifndef __WIREOPT_H__
#define __WIREOPT_H__

#include <cstdint>
#include <omp.h>
#include <vector>

#define MAX_PTS_PER_WIRE 5
#define COST_REPORT_DEPTH 10

/** README(student):

*/
struct validate_wire_t {
  uint8_t num_pts;
  struct {
    uint16_t x;
    uint16_t y;
  } p[MAX_PTS_PER_WIRE];
  validate_wire_t &cleanup(void);
  void print_wire(void) const;
};

/* TODO (student): Define the data structure for wire here.
Feel free to redefine this, notice you also need to change how bend is read in
inside main of wireroute.cpp.

A better readable way of define this is:
int start_x, start_y, end_x, end_y, bend1_x, bend1_y, bend2_x, bend2_y, bend3_x,
bend3_y; but this might not be the most efficient way to define the solution
space for a wire with <= 3 bends.
*/

struct Wire {
  /* Define the data structure for wire here. */
  int start_x, start_y, end_x, end_y;
  int numBends;
  int bendsX[3], bendsY[3];
  validate_wire_t to_validate_format(void) const;
};

struct Triangle {
  int x, y, z;
  vector<int> E; //encroach set
  //neighbors
  //index of the neighbor triangle in the mesh vector M
  //-1 if no neighbor
  int nbr_xy; //neighbor across edge (x,y)
  int nbr_yz;
  int nbr_zx;

  bool active;
};

struct Mesh {
  std::vector<Triangle> triangles;
};

struct Face {
  int a, b;
};

// Definition of the wire checker
struct wr_checker {
  std::vector<Wire> wires;
  std::vector<std::vector<int>> occupancies;
  const int nwires;
  const int dim_x;
  const int dim_y;
  wr_checker(std::vector<Wire> &wires,
             std::vector<std::vector<int>> &occupancies)
      : wires(wires), occupancies(occupancies), nwires(wires.size()),
        dim_x(occupancies[0].size()), dim_y(occupancies.size()) {}
  void validate() const;
};

const char *get_option_string(const char *option_name,
                              const char *default_value);
int get_option_int(const char *option_name, int default_value);
float get_option_float(const char *option_name, float default_value);

#endif
