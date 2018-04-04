#include "computational_geometry.H"

#include "dcel_mesh.H"

class dcel_geometry : public computational_geometry {
public:
  dcel_geometry(const dcel_mesh* const a_mesh);
  dcel_geometry();
  ~dcel_geometry();
};

