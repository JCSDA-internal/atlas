// TL31

#include "atlas/grids/rgg/rgg.h"

namespace atlas {
namespace grids {
namespace rgg {

eckit::ConcreteBuilderT1<Grid,N16> builder_N16 (N16::grid_type_str());

void N16::construct()
{
  int N=16;
  long lon[] = {
      20,
      27,
      32,
      40,
      45,
      48,
      60,
      60,
      64,
      64,
      64,
      64,
      64,
      64,
      64,
      64
  };
  std::vector<double> lats(N);
  gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon,DEG);
}

} // namespace rgg
} // namespace grids
} // namespace atlas
