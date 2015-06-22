// TL1279

#include "atlas/grids/rgg/rgg.h"

namespace atlas {
namespace grids {
namespace rgg {

eckit::ConcreteBuilderT1<Grid,N640> builder_N640 (N640::grid_type_str());

void N640::construct()
{
  int N=640;
  long lon[] = {
    18,
    25,
    32,
    40,
    45,
    50,
    60,
    60,
    72,
    72,
    75,
    81,
    90,
    90,
    96,
    100,
    108,
    120,
    120,
    125,
    135,
    144,
    150,
    160,
    160,
    180,
    180,
    180,
    192,
    192,
    200,
    216,
    216,
    216,
    225,
    240,
    240,
    243,
    250,
    256,
    270,
    270,
    288,
    288,
    288,
    300,
    300,
    320,
    320,
    320,
    360,
    360,
    360,
    360,
    360,
    360,
    375,
    375,
    384,
    384,
    400,
    400,
    400,
    432,
    432,
    432,
    432,
    450,
    450,
    450,
    480,
    480,
    480,
    480,
    480,
    486,
    500,
    500,
    512,
    512,
    540,
    540,
    540,
    540,
    540,
    576,
    576,
    576,
    576,
    576,
    600,
    600,
    600,
    600,
    640,
    640,
    640,
    640,
    640,
    640,
    640,
    648,
    675,
    675,
    675,
    675,
    720,
    720,
    720,
    720,
    720,
    720,
    720,
    720,
    729,
    750,
    750,
    750,
    750,
    768,
    768,
    768,
    800,
    800,
    800,
    800,
    800,
    810,
    810,
    864,
    864,
    864,
    864,
    864,
    864,
    864,
    864,
    900,
    900,
    900,
    900,
    900,
    900,
    960,
    960,
    960,
    960,
    960,
    960,
    960,
    960,
    960,
    960,
    960,
    972,
    972,
    1000,
    1000,
    1000,
    1000,
    1000,
    1024,
    1024,
    1024,
    1024,
    1080,
    1080,
    1080,
    1080,
    1080,
    1080,
    1080,
    1080,
    1080,
    1125,
    1125,
    1125,
    1125,
    1125,
    1125,
    1125,
    1125,
    1152,
    1152,
    1152,
    1152,
    1152,
    1200,
    1200,
    1200,
    1200,
    1200,
    1200,
    1200,
    1200,
    1215,
    1215,
    1215,
    1280,
    1280,
    1280,
    1280,
    1280,
    1280,
    1280,
    1280,
    1280,
    1280,
    1280,
    1280,
    1296,
    1296,
    1350,
    1350,
    1350,
    1350,
    1350,
    1350,
    1350,
    1350,
    1350,
    1350,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1440,
    1458,
    1458,
    1458,
    1458,
    1500,
    1500,
    1500,
    1500,
    1500,
    1500,
    1500,
    1500,
    1536,
    1536,
    1536,
    1536,
    1536,
    1536,
    1536,
    1600,
    1600,
    1600,
    1600,
    1600,
    1600,
    1600,
    1600,
    1600,
    1600,
    1600,
    1600,
    1600,
    1620,
    1620,
    1620,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1728,
    1800,
    1800,
    1800,
    1800,
    1800,
    1800,
    1800,
    1800,
    1800,
    1800,
    1800,
    1800,
    1800,
    1800,
    1800,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1875,
    1920,
    1920,
    1920,
    1920,
    1920,
    1920,
    1920,
    1920,
    1920,
    1920,
    1920,
    1944,
    1944,
    1944,
    1944,
    1944,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2025,
    2025,
    2025,
    2025,
    2025,
    2025,
    2048,
    2048,
    2048,
    2048,
    2048,
    2048,
    2048,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2160,
    2187,
    2187,
    2187,
    2187,
    2187,
    2187,
    2187,
    2187,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2250,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2304,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2400,
    2430,
    2430,
    2430,
    2430,
    2430,
    2430,
    2430,
    2430,
    2430,
    2430,
    2430,
    2430,
    2430,
    2430,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2500,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560,
    2560
  };
  std::vector<double> lats(N);
  gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon,DEG);
}

} // namespace rgg
} // namespace grids
} // namespace atlas
