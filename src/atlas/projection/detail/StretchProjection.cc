/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <functional>
#include <sstream>

#include "eckit/config/Parametrisation.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/StretchProjection.h"
#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/grid.h"
#include "atlas/grid/Grid.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "eckit/testing/Test.h"
#include <iomanip> ///< for setprecision


/**
* Projection for LAM stretching
*
* The origin of the xy-system is at (lon,lat) = (0,0)
*
*/

namespace atlas {
namespace projection {
namespace detail {


// specification parameters
template <typename Rotation>
typename StretchLAM<Rotation>::Spec StretchLAM<Rotation>::spec() const {

    Spec proj_st;
    proj_st.set( "delta_low", delta_low_ ); //resolution of the host model
    proj_st.set( "delta_high", delta_high_ ); //resolution of the regional model (regular grid)
    proj_st.set( "var_ratio", var_ratio_ ); //power used for the stretching
    proj_st.set( "x_reg_start", x_reg_start_ ); //xstart of the internal regional grid
    proj_st.set( "x_reg_end", x_reg_end_ ); //xend of the internal regional grid
    proj_st.set( "y_reg_start", y_reg_start_); //ystart of the internal regional grid
    proj_st.set( "y_reg_end", y_reg_end_ ); //yend of the internal regional grid
    proj_st.set( "startx", startx_ ); //domain startx
    proj_st.set( "endx", endx_ ); //domain endx
    proj_st.set( "starty", starty_ ); //domain starty
    proj_st.set( "endy", endy_ ); //domain endy

    return proj_st;
}

// constructors
template <typename Rotation>
StretchLAM<Rotation>::StretchLAM(const eckit::Parametrisation& proj_st) :
    ProjectionImpl(), rotation_( proj_st ) {


    proj_st.get( "delta_low", delta_low_ = 0.); //resolution of the host model
    proj_st.get( "delta_hi", delta_high_ = 0. ); //resolution of the regional model (regular grid)
    proj_st.get( "var_ratio", var_ratio_ = 0. ); //power used for the stretching
    proj_st.get( "x_reg_start", x_reg_start_ = 0. ); //xstart of the internal regional grid
    proj_st.get( "y_reg_start", y_reg_start_ = 0.); //ystart of the internal regional grid
    proj_st.get( "x_reg_end", x_reg_end_ = 0. ); //xend of the regular part of stretched internal grid
    proj_st.get( "y_reg_end", y_reg_end_ = 0.); //yend of the regular part of stretched internal grid
    proj_st.get( "startx", startx_ = 0. ); //domain startx
    proj_st.get( "endx", endx_ = 0. ); //domain endx
    proj_st.get( "starty", starty_ = 0. ); //domain starty
    proj_st.get( "endy", endy_ = 0. ); //domain endy
}

template<typename Rotation>
void StretchLAM<Rotation>::checkvalue(const double & epsilon,
                                        const double & value_check) const {
    std::string err_message;
    std::string str = std::to_string(value_check);
    err_message = "USER defined limits not in the middle of the area " + str;
    if (value_check > epsilon || value_check < (-1 * epsilon))
          {
             throw eckit::BadValue(err_message, Here());
          }
}
// General stretch from a point in regular grid to the stretched grid
template<typename Rotation>
double StretchLAM<Rotation>::general_stretch(double& lamphi, const bool& L_long,
                                             int& n_int, const int& n_stretched_) const {

       double high_size;
       double lamphi_start;
       double lamphi_end;
       double point = lamphi;

       if ((n_int > 0 ) && (var_ratio_ > 1)){

          double delta_dist;
          double delta_add;

          // input number of points, output intervals
          n_int -= 1;
          high_size = (n_int - n_stretched_ ) * delta_high_;

          // number of new internal regular grid in integer
          int ints_high = (high_size + 0.0001) / delta_high_;

          // number of variable grid
          int var_ints = (n_int - ints_high)/2.;
          // compute ratio
          double var_ints_f = (n_int - ints_high)/2.;
          double logr = log(var_ratio_);
          double log_ratio = (var_ints_f - 0.5) * logr;
          double new_ratio = exp(log_ratio / var_ints);

          /** rim To implement
          * int ints_low = rim_width;
          * total intervals not used here but maybe for later
          * int intervals = ints_high + 2* (var_ints + ints_low);
          */

          // INTERNAL REGULAR GRID the point is mapped to the same point
          if (L_long ) {
              lamphi_start = x_reg_start_;
          } else {
              lamphi_start = y_reg_start_;
          }
          lamphi_end = lamphi_start + high_size;

          if ((point >= lamphi_start) && (point <= lamphi_end)){
                   lamphi = point;
          }
          /* VERY IMPORTANT
          * Stat variable res 'STRETCHED'
          * distance from the regular grid and stretch internal grid delta_high
          */

          if (point < lamphi_start) {
             delta_dist = abs(lamphi_start - point);
          }else if (point > lamphi_end) {
             delta_dist = abs(lamphi_end - point);
          }

          if ((point < lamphi_start) ||  (point > lamphi_end)){
             // number of high resolution points intervals
             // delta_dist is the distance of the actual point to the internal regular grid
             int n_high = delta_dist / delta_high_;
             // part remaining, use modulo
             double p_rem = std:: fmod(delta_dist , delta_high_);
             // computation of the new stretched delta integer part
             double delta = delta_high_;
             // initialization if is not using the cycle (first interval)
             double delta_last = delta;

             // using difference in delta
             for (int i = 0; i < n_high; i+=1){
                delta_last = delta * new_ratio;
                delta_add = delta_last - delta_high_;
                delta = delta_last;
                // recomputation of point for every interval
                if (point > lamphi_start){
                   point = point + delta_add;
                }else {
                   point = point - delta_add;
                }
             }

             // last part, delta last from the cycle before
             double delta_r = p_rem * pow(new_ratio, (n_high+1));
             double delta_addr = delta_r - p_rem;
             if (point > lamphi_start){
                point += delta_addr;
             } else {
                point -= delta_addr;
             }

             lamphi = point;
             lamphi = (lamphi >= 360) ? lamphi -360.0 : lamphi ;
          }
          //Section 5: Reset lat/lons

          /*
          * If longitude domain crosses the meridian (0 degrees) then
          * add 360 degrees to points East of the meridian
          */
          if (L_long){
              lamphi = (lamphi < 180) ? lamphi + 360.0 : lamphi ;
          }

       }
       //  var_ratio = 1.0   regular grid
       else {
           // Section 6: regular grids
            lamphi = point;
            if (L_long){
                lamphi = (lamphi < 180) ? lamphi + 360.0 : lamphi ;
            }
        } //var ratio and number of points


    return lamphi;

}

//xy unstretched, only unrotation
template <typename Rotation>
void StretchLAM<Rotation>::lonlat2xy( double crd[] ) const {

    //unrotate
    rotation_.rotate( crd );

    // PUT the unstretch, I don't have it nows

}

//From unstretched to stretched
template <typename Rotation>
void StretchLAM<Rotation>::xy2lonlat( double crd[] ) const {

    /** PUT here the stretch for a point that come input
    * give number of power as input,work it out using as from the start of regular grid.
    * atlas::PointXY unstretchedXY = crd;
    * atlas::PointXY stretchedXY = stretch_LAM_gen(unstretchedXY);
    * stretch_LAM_gen(crd[]);
    */

    double lam_hires_size;
    double phi_hires_size;
    double lambda_start;
    double phi_start;
    double add_xf_;
    double add_yf_;
    double check_x;
    double check_y;
    int rim_size;
    double epsilon = 0.01;
    int nx_, ny_, n_stretched_;


    double deltax_all = (endx_ - startx_);
    double deltay_all = (endy_ - starty_);

    if (var_ratio_ == 1) {
           lam_hires_size = deltax_all;
           phi_hires_size = deltay_all;
           lambda_start = x_reg_start_ ;
           phi_start = y_reg_start_ ;
         } else {
              lam_hires_size = x_reg_end_ - x_reg_start_;
              phi_hires_size = y_reg_end_ - y_reg_start_;

              /** check for start of the new regular LAM grid x and y
              * in the middle of the previous regular grid
              */

              add_xf_ = (deltax_all - lam_hires_size) / 2.;
              add_yf_ = (deltay_all - phi_hires_size) / 2.;
              //+1 otherwise I have the intervals not the points
              int add_x = ((deltax_all - lam_hires_size) / delta_high_) + 1;
              int add_y = ((deltay_all - phi_hires_size) / delta_high_) + 1;
              lambda_start = x_reg_start_;
              phi_start = y_reg_start_;
              check_x = startx_ + add_xf_ - lambda_start;
              check_y = starty_ + add_yf_ - phi_start;
              double check_st;

              check_st = add_x - add_y;
              checkvalue(epsilon, check_x);
              checkvalue(epsilon, check_y);
              checkvalue(epsilon, check_st);

              n_stretched_ = add_x;
              nx_ = ((deltax_all + 0.0001) / delta_high_) + 1;
              ny_ = ((deltay_all + 0.0001) / delta_high_ ) + 1;
         }

    rim_size = 0;
    crd[0] = general_stretch(crd[0], true, nx_, n_stretched_);
    crd[1] = general_stretch(crd[1], false, ny_, n_stretched_);

    // rotate
    rotation_.unrotate( crd );
}

template <typename Rotation>
ProjectionImpl::Jacobian StretchLAM<Rotation>::jacobian( const PointLonLat& ) const {
    throw_NotImplemented( "StretchProjectionT::jacobian", Here() );
}

template <typename Rotation>
void StretchLAM<Rotation>::hash( eckit::Hash& hsh ) const {
    hsh.add( static_type() );
    rotation_.hash( hsh );
    //hsh.add( radius_ );
}

template class StretchLAM<NotRotated>;
template class StretchLAM<Rotated>;

namespace {
static ProjectionBuilder<StretchProjection> register_1( StretchProjection::static_type() );
static ProjectionBuilder<RotatedStretchProjection> register_2( RotatedStretchProjection::static_type() );
}  // namespace

}  // namespace detail
}  // namespace projection
}  // namespace atlas
