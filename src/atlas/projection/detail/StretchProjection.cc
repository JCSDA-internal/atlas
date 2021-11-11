/**
* (C) Crown copyright 2021, Met Office
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
    proj_st.set( "delta_low", delta_low_ ); ///< resolution of the external regular grid (rim)
    proj_st.set( "delta_high", delta_high_ ); ///< resolution of the regional model (regular grid)
    proj_st.set( "var_ratio", var_ratio_ ); ///< power used for the stretching
    proj_st.set( "x_reg_start", x_reg_start_ ); ///< xstart of the internal regional grid
    proj_st.set( "x_reg_end", x_reg_end_ ); ///< xend of the internal regional grid
    proj_st.set( "y_reg_start", y_reg_start_); ///< ystart of the internal regional grid
    proj_st.set( "y_reg_end", y_reg_end_ ); ///< yend of the internal regional grid
    proj_st.set( "startx", startx_ ); ///< original domain startx
    proj_st.set( "endx", endx_ ); ///< original domain endx
    proj_st.set( "starty", starty_ ); ///< original domain starty
    proj_st.set( "endy", endy_ ); ///< original domain endy
    proj_st.set( "rim_widthx", rim_widthx_ ); ///< xsize of the rim
    proj_st.set( "rim_widthy", rim_widthy_ ); ///< ysize of the rim

    return proj_st;
}

// constructors
template <typename Rotation>
StretchLAM<Rotation>::StretchLAM(const eckit::Parametrisation& proj_st) :
    ProjectionImpl(), rotation_( proj_st ) {


    proj_st.get( "delta_low", delta_low_ = 0.); ///< resolution of the external regular grid (rim)
    proj_st.get( "delta_hi", delta_high_ = 0. ); ///< resolution of the regional model (regular grid)
    proj_st.get( "var_ratio", var_ratio_ = 0. ); ///< power used for the stretching
    proj_st.get( "x_reg_start", x_reg_start_ = 0. ); ///< xstart of the internal regional grid
    proj_st.get( "y_reg_start", y_reg_start_ = 0.); ///< ystart of the internal regional grid
    proj_st.get( "x_reg_end", x_reg_end_ = 0. ); ///< xend of the regular part of stretched internal grid
    proj_st.get( "y_reg_end", y_reg_end_ = 0.); ///< yend of the regular part of stretched internal grid
    proj_st.get( "startx", startx_ = 0. ); ///< original domain startx
    proj_st.get( "endx", endx_ = 0. ); ///< original domain endx
    proj_st.get( "starty", starty_ = 0. ); ///< original domain starty
    proj_st.get( "endy", endy_ = 0. ); ///< original domain endy
    proj_st.get( "rim_widthx", rim_widthx_ ); ///< xsize of the rim
    proj_st.get( "rim_widthy", rim_widthy_ ); ///< ysize of the rim
}

template<typename Rotation>
void StretchLAM<Rotation>::checkvalue(const double & epsilon,
                                        const double & value_check) const {
    std::string err_message;
    std::string str = std::to_string(value_check);
    err_message = "USER defined limits not in the middle of the area " + str;
    if (value_check > epsilon || value_check < (-1. * epsilon))
          {
             throw eckit::BadValue(err_message, Here());
          }
}

/**
 * General stretch from a point in regular grid to the
 * a correspective point in a stretched grid
 */

template<typename Rotation>
double StretchLAM<Rotation>::general_stretch(double& lamphi, const bool& L_long,
                                             int& n_int, const int& n_stretched_, const int& n_rim_) const {

       double high_size; ///< number of new internal regular grid in double
       double lamphi_start; ///< start of the regular grid
       double lamphi_end; ///< end of the regular grid
       double point = lamphi; ///< starting point
       double epstest = 0.00000000001; ///< correction used to change from double to integer

       if ((n_int > 0 ) && (var_ratio_ > 1)){

          double delta_dist; ///< distance from point to reg. grid
          double delta_add; ///< additional part in stretch different from internal high resolution
          int n_high; ///< number of points, from point to reg. grid
          int n_high_st; ///< number of stretched points, from point to reg grid
          int n_high_rim; ///< number of rim points, from point to reg grid
          double p_rem; ///< remaining part in stretch if delta_dist not multiple of delta_high_
          double p_rem_low; ///< remaining part in rim if delta_dist not multiple of delta_high_

          n_int -= 1; ///< input number of points, output intervals
          high_size = (n_int - n_stretched_ - n_rim_) * delta_high_;
          int ints_high = (high_size + epstest) / delta_high_;

          ///< number of variable grid points
          int var_ints = (n_int + epstest - n_rim_ - ints_high)/2.;
          ///< compute ratio
          double var_ints_f = (n_int - n_rim_- ints_high)/2.;
          double logr = log(var_ratio_);
          double log_ratio = (var_ints_f - 0.5) * logr;
          double new_ratio = exp(log_ratio / var_ints);

          /**
           *  SECTION 1
           *  INTERNAL REGULAR GRID the point is mapped to the same point
           */

          if (L_long ) {
              lamphi_start = x_reg_start_;
          } else {
              lamphi_start = y_reg_start_;
          }
          lamphi_end = lamphi_start + high_size;

          if ((point >= lamphi_start) || (point <= lamphi_end)){
                   lamphi = point;
          }

          /** SECTION 2
          * Start variable res 'STRETCHED'
          * The point is mapped in a stretched grid related to
          * distance from the regular grid and stretch internal grid: delta_dist
          */

          if (point < lamphi_start) {
             delta_dist =  lamphi_start - point ;
          }else if (point > lamphi_end) {
             delta_dist = point - lamphi_end;
          }

          if ((point < lamphi_start) ||  (point > lamphi_end)){
             /**
              * number of high resolution points intervals, that are not
              * in the internal regular grid on one side,
              * this work just for the part of the stretch not the rim
              */

              // always the lowest integer
             n_high = (delta_dist + epstest) / delta_high_;

             ///< only for the stretched part take out the rim part
             if (n_high > n_stretched_/2. ) {
                 n_high_st = (n_stretched_/2.);
                 n_high_rim = n_high - n_high_st;
                 p_rem = 0;
                 p_rem_low = std:: fmod((delta_dist + epstest) , delta_high_);
             } else {
                 n_high_st = n_high;
                 n_high_rim = 0;
                 // part remaining, use modulo
                 p_rem = std:: fmod((delta_dist + epstest ) , delta_high_);
                 p_rem_low = 0.;
             }

             ///< computation of the new stretched delta integer part
             double delta = delta_high_;
             ///< initialization if is not using the cycle (first interval)
             double delta_last = delta;
             double deltacheck = 0;
             /**
             * using difference in delta for stretch
             * The point stretched is not in the regular grid and took out the points for the rim
             */
             for (int i = 0; i < n_high_st ; i+=1){
                     delta_last = delta * new_ratio;
                     delta_add = delta_last - delta_high_;
                     delta = delta_last;
                     deltacheck += delta_add   ;
                  }
              ///< recomputation of point for every interval
              if (point > lamphi_start){
                  /**
                  * after the end of the internal high resolution grid
                  * no delta_add in the last points as they are rim
                  */
                 point = point + deltacheck;
              }else {
                 /**
                 * before the begin of the internal high resolution grid
                 * no delta_add in the first points as they are rim
                 */
                 point = point - deltacheck;
             }

             // SECTION 3 last part of stretch adding the remaing non integer part with the same ratio as in the stretching
             double delta_r = p_rem * pow(new_ratio, (n_high_st + 1));
             double delta_addr = delta_r - p_rem;

             if (point > lamphi_start){
                point = point + delta_addr;
             } else {
                point = point - delta_addr;
             }
             // SECTION 4 rim area
             if (n_high > n_stretched_/2. ) {
                 double delta_l_h_ = 0;
                 for (int i = 0; i < n_high_rim; i+=1){
                    delta_l_h_ = delta_l_h_ + (delta_low_ - delta_high_ );
                 }
                 if (point > lamphi_start){
                   point = point + delta_l_h_ + p_rem_low * (delta_low_ - delta_high_ );
                 }else {
                   point = point - delta_l_h_ - p_rem_low * (delta_low_ - delta_high_ );
                 }
             }

             lamphi = point;
             lamphi = (lamphi >= 360) ? lamphi -360.0 : lamphi ;
          }


          /**
          * SECTION 5: Reset lat/lons
          * If longitude domain crosses the meridian (0 degrees) then
          * add 360 degrees to points East of the meridian
          */
          if (L_long){
              lamphi = (lamphi < 180) ? lamphi + 360.0 : lamphi ;
          }

       } 
       else {
           /*
            * SECTION 6: regular grids
            * var_ratio = 1.0
            */
            lamphi = point;
            if (L_long){
                lamphi = (lamphi < 180) ? lamphi + 360.0 : lamphi ;
            }
        }
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

    double lam_hires_size; ///< size regular grid x
    double phi_hires_size; ///< size regular grid y
    double lambda_start; ///< start grid x
    double phi_start; ///< start grid y
    double add_xf_; ///< distance end of the grid and internal regular grid
    double add_yf_; ///< distance end of the grid and internal regular grid
    double check_x; ///< check middle of the previous regular grid
    double check_y; ///< check middle of the previous regular grid
    double epsilon = 0.01; ///< use in check if the same value
    double epstest = 0.00000000001; ///< correction used to change from double to integer
    int nx_, ny_;
    int n_stretchedx_, n_stretchedy_, n_x_rim_, n_y_rim_;

    ///< original domain size includes the points for the rim
    double deltax_all = (endx_ - startx_);
    double deltay_all = (endy_ - starty_);

    if (var_ratio_ == 1) {
           lam_hires_size = deltax_all;
           phi_hires_size = deltay_all;
           lambda_start = x_reg_start_;
           phi_start = y_reg_start_;

         } else {
              lam_hires_size = x_reg_end_ - x_reg_start_;
              phi_hires_size = y_reg_end_ - y_reg_start_;

              /** check for start of the new regular LAM grid x and y
              * in the middle of the previous regular grid
              */

              add_xf_ = (deltax_all + epstest - lam_hires_size) / 2.;
              add_yf_ = (deltay_all + epstest - phi_hires_size) / 2.;
              /**
               *  Compute the number of points for different part of the grid
               *  internal regular grid high resolution
               *  stretched grid
               *  external regular grid low resolution
               *  +1 otherwise I have the intervals not the points
              */
              n_x_rim_ = rim_widthx_/delta_low_ ;
              n_y_rim_ = rim_widthy_/delta_low_ ;
              n_stretchedx_ = ((deltax_all + epstest - lam_hires_size) / delta_high_) - n_x_rim_;
              n_stretchedy_ = ((deltay_all + epstest - phi_hires_size) / delta_high_) - n_y_rim_ ;
              lambda_start = x_reg_start_;
              phi_start = y_reg_start_;
              check_x = startx_ + add_xf_ - lambda_start;
              check_y = starty_ + add_yf_ - phi_start;
              double check_st;

              check_st = n_stretchedx_ - n_stretchedy_;
              checkvalue(epsilon, check_x);
              checkvalue(epsilon, check_y);
              checkvalue(epsilon, check_st);

              nx_ = ((deltax_all + epstest) / delta_high_) + 1;
              ny_ = ((deltay_all + epstest) / delta_high_ ) + 1;

         }

    crd[0] = general_stretch(crd[0], true, nx_, n_stretchedx_, n_x_rim_);
    crd[1] = general_stretch(crd[1], false, ny_, n_stretchedy_, n_y_rim_);

    ///< rotate
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
