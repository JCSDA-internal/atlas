/*
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
// General stretch from a point in regular grid to the stretched grid
template<typename Rotation>
double StretchLAM<Rotation>::general_stretch(double& lamphi, const bool& L_long,
                                             int& n_int, const int& n_stretched_, const int& n_rim_) const {

       double high_size;
       double lamphi_start;
       double lamphi_end;
       double point = lamphi;

       if ((n_int > 0 ) && (var_ratio_ > 1)){

          double delta_dist;
          double delta_add;
          int n_high;
          int n_high_st;
          int n_high_rim;
          double p_rem;
          double p_rem_low;
          double n_high_f;
          double delta_n;

          // input number of points, output intervals
          n_int -= 1;
          high_size = (n_int - n_stretched_ - n_rim_) * delta_high_;

          // number of new internal regular grid in integers
          int ints_high = high_size / delta_high_;

          // number of variable grid
          int var_ints = (n_int - n_rim_ - ints_high)/2.;
          // compute ratio
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

          if ((point >= lamphi_start) && (point <= lamphi_end)){
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
             n_high = (delta_dist + 0.00001) / delta_high_;

             //only for the stretched part
             if (n_high > n_stretched_/2. ) {
                 n_high_st = (n_stretched_/2. + 0.00001);
                 n_high_rim = n_high - n_high_st;
                 n_high_f = (delta_dist - n_high_rim * delta_high_ ) / delta_high_ ;
                 p_rem = 0;
                 p_rem_low = std:: fmod((delta_dist+0.00001) , delta_high_);
                 delta_n = (n_high +1) - n_high_f ;
             } else {
                 n_high_st = n_high;
                 n_high_rim = 0;
                 n_high_f = (delta_dist) / delta_high_;
                 // part remaining, use modulo
                 p_rem = std:: fmod((delta_dist + 0.00001) , delta_high_);
                 p_rem_low = 0.;
                 delta_n = (n_high +1) - n_high_f ;
             }


             /*
             //this is for the case p_rem very near to the high resolution dimension
             if (delta_n < 0.00001 && delta_n > 0.) {
                 p_rem = n_high_f - n_high;
                 p_rem = 0;
                 n_high = n_high + 1;
                 std::cout <<  std::setprecision(14) << "n_high_f: " << n_high_f << std::endl;
                 std::cout << "n_high: " << n_high << std::endl;
             }
             */
             /**
              *  check if we are outside the stretch part,thus in the rim part
              *  p_rem remain the same
              */


             std::cout <<  " n high " << (delta_dist ) / delta_high_ << std::endl;
             std::cout << " n high changed " << (delta_dist + 0.00001) / delta_high_ << std::endl;
             std::cout <<  std::setprecision(18) << " delta dist " << delta_dist << std::endl;
             std::cout <<  std::setprecision(18) << " delta high " << delta_high_ << std::endl;
             std::cout << "number of intervals n_int: " << n_int << std::endl;
             std::cout << "n_high: " << n_high << std::endl;
             std::cout << "n_high_st: " << n_high_st << std::endl;
             std::cout << "n_high_rim: " << n_high_rim << std::endl;
             std::cout <<  std::setprecision(18) << "p_rem: " << p_rem << std::endl;
             std::cout << " " << std::endl;

             // computation of the new stretched delta integer part
             double delta = delta_high_;
             // initialization if is not using the cycle (first interval)
             double delta_last = delta;
             double deltacheck = 0;


             /**
             * using difference in delta for stretch
             * The point stretched is not in the regular grid and took out the points for the rim
             */
             std::cout << "point to change: " << point << std::endl;
             std::cout << "lamphi_start: " << lamphi_start << std::endl;
             std::cout << "n_high_st:  " << n_high_st << std::endl;
             for (int i = 0; i < n_high_st ; i+=1){
                     delta_last = delta * new_ratio;
                     // additional part different from internal high resolution
                     delta_add = delta_last - delta_high_;
                     delta = delta_last;
                     deltacheck += delta_add   ;
                  }
              // recomputation of point for every interval
              if (point > lamphi_start){
                 //after the end of the internal high resolution grid
                 // no delta_add in the last points as they are rim
                 point = point + deltacheck;
              }else {
                 //before the begin of the internal high resolution grid
                 //no delta_add in the first points as they are rim
                 std::cout << "point before adjustment: " << point << std::endl;
                 point = point - deltacheck;
             }


            std::cout << "point after stretch integer: " << point << std::endl;

             // SECTION 3 last part of stretch adding the remaing non integer part with the same ratio as in the stretching
             double delta_r = p_rem * pow(new_ratio, (n_high_st + 1));

             /*
             //this is for the case p_rem very near to the high resolution dimension
             if (abs(delta_n) < 0.00001 && abs(delta_n) > 0.) {
                 // in this case we are not using p_rem for delta
                 //part for correction due to precisioin delta_r =0 until here
                 //delta_r = - delta_n * pow(new_ratio, (n_high_st)) + delta_n;
                 //delta_r = delta_n;
                 //delta_r = - delta_n * pow(new_ratio, (n_high_st-1)) + delta_n;
                 //delta_r = - delta_n;
                 //delta_r = 0. ;
                 //delta_r = delta_n * pow(new_ratio, (n_high_st-1)) - delta_n;
                 //delta_r = - delta_n * pow(new_ratio, (n_high_st)) + delta_n;
                 // n_high >   n_high_f
                 delta_r = (n_high - n_high_f ) * pow(new_ratio, (n_high)) ;
                 //p_rem = - delta_n;

                 std::cout << std::setprecision(14) << "p_rem : " << p_rem_epsilon << std::endl;
                 std::cout << std::setprecision(14) << "delta_r +: " << delta_r << std::endl;
             }
             */

             double delta_addr = delta_r - p_rem;

             if (point > lamphi_start){
                point = point + delta_addr;
             } else {
                point = point - delta_addr;
             }
             std::cout << "point after stretch: " << point << std::endl;


             // SECTION 4 rim area
             if (n_high > n_stretched_/2. ) {
                 std::cout << "n_high_rim: " << n_high_rim << std::endl;
                 double delta_l_h_ = 0;
                 for (int i = 0; i < n_high_rim; i+=1){
                    delta_l_h_ = delta_l_h_ + (delta_low_ - delta_high_ );
                 }
                 std::cout << "delta_l_h_: " << delta_l_h_ << "rest: "<< p_rem_low * (delta_low_ - delta_high_ ) << std::endl;
                 if (point > lamphi_start){
                   point = point + delta_l_h_ + p_rem_low * (delta_low_ - delta_high_ );
                 }else {
                   point = point - delta_l_h_ - p_rem_low * (delta_low_ - delta_high_ );
                 }
             }


             std::cout << "point end: " << point << std::endl;
             std::cout << " " << std::endl;

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
       //  var_ratio = 1.0   regular grid
       else {
           // SECTION 6: regular grids
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

              add_xf_ = (deltax_all+ 0.0001 - lam_hires_size) / 2.;
              add_yf_ = (deltay_all+ 0.0001 - phi_hires_size) / 2.;
              std::cout << "addx " << add_xf_ << " addy " << add_yf_  <<  std::endl;
              /**
               *  Compute the number of points for different part of the grid
               *  internal regular grid high resolution
               *  stretched grid
               *  external regular grid low resolution
               *  +1 otherwise I have the intervals not the points
              */
              //n_stretchedx_ = ((deltax_all + 0.0001 - (rim_widthx_ * delta_low_) - lam_hires_size) / delta_high_) ;
              //n_stretchedy_ = ((deltay_all + 0.0001 - (rim_widthy_ * delta_low_ ) - phi_hires_size) / delta_high_) ;

              n_x_rim_ = rim_widthx_/delta_low_ ;
              n_y_rim_ = rim_widthy_/delta_low_ ;
              n_stretchedx_ = ((deltax_all + 0.0001 - lam_hires_size) / delta_high_) - n_x_rim_;
              n_stretchedy_ = ((deltay_all + 0.0001  - phi_hires_size) / delta_high_) - n_y_rim_ ;
              std::cout << "endx " << endx_ << " startx " << startx_  << " deltax_all " << deltax_all <<  std::endl;
              std::cout << "endy " << endy_ << " starty " << starty_  << " deltay_all " << deltay_all <<  std::endl;
              std::cout << "lam_hires_size: " << lam_hires_size <<  std::endl;
              std::cout << "phi_hires_size: " << phi_hires_size <<  std::endl;
              std::cout << "n_stretchedx +rim " << ((deltax_all+ 0.0001 - lam_hires_size) / delta_high_) << std::endl;
              std::cout << "n_stretchedy +rim  " << ((deltay_all+ 0.0001 - phi_hires_size) / delta_high_) << std::endl;
              std::cout << "rim x  : "<< n_x_rim_ << " y " << n_y_rim_ <<  std::endl;
              lambda_start = x_reg_start_;
              phi_start = y_reg_start_;
              check_x = startx_ + add_xf_ - lambda_start;
              check_y = starty_ + add_yf_ - phi_start;
              double check_st;

              check_st = n_stretchedx_ - n_stretchedy_;
              checkvalue(epsilon, check_x);
              checkvalue(epsilon, check_y);
              checkvalue(epsilon, check_st);

              nx_ = ((deltax_all + 0.0001) / delta_high_) + 1;
              ny_ = ((deltay_all + 0.0001) / delta_high_ ) + 1;
              std::cout << "xall  : "<< nx_ << " yall " << ny_ <<  std::endl;
         }

    rim_size = 0;
    std::cout << "call stretch nx_: " << nx_ << "n_stretchedx_: " << n_stretchedx_ << "n_x_rim_: " << n_x_rim_ << std::endl;
    std::cout << "call stretch ny_: " << ny_ << "n_stretchedy_: " << n_stretchedy_ << "n_y_rim_: " << n_y_rim_ << std::endl;
    crd[0] = general_stretch(crd[0], true, nx_, n_stretchedx_, n_x_rim_);
    crd[1] = general_stretch(crd[1], false, ny_, n_stretchedy_, n_y_rim_);

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
