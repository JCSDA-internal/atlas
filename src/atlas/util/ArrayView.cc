/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <stdexcept>

#include "atlas/mesh/Field.h"
#include "atlas/util/Array.h"
#include "atlas/util/ArrayView.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

#define TEMPLATE_SPECIALIZATION( DATA_TYPE ) \
template<>\
ArrayView <DATA_TYPE, 0 >::ArrayView( const Array<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  rank_ = array.strides().size(); \
  strides_ = array.strides(); \
  shape_ = array.shape(); \
  size_ = std::accumulate(shape_.data(),shape_.data()+rank_,1,std::multiplies<int>()); \
} \
template<>\
ArrayView <DATA_TYPE, 0 >::ArrayView( const Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  rank_ = field.strides().size(); \
  strides_ = field.strides(); \
  shape_ = field.shape(); \
  size_ = std::accumulate(shape_.data(),shape_.data()+rank_,1,std::multiplies<int>()); \
} \
template<>\
ArrayView <DATA_TYPE, 1 >::ArrayView( const Array<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);  shape_[0]=array.shape(0);\
} \
template<> \
ArrayView <DATA_TYPE, 1 >::ArrayView( const Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);  shape_[0]=field.shape(0); \
} \
template<> \
ArrayView<DATA_TYPE,2>::ArrayView( const Array<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);       shape_[0]=array.shape(0); \
  strides_[1]=array.stride(1);       shape_[1]=array.shape(1); \
} \
template<> \
ArrayView<DATA_TYPE,2>::ArrayView( const Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);       shape_[0]=field.shape(0); \
  strides_[1]=field.stride(1);       shape_[1]=field.shape(1); \
} \
template<> \
ArrayView<DATA_TYPE,3>::ArrayView( const Array<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);       shape_[0]=array.shape(0); \
  strides_[1]=array.stride(1);       shape_[1]=array.shape(1); \
  strides_[2]=array.stride(2);       shape_[2]=array.shape(2); \
} \
template<> \
ArrayView<DATA_TYPE,3>::ArrayView( const Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);       shape_[0]=field.shape(0); \
  strides_[1]=field.stride(1);       shape_[1]=field.shape(1); \
  strides_[2]=field.stride(2);       shape_[2]=field.shape(2); \
} \
template<> \
ArrayView<DATA_TYPE,4>::ArrayView( const Array<DATA_TYPE>& array ) : data_( const_cast<DATA_TYPE*>(array.data()) ) \
{ \
  strides_[0]=array.stride(0);       shape_[0]=array.shape(0); \
  strides_[1]=array.stride(1);       shape_[1]=array.shape(1); \
  strides_[2]=array.stride(2);       shape_[2]=array.shape(2); \
  strides_[3]=array.stride(3);       shape_[3]=array.shape(3); \
} \
template<> \
ArrayView<DATA_TYPE,4>::ArrayView( const Field& field ) : data_( const_cast<DATA_TYPE*>(field.data<DATA_TYPE>()) ) \
{ \
  strides_[0]=field.stride(0);       shape_[0]=field.shape(0); \
  strides_[1]=field.stride(1);       shape_[1]=field.shape(1); \
  strides_[2]=field.stride(2);       shape_[2]=field.shape(2); \
  strides_[3]=field.stride(3);       shape_[3]=field.shape(3); \
}\

TEMPLATE_SPECIALIZATION(int);
TEMPLATE_SPECIALIZATION(float);
TEMPLATE_SPECIALIZATION(double);

#undef TEMPLATE_SPECIALIZATION

//------------------------------------------------------------------------------------------------------

} // namespace atlas