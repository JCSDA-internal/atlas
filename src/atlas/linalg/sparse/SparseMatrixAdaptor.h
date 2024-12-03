/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/library/defines.h"
#if ATLAS_HAVE_EIGEN
#include <Eigen/Sparse>
#endif
#include "eckit/linalg/SparseMatrix.h"

#include "atlas/array.h"

namespace atlas {
namespace linalg {

//----------------------------------------------------------------------------------------------------------------------

/// @brief Wrap existing sparse matrix formats into SparseMatrixStorage.
///        This means that the host_data is not owned, but device_data is owned
class SparseMatrixAdaptor : public SparseMatrixStorage {
public:
    SparseMatrixAdaptor(const eckit::linalg::SparseMatrix& m) {
        nnz_   = m.nonZeros();
        rows_  = m.rows();
        cols_  = m.cols();
        outer_.reset(atlas::array::Array::wrap(const_cast<eckit::linalg::Index*>(m.outer()), atlas::array::make_shape(rows_+1)));
        inner_.reset(atlas::array::Array::wrap(const_cast<eckit::linalg::Index*>(m.inner()), atlas::array::make_shape(nnz_)));
        value_.reset(atlas::array::Array::wrap(const_cast<eckit::linalg::Scalar*>(m.data()), atlas::array::make_shape(nnz_)));
    }

#if ATLAS_HAVE_EIGEN
    template <typename Value, typename Index>
    SparseMatrixAdaptor(const Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>& m) {
        nnz_   = m.nonZeros();
        rows_  = m.rows();
        cols_  = m.cols();
        outer_.reset(atlas::array::Array::wrap(const_cast<std::decay_t<Index>*>(m.outerIndexPtr()), atlas::array::make_shape(rows_+1)));
        inner_.reset(atlas::array::Array::wrap(const_cast<std::decay_t<Index>*>(m.innerIndexPtr()), atlas::array::make_shape(nnz_)));
        value_.reset(atlas::array::Array::wrap(const_cast<std::decay_t<Value>*>(m.valuePtr()),      atlas::array::make_shape(nnz_)));
    }
#endif
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace linalg
}  // namespace atlas