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

#include <algorithm>
#include <any>

#include "atlas/library/defines.h"
#if ATLAS_HAVE_EIGEN
#include <Eigen/Sparse>
#endif
#include "eckit/linalg/SparseMatrix.h"

#include "atlas/array.h"

namespace atlas {
namespace linalg {

//----------------------------------------------------------------------------------------------------------------------

template <typename Value, typename Index = eckit::linalg::Index>
struct SparseMatrixConvertor : public SparseMatrixStorage {
    using value_type = Value;
    using index_type = Index;

private:

    template<typename InputT, typename OutputT>
    void host_copy(const InputT* input_data, atlas::array::Array& output) {
        auto size = output.size();
        OutputT* output_data = output.host_data<OutputT>();
        std::copy( input_data, input_data + size, output_data );
    }

    template<typename InputT, typename OutputT>
    void host_copy(const atlas::array::Array& input, atlas::array::Array& output) {
        host_copy<InputT,OutputT>( input.host_data<InputT>(), output );
    }

public:
    virtual ~SparseMatrixConvertor() = default;

    SparseMatrixConvertor() = default;

    explicit SparseMatrixConvertor(eckit::linalg::SparseMatrix&& m) {
        nnz_   = m.nonZeros();
        rows_  = m.rows();
        cols_  = m.cols();

        if constexpr(std::is_same_v<index_type, eckit::linalg::Index> && std::is_same_v<value_type, eckit::linalg::Scalar>) {
            outer_.reset(atlas::array::Array::wrap(const_cast<eckit::linalg::Index*>(m.outer()), atlas::array::make_shape(rows_+1)));
            inner_.reset(atlas::array::Array::wrap(const_cast<eckit::linalg::Index*>(m.inner()), atlas::array::make_shape(nnz_)));
            value_.reset(atlas::array::Array::wrap(const_cast<eckit::linalg::Scalar*>(m.data()), atlas::array::make_shape(nnz_)));

            // We now move the eckit::linalg::SparseMatrix into a generic storage so
            //   the wrapped array data does not go out of scope
            storage_ = std::make_any<eckit::linalg::SparseMatrix>(std::move(m));
        }
        else {
            outer_.reset(atlas::array::Array::create<index_type>(rows_+1));
            inner_.reset(atlas::array::Array::create<index_type>(nnz_));
            value_.reset(atlas::array::Array::create<value_type>(nnz_));
            host_copy<eckit::linalg::Index, index_type>(m.outer(), *outer_);
            host_copy<eckit::linalg::Index, index_type>(m.inner(), *inner_);
            host_copy<eckit::linalg::Scalar,value_type>(m.data(),  *value_);
        }
    }


    /// Create copy of eckit::linalg::SparseMatrix
    explicit SparseMatrixConvertor(const eckit::linalg::SparseMatrix& other) {
        nnz_   = other.nonZeros();
        rows_  = other.rows();
        cols_  = other.cols();
        outer_.reset(atlas::array::Array::create<index_type>(rows_+1));
        inner_.reset(atlas::array::Array::create<index_type>(nnz_));
        value_.reset(atlas::array::Array::create<value_type>(nnz_));
        host_copy<eckit::linalg::Index, index_type>(other.outer(), *outer_);
        host_copy<eckit::linalg::Index, index_type>(other.inner(), *inner_);
        host_copy<eckit::linalg::Scalar,value_type>(other.data(),  *value_);
    }


    template <typename ValueT, typename IndexT>
    SparseMatrixConvertor(const SparseMatrixConvertor<ValueT,IndexT>& other) {
        nnz_   = other.nnz();
        rows_  = other.rows();
        cols_  = other.cols();
        outer_.reset(atlas::array::Array::create<index_type>(rows_+1));
        inner_.reset(atlas::array::Array::create<index_type>(nnz_));
        value_.reset(atlas::array::Array::create<value_type>(nnz_));
        host_copy<IndexT,index_type>(other.outer(), *outer_);
        host_copy<IndexT,index_type>(other.inner(), *inner_);
        host_copy<ValueT,value_type>(other.value(), *value_);
    }

    SparseMatrixConvertor(SparseMatrixConvertor&& other) {
        nnz_   = other.nnz_;
        rows_  = other.rows_;
        cols_  = other.cols_;
        outer_ = std::move(other.outer_);
        inner_ = std::move(other.inner_);
        value_ = std::move(other.value_);
        other.nnz_ = 0;
        other.rows_ = 0;
        other.cols_ = 0;
        storage_ = std::move(other.storage_);
    }

#if ATLAS_HAVE_EIGEN
    template <typename ValueT, typename IndexT>
    SparseMatrixConvertor(const Eigen::SparseMatrix<ValueT, Eigen::RowMajor, IndexT>& m) {
        nnz_   = m.nonZeros();
        rows_  = m.rows();
        cols_  = m.cols();
        outer_.reset(atlas::array::Array::create<Index>(rows_+1));
        inner_.reset(atlas::array::Array::create<Index>(nnz_));
        value_.reset(atlas::array::Array::create<Value>(nnz_));
        host_copy<IndexT,index_type>(m.outerIndexPtr(), *outer_);
        host_copy<IndexT,index_type>(m.innerIndexPtr(), *inner_);
        host_copy<ValueT,value_type>(m.valuePtr(),      *value_);
    }

    template <typename ValueT, typename IndexT>
    SparseMatrixConvertor(Eigen::SparseMatrix<ValueT, Eigen::RowMajor, IndexT>&& m) {
        nnz_   = m.nonZeros();
        rows_  = m.rows();
        cols_  = m.cols();

        if constexpr(std::is_same_v<std::decay_t<index_type>, std::decay_t<IndexT>> && std::is_same_v<std::decay_t<value_type>, std::decay_t<ValueT>>) {

            outer_.reset(atlas::array::Array::wrap(const_cast<index_type*>(m.outerIndexPtr()), atlas::array::make_shape(rows_+1)));
            inner_.reset(atlas::array::Array::wrap(const_cast<index_type*>(m.innerIndexPtr()), atlas::array::make_shape(nnz_)));
            value_.reset(atlas::array::Array::wrap(const_cast<value_type*>(m.valuePtr()),      atlas::array::make_shape(nnz_)));

            // We now move the eckit::linalg::SparseMatrix into a generic storage so
            //   the wrapped array data does not go out of scope
            //   Note: Eigen move constructor only available since 3.5; use swap instead.
            using EigenMatrix = Eigen::SparseMatrix<ValueT, Eigen::RowMajor, IndexT>;
            auto m_ptr = std::make_shared<EigenMatrix>();
            m_ptr->swap(m);
            storage_ = std::make_any<std::shared_ptr<EigenMatrix>>(std::move(m_ptr));
        }
        else {
            outer_.reset(atlas::array::Array::create<index_type>(rows_+1));
            inner_.reset(atlas::array::Array::create<index_type>(nnz_));
            value_.reset(atlas::array::Array::create<value_type>(nnz_));
            host_copy<IndexT, index_type>(m.outerIndexPtr(), *outer_);
            host_copy<IndexT, index_type>(m.innerIndexPtr(), *inner_);
            host_copy<ValueT, value_type>(m.valuePtr(),      *value_);
        }
    }
#endif

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace linalg
}  // namespace atlas