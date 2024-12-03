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
#include <memory>

#include "atlas/library/defines.h"
#if ATLAS_HAVE_EIGEN
#include <Eigen/Sparse>
#endif
#include "eckit/linalg/SparseMatrix.h"

#include "atlas/array.h"
#include "atlas/runtime/Exception.h"
#include "atlas/linalg/sparse/SparseMatrixView.h"

namespace atlas {
namespace linalg {

//----------------------------------------------------------------------------------------------------------------------

/// @brief SparseMatrixStorage
/// Storage class based on atlas::array::Array with GPU offload capability
/// This class contains only the data and no operators, iterators, or special
/// sparse matrix construction from triplets etc.
///
/// The construction is handled via SparseMatrixConvertor class.
/// Examples:
///
/// To construct it, taking ownership of a constructed Eigen::SparseMatrix
///
///     SparseMatrixStorage s {std::move(eigen_matrix)};
///
/// To construct it, taking a copy of a constructed Eigen::SparseMatrix
///
///     SparseMatrixStorage s {eigen_matrix};
///
/// To construct it, taking ownership of a constructed eckit::linalg::SparseMatrix, avoiding copies if data types match
///
///     SparseMatrixStorage s {std::move(eckit_matrix)};
///
///
/// To construct it, taking a copy of a constructed eckit::linalg::SparseMatrix (no std::move)
///
///     SparseMatrixStorage s {eckit_matrix};
///
///
/// It is also possible to initialise empty and move into it at later stage:
///
///     SparseMatrixStorage s;
///     s = SparseMatrixStorage{eckit_matrix};
///
///
/// To construct it, taking a single precision copy of a constructed eckit::linalg::SparseMatrix
///
///     SparseMatrixStorage s {SparseMatrixConvertor<float>{eckit_matrix}};
///

class SparseMatrixStorage {
public:

    // Virtual destructor
    virtual ~SparseMatrixStorage() = default;

    /// Default constructor
    SparseMatrixStorage() = default;

    /// Move constructor, takes ownership!
    SparseMatrixStorage(SparseMatrixStorage&& other);

    /// Copy constructor, makes copy!
    SparseMatrixStorage(const SparseMatrixStorage& other);

    /// Copy from a SparseMatrixView
    template<typename Value, typename Index>
    SparseMatrixStorage(const SparseMatrixView<Value,Index>& host_view);

    /// Move constructor from eckit::linalg::SparseMatrix, takes ownership!
    explicit SparseMatrixStorage(eckit::linalg::SparseMatrix&& eckit_matrix);

    /// Copy constructor from eckit::linalg::SparseMatrix, makes copy!
    explicit SparseMatrixStorage(const eckit::linalg::SparseMatrix& eckit_matrix);

    /// Move constructor from Eigen::SparseMatrix, takes ownership!
    template<typename Value, typename Index>
    explicit SparseMatrixStorage(Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>&& eigen_matrix);

    /// Copy constructor from Eigen::SparseMatrix, makes copy!
    template<typename Value, typename Index>
    explicit SparseMatrixStorage(const Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>& eigen_matrix);

    /// Move assign from other SparseMatrixStorage, takes ownership!
    SparseMatrixStorage& operator=(SparseMatrixStorage&& other);

    /// Copy assign from other SparseMatrixStorage, takes ownership!
    SparseMatrixStorage& operator=(const SparseMatrixStorage& other);

    /// Empty if rows and cols are zero.
    bool empty() const { return rows_ == 0 && cols_ == 0; }

    /// Footprint in bytes in host memory space
    std::size_t footprint() const { return value_->footprint() + outer_->footprint() + inner_->footprint(); }

    const atlas::array::Array& value() const { return *value_; }
    const atlas::array::Array& outer() const { return *outer_; }
    const atlas::array::Array& inner() const { return *inner_; }
    std::size_t rows() const { return rows_;}
    std::size_t cols() const { return cols_;}
    std::size_t nnz()  const { return nnz_;}

    void updateDevice() const;

    void updateHost() const;

    bool hostNeedsUpdate() const;

    bool deviceNeedsUpdate() const;

    void setHostNeedsUpdate(bool v) const;

    void setDeviceNeedsUpdate(bool v) const;

    bool deviceAllocated() const;

    void allocateDevice() const;

    void deallocateDevice() const;
    
protected:
    std::size_t nnz_{0};
    std::size_t rows_{0};
    std::size_t cols_{0};
    std::unique_ptr<atlas::array::Array> outer_;
    std::unique_ptr<atlas::array::Array> inner_;
    std::unique_ptr<atlas::array::Array> value_;

    std::any storage_;
        // This storage is usually empty.
        // It is used to allow to move alternative sparse matrix formats into it,
        // and then wrap this data using the atlas Arrays

private:
    template<typename ValueT>
    void host_copy(const ValueT* input_data, atlas::array::Array& output) {
        auto size = output.size();
        ValueT* output_data = output.host_data<ValueT>();
        std::copy( input_data, input_data + size, output_data );
    }
};

//----------------------------------------------------------------------------------------------------------------------

template<typename Value, typename Index>
SparseMatrixStorage::SparseMatrixStorage(const SparseMatrixView<Value,Index>& host_view) {
    nnz_   = host_view.nnz_;
    rows_  = host_view.rows_;
    cols_  = host_view.cols_;
    outer_.reset(atlas::array::Array::create<Index>(host_view.outer_size()));
    inner_.reset(atlas::array::Array::create<Index>(host_view.inner_size()));
    value_.reset(atlas::array::Array::create<Value>(host_view.value_size()));
    host_copy(host_view.outer(),*outer_);
    host_copy(host_view.inner(),*inner_);
    host_copy(host_view.value(),*value_);
}

//----------------------------------------------------------------------------------------------------------------------

#if ATLAS_HAVE_EIGEN
template <typename Value, typename Index>
SparseMatrixStorage::SparseMatrixStorage(Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>&& m) {
    nnz_   = m.nonZeros();
    rows_  = m.rows();
    cols_  = m.cols();

    outer_.reset(atlas::array::Array::wrap(const_cast<Index*>(m.outerIndexPtr()), atlas::array::make_shape(rows_+1)));
    inner_.reset(atlas::array::Array::wrap(const_cast<Index*>(m.innerIndexPtr()), atlas::array::make_shape(nnz_)));
    value_.reset(atlas::array::Array::wrap(const_cast<Value*>(m.valuePtr()),      atlas::array::make_shape(nnz_)));

    // We now move the eckit::linalg::SparseMatrix into a generic storage so
    //   the wrapped array data does not go out of scope
    //   Note: Eigen move constructor only available since 3.5; use swap instead.
    using EigenMatrix = Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>;
    auto m_ptr = std::make_shared<EigenMatrix>();
    m_ptr->swap(m);
    storage_ = std::make_any<std::shared_ptr<EigenMatrix>>(std::move(m_ptr));
}

template <typename Value, typename Index>
SparseMatrixStorage::SparseMatrixStorage(const Eigen::SparseMatrix<Value, Eigen::RowMajor, Index>& m) {
    nnz_   = m.nonZeros();
    rows_  = m.rows();
    cols_  = m.cols();

    outer_.reset(atlas::array::Array::create<Index>(rows_+1));
    inner_.reset(atlas::array::Array::create<Index>(nnz_));
    value_.reset(atlas::array::Array::create<Value>(nnz_));
    host_copy(m.outerIndexPtr(), *outer_);
    host_copy(m.innerIndexPtr(), *inner_);
    host_copy(m.valuePtr(),      *value_);
}
#endif

//----------------------------------------------------------------------------------------------------------------------

template<typename Value, typename Index = eckit::linalg::Index>
inline SparseMatrixView<Value,Index> make_host_view(const SparseMatrixStorage& m) {
    if( m.value().datatype().kind() != DataType::kind<Value>() || m.outer().datatype().kind() != DataType::kind<Index>() ) {
        ATLAS_THROW_EXCEPTION("Cannot make_host_view<" + DataType::str<Value>() + "," << DataType::str<Index>() +
            ">(const SparseMatrixStorage&) from SparseMatrixStorage containing values of type <" + m.value().datatype().str() + "> and indices of type <" + m.outer().datatype().str() +">" );
    }
    return SparseMatrixView<Value,Index> {
        m.rows(),
        m.cols(),
        m.nnz(),
        m.value().host_data<Value>(),
        m.inner().host_data<Index>(),
        m.outer().host_data<Index>()
    };
}

//----------------------------------------------------------------------------------------------------------------------

template<typename Value, typename Index = eckit::linalg::Index>
inline SparseMatrixView<Value,Index> make_device_view(const SparseMatrixStorage& m) {
    if( m.value().datatype().kind() != DataType::kind<Value>() || m.outer().datatype().kind() != DataType::kind<Index>() ) {
        ATLAS_THROW_EXCEPTION("Cannot make_device_view<" + DataType::str<Value>() + "," << DataType::str<Index>() +
            ">(const SparseMatrixStorage&) from SparseMatrixStorage containing values of type <" + m.value().datatype().str() + "> and indices of type <" + m.outer().datatype().str() +">" );
    }
    if( ! m.deviceAllocated() ) {
        ATLAS_THROW_EXCEPTION("Cannot make_device_view(const SparseMatrixStorage&) as the device data is not allocated");
    }
    return SparseMatrixView<Value,Index>{
        m.rows(),
        m.cols(),
        m.nnz(),
        m.value().device_data<Value>(),
        m.inner().device_data<Index>(),
        m.outer().device_data<Index>()
    };
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace linalg
}  // namespace atlas