/*
 * (C) Copyright 2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

#include "atlas/runtime/Exception.h"
#include "atlas/util/detail/KDTree.h"
#include "atlas/util/detail/nanoflann.hpp"

namespace atlas {
namespace util {
namespace detail {

//------------------------------------------------------------------------------------------------------

namespace nanoflann_detail {

/// Like nanoflann::RadiusResultSet, but includes points at exactly the search radius (as eckit::KDTree does),
/// and stores (index, squared distance) pairs. nanoflann only offers points with `dist < worstDist()` to the
/// result set, so worstDist() is nudged up by one ulp.
class InclusiveRadiusResultSet {
public:
    using DistanceType = double;
    using IndexType    = size_t;

    InclusiveRadiusResultSet(DistanceType radius_squared, std::vector<std::pair<IndexType, DistanceType>>& out):
        radius_squared_(radius_squared),
        worst_(std::nextafter(radius_squared, std::numeric_limits<DistanceType>::infinity())),
        out_(out) {
        out_.clear();
    }

    size_t size() const { return out_.size(); }
    bool empty() const { return out_.empty(); }
    bool full() const { return true; }

    bool addPoint(DistanceType dist, IndexType index) {
        if (dist <= radius_squared_) {
            out_.emplace_back(index, dist);
        }
        return true;
    }

    DistanceType worstDist() const { return worst_; }

    void sort() {
        std::sort(out_.begin(), out_.end(),
                  [](const std::pair<IndexType, DistanceType>& a, const std::pair<IndexType, DistanceType>& b) {
                      return a.second < b.second;
                  });
    }

private:
    DistanceType radius_squared_;
    DistanceType worst_;
    std::vector<std::pair<IndexType, DistanceType>>& out_;
};

}  // namespace nanoflann_detail

//------------------------------------------------------------------------------------------------------

/// @brief KDTreeBase implementation backed by nanoflann (https://github.com/jlblancoc/nanoflann)
///
/// Behaves like KDTree_eckit (eckit::KDTreeMemory):
///   - After reserve(), inserted values are deferred: they are not counted by size(), and queries throw, until
///     build() is called.
///   - Without reserve(), inserted values are immediately visible to queries.
///   - closestPointsWithinRadius() includes points at exactly the given radius.
///   - Queries on an empty tree throw eckit::AssertionFailed.
///
/// Implementation notes and differences with KDTree_eckit:
///   - Points and payloads are stored once, compactly; nanoflann reads coordinates in place.
///     There is no temporary staging buffer, so peak memory during build() equals the final footprint.
///   - The index is (re)built lazily, at most once, on the first query after insertions (or by build()).
///     This makes "insert all, then query" cheap, but alternating single inserts and queries costs
///     O(n^2 log n) in total, compared to O(n log n) for eckit.
///   - A query may trigger that rebuild, so a const query can be expensive and briefly takes a lock.
///     Queries are otherwise const and thread-safe.
///   - Points at equal distance from the query point may be returned in a different order than with eckit.
template <typename PayloadT, typename PointT = Point3>
class KDTree_nanoflann : public KDTreeBase<PayloadT, PointT> {
    using Base = KDTreeBase<PayloadT, PointT>;

public:
    using Point       = typename Base::Point;
    using Payload     = typename Base::Payload;
    using PayloadList = typename Base::PayloadList;
    using Value       = typename Base::Value;
    using ValueList   = typename Base::ValueList;

    using Base::build;
    using Base::closestPoint;
    using Base::closestPoints;
    using Base::closestPointsWithinRadius;
    using Base::insert;
    using Base::reserve;

    static constexpr int DIMS = static_cast<int>(Point::DIMS);

public:
    KDTree_nanoflann() = default;
    explicit KDTree_nanoflann(const Geometry& geometry): Base(geometry) {}

    // The nanoflann index refers to this object's storage
    KDTree_nanoflann(const KDTree_nanoflann&)            = delete;
    KDTree_nanoflann& operator=(const KDTree_nanoflann&) = delete;

    idx_t size() const override { return static_cast<idx_t>(values_.size() - pending_count_); }

    /// @brief Memory used by stored points/payloads and by the index
    size_t footprint() const override {
        ensure_index();
        return values_.capacity() * sizeof(Entry) + (index_ ? index_->usedMemory(*index_) : 0);
    }

    void reserve(idx_t size) override {
        values_.reserve(values_.size() + static_cast<size_t>(size));
        pending_build_ = pending_build_ || size > 0;
    }

    void insert(const Value& value) override {
        values_.push_back(Entry{value.point(), value.payload()});
        if (pending_build_) {
            ++pending_count_;
        }
        dirty_.store(true, std::memory_order_release);
    }

    void build() override {
        pending_build_ = false;
        pending_count_ = 0;
        std::lock_guard<std::mutex> lock(mutex_);
        rebuild_index();
    }

    /// @brief Build with given values; previously inserted values are discarded
    void build(std::vector<Value>& values) override {
        values_.clear();
        values_.reserve(values.size());
        for (const auto& value : values) {
            values_.push_back(Entry{value.point(), value.payload()});
        }
        build();
    }

private:
    ValueList do_closestPoints(const Point& p, size_t k) const override {
        assert_searchable();
        k = std::min(k, values_.size());
        ValueList result;
        if (k == 0) {
            return result;
        }
        // Avoid heap allocations for typical small k
        constexpr size_t stack_size = 32;
        size_t stack_indices[stack_size];
        double stack_distances[stack_size];
        std::vector<size_t> heap_indices;
        std::vector<double> heap_distances;
        size_t* indices   = stack_indices;
        double* distances = stack_distances;
        if (k > stack_size) {
            heap_indices.resize(k);
            heap_distances.resize(k);
            indices   = heap_indices.data();
            distances = heap_distances.data();
        }
        nanoflann::KNNResultSet<double, size_t, size_t> result_set(k);
        result_set.init(indices, distances);
        index_->findNeighbors(result_set, p.data());
        const size_t found = result_set.size();
        result.reserve(found);
        for (size_t i = 0; i < found; ++i) {
            const Entry& e = values_[indices[i]];
            result.emplace_back(e.point, e.payload, std::sqrt(distances[i]));
        }
        return result;
    }

    Value do_closestPoint(const Point& p) const override {
        auto result = do_closestPoints(p, 1);
        return result.front();
    }

    ValueList do_closestPointsWithinRadius(const Point& p, double radius) const override {
        assert_searchable();
        ValueList result;
        if (radius < 0.) {
            return result;
        }
        std::vector<std::pair<size_t, double>> found;
        nanoflann_detail::InclusiveRadiusResultSet result_set(radius * radius, found);
        index_->findNeighbors(result_set, p.data(), nanoflann::SearchParameters(0.f, /*sorted=*/true));
        result.reserve(found.size());
        for (const auto& item : found) {
            const Entry& e = values_[item.first];
            result.emplace_back(e.point, e.payload, std::sqrt(item.second));
        }
        return result;
    }

private:
    struct Entry {
        Point point;
        Payload payload;
    };

    /// nanoflann dataset adaptor, reading coordinates in place from the stored entries
    struct Dataset {
        const Entry* data{nullptr};
        size_t size{0};

        size_t kdtree_get_point_count() const { return size; }
        double kdtree_get_pt(const size_t idx, const size_t dim) const { return data[idx].point.data()[dim]; }
        template <class BBOX>
        bool kdtree_get_bbox(BBOX&) const {
            return false;
        }
    };

    using Metric = nanoflann::L2_Simple_Adaptor<double, Dataset, double, size_t>;
    using Index  = nanoflann::KDTreeSingleIndexAdaptor<Metric, Dataset, DIMS, size_t>;

    static constexpr size_t leaf_max_size = 10;

    /// Throws if reserve() was called without a subsequent build(), as KDTree_eckit does; also checks non-empty.
    /// Ensures the index is up to date.
    void assert_searchable() const {
        if (pending_build_) {
            throw_AssertionFailed("KDTree was used before calling build()");
        }
        ATLAS_ASSERT(!values_.empty(), "Cannot search an empty k-d tree");
        ensure_index();
    }

    /// Rebuild index if values were inserted since the last build. Thread-safe (double-checked locking).
    void ensure_index() const {
        if (dirty_.load(std::memory_order_acquire)) {
            std::lock_guard<std::mutex> lock(mutex_);
            if (dirty_.load(std::memory_order_relaxed)) {
                rebuild_index();
            }
        }
    }

    /// Caller must hold mutex_ or have exclusive access
    void rebuild_index() const {
        index_.reset();
        dataset_.data = values_.data();
        dataset_.size = values_.size();
        if (!values_.empty()) {
            // The constructor builds the index (single-threaded)
            index_.reset(new Index(DIMS, dataset_, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size)));
        }
        dirty_.store(false, std::memory_order_release);
    }

private:
    std::vector<Entry> values_;  // the index reads coordinates in place from these entries
    bool pending_build_{false};  // reserve() was called, build() not yet
    size_t pending_count_{0};    // number of values inserted while pending_build_

    mutable Dataset dataset_;
    mutable std::unique_ptr<Index> index_;
    mutable std::mutex mutex_;
    mutable std::atomic<bool> dirty_{false};
};

//------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace util
}  // namespace atlas
