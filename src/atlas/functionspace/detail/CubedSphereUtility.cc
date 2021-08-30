/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereUtility.h"

namespace atlas {
namespace functionspace {
namespace detail {



CubedSphereUtility::CubedSphereUtility(
  idx_t halo, const Field& ijt, const Field& ghost) :
  halo_(halo), ijtView_(array::make_view<idx_t, 2>(ijt)),
  ghostView_(array::make_view<idx_t, 1>(ghost))
{

  // loop over ijt_ and find min and max ij bounds.
  for (idx_t index = 0; index < ijtView_.shape(0); ++index) {
    const auto i = ijtView_(index, CSIndex::I);
    const auto j = ijtView_(index, CSIndex::J);
    const auto t = idx2st(ijtView_(index, CSIndex::T));

    ijBounds_[t].iBegin = std::min(i    , ijBounds_[t].iBegin);
    ijBounds_[t].jBegin = std::min(j    , ijBounds_[t].jBegin);
    ijBounds_[t].iEnd   = std::max(i + 1, ijBounds_[t].iEnd  );
    ijBounds_[t].jEnd   = std::max(j + 1, ijBounds_[t].jEnd  );

  }

  // Set ijtToIdx vectors
  for (idx_t t = 0; t < 6; ++t) {

    // Set data array.
    const auto vecSize = (j_end_halo(t) - j_begin_halo(t))
                       * (i_end_halo(t) - i_begin_halo(t));
    ijtToIdx_.push_back(std::vector<idx_t>(idx2st(vecSize), invalid_index()));

  }

  // loop over ijt_ and set ijtToIdx
  for (idx_t index = 0; index < ijtView_.shape(0); ++index) {
    const auto i = ijtView_(index, CSIndex::I);
    const auto j = ijtView_(index, CSIndex::J);
    const auto t = ijtView_(index, CSIndex::T);

   ijtToIdx_[idx2st(t)][vecIndex(i, j, t)] = index;
  }
}

CubedSphereUtility::~CubedSphereUtility() {};

idx_t CubedSphereUtility::i_begin(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[idx2st(t)].iBegin + halo_;
}

idx_t CubedSphereUtility::i_end(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[idx2st(t)].iEnd - halo_;
}

idx_t CubedSphereUtility::i_begin_halo(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[idx2st(t)].iBegin;
}

idx_t CubedSphereUtility::i_end_halo(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[idx2st(t)].iEnd;
}

idx_t CubedSphereUtility::j_begin(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[idx2st(t)].jBegin + halo_;
}

idx_t CubedSphereUtility::j_end(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[idx2st(t)].jEnd - halo_;
}

idx_t CubedSphereUtility::j_begin_halo(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[idx2st(t)].jBegin;
}

idx_t CubedSphereUtility::j_end_halo(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[idx2st(t)].jEnd;
}

idx_t CubedSphereUtility::index(idx_t i, idx_t j, idx_t t) const {

  // Check bounds.
  iBoundsCheck(i, t);
  jBoundsCheck(j, t);

  return ijtToIdx_[idx2st(t)][vecIndex(i, j, t)];
}

void CubedSphereUtility::tBoundsCheck(idx_t t) const {
  if (t < 0 or t > 5) throw_OutOfRange("t", t, 6);
}

void CubedSphereUtility::jBoundsCheck(idx_t j, idx_t t) const {
  const auto jSize = j_end_halo(t) - j_begin_halo(t);
  j -= j_begin_halo(t);
  if (j < 0 or j >= jSize) throw_OutOfRange("j - jMin", j, jSize);
}

void CubedSphereUtility::iBoundsCheck(idx_t i, idx_t t) const {
  const auto iSize = i_end_halo(t) - i_begin_halo(t);
  i -= i_begin_halo(t);
  if (i < 0 or i >= iSize) throw_OutOfRange("i - iMin", i, iSize);
}

size_t CubedSphereUtility::vecIndex(idx_t i, idx_t j, idx_t t) const {
  return idx2st((j - j_begin_halo(t)) * (i_end_halo(t) - i_begin_halo(t))
               + i - i_begin_halo(t));
}



} // detail
} // functionspace
} // atlas
