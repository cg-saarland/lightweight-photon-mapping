#ifndef __SPARSE_HISTOGRAM_H
#define __SPARSE_HISTOGRAM_H

#include <unordered_map>
#include <atomic>
#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

namespace lpm {

/// A datatype that allows copying std::atomic
/// The copy procedure itself is not atomic
template <typename T>
struct CopyableAtomic : public std::atomic<T> {
    CopyableAtomic() { }
    CopyableAtomic& operator= (const CopyableAtomic& a) { this->store(a.load()); return *this; }
    CopyableAtomic(const CopyableAtomic& a) { this->store(a.load()); }
    CopyableAtomic& operator= (T i) { this->store(i); return *this; }
};

/// Histogram on the 4D unit hypercube.
/// Uses hashing (std::unordered_map) to reduce memory and allow a very high resolution for sparse histograms.
/// This class can also be used to represent lower-dimensional histograms, by setting the resolution of some dimensions to one.
/// A simple box filter is applied to the data to reduce noise.
class SparseHistogram4D { // TODO rename class (no longer sparse)
public:
    /// Constructs a sparse histogram with the given resolution over the 4D unit hypercube.
    SparseHistogram4D(uint32_t res_u_prim, uint32_t res_v_prim, uint32_t res_u_sec, uint32_t res_v_sec, uint32_t filter_radius);

    void normalize();

    /// Subdivides the grid if it received a sufficient number of samples during the last iteration,
    /// has a variance above a threshold, and a resolution below the maximum.
    /// Invalidates the CDF and must be called before \c normalize()
    void subdiv();

    /// Transforms the given uniformly distributed samples to a distribution proportional to this histogram
    /// \return
    ///     The jacobian determinant of the transformation (i.e., the effective pdf on the unit hypercube)
    Float transform(Point2& prim, Point2& sec, Float bin) const; // TODO this should not need a 5th sample!

    /// Returns the pdf on the unit hypercube that corresponds to the pmf represented by this histogram.
    Float pdf(const Point2& prim, const Point2& sec) const;

    /// Adds the given weight to the corresponding bins (with filtering)
    void put(const Point2& prim, const Point2& sec, Float weight);

    Float& operator()(const Point2& prim, const Point2& sec);
    Float  operator()(const Point2& prim, const Point2& sec) const;

    Float& operator()(uint32_t u_prim, uint32_t v_prim, uint32_t u_sec, uint32_t v_sec);
    Float  operator()(uint32_t u_prim, uint32_t v_prim, uint32_t u_sec, uint32_t v_sec) const;

    void dump(std::ostream&) const;

    CopyableAtomic<uint64_t> num_samples;

    uint32_t bin(const Point2& prim, const Point2& sec) const;
    void coords(uint32_t bin_idx, Point2& prim_min, Point2& sec_min) const;

private:
    std::vector<Float> pdf_; // TODO rename to "weights" (its not a distribution)
    std::vector<Float> cdf_;

    std::vector<Float> pos_marginal_;
    std::vector<Float> dir_marginal_;

    Point2u  res_prim_;
    Point2u  res_sec_;
    Point2   inv_res_prim_;
    Point2   inv_res_sec_;
    Float    inv_normalization_;
    uint32_t filter_radius_;
    uint32_t num_bins_;
};

/// A one dimensional histogram in primary sample space.
class Histogram1D {
public:
    Histogram1D(uint32_t count);
    void put(uint32_t idx, Float w);
    void normalize();

    uint32_t sample(Float u, Float& pdf) const;

    Float pdf(uint32_t idx) const;

    void dump(std::ostream&) const;

private:
    std::vector<Float> weights_;
    std::vector<Float> cdf_;
    Float inv_normalization_;
};

} // namespace lpm

MTS_NAMESPACE_END

#endif // __SPARSE_HISTOGRAM_H