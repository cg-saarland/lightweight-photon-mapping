#ifndef __LPM_PHOTONS_H
#define __LPM_PHOTONS_H

#include <mutex>
#include <atomic>
#include <mitsuba/mitsuba.h>
#include "emission_sampler.h"

MTS_NAMESPACE_BEGIN

namespace lpm {

template <typename T>
static T atomic_add(std::atomic<T>& a, T b) {
    T old_val = a.load();
    T desired_val = old_val + b;
    while(!a.compare_exchange_weak(old_val, desired_val))
        desired_val = old_val + b;
    return desired_val;
}

struct Photon {
    /// Weight carried by the photon, used in density estimation.
    Spectrum weight;

    /// Direction towards the previous vertex along the path.
    Vector dir;

    /// Position of the photon in world space.
    Point pos;

    /// Index within the photon map of the previous vertex along the path (or -1).
    /// The ancestor is not necesarily directly connected to this photon as not all
    /// intermediate photons are stored (e.g. on specular surfaces)
    size_t ancestor_idx;

    /// Number of edges between this photon and the light source
    uint32_t depth;

    /// Partial MIS weights for fast weight computation in the simple case
    Float partial_mis_sum;

    /// Whether this photon is more efficient than a path tracer sample
    /// \see{is_useful_photon}
    bool is_useful;

    /// Random numbers used to sample the emission of the path that generated this photon.
    PrimarySample emission_sample;

    /// The weight carried by the photon without the Le term.
    Float throughput;

    /// Contribution (e.g. luminance of image contribution) accumulated to be used with guiding approaches.
    CopyableAtomic<Float> guiding_contrib;
};

constexpr size_t NO_ANCESTOR = (size_t)(-1);

class PhotonContainer {
public:
    PhotonContainer(size_t capacity = 1000000) : next_(0), photons_(capacity) { }

    template <typename ConstIter>
    void add(ConstIter first, ConstIter last) {
        // TODO after the first (or second) iteration, the photon map will most likely not have to grow anymore,
        //      we can fix the size and replace the locks by atomic operations!
        //      In the unlikely event of having too many photons, those additional ones can be discarded with negligible bias

        // reserve memory, allocate new memory if necessary
        std::lock_guard<std::mutex> guard(write_lock_);
        auto n = last - first;
        auto start = next_;
        next_ += n;
        if (photons_.size() < next_) {
            SLog(EDebug, "Growing photon container to store %i photons...", next_);
            photons_.resize(2 * photons_.size());
            SAssert(next_ < photons_.size()); // the container should be large enough
        }

        // copy the photons and adjust the ancestor indices
        auto idx = start;
        for (auto i = first; i != last; ++i, ++idx) {
            photons_[idx] = *i;
            if (photons_[idx].ancestor_idx != NO_ANCESTOR)
                photons_[idx].ancestor_idx += start;
        }
    }

    /// Delets all photons in the container. No memory will be freed, it will instead be available for
    /// re-use during the next iteration.
    void reset() { next_ = 0; }

    Photon& operator[] (size_t idx) {
        SAssert(idx < next_);
        return photons_[idx];
    }

    const Photon& operator[] (size_t idx) const {
        SAssert(idx < next_);
        return photons_[idx];
    }

    size_t size() const { return next_; }

    Photon* begin() { return photons_.data(); }
    const Photon* begin() const { return photons_.data(); }

    Photon* end() { return photons_.data() + next_; }
    const Photon* end() const { return photons_.data() + next_; }

private:
    size_t next_;
    std::vector<Photon> photons_;
    std::mutex write_lock_;
};

} // namespace lpm

MTS_NAMESPACE_END

#endif // __LPM_PHOTONS_H