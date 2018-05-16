#ifndef __EMISSION_SAMPLER_H
#define __EMISSION_SAMPLER_H

#include <memory>
#include <unordered_map>
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/pmf.h>
#include "path_sampler.h"
#include "sparse_histogram.h"

MTS_NAMESPACE_BEGIN

namespace lpm {

/// Random numbers on the 4D hypercube used to sample a position and direction from an emitter.
struct PrimarySample {
    uint32_t em_id;
    Point2 pos_uv;
    Point2 dir_uv;

    PrimarySample() : em_id(-1) {}
};

/// Stores the emission guiding information
class EmissionDistribution {
public:
    EmissionDistribution(const Scene*, uint32_t guiding_start, uint32_t useful_start);

    Float sample_primaries(Sampler*, uint32_t& em, Point2& pos_uv, Point2& dir_uv);
    Float pdf_primaries(const Emitter* em, Point2& pos_uv, Point2& dir_uv) const;

    void add_record(const PrimarySample&, Float w, bool useful);
    void normalize();

    void dump(const std::string& filename, bool useful) const;

private:
    uint32_t iterations_;
    const uint32_t guiding_start_;
    const uint32_t useful_start_;

    /// total power based importance sampling of emitters
    DiscreteDistribution emitter_pdf_;

    std::unordered_map<const Emitter*, uint32_t> emitter_obj_to_id_;

    std::vector<SparseHistogram4D> emitter_histograms_all_;
    std::vector<SparseHistogram4D> emitter_histograms_useful_;
    Histogram1D emitter_dist_all_;
    Histogram1D emitter_dist_useful_;
};

using EmitDistPtr = std::shared_ptr<EmissionDistribution>;

/// Importance samples the emission from the light sources
class EmissionSampler : public PathSampler {
public:
    Spectrum sample_emitter(Float time, PositionSamplingRecord&, DirectionSamplingRecord&, Float& uvpdf) override;
    Float pdf_emit(const Emitter* em, const PositionSamplingRecord& pos, const Vector& dir) override;

    /// The random numbers that were used to sample the last path traced with this instance (if any).
    PrimarySample emission_sample;

protected:
    template <typename... Args>
    EmissionSampler(EmitDistPtr d, Args&&... args)
    : dist_(d), PathSampler(std::forward<Args>(args)...) {}

private:
    EmitDistPtr dist_;
};

} // namespace lpm

MTS_NAMESPACE_END

#endif // __EMISSION_SAMPLER_H