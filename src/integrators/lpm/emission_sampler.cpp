#include <fstream>
#include <mitsuba/render/scene.h>
#include "emission_sampler.h"

MTS_NAMESPACE_BEGIN

namespace lpm {

uint32_t next_pow2(uint32_t x) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

EmissionDistribution::EmissionDistribution(const Scene* sc, uint32_t guiding_start, uint32_t useful_start)
: iterations_(0)
, guiding_start_(guiding_start)
, useful_start_(useful_start)
, emitter_dist_all_(sc->getEmitters().size())
, emitter_dist_useful_(sc->getEmitters().size())
{
    // initialize the discrete sampling distribution for emitters
    uint32_t idx = 0;
    for (auto& em : sc->getEmitters()) {
        emitter_pdf_.append(em.get()->getSamplingWeight());
        emitter_obj_to_id_[em.get()] = idx++;

        // TODO degenerate lights can eliminate some (or all) dimensions

        // initial resolution is chosen s.t. every bin receives eight samples on average
        uint32_t num_lp = 1000000; // TODO determine correct value of number of light paths per iteration!
        uint32_t num_ls = sc->getEmitters().size();
        auto target_samples_per_bin = 8;
        auto blur_radius = 4;
        auto init_res = static_cast<uint32_t>(/*std::sqrt*/(std::sqrt((double)num_lp / (target_samples_per_bin * num_ls))));
        // ensure that the histogram is neither too small nor too large
        init_res = std::max(4u, std::min(256u, next_pow2(init_res)));

        emitter_histograms_all_.emplace_back(1, 1, /*init_res, init_res, */init_res, init_res, blur_radius);
        emitter_histograms_useful_.emplace_back(1, 1, /*init_res, init_res, */init_res, init_res, blur_radius);
    }
    emitter_pdf_.normalize();
}

Float EmissionDistribution::sample_primaries(Sampler* sampler, uint32_t& em_idx, Point2& pos_uv, Point2& dir_uv) {
    bool contrib_guiding = iterations_ >= guiding_start_ && iterations_ < useful_start_;
    bool useful_guiding = iterations_ >= useful_start_;

    Float em_pdf = 1.0f;

    auto em_u = sampler->next1D();
    if (contrib_guiding)
        em_idx = emitter_dist_all_.sample(em_u, em_pdf);
    else if (useful_guiding)
        em_idx = emitter_dist_useful_.sample(em_u, em_pdf);
    else
        em_idx = emitter_pdf_.sample(em_u, em_pdf);

    pos_uv = sampler->next2D();
    dir_uv = sampler->next2D();

    emitter_histograms_all_[em_idx].num_samples++;
    emitter_histograms_useful_[em_idx].num_samples++;

    auto bin_u = sampler->next1D();
    Float uv_pdf = 1.0f;
    if (useful_guiding)
        uv_pdf *= emitter_histograms_useful_[em_idx].transform(pos_uv, dir_uv, bin_u);
    else if (contrib_guiding)
        uv_pdf *= emitter_histograms_all_[em_idx].transform(pos_uv, dir_uv, bin_u);

    return em_pdf * uv_pdf;
}

Float EmissionDistribution::pdf_primaries(const Emitter* em, Point2& pos_uv, Point2& dir_uv) const {
    auto iter = emitter_obj_to_id_.find(em);
    SAssert(iter != emitter_obj_to_id_.end());
    auto idx = iter->second;

    bool contrib_guiding = iterations_ >= guiding_start_ && iterations_ < useful_start_;
    bool useful_guiding = iterations_ >= useful_start_;

    Float uvpdf = 1.0f;
    if (useful_guiding) {
        uvpdf = emitter_histograms_useful_[idx].pdf(pos_uv, dir_uv)
              * emitter_dist_useful_.pdf(idx);
    } else if (contrib_guiding) {
        uvpdf = emitter_histograms_all_[idx].pdf(pos_uv, dir_uv)
              * emitter_dist_all_.pdf(idx);
    } else {
        uvpdf = emitter_pdf_[idx];
    }

    return uvpdf;
}

void EmissionDistribution::add_record(const PrimarySample& sample, Float weight, bool useful) {
    // update the histogram on the emitter
    if (useful) emitter_histograms_useful_[sample.em_id].put(sample.pos_uv, sample.dir_uv, weight);
    emitter_histograms_all_[sample.em_id].put(sample.pos_uv, sample.dir_uv, weight);

    // update the emitter sampling distribution
    if (useful) emitter_dist_useful_.put(sample.em_id, weight);
    emitter_dist_all_.put(sample.em_id, weight);
}

void EmissionDistribution::normalize() {
    for (auto& h : emitter_histograms_all_) {
        h.subdiv();
        h.normalize();
    }
    for (auto& h : emitter_histograms_useful_) {
        h.subdiv();
        h.normalize();
    }
    emitter_dist_all_.normalize();
    emitter_dist_useful_.normalize();

    iterations_++;
}

void EmissionDistribution::dump(const std::string& filename, bool useful) const {
    std::ofstream str(filename, std::ios_base::binary);

    auto num_lights = emitter_histograms_all_.size();
    str.write(reinterpret_cast<const char*>(&num_lights), sizeof(num_lights));

    (useful ? emitter_dist_useful_ : emitter_dist_all_).dump(str);

    for (auto& h : (useful ? emitter_histograms_useful_ : emitter_histograms_all_))
        h.dump(str);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

Spectrum EmissionSampler::sample_emitter(Float time, PositionSamplingRecord& pos, DirectionSamplingRecord& dir, Float& uvpdf) {
    uvpdf = dist_->sample_primaries(sampler_, emission_sample.em_id, emission_sample.pos_uv, emission_sample.dir_uv);

    if (uvpdf == 0.0f) return Spectrum(0.0f);

    // pick an emitter
    auto em = scene_->getEmitters()[emission_sample.em_id].get();
    pos.object = em;

    // sample a position
    auto power = em->samplePosition(pos, emission_sample.pos_uv);

    // allow the integrator to take action
    auto medium  = em->getMedium();
    handle_emitter(medium, pos, power);

    // sample a direction
    power *= em->sampleDirection(dir, pos, emission_sample.dir_uv) / uvpdf;

    return power;
}

Float EmissionSampler::pdf_emit(const Emitter* em, const PositionSamplingRecord& pos, const Vector& dir) {
    auto pdf = em->pdfPosition(pos);
    auto posuv = em->samplePositionInv(pos);

    DirectionSamplingRecord dir_sample;
    dir_sample.measure = ESolidAngle;
    dir_sample.d = dir;
    pdf *= em->pdfDirection(dir_sample, pos);
    auto diruv = em->sampleDirectionInv(dir_sample, pos);

    auto uvpdf = dist_->pdf_primaries(em, posuv, diruv);

    return pdf * uvpdf;
}

} // namespace lpm

MTS_NAMESPACE_END