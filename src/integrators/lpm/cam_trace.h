#ifndef __CAM_TRACE_H
#define __CAM_TRACE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/imageblock.h>
#include <mitsuba/render/renderproc.h>
#include <mitsuba/core/sfcurve.h>
#include "emission_sampler.h"
#include "lpm.h"

MTS_NAMESPACE_BEGIN

struct CTWorkResult : public WorkResult {
    CTWorkResult(const Vector2i& size, const ReconstructionFilter* filter)
        : img(new ImageBlock(Bitmap::ESpectrum, size, filter))
        , img_useful(new ImageBlock(Bitmap::ESpectrum, size, filter))
    {}

    ref<ImageBlock> img;
    ref<ImageBlock> img_useful;

    void load(Stream *stream) override {
        img->load(stream);
        img_useful->load(stream);
    }

    /// Serialize a work result to a binary data stream
    void save(Stream *stream) const override {
        img->save(stream);
        img_useful->save(stream);
    }

    /// Return a string representation
    std::string toString() const override {
        return img->toString() + img_useful->toString();
    }

    MTS_DECLARE_CLASS()
};

namespace lpm {

class CamPathTracer : public EmissionSampler {
public:
    CamPathTracer(EmitDistPtr, int rr_depth, int max_depth, Scene*, Sampler*, Sensor*, VertexMemPool*, CTWorkResult*, const lpm::MisComputer*,
                  const lpm::PhotonContainer*, const lpm::HashGrid*, Float radius, int lp_count, int min_lp, bool vis_photons);

    void handle_surface(int depth, int nullInteractions, bool delta,
                        const Intersection &its, const Medium *medium,
                        const Spectrum &emitted, const Spectrum& throughput) override;

    void handle_medium(int depth, int nullInteractions, bool delta,
                       const MediumSamplingRecord &mRec, const Medium *medium,
                       const Vector &wi, const Spectrum &weight) override;

    void handle_miss(int depth, int nullInteractions, bool delta, const Spectrum& weight) override;

    Spectrum contrib;
    Spectrum contrib_useful;

protected:
    ref<CTWorkResult> result_;
    const lpm::MisComputer* mis_;

    const lpm::PhotonContainer* photons_;
    const lpm::HashGrid* photon_accel_;
    Float radius_;
    Float radius_sqr_;
    int   lp_count_;
    int   min_lp_;
    bool  vis_photons_;

    void merge(int depth, const BSDF*, const Intersection&, const Spectrum& throughput);
    void hit(int depth, const BSDF*, const Intersection&, const Spectrum& throughput);
    void nee(int depth, const BSDF*, const Intersection&, const Spectrum& throughput);
};

} // namespace lpm

/// Traces the requested number of light paths. Renders an image and builds a set of photons.
class CameraTracer : public WorkProcessor {
public:
    CameraTracer(const LPMConfig&, const lpm::PhotonContainer*, const lpm::HashGrid*, lpm::EmitDistPtr);
    CameraTracer(Stream*, InstanceManager*);

    void serialize(Stream*, InstanceManager*) const final;
    ref<WorkResult> createWorkResult() const final;
    ref<WorkProcessor> clone() const final;
    ref<WorkUnit> createWorkUnit() const final;
    void prepare() final;
    void process(const WorkUnit*, WorkResult*, const bool& stop) final;

    MTS_DECLARE_CLASS()

private:
    LPMConfig config_;
    ref<Scene>                scene_;
	ref<Sensor>               sensor_;
	ref<Sampler>              sampler_;
	ref<ReconstructionFilter> filter_;
    ref<CTWorkResult>         result_;
    HilbertCurve2D<uint8_t>   hilbert_curve_;

    // TODO this will not work in a distributed setting
    lpm::EmitDistPtr emit_dist_;
    const lpm::PhotonContainer* photons_;
    const lpm::HashGrid*        photon_accel_;
};

class CTProcess : public BlockedRenderProcess {
public:
    CTProcess(const RenderJob *job, RenderQueue *queue, const LPMConfig& cfg,
              const lpm::PhotonContainer*, const lpm::HashGrid*, lpm::EmitDistPtr);

    void processResult(const WorkResult *wr, bool cancelled) override;
	ref<WorkProcessor> createWorkProcessor() const override;
    void bindResource(const std::string &name, int id) override;

    MTS_DECLARE_CLASS()

    ref<ImageBlock> img;
    ref<ImageBlock> img_useful;

private:
    ref<Timer> refresh_timer_;
    LPMConfig config_;
    ref<Film> film_;

    lpm::EmitDistPtr emit_dist_;
    const lpm::PhotonContainer* photons_;
    const lpm::HashGrid* photon_accel_;
};

MTS_NAMESPACE_END

#endif // __CAM_TRACE_H