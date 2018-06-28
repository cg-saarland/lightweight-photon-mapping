#ifndef __LPM_LIGHT_TRACE_H
#define __LPM_LIGHT_TRACE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/render/imageblock.h>
#include <mitsuba/render/range.h>

#include "lpm.h"
#include "photon_map.h"
#include "emission_sampler.h"

MTS_NAMESPACE_BEGIN

/// The LT result is an (MIS weighted) image rendered by connecting the vertices to the camera and a set of photons.
struct LTWorkResult : public WorkResult {
    LTWorkResult(const Vector2i& size, const ReconstructionFilter* filter)
        : img(new ImageBlock(Bitmap::ESpectrum, size, filter))
        , img_useful(new ImageBlock(Bitmap::ESpectrum, size, filter))
        , range(new RangeWorkUnit)
    {}

    ref<RangeWorkUnit> range;

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

/// Samples a path from an emitter into the scene. Every vertex is connected to the sensor.
/// Hitting the sensor is only accepted if the last vertex was specular.
class LightPathTracer : public EmissionSampler {
public:
    LightPathTracer(EmitDistPtr, int rr_depth, int max_depth, Float pm_radius, int min_lp,
                    Scene*, Sampler*, Sensor*,
                    VertexMemPool*, LTWorkResult*, const lpm::MisComputer*);

    void handle_surface(int depth, int nullInteractions, bool delta,
                        const Intersection &its, const Medium *medium,
                        const Spectrum &emitted, const Spectrum& throughput) override;

    void handle_medium(int depth, int nullInteractions, bool delta,
                       const MediumSamplingRecord &mRec, const Medium *medium,
                       const Vector &wi, const Spectrum &weight) override;

protected:
    ref<LTWorkResult> result_;
    const lpm::MisComputer* mis_;
    Float radius_;
    int min_lp_;
};

} // namespace lpm

/// Traces the requested number of light paths. Renders an image and builds a set of photons.
class LightTracer : public WorkProcessor {
public:
    // TODO for this to work in a distributed system, the photons would have to be serialized

    LightTracer(const LPMConfig&, lpm::PhotonContainer*, lpm::EmitDistPtr);
    LightTracer(Stream*, InstanceManager*);

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
    ref<LTWorkResult>         result_;

    lpm::PhotonContainer* photons_;
    lpm::EmitDistPtr emit_dist_;
};

class LTProcess : public ParallelProcess {
public:
    LTProcess(const RenderJob *job, RenderQueue *queue, const LPMConfig& cfg, lpm::PhotonContainer*, lpm::EmitDistPtr);

    EStatus generateWork(WorkUnit*, int worker) override;
    ref<WorkProcessor> createWorkProcessor() const final;
    void processResult(const WorkResult*, bool cancelled) final;
    void bindResource(const std::string &name, int id);

    MTS_DECLARE_CLASS()

    ref<ImageBlock> img;
    ref<ImageBlock> img_useful;
    Float weight;

private:
    lpm::PhotonContainer* photons_;
    lpm::EmitDistPtr emit_dist_;

    ref<const RenderJob> job_;
    ref<RenderQueue> queue_;
    size_t complete_; ///< number of finished paths.
    size_t generated_;
    size_t granularity_; ///< number of paths in each work unit
    ref<ProgressReporter> progress_;
    ref<Mutex> result_mutex_;

    ref<Film> film_;

    LPMConfig config_;

    /// Notifies the scheduler that the given number of paths just finished.
    /// \returns the total number of paths that are finished.
    size_t notify_completed(size_t count);
};

MTS_NAMESPACE_END

#endif // __LPM_LIGHT_TRACE_H