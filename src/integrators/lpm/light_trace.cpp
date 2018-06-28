#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/core/statistics.h>
#include "light_trace.h"
#include "mis.h"
#include "hash_grid.h"
#include "usefulness.h"

MTS_NAMESPACE_BEGIN

MTS_IMPLEMENT_CLASS(LTWorkResult, false, WorkResult);
MTS_IMPLEMENT_CLASS_S(LightTracer, false, WorkProcessor);
MTS_IMPLEMENT_CLASS(LTProcess, false, ParallelProcess);

namespace lpm {

LightPathTracer::LightPathTracer(EmitDistPtr dist_ptr, int rr_depth, int max_depth, Float pm_radius, int min_lp,
                                 Scene* scene, Sampler* sampler, Sensor* sensor,
                                 VertexMemPool* pool, LTWorkResult* res, const lpm::MisComputer* mis)
: EmissionSampler(dist_ptr, rr_depth, max_depth, scene, sampler, sensor, pool)
, result_(res), mis_(mis), radius_(pm_radius), min_lp_(min_lp)
{}

void LightPathTracer::handle_surface(int depth, int nullInteractions, bool delta,
                                     const Intersection &its, const Medium *medium,
                                     const Spectrum &emitted, const Spectrum& throughput)
{
    if (its.isSensor())
        return;

    if (depth >= max_depth_ && max_depth_ > 0) return;

    path_.last->useful_for_merging = is_useful_photon(&path_, min_lp_, radius_);

    DirectSamplingRecord dRec(its);
    Spectrum value = scene_->sampleSensorDirect(dRec, sampler_->next2D(), true);

    if (value.isZero()) return;

    SAssert(!value.isNaN());

    const BSDF *bsdf = its.getBSDF();

    Vector wo = dRec.d;
    BSDFSamplingRecord bRec(its, its.toLocal(wo), EImportance);

    /* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
    Vector wi = its.toWorld(its.wi);
    Float wiDotGeoN = dot(its.geoFrame.n, wi),
          woDotGeoN = dot(its.geoFrame.n, wo);
    // if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
    //     woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
    //     return;

    /* Adjoint BSDF for shading normals -- [Veach, p. 155] */
    Float correction = std::abs(
        (Frame::cosTheta(bRec.wi) * woDotGeoN)/
        (Frame::cosTheta(bRec.wo) * wiDotGeoN));
    value *= bsdf->eval(bRec) * correction;

    SAssert(!value.isNaN());

    auto rev_bsdf_sample = bRec;
    rev_bsdf_sample.reverse();

    DirectionSamplingRecord cam_dir;
    cam_dir.d = -wo;
    cam_dir.measure = ESolidAngle;
    PositionSamplingRecord cam_pos = dRec;

    // compute pdfs and MIS weight
    auto pdf_camera = sensor_->pdfDirection(cam_dir, cam_pos);
    pdf_camera *= sensor_->pdfPosition(cam_pos);

    SAssert(cam_dir.measure == ESolidAngle);
    auto cos_to_cam = std::abs(dot(its.geoFrame.n, wo));
    pdf_camera *= cos_to_cam / (dRec.dist * dRec.dist);
    auto pdf_reverse = bsdf->pdf(rev_bsdf_sample);
    auto mis_weight = mis_->weight_lt_nee(path_.last, pdf_camera, pdf_reverse);

    value *= mis_weight;

    // log the image contribution for use with guiding methods
    path_.last->guiding_contrib = value * path_.last->throughput;

    value *= emitted * throughput;

    // Splat onto the accumulation buffer
    result_->img->put(dRec.uv, (Float *) &value[0]);

    if (path_.last->useful_for_merging)
        result_->img_useful->put(dRec.uv, (Float *) &value[0]);
}

void LightPathTracer::handle_medium(int depth, int nullInteractions, bool delta,
                                    const MediumSamplingRecord &mRec, const Medium *medium,
                                    const Vector &wi, const Spectrum &weight)
{
}

} // namespace lpm

/////////////////////////////////////////////////////////////////////
// Light Tracer
/////////////////////////////////////////////////////////////////////

LightTracer::LightTracer(const LPMConfig& cfg, lpm::PhotonContainer* p, lpm::EmitDistPtr ed)
: config_(cfg), photons_(p), emit_dist_(ed)
{
}

LightTracer::LightTracer(Stream* stream, InstanceManager* manager)
: config_(stream)
{
}

void LightTracer::serialize(Stream* stream, InstanceManager* manager) const {
    config_.serialize(stream);
}

ref<WorkResult> LightTracer::createWorkResult() const {
    return new LTWorkResult(config_.img_size, filter_.get());
}

ref<WorkProcessor> LightTracer::clone() const {
    return new LightTracer(config_, photons_, emit_dist_);
}

ref<WorkUnit> LightTracer::createWorkUnit() const {
    return new RangeWorkUnit();
}

void LightTracer::prepare() {
    Scene* scene = static_cast<Scene *>(getResource("scene"));
    scene_   = new Scene(scene);
    sampler_ = static_cast<Sampler *>(getResource("sampler"));
    sensor_  = static_cast<Sensor *>(getResource("sensor"));
    scene_->removeSensor(scene->getSensor());
    scene_->addSensor(sensor_);
    scene_->setSensor(sensor_);
    scene_->initializeBidirectional();
    filter_ = sensor_->getFilm()->getReconstructionFilter();
}

void LightTracer::process(const WorkUnit* wu, WorkResult* res, const bool& stop) {
    auto range  = static_cast<const RangeWorkUnit*>(wu);
    result_ = static_cast<LTWorkResult*>(res);
    result_->range->set(range);
    result_->img->clear();
    result_->img_useful->clear();

    // generate some samples
    sampler_->generate(Point2i(0));

    lpm::VmMis mis(config_.lp_per_cp, config_.pm_radius, config_.lp_count,
                   config_.light_trace_di, config_.enable_merge, config_.merge_primary, config_.merge_di);
    lpm::VertexMemPool pool;

    auto lp_tracer = lpm::LightPathTracer(emit_dist_, config_.rr_depth, config_.max_depth, config_.pm_radius, config_.min_lp,
                                          scene_, sampler_, sensor_, &pool, result_, &mis);

    for (size_t index = range->getRangeStart(); index <= range->getRangeEnd() && !stop; ++index) {
        pool.reset();
        lp_tracer.trace(lpm::PathSampler::FROM_EMITTER, Point2i(0), &mis);

        // TODO photon generation could also be done on the fly to improve memory performance

        // generate photons along the path
        // auto path_photons = STACK_ARRAY(lpm::Photon, lp_tracer.path().len);
        std::vector<lpm::Photon> cont(lp_tracer.path().len); // TODO reuse / reduce number of allocations
        auto path_photons = cont.data();

        if (lp_tracer.path().len < 1) continue; // never intersected anything

        uint32_t generated = 0;
        uint32_t depth     = 1;
        for (auto vert = lp_tracer.path().first->succ; vert; vert = vert->succ, ++depth) {
            if (!vert->specular && vert->on_surface) {
                auto& new_photon = path_photons[generated++];
                new_photon.ancestor_idx    = generated - 1;
                new_photon.pos             = vert->pos;
                new_photon.dir             = normalize(vert->predec->pos - vert->pos);
                new_photon.weight          = vert->weight;
                new_photon.depth           = depth;
                new_photon.is_useful       = vert->useful_for_merging;
                new_photon.emission_sample = lp_tracer.emission_sample;
                new_photon.guiding_contrib = vert->guiding_contrib.getLuminance();
                new_photon.throughput      = vert->throughput.getLuminance();

                // multiply MIS partials with geometry terms
                new_photon.partial_mis_sum = vert->partial_mis_sum / vert->cos_to_prev;

                // include weight for path tracer next event estimation if not yet part of partial
                if (vert->sub_path_len == 1) {
                    Float pdf_next_event_surf_area = vert->predec->pdf_reverse;
                    Float pdf_emission_surf_area   = vert->pdf_forward * vert->cos_to_prev / vert->dist_sqr;
                    new_photon.partial_mis_sum += pdf_next_event_surf_area / pdf_emission_surf_area;
                }
            }
        }

        photons_->add(path_photons, path_photons + generated);
    }

    result_ = nullptr;
}

/////////////////////////////////////////////////////////////////////
// LT Work Processor
/////////////////////////////////////////////////////////////////////

LTProcess::LTProcess(const RenderJob *job, RenderQueue *queue, const LPMConfig& cfg, lpm::PhotonContainer* p, lpm::EmitDistPtr emit_dist)
: photons_(p), emit_dist_(emit_dist), job_(job), queue_(queue)
, complete_(0), generated_(0), config_(cfg)
{
    // Choose a suitable work unit granularity if none was specified
    granularity_ = std::max((size_t) 1, cfg.lp_count / (16 * Scheduler::getInstance()->getWorkerCount()));

    //  Create a visual progress reporter
    progress_ = new ProgressReporter("Tracing light paths", cfg.lp_count, job);
    result_mutex_ = new Mutex();
}

ref<WorkProcessor> LTProcess::createWorkProcessor() const {
    return new LightTracer(config_, photons_, emit_dist_);
}

void LTProcess::bindResource(const std::string &name, int id) {
    if (name == "sensor") {
        auto sensor = static_cast<Sensor*>(Scheduler::getInstance()->getResource(id));
        film_  = sensor->getFilm();
        img = new ImageBlock(Bitmap::ESpectrum, sensor->getFilm()->getCropSize(), nullptr);
        img->clear();

        img_useful = new ImageBlock(Bitmap::ESpectrum, sensor->getFilm()->getCropSize(), nullptr);
        img_useful->clear();
    }
    ParallelProcess::bindResource(name, id);
}

LTProcess::EStatus LTProcess::generateWork(WorkUnit* unit, int worker) {
    if (generated_ == static_cast<size_t>(config_.lp_count))
        return EFailure; // There is no more work

    auto work_size = std::min(granularity_, config_.lp_count - generated_);

    auto range = static_cast<RangeWorkUnit *>(unit);
    range->setRange(generated_, generated_ + work_size - 1);
    generated_ += work_size;

    return ESuccess;
}

void LTProcess::processResult(const WorkResult* wr, bool cancelled) {
    if (cancelled) return;

    auto result = static_cast<const LTWorkResult*>(wr);
    auto range = result->range.get();

    LockGuard lock(result_mutex_);
    complete_ += range->getSize();
    progress_->update(complete_);
    img->put(result->img.get());
    img_useful->put(result->img_useful.get());

    if (complete_ == static_cast<size_t>(config_.lp_count)) {
        weight = (img->getWidth() * img->getHeight()) / (Float)config_.lp_count;
        queue_->signalRefresh(job_);
    }
}

MTS_NAMESPACE_END