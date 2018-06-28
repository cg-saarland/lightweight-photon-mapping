#include <mitsuba/core/statistics.h>
#include "cam_trace.h"
#include "mis.h"

MTS_NAMESPACE_BEGIN

MTS_IMPLEMENT_CLASS(CTWorkResult, false, ImageBlock);
MTS_IMPLEMENT_CLASS_S(CameraTracer, false, WorkProcessor);
MTS_IMPLEMENT_CLASS(CTProcess, false, BlockedRenderProcess);

namespace lpm {

CamPathTracer::CamPathTracer(EmitDistPtr d, int rr_depth, int max_depth, Scene* scene, Sampler* sampler, Sensor* sensor,
                             VertexMemPool* pool, CTWorkResult* res, const lpm::MisComputer* mis,
                             const lpm::PhotonContainer* p, const lpm::HashGrid* h, Float radius, int lp_count, int min_lp, bool visp)
: EmissionSampler(d, rr_depth, max_depth, scene, sampler, sensor, pool)
, contrib(0.0f), contrib_useful(0.0f), result_(res), mis_(mis)
, photons_(p), photon_accel_(h), radius_(radius), radius_sqr_(radius * radius), lp_count_(lp_count), min_lp_(min_lp)
, vis_photons_(visp)
{
}

void CamPathTracer::handle_surface(int depth, int nullInteractions, bool delta,
                                   const Intersection &its, const Medium *medium,
                                   const Spectrum &emitted, const Spectrum& throughput)
{
    PathSampler::handle_surface(depth, nullInteractions, delta, its, medium, emitted, throughput);

    auto weight = emitted * throughput;

    const BSDF *bsdf = its.getBSDF();

    if (its.isEmitter() && !vis_photons_) {
        hit(depth, bsdf, its, weight);
        return;
    }

    // TODO do not perform merging if it is disabled!
    merge(depth, bsdf, its, weight);
    if (vis_photons_) return;

    nee(depth, bsdf, its, weight);
}

void CamPathTracer::handle_medium(int depth, int nullInteractions, bool delta,
                                  const MediumSamplingRecord &mRec, const Medium *medium,
                                  const Vector &wi, const Spectrum &weight)
{
    PathSampler::handle_medium(depth, nullInteractions, delta, mRec, medium, wi, weight);
}

void CamPathTracer::handle_miss(int depth, int nullInteractions, bool delta, const Spectrum& weight) {
    // TODO add support for environment maps
}

void CamPathTracer::nee(int depth, const BSDF* bsdf, const Intersection& its, const Spectrum& throughput) {
    if (depth >= max_depth_ && max_depth_ > 0)
        return;

    DirectSamplingRecord dRec(its);
    auto value = throughput * scene_->sampleEmitterDirect(dRec, sampler_->next2D(), true);

    BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
    value *= bsdf->eval(bRec);

    if (value.isZero()) return;

    Vector wi = its.toWorld(its.wi);
    auto wiDotGeoN = dot(its.geoFrame.n, wi);
    auto woDotGeoN = dot(its.geoFrame.n, dRec.d);
    // if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
    //     woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
    //     return;

    auto rev_bsdf_sample = bRec;
    rev_bsdf_sample.reverse();

    auto emitter     = static_cast<const Emitter*>(dRec.object);
    auto pdf_em      = pdf_emit(emitter, dRec, -dRec.d);
    auto pdf_nee     = dRec.pdf;
    auto cos_surf    = std::abs(woDotGeoN);
    auto cos_emitter = std::abs(dot(dRec.n, -dRec.d));
    auto pdf_hit     = bsdf->pdf(bRec);
    auto pdf_reverse = bsdf->pdf(rev_bsdf_sample);
    auto mis_weight  = mis_->weight_pt_nee(path_.last, cos_surf, cos_emitter, pdf_nee, pdf_hit, pdf_em, pdf_reverse);

    contrib += mis_weight * value;
}

void CamPathTracer::hit(int depth, const BSDF* bsdf, const Intersection& its, const Spectrum& throughput) {
    Vector wi = its.toWorld(its.wi);
    auto emitter = its.shape->getEmitter();
    auto radiance = emitter->eval(its, wi);

    if (radiance.isZero()) return;

    // Compute pdfs and MIS weight
    auto pos_prev    = path_.last->predec->pos;
    auto normal_prev = path_.last->predec->geom_normal;

    PositionSamplingRecord pos_sample(its);
    auto pdf_em  = pdf_emit(emitter, pos_sample, wi);
    auto pdf_nee = pdf_next_event(emitter, pos_sample, pos_prev, normal_prev);

    auto mis_weight = mis_->weight_pt_hit(path_.last, pdf_em, pdf_nee);

    contrib += mis_weight * radiance * throughput;
}

void CamPathTracer::merge(int depth, const BSDF* bsdf, const Intersection& its, const Spectrum& throughput) {
    // compute a contribution with density estimation
    Spectrum photon_contrib(0.0f);

    constexpr uint32_t k = 10; // TODO make configurable
    auto photons = STACK_ARRAY(Photon*, k);
    auto count = photon_accel_->query(its.p, photons, k);
    const Float r2 = (count == k) ? (photons[k - 1]->pos - its.p).lengthSquared() : radius_sqr_;
    const Float ep_kernel_factor = 2.0f / (r2 * M_PI);
    for (int i = 0; i < count; ++i) {
        auto& p = *(photons[i]);

        if (p.depth + depth > max_depth_ && max_depth_ > 0) continue;

        BSDFSamplingRecord bsdf_sample(its, its.toLocal(p.dir));
        auto bsdf_value = bsdf->eval(bsdf_sample);

        if (bsdf_value.isZero() && !vis_photons_) continue;

        // compute the MIS weight
        auto pdf_rev_light = bsdf->pdf(bsdf_sample);
        auto rev_sample    = bsdf_sample; rev_sample.reverse();
        auto pdf_rev_cam   = bsdf->pdf(rev_sample);
        auto mis_weight    = mis_->weight_merge(path_.last, &p, pdf_rev_cam, pdf_rev_light);

        const float d2 = (p.pos - its.p).lengthSquared();
        const float kernel = (1.0f - d2 / r2) * ep_kernel_factor;

        auto adjoint_correction = 1.0f / fabsf(Frame::cosTheta(its.toLocal(p.dir)));
        Spectrum value = mis_weight * kernel * throughput;
        photon_contrib += value * (vis_photons_ ? Spectrum(1.0f) : bsdf_value * p.weight * adjoint_correction);

        if (p.is_useful) contrib_useful += value * bsdf_value * p.weight * adjoint_correction * (1.0f / lp_count_);

        // Log the luminance of the photon's contribution to the image
        // Do not divide by the number of light paths because the light tracer doesn't
        Spectrum guiding_weight = value * bsdf_value * p.throughput;
        atomic_add(p.guiding_contrib, guiding_weight.getLuminance());
    }

    // Divide by number of light paths (each light path is a Monte Carlo estimate for the entire contribution of that light to the scene)
    photon_contrib *= 1.0f / lp_count_;

    if (!vis_photons_ || depth == 1)
        contrib += photon_contrib;
    else if (vis_photons_) {
        bool all_specular = true;
        for (auto v = path_.first->succ->succ; v; v = v->succ)
            if (!v->last_specular) { all_specular = false; break; }
        if (all_specular)
            contrib += photon_contrib;
    }
}

} // namespace lpm

CameraTracer::CameraTracer(const LPMConfig& cfg, const lpm::PhotonContainer* p, const lpm::HashGrid* h, lpm::EmitDistPtr ed)
: config_(cfg), emit_dist_(ed), photons_(p), photon_accel_(h)
{
}

CameraTracer::CameraTracer(Stream* stream, InstanceManager* manager)
: config_(stream)
{
}

void CameraTracer::serialize(Stream* stream, InstanceManager* manager) const {
    config_.serialize(stream);
}

ref<WorkResult> CameraTracer::createWorkResult() const {
    // TODO adapt block size
    auto blockSize = Vector2i(64);
    auto block = new CTWorkResult(blockSize, filter_.get());
    block->img->setOffset(Point2i(0, 0));
    block->img->setSize(blockSize);
    block->img_useful->setOffset(Point2i(0, 0));
    block->img_useful->setSize(blockSize);
    return block;
}

ref<WorkProcessor> CameraTracer::clone() const {
    return new CameraTracer(config_, photons_, photon_accel_, emit_dist_);
}

ref<WorkUnit> CameraTracer::createWorkUnit() const {
    return new RectangularWorkUnit();
}

void CameraTracer::prepare() {
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

void CameraTracer::process(const WorkUnit* wu, WorkResult* res, const bool& stop) {
    auto rect  = static_cast<const RectangularWorkUnit*>(wu);
    result_ = static_cast<CTWorkResult*>(res);
    result_->img->setOffset(rect->getOffset());
    result_->img->setSize(rect->getSize());;
    result_->img->clear();
    result_->img_useful->setOffset(rect->getOffset());
    result_->img_useful->setSize(rect->getSize());;
    result_->img_useful->clear();

    lpm::VmMis mis(config_.lp_per_cp, config_.pm_radius, config_.lp_count,
                   config_.light_trace_di, config_.enable_merge, config_.merge_primary, config_.merge_di);
    lpm::VertexMemPool pool;

    hilbert_curve_.initialize(TVector2<uint8_t>(rect->getSize()));
    for (size_t i=0; i<hilbert_curve_.getPointCount(); ++i) {
        Point2i pixel = Point2i(hilbert_curve_[i]) + Vector2i(rect->getOffset());
        sampler_->generate(pixel);

        for (size_t j = 0; j < config_.cp_per_iter; j++) {
            if (stop) break;

            auto cp_tracer = lpm::CamPathTracer(emit_dist_, config_.rr_depth, config_.max_depth, scene_, sampler_, sensor_,
                                                &pool, result_, &mis, photons_, photon_accel_, config_.pm_radius,
                                                config_.lp_count, config_.min_lp, config_.vis_photons);
            cp_tracer.trace(lpm::PathSampler::FROM_PIXEL, pixel, &mis);
            result_->img->put(cp_tracer.pixel, cp_tracer.contrib, 1.0f);
            result_->img_useful->put(cp_tracer.pixel, cp_tracer.contrib_useful, 1.0f);
            pool.reset();
        }
    }

    result_ = nullptr;
}

CTProcess::CTProcess(const RenderJob *job, RenderQueue *queue, const LPMConfig& cfg,
                     const lpm::PhotonContainer* p, const lpm::HashGrid* h, lpm::EmitDistPtr ed)
// TODO configurable block size
: BlockedRenderProcess(job, queue, 64), refresh_timer_(new Timer), config_(cfg)
, emit_dist_(ed), photons_(p), photon_accel_(h)
{
}

void CTProcess::processResult(const WorkResult *wr, bool cancelled) {
    if (cancelled) return;

    auto result = static_cast<const CTWorkResult*>(wr);
    LockGuard lock(m_resultMutex);
    m_progress->update(++m_resultCount);

    img->put(static_cast<const ImageBlock*>(result->img.get()));
    img_useful->put(static_cast<const ImageBlock*>(result->img_useful.get()));

    m_queue->signalWorkEnd(m_parent, result->img.get(), false);
}

void CTProcess::bindResource(const std::string &name, int id) {
    BlockedRenderProcess::bindResource(name, id);
    if (name == "sensor") {
        auto sensor = static_cast<Sensor*>(Scheduler::getInstance()->getResource(id));
        film_  = sensor->getFilm();
        img = new ImageBlock(Bitmap::ESpectrum, sensor->getFilm()->getCropSize(), nullptr);
        img->clear();

        img_useful = new ImageBlock(Bitmap::ESpectrum, sensor->getFilm()->getCropSize(), nullptr);
        img_useful->clear();
    }
}

ref<WorkProcessor> CTProcess::createWorkProcessor() const {
    return new CameraTracer(config_, photons_, photon_accel_, emit_dist_);
}


MTS_NAMESPACE_END