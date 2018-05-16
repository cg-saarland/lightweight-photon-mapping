#include <mitsuba/render/scene.h>
#include <mitsuba/render/sensor.h>

#include "mis.h"

MTS_NAMESPACE_BEGIN

namespace lpm {

void Path::append(VertexMemPool* pool) {
    last->succ = pool->alloc();
    last->succ->predec = last;
    last = last->succ;
    last->succ = nullptr;
    len++;
}

void Path::remove_last() {
    last = last->predec;
    last->succ = nullptr;
    len--;
}

PathSampler::PathSampler(int rr_depth, int max_depth, Scene* scene, Sampler* sampler, Sensor* sensor, VertexMemPool* pool)
: pool_(pool), sampler_(sampler), scene_(scene), sensor_(sensor)
, rr_depth_(rr_depth), max_depth_(max_depth)
{
}

void PathSampler::trace(PathOrigin org, const Point2i& px, MisComputer* mis) {
    origin_ = org;
    pixel   = Point2(px);

    // initialize the path
    if (!pool_) path_.first = nullptr;
    else {
        path_.first = pool_->alloc();
        path_.first->predec = nullptr;
        path_.last  = path_.first;
        path_.len   = 0;
    }

    // if required: sample a time from the sensor
    ref<Sensor> sensor = scene_->getSensor();
    auto t = sensor->needsTimeSample()
           ? sensor->sampleTime(sampler_->next1D())
           : sensor->getShutterOpen();

    // sample a primary ray from an emitter in the scene
    PositionSamplingRecord  pos_sample;
    DirectionSamplingRecord dir_sample;
    Spectrum emitted;
    Float emission_uvpdf = 1.0f;
    if (origin_ == FROM_PIXEL) {
        emitted = sample_sensor(t, pos_sample, dir_sample, Point2(px), sensor);
        sensor->getSamplePosition(pos_sample, dir_sample, pixel);
    } else {
        emitted = sample_emitter(t, pos_sample, dir_sample, emission_uvpdf);
        if (emitted.isZero() || emission_uvpdf == 0.0f) return;
    }

    // TODO determine this properly to allow point light sources and directional lights
    bool delta_emitter = false;

    if (pool_) {
        // initialize the first vertex along the path
        path_.first->pdf_forward        = 1.0f;
        path_.first->pdf_reverse        = 1.0f;
        path_.first->on_surface         = true;
        path_.first->pos                = pos_sample.p;
        path_.first->sub_path_len       = 0;
        path_.first->last_specular      = false;
        path_.first->geom_normal        = pos_sample.n;
        path_.first->specular           = delta_emitter;
        path_.first->weight             = Spectrum(0.0f);
        path_.first->useful_for_merging = false;
        path_.first->guiding_contrib    = Spectrum(0.0f);
        path_.first->throughput         = Spectrum(1.0f) / emission_uvpdf;
        SAssert(dir_sample.measure == ESolidAngle);

        // prepare the next vertex
        path_.append(pool_);
        Float cos_theta_o = std::abs(dot(pos_sample.n, dir_sample.d));
        path_.last->pdf_forward        = pos_sample.pdf * dir_sample.pdf * emission_uvpdf;
        path_.last->cos_from_prev      = cos_theta_o;
        path_.last->on_surface         = false;
        path_.last->sub_path_len       = 1;
        path_.first->last_specular     = false;
        path_.last->useful_for_merging = false;
        path_.last->guiding_contrib    = Spectrum(0.0f);
        path_.first->throughput        = Spectrum(1.0f) / emission_uvpdf;

        if (origin_ == FROM_EMITTER)
            mis->init_emitter(path_.last, origin_);
        else
            mis->init_sensor(path_.last, origin_);
    }

    Ray ray;
    ray.setTime(pos_sample.time);
    ray.setOrigin(pos_sample.p);
    ray.setDirection(normalize(dir_sample.d));

    auto emitter = static_cast<const Emitter*>(pos_sample.object);
    auto medium  = emitter->getMedium();

    // incrementally build the full path
    int depth = 1;  // number of vertices
    int null_interactions = 0; // number of medium-matched transitions
    Spectrum throughput(1.0f);
    Intersection isect;
    bool last_delta = false; // was the last BSDF specular?
    while (!throughput.isZero() && (depth <= max_depth_ || max_depth_ < 0)) {
        scene_->rayIntersectAll(ray, isect);

        auto mode = origin_ == FROM_EMITTER ? EImportance : ERadiance;

        MediumSamplingRecord medium_sample;
        if (medium && sample_medium(medium, ray, isect, medium_sample)) {
            // we sampled a point inside a medium
            throughput *= medium_sample.sigmaS * medium_sample.transmittance / medium_sample.pdfSuccess;

            handle_medium(depth, null_interactions, last_delta, medium_sample, medium, -ray.d, throughput * emitted);

            PhaseFunctionSamplingRecord phase_sample(medium_sample, -ray.d, mode);
            throughput *= sample_phase_fn(medium, phase_sample);

            ray = Ray(medium_sample.p, phase_sample.wo, ray.time);
            ray.mint = 0;

            // TODO add vertex inside medium to compute MIS and perform photon mapping!
            // TODO compute MIS for vertices inside the medium
            // TODO compute next event pdf inside the medium

            last_delta = false; // participating media cannot be described by a delta function
        } else if (isect.t != std::numeric_limits<Float>::infinity()) {
            // a surface was intersected
            if (medium) throughput *= medium_sample.transmittance / medium_sample.pdfFailure;

            Vector wi = -ray.d;
            Float wiDotGeoN = dot(isect.geoFrame.n, wi);
            auto bsdf = isect.getBSDF();

            if (pool_) {
                path_.last->cos_to_prev = std::abs(wiDotGeoN);
                path_.last->on_surface  = true;
                path_.last->pos         = isect.p;
                path_.last->dist_sqr    = (ray.o - isect.p).lengthSquared();
                path_.last->weight      = emitted * throughput;
                path_.last->geom_normal = isect.geoFrame.n;
                path_.last->specular    = !(bsdf->getType() & BSDF::ESmooth); // are all components specular?
                path_.last->throughput  = throughput / emission_uvpdf;

                if (depth == 1 && origin_ == FROM_EMITTER)
                    path_.first->pdf_reverse = pdf_next_event(emitter, pos_sample, isect.p, isect.geoFrame.n);

                mis->update_hit(path_.last, origin_);
            }

            handle_surface(depth, null_interactions, last_delta, isect, medium, emitted, throughput);

            BSDFSamplingRecord bsdf_sample(isect, sampler_, mode);
            Float pdf_forward, pdf_reverse;
            Spectrum bsdf_weight = sample_surface(bsdf, bsdf_sample, pdf_forward, pdf_reverse);

            if (bsdf_weight.isZero()) break;

            /* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
            auto wo = isect.toWorld(bsdf_sample.wo);
            auto woDotGeoN = dot(isect.geoFrame.n, wo);
            if (wiDotGeoN * Frame::cosTheta(bsdf_sample.wi) <= 0 ||
                woDotGeoN * Frame::cosTheta(bsdf_sample.wo) <= 0)
                break;

            /* Keep track of the weight, medium and relative
                refractive index along the path */
            throughput *= bsdf_weight;
            if (isect.isMediumTransition())
                medium = isect.getTargetMedium(woDotGeoN);

            if (bsdf_sample.sampledType & BSDF::ENull) {
                ++null_interactions;
                last_delta = true;
            } else
                last_delta = bsdf_sample.sampledType & BSDF::EDelta;

            // Add the vertex to the path, complete pdf computations of the predecessor
            if (pool_) {
                path_.last->pdf_reverse = pdf_reverse;

                path_.append(pool_);
                path_.last->pdf_forward        = pdf_forward;
                path_.last->cos_from_prev      = std::abs(woDotGeoN);
                path_.last->on_surface         = false;
                path_.last->sub_path_len       = path_.last->predec->sub_path_len + 1;
                path_.last->last_specular      = last_delta;
                path_.last->useful_for_merging = false;
                path_.last->guiding_contrib    = Spectrum(0.0f);

                mis->update_bounce(path_.last, origin_);
            }

            ray.setOrigin(isect.p);
            ray.setDirection(wo);
            ray.mint = Epsilon;
        } else {
            // neither a medium nor a surface was intersected
            handle_miss(depth, null_interactions, last_delta, throughput * emitted);
            break;
        }

        if (depth++ >= rr_depth_) {
            /* Russian roulette: try to keep path weights equal to one,
                Stop with at least some probability to avoid
                getting stuck (e.g. due to total internal reflection) */
            Float q = std::min(throughput.max(), (Float) 0.95f);
            if (sampler_->next1D() >= q) {
                break;
            }
            throughput /= q;
        }
    }
}

Float PathSampler::pdf_next_event(const Emitter* em, const PositionSamplingRecord& pos, const Point& from, const Vector& from_normal) {
    DirectSamplingRecord sample;
    sample.object  = em;
    sample.p       = pos.p;
    const auto d   = pos.p - from;
    sample.dist    = d.length();
    sample.d       = d / sample.dist;
    sample.measure = EArea;
    sample.uv      = pos.uv;
    sample.n       = pos.n;
    sample.ref     = from;
    sample.refN    = from_normal;

    return scene_->pdfEmitterDirect(sample);
}

Float PathSampler::pdf_emit(const Emitter* em, const PositionSamplingRecord& pos, const Vector& dir) {
    PositionSamplingRecord p = pos;
    p.object = em;
    auto pdf = scene_->pdfEmitterPosition(p);
    SAssert(pdf > 0.0f);
    DirectionSamplingRecord dir_sample;
    dir_sample.measure = ESolidAngle;
    dir_sample.d = dir;
    pdf *= em->pdfDirection(dir_sample, p);
    return pdf;
}

Spectrum PathSampler::sample_emitter(Float time, PositionSamplingRecord& pos, DirectionSamplingRecord& dir, Float& uvpdf) {
    uvpdf = 1.0f;

    auto power   = scene_->sampleEmitterPosition(pos, sampler_->next2D());
    auto emitter = static_cast<const Emitter*>(pos.object);
    auto medium  = emitter->getMedium();

    handle_emitter(medium, pos, power);

    auto diruv = emitter->needsDirectionSample() ? sampler_->next2D() : Point2(0.5f);
    return power * emitter->sampleDirection(dir, pos, diruv);
}

Spectrum PathSampler::sample_sensor(Float time, PositionSamplingRecord& pos_sample, DirectionSamplingRecord& dir_sample,
                                    const Point2& pixel_pos, Sensor* sensor) {
    Point2 pixel_sample = sampler_->next2D();
    Point2 aperture_sample = sensor->needsApertureSample() ? sampler_->next2D() : Point2(0.5f);

    pos_sample = PositionSamplingRecord(time);
    Spectrum importance = scene_->sampleSensorPosition(pos_sample,
		(sensor->getType() & Sensor::EPositionSampleMapsToPixels) ? pixel_sample : aperture_sample,
        &pixel_pos);

	importance *= sensor->sampleDirection(dir_sample, pos_sample,
		(sensor->getType() & Sensor::EPositionSampleMapsToPixels) ? aperture_sample : pixel_sample,
        &pixel_pos);

    return importance;
}

bool PathSampler::sample_medium(const Medium* medium, Ray& ray, Intersection& isect, MediumSamplingRecord& medium_sample) {
    return medium->sampleDistance(Ray(ray, 0, isect.t), medium_sample, sampler_);
}

Float PathSampler::sample_phase_fn(const Medium* medium, PhaseFunctionSamplingRecord& phase_sample) {
    return medium->getPhaseFunction()->sample(phase_sample, sampler_);
}

Spectrum PathSampler::sample_surface(const BSDF* bsdf, BSDFSamplingRecord& sample, Float& pdf_forward, Float& pdf_reverse) {
    auto val = bsdf->sample(sample, pdf_forward, sampler_->next2D());

    auto rev_sample = sample;
    rev_sample.reverse();
    pdf_reverse = bsdf->pdf(rev_sample);

    return val;
}

} // namespace lpm

MTS_NAMESPACE_END