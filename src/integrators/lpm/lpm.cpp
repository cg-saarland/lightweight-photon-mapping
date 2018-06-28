#include <mitsuba/render/scene.h>

#include "lpm.h"
#include "light_trace.h"
#include "cam_trace.h"
#include "pixel_classification.h"

MTS_NAMESPACE_BEGIN

MTS_IMPLEMENT_CLASS_S(LPMIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(LPMIntegrator, "Lightweight photon mapping integrator");

LPMIntegrator::LPMIntegrator(const Properties &props) : Integrator(props) {
    config_.max_depth             = props.getInteger("maxDepth"        ,    -1);
    config_.rr_depth              = props.getInteger("rrDepth"         ,     5);
    config_.lp_per_cp             = props.getFloat  ("lpRatio"         ,  1.0f);

    config_.initial_radius        = props.getFloat  ("initialRadius"   ,  0.0f);
    config_.iterations            = props.getInteger("iterations"      ,     1);
    config_.cp_per_iter           = props.getInteger("sppPerIter"      ,     1);

    config_.vis_photons           = props.getBoolean("visPhotons"      , false);
    config_.guiding_start         = props.getInteger("guideStart"      ,     1);
    config_.guiding_useful_start  = props.getInteger("guideUsefulStart",     2);
    config_.min_lp                = props.getInteger("minLp"           ,  5000);

    config_.enable_merge          = props.getBoolean("enableMerge"     , true);
    config_.merge_primary         = props.getBoolean("mergePrimary"    , false);
    config_.light_trace_di        = props.getBoolean("lightTraceDi"    , false);
    config_.merge_di              = props.getBoolean("mergeDi"         , false);

    config_.vis_lpcontrib = props.getBoolean("visLpContrib", false);

    if (config_.rr_depth <= 0)
        Log(EError, "'rr_depth' must be set to a value greater than zero!");

    if (config_.max_depth <= 0 && config_.max_depth != -1)
        Log(EError, "'max_depth' must be set to -1 (infinite) or a value greater than zero!");
}

LPMIntegrator::LPMIntegrator(Stream *stream, InstanceManager *manager)
    : Integrator(stream, manager) {
    config_ = LPMConfig(stream);
}

void LPMIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
    Integrator::serialize(stream, manager);
    config_.serialize(stream);
}

bool LPMIntegrator::preprocess(const Scene *scene, RenderQueue *queue,
                               const RenderJob *job, int sceneResID, int sensorResID,
                               int samplerResID)
{
    Integrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);

    canceled_ = false;

    if (scene->getSampler()->getClass()->getName() != "IndependentSampler")
        Log(EError, "The LPM integrator only works with an independent sampler.");

    if (scene->getSubsurfaceIntegrators().size() > 0)
        Log(EError, "Subsurface integrators are not supported "
            "by the lightweight photon mapper!");

    if (config_.initial_radius == 0) {
        /* Guess an initial radius if not provided
            (scene width / horizontal or vertical pixel count) * 5 */
        Float rad = scene->getBSphere().radius;
        Vector2i filmSize = scene->getSensor()->getFilm()->getSize();

        config_.initial_radius = std::min(rad / filmSize.x, rad / filmSize.y) * 5;
    }

    return true;
}

bool LPMIntegrator::render(Scene *scene, RenderQueue *queue, const RenderJob *job,
    int sceneResID, int sensorResID, int samplerResID)
{
    lpm::EmitDistPtr emit_dist(new lpm::EmissionDistribution(scene, config_.guiding_start, config_.guiding_useful_start));

    ref<Sensor> sensor = scene->getSensor();
    Film *film         = sensor->getFilm();
    config_.img_size   = film->getCropSize();
    config_.lp_count   = config_.img_size.x * config_.img_size.y * config_.lp_per_cp;
    config_.pm_radius  = config_.initial_radius;

    ref<Bitmap> img = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize());
    img->clear();

    ref<Bitmap> img_useful = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getCropSize());
    img_useful->clear();

    ref<Timer> timer_build = new Timer(false);
    ref<Timer> timer_norm  = new Timer(false);

    int num_discarded = 0;
    for (int i = 0; i < config_.iterations && !canceled_; ++i) {
        SLog(EInfo, "rendering iteration %i out of %i", i + 1, config_.iterations);

        // trace light paths, store vertices, connect to the camera
        ref<LTProcess> lt_proc = new LTProcess(job, queue, config_, &photons_, emit_dist);
        ref<Scheduler> sched = Scheduler::getInstance();
        lt_proc->bindResource("scene"  , sceneResID);
        lt_proc->bindResource("sensor" , sensorResID);
        lt_proc->bindResource("sampler", samplerResID);
        sched->schedule(lt_proc);
        sched->wait(lt_proc);

        // build the photon map
        photon_accel_.build(photons_.begin(), photons_.end(), config_.pm_radius);

        SLog(EInfo, "%i photons", photons_.size());

        // trace camera paths, perform next event, perform density estimation
        ref<CTProcess> ct_proc = new CTProcess(job, queue, config_, &photons_, &photon_accel_, emit_dist);
        ct_proc->bindResource("scene"  , sceneResID);
        ct_proc->bindResource("sensor" , sensorResID);
        ct_proc->bindResource("sampler", samplerResID);
        sched->schedule(ct_proc);
        sched->wait(ct_proc);

        // update the guiding distribution with the weights (stored inside the photons)
        // TODO this could be done much faster by sorting and parallelizing over the light sources
        //      also, there is no need to refine the histograms forever
        timer_build->reset();
        for (const auto& p : photons_) {
            emit_dist->add_record(p.emission_sample, p.guiding_contrib, p.is_useful);
        }
        auto secs = timer_build->stop();
        SLog(EInfo, "updating distribution done after %f seconds", secs);

        // TODO this only needs to be done when / if the histograms are acutally used
        timer_norm->reset();
        emit_dist->normalize();
        secs = timer_norm->stop();
        SLog(EInfo, "normalizing emission distribution done after %f seconds", secs);

        photons_.reset();

        if (i == config_.guiding_useful_start || i == config_.guiding_start) {
            // reset the image
            img->clear();
            img_useful->clear();
            num_discarded = i;
        }

        // combine the images of this iteration with the one of the previous iterations
        auto lp_img = lt_proc->img->getBitmap()->convert(Bitmap::ESpectrum, Bitmap::EFloat, 1.0f, lt_proc->weight);
        auto lp_img_useful = lt_proc->img_useful->getBitmap()->convert(Bitmap::ESpectrum, Bitmap::EFloat, 1.0f, lt_proc->weight);
        img->accumulate(lp_img);
        img_useful->accumulate(lp_img_useful);

        auto cp_img = ct_proc->img->getBitmap()->convert(Bitmap::ESpectrum, Bitmap::EFloat, 1.0f, 1.0f);
        auto cp_img_useful = ct_proc->img_useful->getBitmap()->convert(Bitmap::ESpectrum, Bitmap::EFloat, 1.0f, 1.0f);
        img->accumulate(cp_img);
        img_useful->accumulate(cp_img_useful);

        if (!config_.vis_photons)
            film->setBitmap(img, 1.0f / (i - num_discarded + 1));
        else
            film->setBitmap(cp_img);

        if (config_.vis_lpcontrib)
            film->setBitmap(img_useful, 1.0f / (i - num_discarded + 1));

        // After the initial training is over, adapt the number of light paths
        if (i >= config_.guiding_useful_start) {
            config_.lp_count = lpm::num_pixels_requiring_light_paths(img_useful.get(), img.get()) * config_.lp_per_cp;
            SLog(EInfo, "Reducing number of light paths to: %d", config_.lp_count);
        }
    }

    SLog(EInfo, "writing emission distributions to file ... ");
    emit_dist->dump("emdist.hist", true);
    emit_dist->dump("emdist_contrib.hist", false);

    return true;
}

void LPMIntegrator::cancel() {
    Scheduler::getInstance()->cancel(process_);
    canceled_ = true;
}

MTS_NAMESPACE_END