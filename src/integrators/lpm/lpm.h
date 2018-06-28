#ifndef __LPM_H
#define __LPM_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/integrator.h>
#include "hash_grid.h"
#include "emission_sampler.h"

MTS_NAMESPACE_BEGIN

/// Stores the parameters of the lightweight photon mapper
struct LPMConfig {
    int max_depth;
    int rr_depth;
    Float lp_per_cp;
    int lp_count;
    Vector2i img_size;

    Float initial_radius;
    Float pm_radius;
    int iterations;
    int cp_per_iter;

    bool vis_photons;
    int guiding_start;
    int guiding_useful_start;
    int min_lp;

    bool enable_merge;
    bool merge_primary;
    bool light_trace_di;
    bool merge_di;

    bool vis_lpcontrib;

    inline LPMConfig() { }

    inline LPMConfig(Stream *stream) {
        max_depth            = stream->readInt();
        rr_depth             = stream->readInt();
        lp_per_cp            = stream->readFloat();
        lp_count             = stream->readInt();
        img_size             = Vector2i(stream);

        initial_radius       = stream->readFloat();
        pm_radius            = stream->readFloat();
        iterations           = stream->readInt();
        cp_per_iter          = stream->readInt();

        vis_photons          = stream->readBool();
        guiding_start        = stream->readInt();
        guiding_useful_start = stream->readInt();
        min_lp               = stream->readInt();

        enable_merge         = stream->readBool();
        merge_primary        = stream->readBool();
        light_trace_di       = stream->readBool();
        merge_di             = stream->readBool();

        vis_lpcontrib = stream->readBool();
    }

    inline void serialize(Stream *stream) const {
        stream->writeInt(max_depth);
        stream->writeInt(rr_depth);
        stream->writeFloat(lp_per_cp);
        stream->writeInt(lp_count);
        img_size.serialize(stream);

        stream->writeFloat(initial_radius);
        stream->writeFloat(pm_radius);
        stream->writeInt(iterations);
        stream->writeInt(cp_per_iter);

        stream->writeBool(vis_photons);
        stream->writeInt(guiding_start);
        stream->writeInt(guiding_useful_start);
        stream->writeInt(min_lp);

        stream->writeBool(enable_merge);
        stream->writeBool(merge_primary);
        stream->writeBool(light_trace_di);
        stream->writeBool(merge_di);

        stream->writeBool(vis_lpcontrib);
    }

    void dump() const {
        SLog(EDebug, "Bidirectional path tracer configuration:");
        SLog(EDebug, "   Maximum path depth            : %i", max_depth);
        SLog(EDebug, "   Russian roulette depth        : %i", rr_depth);
        SLog(EDebug, "   Light path ratio              : %f", lp_per_cp);
        SLog(EDebug, "   Number of light paths         : %i", lp_count);
        SLog(EDebug, "   Image size                    : %ix%i", img_size.x, img_size.y);
        SLog(EDebug, "   Initial radius                : %f", initial_radius);
        SLog(EDebug, "   Current photon radius         : %f", pm_radius);
        SLog(EDebug, "   Minimum number of light paths : %i", min_lp);
        SLog(EDebug, "   Number of iterations          : %i", iterations);
        SLog(EDebug, "   Camera paths per iteration    : %i", cp_per_iter);
        SLog(EDebug, "   Visualize photons only?       : %i", vis_photons);
        // TODO
    }
};

/*!\plugin{lpm}{Lightweight photon mapping}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *	       A value of \code{1} will only render directly visible light sources.
 *	       \code{2} will lead to single-bounce (direct-only) illumination,
 *	       and so on. \default{\code{-1}}
 *	   }
 *	   \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *	      which the implementation will start to use the ``russian roulette''
 *	      path termination criterion. \default{\code{5}}
 *	   }
 *     \parameter{lpCount}{\Integer}{Specifies the (intial) number of light paths
 *        to trace during each iteration.
 *     }
 * }
 *
 * This plugin implements lightweight photon mapping (short: LPM).
 *
 * \remarks{
 *    \item This integrator does not work with dipole-style subsurface
 *    scattering models.
 * }
 */
class LPMIntegrator : public Integrator {
public:
    LPMIntegrator(const Properties &props);

    /// Unserialize from a binary data stream
    LPMIntegrator(Stream *stream, InstanceManager *manager);

    void serialize(Stream *stream, InstanceManager *manager) const final;

    bool preprocess(const Scene *scene, RenderQueue *queue,
        const RenderJob *job, int sceneResID, int sensorResID,
        int samplerResID) final;

    bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID) final;

    void cancel() final;

    MTS_DECLARE_CLASS()

private:
    LPMConfig config_;
    ref<ParallelProcess> process_;
    bool canceled_;

    lpm::PhotonContainer photons_;
    lpm::HashGrid photon_accel_;
};

MTS_NAMESPACE_END

#endif // __LPM_H