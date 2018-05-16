#ifndef __MIS_H
#define __MIS_H

#include "path_sampler.h"
#include "photon_map.h"

MTS_NAMESPACE_BEGIN

namespace lpm {

/// Updates the partial sum of MIS techniques along a sub-path.
struct MisComputer {
    virtual void update_hit   (Vertex*, PathSampler::PathOrigin) = 0;
    virtual void update_bounce(Vertex*, PathSampler::PathOrigin) = 0;
    virtual void init_sensor  (Vertex*, PathSampler::PathOrigin) = 0;
    virtual void init_emitter (Vertex*, PathSampler::PathOrigin) = 0;

    virtual Float weight_pt_hit(const Vertex*, Float pdf_emit_solid_angle, Float pdf_next_event_surf_area) const = 0;
    virtual Float weight_pt_nee(const Vertex*, Float cos_surf, Float cos_emitter, Float pdf_next_event_solid_angle,
                                Float pdf_hit_solid_angle, Float pdf_emit_solid_angle, Float pdf_reverse) const = 0;

    virtual Float weight_lt_hit(const Vertex*) const = 0;
    virtual Float weight_lt_nee(const Vertex*, Float pdf_camera_surf_area, Float pdf_reverse) const = 0;

    virtual Float weight_merge(const Vertex* cam_v, const Photon* light_v, Float pdf_reverse_cam, Float pdf_reverse_light) const = 0;
};

/// MIS weights for two-way path tracing: light tracing + path tracing (aka bdpt without connections)
struct VmMis : public MisComputer {
    /// The number of samples from the light tracer for the whole (current) image, divided by the
    /// number of samples from the path tracer for the whole (current) image.
    Float light_paths_per_camera_path;

    Float merge_accept_weight;
    Float inv_merge_accept_weight;

    /// Set to false to disregard camera paths of length 2
    /// This can significantly reduce the variance due to the sub-optimal MIS weights!
    bool light_trace_di;

    bool enable_merging;
    bool merge_primary;
    bool merge_di;

    VmMis(Float lp_per_cp, Float radius, int num_light_paths, bool lt_di, bool merge, bool merge_prim, bool merge_di)
    : light_paths_per_camera_path(lp_per_cp), light_trace_di(lt_di)
    , merge_accept_weight(M_PI * radius * radius * num_light_paths)
    , inv_merge_accept_weight(1.0f / merge_accept_weight)
    , enable_merging(merge)
    , merge_primary(merge_prim)
    , merge_di(merge_di)
    {}

    void update_hit   (Vertex*, PathSampler::PathOrigin) final;
    void update_bounce(Vertex*, PathSampler::PathOrigin) final;
    void init_sensor  (Vertex*, PathSampler::PathOrigin) final;
    void init_emitter (Vertex*, PathSampler::PathOrigin) final;

    Float weight_pt_hit(const Vertex*, Float pdf_emit_solid_angle, Float pdf_next_event_surf_area) const final;
    Float weight_pt_nee(const Vertex*, Float cos_surf, Float cos_emitter, Float pdf_next_event_solid_angle,
                        Float pdf_hit_solid_angle, Float pdf_emit_solid_angle, Float pdf_reverse) const final;

    Float weight_lt_hit(const Vertex*) const final;
    Float weight_lt_nee(const Vertex*, Float pdf_sensor_surf_area, Float pdf_reverse) const final;

    Float weight_merge(const Vertex* cam_v, const Photon* light_v, Float pdf_reverse_cam, Float pdf_reverse_light) const final;
};

} // namespace lpm

MTS_NAMESPACE_END

#endif // __MIS_H