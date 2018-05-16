#include "mis.h"

MTS_NAMESPACE_BEGIN

namespace lpm {

void VmMis::update_hit(Vertex* v, PathSampler::PathOrigin type) {
    v->partial_mis_sum /= v->cos_to_prev;
}

void VmMis::update_bounce(Vertex* next, PathSampler::PathOrigin type) {
    auto current = next->predec;

    // propagate the partial mis weights to the new vertex
    next->partial_mis_sum = current->partial_mis_sum;

    if (next->last_specular) {
        next->partial_mis_sum *= next->cos_from_prev;
        return;
    }

    next->partial_mis_sum *= current->pdf_reverse;

    if (type == PathSampler::FROM_PIXEL && current->sub_path_len == 1) { // first bounce after the camera
        // add weight for connecting to the camera VS all techniques requiring a camera sub-path
        auto pdf_sensor_surf_area = current->pdf_forward * current->cos_to_prev / current->dist_sqr;
        next->partial_mis_sum += light_paths_per_camera_path / pdf_sensor_surf_area;
    }

    if (type == PathSampler::FROM_EMITTER && current->sub_path_len == 1) { // first bounce after the light source
        // add weight for next event estimation from the path tracer VS all emission-based techniques
        Float pdf_next_event_surf_area = current->predec->pdf_reverse;
        Float pdf_emission_surf_area   = current->pdf_forward * current->cos_to_prev / current->dist_sqr;
        next->partial_mis_sum += pdf_next_event_surf_area / pdf_emission_surf_area;
    }

    // add weight for merging here instead
    if (enable_merging && type == PathSampler::FROM_EMITTER && (current->sub_path_len > 1 || merge_di))
        next->partial_mis_sum += merge_accept_weight;
    else if (enable_merging && type == PathSampler::FROM_PIXEL && (current->sub_path_len > 1 || merge_primary))
        next->partial_mis_sum += merge_accept_weight;

    next->partial_mis_sum *= next->cos_from_prev / next->pdf_forward;
}

void VmMis::init_sensor(Vertex* v, PathSampler::PathOrigin type) {
    v->partial_mis_sum = 0.0f;
}

void VmMis::init_emitter(Vertex* next, PathSampler::PathOrigin type) {
    // TODO special case for infinite light sources (aka those without meaningful surface area pdf)
    // TODO special case (do not add) for emitters that cannot be intersected

    // include weight for the path tracer hitting the light source VS all emission-based techniques
    next->partial_mis_sum = next->cos_from_prev / next->pdf_forward;
}

Float VmMis::weight_pt_nee(const Vertex* v, Float cos_surf, Float cos_emitter, Float pdf_next_event_solid_angle,
                           Float pdf_hit_solid_angle, Float pdf_emit_solid_angle, Float pdf_reverse) const {
    // compute weight of hitting the light VS next event estimation
    // TODO special case for delta light source!
    auto tech_on_light_sub_path = pdf_hit_solid_angle / pdf_next_event_solid_angle; // the geometry term cancels out

    auto partial_cam = pdf_reverse * v->partial_mis_sum;

    // add the weight for tracing the full path via light tracing, if "bounce()" was not called
    if (light_trace_di && v->sub_path_len == 1) {
        auto pdf_sensor_surf_area = v->pdf_forward * v->cos_to_prev / v->dist_sqr;
        partial_cam += light_paths_per_camera_path / pdf_sensor_surf_area;
    }

    // add the weight for merging
    if (enable_merging && merge_di && (v->sub_path_len > 1 || merge_primary)) {
        partial_cam += merge_accept_weight;
    }

    auto tech_on_cam_sub_path = (pdf_emit_solid_angle * cos_surf) / (pdf_next_event_solid_angle * cos_emitter)
                              * partial_cam;

    return 1.0f / (tech_on_light_sub_path + 1.0f + tech_on_cam_sub_path);
}

Float VmMis::weight_pt_hit(const Vertex* v, Float pdf_emit_solid_angle, Float pdf_next_event_surf_area) const {
    if (v->sub_path_len == 1) {
        // this assumes only the path tracer renders directly visible lights (reasonable method)
        return 1.0f;
    }

    // if DI light tracing is disabled, do not consider that technique
    auto partial = v->partial_mis_sum;
    if (v->sub_path_len == 2 && !v->last_specular && !light_trace_di) {
        if (enable_merging && merge_primary)
            partial = merge_accept_weight * v->cos_from_prev / (v->pdf_forward * v->cos_to_prev);
        else partial = 0.0f;
    }

    // if DI merging is disabled, do not consider merging at the previous vertex along the path
    if (!v->last_specular && enable_merging && !merge_di && (v->sub_path_len > 2 || merge_primary))
        partial -= merge_accept_weight * v->cos_from_prev / (v->pdf_forward * v->cos_to_prev);
    partial = std::max(partial, 0.0f);

    // add the weight for next event estimation, if the last bounce was not specular
    Float nee_weight = 0.0f;
    if (!v->last_specular)
        nee_weight = pdf_next_event_surf_area / (v->pdf_forward * v->cos_to_prev / v->dist_sqr);

    Float other_techs = pdf_emit_solid_angle * partial + nee_weight;

    return 1.0f / (other_techs + 1.0f);
}

Float VmMis::weight_lt_nee(const Vertex* v, Float pdf_sensor_surf_area, Float pdf_reverse) const {
    if (!light_trace_di && v->sub_path_len == 1) return 0.0f;

    // include weight for path tracer next event estimation if not yet part of partial
    Float pt_nee = 0.0f;
    if (v->sub_path_len == 1) {
        Float pdf_next_event_surf_area = v->predec->pdf_reverse;
        Float pdf_emission_surf_area   = v->pdf_forward * v->cos_to_prev / v->dist_sqr;
        pt_nee += pdf_next_event_surf_area / pdf_emission_surf_area;
    }

    // consider merging if this is not a direct illumination path
    bool consider_merge = (v->sub_path_len > 1 || merge_di) && enable_merging && merge_primary;
    Float merge_w = consider_merge ? merge_accept_weight : 0.0f;

    auto partial = pdf_reverse * v->partial_mis_sum + pt_nee + merge_w;
    auto tech_on_light_sub_path = pdf_sensor_surf_area / light_paths_per_camera_path * partial;

    return 1.0f / (tech_on_light_sub_path + 1.0f);
}

Float VmMis::weight_lt_hit(const Vertex* v) const {
    // this technique is not sensible and always disabled
    return 0.0f;
}

Float VmMis::weight_merge(const Vertex* cam_v, const Photon* light_v, Float pdf_reverse_cam, Float pdf_reverse_light) const {
    // do not use merging for direct illumination, path tracing is always more efficient at that (unless specified otherwise)
    if ((light_v->depth == 1 && !merge_di) || !enable_merging || (cam_v->sub_path_len == 1 && !merge_primary))
        return 0.0f;

    // include the weight for tracing the full path via light tracing, if "bounce()" was not called
    Float lt_nee = 0.0f;
    if (light_trace_di && cam_v->sub_path_len == 1) {
        auto pdf_sensor_surf_area = cam_v->pdf_forward * cam_v->cos_to_prev / cam_v->dist_sqr;
        lt_nee = light_paths_per_camera_path / pdf_sensor_surf_area;
    }

    auto partial_cam   = pdf_reverse_cam   *   cam_v->partial_mis_sum + lt_nee;
    auto partial_light = pdf_reverse_light * light_v->partial_mis_sum;

    auto tech_on_cam_sub_path   = inv_merge_accept_weight * partial_cam;
    auto tech_on_light_sub_path = inv_merge_accept_weight * partial_light;
    return 1.0f / (tech_on_cam_sub_path + tech_on_light_sub_path + 1.0f);
}

} // namespace lpm

MTS_NAMESPACE_END