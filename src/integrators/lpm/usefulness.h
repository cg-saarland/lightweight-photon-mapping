#ifndef __USEFULNESS_H
#define __USEFULNESS_H

#include <mitsuba/mitsuba.h>
#include "path_sampler.h"

MTS_NAMESPACE_BEGIN
namespace lpm {

bool is_useful_photon(const Path* p, Float min_num_lp, Float r) {
    if (p->len < 1) return false; // the path never left the light source
    if (p->last->specular) return false; // the photon is on a purely specular surface

    // compute the ratio of pdf(path tracing) / pdf(light tracing)
    // this prevents division by zero (the light tracer actually sampled this path, so no pdf can be zero)
    Float pdf_ratio = 0.0f;

    auto vert_on_emitter = p->first;
    auto first_surf      = p->first->succ;

    // path tracer next event VS emission
    if (!first_surf->specular) {
        Float pdf_next_event_surf_area = vert_on_emitter->pdf_reverse;
        Float pdf_emission_surf_area   = first_surf->pdf_forward * first_surf->cos_to_prev / first_surf->dist_sqr;
        pdf_ratio += pdf_next_event_surf_area / pdf_emission_surf_area;
    }

    if (p->len > 1) {
        // path tracer hitting the emitter VS emission
        // the squared distance from the geometry terms cancels out
        pdf_ratio += first_surf->cos_from_prev * first_surf->pdf_reverse / (first_surf->pdf_forward * first_surf->cos_to_prev);

        // going forward VS going backwards along the inner path
        // the squared distances cancel out
        for (auto v = first_surf->succ; v && v->succ && v != p->last; v = v->succ) {
            // multiply by the pdfs for each direction, if the sampled BSDF component was not specular
            if (!v->succ->last_specular) pdf_ratio *= v->pdf_reverse;
            if (!v->last_specular) pdf_ratio /= v->pdf_forward;

            // ratio of geometry terms
            pdf_ratio *= v->cos_from_prev / v->cos_to_prev;
        }
    }

    // The last vertex (the photon in question) has to be treated separately:
    // the reverse pdf is not known and depends on the sensor sub-path anyway.
    // We assume a diffuse surface to compute an upper bound of the usefulness.
    // The cosine term cancels out.
    Float pdf_revserse_diffuse = /*p->last->cos_to_prev * */M_1_PI;
    pdf_ratio *= pdf_revserse_diffuse * p->last->cos_from_prev / (p->last->pdf_forward/* * p->last->cos_to_prev*/);

    // The denominator pdf (light tracer) now has measure m^(2 * num_verts)
    // whereas the numerator (path tracer) does not consider sampling the last vertex
    // and therefore has measure m^(2 * (num_verts - 1))

    // In order to achieve a meaningful comparison, we assume a uniform distribution of vertices
    // from the path tracer within the photon mapping radius around the photon in question.
    Float uniform_disc_pdf = 1.0f / (M_PI * r * r);

    if (pdf_ratio == 0.0f) {
        // the path tracer cannot sample this path at all, probably due to a delta distribution in the light source.
        return true;
    }

    auto usefulness = 1.0f / (uniform_disc_pdf * pdf_ratio);

    if (min_num_lp * usefulness > 1.0f) {
        return true; // even with only "min_num_lp" paths, the photon mapper has a higher probability
    } else
        return false; // the path tracer can capture this path efficiently
}

} // namespace lpm
MTS_NAMESPACE_END

#endif // __USEFULNESS_H