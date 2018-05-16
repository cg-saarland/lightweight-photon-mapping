#ifndef __PATH_SAMPLER_H
#define __PATH_SAMPLER_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/mempool.h>

MTS_NAMESPACE_BEGIN

namespace lpm {

// TODO many reduntant values here - recomputation could be faster than storing everything.

struct Vertex {
    Vertex* predec;
    Vertex* succ;

    /// pdf of generating this vertex from the previous one - solid angle measure.
    Float pdf_forward;

    /// pdf of generating the previous vertex when starting from this one - solid angle measure.
    /// For the first vertex of a light sub-path: the pdf of next event estimation sampling this vertex (surface area measure)
    Float pdf_reverse;

    /// Sum of the pdf ratios of the techniques along the current sub-path.
    /// Most combined estimators can use this to efficiently compute MIS without iterating over the path or storing it.
    Float partial_mis_sum;

    /// False if the ray left the scene
    bool on_surface;

    /// True if the ray that generated this vertex was sampled from a discrete distribution
    bool last_specular;

    /// True, if the BSDF at this vertex is purely specular (i.e. there is no point in performing density estimation and connections)
    bool specular;

    /// Positon of the vertex in the world
    Point pos;

    /// Cosine between the normal and the vector towards the predecessor
    Float cos_to_prev;

    /// Cosine at the predecessor between the normal and direction towards this vertex
    Float cos_from_prev;

    /// Squared distance between this vertex and its predecessor
    Float dist_sqr;

    /// Weigth carried by the sub-path up until this vertex.
    Spectrum weight;

    /// Geometric (i.e. "correct") normal of the surface at the vertex (if it is not inside a medium)
    Vector geom_normal;

    /// Length of the sub-path up to and including this vertex
    int sub_path_len;

    /// The estimated (and MIS weighted) image contribution of this vertex, excluding the emission.
    /// Used for guiding the emission in primary sample space (the Le term is already importance sampled)
    Spectrum guiding_contrib;

    /// The throughput of the path (weight without emission) up until this vertex.
    Spectrum throughput;

    /// Is this vertex useful when performing density estimation on the other end? (currently only used for photons)
    bool useful_for_merging;
};

template <typename T>
class SimpleBuffer {
public:
    SimpleBuffer(size_t count = 128) : next_(0), chunk_sz_(count) {
        chunks_.emplace_back(chunk_sz_);
    }

    T* alloc() {
        // ensure that there is enough memory for the full path!
        if (chunks_.size() * chunk_sz_ <= next_)
            chunks_.emplace_back(chunk_sz_);
        T* e = &(chunks_[next_ / 128][next_ % 128]);
        next_++;
        return e;
    }

    void reset() { next_ = 0; }

private:
    std::vector<std::vector<T>> chunks_;
    size_t chunk_sz_;
    size_t next_;
};

using VertexMemPool = SimpleBuffer<Vertex>;

// TODO / REFACTOR: use std::list with a custom allocator instead of re-inventing the wheel?
struct Path {
    Vertex* first;
    Vertex* last;

    uint32_t len;

    Path() : first(nullptr), last(nullptr), len(0) {}

    void append(VertexMemPool* pool);
    void remove_last();
};

struct MisComputer;

/// Incrementally constructs a path from either the emitters or the sensor.
/// Virtual callbacks are triggered for every event along the path.
/// The importance sampling at every vertex can also be controlled by derived classes.
class PathSampler {
public:
    enum PathOrigin {
        FROM_PIXEL,
        FROM_EMITTER
    };

    Point2 pixel;

    /// Initializes the sampler.
    /// \param pool     the memory pool used to store the vertices of the path, or nullptr to not store any vertices.
    PathSampler(int rr_depth, int max_depth, Scene*, Sampler*, Sensor*, VertexMemPool* pool = nullptr);

    virtual ~PathSampler() {}

    void trace(PathOrigin, const Point2i& px, MisComputer*);

    /// Returns the path.
    /// Only valid as long as the memory pool is valid.
    Path path() { return path_; }

protected:
    PathOrigin origin_;

    /// If this is not null, it is used to store the vertices of the path
    VertexMemPool* pool_;
    Path path_;

    ref<Sampler> sampler_;
    ref<Scene>   scene_;
    ref<Sensor>  sensor_;

    int rr_depth_;
    int max_depth_;

    /// Samples an outgoing ray from an emitter in the scene.
    /// The default implementation uses the \see{sampleEmitter...} methods of \see{Scene}.
    /// Calls the \see{handle_emitter} callback function after sampling a position on the light source
    ///
    /// \returns the MC estimate of the emission along the ray
    virtual Spectrum sample_emitter(Float time, PositionSamplingRecord&, DirectionSamplingRecord&, Float& uvpdf);

    virtual Spectrum sample_sensor(Float time, PositionSamplingRecord&, DirectionSamplingRecord&,
                                   const Point2& pixel, Sensor* sensor);

    /// Samples a scattering position inside a medium
    virtual bool sample_medium(const Medium*, Ray&, Intersection&, MediumSamplingRecord&);

    /// Samples a continuation ray from the phase function at a point in a medium
    virtual Float sample_phase_fn(const Medium*, PhaseFunctionSamplingRecord&);

    /// Samples a continuation ray from a point on a surface
    virtual Spectrum sample_surface(const BSDF*, BSDFSamplingRecord&, Float& pdf_forward, Float& pdf_reverse);

    virtual Float pdf_next_event(const Emitter* em, const PositionSamplingRecord& pos, const Point& from, const Vector& from_normal);

    virtual Float pdf_emit(const Emitter* em, const PositionSamplingRecord& pos, const Vector& dir);

    virtual void handle_surface(int depth, int nullInteractions, bool delta,
                                const Intersection &its, const Medium *medium,
                                const Spectrum &emitted, const Spectrum& throughput) {}

    virtual void handle_medium(int depth, int nullInteractions, bool delta,
                               const MediumSamplingRecord &mRec, const Medium *medium,
                               const Vector &wi, const Spectrum &weight) {}

    /// Called by \see{sample_emitter} after a point on an emitter was sampled.
    /// Useful to connect the emitter to the sensor if path tracing is disabled and the sensor has no area.
    /// Default implementation adds the vertex to the path.
    virtual void handle_emitter(const Medium *medium, const PositionSamplingRecord&, const Spectrum& power) {}

    /// Called when nothing was intersected.
    /// Useful to add the contribution of the environment map, etc.
    virtual void handle_miss(int depth, int nullInteractions, bool delta, const Spectrum& weight) {}
};

} // namespace lpm

MTS_NAMESPACE_END

#endif // __PATH_SAMPLER_H