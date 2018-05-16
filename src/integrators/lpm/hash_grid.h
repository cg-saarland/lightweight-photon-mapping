#ifndef __HASH_GRID_H
#define __HASH_GRID_H

#include <numeric>
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>
#include "photon_map.h"

MTS_NAMESPACE_BEGIN

namespace lpm {

/// Returns the initializer for Bernstein's hash function
inline uint32_t bernstein_init() { return 5381; }

/// Hashes 4 bytes using Bernstein's hash
inline uint32_t bernstein_hash(uint32_t h, uint32_t d) {
    h = (h * 33) ^ ( d        & 0xFF);
    h = (h * 33) ^ ((d >>  8) & 0xFF);
    h = (h * 33) ^ ((d >> 16) & 0xFF);
    h = (h * 33) ^ ((d >> 24) & 0xFF);
    return h;
}

/// Returns the integer that is greater or equal to the logarithm base 2 of the argument.
template <typename T>
inline T closest_log2(T i) {
    T p = 1, q = 0;
    while (i > p) p <<= 1, q++;
    return q;
}

#define STACK_ARRAY(T, N) static_cast<T*>(alloca((N) * sizeof(T)));

struct CellIdx {
    int x;
    int y;

    CellIdx(int x, int y) : x(x), y(y) {}
    CellIdx() : x(0), y(0) {}
};

class HashGrid {
    typedef unsigned int uint;
public:
    void build(Photon* photons_begin, Photon* photons_end, float radius) {
        constexpr uint32_t inv_load_factor = 2;
        radius_        = radius;
        radius_sqr_    = radius * radius;
        cell_size_     = radius_ * 2.f;
        inv_cell_size_ = 1.f / cell_size_;

        auto photon_count = photons_end - photons_begin;
        if (cell_ends_.size() < size_t(photon_count * inv_load_factor))
            cell_ends_ = std::vector<int>(photon_count * inv_load_factor);

        // Compute the extents of the bounding box.
        bbox_ = AABB();
        for (Photon* it = photons_begin; it != photons_end; ++it)
            bbox_.expandBy(it->pos);

        auto extents = bbox_.max - bbox_.min;
        bbox_.max += extents * 0.001f;
        bbox_.min -= extents * 0.001f;

        // Distribute the photons to the HashGrid cells using Counting Sort.
        photons_.resize(photon_count);
        std::fill(cell_ends_.begin(), cell_ends_.end(), 0);

        // Count the number of photons in each cell.
        for (Photon* it = photons_begin; it != photons_end; ++it) {
            cell_ends_[cell_index(it->pos)]++;
        }

        // Set the cell_ends_[x] to the first index that belongs to the respective cell.
        int sum = 0;
        for(size_t i = 0; i < cell_ends_.size(); i++) {
            int temp = cell_ends_[i];
            cell_ends_[i] = sum;
            sum += temp;
        }

        // Assign the photons to the cells.
        for (Photon* it = photons_begin; it != photons_end; ++it) {
            const Point &pos = it->pos;
            const int target_idx = cell_ends_[cell_index(pos)]++;
            photons_[target_idx] = it;
        }
    }

    /// Writes pointers to the k nearest photons around a point into the given buffer.
    int query(const Point& query_pos, Photon** out, int k) const {
        // Check if the position is outside the bounding box.
        if (!bbox_.contains(query_pos)) return 0;

        auto cell = inv_cell_size_ * (query_pos - bbox_.min);
        const Point coord(
            std::floor(cell.x),
            std::floor(cell.y),
            std::floor(cell.z));

        const int px = int(coord.x);
        const int py = int(coord.y);
        const int pz = int(coord.z);

        auto fract_coord = -coord + cell;

        const int pxo = px + (fract_coord.x < 0.5f ? -1 : 1);
        const int pyo = py + (fract_coord.y < 0.5f ? -1 : 1);
        const int pzo = pz + (fract_coord.z < 0.5f ? -1 : 1);

        auto distances = STACK_ARRAY(float, k);
        int count = 0;
        for (int j = 0; j < 8; j++) {
            const int x = j & 4 ? pxo : px;
            const int y = j & 2 ? pyo : py;
            const int z = j & 1 ? pzo : pz;
            CellIdx active_range = cell_range(cell_index(x , y , z ));

            for (; active_range.x < active_range.y; active_range.x++) {
                auto& photon = photons_[active_range.x];
                const float dist_sqr = (query_pos - photon->pos).lengthSquared();

                if (dist_sqr <= radius_sqr_) {
                    if (count == k) {
                        if (distances[count - 1] < dist_sqr) continue;
                    } else count++;

                    // Insertion sort
                    distances[count - 1] = dist_sqr;
                    out[count - 1] = photon;
                    for (int l = count - 2; l >= 0; --l) {
                        if (distances[l] > dist_sqr) {
                            std::swap(distances[l], distances[l + 1]);
                            std::swap(out[l], out[l + 1]);
                        } else break;
                    }
                }
            }
        }

        return count;
    }

private:
    CellIdx cell_range(int cell_idx) const {
        if(cell_idx == 0) return CellIdx(0, cell_ends_[0]);
        return CellIdx(cell_ends_[cell_idx-1], cell_ends_[cell_idx]);
    }

    int cell_index(uint x, uint y, uint z) const {
        return int(((x * 73856093) ^ (y * 19349663) ^
            (z * 83492791)) % uint(cell_ends_.size()));
    }

    int cell_index(const Point &point) const {
        auto dist_min = inv_cell_size_ * (point - bbox_.min);
        int coord_x = std::floor(dist_min.x);
        int coord_y = std::floor(dist_min.y);
        int coord_z = std::floor(dist_min.z);
        return cell_index(coord_x, coord_y, coord_z);
    }

    AABB bbox_;
    std::vector<Photon*> photons_;
    std::vector<int> cell_ends_;

    float radius_;
    float radius_sqr_;
    float cell_size_;
    float inv_cell_size_;
};

} // namespace lpm

MTS_NAMESPACE_END

#endif // __HASH_GRID_H