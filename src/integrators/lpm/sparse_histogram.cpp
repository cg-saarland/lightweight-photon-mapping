#include <numeric>
#include "sparse_histogram.h"

MTS_NAMESPACE_BEGIN

namespace lpm {

SparseHistogram4D::SparseHistogram4D(uint32_t res_u_prim, uint32_t res_v_prim, uint32_t res_u_sec, uint32_t res_v_sec,
                                     uint32_t filter_radius)
: res_prim_(res_u_prim, res_v_prim)
, res_sec_(res_u_sec, res_v_sec)
, inv_res_prim_(1.0f / res_u_prim, 1.0f / res_v_prim)
, inv_res_sec_ (1.0f / res_u_sec , 1.0f / res_v_sec)
, inv_normalization_(1.0f)
, filter_radius_(filter_radius)
, pdf_(res_u_prim * res_v_prim * res_u_sec * res_v_sec)
, cdf_(res_u_prim * res_v_prim * res_u_sec * res_v_sec)
, dir_marginal_(res_u_sec * res_v_sec)
, pos_marginal_(res_u_prim * res_v_prim)
, num_bins_(res_u_prim * res_v_prim * res_u_sec * res_v_sec)
{
    num_samples.store(0);
}

void SparseHistogram4D::normalize() {
    Float sum = std::accumulate(pdf_.begin(), pdf_.end(), 0.0f);
    if (sum == 0.0f) return;
    inv_normalization_ = 1.0f / sum;

    SAssert(cdf_.size() == pdf_.size());

    // compute the CDF
    Float accum = 0.0f;
    for (uint32_t i = 0; i < num_bins_; ++i) {
        accum += pdf_[i] * inv_normalization_;
        cdf_[i] = accum;
    }
}

void SparseHistogram4D::subdiv() {
    if (num_samples < num_bins_ * 8) // TODO make configurable
        return;

    // only split if the maximum has not been reached, and only if the dimension is used
    auto split_pos = res_prim_.x < 256 && res_prim_.y < 256 && res_prim_.x != 1;
    auto split_dir = res_sec_.x  < 256 && res_sec_.y  < 256 && res_sec_.x  != 1;

    SLog(EInfo, "splitting dir %i (res=%i), pos %i (res=%i)",
         split_dir, res_sec_.x * 4, split_pos, res_prim_.x * 4);

    if (!split_dir && !split_pos) return;

    // use the cdf as a buffer
    std::copy(pdf_.begin(), pdf_.end(), cdf_.begin());
    auto& buffer = cdf_;

    // increase resolution and allocate memory
    num_bins_ *= (split_pos ? 4 : 1) * (split_dir ? 4 : 1);

    pdf_.resize(num_bins_);

    // the weight is not adjusted which effectively reduces the weight of the lower resultion by a factor of 4
    Float factor = 1.0f;// / ((split_dir ? 4.0f : 1.0f) * (split_pos ? 4.0f : 1.0f));

    // adapt the histogram by duplicating the cells
    uint32_t idx = 0;
    for (uint32_t l = 0; l < res_prim_.y; ++l) {
       for (uint32_t k = 0; k < res_prim_.x; ++k) {
            for (uint32_t j = 0; j < res_sec_.y; ++j) {
                for (uint32_t i = 0; i < res_sec_.x; ++i) {
                    auto binidx = res_sec_.x * (res_sec_.y * (res_prim_.x * l + k) + j) + i;
                    pdf_[idx++] = buffer[binidx] * factor;
                    if (split_dir) pdf_[idx++] = buffer[binidx] * factor;
                }
                if (split_dir) {
                    std::copy(pdf_.begin() + idx - 2 * res_sec_.x, pdf_.begin() + idx, pdf_.begin() + idx);
                    idx += 2 * res_sec_.x;
                }
            }
            if (split_pos) {
                int stride = split_dir ? 4 : 1;
                std::copy(pdf_.begin() + idx - stride * res_sec_.x * res_sec_.y, pdf_.begin() + idx, pdf_.begin() + idx);
                idx += stride * res_sec_.x * res_sec_.y;
            }
        }
        if (split_pos) {
            auto size = res_sec_.x * res_sec_.y * res_prim_.x * 2 * (split_dir ? 4 : 1);
            std::copy(pdf_.begin() + idx - size, pdf_.begin() + idx, pdf_.begin() + idx);
            idx += size;
        }
    }

    SAssert(idx == num_bins_);

    // adapt the marginals
    if (split_pos) {
        int64_t idx = res_prim_.x * res_prim_.y * 4 - 1;
        pos_marginal_.resize(idx + 1);
        for (int64_t j = res_prim_.y - 1; j >= 0; --j) {
            for (int64_t i = res_prim_.x - 1; i >= 0; --i) {
                pos_marginal_[idx--] = pos_marginal_[i + res_prim_.x * j] * factor;
                pos_marginal_[idx--] = pos_marginal_[i + res_prim_.x * j] * factor;
            }
            std::copy(pos_marginal_.begin() + idx + 1, pos_marginal_.begin() + idx + 1 + 2 * res_prim_.x,
                      pos_marginal_.begin() + idx + 1 - 2 * res_prim_.x);
            idx -= 2 * res_prim_.x;
        }
        SAssert(idx == -1);
        res_prim_ *= 2;
        inv_res_prim_ *= 0.5f;
    }

    if (split_dir) {
        int64_t idx = res_sec_.x * res_sec_.y * 4 - 1;
        dir_marginal_.resize(idx + 1);
        for (int64_t j = res_sec_.y - 1; j >= 0; --j) {
            for(int64_t i = res_sec_.x - 1; i >= 0; --i) {
                dir_marginal_[idx--] = dir_marginal_[i + res_sec_.x * j] * factor;
                dir_marginal_[idx--] = dir_marginal_[i + res_sec_.x * j] * factor;
            }
            std::copy(dir_marginal_.begin() + idx + 1, dir_marginal_.begin() + idx + 2 * res_sec_.x,
                      dir_marginal_.begin() + idx + 1 - 2 * res_sec_.x);
            idx -= 2 * res_sec_.x;
        }
        SAssert(idx == -1);
        res_sec_ *= 2;
        inv_res_sec_ *= 0.5f;
    }

    SAssert(num_bins_ == res_prim_.x * res_prim_.y * res_sec_.x * res_sec_.y);

    cdf_.resize(num_bins_);
}

Float SparseHistogram4D::transform(Point2& prim, Point2& sec, Float bin_sample) const {
    auto i = std::lower_bound(cdf_.begin(), cdf_.end(), bin_sample);

    // handle corner cases
    auto bin_idx = i - cdf_.begin();
    bin_idx = std::min(num_bins_ - 1, (uint32_t)std::max((ptrdiff_t)0, bin_idx));

    auto pdf = pdf_[bin_idx] * inv_normalization_ * num_bins_;

    // if the pdf is still zero, quit
    if (pdf == 0.0f) return 0.0f;

    // sample a point within the bin
    Point2 prim_min, sec_min;
    coords(bin_idx, prim_min, sec_min);

    prim.x = std::min(prim.x, 1.0f - Epsilon);
    prim.y = std::min(prim.y, 1.0f - Epsilon);
    sec.x  = std::min(sec.x , 1.0f - Epsilon);
    sec.y  = std::min(sec.y , 1.0f - Epsilon);

    prim.x = prim_min.x + prim.x * inv_res_prim_.x;
    prim.y = prim_min.y + prim.y * inv_res_prim_.y;
    sec.x  = sec_min.x  + sec.x  * inv_res_sec_.x;
    sec.y  = sec_min.y  + sec.y  * inv_res_sec_.y;

    return pdf;
}

void SparseHistogram4D::put(const Point2& prim, const Point2& sec, Float weight) {
    auto px = std::min(static_cast<const int32_t>(prim.x * res_prim_.x), int32_t(res_prim_.x) - 1);
    auto py = std::min(static_cast<const int32_t>(prim.y * res_prim_.y), int32_t(res_prim_.y) - 1);
    auto sx = std::min(static_cast<const int32_t>(sec.x  * res_sec_.x ), int32_t(res_sec_.x)  - 1);
    auto sy = std::min(static_cast<const int32_t>(sec.y  * res_sec_.y ), int32_t(res_sec_.y)  - 1);

    if (filter_radius_ == 0) {
        (*this)(px, py, sx, sy) += weight;
        pos_marginal_[px + res_prim_.x * py] += weight;
        dir_marginal_[sx + res_sec_.x  * sy] += weight;
        return;
    }

    const int32_t fr = filter_radius_;

    auto w = weight / (fr * fr * fr * fr);

    for (auto up = std::max(0, px - fr);
              up < std::min(int32_t(res_prim_.x), px + fr);
              ++up)
    for (auto vp = std::max(0, py - fr);
              vp < std::min(int32_t(res_prim_.y), py + fr);
              ++vp)
    for (auto us = std::max(0, sx - fr);
              us < std::min(int32_t(res_sec_.x), sx + fr);
              ++us)
    for (auto vs = std::max(0, sy - fr);
              vs < std::min(int32_t(res_sec_.y), sy + fr);
              ++vs)
    {
        (*this)(up, vp, us, vs) += w;
        pos_marginal_[up + res_prim_.x * vp] += w;
        dir_marginal_[us + res_sec_.x  * vs] += w;
    }
}

Float& SparseHistogram4D::operator()(const Point2& prim, const Point2& sec) {
    return pdf_[bin(prim, sec)];
}

Float SparseHistogram4D::operator()(const Point2& prim, const Point2& sec) const {
    return pdf_[bin(prim, sec)];
}

Float& SparseHistogram4D::operator()(uint32_t u_prim, uint32_t v_prim, uint32_t u_sec, uint32_t v_sec) {
    return pdf_[res_sec_.x * (res_sec_.y * (res_prim_.x * v_prim + u_prim) + v_sec) + u_sec];
}

Float SparseHistogram4D::operator()(uint32_t u_prim, uint32_t v_prim, uint32_t u_sec, uint32_t v_sec) const {
    return pdf_[res_sec_.x * (res_sec_.y * (res_prim_.x * v_prim + u_prim) + v_sec) + u_sec];
}

Float SparseHistogram4D::pdf(const Point2& prim, const Point2& sec) const {
    return (*this)(prim, sec) * inv_normalization_ * num_bins_;
}

uint32_t SparseHistogram4D::bin(const Point2& prim, const Point2& sec) const {
    auto prim_x_bin = std::min(static_cast<const uint32_t>(prim.x * res_prim_.x), res_prim_.x - 1);
    auto prim_y_bin = std::min(static_cast<const uint32_t>(prim.y * res_prim_.y), res_prim_.y - 1);
    auto  sec_x_bin = std::min(static_cast<const uint32_t>( sec.x *  res_sec_.x),  res_sec_.x - 1);
    auto  sec_y_bin = std::min(static_cast<const uint32_t>( sec.y *  res_sec_.y),  res_sec_.y - 1);

    return res_sec_.x * (res_sec_.y * (res_prim_.x * prim_y_bin + prim_x_bin) + sec_y_bin) + sec_x_bin;
}

void SparseHistogram4D::coords(uint32_t bin_idx, Point2& prim_min, Point2& sec_min) const {
    // extract the 4D integer coordinates
    auto sec_x_bin = bin_idx % res_sec_.x;
    bin_idx /= res_sec_.x;
    auto sec_y_bin = bin_idx % res_sec_.y;
    bin_idx /= res_sec_.y;
    auto prim_x_bin = bin_idx % res_prim_.x;
    auto prim_y_bin = bin_idx / res_prim_.y;

    // compute position on the unit hypercube
    prim_min.x = prim_x_bin / (Float)res_prim_.x;
    prim_min.y = prim_y_bin / (Float)res_prim_.y;
    sec_min.x  =  sec_x_bin / (Float) res_sec_.x;
    sec_min.y  =  sec_y_bin / (Float) res_sec_.y;
}

void SparseHistogram4D::dump(std::ostream& str) const {
    str.write(reinterpret_cast<const char*>(&res_prim_.x), sizeof(res_prim_.x));
    str.write(reinterpret_cast<const char*>(&res_prim_.y), sizeof(res_prim_.y));
    str.write(reinterpret_cast<const char*>(&res_sec_.x ), sizeof(res_sec_.x ));
    str.write(reinterpret_cast<const char*>(&res_sec_.y ), sizeof(res_sec_.y ));

    int index = 0;
    for (uint32_t pv = 0; pv < res_prim_.y; ++pv)
    for (uint32_t pu = 0; pu < res_prim_.x; ++pu)
    for (uint32_t sv = 0; sv < res_sec_.y; ++sv)
    for (uint32_t su = 0; su < res_sec_.x; ++su)
    {
        auto val = (*this)(pu, pv, su, sv);
        str.write(reinterpret_cast<const char*>(&val), sizeof(val));
    }
}


Histogram1D::Histogram1D(uint32_t count)
: weights_(count), cdf_(count)
{
}

void Histogram1D::put(uint32_t idx, Float w) {
    weights_[idx] += w;
}

void Histogram1D::normalize() {
    Float sum = std::accumulate(weights_.begin(), weights_.end(), 0.0f);
    if (sum == 0.0f) return;
    inv_normalization_ = 1.0f / sum;

    // compute the CDF
    Float accum = 0.0f;
    for (uint32_t i = 0; i < weights_.size(); ++i) {
        accum += weights_[i] * inv_normalization_;
        cdf_[i] = accum;
    }
}

uint32_t Histogram1D::sample(Float u, Float& pdf) const {
    if (cdf_.size() == 1) return 0;

    auto i = std::lower_bound(cdf_.begin(), cdf_.end(), u);

    // handle corner cases
    auto bin_idx = i - cdf_.begin();
    bin_idx = std::min(cdf_.size() - 1, (size_t)std::max((ptrdiff_t)0, bin_idx));

    pdf = weights_[bin_idx] * inv_normalization_;

    return bin_idx;
}

Float Histogram1D::pdf(uint32_t idx) const {
    return weights_[idx] * inv_normalization_;
}

void Histogram1D::dump(std::ostream& str) const {
    str.write(reinterpret_cast<const char*>(weights_.data()), sizeof(Float) * weights_.size());
}

} // namespace lpm

MTS_NAMESPACE_END
