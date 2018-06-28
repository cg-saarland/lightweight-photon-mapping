#ifndef __PIXEL_CLASSIFICATION_H
#define __PIXEL_CLASSIFICATION_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/bitmap.h>

MTS_NAMESPACE_BEGIN
namespace lpm {

int num_pixels_requiring_light_paths(const Bitmap* img_lpcontrib, const Bitmap* img_full) {
    int count = 0;
    for (int i = 0; i < img_full->getSize().x; ++i)
    for (int j = 0; j < img_full->getSize().y; ++j) {
        Float luminance_lpcontrib = img_lpcontrib->getPixel(Point2i(i, j)).getLuminance();
        Float luminance_full = img_full->getPixel(Point2i(i, j)).getLuminance();
        if (luminance_lpcontrib / luminance_full > 0.01f) // TODO turn hard-coded constant into parameter
            count++;
    }
    return count;
}

} // namespace lpm
MTS_NAMESPACE_END

#endif // __PIXEL_CLASSIFICATION_H