#ifndef _EDGEDETECTOR_H
#define _EDGEDETECTOR_H

#include <image.h>

// Marr-Hildreth edge detector
VISImage<float> marrHildreth(VISImage<float> I, float sigma, float threshold);
// Finds the zero crossings in the given image
// Marks it black on a white background
VISImage<float> zeroCrossings(VISImage<float> I);
// Finds the zero crossings for the given image with a threshold on the gradient values
VISImage<float> zeroCrossings(VISImage<float> I, float threshold);

#endif /* EDGEDETECTOR_H */

