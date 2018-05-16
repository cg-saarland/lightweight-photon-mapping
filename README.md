# Lightweight Photon Mapping (Mitsuba)

Simplified Mitsuba implementation of the "Lightweight Photon Mapping" algorithm. This is not the implementation that was used in the paper.
The implementation here lacks many optimizations and only works for area light sources (triangle meshes or rectangles).
Not every feature is implemented, and by far not every corner case has been tested.
In particular, there is no support for volume rendering (although adding it should be straighforward),
no support for subsurface scattering, no support for low-discrepancy samplers, and some missing special cases in the
MIS weight computations, like delta light sources and infinitely far away light sources.
Also, the performance is not optimized at all. Some parts (in particular the MIS computations and the usefulness estimation) are deliberately implemented in a less than optimal way. The focus was on the flexibility and readability of the code.

The implementation is based on Mitsuba 0.6. See the [Mitsuba](https://www.mitsuba-renderer.org/) webpage for build instructions, dependencies, and other documentation.

## Overview of the Source Code

Our extensions to the original Mitsuba source code consist of the following major components distributed across multiple files:
- An integrator implementing the vertex merging (VM) algorithm (vertex connection and merging without bidirectional connections)
- A function to compute the usefulness of a photon
- A simple emission guiding approach using a 4D histogram for each light source in the scene
- Functions in the emitter classes to compute the inverse mapping of the position and direction sampling

The core contribution of the lightweight photon mapping paper, determining the usefulness of a photon, is implemented in the file *usefulness.h*.
Histograms of the image contribution of (useful) photons are built and used for the importance sampling of emission in the files *emission_sampler.h* and *emission_sampler.cpp*.

The integrator itself is defined in the files *lpm.h* and *lpm.cpp*.
The files *path_sampler.h* and *path_sampler.cpp* provide functionality to trace paths from either the camera or the lights into the scene.
They are used in *cam_trace.h*, *cam_trace.cpp*, *light_trace.h*, and *light_trace.cpp* to compute image contributions from all path sampling techniques.
The files *mis.h* and *mis.cpp* provide code to combine the estimates via multiple importance sampling (MIS).
The derivation of the MIS weights can be found in the [VCM paper](http://www.iliyan.com/publications/VertexMerging), the method used to compute the weights is inspired by the approach described in the [tech report](http://www.iliyan.com/publications/ImplementingVCM) on implementing VCM.
The remaining files define datastructures and utility functions required for the implementation, like the hash grid used for photon mapping and the 4D histograms used for emission guiding.


### Remarks and Known Issues

- Our code requires compiler support for C++11
