/*
 ==============================================================================

 RPGPUForceAtlas2.hpp
 Copyright © 2016, 2017, 2018  G. Brinkmann

 This file is part of graph_viewer.

 graph_viewer is free software: you can redistribute it and/or modify
 it under the terms of version 3 of the GNU Affero General Public License as
 published by the Free Software Foundation.

 graph_viewer is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with graph_viewer.  If not, see <https://www.gnu.org/licenses/>.

 ==============================================================================
*/

#ifndef RPGPUForceAtlas2_hpp
#define RPGPUForceAtlas2_hpp
#include "RPForceAtlas2.hpp"

__global__
void PinPolesKernel(int nbodiesd, int npoles, float radius, 
                        volatile float2 * __restrict body_posd, 
                        volatile int * __restrict poleid);

namespace RPGraph
{
    class CUDAForceAtlas2: public ForceAtlas2
    {
    public:
        CUDAForceAtlas2(GraphLayout &layout, bool use_barneshut,
                        bool strong_gravity, float gravity, float scale);
        ~CUDAForceAtlas2();
        void doStep(uint32_t *nodemap) override;
        void sync_layout() override;

        float count_misaligned_edges(float threshold) override;

    private:
        /// CUDA Specific stuff.
        // Host storage.
        float *body_mass;
        float2 *body_pos;
        float *fx, *fy, *fx_prev, *fy_prev;

        // Quick way to represent a graph on the GPU
        int *sources, *targets;

        // Pointers to device memory (all suffixed with 'l').
        int   *errl,  *sortl, *childl, *countl, *startl;
        int   *sourcesl, *targetsl;
        float *body_massl, *node_massl;
        float2 *body_posl, *node_posl;
        float *minxl, *minyl, *maxxl, *maxyl;
        float *fxl, *fyl, *fx_prevl, *fy_prevl;
        float *swgl, *etral;

        int *poleid, *poleidl;

        int mp_count; // Number of multiprocessors on GPU.
        int max_threads_per_block;
        int nnodes;
        int nbodies;
        int nedges;

        int npoles;

        void sendGraphToGPU();
        void sendLayoutToGPU();
        void retrieveLayoutFromGPU();
        void freeGPUMemory();

        // New helper methods
        void pinRoots(float2 *body_posl);
    };
};


#endif /* RPGPUForceAtlas2_hpp */
