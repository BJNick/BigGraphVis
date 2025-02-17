/*
 ==============================================================================

 RPGPUForceAtlas2.cu
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

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <fstream>
#include <chrono>
#include <algorithm>
#include "time.h"
#include <cstdio>
#include <iostream>
#include <math.h>

#include "RPGPUForceAtlas2.hpp"
#include "RPBHFA2LaunchParameters.cuh"
#include "RPBHKernels.cuh"
#include "RPFA2Kernels.cuh"
#include <math_constants.h>

__global__
void PinPolesKernel(int nbodiesd, int npoles, float radius, volatile float2 * __restrict body_posd, volatile int * __restrict poleid)
{
    int i = threadIdx.x;
    if (i >= 0 && i < npoles) {
        body_posd[poleid[i]].x = cosf(2*CUDART_PI_F*i/npoles)*radius;
        body_posd[poleid[i]].y = sinf(2*CUDART_PI_F*i/npoles)*radius;
    }
}

namespace RPGraph
{
    CUDAForceAtlas2::CUDAForceAtlas2(GraphLayout &layout, bool use_barneshut,
                                     bool strong_gravity, float gravity,
                                     float scale)
    : ForceAtlas2(layout, use_barneshut, strong_gravity, gravity, scale)
    {
        int deviceCount;
        cudaGetDeviceCount(&deviceCount);
        if (deviceCount == 0)
        {
            fprintf(stderr, "error: No CUDA devices found.\n");
            exit(EXIT_FAILURE);
        }

        // Host initialization and setup //
        nbodies = layout.graph.num_nodes();
        nedges  = layout.graph.num_edges();
        npoles = layout.pole_list_size;

        body_pos = (float2 *)malloc(sizeof(float2) * layout.graph.num_nodes());
        body_mass = (float *)malloc(sizeof(float) * layout.graph.num_nodes());
        sources  = (int *)  malloc(sizeof(int)   * layout.graph.num_edges());
        targets  = (int *)  malloc(sizeof(int)   * layout.graph.num_edges());
        fx       = (float *)malloc(sizeof(float) * layout.graph.num_nodes());
        fy       = (float *)malloc(sizeof(float) * layout.graph.num_nodes());
        fx_prev  = (float *)malloc(sizeof(float) * layout.graph.num_nodes());
        fy_prev  = (float *)malloc(sizeof(float) * layout.graph.num_nodes());

        poleid = (int *)malloc(sizeof(int) * layout.pole_list_size);

        for (nid_t n = 0; n < layout.graph.num_nodes(); ++n)
        {
            body_pos[n] = {layout.getX(n), layout.getY(n)};
            body_mass[n] = ForceAtlas2::mass(n);
            fx[n] = 0.0;
            fy[n] = 0.0;
            fx_prev[n] = 0.0;
            fy_prev[n] = 0.0;
        }

        for (int i = 0; i < layout.pole_list_size; i++) {
            poleid[i] = layout.graph.node_map[layout.pole_list[i]];
        }

        int cur_sources_idx = 0;
        int cur_targets_idx = 0;

        // Initialize the sources and targets arrays with edge-data.
        for (nid_t source_id = 0; source_id < layout.graph.num_nodes(); ++source_id)
        {
            for (nid_t target_id : layout.graph.neighbors_with_geq_id(source_id))
            {
                sources[cur_sources_idx++] = source_id;
                targets[cur_targets_idx++] = target_id;
            }
        }

        // GPU initialization and setup //
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, 0);

        if (deviceProp.warpSize != WARPSIZE)
        {
            printf("Warpsize of device is %d, but we anticipated %d\n", deviceProp.warpSize, WARPSIZE);
            exit(EXIT_FAILURE);

        }
        cudaFuncSetCacheConfig(BoundingBoxKernel, cudaFuncCachePreferShared);
        cudaFuncSetCacheConfig(TreeBuildingKernel, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(ClearKernel1, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(ClearKernel2, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(SummarizationKernel, cudaFuncCachePreferShared);
        cudaFuncSetCacheConfig(SortKernel, cudaFuncCachePreferL1);
#if __CUDA_ARCH__ < 300
        cudaFuncSetCacheConfig(ForceCalculationKernel, cudaFuncCachePreferL1);
#endif
        cudaFuncSetCacheConfig(DisplacementKernel, cudaFuncCachePreferL1);

        cudaGetLastError();  // reset error value

        // Allocate space on device.
        mp_count = deviceProp.multiProcessorCount;
        max_threads_per_block = deviceProp.maxThreadsPerBlock;

        nnodes = std::max(2 * nbodies, mp_count * max_threads_per_block);

        // Round up to next multiple of WARPSIZE
        while ((nnodes & (WARPSIZE-1)) != 0) nnodes++;
        nnodes--;

        // child stores structure of the quadtree. values point to IDs.
        cudaCatchError(cudaMalloc((void **)&childl,  sizeof(int)   * (nnodes+1) * 4));

        // the following properties, for each node in the quadtree (both internal and leaf)
        cudaCatchError(cudaMalloc((void **)&body_massl,   sizeof(float) * nbodies));
        cudaCatchError(cudaMalloc((void **)&node_massl,   sizeof(float) * (nnodes+1)));
        cudaCatchError(cudaMalloc((void **)&body_posl,sizeof(float2) * nbodies));
        cudaCatchError(cudaMalloc((void **)&node_posl,    sizeof(float2) * (nnodes+1)));
        // count contains the number of nested nodes for each node in quadtree
        cudaCatchError(cudaMalloc((void **)&countl,  sizeof(int)   * (nnodes+1)));
        // start contains ...
        cudaCatchError(cudaMalloc((void **)&startl,  sizeof(int)   * (nnodes+1)));
        cudaCatchError(cudaMalloc((void **)&sortl,   sizeof(int)   * (nnodes+1)));


        cudaCatchError(cudaMalloc((void **)&sourcesl,sizeof(int)   * (nedges)));
        cudaCatchError(cudaMalloc((void **)&targetsl,sizeof(int)   * (nedges)));
        cudaCatchError(cudaMalloc((void **)&fxl,     sizeof(float) * (nbodies)));
        cudaCatchError(cudaMalloc((void **)&fyl,     sizeof(float) * (nbodies)));
        cudaCatchError(cudaMalloc((void **)&fx_prevl,sizeof(float) * (nbodies)));
        cudaCatchError(cudaMalloc((void **)&fy_prevl,sizeof(float) * (nbodies)));
        
        cudaCatchError(cudaMalloc((void **)&poleidl, sizeof(int) * npoles));

        // Used for reduction in BoundingBoxKernel
        cudaCatchError(cudaMalloc((void **)&maxxl,   sizeof(float) * mp_count * FACTOR1));
        cudaCatchError(cudaMalloc((void **)&maxyl,   sizeof(float) * mp_count * FACTOR1));
        cudaCatchError(cudaMalloc((void **)&minxl,   sizeof(float) * mp_count * FACTOR1));
        cudaCatchError(cudaMalloc((void **)&minyl,   sizeof(float) * mp_count * FACTOR1));

        // Used for reduction in SpeedKernel
        cudaCatchError(cudaMalloc((void **)&swgl,    sizeof(float) * mp_count * FACTOR1));
        cudaCatchError(cudaMalloc((void **)&etral,   sizeof(float) * mp_count * FACTOR1));

        // Copy host data to device.
        cudaCatchError(cudaMemcpy(body_massl, body_mass, sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(body_posl,  body_pos,  sizeof(float2) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(sourcesl, sources, sizeof(int) * nedges, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(targetsl, targets, sizeof(int) * nedges, cudaMemcpyHostToDevice));

        // cpy fx, fy , fx_prevl, fy_prevl so they are all initialized to 0 in device memory.
        cudaCatchError(cudaMemcpy(fxl, fx,           sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(fyl, fy,           sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(fx_prevl, fx_prev, sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(fy_prevl, fy_prev, sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        
        cudaCatchError(cudaMemcpy(poleidl, poleid, sizeof(int) * npoles, cudaMemcpyHostToDevice));
    }

    void CUDAForceAtlas2::freeGPUMemory()
    {
        cudaFree(childl);

        cudaFree(body_massl);
        cudaFree(node_massl);
        cudaFree(body_posl);
        cudaFree(node_posl);
        cudaFree(sourcesl);
        cudaFree(targetsl);
        cudaFree(countl);
        cudaFree(startl);
        cudaFree(sortl);

        cudaFree(poleidl);

        cudaFree(fxl);
        cudaFree(fx_prevl);
        cudaFree(fyl);
        cudaFree(fy_prevl);

        cudaFree(maxxl);
        cudaFree(maxyl);
        cudaFree(minxl);
        cudaFree(minyl);

        cudaFree(swgl);
        cudaFree(etral);
    }

    CUDAForceAtlas2::~CUDAForceAtlas2()
    {
        free(body_mass);
        free(body_pos);
        free(sources);
        free(targets);
        free(fx);
        free(fy);
        free(fx_prev);
        free(fy_prev);
        free(poleid);

        freeGPUMemory();
    }

    float CUDAForceAtlas2::count_misaligned_edges(float threshold) {
        return 0;
    }

    void CUDAForceAtlas2::doStep(uint32_t *nodemap)
    {
        float pin_radius;

        if (pin_poles && pole_list_size >= 2) {
            pin_radius = magetic_pole_separation / (2*sin(M_PI/pole_list_size));
            PinPolesKernel<<<1, pole_list_size>>>(nbodies, pole_list_size, pin_radius, body_posl, poleidl);
        }

        GravityKernel<<<mp_count * FACTOR6, THREADS6>>>(nbodies, k_g, strong_gravity, body_massl, body_posl, fxl, fyl);

        AttractiveForceKernel<<<mp_count * FACTOR6, THREADS6>>>(nedges, body_posl, fxl, fyl, sourcesl, targetsl,nodemap);

        BoundingBoxKernel<<<mp_count * FACTOR1, THREADS1>>>(nnodes, nbodies, startl, childl, node_massl, body_posl, node_posl, maxxl, maxyl, minxl, minyl);

        // Build Barnes-Hut Tree
        // 1.) Set all child pointers of internal nodes (in childl) to null (-1)
        ClearKernel1<<<mp_count, 1024>>>(nnodes, nbodies, childl);
        // 2.) Build the tree
        TreeBuildingKernel<<<mp_count * FACTOR2, THREADS2>>>(nnodes, nbodies, childl, body_posl, node_posl);
        // 3.) Set all cell mass values to -1.0, set all startd to null (-1)
        ClearKernel2<<<mp_count, 1024>>>(nnodes, startl, node_massl);

        // Recursively compute mass for each BH. cell.
        SummarizationKernel<<<mp_count * FACTOR3, THREADS3>>>(nnodes, nbodies, countl, childl, body_massl, node_massl, body_posl, node_posl);

        SortKernel<<<mp_count * FACTOR4, THREADS4>>>(nnodes, nbodies, sortl, countl, startl, childl);

        // Compute repulsive forces between nodes using BH. tree.
        ForceCalculationKernel<<<mp_count * FACTOR5, THREADS5>>>(nnodes, nbodies, itolsq, epssq, sortl, childl, body_massl, node_massl, body_posl, node_posl, fxl, fyl, k_r);

        SpeedKernel<<<mp_count * FACTOR1, THREADS1>>>(nbodies, fxl, fyl, fx_prevl, fy_prevl, body_massl, swgl, etral);

        DisplacementKernel<<<mp_count * FACTOR6, THREADS6>>>(nbodies, body_posl, fxl, fyl, fx_prevl, fy_prevl);

        if (pin_poles && pole_list_size >= 2) {
            PinPolesKernel<<<1, pole_list_size>>>(nbodies, pole_list_size, pin_radius, body_posl, poleidl);
        }

        cudaCatchError(cudaDeviceSynchronize());
        iteration++;
    }

    void CUDAForceAtlas2::retrieveLayoutFromGPU()
    {
        cudaCatchError(cudaMemcpy(body_pos, body_posl, sizeof(float2) * nbodies, cudaMemcpyDeviceToHost));
        cudaDeviceSynchronize();
    }

    void CUDAForceAtlas2::sendLayoutToGPU()
    {
        cudaCatchError(cudaMemcpy(body_posl, body_pos, sizeof(float2) * nbodies, cudaMemcpyHostToDevice));
        cudaDeviceSynchronize();
    }

    void CUDAForceAtlas2::sendGraphToGPU()
    {
        cudaCatchError(cudaMemcpy(body_massl, body_mass, sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(sourcesl, sources, sizeof(int) * nedges, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(targetsl, targets, sizeof(int) * nedges, cudaMemcpyHostToDevice));
        cudaDeviceSynchronize();
    }

    void CUDAForceAtlas2::sync_layout()
    {
        retrieveLayoutFromGPU();
        for(nid_t n = 0; n < layout.graph.num_nodes(); ++n)
        {
            layout.setX(n, body_pos[n].x);
            layout.setY(n, body_pos[n].y);
        }
    }
}
