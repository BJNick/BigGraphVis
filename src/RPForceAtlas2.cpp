/*
 ==============================================================================

 RPForceAtlas2.cpp
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

#include "RPForceAtlas2.hpp"

namespace RPGraph
{
    ForceAtlas2::ForceAtlas2(GraphLayout &layout, bool use_barneshut,
                             bool strong_gravity, float gravity, float scale)
    : LayoutAlgorithm(layout), use_barneshut{use_barneshut},
      strong_gravity{strong_gravity}
    {
        iteration = 0;

        setGravity(gravity);
        setScale(scale);

        global_speed = 1.0;
        speed_efficiency = 1.0;
        jitter_tolerance = 1.0;

        k_s = 0.1;
        k_s_max = 10.0;

        theta = 1.0;
        epssq  = 0.05 * 0.05;
        itolsq = 1.0f / (theta * theta);

        delta = 0.0;

        prevent_overlap = false;
        use_linlog = false;

        layout.randomizePositions();

        // Additional parameters
        use_magnetic_field = false;
        k_attraction = 1;
        attraction_exponent = 1;
        pin_2_roots = false;
        magetic_pole_separation = 10000;
        max_force = 0;
        repulsion_d_squared = false;
        pole_list = {};
        pole_list_size = 0;
    }

    ForceAtlas2::~ForceAtlas2(){};

    void ForceAtlas2::doSteps(int n,uint32_t *nodemap)
    {
        for (int i = 0; i < n; ++i) doStep( nodemap);
    }

    void ForceAtlas2::setScale(float s)
    {
        k_r = s;
    }

    void ForceAtlas2::setGravity(float g)
    {
        k_g = g;
    }

    void ForceAtlas2::setMagneticParameters(std::string field_type, bool bi_directional, float field_strength, float c_m, float alpha, float beta) 
    {
        this->use_magnetic_field = true;
        this->field_type = field_type;
        this->bi_directional = bi_directional;
        this->field_strength = field_strength;
        this->c_m = c_m;
        this->alpha = alpha;
        this->beta = beta;
    }

    float ForceAtlas2::mass(nid_t n)
    {
        return layout.graph.degree(n) + 1.0;
    }
}
