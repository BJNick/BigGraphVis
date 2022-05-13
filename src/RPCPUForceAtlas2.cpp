/*
 ==============================================================================

 RPCPUForceAtlas2.cpp
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

#include "RPCPUForceAtlas2.hpp"
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <cmath>
#include <chrono>

#include <stdexcept>
#include <iostream>

namespace RPGraph
{
    // CPUForceAtlas2 definitions.
    CPUForceAtlas2::CPUForceAtlas2(GraphLayout &layout, bool use_barneshut,
                                   bool strong_gravity, float gravity,
                                   float scale)
    :  ForceAtlas2(layout, use_barneshut, strong_gravity, gravity, scale),
       BH_Approximator{layout.getCenter(), layout.getSpan()+10, theta}
    {
        forces      = (Real2DVector *)malloc(sizeof(Real2DVector) * layout.graph.num_nodes());
        prev_forces = (Real2DVector *)malloc(sizeof(Real2DVector) * layout.graph.num_nodes());
        for (nid_t n = 0; n < layout.graph.num_nodes(); ++n)
        {
            forces[n]      = Real2DVector(0.0f, 0.0f);
            prev_forces[n] = Real2DVector(0.0f, 0.0f);
        }
    }

    CPUForceAtlas2::~CPUForceAtlas2()
    {
        free(forces);
        free(prev_forces);
    }

    void CPUForceAtlas2::apply_attract(nid_t n)
    {
        Real2DVector f = Real2DVector(0.0, 0.0);
        for (nid_t t : layout.graph.neighbors_with_geq_id(n))
        {
            // Here we define the magnitude of the attractive force `f_a'
            // *divided* by the length distance between `n' and `t', i.e. `f_a_over_d'
            float f_a_over_d;
            float dist = layout.getDistance(n, t);

            if (use_linlog)
            {
                f_a_over_d = dist == 0.0 ? std::numeric_limits<float>::max() : logf(1+dist) / dist;
            }
            else
            {
                f_a_over_d = 1.0;
            }

            // ADDITIONAL FACTOR: attraction_exponent and k_attraction
            float exp_factor = powf(dist, attraction_exponent-1.0) * k_attraction;
            // Check if not finite
            if (!std::isfinite(exp_factor))
            {
                fprintf(stderr, "exp_factor is not finite, dist = %f\n", dist);
                exit(EXIT_FAILURE);
            }

            f_a_over_d *= exp_factor;

            f += layout.getDistanceVector(n, t) * f_a_over_d;

            //TODO: this is temporary, but required due to
            //      iteration over neighbors_with_geq_id
            forces[t] += layout.getDistanceVector(n, t) * (-f_a_over_d);

    //            forces[n] += getNormalizedDistanceVector(n, t) * f_a(n, t);
        }
        forces[n] += f;
    }

    //====== New magnetic force code ======

    void CPUForceAtlas2::apply_magnetic(nid_t n)
    {
        Real2DVector f = Real2DVector(0.0, 0.0);
        for (nid_t t : layout.graph.neighbors_with_geq_id(n))
        {
            // Here we define the magnitude of the attractive force `f_a'
            // *divided* by the length distance between `n' and `t', i.e. `f_a_over_d'
            Real2DVector pos_n = layout.getCoordinate(n).toVector();
            Real2DVector pos_t = layout.getCoordinate(t).toVector();
            Real2DVector disp = layout.getDistanceVector(n, t);
            float dist = std::sqrt(disp.magnitude());

            Real2DVector field_direction = get_magnetic_field(center_of_mass(n, t));

            // Assume the graph is directed from lower to higher node ids
            int edge_dir = layout.graph.node_map_r[n] < layout.graph.node_map_r[t] ? 1 : -1;
            field_direction = field_direction * edge_dir;

            if (dist == 0.0 || field_direction.magnitude() == 0.0)
                continue; // Cannot compute the angle when either is zero

            Real2DVector force_on_n = magnetic_equation(field_direction, disp, field_strength, c_m, alpha, beta);
            
            Real2DVector force_on_t = force_on_n * -1;

            Real2DVector accel_on_n = force_on_n / mass(n);
            Real2DVector accel_on_t = force_on_t / mass(t);

            f += accel_on_n;
            forces[t] += accel_on_t;

            // Error case, something went wrong with the above computation
            if (!std::isfinite(accel_on_n.magnitude()) || !std::isfinite(accel_on_t.magnitude())) {
                std::cout << "Magnetic acceleration is not finite: " << accel_on_n.magnitude() << "\n";
                exit(EXIT_FAILURE);
            }
        }
        forces[n] += f;
    }

    Real2DVector CPUForceAtlas2::center_of_mass(nid_t n, nid_t t)
    {
        Real2DVector pos_n = layout.getCoordinate(n).toVector();
        Real2DVector pos_t = layout.getCoordinate(t).toVector();
        float mass_n = mass(n) == 0 ? 1 : mass(n);
        float mass_t = mass(t) == 0 ? 1 : mass(t);
        float one_over_total_mass = 1 / (mass_n + mass_t);
        Real2DVector center_of_mass = (pos_n*mass_n + pos_t*mass_t) * one_over_total_mass;
        return center_of_mass;
    }

    Real2DVector CPUForceAtlas2::magnetic_equation(Real2DVector m, Real2DVector d, float b, float c, float alpha, float beta) 
    {
        float dist = std::sqrt(d.magnitude());
        Real2DVector force_on_n = Real2DVector(0.0, 0.0);
        if (!bi_directional)
            force_on_n = d.rotate90clockwise() * -(b * c * powf(dist, alpha-1) * powf(m.angleCos(d), beta) * sign(m.cross(d)));
        else
            force_on_n = d.rotate90clockwise() * -(b * c * powf(dist, alpha-1) * powf(std::abs(m.angleSin(d)), beta) * sign(m.cross(d) * m.dot(d)));
        return force_on_n;
    }
    
    Real2DVector CPUForceAtlas2::get_magnetic_field(Real2DVector pos) 
    {
        // TODO: Make coordinates configurable
        Real2DVector pole1 = Real2DVector(-magetic_pole_separation/2, 0); // positive (from)
        Real2DVector pole2 = Real2DVector(magetic_pole_separation/2, 0);  // negative (to)
        if (field_type == "none") {
            return Real2DVector(0, 0);
        } else if (field_type == "parallel" || field_type == "linear") {
            return Real2DVector(1, 0);
        } else if (field_type == "polar") {
            return pos.getNormalizedFinite();
        } else if (field_type == "polar-reversed") {
            return pos.getNormalizedFinite() * (-1);
        } else if (field_type == "concentric") {
            return pos.getNormalizedFinite().rotate90clockwise();
        } else if (field_type == "dipole-1" || field_type == "dipole-2") {
            // Dipole 1 has one the denominator with power of 1, and the other with power of 2
            Real2DVector v1 = pos - pole1;
            Real2DVector v2 = pole2 - pos;
            float dist1 = v1.magnitude();
            float dist2 = v2.magnitude();
            if (field_type == "dipole-1") {
                float m1 = 1 / (dist1);
                float m2 = 1 / (dist2);
                return (v1 * m1 + v2 * m2).getNormalizedFinite();
            } else { // field_type == "dipole-2"
                float m1 = 1 / (dist1 * dist1);
                float m2 = 1 / (dist2 * dist2);
                return (v1 * m1 + v2 * m2).getNormalizedFinite();
            }
        } else if (field_type == "negative-charges") {
            // Two negative charges at the dipole positions
            Real2DVector v1 = pole1 - pos; // Reversed from dipole
            Real2DVector v2 = pole2 - pos;
            float dist1 = v1.magnitude();
            float dist2 = v2.magnitude();
            float m1 = 1 / (dist1 * dist1);
            float m2 = 1 / (dist2 * dist2);
            return (v1 * m1 + v2 * m2).getNormalizedFinite();
        } else if (field_type == "linear-dipole") {
            // Looks like: <-<-<- | ->->-> | <-<-<- 
            if (pos.x < pole1.x) {
                return Real2DVector(-1, 0);
            } else if (pos.x < pole2.x) {
                return Real2DVector(1, 0);
            } else {
                return Real2DVector(-1, 0);
            }
        } else {
            throw std::invalid_argument("Invalid magnetic field type: " + field_type);
        }
    }

    //====================================

    void CPUForceAtlas2::apply_repulsion(nid_t n)
    {
        if (use_barneshut)
        {
            forces[n] += (BH_Approximator.approximateForce(layout.getCoordinate(n), mass(n), theta) * k_r);
        }

        else
        {
            for (nid_t t = 0; t < layout.graph.num_nodes(); ++t)
            {
                if (n == t) continue;
                float  distance = layout.getDistance(n, t);
                float f_r = distance == 0.0 ? std::numeric_limits<float>::max() : k_r * mass(n) * mass(t) / distance / distance;
                forces[n] += layout.getDistanceVector(n, t) * f_r;
            }
        }
    }

    void CPUForceAtlas2::apply_gravity(nid_t n)
    {
        float f_g, d;

        // `d' is the distance from `n' to the center (0.0, 0.0)
        d = std::sqrt(layout.getX(n)*layout.getX(n) + layout.getY(n)*layout.getY(n));
        if(d == 0.0) return;

        // Here we define the magnitude of the gravitational force `f_g'.
        if (strong_gravity)
        {
            f_g = k_g*mass(n);
        }

        else
        {
            f_g = k_g*mass(n) / d; 
        }

        forces[n] += (Real2DVector(-layout.getX(n), -layout.getY(n)) * f_g);
    }

    // Eq. (8)
    float CPUForceAtlas2::swg(nid_t n)
    {
        return (forces[n] - prev_forces[n]).magnitude();
    }

    // Eq. (9)
    float CPUForceAtlas2::s(nid_t n)
    {
        return (k_s * global_speed)/(1.0f+global_speed*std::sqrt(swg(n)));
    }

    // Eq. (12)
    float CPUForceAtlas2::tra(nid_t n)
    {
        return (forces[n] + prev_forces[n]).magnitude() / 2.0;
    }

    void CPUForceAtlas2::updateSpeeds()
    {
        // The following speed-update procedure for ForceAtlas2 follows
        // the one by Gephi:
        // https://github.com/gephi/gephi/blob/6efb108718fa67d1055160f3a18b63edb4ca7be2/modules/LayoutPlugin/src/main/java/org/gephi/layout/plugin/forceAtlas2/ForceAtlas2.java

        // `Auto adjust speeds'
        float total_swinging = 0.0;
        float total_effective_traction = 0.0;
        for (nid_t nid = 0; nid < layout.graph.num_nodes(); ++nid)
        {
            total_swinging += mass(nid) * swg(nid); // Eq. (11)
            total_effective_traction += mass(nid) * tra(nid); // Eq. (13)
        }

        // We want to find the right jitter tollerance for this graph,
        // such that totalSwinging < tolerance * totalEffectiveTraction

        float estimated_optimal_jitter_tollerance = 0.05 * std::sqrt(layout.graph.num_nodes());
        float minJT = std::sqrt(estimated_optimal_jitter_tollerance);
        float jt = jitter_tolerance * fmaxf(minJT,
                                           fminf(k_s_max,
                                                 estimated_optimal_jitter_tollerance * total_effective_traction / powf(layout.graph.num_nodes(), 2.0)
                                                 )
                                           );
        float min_speed_efficiency = 0.05;

        // `Protect against erratic behavior'
        if (total_swinging / total_effective_traction > 2.0)
        {
            if (speed_efficiency > min_speed_efficiency) speed_efficiency *= 0.5;
            jt = fmaxf(jt, jitter_tolerance);
        }

        // `Speed efficiency is how the speed really corrosponds to the swinging vs. convergence tradeoff.'
        // `We adjust it slowly and carefully'
        float targetSpeed = jt * speed_efficiency * total_effective_traction / total_swinging;

        if (total_swinging > jt * total_effective_traction)
        {
            if (speed_efficiency > min_speed_efficiency)
            {
                speed_efficiency *= 0.7;
            }
        }
        else if (global_speed < 1000)
        {
            speed_efficiency *= 1.3;
        }

        // `But the speed shouldn't rise much too quickly, ... would make convergence drop dramatically'.
        float max_rise = 0.5;
        global_speed += fminf(targetSpeed - global_speed, max_rise * global_speed);
    }

    void CPUForceAtlas2::apply_displacement(nid_t n)
    {
        if (prevent_overlap)
        {
            // Not yet implemented
            exit(EXIT_FAILURE);
        }

        else
        {

            float factor = global_speed / (1.0 + std::sqrt(global_speed * swg(n)));
            layout.moveNode(n, forces[n] * factor);
        }
    }

    void CPUForceAtlas2::rebuild_bh()
    {
        BH_Approximator.repulsion_d_squared = repulsion_d_squared;
        BH_Approximator.reset(layout.getCenter(), layout.getSpan()+10);

        for (nid_t n = 0; n < layout.graph.num_nodes(); ++n)
        {
            BH_Approximator.insertParticle(layout.getCoordinate(n),
                                           layout.graph.degree(n)+1);
        }
    }

    void CPUForceAtlas2::doStep(uint32_t *nodemap)
    {
        // Context
        // std::unordered_map<nid_t, nid_t> node_map; // el id -> UGraph id
        // std::unordered_map<nid_t, nid_t> node_map_r; // UGraph id -> el id

        if (pin_2_roots)
        {
            nid_t r1 = layout.graph.node_map[1];
            nid_t r2 = layout.graph.node_map[layout.graph.num_nodes()];
            layout.setX(r1, -magetic_pole_separation/2);
            layout.setY(r1, 0);
            layout.setX(r2, magetic_pole_separation/2);
            layout.setY(r2, 0);
        }

        if (use_barneshut) rebuild_bh();

        max_force = 0.0;

        for (nid_t n = 0; n < layout.graph.num_nodes(); ++n)
        {
            apply_gravity(n);
            apply_attract(n);
            apply_repulsion(n);
            if (use_magnetic_field)
                apply_magnetic(n);
            max_force = std::max(max_force, forces[n].magnitude());
        }

        updateSpeeds();

        for (nid_t n = 0; n < layout.graph.num_nodes(); ++n)
        {
            apply_displacement(n);
            prev_forces[n]  = forces[n];
            forces[n]       = Real2DVector(0.0f, 0.0f);
        }
        iteration++;

        if (pin_2_roots)
        {
            nid_t r1 = layout.graph.node_map[1];
            nid_t r2 = layout.graph.node_map[layout.graph.num_nodes()];
            layout.setX(r1, -magetic_pole_separation/2);
            layout.setY(r1, 0);
            layout.setX(r2, magetic_pole_separation/2);
            layout.setY(r2, 0);
        }
    }

    void CPUForceAtlas2::sync_layout() {}

}
