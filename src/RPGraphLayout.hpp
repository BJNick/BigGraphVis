/*
 ==============================================================================

 RPGraphLayout.cpp
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

#ifndef RPGraphLayout_hpp
#define RPGraphLayout_hpp

#include "RPGraph.hpp"
#include "RPCommon.hpp"
#include <string>
#include <unordered_set>

namespace RPGraph
{
    class GraphLayout
    {
    private:
        Coordinate *coordinates;

    protected:
        float width, height;
        float minX(), minY(), maxX(), maxY();

        float node_alpha, edge_alpha; // Additional cosmetic parameters

    public:
        GraphLayout(RPGraph::UGraph &graph,
                    float width = 10000, float height = 10000);
        ~GraphLayout();

        UGraph &graph; // to lay-out

        // randomize the layout position of all nodes.
        void randomizePositions();

        float getX(nid_t node_id), getY(nid_t node_id);
        float getXRange(), getYRange(), getSpan();
        float getDistance(nid_t n1, nid_t n2);
        Real2DVector getDistanceVector(nid_t n1, nid_t n2);
        Real2DVector getNormalizedDistanceVector(nid_t n1, nid_t n2);
        Coordinate getCoordinate(nid_t node_id);
        Coordinate getCenter();


        void setX(nid_t node_id, float x_value), setY(nid_t node_id, float y_value);
        void moveNode(nid_t, Real2DVector v);
        void setCoordinates(nid_t node_id, Coordinate c);
        void writeToPNG(const int image_w, const int image_h, std::string path);
        void writeToCSV(std::string path);
        void writeToBin(std::string path);

        // Additional cosmetic parameters
        bool square_coordinates;
        bool draw_arrows;
        float min_arrow_length;
        void setAlphaParameters(float node_alpha, float edge_alpha);
        int* pole_list;
        int pole_list_size;

        //std::unordered_set<nid_t> connected_to_poles = std::unordered_set<nid_t>();
        std::unordered_set<nid_t>* getConnectedToList();
        void addConnectedNodes(std::unordered_set<nid_t> &connected_nodes, nid_t node);
        void getNodeColor(nid_t n, double &r, double &g, double &b);
        bool isConnectedTo(nid_t node, int pole);
        bool isDisconnected(nid_t node);
        bool isConnectedToTwoPoles(nid_t node);
    };
}

#endif /* RPGraphLayout_hpp */
