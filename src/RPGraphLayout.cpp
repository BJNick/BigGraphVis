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

#include <vector>

#include "RPGraphLayout.hpp"
#include "../lib/pngwriter/src/pngwriter.h"
#include <algorithm>
#include <fstream>
#include <cmath>
#include <limits>
#include <iostream>
#include <bits/stdc++.h>

#include <queue>
#include <unordered_set>
#include <unordered_map>

namespace RPGraph
{
    GraphLayout::GraphLayout(UGraph &graph, float width, float height)
        : graph(graph), width(width), height(height)
    {
        coordinates = (Coordinate *)malloc(graph.num_nodes() * sizeof(Coordinate));
        node_alpha = 0.8;
        edge_alpha = 0.005;
        square_coordinates = false;
        draw_arrows = false;
        min_arrow_length = 50;
        draw_common_edges = true;
        use_distance_based_edge_direction = false;
        max_influence_distance = -1; // -1 means no limit
        pole_size_factor = 3;
        colored_fraction = 1;
        predraw_edges = false;
    }

    GraphLayout::~GraphLayout()
    {
        free(coordinates);
    }

    void GraphLayout::randomizePositions()
    {
        for (nid_t i = 0; i < graph.num_nodes(); ++i)
        {
            setX(i, get_random(-width / 2.0, width / 2.0));
            setY(i, get_random(-height / 2.0, height / 2.0));
        }
    }

    float GraphLayout::getX(nid_t node_id)
    {
        return coordinates[node_id].x;
    }

    float GraphLayout::getY(nid_t node_id)
    {
        return coordinates[node_id].y;
    }

    float GraphLayout::minX()
    {
        float minX = std::numeric_limits<float>::max();
        for (nid_t n = 0; n < graph.num_nodes(); ++n)
            if (getX(n) < minX)
                minX = getX(n);
        return minX;
    }

    float GraphLayout::maxX()
    {
        float maxX = std::numeric_limits<float>::lowest();
        for (nid_t n = 0; n < graph.num_nodes(); ++n)
            if (getX(n) > maxX)
                maxX = getX(n);
        return maxX;
    }

    float GraphLayout::minY()
    {
        float minY = std::numeric_limits<float>::max();
        for (nid_t n = 0; n < graph.num_nodes(); ++n)
            if (getY(n) < minY)
                minY = getY(n);
        return minY;
    }

    float GraphLayout::maxY()
    {
        float maxY = std::numeric_limits<float>::lowest();
        for (nid_t n = 0; n < graph.num_nodes(); ++n)
            if (getY(n) > maxY)
                maxY = getY(n);
        return maxY;
    }

    float GraphLayout::getXRange()
    {
        return maxX() - minX();
    }

    float GraphLayout::getYRange()
    {
        return maxY() - minY();
    }

    float GraphLayout::getSpan()
    {
        return ceil(fmaxf(getXRange(), getYRange()));
    }

    float GraphLayout::getDistance(nid_t n1, nid_t n2)
    {
        const float dx = getX(n1) - getX(n2);
        const float dy = getY(n1) - getY(n2);
        return std::sqrt(dx * dx + dy * dy);
    }

    Real2DVector GraphLayout::getDistanceVector(nid_t n1, nid_t n2)
    {
        return Real2DVector(getX(n2) - getX(n1), getY(n2) - getY(n1));
    }

    Real2DVector GraphLayout::getNormalizedDistanceVector(nid_t n1, nid_t n2)
    {
        const float x1 = getX(n1);
        const float x2 = getX(n2);
        const float y1 = getY(n1);
        const float y2 = getY(n2);
        const float dx = x2 - x1;
        const float dy = y2 - y1;
        const float len = std::sqrt(dx * dx + dy * dy);

        return Real2DVector(dx / len, dy / len);
    }

    Coordinate GraphLayout::getCoordinate(nid_t node_id)
    {
        return coordinates[node_id];
    }

    Coordinate GraphLayout::getCenter()
    {
        float x = minX() + getXRange() / 2.0;
        float y = minY() + getYRange() / 2.0;
        return Coordinate(x, y);
    }

    void GraphLayout::setX(nid_t node_id, float x_value)
    {

        coordinates[node_id].x = x_value;
    }

    void GraphLayout::setY(nid_t node_id, float y_value)
    {
        coordinates[node_id].y = y_value;
    }

    void GraphLayout::moveNode(nid_t n, RPGraph::Real2DVector v)
    {
        setX(n, getX(n) + v.x);
        setY(n, getY(n) + v.y);
    }

    void GraphLayout::setCoordinates(nid_t node_id, Coordinate c)
    {
        setX(node_id, c.x);
        setY(node_id, c.y);
    }

    void GraphLayout::writeToPNG(const int image_w, const int image_h, std::string path)
    {

        float xRange = getXRange();
        float yRange = getYRange();
        if (square_coordinates) {
            xRange = yRange = std::max(xRange, yRange);
        }
        const RPGraph::Coordinate center = getCenter();
        const float xCenter = center.x;
        const float yCenter = center.y;
        const float minX = xCenter - xRange / 2.0;
        const float minY = yCenter - yRange / 2.0;
        const float xScale = image_w / xRange;
        const float yScale = image_h / yRange;

        // Here we need to do some guessing as to what the optimal
        // opacity of nodes and edges might be, given network size.
        float node_opacity = 0.5; // 10000.0  / graph.num_nodes();
        const float edge_opacity = 1000.0 / graph.num_edges();

        // Write to file.
        pngwriter layout_png(image_w, image_h, 0, path.c_str());
        layout_png.invert(); // set bg. to white

        int *sorteddegree = new int[graph.num_nodes()];
        int *largest = new int[11];
        int *sortof12 = new int[graph.num_nodes()];
        double r, g, b;

        for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1)
            sorteddegree[n1] = graph.degree(n1);

        int n = graph.num_nodes();

        std::ofstream out_color("../../../files/outcolor.txt");
        int largest_counter = 0;

        std::vector<size_t> idx(graph.num_nodes());
        std::iota(idx.begin(), idx.end(), 0);

        stable_sort(idx.begin(), idx.end(),
                    [&sorteddegree](size_t i1, size_t i2)
                    { return sorteddegree[i1] < sorteddegree[i2]; });
        int sum_degre = 0;
        int ct = 0;
        for (auto i : idx)
        {
            if (ct >= graph.num_nodes() / 2)
                break;
            sum_degre = sqrt(graph.degree(i)) + sum_degre;
            ct++;
        }
        int space = sum_degre;
        space = (space * 0.1);

        int c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0, c7 = 0, c8 = 0, c9 = 0, c10 = 0, c11 = 0;

        ct = 0;
        int sum_track = 0;
        space = sum_degre * .0001;
        // Do not draw edges twice if not necessary
        if (predraw_edges)
            for (auto i : idx)
            {

                for (nid_t n2 : graph.neighbors_with_geq_id(i))
                {
                    // ... and edge.
                    layout_png.line_blend((getX(i) - minX) * xScale, (getY(i) - minY) * yScale,
                                        (getX(n2) - minX) * xScale, (getY(n2) - minY) * yScale,
                                        0.08, 0.3, 0.3, 0.3);
                }
            }
        for (auto i : idx)
        {
            r = double(2) / double(255);
            g = double(10) / double(255);
            b = double(40) / double(255);
            // std::cout<<sorteddegree[i]<<'\n';
            sum_track = sqrt(graph.degree(i));
            if (sum_track <= space)
            {

                /* //brown
                 r=double(177)/double(255);
                 g=double(89)/double(255);
                 b=double(40)/double(255);
                 c1++;
              }
                  else if(sum_track>space&&sum_track<=space*2)
              {
                   //light purpile
                 r=double(202)/double(255);
                 g=double(178)/double(255);
                 b=double(214)/double(255);
                      c2++;
                                    //   std::cout<<r<<" "<<g<<" "<<b<<" 1 "<<'\n';

              }
                  else if(sum_track>space*2&&sum_track<=space*3)
              {


                  //purpule
                 r=double(106)/double(255);
                 g=double(61)/double(255);
                 b=double(254)/double(255);
                                    //  std::cout<<r<<" "<<g<<" "<<b<<" 2 "<<'\n';
 c3++;
              }
                else if(sum_track>space*3&&sum_track<=space*4)
              {
                    //light orange

                 r=double(253)/double(255);
                 g=double(191)/double(255);
                 b=double(111)/double(255);           //     std::cout<<r<<" "<<g<<" "<<b<<" 3 "<<'\n';
 c4++;
              }
          else if(sum_track>space*4&&sum_track<=space*5)
              {


                     //orange
                     r=double(255)/double(255);
                 g=double(127)/double(255);
                 b=double(0)/double(255);
                    c5++;             //  std::cout<<r<<" "<<g<<" "<<b<<" 4 "<<'\n';

              }
          else if(sum_track>space*5&&sum_track<=space*6)
              {
              //light red
                  r=double(251)/double(255);
                 g=double(154)/double(255);
                 b=double(153)/double(255);
                        c6++;           //  std::cout<<r<<" "<<g<<" "<<b<<" 5 "<<'\n';

              }
             else if(sum_track>space*6&&sum_track<=space*7)
              {


                 //red
                 r=double(227)/double(255);
                 g=double(26)/double(255);
                 b=double(28)/double(255);
                   c7++;            //      std::cout<<r<<" "<<g<<" "<<b<<" 6 "<<'\n';

              }
              else if(sum_track>space*7&&sum_track<=space*8)
              {
                  //light green
                  r=double(178)/double(255);
                 g=double(233)/double(255);
                 b=double(138)/double(255);
                              //        std::cout<<r<<" "<<g<<" "<<b<<" 7 "<<'\n';
 c8++;
              }
            else if(sum_track>space*8&&sum_track<=space*9)
              {



                 //green
                 r=double(52)/double(255);
                 g=double(160)/double(255);
                 b=double(44)/double(255);
      c9++;                          //      std::cout<<r<<" "<<g<<" "<<b<<" 8 "<<'\n';

              }
                  else if(sum_track>space*9&&sum_track<=space*10)
              {
                  //light blue
                 r=double(166)/double(255);
                 g=double(206)/double(255);
                 b=double(227)/double(255);
                               //       std::cout<<r<<" "<<g<<" "<<b<<" 9 "<<'\n';
 c10++;
              }


                  else
              {
                  c11++;
                if(sum_track>space*10&&sum_track<=space*20)  {r=double(31)/double(255);
                 g=double(120)/double(255);
                 b=double(180)/double(255);
                 std::cout<<"31";
                 }    //dark blue


                 else if(sum_track>space*20&&sum_track<=space*30) {r=double(60)/double(255);
                 g=double(120)/double(255);
                 b=double(180)/double(255);
                 std::cout<<"60";}

                else  if(sum_track>space*30&&sum_track<=space*1000) {r=double(90)/double(255);
                 g=double(120)/double(255);
                 b=double(180)/double(255);
                 std::cout<<"90";}

                 else{r=double(120)/double(255);
                 g=double(120)/double(255);
                 b=double(180)/double(255);
                 std::cout<<"120";}*/
                // brown
                r = double(247) / double(255);
                g = double(251) / double(255);
                b = double(255) / double(255);
                c1++;
            }
            else if (sum_track > space && sum_track <= space * 2)
            {
                // light purpile
                r = double(222) / double(255);
                g = double(235) / double(255);
                b = double(247) / double(255);
                c2++;
                //   std::cout<<r<<" "<<g<<" "<<b<<" 1 "<<'\n';
            }
            else if (sum_track > space * 2 && sum_track <= space * 3)
            {

                // purpule
                r = double(198) / double(255);
                g = double(219) / double(255);
                b = double(239) / double(255);
                //  std::cout<<r<<" "<<g<<" "<<b<<" 2 "<<'\n';
                c3++;
            }
            else if (sum_track > space * 3 && sum_track <= space * 4)
            {
                // light orange

                r = double(158) / double(255);
                g = double(202) / double(255);
                b = double(225) / double(255); //     std::cout<<r<<" "<<g<<" "<<b<<" 3 "<<'\n';
                c4++;
            }
            else if (sum_track > space * 4 && sum_track <= space * 5)
            {

                // orange
                r = double(107) / double(255);
                g = double(174) / double(255);
                b = double(214) / double(255);
                c5++; //  std::cout<<r<<" "<<g<<" "<<b<<" 4 "<<'\n';
            }
            else if (sum_track > space * 5 && sum_track <= space * 6)
            {
                // light red
                r = double(66) / double(255);
                g = double(146) / double(255);
                b = double(198) / double(255);
                c6++; //  std::cout<<r<<" "<<g<<" "<<b<<" 5 "<<'\n';
            }
            else if (sum_track > space * 6 && sum_track <= space * 7)
            {

                // red
                r = double(33) / double(255);
                g = double(113) / double(255);
                b = double(181) / double(255);
                c7++; //      std::cout<<r<<" "<<g<<" "<<b<<" 6 "<<'\n';
            }
            else if (sum_track > space * 7 && sum_track <= space * 8)
            {
                // light green
                r = double(8) / double(255);
                g = double(81) / double(255);
                b = double(156) / double(255);
                //        std::cout<<r<<" "<<g<<" "<<b<<" 7 "<<'\n';
                c8++;
            }
            else if (sum_track > space * 8 && sum_track <= space * 9)
            {

                // green
                r = double(8) / double(255);
                g = double(48) / double(255);
                b = double(107) / double(255);
                c9++; //      std::cout<<r<<" "<<g<<" "<<b<<" 8 "<<'\n';
            }
            else if (sum_track > space * 9 && sum_track <= space * 10)
            {
                // light blue
                r = double(5) / double(255);
                g = double(24) / double(255);
                b = double(80) / double(255);
                //       std::cout<<r<<" "<<g<<" "<<b<<" 9 "<<'\n';
                c10++;
            }

            int radian = (sqrt(graph.degree(i) + 1)) * 0.6; // sqrt(graph.degree(i));
            
            // TODO: Make only for trees
            // Make the root(s) of a tree appear bigger
            // The first node is the root and nodes with the same degree are also roots

            //std::cout << graph.degree(i) << "\n";
            /*if (i == 0 || graph.degree(i) == graph.degree(0))
            {
                int adj_degree = sorteddegree[0];
                radian = ((sqrt(adj_degree + 1)) * 0.6) * 2;
            }*/

            int threash = 100;
            if (r == double(177) / double(255))
                radian = 1;

            if (pole_list_size > 0)
                getNodeColor(i, r, g, b);

            // Paint roots/poles red
            for (int j = 0; j < pole_list_size; j++)
                if (i == graph.node_map[pole_list[j]])
                {
                    r = 1.0, g = 0.0, b = 0.0; 
                    radian = radian * pole_size_factor;
                    break;
                }
            

            // if(radian<=threash)//sqrt(graph.degree(i)<10))
            /*char  c_string[]= "~/OpenSans-Bold.ttf";
            char c_string2[]= "7";*/
            //    if(graph.node_map_r[i]==748264||graph.node_map_r[i]==72606)
            layout_png.filledcircle_blend((getX(i) - minX) * xScale,
                                          (getY(i) - minY) * yScale,
                                          std::min(threash, radian), node_alpha, r, g, b);
            /*layout_png.plot_text_utf8_blend(c_string, 5,

                                (getX(i) - minX)*xScale, (getY(i) - minY)*yScale, 0.0, c_string2,
                                1.0,

                                0.0,0.0 ,0.0 );*/

            for (nid_t n2 : graph.neighbors_with_geq_id(i))
            {
                if (pole_list_size > 0)
                    getNodeColor(primary(i, n2), r, g, b);

                if (pole_list_size > 0 && getEdgeDirection(i, n2) == 0 && !sameColor(i, n2))
                    r = 0.5, g = 0.5, b = 0.5;

                // Do not paint if the node is not connected to the root
                if (!draw_common_edges && (isDisconnected(primary(i, n2)) || isConnectedToTwoPoles(primary(i, n2))))
                    continue;

                // Draw a line from node to node
                if (colored_fraction <= 0.5) {
                    float distance = sqrt(pow(getX(i) - getX(n2), 2) + pow(getY(i) - getY(n2), 2))*std::min(xScale, yScale);
                    int edge_dir = getEdgeDirection(i, n2);

                    float third1X = (getX(i)*(1-colored_fraction) + getX(n2)*colored_fraction);
                    float third1Y = (getY(i)*(1-colored_fraction) + getY(n2)*colored_fraction);
                    float third2X = (getX(i)*colored_fraction + getX(n2)*(1-colored_fraction));
                    float third2Y = (getY(i)*colored_fraction + getY(n2)*(1-colored_fraction));

                    double r1 = r, g1 = g, b1 = b;
                    double r2 = r, g2 = g, b2 = b;

                    if (edge_dir != 0) {
                        if (edge_dir < 0)
                            r2 = 0.5, g2 = 0.5, b2 = 0.5;
                        if (edge_dir > 0)
                            r1 = 0.5, g1 = 0.5, b1 = 0.5;
                    }
                    layout_png.line_blend((getX(i) - minX) * xScale, (getY(i) - minY) * yScale, (third1X - minX) * xScale, (third1Y - minY) * yScale, edge_alpha, r1, g1, b1);
                    layout_png.line_blend((getX(n2) - minX) * xScale, (getY(n2) - minY) * yScale, (third2X - minX) * xScale, (third2Y - minY) * yScale, edge_alpha, r2, g2, b2);
                    if (colored_fraction < 0.5)
                        layout_png.line_blend((third1X - minX) * xScale, (third1Y- minY) * yScale, (third2X - minX) * xScale, (third2Y - minY) * yScale, edge_alpha, 0.5, 0.5, 0.5);

                } else {
                    layout_png.line_blend((getX(i) - minX) * xScale, (getY(i) - minY) * yScale, (getX(n2) - minX) * xScale, (getY(n2) - minY) * yScale, edge_alpha, r, g, b);
                }

                if (draw_arrows) {
                    // Draw arrows if the edge is long enough
                    float distance = sqrt(pow(getX(i) - getX(n2), 2) + pow(getY(i) - getY(n2), 2))*std::min(xScale, yScale);
                    if (distance > min_arrow_length)
                    {
                        // Draw an arrow head halway between the nodes
                        int edge_dir = getEdgeDirection(i, n2);
                        if (edge_dir == 0) continue;

                        float nudge = edge_dir * 7.0 / distance + 0.5;
                        float midpoint1X = (getX(i)*nudge + getX(n2)*(1-nudge));
                        float midpoint1Y = (getY(i)*nudge + getY(n2)*(1-nudge));
                        float midpoint2X = (getX(i)*(1-nudge) + getX(n2)*nudge);
                        float midpoint2Y = (getY(i)*(1-nudge) + getY(n2)*nudge);

                        layout_png.filledarrow((midpoint1X - minX)*xScale, (midpoint1Y - minY)*yScale, (midpoint2X - minX)*xScale, (midpoint2Y - minY)*yScale, 15, 0.2, r, g, b);
                    }
                }
            }

            nid_t id = graph.node_map_r[i]; // id as found in edgelist
            out_color << id << " " << r << " " << g << " " << b << "\n";

            // ct++;
        }

        // Write it to disk.
        layout_png.write_png();
        out_color.close();
    }

    void GraphLayout::writeToCSV(std::string path)
    {
        if (is_file_exists(path.c_str()))
        {
            printf("Error: File exists at %s\n", path.c_str());
            exit(EXIT_FAILURE);
        }

        std::ofstream out_file(path);

        for (nid_t n = 0; n < graph.num_nodes(); ++n)
        {
            nid_t id = graph.node_map_r[n]; // id as found in edgelist
            out_file << id << "," << getX(n) << "," << getY(n) << "\n";
        }

        out_file.close();
    }

    void GraphLayout::writeToNET(std::string path, std::unordered_map<long, std::string> node_labels)
    {
        if (is_file_exists(path.c_str()))
        {
            printf("Error: File exists at %s\n", path.c_str());
            exit(EXIT_FAILURE);
        }

        std::ofstream out_file(path);

        // Find min/max reverse id
        nid_t max_id = 0;
        nid_t min_id = INT_MAX;
        for (nid_t n = 0; n < graph.num_nodes(); ++n)
        {
            nid_t id = graph.node_map_r[n]; // id as found in edgelist
            if (id > max_id)
                max_id = id;
            if (id < min_id)
                min_id = id;
        }
        int vertex_count = vertex_count = max_id - min_id + 1;

        out_file << "*Vertices " << vertex_count << "\n";

        for (nid_t n = 0; n < graph.num_nodes(); ++n)
        {
            nid_t id = graph.node_map_r[n]; // id as found in edgelist
            std::string label = "";
            if (node_labels.find(id) != node_labels.end())
                label = '\"' + node_labels[id] + '\"';
            else 
                label = '\"' + std::to_string(id) + '\"';
            out_file << id << " " << label << " " << getX(n) << " " << getY(n) << "\n";
        }

        // Now output the arcs (directed)

        out_file << "*Arcs\n";

        for (nid_t n = 0; n < graph.num_nodes(); ++n)
        {
            for (nid_t n2 : graph.neighbors_with_geq_id(n))
            {
                if (getEdgeDirection(n, n2) > 0)
                    out_file << graph.node_map_r[n] << " " << graph.node_map_r[n2] << " " << "\n";
                else if (getEdgeDirection(n, n2) < 0)
                    out_file << graph.node_map_r[n2] << " " << graph.node_map_r[n] << " " << "\n";
            }
        }

        // Finally output the edges (undirected)

        out_file << "*Edges\n";

        for (nid_t n = 0; n < graph.num_nodes(); ++n)
        {
            for (nid_t n2 : graph.neighbors_with_geq_id(n))
            {
                if (getEdgeDirection(n, n2) == 0) {
                    // Print both directions
                    out_file << graph.node_map_r[n] << " " << graph.node_map_r[n2] << " " << "\n";
                    out_file << graph.node_map_r[n2] << " " << graph.node_map_r[n] << " " << "\n";
                }
            }
        }
        
        out_file.close();
    }

    void GraphLayout::writeToBin(std::string path)
    {
        if (is_file_exists(path.c_str()))
        {
            printf("Error: File exists at %s\n", path.c_str());
            exit(EXIT_FAILURE);
        }

        std::ofstream out_file(path, std::ofstream::binary);

        for (nid_t n = 0; n < graph.num_nodes(); ++n)
        {
            nid_t id = graph.node_map_r[n]; // id as found in edgelist
            float x = getX(n);
            float y = getY(n);

            out_file.write(reinterpret_cast<const char *>(&id), sizeof(id));
            out_file.write(reinterpret_cast<const char *>(&x), sizeof(x));
            out_file.write(reinterpret_cast<const char *>(&y), sizeof(y));
        }

        out_file.close();
    }

    // Additional cosmetic parameters
    void GraphLayout::setAlphaParameters(float node_alpha, float edge_alpha)
    {
        this->node_alpha = node_alpha;
        this->edge_alpha = edge_alpha;
    }
    
    // Contains a singleton instance of the connected sets
    // No longer used
    /*std::unordered_set<nid_t>* GraphLayout::getConnectedToList() {
        static std::unordered_set<nid_t>* connected_to = new std::unordered_set<nid_t>[pole_list_size];
        static bool initialized = false;
        if (!initialized) {
            for (int i = 0; i < pole_list_size; i++)
            {
                addConnectedNodes(connected_to[i], graph.node_map[pole_list[i]]);
            }
            initialized = true;
        }
        return connected_to;
    }*/

    nid_t GraphLayout::primary(nid_t n, nid_t t)
    {
        if (getEdgeDirection(n, t) >= 0)
            return t;
        else
            return n;
    }

    nid_t GraphLayout::secondary(nid_t n, nid_t t)
    {
        if (primary(n, t) == t)
            return n;
        else
            return t;
    }

    int GraphLayout::getEdgeDirection(nid_t n, nid_t t)
    {
        if (use_distance_based_edge_direction)
            return getDistanceBasedEdgeDirection(n, t);
        else
            return getInitialEdgeDirection(n, t);
    }

    int GraphLayout::getInitialEdgeDirection(nid_t n, nid_t t)
    {
        if (!graph.is_edge_directed[n][t])
            return 0;
        else if (graph.initial_edge_direction[n][t])
            return 1;
        else
            return -1;
    }

    // Is the given node connected to the given pole? uses connected_to to determine
    bool GraphLayout::isConnectedTo(nid_t node, int pole)
    {
        // Rewrite:
        int closest, distance;
        getClosestPole(node, closest, distance);
        if (closest == -1) return false;
        if (closest == -2) return true;
        return pole == closest;
    }

    // Is the given node completely disconnected from all poles?
    bool GraphLayout::isDisconnected(nid_t node)
    {
        int closest, distance;
        getClosestPole(node, closest, distance);
        return closest == -1;
    }

    // Is the given node connected to two or more poles?
    bool GraphLayout::isConnectedToTwoPoles(nid_t node)
    {
        int closest, distance;
        getClosestPole(node, closest, distance);
        return closest == -2;
    }

    // Is the given node connected to ONLY one pole?
    bool GraphLayout::isConnectedToOneOnly(nid_t node)
    {
        int closest, distance;
        getClosestPole(node, closest, distance);
        return closest >= 0;
    }

    // Get color of a node depending on which pole it is connected to
    void GraphLayout::getNodeColor(nid_t n, double &r, double &g, double &b)
    { 
        if (pole_list_size == 0)
            r = 0.0, g = 0.0, b = 0.0;
        // If connected to multiple, set to gray
        else if (isConnectedToTwoPoles(n)) 
            r = 0.5, g = 0.5, b = 0.5;
        // If connected to a particular one, set to one from a pallette
        else if (isConnectedTo(n, 0))
            r = 1.0, g = 0.5, b = 0.0; // orange
        else if (isConnectedTo(n, 1))
            r = 0.0, g = 0.65, b = 0.93; // sky blue
        else if (isConnectedTo(n, 2))
            r = 0.85, g = 0.39, b = 0.6; // pink
        else if (isConnectedTo(n, 3))
            r = 0.5, g = 0.72, b = 0;  // lime
        else if (isConnectedTo(n, 4))
            r = 0.14, g = 0, b = 0.68;  // navy
        else if (isConnectedTo(n, 5))
            r = 0.53, g = 0.36, b = 0.0; // brown
        else if (isConnectedTo(n, 6))
            r = 0.2, g = 0.75, b = 0.58; // teal
        else if (isConnectedTo(n, 7))
            r = 0.65, g = 0, b = 0; // dark red
        else if (isConnectedTo(n, 8))
            r = 0, g = 0.41, b = 0.07;  // dark green
        else if (isConnectedTo(n, 9))
            r = 0.67, g = 0, b = 0.64;  // purple
        else if (isConnectedToOneOnly(n))
            r = 1.0, g = 0.0, b = 0.0; // red for all other nodes
        // If not connected to either, set to black
        else
            r = 0.0, g = 0.0, b = 0.0;
    }

    // Add nodes to the set that are connected to the given node through any path
    // No longer used
    /*void GraphLayout::addConnectedNodes(std::unordered_set<nid_t> &connected_nodes, nid_t node)
    {
        std::queue<nid_t> q;
        q.push(node);
        connected_nodes.insert(node);

        while (!q.empty())
        {
            nid_t n = q.front();
            q.pop();

            for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1)
            {
                for (nid_t n2 : graph.neighbors_with_geq_id(n1))
                {
                    // INEFFICIENT ITERATION; connects from higher to lower id
                    if (n1 == n && connected_nodes.find(n2) == connected_nodes.end() && (getEdgeDirection(n1, n2) < 0 
                        || (!use_distance_based_edge_direction && getEdgeDirection(n1, n2) == 0)))
                    {
                        connected_nodes.insert(n2);
                        q.push(n2);
                    }
                    if (n2 == n && connected_nodes.find(n1) == connected_nodes.end() && (getEdgeDirection(n1, n2) > 0
                        || (!use_distance_based_edge_direction && getEdgeDirection(n1, n2) == 0)))
                    {
                        connected_nodes.insert(n1);
                        q.push(n1);
                    }
                }
            }
        }
    }*/

    bool GraphLayout::sameColor(nid_t n1, nid_t n2)
    {
        double r1, g1, b1;
        double r2, g2, b2;
        getNodeColor(n1, r1, g1, b1);
        getNodeColor(n2, r2, g2, b2);
        return r1 == r2 && g1 == g2 && b1 == b2;
    }

    // Performs BFS to get a list of shortest distances from the given node to all other nodes
    void GraphLayout::getShortestDistances(nid_t node, std::vector<int> &distances)
    {
        std::queue<nid_t> q;
        std::unordered_set<nid_t> visited;
        q.push(node);
        visited.insert(node);
        distances.resize(graph.num_nodes());
        // Fill in distances with -1
        for (int i = 0; i < graph.num_nodes(); ++i)
            distances[i] = -1;

        distances[node] = 0;

        while (!q.empty())
        {
            nid_t n = q.front();
            q.pop();

            for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1)
            {
                for (nid_t n2 : graph.neighbors_with_geq_id(n1))
                {
                    // INEFFICIENT ITERATION
                    if (n1 == n && visited.find(n2) == visited.end() && (getInitialEdgeDirection(n1, n2) <= 0 
                        || use_distance_based_edge_direction))
                    {
                        if (max_influence_distance != -1 && distances[n1] >= max_influence_distance) continue;
                        visited.insert(n2);
                        q.push(n2);
                        distances[n2] = distances[n] + 1;
                    }
                    if (n2 == n && visited.find(n1) == visited.end() && (getInitialEdgeDirection(n1, n2) >= 0 
                        || use_distance_based_edge_direction))
                    {
                        if (max_influence_distance != -1 && distances[n2] >= max_influence_distance) continue;
                        visited.insert(n1);
                        q.push(n1);
                        distances[n1] = distances[n] + 1;
                    }
                }
            }
        }
    }

    // Contains a singleton instance of the shortest distances
    std::vector<int>* GraphLayout::getShortestDistancesList() {
        static std::vector<int>* distances = new std::vector<int>[pole_list_size];
        static bool initialized = false;
        if (!initialized) {
            for (int i = 0; i < pole_list_size; i++)
            {
                getShortestDistances(graph.node_map[pole_list[i]], distances[i]);
            }
            initialized = true;
        }
        return distances;
    }

    // Returns the closest pole to the given node and the distance to it
    void GraphLayout::getClosestPole(nid_t node, int &pole, int &distance)
    {
        int min_distance = INT_MAX;
        bool the_same_distance = false;
        for (int i = 0; i < pole_list_size; i++)
        {
            int d = getShortestDistancesList()[i][node];
            if (d < min_distance && d >= 0)
            {
                min_distance = d;
                pole = i;
                the_same_distance = false;
            } else if (d == min_distance)
            {
                the_same_distance = true;
            }
        }
        if (the_same_distance) { 
            distance = min_distance;
            pole = -2;
        } else if (min_distance != INT_MAX) {
            distance = min_distance;
        } else {
            distance = 0;
            pole = -1;
        }
    }

    int GraphLayout::getDistanceBasedEdgeDirection(nid_t n, nid_t t)
    {
        int n_dist, t_dist;
        int n_pole, t_pole;
        getClosestPole(n, n_pole, n_dist);
        getClosestPole(t, t_pole, t_dist);
        if (max_influence_distance != -1 && (n_dist > max_influence_distance || t_dist > max_influence_distance))
            return 0;
        if (n_pole == -2 && t_pole >= 0)
            return 1;
        if (t_pole == -2 && n_pole >= 0)
            return -1;
        if (n_pole != t_pole)
            return 0;
        else if (n_dist > t_dist)
            return 1;
        else if (n_dist == t_dist)
            return 0;
        else
            return -1;
    }


}
