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
namespace RPGraph
{
    GraphLayout::GraphLayout(UGraph &graph, float width, float height)
        : graph(graph), width(width), height(height)
    {
        coordinates = (Coordinate *)malloc(graph.num_nodes() * sizeof(Coordinate));
    }

    GraphLayout::~GraphLayout()
    {
        free(coordinates);
    }

    void GraphLayout::randomizePositions()
    {
        /*if (is_file_exists("communities.txt")&&is_file_exists("coordinate.txt"))
        {
         uint32_t node,comm;
            int xc,yc;
          uint32_t *communities2 = new uint32_t[graph.num_nodes()];
       int *xcoordinate = new int[graph.num_nodes()];
       int *ycoordinate = new int[graph.num_nodes()];

         std::ifstream inFile;
         inFile.open("communities.txt");

          for(int i=0;inFile >>node >>comm != NULL ;i++)
          {
              communities2[i]=comm;
          }
            inFile.close();
          inFile.open("coordinate.txt");
          for(int i=0;inFile >>comm >>xc >>yc != NULL ;i++)
          {
              xcoordinate[comm]=xc;
              ycoordinate[comm]=yc;
          }

           for (nid_t i = 0; i <  graph.num_nodes(); ++i)
          {
              setX(i,xcoordinate[communities2[i]]);
              setY(i,ycoordinate[communities2[i]]);
          }
         }
          else*/
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
        float maxX = std::numeric_limits<float>::min();
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
        float maxY = std::numeric_limits<float>::min();
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

        const float xRange = getXRange();
        const float yRange = getYRange();
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

            int threash = 100;
            if (r == double(177) / double(255))
                radian = 1;
            // if(radian<=threash)//sqrt(graph.degree(i)<10))
            /*char  c_string[]= "~/OpenSans-Bold.ttf";
            char c_string2[]= "7";*/
            //    if(graph.node_map_r[i]==748264||graph.node_map_r[i]==72606)
            layout_png.filledcircle_blend((getX(i) - minX) * xScale,
                                          (getY(i) - minY) * yScale,
                                          std::min(threash, radian), 0.8, r, g, b);
            /*layout_png.plot_text_utf8_blend(c_string, 5,

                                (getX(i) - minX)*xScale, (getY(i) - minY)*yScale, 0.0, c_string2,
                                1.0,

                                0.0,0.0 ,0.0 );*/

            for (nid_t n2 : graph.neighbors_with_geq_id(i))
            {
                // ... and edge.
                layout_png.line_blend((getX(i) - minX) * xScale, (getY(i) - minY) * yScale,
                                      (getX(n2) - minX) * xScale, (getY(n2) - minY) * yScale,
                                      0.005, r, g, b);
            }

            nid_t id = graph.node_map_r[i]; // id as found in edgelist
            out_color << id << " " << r << " " << g << " " << b << "\n";

            // ct++;
        }

        // Write it to disk.
        layout_png.write_png();
        out_color.close();
        //********************************************
        //********************************************
        //********************************************
        //********************************************
        //********************************************

        /* std::ofstream out_file("../../../files/coordinate.txt");
         std::ofstream out_file2("../../../files/coordinate2.txt");
     for (nid_t n = 0; n < graph.num_nodes(); ++n)
     {
         nid_t id = graph.node_map_r[n]; // id as found in edgelist
         out_file << id << " " << (getX(n) - minX)*xScale << " " <<  (getY(n) - minY)*yScale<<" " << graph.degree(n)<<"\n";
         out_file2 << id << " " << getX(n) << " " <<  getY(n)<<"\n";
     }

     out_file.close();
         out_file2.close();*/
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

}
