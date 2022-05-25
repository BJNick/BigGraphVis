/*
 ==============================================================================

 RPGraph.cpp
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

#include <iostream>
using namespace std; 
#include <string>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "RPGraph.hpp"

namespace RPGraph
{
    /* Definitions for UGraph */
   // UGraph::UGraph(std::string edgelist_path)
    UGraph::UGraph(nid_t *s, nid_t *t,nid_t n,nid_t m,nid_t *communities,nid_t *degree_S)
    {
        node_count = 0;
        edge_count = 0;
        cout<<"infos "<<n<<"```"<<m<<'\n';
/*
        std::fstream edgelist_file(edgelist_path, std::ifstream::in);

        std::string line;
        while(std::getline(edgelist_file, line))
        {
            // Skip any comments
            if(line[0] == '#') continue;

            // Read source and target from file
            nid_t s, t;
            std::istringstream(line) >> s >> t;

            if(s != t and !has_edge(s, t)) add_edge(s, t);
        }

        edgelist_file.close();*/
        //========================
        for(int i=0;i<n;i++)
        {
            if(s[i]<=m&&t[i]<=m)
                if(communities[s[i]]>0&&communities[s[i]]<m&&communities[t[i]]>0&&communities[t[i]]<m)
                    if(communities[s[i]] != communities[t[i]]) 
                        add_edge(communities[s[i]], communities[t[i]],degree_S[communities[s[i]]],degree_S[communities[t[i]]]); 
           
        }
       
        
        
        
        
        
        //*********************************************************************************************************************************************************
    }

    bool UGraph::has_node(nid_t nid)
    {
        return node_map.count(nid) > 0;
    }

    bool UGraph::has_edge(nid_t s, nid_t t)
    {
        if(!has_node(s) or !has_node(t)) return false;

        nid_t s_mapped = node_map[s];
        nid_t t_mapped = node_map[t];

        if(adjacency_list.count(std::min(s_mapped, t_mapped)) == 0) return false;

        std::vector<nid_t> neighbors = adjacency_list[std::min(s_mapped, t_mapped)];
        if(std::find(neighbors.begin(), neighbors.end(), std::max(s_mapped, t_mapped)) == neighbors.end())
            return false;
        else
            return true;
    }

    void UGraph::add_node(nid_t nid,nid_t degree_s)
    {
        if(!has_node(nid))
        {
            node_map[nid] = node_count;
            node_map_r[node_count] = nid;
            //degrees[node_map[nid]] = degree_s;
            node_count++;
            //cout<<node_count<<'\n';
        }
    }

    void UGraph::add_edge(nid_t s, nid_t t,nid_t degree_s,nid_t degree_t)
    {
        
        if(!has_node(s)) add_node(s,degree_s);
        if(!has_node(t)) add_node(t,degree_t);
        nid_t s_mapped = node_map[s];
        nid_t t_mapped = node_map[t];
      /*  if(degree_s>0||degree_t>0)
           cout<<degree_s<<" "<<degree_t<<'\n';*/
        // Insert edge into adjacency_list
        if (!has_edge(s, t)) {
            adjacency_list[std::min(s_mapped, t_mapped)].push_back(std::max(s_mapped, t_mapped));
            degrees[s_mapped] += degree_s+1;
            degrees[t_mapped] +=degree_t+ 1;
            edge_count++;
        }

        // Check if connection already exists in the map
        if (std::count(in_adj_list[t_mapped].begin(), in_adj_list[t_mapped].end(), s_mapped) != 0) {
            return;
        }

        // Account for directed edges
        in_adj_list[t_mapped].push_back(s_mapped);
        in_degrees[t_mapped] += 1;
        initial_edge_direction[s_mapped][t_mapped] = true;
        initial_edge_direction[t_mapped][s_mapped] = false;
        is_edge_directed[s_mapped][t_mapped] = !is_edge_directed[s_mapped][t_mapped];
        is_edge_directed[t_mapped][s_mapped] = !is_edge_directed[t_mapped][s_mapped];
    }

    nid_t UGraph::num_nodes()
    {
        return node_count;
    }

    nid_t UGraph::num_edges()
    {
        return edge_count;
    }

    nid_t UGraph::degree(nid_t nid)
    {
        return degrees[nid];
    }
  /*  nid_t UGraph::weights(nid_t nid)
    {
        return weight_S[nid];
    }*///**

    nid_t UGraph::in_degree(nid_t nid)
    {
        return in_degrees[nid];
    }

    nid_t UGraph::out_degree(nid_t nid)
    {
        return degree(nid);
    }
    std::vector<nid_t> UGraph::neighbors_with_geq_id(nid_t nid)
    {
        return adjacency_list[nid];
    }

    // Get OG indices of top N nodes with highest degree without sorting 
    int* UGraph::get_top_nodes(int N)
    {
        std::unordered_map<nid_t, nid_t> degrees = this->in_degrees;
	    int* top_nodes = new int[N];
        int* top_degrees = new int[N];
        int i = 0;
        for(auto it = degrees.begin(); it != degrees.end(); ++it)
        {
            if(i < N)
            {
                top_nodes[i] = it->first;
                top_degrees[i] = it->second;
                i++;
            }
            else
            {
                if(it->second > top_degrees[N-1])
                {
                    top_nodes[N-1] = it->first;
                    top_degrees[N-1] = it->second;
                    int indices[N];
                    for(int j = 0; j < N; j++)
                        indices[j] = j;
                    std::sort(indices, indices + N, [top_degrees](int a, int b) { return top_degrees[a] > top_degrees[b]; });
                    int new_top_nodes[N];
                    int new_top_degrees[N];
                    for(int j = 0; j < N; j++)
                    {
                        new_top_nodes[j] = top_nodes[indices[j]];
                        new_top_degrees[j] = top_degrees[indices[j]];
                    }
                    for(int j = 0; j < N; j++)
                    {
                        top_nodes[j] = new_top_nodes[j];
                        top_degrees[j] = new_top_degrees[j];
                    }
                }
            }
        }
        // Convert all of them using node_map_r
        for(int i = 0; i < N; i++)
            top_nodes[i] = node_map_r[top_nodes[i]];
        /*std::cout << "Top nodes: {";
        for (int i = 0; i < 10; i++)
            cout << top_nodes[i] << "," << top_degrees[i] << (i == 10 - 1 ? "" : "|");
        std::cout << "}\n";*/
        return top_nodes;
    }

    /* Definitions for CSRUGraph */

// CSRUGraph represents an undirected graph using a
// compressed sparse row (CSR) datastructure.
    CSRUGraph::CSRUGraph(nid_t num_nodes, nid_t num_edges)
    {
        // `edges' is a concatenation of all edgelists
        // `offsets' contains offset (in `edges`) for each nodes' edgelist.
        // `nid_to_offset` maps nid to index to be used in `offset'

        // e.g. the edgelist of node with id `nid' starts at
        // edges[offsets[nid_to_offset[nid]]] and ends at edges[offsets[nid_to_offset[nid]] + 1]
        // (left bound inclusive right bound exclusive)

        edge_count = num_edges; // num_edges counts each bi-directional edge once.
        node_count = num_nodes;
        edges =   (nid_t *) malloc(sizeof(nid_t) * 2 * edge_count);
        offsets = (nid_t *) malloc(sizeof(nid_t) * node_count);
        offset_to_nid = (nid_t *) malloc(sizeof(nid_t) * node_count);

        // Create a map from original ids to ids used throughout CSRUGraph
        nid_to_offset.reserve(node_count);

        first_free_id = 0;
        edges_seen = 0;
    }

    CSRUGraph::~CSRUGraph()
    {
        free(edges);
        free(offsets);
        free(offset_to_nid);
    }

    void CSRUGraph::insert_node(nid_t node_id, std::vector<nid_t> nbr_ids)
    {
        nid_t source_id_old = node_id;
        nid_t source_id_new = first_free_id;
        nid_to_offset[source_id_old] = first_free_id;
        offset_to_nid[first_free_id] = source_id_old;
        first_free_id++;

        offsets[source_id_new] = edges_seen;
        for (auto nbr_id : nbr_ids)
        {
            nid_t dest_id_old = nbr_id;
            edges[edges_seen] = dest_id_old;
            edges_seen++;
        }
    }

    void CSRUGraph::fix_edge_ids()
    {
        for (eid_t ei = 0; ei < 2*edge_count; ei++)
        {
            edges[ei] = nid_to_offset[edges[ei]];
        }
    }

    nid_t CSRUGraph::degree(nid_t nid)
    {
        // If nid is last element of `offsets'... we prevent out of bounds.
        nid_t r_bound;
        if (nid < node_count - 1) r_bound = offsets[nid+1];
        else r_bound = edge_count * 2;
        nid_t l_bound = offsets[nid];
        return (r_bound - l_bound);
    }

    nid_t CSRUGraph::out_degree(nid_t nid)
    {
        return degree(nid);
    }

    nid_t CSRUGraph::in_degree(nid_t nid)
    {
        return degree(nid);
    }

    nid_t CSRUGraph::nbr_id_for_node(nid_t nid, nid_t edge_no)
    {
        return edges[offsets[nid] + edge_no];
    }
    nid_t CSRUGraph::num_nodes()
    {
        return node_count;
    }

    nid_t CSRUGraph::num_edges()
    {
        return edge_count;
    }
};
