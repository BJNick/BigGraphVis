
using namespace std;

#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <chrono> 
#include <curand.h>
#include <curand_kernel.h>
#include <ctime>
#include <functional>
#include <cctype>
#include <algorithm> 
#include <sstream>
#include <iomanip>
#include <map>
#include <ctype.h>
#include <unordered_map>

using namespace std::chrono;

#include "RPCommon.hpp"
#include "RPGraph.hpp"
#include "RPGraphLayout.hpp"
#include "RPCPUForceAtlas2.hpp"

#ifdef __NVCC__
#include <cuda_runtime_api.h>
#include "RPGPUForceAtlas2.hpp"
#endif

/* 
* The following are global CUDA kernels used by the modified SCoDA algorithm
* for running operations in parallel on the GPU.
*/

__global__
void initial(int num_of_edges, int num_of_nodes, uint32_t* communities, uint32_t* degree, uint32_t* degree_cmt, uint32_t* degree_S, int heunumber, int* sketch)
{
	// Get the index of a NODE that is being processed by the thread
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index < num_of_nodes)
	{
		// Set default values for each node in the graph
		communities[index] = index;
		degree[index] = 0;
		degree_S[index] = 1;
		degree_cmt[index] = 0;
		sketch[index] = 1;
	}

}

__global__
void mainfor(uint32_t num_of_edges, uint32_t num_of_nodes, int degree_threshold, uint32_t* degree_cmt, uint32_t* src, uint32_t* dst, uint32_t* communities)
{
	// Get the index of an EDGE that is being processed by the thread
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (src[index] > 0 && src[index] <= num_of_nodes && dst[index] > 0 && dst[index] <= num_of_nodes)
	{
		if (degree_cmt[src[index]] <= degree_threshold && degree_cmt[dst[index]] <= degree_threshold)
		{
			// If degrees are less than the threshold, join the two communities
			if (degree_cmt[src[index]] > degree_cmt[dst[index]])
			{
				communities[dst[index]] = communities[src[index]];
				atomicAdd(&degree_cmt[src[index]], 1);
			}
			else
			{
				communities[src[index]] = communities[dst[index]];
				atomicAdd(&degree_cmt[dst[index]], 1);
			}
		}
	}
}

__global__
void communities_hashing(int num_of_edges, int num_of_nodes, uint32_t* communities, uint32_t* degree, int huenumber, uint32_t* dst, uint32_t* src, uint32_t degree_thresholdS, int* sketch, int* Degree_done)
{
	// Get the index of an EDGE that is being processed by the thread
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	// Creates hashes for finding the size of the communities
	if (index < num_of_edges)
	{
		//hash 1 FVN-- hash2 ELF --hash3 BKDR--hash4 DJB
		const unsigned int fnv_prime = 0x811C9DC5;
		unsigned int hash1 = 0, hash2 = 0, hash3 = 0, hash4 = 5381;
		uint32_t community = 0;
		unsigned int x = 0;
		int deg = 0;
		unsigned int seed = 131;

		if (Degree_done[src[index]] == 0 || Degree_done[dst[index]] == 0)// check if the node is visited before
		{
			if (Degree_done[src[index]] == 0)
			{
				community = communities[src[index]] + 1;
				deg = int(degree[src[index]] / degree_thresholdS);
				Degree_done[src[index]] = 1;
			}
			else
			{
				community = communities[dst[index]] + 2;
				deg = int(degree[dst[index]] / degree_thresholdS);
				Degree_done[dst[index]] = 1;
			}
			while (community > 10)
			{
				hash1 *= fnv_prime;
				hash1 ^= (community % 10);

				hash2 = (hash2 << 4) + (community % 10);
				if ((x = hash2 & 0xF0000000) != 0)
				{
					hash2 ^= (x >> 24);
				}
				hash2 &= ~x;

				hash3 = (hash3 * seed) + (community);

				hash4 = ((hash4 << 5) + hash4) + (community);

				community = uint32_t(community / 10);
			}
			sketch[int((hash1 % 1000) % int(huenumber / 4))] += deg + 1;
			sketch[int((hash2 % 1000) % (int(huenumber / 4) + int(huenumber / 4)))] += deg + 1;
			sketch[int((hash3 % 1000) % (int(huenumber / 4) + int(huenumber / 2)))] += deg + 1;
			sketch[int((hash4 % 1000) % (int(huenumber / 4) + int((3 * huenumber) / 4)))] += deg + 1;
		}
	}
}

__global__
void find_degree_S(int num_of_edges, int num_of_nodes, uint32_t* communities, uint32_t* degree, uint32_t* degree_S, int* sketch, int huenumber, uint32_t* dst, uint32_t* src)
{
	// Get the index of a COMMUNITY that is being processed by the thread
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	// Computes the size (weight) of each community for the supergraph
	if (index < num_of_nodes)
	{
		//hash 1 FVN-- hash2 ELF --hash3 BKDR--hash4 DJB
		const unsigned int fnv_prime = 0x811C9DC5;
		unsigned int hash1 = 0, hash2 = 0, hash3 = 0, hash4 = 5381;
		uint32_t community = communities[index] + 1;
		unsigned int x = 0;
		unsigned int seed = 131;

		while (community > 10)
		{
			hash1 *= fnv_prime;
			hash1 ^= (community % 10);

			hash2 = (hash2 << 4) + (community % 10);
			if ((x = hash2 & 0xF0000000) != 0)
			{
				hash2 ^= (x >> 24);
			}
			hash2 &= ~x;

			hash3 = (hash3 * seed) + (community);

			hash4 = ((hash4 << 5) + hash4) + (community);

			community = uint32_t(community / 10);
		}
		degree_S[communities[index]] = fminf
		(
			fminf
			(
				sketch[int(hash1 % int(huenumber / 4))], sketch[int(hash2 % (int(huenumber / 4) + int(huenumber / 4)))]
			)
			, fminf
			(
				sketch[int(hash3 % (int(huenumber / 4) + int(huenumber / 2)))], sketch[int(hash4 % (int(huenumber / 4) + int((3 * huenumber) / 4)))]
			)
		);

		if (degree_S[communities[index]] == 0)
			degree_S[communities[index]] = 1;
	}

}

//============================================================

const int num_of_parameters = 60; // arbitrary number

std::string parameter_keys[num_of_parameters] = {
	// ForceAtlas2 parameters:
	"program_call", "cuda_requested", "max_iterations", "num_screenshots", "strong_gravity", "scale", "gravity", "approximate",
	"in_path", "out_path", "out_format", "image_w", "image_h", "degree_threshold", "rounds", "huenumber", "use_linlog",
	// Configuration file parameters:
	"config_folder", "config_chain", "chain_output_name", "chain_separator", "include_timestamp",
	// Extra parameters:
	"community_detection", "attraction_exponent", "attraction", "random_seed", "pin_2_roots", "repulsion_d_squared",
	"stop_on_divergence", "divergence_factor", "divergence_threshold", "avoid_cpu_allocation",
	// Pole parameters:
	"use_distance_based_edge_direction", "magnetic_pole_separation", "draw_common_edges",
	"max_influence_distance", "pin_poles", "extra_pole_attraction", "use_pole_segmentation", "pole_list",
	"pole_size_factor", "top_N_nodes", "pole_gravity_factor",
	// Magnetic field parameters:
	"use_magnetic_field", "field_type", "bi_directional", "field_strength", "magnetic_constant", "magnetic_alpha", "magnetic_beta",
	"legacy_segmentation", "simple_center_of_mass",
	// Cosmetic parameters:
	"node_alpha", "edge_alpha", "square_coordinates", "draw_arrows", "min_arrow_length", "colored_fraction",
	"hide_iteration_output",
}; 

// A helpful method for naming the output files
std::string fill_zeros(int number, int digits)
{
	std::stringstream ss;
	ss << std::setw(digits) << std::setfill('0') << number;
	return ss.str();
}

// Store command line arguments to a map
void store_argv(int argc, const char** argv, map<string, string>& map)
{
	if (argc == 16) {
		for (int i = 0; i < argc; i++)
		{
			string key = parameter_keys[i];
			string value = argv[i];
			map[key] = value;
		}
	}
}

// Read the input file and store the fields to a map
void read_args_from_file(string file_path, map<string, string>& map)
{
	// Check if this is a file name not complete path
	bool name_only = file_path.find("/") == string::npos;
	if (name_only) 
	{
		// Check if file extension is ommitted
		if (file_path.find(".") == string::npos)
			file_path += ".txt";
		file_path = map["config_folder"] + file_path;
	}

	if (!is_file_exists(file_path))
	{
		cout << "File does not exist " << file_path << "\n";
		exit(EXIT_FAILURE);
	}

	std::ifstream file(file_path);
	string line;
	while (std::getline(file, line))
	{
		// Skip lines that start with #, they are comments
		if (line[0] == '#')
			continue;
		std::istringstream ss(line);
		string key, value;
		if (std::getline(ss, key, '='))
		{
			// If the key is whitespace, skip the line
			key.erase(remove_if(key.begin(), key.end(), ::isspace), key.end());
			if (key.empty())
				continue;
			if (key == "BREAK" || key == "break")
				break;
			// Make sure the key is in the list of keys
			bool found = (std::find(parameter_keys, parameter_keys + num_of_parameters, key) != parameter_keys + num_of_parameters);
			if (!found)
			{
				cout << "Config key \"" << key << "\" is not a valid key\n";
				cout << "Valid keys are: ";
				for (int i = 0; i < num_of_parameters; i++)
				{
					cout << parameter_keys[i] << ", ";
				}
				cout << "\n";
				exit(EXIT_FAILURE);
			}
			if (std::getline(ss, value))
			{
				// Erase unnecessary spaces
				value.erase(value.find_last_not_of("\n\r\t")+1);
				map[key] = value;
			}
		}
	}

	// Store all config names in a chain
	std::string just_name = file_path.substr(file_path.find_last_of("/") + 1);
	just_name = just_name.substr(0, just_name.find_last_of("."));
	map["config_chain"] += just_name + map["chain_separator"];
}

// Set the default values for the parameters
void set_default_args(map<string, string>& map)
{
	// Set default values for the parameters, following
	// ./graph_viewer gpu 500 1 sg 80 1 approximate ~/net/web-BerkStan.txt  ~/output png 1024 1024 11 5 6500 
	map["program_call"] = "";
	map["cuda_requested"] = "gpu";
	map["max_iterations"] = "500";
	map["num_screenshots"] = "1";
	map["strong_gravity"] = "sg";
	map["scale"] = "80";
	map["gravity"] = "1";
	map["approximate"] = "approximate";
	map["in_path"] = "../../../net/web-BerkStan.txt";
	map["out_path"] = "../../../output/";
	map["out_format"] = "png";
	map["image_w"] = "1024";
	map["image_h"] = "1024";
	map["degree_threshold"] = "11";
	map["rounds"] = "5";
	map["huenumber"] = "6500";
	map["use_linlog"] = "false";
	// Configuration file parameters
	map["config_folder"] = "../../../config/";
	map["config_chain"] = "";
	map["chain_output_name"] = "false";
	map["chain_separator"] = " ";
	map["include_timestamp"] = "true";
	// Extra parameters
	map["community_detection"] = "SCoDA";
	map["attraction_exponent"] = "1";
	map["attraction"] = "1";
	map["random_seed"] = "1234";
	map["pin_2_roots"] = "false";
	map["repulsion_d_squared"] = "false";
	map["stop_on_divergence"] = "false";
	map["divergence_factor"] = "1.75";
	map["divergence_threshold"] = "1e+8";
	map["avoid_cpu_allocation"] = "false";
	// Pole parameters
	map["pole_list"] = "";
	map["pin_poles"] = "false";
	map["use_distance_based_edge_direction"] = "false";
	map["max_influence_distance"] = "-1";	
	map["extra_pole_attraction"] = "1";
	map["use_pole_segmentation"] = "false";
	map["magnetic_pole_separation"] = "10000";
	map["draw_common_edges"] = "true";
	map["pole_size_factor"] = "3";
	map["top_N_nodes"] = "0";
	map["pole_gravity_factor"] = "0";
	map["legacy_segmentation"] = "false";
	// Magnetic force parameters
	map["use_magnetic_field"] = "false";
	map["field_type"] = "linear";
	map["bi_directional"] = "false";
	map["field_strength"] = "16";
	map["magnetic_constant"] = "1";
	map["magnetic_alpha"] = "1";
	map["magnetic_beta"] = "1";
	map["simple_center_of_mass"] = "false";
	// Cosmetic parameters
	map["node_alpha"] = "0.8";
	map["edge_alpha"] = "0.005";
	map["square_coordinates"] = "false";
	map["draw_arrows"] = "false";
	map["min_arrow_length"] = "50";
	map["colored_fraction"] = "1";
	map["hide_iteration_output"] = "false";
}

// Split a string into an array of integers
int* split_to_int(string& s, char delim, int* arr, int max_size)
{
	int i = 0;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim) && i < max_size)
	{
		arr[i] = atoi(item.c_str());
		i++;
	}
	while (i < max_size)
	{
		arr[i] = -1;
		i++;
	}
	return arr;
}

/*
	uint32_t *communities, *src, *dst, num_of_edges, num_of_nodes;
	uint32_t *degree, *degree_cmt, *degree_S;
	int src_id, dst_id, *sketch, *h1, *h2, *h3, *h4, *Degree_done, *weight_S;
*/

void read_data_from_file(string in_path, uint32_t* src, uint32_t* dst, uint32_t* degree, unordered_map<long, string>& node_labels) {
	// Read the input file again, now storing the edge connections
	// in the CUDA allocated variables
	// Read the file line by line instead of tokens
	ifstream inFile;
	if (!is_file_exists(in_path)) {
		cout << "error: File not found: " << in_path << "\n";
		exit(EXIT_FAILURE);
	}
	inFile.open(in_path);
	std::string line;
	int edge_id = 0;
	bool edge_mode = true;
	bool vertex_mode = false;
	while (getline(inFile, line)) {
		// if line starts with a *, check if it says *Edges or *Arcs
		if (line[0] == '*') {
			if (line.find("Edges") != std::string::npos || line.find("Arcs") != std::string::npos) {
				edge_mode = true;
				vertex_mode = false;
				continue;
			} else if (line.find("Vertices") != std::string::npos) {
				edge_mode = false;
				vertex_mode = true;
				continue;
			} else {
				cout << "Error reading input file" << '\n';
				exit(EXIT_FAILURE);
			}
		}
		if (edge_mode) {
			std::string src_id_s, dst_id_s;
			long src_id, dst_id;
			istringstream iss(line);
			iss >> src_id_s >> dst_id_s;
			src_id = stol(src_id_s);
			dst_id = stol(dst_id_s);
			src[edge_id] = src_id;
			dst[edge_id] = dst_id;
			// Compute the degree of each node
			degree[src_id]++;
			if (dst_id != src_id)
				degree[dst_id]++;
			edge_id++;
		} else if (vertex_mode) {
			// The string is of format 
			//  1 "S. Kambhampati"    0.0000    0.0000    0.5000 metadata
			long node_id = 0;
			std::string node_label;
			// Split the string using quotes as delimiter
			std::stringstream ss(line);
			std::string token;
			// if line contains double quotes, split using double quotes
			if (line.find("\"") != std::string::npos) {
				int i = 0;
				while (getline(ss, token, '\"')) {
					if (i == 0) {
						node_id = stol(token);
					} else if (i == 1) {
						node_label = token;
						break;
					}
					i++;
				}
			} else {
				// if line does not contain double quotes, split using whitespace
				istringstream iss(line);
				std::string node_id_s;
				iss >> node_id_s >> node_label;
				node_id = stol(node_id_s);
			}
			node_labels[node_id] = node_label;
		}
	}
	inFile.close();
}

void read_node_edge_counts(string in_path, uint32_t& num_of_nodes, uint32_t& num_of_edges) {
	// Read the file line by line instead of tokens
	ifstream inFile;
	if (!is_file_exists(in_path)) {
		cout << "error: File not found: " << in_path << "\n";
		exit(EXIT_FAILURE);
	}
	inFile.open(in_path);
	std::string line;
	int edge_id = 0;
	bool edge_mode = true;
	bool vertex_mode = false;
	num_of_nodes = 0;
	num_of_edges = 0;
	uint16_t max_node_id = 0;
	while (getline(inFile, line)) {
		// if line starts with a *, check if it says *Edges or *Arcs
		if (line[0] == '*') {
			if (line.find("Edges") != std::string::npos || line.find("Arcs") != std::string::npos) {
				edge_mode = true;
				vertex_mode = false;
				continue;
			} else if (line.find("Vertices") != std::string::npos) {
				edge_mode = false;
				vertex_mode = true;
				continue;
			} else {
				cout << "Error reading input file" << '\n';
				exit(EXIT_FAILURE);
			}
		}
		if (edge_mode) {
			std::string src_id_s, dst_id_s;
			long src_id, dst_id;
			istringstream iss(line);
			iss >> src_id_s >> dst_id_s;
			src_id = stol(src_id_s);
			dst_id = stol(dst_id_s);
			if (src_id > max_node_id)
				max_node_id = src_id;
			if (dst_id > max_node_id)
				max_node_id = dst_id;
			num_of_edges++;
			edge_id++;
		} else if (vertex_mode) {
			num_of_nodes++;
		}
	}
	num_of_nodes = max_node_id;
	num_of_edges = edge_id;
	inFile.close();
}

//============================================================

int main(int argc, const char** argv)
{
	// --- Parse all the arguments ---

	// Check command-line usage
	if (argc < 2 || (argc > 8 && argc < 16))
	{
		fprintf(stderr, "Usage: graph_viewer gpu|cpu max_iterations num_snaps sg|wg scale gravity exact|approximate edgelist_path out_path png|csv|bin image_w image_h degree_threshold rounds heuristic\n");
		fprintf(stderr, "Or:    graph_viewer -c config_file ...\n");
		exit(EXIT_FAILURE);
	}

	// Either use the command-line arguments or the config file
	map<string, string> arg_map;
	set_default_args(arg_map);

	if (argc < 16) {
		// Use listed config files
		for (int i = 1; i < argc; i++)
		{
			string arg = argv[i];
			if (arg == "-c")
				continue;
			read_args_from_file(arg, arg_map);
		}
	} else {
		// Use command-line arguments
		store_argv(argc, argv, arg_map);
	}

	// --- Parsing finished ---

	// Measure the execution time
	auto start = high_resolution_clock::now();
	auto end = high_resolution_clock::now();
	auto end_cmt = high_resolution_clock::now();

	cout << "Initialization" << "\n";

	// Use a consistent random seed for reproducibility 
	// (note: this may have no effect on CUDA's random generation)
	srandom(stoi(arg_map["random_seed"]));
	
	// Allocate variables
	curandState_t* states; // CUDA random number generation states
	uint32_t *communities, *src, *dst, num_of_edges, num_of_nodes;
	uint32_t *degree, *degree_cmt, *degree_S;
	int src_id, dst_id, *sketch, *h1, *h2, *h3, *h4, *Degree_done, *weight_S;

	std::unordered_map<long, string> node_labels;

	// SCoDA Community Detection parameters
	uint32_t degree_threshold = std::stoul(arg_map["degree_threshold"]);
	uint32_t degree_thresholdS = degree_threshold;

	int rounds = std::stoi(arg_map["rounds"]);  // Number of  SCoDA rounds
	int huenumber = std::stoul(arg_map["huenumber"]);   // Heuristic number
	
	// Filepaths for I/O
	std::string in_path = arg_map["in_path"];
	std::string out_path = arg_map["out_path"];

	read_node_edge_counts(in_path, num_of_nodes, num_of_edges);

	// TODO: Better solution for iterating nodes
	num_of_nodes += 1; // Count the last node id as well

	cout << "Number of nodes: " << num_of_nodes << '\n';
	cout << "Number of edges: " << num_of_edges << '\n';

	const bool cuda_requested = std::string(arg_map["cuda_requested"]) == "gpu" or std::string(arg_map["cuda_requested"]) == "cuda";

	const bool avoid_cpu_allocation = std::string(arg_map["avoid_cpu_allocation"]) == "true";

	if (avoid_cpu_allocation || cuda_requested) {

		if (avoid_cpu_allocation)
			cout << "NOTICE: CUDA is still used for memory allocation due to errors with the CPU version.\n";

		// Sets the CUDA computing device
		int device_count;
		cudaGetDeviceCount(&device_count);

		if (device_count == 0) 
		{
			fprintf(stderr, "No CUDA-capable devices found.\n");
			exit(EXIT_FAILURE);
		} else {
			cout << device_count << " CUDA-capable devices found" << '\n';
		}

		// Just use the first one found
		cudaSetDevice(0);

		cout << "Allocating CUDA memory" << '\n';

		// Allocate GPU memory
		cudaMallocManaged(&src, num_of_edges * sizeof(uint32_t));
		cudaMallocManaged(&dst, num_of_edges * sizeof(uint32_t));
		cudaMallocManaged(&Degree_done, num_of_edges * sizeof(int));
		cudaMallocManaged(&degree, num_of_nodes * sizeof(uint32_t));
		cudaMallocManaged(&degree_cmt, num_of_nodes * sizeof(uint32_t));
		cudaMallocManaged(&degree_S, num_of_nodes * sizeof(uint32_t));
		cudaMallocManaged(&weight_S, num_of_edges * sizeof(int));

		cudaMallocManaged(&communities, num_of_nodes * sizeof(uint32_t));
		cudaMallocManaged(&h1, num_of_nodes * sizeof(int));
		cudaMallocManaged(&h2, num_of_nodes * sizeof(int));
		cudaMallocManaged(&h3, num_of_edges * sizeof(int));
		cudaMallocManaged(&h4, num_of_edges * sizeof(int));
		cudaMallocManaged(&sketch, huenumber * sizeof(int));

		cudaMallocManaged((void**)&states, num_of_edges * sizeof(curandState_t));
		
		// Compute the number of CUDA blocks based on total NODES
		int blockSize = 256;
		int remain = num_of_nodes % blockSize;
		int numBlocks = num_of_nodes / blockSize + (remain > 0 ? 1 : 0);

		// Run kernel on the GPU that sets the initial values for the thread variables
		initial<<<numBlocks, blockSize>>>(num_of_edges, num_of_nodes, communities, degree, degree_cmt, degree_S, huenumber, sketch);

		cudaDeviceSynchronize();

		cudaError_t error = cudaGetLastError();
		if (error != cudaSuccess)
		{
			printf("CUDA error in initialization: %s\n", cudaGetErrorString(error));
			exit(-1);
		}

		cout << "CUDA initialization complete" << '\n';

	} else {

		// THIS ALLOCATION DOES NOT WORK AS INTENDED

		// Allocate CPU memory for the arrays
		src = (uint32_t*)malloc(num_of_edges * sizeof(uint32_t));
		dst = (uint32_t*)malloc(num_of_edges * sizeof(uint32_t));
		Degree_done = (int*)malloc(num_of_edges * sizeof(int));
		degree = (uint32_t*)malloc(num_of_nodes * sizeof(uint32_t));
		degree_cmt = (uint32_t*)malloc(num_of_nodes * sizeof(uint32_t));
		degree_S = (uint32_t*)malloc(num_of_nodes * sizeof(uint32_t));
		weight_S = (int*)malloc(num_of_edges * sizeof(int));

		// For SCoDA detection
		communities = (uint32_t*)malloc(num_of_nodes * sizeof(uint32_t));
		h1 = (int*)malloc(num_of_nodes * sizeof(int));
		h2 = (int*)malloc(num_of_nodes * sizeof(int));
		h3 = (int*)malloc(num_of_edges * sizeof(int));
		h4 = (int*)malloc(num_of_edges * sizeof(int));
		sketch = (int*)malloc(huenumber * sizeof(int));
		states = (curandState_t*)malloc(num_of_edges * sizeof(curandState_t));

		// Set default values for each node in the graph
		for (int i = 0; i < num_of_nodes; i++) {
			communities[i] = i;
			degree[i] = 0;
			degree_S[i] = 1;
			degree_cmt[i] = 0;
			sketch[i] = 1;
		}

	}
	
	std::cout << "Started file reading" << '\n';

	read_data_from_file(in_path, src, dst, degree, node_labels);
	
	cout << "Graph Data Loaded" << "\n";
	
	// ######################################## SCoDA 

	/*
	* The following is a modified version of the SCoDA algorithm for 
	* community detection, by Hollocou et al.
	*/

	if (arg_map["community_detection"] == "SCoDA") {
		
		// Compute the number of CUDA blocks based on total NODES
		int blockSize = 256;
		int remain = num_of_nodes % blockSize;
		int numBlocks = num_of_nodes / blockSize + (remain > 0 ? 1 : 0);

		cout << "--- SCoDA Community Detection ---" << "\n";
		
		// Recompute the number of CUDA blocks based on total EDGES
		remain = num_of_edges % blockSize;
		numBlocks = num_of_edges / blockSize + (remain > 0 ? 1 : 0);

		if (rounds == 0)
			rounds = num_of_edges * 0.01;
		
		for (int i = 0; i < rounds; i++)
		{
			cout << "Community Detection Round " << i + 1 << "\n";

			degree_threshold = degree_threshold * degree_threshold;

			// Run kernel on the GPU that joins communities together according to SCoDA procedure
			mainfor<<<numBlocks, blockSize>>>(num_of_edges, num_of_nodes, degree_threshold, degree_cmt, src, dst, communities); 
			
			cudaDeviceSynchronize(); 
			cudaError_t error = cudaGetLastError(); 
			if (error != cudaSuccess) { 
				printf("CUDA error in SCoDA main loop: %s\n", cudaGetErrorString(error)); 
				exit(-1); 
			}
		}

		cout << "Hashing . . . " << "\n";

		communities_hashing<<<numBlocks, blockSize>>>(num_of_edges, num_of_nodes, communities, degree, huenumber, dst, src, degree_thresholdS, sketch, Degree_done); 
		
		cudaDeviceSynchronize(); 
		cudaError_t error = cudaGetLastError(); 
		if (error != cudaSuccess) { 
			printf("CUDA error in community hashing: %s\n", cudaGetErrorString(error)); 
			exit(-1); 
		}
		
		cout << "Counting . . . " << "\n";

		// Recompute the number of CUDA blocks based on COMMUNITIES
		remain = (num_of_edges / huenumber) % blockSize;
		numBlocks = (num_of_edges / huenumber + blockSize - remain) / blockSize;

		// Finds the degrees of nodes of the supergraph
		find_degree_S<<<numBlocks, blockSize>>>(num_of_edges, num_of_nodes, communities, degree, degree_S, sketch, huenumber, dst, src); 
		
		cudaDeviceSynchronize(); 
		error = cudaGetLastError(); 
		if (error != cudaSuccess) { 
			printf("CUDA error in finding degree_S: %s\n", cudaGetErrorString(error)); 
			exit(-1); 
		}
		
		// Write nodes categorized by community to file (for debugging)

		/*std::ofstream out_file("../../../files/communities.txt");
		for (uint32_t n = 0; n < num_of_nodes; ++n)
		{
			out_file << n << " " << communities[n] << "\n";
		}
		out_file.close();*/

		// Measure passed execution time
		end_cmt = high_resolution_clock::now();
		long cmt_time = chrono::duration_cast<chrono::milliseconds>(end_cmt - start).count();
		cout << "SCoDA algorithm finished in " << cmt_time << "ms \n";

	}
	
	// ######################################## ForceAtlas2

	/* 
	* The following is an applied version of the ForceAtlas2 algorithm
	*/

	// ForceAtlas2 parameters
	const int max_iterations = std::stoi(arg_map["max_iterations"]);
	const int num_screenshots = std::stoi(arg_map["num_screenshots"]);
	const bool strong_gravity = std::string(arg_map["strong_gravity"]) == "sg" or std::string(arg_map["strong_gravity"]) == "strong";
	const float scale = std::stof(arg_map["scale"]);
	const float gravity = std::stof(arg_map["gravity"]);
	const bool approximate = std::string(arg_map["approximate"]) == "approximate" or std::string(arg_map["approximate"]) == "yes";

	// Output format and parameters
	std::string out_format = arg_map["out_format"];
	int image_w = std::stoi(arg_map["image_w"]);
	int image_h = std::stoi(arg_map["image_h"]);

	cout << "--- Force Atlas 2 ---" << "\n";

	if (cuda_requested and not approximate)
	{
		fprintf(stderr, "error: The CUDA implementation (currently) requires Barnes-Hut approximation.\n");
		exit(EXIT_FAILURE);
	}

	// Check in_path and out_path
	if (!is_file_exists(in_path))
	{
		fprintf(stderr, "error: No edgelist at %s\n", in_path.c_str());
		exit(EXIT_FAILURE);
	}

	if (!is_file_exists(out_path))
	{
		fprintf(stderr, "error: No output folder at %s\n", out_path.c_str());
		exit(EXIT_FAILURE);
	}

	// If not compiled with cuda support, check if cuda is requested.
#ifndef __NVCC__
	if (cuda_requested)
	{
		fprintf(stderr, "error: CUDA was requested, but not compiled for.\n");
		exit(EXIT_FAILURE);
	}
#endif

	// Load the graph into ForceAtlas2 data structure
	RPGraph::UGraph graph = RPGraph::UGraph(src, dst, num_of_edges, num_of_nodes, communities, degree_S);
	cout << "Graph converted to force atlas" << "\n";

  	// Find graph modularity
	uint32_t qm = 0;
	int sigma = 0;
	for (uint32_t i = 0; i < num_of_edges; i++)
	{
		sigma = communities[src[i]] == communities[dst[i]] ? 1 : 0;
		qm = qm + (1 - ((degree[src[i]] * degree[dst[i]]) / (2 * num_of_edges))) * (sigma);
	}
	cout << "The modularity is: " << qm << "/" << (2 * num_of_edges) << "\n";

	if (avoid_cpu_allocation || cuda_requested) {
		// Free CUDA memory from the SCoDA algorithm for it is no longer needed
		cudaFree(src);
		cudaFree(dst);
		cudaFree(degree);
		cudaFree(Degree_done);
		cudaFree(degree_cmt);
		cudaFree(degree_S);
	} else {
		// Free memory from the SCoDA algorithm for it is no longer needed
		free(src);
		free(dst);
		free(degree);
		free(Degree_done);
		free(degree_cmt);
		free(degree_S);
	}

	// POLES FOR COLORING

	// Convert arg_map["pole_list"] = "12,144,2" to array of ints pole_list[3] = {12, 144, 2}
	int pole_list_size = std::count(arg_map["pole_list"].begin(), arg_map["pole_list"].end(), ',') + 1;
	if (arg_map["pole_list"].size() == 0)
		pole_list_size = 0;
	int* pole_list = new int[pole_list_size];
	split_to_int(arg_map["pole_list"], ',', pole_list, pole_list_size);
	
	int top_N_nodes = stoi(arg_map["top_N_nodes"]);
	if (top_N_nodes > 0) {
		// Replace the pole_list with the top N nodes
		int* top_nodes = graph.get_top_nodes(top_N_nodes);
		delete[] pole_list;
		pole_list = top_nodes;
		pole_list_size = top_N_nodes;
		cout << "Top " << top_N_nodes << "in-degree nodes were selected" << "\n";
	}

	cout << "Pole list: {";
	for (int i = 0; i < pole_list_size; i++)
		cout << pole_list[i] << (i == pole_list_size - 1 ? "" : ",");
	cout << "}\n";
	
	// Create the GraphLayout and ForceAtlas2 objects.
	RPGraph::GraphLayout layout(graph);

	layout.pole_list = pole_list;
	layout.pole_list_size = pole_list_size;
	
	RPGraph::ForceAtlas2* fa2;

	// Choose an appropriate ForceAtlas2 implementation
	// NOTE: Magnets are not supported in the CUDA implementation.
#ifdef __NVCC__
	if(cuda_requested)
		fa2 = new RPGraph::CUDAForceAtlas2(layout, approximate, 
			strong_gravity, gravity, scale);
	else
#endif
	fa2 = new RPGraph::CPUForceAtlas2(layout, approximate,
									  strong_gravity, gravity, scale);

	// FA 2 parameters
	fa2->use_linlog = std::string(arg_map["use_linlog"]) == "true";

	fa2->attraction_exponent = std::stof(arg_map["attraction_exponent"]);
	fa2->k_attraction = std::stof(arg_map["attraction"]);
	fa2->pin_2_roots = std::string(arg_map["pin_2_roots"]) == "true";
	fa2->magetic_pole_separation = std::stof(arg_map["magnetic_pole_separation"]);
	fa2->repulsion_d_squared = std::string(arg_map["repulsion_d_squared"]) == "true";
	fa2->use_pole_segmentation = std::string(arg_map["use_pole_segmentation"]) == "true";
	fa2->pin_poles = std::string(arg_map["pin_poles"]) == "true";
	fa2->extra_pole_attraction = std::stof(arg_map["extra_pole_attraction"]);
	fa2->pole_gravity_factor = std::stof(arg_map["pole_gravity_factor"]);
	fa2->legacy_segmentation = std::string(arg_map["legacy_segmentation"]) == "true";
	fa2->simple_center_of_mass = std::string(arg_map["simple_center_of_mass"]) == "true";

	// COSMETIC PARAMETERS
	layout.setAlphaParameters(std::stof(arg_map["node_alpha"]), std::stof(arg_map["edge_alpha"]));
	layout.square_coordinates = std::string(arg_map["square_coordinates"]) == "true";
	layout.draw_arrows = std::string(arg_map["draw_arrows"]) == "true";
	layout.min_arrow_length = std::stoi(arg_map["min_arrow_length"]);
	layout.draw_common_edges = std::string(arg_map["draw_common_edges"]) == "true";
	layout.pole_size_factor = std::stof(arg_map["pole_size_factor"]);
	layout.colored_fraction = std::stof(arg_map["colored_fraction"]);

	fa2->pole_list = pole_list;
	fa2->pole_list_size = pole_list_size;

	layout.use_distance_based_edge_direction = std::string(arg_map["use_distance_based_edge_direction"]) == "true";
	layout.max_influence_distance = std::stoi(arg_map["max_influence_distance"]);

	// Print the degree of each pole in pole_list
	cout << "Pole in-degrees: {";
	for (int i = 0; i < pole_list_size; i++)
		cout << layout.graph.in_degree(layout.graph.node_map[pole_list[i]]) << (i == pole_list_size - 1 ? "" : ",");
	cout << "}\n";

	// Print pole labels
	if (node_labels.size() > 0) {
		cout << "Pole labels: {";
		for (int i = 0; i < pole_list_size; i++)
			cout << node_labels[pole_list[i]] << (i == pole_list_size - 1 ? "" : ", ");
		cout << "}\n";
	}

	// MAGNETIC FORCE PARAMETERS
	bool use_magnetic_field = std::string(arg_map["use_magnetic_field"]) == "yes" or std::string(arg_map["use_magnetic_field"]) == "true";

	if (use_magnetic_field) {
		std::string field_type = arg_map["field_type"];
		bool bi_directional = std::string(arg_map["bi_directional"]) == "true";
		float field_strength = std::stof(arg_map["field_strength"]);
		float magnetic_constant = std::stof(arg_map["magnetic_constant"]);
		float magnetic_alpha = std::stof(arg_map["magnetic_alpha"]);
		float magnetic_beta = std::stof(arg_map["magnetic_beta"]);

		fa2->setMagneticParameters(field_type, bi_directional, field_strength, magnetic_constant, magnetic_alpha, magnetic_beta);
	}

	cout << "Started ForceAtlas2 layout algorithm..." << "\n";

	const int snap_period = ceil((float) max_iterations / num_screenshots);
	const int print_period = ceil((float) max_iterations * 0.05);

	// Load the working map into CUDA memory
	uint32_t* nodemap;
	cudaMallocManaged(&nodemap, graph.num_nodes() * sizeof(uint32_t));

	for (uint32_t i = 0; i < graph.num_nodes(); i++)
		nodemap[i] = layout.graph.node_map_r[i];

	float cumulative_max_force = 0;
	float last_avg_max_force = std::numeric_limits<float>::max();

	// The procedure for saving screenshots of the layout
	auto saveScreenshot = [&] (int iteration) {
        // Determine output filename
		std::string edgelist_basename = basename(in_path);
		time_t now = time(0);

		// Use a simplified output name in lexico-graphical order
		std::string timestamp = std::to_string(now) + "_";

		if (arg_map["include_timestamp"] != "true")
			timestamp = "";

		std::string out_filename = edgelist_basename + "_" + timestamp + fill_zeros(iteration, 4) + ".";	
		// Example: "tree_edgelist.txt_123456789_0500.png"

		// Or use a config chain name instead
		if (arg_map["chain_output_name"] == "true")
			out_filename = arg_map["config_chain"] + timestamp + fill_zeros(iteration, 4) + ".";
		// Example: "tree m-polar 123456789_0500.png"

		std::string out_filepath = out_path + "/" + out_filename;

		printf("writing %s... ", out_format.c_str());
		fflush(stdout);

		// Make sure to synchronize the layout and retrieve it from the GPU first
		fa2->sync_layout();

		if (out_format == "png")
			layout.writeToPNG(image_w, image_h, out_filepath + "png");
		else if (out_format == "csv")
			layout.writeToCSV(out_filepath + "csv");
		else if (out_format == "bin")
			layout.writeToBin(out_filepath + "bin");
		else if (out_format == "net")
			layout.writeToNET(out_filepath + "net", node_labels);
		else if (out_format == "png+net") {
			layout.writeToPNG(image_w, image_h, out_filepath + "png");
			layout.writeToNET(out_filepath + "net", node_labels);
		}

		printf("done.\n");
    };
	
	for (int iteration = 1; iteration <= max_iterations; iteration++)
	{
		// Run a CUDA kernel to compute the next step of the force layout
		fa2->doStep(nodemap);
		cumulative_max_force += fa2->max_force;

		// If we need to, write the result to a png
		if (num_screenshots > 0 && (iteration % snap_period == 0 || iteration == max_iterations))
		{
			if (arg_map["hide_iteration_output"] != "true")
				printf("Starting iteration %3d (%3.0f%%), ", iteration, 100 * (float)iteration / max_iterations);
			saveScreenshot(iteration);
		}
		// Print out the iteration progress
		else if (iteration % print_period == 0 && arg_map["hide_iteration_output"] != "true")
		{
			printf("Starting iteration %3d (%3.0f%%), ", iteration, 100 * (float)iteration / max_iterations);
			// Print out max force during the last iteration
			float past_avg_max_force =  cumulative_max_force / print_period;
			// If the exponent is larger than 5 render using scientific notation
			if (to_string(past_avg_max_force).length() > 10)
				printf("f = %.1e", past_avg_max_force);
			else
				printf("f = %7.2f", past_avg_max_force);
			// Warning about not converging
			if (past_avg_max_force >= 1e5)
				printf(" (!)");
			
			if (pole_list_size >= 1) {
				float misaligned = fa2->count_misaligned_edges(M_PI/3);
				printf(" err = %2.2f%%", 100.0*misaligned);
			}
			
			printf("\n");
			// If the force is getting closer to 0, it converges
			if (past_avg_max_force > last_avg_max_force * stof(arg_map["divergence_factor"]) && arg_map["stop_on_divergence"] == "true") {
				std::cout << "DIVERGENCE DETECTED, ";
				saveScreenshot(iteration);
				break;
			}
			if (past_avg_max_force > stof(arg_map["divergence_threshold"]) && iteration >= 500 && arg_map["stop_on_divergence"] == "true") {
				std::cout << "DIVERGENCE DETECTED, ";
				saveScreenshot(iteration);
				break;	
			}
				
			last_avg_max_force = past_avg_max_force;
			cumulative_max_force = 0;
		}
	}

	if (pole_list_size >= 1) {
		float misaligned = fa2->count_misaligned_edges(M_PI/3);
		printf("Final err = %2.2f%%\n", 100.0*misaligned);
	}

	delete fa2;

	// Print out the total execution time
	end = high_resolution_clock::now();
	long exec_time = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	cout << "Execution finished in " << exec_time << "ms\n";

	exit(EXIT_SUCCESS);
}
