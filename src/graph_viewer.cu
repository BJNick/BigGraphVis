//in code ghable taghir be code scoda hast
//in kar mikone ghable dastkarie 30-08-2020
#include <iostream>
using namespace std; 
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
using namespace std::chrono; 
#include "RPCommon.hpp"
#include "RPGraph.hpp"
#include "RPGraphLayout.hpp"
#include "RPCPUForceAtlas2.hpp"
#include <algorithm> 
#ifdef __NVCC__
#include <cuda_runtime_api.h>
#include "RPGPUForceAtlas2.hpp"
#endif

__global__
void initial(int numb_of_links,int numb_of_nodes,uint32_t *communities,uint32_t *degree,uint32_t *degree_cmt,uint32_t *degree_S, int heunumber,int *sketch)
{
 int index = blockIdx.x * blockDim.x + threadIdx.x;


          if(index<numb_of_nodes)
        {
          
          communities[index] = index;
          degree[index]=0;
		  degree_S[index]=1;
		  degree_cmt[index]=0;
		  sketch[index]=1;
		}
          
}
//=========================================================================
__global__
void communities_hashing(int numb_of_links,int numb_of_nodes,uint32_t *communities,uint32_t *degree,int huenumber,uint32_t *dst,uint32_t *src,uint32_t degree_thresholdS,int *sketch,int *Degree_done)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index<numb_of_nodes)
     {
	//hash 1 FVN-- hash2 ELF --hash3 BKDR--hash4 DJB
	const unsigned int fnv_prime = 0x811C9DC5;
	unsigned int hash1=0,hash2 = 0,hash3=0,hash4 = 5381;
	uint32_t community=0;
	unsigned int x = 0;
    int deg = 0;
	unsigned int seed = 131;

	if(Degree_done[src[index]]==0||Degree_done[dst[index]]==0)// check if the node is visited before
         {
			if(Degree_done[src[index]]==0)
			{
			community=communities[src[index]]+1;
			deg=int(degree[src[index]]/degree_thresholdS);
			Degree_done[src[index]]=1;
			}
			else
			{
				community=communities[dst[index]]+2;
				deg=int(degree[dst[index]]/degree_thresholdS);
				Degree_done[dst[index]]=1;
			}
	while(community>10)
	{
		hash1 *= fnv_prime;
		hash1 ^= (community%10);
		
		hash2 = (hash2 << 4) + (community%10);
		if ((x = hash2 & 0xF0000000) != 0)
		{
			hash2 ^= (x >> 24);
		}
		hash2 &= ~x;
		
		hash3 = (hash3 * seed) + (community);
		
		hash4 = ((hash4 << 5) + hash4) + (community);
		
		community=uint32_t(community/10);
	}
	sketch[int((hash1%1000)%int(huenumber/4))]+=deg+1;
	sketch[int((hash2%1000)%(int(huenumber/4)+int(huenumber/4)))]+=deg+1;
	sketch[int((hash3%1000)%(int(huenumber/4)+int(huenumber/2)))]+=deg+1;
	sketch[int((hash4%1000)%(int(huenumber/4)+int((3*huenumber)/4)))]+=deg+1;

		
     }
          
}}
//***************************************
__global__
void Fdegree_S(int numb_of_links,int numb_of_nodes,uint32_t *communities,uint32_t *degree,uint32_t *degree_S,int *sketch,int huenumber,uint32_t *dst,uint32_t *src)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if(index<numb_of_nodes)
     {
        	//hash 1 FVN-- hash2 ELF --hash3 BKDR--hash4 DJB
	const unsigned int fnv_prime = 0x811C9DC5;
	unsigned int hash1=0,hash2 = 0,hash3=0,hash4 = 5381;
	uint32_t community=communities[index]+1;
	unsigned int x = 0,deg=0;
	unsigned int seed = 131;

	while(community>10)
	{
		hash1 *= fnv_prime;
		hash1 ^= (community%10);
		
		hash2 = (hash2 << 4) + (community%10);
		if ((x = hash2 & 0xF0000000) != 0)
		{
			hash2 ^= (x >> 24);
		}
		hash2 &= ~x;
		
		hash3 = (hash3 * seed) + (community);
		
		hash4 = ((hash4 << 5) + hash4) + (community);
		
		community=uint32_t(community/10);
	}
	degree_S[communities[index]]= fminf
						(
							fminf
								(
									sketch[int(hash1%int(huenumber/4))],sketch[int(hash2%(int(huenumber/4)+int(huenumber/4)))]
								 )
						,fminf
								(
									sketch[int(hash3%(int(huenumber/4)+int(huenumber/2)))],sketch[int(hash4%(int(huenumber/4)+int((3*huenumber)/4)))]
								)
						);
						
          if(degree_S[communities[index]]==0)
			  degree_S[communities[index]]=1;
     }

}
//=========================================================================
__global__
void mainfor(uint32_t num_of_links,uint32_t num_of_nodes,int degree_threshold,uint32_t *degree_cmt,uint32_t *src,uint32_t *dst,uint32_t *communities)
{
 int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index>0)
    if(src[index]>0&&src[index]<=num_of_nodes&&dst[index]>0&&dst[index]<=num_of_nodes)
        if( degree_cmt[src[index]] <= degree_threshold && degree_cmt[dst[index]]<= degree_threshold )
        {
            
            if( degree_cmt[src[index]] > degree_cmt[dst[index]] )
            {
                communities[dst[index]] = communities[src[index]];
               atomicAdd(&degree_cmt[src[index]],1);
                
            } 
            else 
            {
                communities[src[index]] = communities[dst[index]]; 
               atomicAdd(&degree_cmt[dst[index]],1);
            }
                
               
        }
 
}

//=========================================================================


//***********************************************************
//============================================================
//***********************************************************

int main(int argc, const char **argv)
{

            cudaSetDevice(1);

    cout<<"The Algorithm started "<< "\n"; 
    // For reproducibility.
    srandom(1234);
    // Parse commandline arguments
    if (argc < 15)
    {
        fprintf(stderr, "Usage: graph_viewer gpu|cpu max_iterations num_snaps sg|wg scale gravity exact|approximate edgelist_path out_path [png image_w image_h|csv|bin]\n");
        exit(EXIT_FAILURE);
    }
auto start = high_resolution_clock::now();     
auto end = high_resolution_clock::now(); 
    auto end_cmt = high_resolution_clock::now(); 
curandState_t* states;
uint32_t *communities,*src,*dst,degree_threshold,degree_thresholdS,num_of_links,num_of_nodes,*degree,*degree_cmt,*degree_S,degree_name;
int huenumber=6000,src_id, dst_id,rounds,*sketch,*h1,*h2,*h3,*h4,*Degree_done,*weight_S;
std::string sss=argv[15];
std::string src_id_s,dst_id_s;
huenumber=std::stoul(sss.c_str());
std::string s = argv[13];
degree_threshold = std::stoul(s.c_str());
degree_thresholdS=degree_threshold;
    degree_name=degree_threshold;
num_of_nodes = 0; 
num_of_links = 0;
s=argv[14];
rounds = std::atoi(s.c_str());
ifstream inFile; 
inFile.open(argv[8]);
cout<<"Initialization"<< "\n"; 
for(int i=0;inFile >>src_id_s >>dst_id_s != NULL ;i++)//reading the file for counting nodes and links
 {
	 if(isdigit(src_id_s[0]))
     {
		 src_id=stol(src_id_s); 
         dst_id=stol(dst_id_s);
         num_of_links++; 
         if(src_id>num_of_nodes)
             num_of_nodes=src_id; 
         if(dst_id>num_of_nodes) 
             num_of_nodes=dst_id;
     } 
}

inFile.close(); 
cout<<"Number of nodes: "<<num_of_nodes<<" number of Links: "<<num_of_links<<'\n';
int blockSize = 256;
int remain=num_of_nodes%blockSize;
int numBlocks = (num_of_nodes + blockSize - remain) / blockSize;
cout<<"line 251"<<'\n';    
cudaMallocManaged(&src, num_of_links*sizeof(uint32_t));
cudaMallocManaged(&dst, num_of_links*sizeof(uint32_t));
cudaMallocManaged(&Degree_done, num_of_links*sizeof(int));
cudaMallocManaged(&degree, num_of_nodes*sizeof(uint32_t));  
cudaMallocManaged(&degree_cmt, num_of_nodes*sizeof(uint32_t));  
cudaMallocManaged(&degree_S, num_of_nodes*sizeof(uint32_t));  
cudaMallocManaged(&weight_S, num_of_links*sizeof(int));  
cout<<"line 263"<<'\n';    
cudaMallocManaged(&communities, num_of_nodes*sizeof(uint32_t));
cudaMallocManaged(&h1, num_of_nodes*sizeof(int));
cudaMallocManaged(&h2, num_of_nodes*sizeof(int));
cudaMallocManaged(&h3, num_of_links*sizeof(int));
cudaMallocManaged(&h4, num_of_links*sizeof(int));
cudaMallocManaged(&sketch, huenumber*sizeof(int));

cudaMallocManaged((void**) &states, num_of_links*sizeof(curandState_t));
cout<<"line 276"<<'\n';    
//################################################################# INITIALIZATION
remain=num_of_links%blockSize;
numBlocks = (num_of_links + blockSize - remain) / blockSize;
initial<<<numBlocks, blockSize>>>( num_of_links,num_of_nodes,communities,degree,degree_cmt,degree_S,  huenumber,sketch); 
    cout<<"LINE 281"<<'\n';
cudaDeviceSynchronize();
cudaError_t error = cudaGetLastError();  
    if(error != cudaSuccess) 
    { 
        printf("CUDA error in LINE 282: %s\n", cudaGetErrorString(error)); 
     exit(-1); 
    }  
    
//************************************************************** READING NETWORK
    cout<<"line 284"<<'\n';    
inFile.open(argv[8]); 
int adjmat=0;
for(int i=1;inFile >>src_id_s >>dst_id_s != NULL ;i++)
	{
        if(isdigit(src_id_s[0]))
        {
                src_id=stol(src_id_s);
                dst_id=stol(dst_id_s);
          
                src[adjmat]=src_id;
                dst[adjmat]=dst_id; 
                adjmat++;
                degree[src_id]++;
                if(dst_id!=src_id) 
                    degree[dst_id]++;          
        }
    }
//start = high_resolution_clock::now();   

    cout<<"Community Detection . . . "<< "\n"; 
//################################################################# SCODA MAIN LOOP
if(rounds==0)
 rounds=num_of_links*0.01;
remain=(num_of_links)%blockSize;
numBlocks = (num_of_links + blockSize - remain) / blockSize;
for(int i=1;i<rounds;i++)
 {
    //if(degree_threshold<num_of_links*.05)
    degree_threshold=degree_threshold*degree_threshold;
   mainfor<<<numBlocks, blockSize>>>(num_of_links,num_of_nodes,degree_threshold,degree_cmt,src,dst,communities); cudaDeviceSynchronize(); error = cudaGetLastError(); if(error != cudaSuccess) { printf("CUDA error: %s\n", cudaGetErrorString(error)); exit(-1); }
cout<<"Community Detection Round  "<<i+1<<"\n"; 
}
//communities[0]=1;
//#################################################################Hash Function
cout<<"   Hashing . . . "<< "\n"; 
communities_hashing<<<numBlocks, blockSize>>>(  num_of_links, num_of_nodes,communities,degree, huenumber,dst,src,degree_thresholdS,sketch,Degree_done); cudaDeviceSynchronize(); error = cudaGetLastError(); if(error != cudaSuccess) { printf("CUDA error in LINE 218: %s\n", cudaGetErrorString(error)); exit(-1); }

//################################################################# min  Sketch
cout<<"Counting . . . "<< "\n"; 
remain=(num_of_links/huenumber)%blockSize;
numBlocks = (num_of_links/huenumber + blockSize - remain) / blockSize;
    Fdegree_S<<<numBlocks, blockSize>>>( num_of_links, num_of_nodes,communities,degree,degree_S,sketch, huenumber,dst,src); cudaDeviceSynchronize(); error = cudaGetLastError(); if(error != cudaSuccess) { printf("CUDA error in LINE 248: %s\n", cudaGetErrorString(error)); exit(-1); }
//#######################################################################
 //#################################### Writing of communities
    std::ofstream out_file("../../../files/communities.txt");

        for ( uint32_t n = 0; n < num_of_nodes; ++n)
        {
            
            out_file << n << " " << communities[n] << "\n";
        }

        out_file.close();
        end_cmt = high_resolution_clock::now();   
 cout<<"Time of the  CMT algorithm"<<chrono::duration_cast<chrono::milliseconds>(end_cmt - start).count()<<"\n";
//####################################### FORCEATLAS2
cout<<"The Force Atlas 2"<< "\n"; 

    const bool cuda_requested = std::string(argv[1]) == "gpu" or std::string(argv[1]) == "cuda";
    const int max_iterations = std::stoi(argv[2]);
    const int num_screenshots = std::stoi(argv[3]);
    const bool strong_gravity = std::string(argv[4]) == "sg";
    const float scale = std::stof(argv[5]);
    const float gravity = std::stof(argv[6]);
    const bool approximate = std::string(argv[7]) == "approximate";
    std::string edgelist_path = argv[8];
    std::string out_path = argv[9];
    std::string out_format = "png";
    int image_w = 1250;
    int image_h = 1250;

    for (int arg_no = 10; arg_no < argc; arg_no++)
    {
        if(std::string(argv[arg_no]) == "png")
        {
            out_format = "png";
            image_w = std::stoi(argv[arg_no+1]);
            image_h = std::stoi(argv[arg_no+2]);
            arg_no += 2;
        }

        else if(std::string(argv[arg_no]) == "csv")
        {
            out_format = "csv";
        }

        else if(std::string(argv[arg_no]) == "bin")
        {
            out_format = "bin";
        }
    }


    if(cuda_requested and not approximate)
    {
        fprintf(stderr, "error: The CUDA implementation (currently) requires Barnes-Hut approximation.\n");
        exit(EXIT_FAILURE);
    }

    // Check in_path and out_path
    if (!is_file_exists(edgelist_path))
    {
        fprintf(stderr, "error: No edgelist at %s\n", edgelist_path.c_str());
        exit(EXIT_FAILURE);
    }
    if (!is_file_exists(out_path))
    {
        fprintf(stderr, "error: No output folder at %s\n", out_path.c_str());
        exit(EXIT_FAILURE);
    }

    // If not compiled with cuda support, check if cuda is requested.
    #ifndef __NVCC__
    if(cuda_requested)
    {
        fprintf(stderr, "error: CUDA was requested, but not compiled for.\n");
        exit(EXIT_FAILURE);
    }
    #endif

    // Load graph.
//    int control=0;
  
//int32_t num_superlink=0;
ofstream myfile;
myfile.open (argv[4]);

RPGraph::UGraph graph = RPGraph::UGraph(src,dst,num_of_links,num_of_nodes,communities,degree_S);
  // RPGraph::UGraph graph = RPGraph::UGraph(src,dst,num_of_links,num_of_nodes,communities,degree_cmt);
//***********************
//******************************
//Moudularity
uint32_t qm=0;
int sigm=0;
for(uint32_t i=0;i<num_of_links;i++)
{
    if(communities[src[i]]==communities[dst[i]])
    sigm=1;
    else
    sigm=0;
qm=qm+(1-((degree[src[i]]*degree[dst[i]])/(2*num_of_links)))*(sigm);
}
cout<<"The moudularity is : "<<qm<<"/"<<(2*num_of_links)<<"\n";

//end of modularity




  cudaFree(src);
  cudaFree(dst);
  cudaFree(degree);
      cudaFree(Degree_done);
    cudaFree(degree_cmt);
    cudaFree(degree_S);

         

    
    
    
    
    
    
    
    
    
    printf("Loading edgelist at '%s'...", edgelist_path.c_str());
    fflush(stdout);
    //RPGraph::UGraph graph = RPGraph::UGraph(edgelist_path);
    printf("done.\n");
    printf("    fetched %d nodes and %d edges.\n", graph.num_nodes(), graph.num_edges());
if(graph.num_nodes()<0||graph.num_nodes()>num_of_nodes)
     exit(EXIT_SUCCESS);
    // Create the GraphLayout and ForceAtlas2 objects.
    RPGraph::GraphLayout layout(graph);
    RPGraph::ForceAtlas2 *fa2;
  //  #ifdef __NVCC__
    //if(cuda_requested)//commented
        fa2 = new RPGraph::CUDAForceAtlas2(layout, approximate,
                                           strong_gravity, gravity, scale);
    /*else//commented
    #endif
        fa2 = new RPGraph::CPUForceAtlas2(layout, approximate,
                                          strong_gravity, gravity, scale);*/

    printf("Started Layout algorithm...\n");
    const int snap_period = ceil((float)max_iterations/num_screenshots);
    const int print_period = ceil((float)max_iterations*0.05);
    //cout<<"*************"<< layout.graph.node_map_r[100]<<"\n";
    uint32_t *nodemap;
    cudaMallocManaged(&nodemap, graph.num_nodes()*sizeof(uint32_t));
    for(uint32_t i=0;i<graph.num_nodes();i++)
        nodemap[i]=layout.graph.node_map_r[i];
    for (int iteration = 1; iteration <= max_iterations; ++iteration)
    {
      
        fa2->doStep(nodemap);
        // If we need to, write the result to a png
        if (num_screenshots > 0 && (iteration % snap_period == 0 || iteration == max_iterations))
        {
	    // Determine output filename
	    std::string edgelist_basename = basename(edgelist_path);
		time_t now = time(0);
	    std::string out_filename = edgelist_basename+ to_string(rounds)+"_"+to_string(degree_name)+ "_" + std::to_string(iteration)+"_"+std::to_string(huenumber)+"_" +std::to_string(now) +"." + out_format;
            std::string out_filepath = out_path + "/" + out_filename;
            printf("Starting iteration %d (%.2f%%), writing %s...", iteration, 100*(float)iteration/max_iterations, out_format.c_str());
            fflush(stdout);
            fa2->sync_layout();

            if (out_format == "png")
                layout.writeToPNG(image_w, image_h, out_filepath);
            else if (out_format == "csv")
                layout.writeToCSV(out_filepath);
            else if (out_format == "bin")
                layout.writeToBin(out_filepath);

            printf("done.\n");
        }

        // Else we print (if we need to)
        else if (iteration % print_period == 0)
        {
            printf("Starting iteration %d (%.2f%%).\n", iteration, 100*(float)iteration/max_iterations);
        }
    }
    
    
    
    

    delete fa2;
    end = high_resolution_clock::now();   
 cout<<"Time of the  algorithm"<<chrono::duration_cast<chrono::milliseconds>(end - start).count()<<"\n";
    

    exit(EXIT_SUCCESS);
}
