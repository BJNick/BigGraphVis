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
void initial(int numb_of_links,int numb_of_nodes,uint32_t *communities,int *degree,int *degree_cmt,int *degree_S,int *hashtable1,int *hashtable2,int *hashtable3,int *hashtable4,int seed,curandState_t* states, int heunumber)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;


        if(index<numb_of_nodes&&index>0)
        {
            
          communities[index] = index;
          degree[index]=0;
            degree_S[index]=0;
           degree_cmt[index]=0;
       /*  if(index<=numb_of_nodes/heunumber)
         {
              curand_init(seed, index,0,&states[index]);
              hashtable1[index]=curand(&states[index])%2; 
              hashtable2[index]=curand(&states[index])%2; 
              hashtable3[index]=curand(&states[index])%2; 
              hashtable4[index]=curand(&states[index])%2; 
         }*/
             
      
        }
    if(index<numb_of_links&&index>0)
        {
        curand_init(seed, index,0,&states[index]);
              hashtable1[index]=curand(&states[index])%2; 
              hashtable2[index]=curand(&states[index])%2; 
              hashtable3[index]=curand(&states[index])%2; 
              hashtable4[index]=curand(&states[index])%2; 
    }

}
//=========================================================================
__global__
void communities_hashing(int numb_of_links,int numb_of_nodes,uint32_t *communities,int *degree,int *hashtable1,int *hashtable2,int *hashtable3,int *hashtable4,int seed,curandState_t* states,int *h1,int *h2,int *h3,int *h4,int huenumber,uint32_t *dst,uint32_t *src)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index<numb_of_nodes)
     {
        
             h1[communities[index]]=((((communities[index]+7)/3)*2)+1+hashtable1[index])%huenumber;
             h2[communities[index]]=((((communities[index]+16)/4)*3)+1+hashtable2[index])%huenumber;

          
     }
            h3[index]=((((communities[src[index]]+communities[dst[index]])/53)*4)+1+hashtable3[src[index]])%huenumber;
            h4[index]=((((communities[src[index]]+communities[dst[index]])/74)*3)+1+hashtable4[dst[index]])%huenumber;
}
//====================================================================================
/*
__global__
void links_remove(int numb_of_links,int numb_of_nodes,uint32_t *communities,int *degree,int *hashtable1,int *hashtable2,int *hashtable3,int *hashtable4,int seed,curandState_t* states,int *h1,int *h2,int *h3,int *h4,int huenumber,uint32_t *dst,uint32_t *src,int *weight_temp,int *weight_S)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
  if(weight_temp[index]==0)
  {
      weight_temp[index]=weight_S[index];
  }
    else
    {
        src[index]=numb_of_links+10;
        dst[index]=numb_of_links+10;
    }
} */
//=========================================================================

__global__
void Fdegree_S(int numb_of_links,int numb_of_nodes,uint32_t *communities,int *degree,int *degree_S,int *sketch,int *sketch2,int *hashtable1,int *hashtable2,int *hashtable3,int *hashtable4,int seed,curandState_t* states,int *h1,int *h2,int *h3,int *h4,int huenumber,uint32_t *dst,uint32_t *src,int *weight_S)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int a,b,c,d;
    if(index<numb_of_nodes)
     {
        
           a=((((communities[index]+7)/3)*2)+1+hashtable1[index])%huenumber;
            b=((((communities[index]+16)/4)*3)+1+hashtable2[index])%huenumber;
             degree_S[index]= fminf(sketch[a],sketch[b]);
          
     }
    c=((((communities[src[index]]+communities[dst[index]])/53)*4)+1+hashtable3[communities[src[index]]])%huenumber;
    d=((((communities[src[index]]+communities[dst[index]])/74)*3)+1+hashtable4[index])%huenumber;
    weight_S[index]=fminf(sketch2[a],sketch2[b]);
}
//=========================================================================
__global__
void mainfor(uint32_t num_of_links,uint32_t num_of_nodes,int degree_threshold,int *degree_cmt,uint32_t *src,uint32_t *dst,uint32_t *communities)
{
 int index = blockIdx.x * blockDim.x + threadIdx.x;
    if(index>0)
    if(src[index]>0&&src[index]<=num_of_nodes&&dst[index]>0&&dst[index]<=num_of_nodes)
        if( degree_cmt[src[index]] <= degree_threshold && degree_cmt[dst[index]]<= degree_threshold )
        {
            
            if( degree_cmt[src[index]] > degree_cmt[dst[index]] )
            {
                communities[dst[index]] = communities[src[index]];
               atomicAdd(&degree_cmt[dst[index]],1);
               // src[index]=0;
                
            } 
            else 
            {
                communities[src[index]] = communities[dst[index]]; 
                atomicAdd(&degree_cmt[src[index]],1);
                 //dst[index]=0;
            }
                
               
        }
 
}

//=========================================================================
__global__
void communities_sketch(int num_of_links,int num_of_nodes,uint32_t *communities,int *degree,int *hashtable1,int *hashtable2,int *hashtable3,int *hashtable4,int seed,curandState_t* states,int *h1,int *h2,int *h3,int *h4,int huenumber, int *Degree_done,int *sketch,int *sketch2,uint32_t *src,uint32_t *dst,int *weight_src,int *weight_dst,int *weight,int num_of_threads)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int weight_src_flag=0;
    int weight_dst_flag=0;
    int x=num_of_links/num_of_threads;
//************************************* each threads will cout only its section
//************************************* each threads will only count the data for it self
for(int i=index*x;i<=(index*x)+x;i++)
{
//************************************* counting outside links ignoring intra links
   

//************************************* Sketching the source of links
    
     if( h1[communities[src[i]]]<=(index*x)+x && h1[communities[src[i]]]>=(index*x) )
    //if( h1[i]<=(index*x)+x && h1[i]>=(index*x) )
     {
         if(Degree_done[src[i]]==0)// check if the node is visited before
         {
             sketch[h1[communities[src[i]]]]+=degree[src[i]];// instrad of increasing by one! we increased it by the degree of the node
             Degree_done[src[i]]=1; // the node is visited 
         }
     }
     if( h2[communities[src[i]]]<=(index*x)+x && h2[communities[src[i]]]>=(index*x) )
    // if( h2[i]<=(index*x)+x && h2[i]>=(index*x) )
     {
         if(Degree_done[src[i]]==0) 
         {
             sketch[h2[communities[src[i]]]]+=degree[src[i]];
             Degree_done[src[i]]=1;
         }
     }        
        
           if( h1[communities[src[i]]]<=(index*x)+x && h1[communities[src[i]]]>=(index*x) )
   // if( h1[i]<=(index*x)+x && h1[i]>=(index*x) )
     {
         if(Degree_done[dst[i]]==0)
         {
             sketch[h1[communities[dst[i]]]]+=degree[dst[i]];
             Degree_done[dst[i]]=1;
         }
     }
    if( h2[communities[src[i]]]<=(index*x)+x && h2[communities[src[i]]]>=(index*x) )
    // if( h2[i]<=(index*x)+x && h2[i]>=(index*x) )
     {
         if(Degree_done[dst[i]]==0)
         {
             sketch[h2[communities[dst[i]]]]+=degree[dst[i]];
             Degree_done[dst[i]]=1;
         }
     }   
//************************************* Sketching the destination of links
       if( communities[src[i]] != communities[dst[i]] )
    {
     if( h3[src[i]]<=(index*x)+x && h3[src[i]]>=(index*x) || h3[dst[i]]<=(index*x)+x && h3[dst[i]]>=(index*x))
     {
             sketch2[h3[i]]++;
     }
     if( h4[src[i]]<=(index*x)+x && h4[src[i]]>=(index*x) || h4[dst[i]]<=(index*x)+x && h4[dst[i]]>=(index*x))
     {
             sketch2[h4[i]]++;
     }

    }

}
    
    
}

//***********************************************************
//============================================================
//***********************************************************

int main(int argc, const char **argv)
{
    cout<<"The Algorithm started "<< "\n"; 
    // For reproducibility.
    srandom(1234);
    // Parse commandline arguments
    if (argc < 14)
    {
        fprintf(stderr, "Usage: graph_viewer gpu|cpu max_iterations num_snaps sg|wg scale gravity exact|approximate edgelist_path out_path [png image_w image_h|csv|bin]\n");
        exit(EXIT_FAILURE);
    }auto start = high_resolution_clock::now();     auto end = high_resolution_clock::now();
curandState_t* states;
uint32_t *communities,*src,*dst,*counters_matrix_weight,degree_threshold,num_of_links,num_of_nodes;
int huenumber=100,src_id, dst_id,weight_id,rounds,*hashtable1,*hashtable2,*hashtable3,*hashtable4,*sketch,*sketch2,*h1,*h2,*h3,*h4,*weight_dst,*weight_src,*weight,*Degree_done,num_of_threads,*degree,*degree_cmt,*degree_S,*weight_S,*weight_temp;
std::string s = argv[13];
degree_threshold = std::atoi(s.c_str());
num_of_nodes = 0; num_of_links = 0;
s=argv[14];rounds = std::atoi(s.c_str());
ifstream inFile; inFile.open(argv[8]);
cout<<"Initialization"<< "\n"; 
for(int i=0;inFile >>src_id >>dst_id != NULL ;i++)//reading the file for counting nodes and links
 {
   num_of_links++; if(src_id>num_of_nodes) num_of_nodes=src_id; if(dst_id>num_of_nodes) num_of_nodes=dst_id;
 } 
inFile.close(); cout<<"Number of nodes: "<<num_of_nodes<<" number of Links: "<<num_of_links<<'\n';
         
int blockSize = 256;
int remain=num_of_nodes%blockSize;
int numBlocks = (num_of_nodes + blockSize - remain) / blockSize;
int number_of_threads=(numBlocks*blockSize);  
    
cudaMallocManaged(&src, num_of_links*sizeof(uint32_t));
cudaMallocManaged(&dst, num_of_links*sizeof(uint32_t));
cudaMallocManaged(&weight_dst, num_of_links*sizeof(int));
cudaMallocManaged(&weight_src, num_of_links*sizeof(int));
cudaMallocManaged(&Degree_done, num_of_links*sizeof(int));
cudaMallocManaged(&weight, num_of_links*sizeof(int));
cudaMallocManaged(&degree, num_of_nodes*sizeof(int));  
cudaMallocManaged(&degree_cmt, num_of_nodes*sizeof(int));  
cudaMallocManaged(&degree_S, num_of_nodes*sizeof(int));  
cudaMallocManaged(&weight_S, num_of_links*sizeof(int));  
cudaMallocManaged(&weight_temp, num_of_links*sizeof(int));  

cudaMallocManaged(&communities, num_of_nodes*sizeof(uint32_t));
cudaMallocManaged(&h1, num_of_nodes*sizeof(int));
cudaMallocManaged(&h2, num_of_nodes*sizeof(int));
cudaMallocManaged(&h3, num_of_links*sizeof(int));
cudaMallocManaged(&h4, num_of_links*sizeof(int));
cudaMallocManaged(&sketch, huenumber*sizeof(int));
cudaMallocManaged(&sketch2, huenumber*sizeof(int));
cudaMallocManaged(&hashtable1, num_of_nodes*sizeof(int));  
cudaMallocManaged(&hashtable2, num_of_nodes*sizeof(int));  
cudaMallocManaged(&hashtable3, num_of_links*sizeof(int));  
cudaMallocManaged(&hashtable4, num_of_links*sizeof(int));  
cudaMallocManaged((void**) &states, num_of_nodes*sizeof(curandState_t));
//################################################################# INITIALIZATION
remain=num_of_links%blockSize;
numBlocks = (num_of_links + blockSize - remain) / blockSize;
initial<<<numBlocks, blockSize>>>(num_of_links,num_of_nodes,communities,degree,degree_cmt,degree_S,hashtable1,hashtable2,hashtable3,hashtable4,time(0),states,huenumber); 
cudaDeviceSynchronize();
cudaError_t error = cudaGetLastError();  if(error != cudaSuccess) { printf("CUDA error in LINE 201: %s\n", cudaGetErrorString(error)); exit(-1); }   
//************************************************************** READING NETWORK
inFile.open(argv[8]); for(int i=1;inFile >>src_id >>dst_id != NULL ;i++) {src[i]=src_id; dst[i]=dst_id; degree[src_id]++; if(dst_id!=src_id) degree[dst_id]++; }
start = high_resolution_clock::now();   
    cout<<"Community Detection . . . "<< "\n"; 
//################################################################# SCODA MAIN LOOP
if(rounds==0)
 rounds=num_of_links*0.01;
remain=(num_of_links)%blockSize;
numBlocks = (num_of_links + blockSize - remain) / blockSize;
for(int i=0;i<rounds;i++)
 {
   mainfor<<<numBlocks, blockSize>>>(num_of_links,num_of_nodes,degree_threshold,degree_cmt,src,dst,communities); cudaDeviceSynchronize(); error = cudaGetLastError(); if(error != cudaSuccess) { printf("CUDA error: %s\n", cudaGetErrorString(error)); exit(-1); }
cout<<"Community Detection Round  "<<i+1<<"\n"; 
}
communities[0]=1;
//#################################################################Hash Function
cout<<"   Hashing . . . "<< "\n"; 
communities_hashing<<<numBlocks, blockSize>>>( num_of_links, num_of_nodes, communities, degree, hashtable1, hashtable2, hashtable3, hashtable4, time(0), states, h1, h2, h3, h4, huenumber, dst, src); cudaDeviceSynchronize(); error = cudaGetLastError(); if(error != cudaSuccess) { printf("CUDA error in LINE 218: %s\n", cudaGetErrorString(error)); exit(-1); }

//################################################################# min  Sketch
cout<<"   min sketch . . . "<< "\n"; 

remain=(num_of_links/huenumber)%blockSize;
numBlocks = (num_of_links/huenumber + blockSize - remain) / blockSize;
numBlocks=8;    
num_of_threads=numBlocks*blockSize;
communities_sketch<<<numBlocks, blockSize>>>( num_of_links, num_of_nodes, communities, degree, hashtable1, hashtable2, hashtable3, hashtable4, time(0),states, h1, h2, h3, h4, huenumber,  Degree_done, sketch, sketch2, src, dst, weight_src, weight_dst, weight, num_of_threads); cudaDeviceSynchronize(); error = cudaGetLastError(); if(error != cudaSuccess) { printf("CUDA error in LINE 224: %s\n", cudaGetErrorString(error)); exit(-1); }
//cout<<numBlocks<<" "<<blockSize;
  
        //cout<<weight[i]<<"***SRC"<<weight_src[i]<<"***DST"<<weight_dst[i]<<"***Wsrc"<<sketch[weight_src[i]]<<"***Wdst"<<sketch[weight_dst[i]]<<'\n';
//#######################################################
cout<<"Counting . . . "<< "\n"; 
remain=(num_of_links/huenumber)%blockSize;
numBlocks = (num_of_links/huenumber + blockSize - remain) / blockSize;
    Fdegree_S<<<numBlocks, blockSize>>>(num_of_links,num_of_nodes,communities,degree,degree_S,sketch,sketch2,hashtable1,hashtable2,hashtable3,hashtable4,time(0),states,h1,h2,h3,h4,huenumber,dst,src,weight_S); cudaDeviceSynchronize(); error = cudaGetLastError(); if(error != cudaSuccess) { printf("CUDA error in LINE 248: %s\n", cudaGetErrorString(error)); exit(-1); }
//#######################################################################
  /*  cout<<"Removing Duplicate links . . . "<< "\n"; 

 remain=(num_of_links/huenumber)%blockSize;
numBlocks = (num_of_links/huenumber + blockSize - remain) / blockSize;
    links_remove<<<numBlocks, blockSize>>>(num_of_links,num_of_nodes,communities,degree,hashtable1,hashtable2,
                                           hashtable3,hashtable4,time(0),states,h1, h2,h3,h4,huenumber,dst,src,weight_temp,weight_S); cudaDeviceSynchronize(); error = cudaGetLastError(); if(error != cudaSuccess) { printf("CUDA error in LINE 248: %s\n", cudaGetErrorString(error)); exit(-1); }
*/
//####################################### FORCEATLAS2
cout<<"The Force Atlas 2"<< "\n"; 

 /*  for(int i=0;i<num_of_nodes;i++){
       if(degree_S[i]>0){
       cout<<"degree "<<degree_S[i]<<'\n';
       cout<<"h1 "<<h1[i]<<" "<<h2[i]<<'\n';}
   }*/
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
    int control=0;
  
int32_t num_superlink=0;
ofstream myfile;
myfile.open (argv[4]);
RPGraph::UGraph graph = RPGraph::UGraph(src,dst,num_of_links,communities,degree_S);
   
  cudaFree(src);
  cudaFree(dst);
  cudaFree(degree);
  cudaFree(communities);
        

    
    
    
    
    
    
    
    
    
    printf("Loading edgelist at '%s'...", edgelist_path.c_str());
    fflush(stdout);
    //RPGraph::UGraph graph = RPGraph::UGraph(edgelist_path);
    printf("done.\n");
    printf("    fetched %d nodes and %d edges.\n", graph.num_nodes(), graph.num_edges());

    // Create the GraphLayout and ForceAtlas2 objects.
    RPGraph::GraphLayout layout(graph);
    RPGraph::ForceAtlas2 *fa2;
    #ifdef __NVCC__
    if(cuda_requested)
        fa2 = new RPGraph::CUDAForceAtlas2(layout, approximate,
                                           strong_gravity, gravity, scale);
    else
    #endif
        fa2 = new RPGraph::CPUForceAtlas2(layout, approximate,
                                          strong_gravity, gravity, scale);

    printf("Started Layout algorithm...\n");
    const int snap_period = ceil((float)max_iterations/num_screenshots);
    const int print_period = ceil((float)max_iterations*0.05);

    for (int iteration = 1; iteration <= max_iterations; ++iteration)
    {
        fa2->doStep();
        // If we need to, write the result to a png
        if (num_screenshots > 0 && (iteration % snap_period == 0 || iteration == max_iterations))
        {
	    // Determine output filename
	    std::string edgelist_basename = basename(edgelist_path);
	    std::string out_filename = edgelist_basename + "_" + std::to_string(iteration) + "." + out_format;
            std::string out_filepath = out_path + "/" + out_filename;
            printf("Starting iteration %d (%.2f%%), writing %s...", iteration, 100*(float)iteration/max_iterations, out_format.c_str());
            fflush(stdout);
            fa2->sync_layout();

            if (out_format == "png")
                layout.writeToPNG(image_w, image_h, out_filepath,hashtable3,hashtable4,sketch2, huenumber);
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
 cout<<"Time of the  algorithm"<<chrono::duration_cast<chrono::seconds>(end - start).count()<<"\n";
    
    
    //*****************
    exit(EXIT_SUCCESS);
}
