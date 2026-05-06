#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

// Macros:
// Convert coordinates to index
#define xyz2index(x, y, z, DIMS) ( (z-1)*DIMS[0]*DIMS[1] + (y-1)*DIMS[0] + (x-1) )

// index (0-based) -> [x, y, z] (1-based)
inline std::vector<int> index2xyz(int index,
                                  Rcpp::IntegerVector& DIMS)
{
    std::vector<int> XYZ;
    XYZ.reserve(3);
    XYZ.push_back( index % DIMS[0] + 1 );
    XYZ.push_back( (index/DIMS[0]) % DIMS[1] + 1 );
    XYZ.push_back( index/(DIMS[0]*DIMS[1]) + 1 );
    
    return XYZ;
}

// check voxel validity
inline bool xyz_check(int x, int y, int z, int index,
                      Rcpp::IntegerVector& DIMS,
                      Rcpp::IntegerVector& MASK)
{
    return (x >= 1 && x <= DIMS[0] && \
            y >= 1 && y <= DIMS[1] && \
            z >= 1 && z <= DIMS[2] && \
            MASK[index] != 0);
}

// Find valid neighbours of a voxel
Rcpp::IntegerVector findNeighbours(Rcpp::IntegerVector& MASK,  // 3D mask of original orders (1:m)
                                   Rcpp::IntegerVector& DIMS,  // image dimensions
                                   int index, int conn)
{
    // compute [x y z] coordinates based on voxel index
    std::vector<int> XYZ = index2xyz(index, DIMS);
    
    // xyz coordinate adjustment vectors
    int DX[26] = {1,-1,0, 0,0, 0,  1,-1, 1,-1,1,-1, 1,-1,0, 0, 0, 0,  1,-1, 1,-1, 1,-1, 1,-1};
    int DY[26] = {0, 0,1,-1,0, 0,  1, 1,-1,-1,0, 0, 0, 0,1,-1, 1,-1,  1, 1,-1,-1, 1, 1,-1,-1};
    int DZ[26] = {0, 0,0, 0,1,-1,  0, 0, 0, 0,1, 1,-1,-1,1, 1,-1,-1,  1, 1, 1, 1,-1,-1,-1,-1};
    
    // find all valid neighbours of a voxel
    Rcpp::IntegerVector IDS(conn);
    int len = 0;
    for (int i = 0; i < conn; i++)
    {
        int      id = xyz2index(XYZ[0]+DX[i], XYZ[1]+DY[i], XYZ[2]+DZ[i], DIMS);
        bool inmask = xyz_check(XYZ[0]+DX[i], XYZ[1]+DY[i], XYZ[2]+DZ[i], id, DIMS, MASK);
        if (inmask)
        {
            IDS[len] = MASK[id];
            len++;
        }
    }
    
    if (len < conn) IDS.erase(len, conn);
    
    return IDS;
}

// Find the adjacency list for all in-mask voxels
// [[Rcpp::export]]
Rcpp::List findAdjList(Rcpp::IntegerVector& MASK,    // 3D mask of original orders (1:m)
                       Rcpp::IntegerVector& INDEXP,  // voxel indices of unsorted p-values
                       Rcpp::IntegerVector& DIMS,    // image dimensions
                       int m, int conn)
{
    Rcpp::List ADJ(m);
    for (int i = 0; i < m; i++)
    {
        Rcpp::IntegerVector IDS = findNeighbours(MASK, DIMS, INDEXP[i], conn);
        ADJ[i] = IDS;
    }
    
    return ADJ;
}

// ---------------- FIND ----------------
// Find with path halving
int Find(int x, std::vector<int>& PARENT)
{
    while (PARENT[x] != x)
    {
        PARENT[x] = PARENT[PARENT[x]];
        x         = PARENT[x];
    }
    
    return x;
}

// ---------------- UNION ----------------
// Union by size
void UnionBySize(int i, int j,
                 std::vector<int>& PARENT,
                 std::vector<int>& FORESTROOT,
                 std::vector<int>& SIZE)
{
    int irep = Find(i, PARENT);
    int jrep = Find(j, PARENT);
    
    // if i and j are already in the same set
    if (irep == jrep) return;
    
    // i and j are not in same set, so we merge
    int iroot = FORESTROOT[irep];
    int jroot = FORESTROOT[jrep];
    if (SIZE[iroot] < SIZE[jroot])
    {
        PARENT[irep] = jrep;
        FORESTROOT[jrep] = iroot;
    }
    else
    {
        PARENT[jrep] = irep;
    }
    SIZE[iroot] += SIZE[jroot];
}

// =====================================================
// Find voxel at which the formed STC has extent > k
// =====================================================
// [[Rcpp::export]]
int FindIndex(int k, int m, Rcpp::List& ADJ,  // a list of neighbours for all nodes (unsorted!)
              Rcpp::IntegerVector& INDEXP,    // voxel indices of unsorted p-values
              Rcpp::IntegerVector& ORD,       // sorted orders for non-decreasing p-values
              Rcpp::IntegerVector& RANK)      // sorting ranks for all p-values
{
    std::vector<int> SIZE(m, 1);
    std::vector<int> PARENT(m);
    std::vector<int> FORESTROOT(m);
    
    for (int i = 0; i < m; i++)
    {
        PARENT[i] = i;
        FORESTROOT[i] = i;
    }
    
    for (int i = 0; i < m; i++)
    {
        int v = ORD[i]-1;
        
        // find neighbours for node with the ith smallest p-value
        Rcpp::IntegerVector IDS = ADJ[v];
        
        // loop through all its neighbours
        for (int j = 0; j < IDS.size(); j++)
        {
            if (RANK[IDS[j]-1] < i+1)  // check if the neighbour has a smaller rank
            {
                int jrep = Find(IDS[j]-1, PARENT);  // representative of the tree
                int    w = FORESTROOT[jrep];        // forest root of the tree
                
                if (v != w) UnionBySize(v, jrep, PARENT, FORESTROOT, SIZE);
            }
        }
        
        if (SIZE[v] > k) return INDEXP[v] + 1;
    }
    
    return 0;
}

// =====================================================
// Find max(size) for STCs formed with a z-threshold
// =====================================================
// [[Rcpp::export]]
int FindSize(int rnk, int m, Rcpp::List& ADJ,  // rnk: rank corresponding to the given p- or z-threshold
             Rcpp::IntegerVector& ORD,
             Rcpp::IntegerVector& RANK)
{
    std::vector<int> SIZE(m, 1);
    std::vector<int> PARENT(m);
    std::vector<int> FORESTROOT(m);

    for (int i = 0; i < m; i++)
    {
        PARENT[i] = i;
        FORESTROOT[i] = i;
    }

    int maxsize = 0;
    
    // loop through voxels until rnk
    for (int i = 0; i < rnk; i++)
    {
        int v = ORD[i]-1;
        
        // find neighbours for node with the ith smallest p-value
        Rcpp::IntegerVector IDS = ADJ[v];
        
        // loop through all its neighbours
        for (int j = 0; j < IDS.size(); j++)
        {
            if (RANK[IDS[j]-1] < i+1)  // check if the neighbour has a smaller rank
            {
                int jrep = Find(IDS[j]-1, PARENT);  // representative of the tree
                int    w = FORESTROOT[jrep];        // forest root of the tree
                
                if (v != w) UnionBySize(v, jrep, PARENT, FORESTROOT, SIZE);
            }
        }
        
        if (SIZE[v] > maxsize) maxsize = SIZE[v];
    }
    
    return maxsize;
}
