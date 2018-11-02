#ifndef PATH_COVER_H
#define PATH_COVER_H

#include <iostream>
#include <vector>
#include <algorithm>

#include <lemon/smart_graph.h>
#include <lemon/matching.h>

/*
input:
n m
x y w
...
*/

#define DERR \
    cerr << "[" << __FILE__ << " " << __LINE__ << "]: "

namespace loon
{

//typedef long long Weight_T;
typedef double Weight_T;

class PathCover
{
public:
    class Path
    {
    public:
        int path_type;
        std::vector<int> edges;
    public:
        int getPathType() const
        {
            /*
            1: singleton
            2: path
            3: cycle
            */
            return path_type;
        }
        void swap(Path& p)
        {
            std::swap(this->path_type, p.path_type);
            this->edges.swap( p.edges );
        }
    };

    int n;
    lemon::SmartGraph g;
    lemon::SmartGraph::EdgeMap<Weight_T> w;
    std::vector<Path> paths;
    
public:
    PathCover();
    PathCover(int n, int m);
    void reserveNode(int n);
    void reserveEdge(int m);
    void genNodes();
    void addEdge(int i, int j, Weight_T weight);
    void run();
};

}
#endif
