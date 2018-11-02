#include "path_cover.h"

namespace loon
{
PathCover::PathCover(): w(g)
{
}

PathCover::PathCover(int n, int m): w(g)
{
    this->n = n;
    g.reserveNode(n);
    g.reserveEdge(m);
}

void PathCover::reserveNode(int n)
{
    this->n = n;
    g.reserveNode(n);
}

void PathCover::reserveEdge(int m)
{
    g.reserveEdge(m);
}

void PathCover::genNodes()
{
    for(int i = 0; i < n; ++i)
        g.addNode();
}

void PathCover::addEdge(int i, int j, Weight_T weight)
{
    w[ g.addEdge( g.nodeFromId(i), g.nodeFromId(j) )] = weight;
}

void PathCover::run()
{
    lemon::MaxWeightedMatching<lemon::SmartGraph, lemon::SmartGraph::EdgeMap<Weight_T> > mc(g, w);
    mc.run();

#ifdef DEBUG_MC
    std::cerr << std::endl;
    std::DERR << "output matching:" << std::endl;
    std::vector<bool> matched(n, false);
    for(lemon::SmartGraph::NodeIt u(g); u != lemon::INVALID; ++u)
        if(! matched[ g.id(u) ])
        {
            lemon::SmartGraph::Arc arc = mc.matching(u);
            if(arc == lemon::INVALID)
                continue;
            int id1 = g.id( g.source(arc) ), id2 = g.id( g.target(arc));
            matched[ id1 ] = matched[ id2 ] = true;
            std::DERR << id1 << " --- " << id2 << ": " << w[arc] << std::endl;
        }
    std::cerr << std::endl;
#endif
    // compute paths

    std::vector<int> tmp_path_fwd, tmp_path_bwd;
    int min_fwd;
    Weight_T min_weight;
    std::vector<bool> visited(n, false);

    for(lemon::SmartGraph::NodeIt u(g); u != lemon::INVALID; ++u)
    {
        int u_id = g.id(u);
        int u_pair_id = u_id ^ 1;
        if(! visited[ u_id ])
        {
            bool found_cycle = false;
            visited[ u_id ] = visited[ u_pair_id ] = true;
            // forward
            lemon::SmartGraph::Node v = mc.mate(u);
            while(v != lemon::INVALID)
            {
                lemon::SmartGraph::Arc arc = mc.matching( g.nodeFromId(u_id) );
                if(tmp_path_fwd.empty())
                {
                    min_fwd = 0;
                    min_weight = w[arc];
                }
                else if(w[arc] < min_weight)
                {
                    min_fwd = tmp_path_fwd.size();
                    min_weight = w[arc];
                }


                int v_id = g.id(v);
                tmp_path_fwd.push_back(u_id);
                tmp_path_fwd.push_back(v_id);

                u_id = v_id ^ 1;
                visited[u_id] = visited[v_id] = true;
                if(u_id == tmp_path_fwd[0])
                {
                    found_cycle = true;
                    break;
                }
                v = mc.mate( g.nodeFromId(u_id) );
            }

            if(found_cycle)
            {
                paths.push_back( Path() );
                paths.back().path_type = 3;
                // handle cycle
                std::vector<int>& p = paths.back().edges;
                p.insert(p.end(), tmp_path_fwd.begin() + (min_fwd+2), tmp_path_fwd.end());
                p.insert(p.end(), tmp_path_fwd.begin(), tmp_path_fwd.begin() + min_fwd);
                tmp_path_fwd.clear();
                continue;
            }

            // backward
            v = mc.mate( g.nodeFromId(u_pair_id) );
            while(v != lemon::INVALID)
            {
                int v_pair_id = g.id(v);
                tmp_path_bwd.push_back( u_pair_id );
                tmp_path_bwd.push_back( v_pair_id );

                u_pair_id = v_pair_id ^ 1;
                visited[ u_pair_id ] = visited[ v_pair_id ] = true;
                v = mc.mate( g.nodeFromId(u_pair_id) );
            }

            // merge paths
            if(tmp_path_fwd.empty() && tmp_path_bwd.empty())
            {
                paths.push_back( Path() );
                paths.back().path_type = 1;
                paths.back().edges.push_back( u_id );
            }
            else
            {
                paths.push_back( Path() );
                paths.back().path_type = 2;
                if(!tmp_path_bwd.empty())
                {
                    std::reverse(tmp_path_bwd.begin(), tmp_path_bwd.end());
                    paths.back().edges.swap( tmp_path_bwd );
                }
                if(!tmp_path_fwd.empty())
                {
                    std::vector<int>& p = paths.back().edges;
                    p.insert(p.end(), tmp_path_fwd.begin(), tmp_path_fwd.end());
                    tmp_path_fwd.clear();
                }
            }
        }
    }

#ifdef DEBUG_MC
    std::cerr << std::endl;
    std::DERR << "output paths" << std::endl;
    for(size_t i = 0; i < paths.size(); ++i)
    {
        
        std::DERR << "path[" << i << "]: type = " << paths[i].path_type << std::endl;
        if(paths[i].path_type  == 1)
        {
            int u_id = paths[i].edges[0];
            std::cerr << "(" << u_id << ", " << (u_id ^ 1) << ")" << std::endl;
        }
        else
        {
            std::cerr << "(" << (paths[i].edges[0]^1) << ", ";
            for(size_t j = 0; j < paths[i].edges.size(); j+=2)
            {
                std::cerr << paths[i].edges[j] << ") ---- (" << paths[i].edges[j+1] << ", ";
            }
            std::cerr << (paths[i].edges.back() ^ 1) << ")" << std::endl;
        }
    }
    std::cerr << std::endl;
#endif
}
}
