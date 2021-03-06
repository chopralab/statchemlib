/* This is graph.hpp and is part of StatChemLIB
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#ifndef GRAPH_H
#define GRAPH_H
#include <algorithm>
#include <functional>
#include <memory>
#include <queue>
#include <queue>
#include <set>
#include "statchem/graph/mnts.hpp"
#include "statchem/graph/ullsubstate.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/helper/error.hpp"
#include "statchem/helper/smiles.hpp"
#include "statchem/molib/it.hpp"

namespace statchem {

namespace graph {
typedef std::vector<std::vector<bool>> AdjacencyMatrix;

template <class Vertex>
class Graph : public molib::template_vector_container<Vertex*, Vertex> {
   public:
    typedef std::set<Vertex*> VertexSet;
    typedef std::vector<Vertex*> Path;
    typedef std::vector<Path> Cliques;
    typedef std::set<VertexSet> Cycles;
    typedef std::map<const Vertex*, std::vector<size_t>> VertexRingMap;
    typedef std::pair<std::vector<node_id>, std::vector<node_id>>
        MatchedVertices;
    typedef std::vector<MatchedVertices> Matches;

   private:
    std::vector<std::unique_ptr<Vertex>> __vertices;  // if graph owns the
                                                      // vertices, they (the
                                                      // unique_ptrs) are stored
                                                      // here
    AdjacencyMatrix __conn;
    std::vector<int> __num_edges;
    void __expand(Vertex&, Path&, Cycles&, VertexSet&);
    template <class Vertex2>
    bool __match(std::vector<node_id>&, std::vector<node_id>&, Matches&,
                 UllSubState<Graph<Vertex>, Graph<Vertex2>>*) const;
    template <typename T>
    void __init(const T&, const bool);

   public:
    Graph() {}
    template <typename T>
    Graph(const T& vertices, const bool ict = false) {
        __init(vertices, ict);
    }
    Graph(std::vector<std::unique_ptr<Vertex>> vertices, const bool ict,
          const bool)
        : __vertices(std::move(vertices)) {  // here graph owns the vertices
        __init(__vertices, ict);
    }
    void init_conn() {
        __conn.resize(this->size());
        for (auto& row : __conn) row.resize(this->size(), false);
    }
    void set_conn(node_id i, node_id j) {
        __conn[i][j] = true;
        __conn[j][i] = true;
    }
    bool get_conn(node_id i, node_id j) const { return __conn[i][j]; }
    const Vertex& vertex(size_t i) const { return *__vertices[i]; }
    int get_num_edges(const node_id i) const {
        return __num_edges[i];
    }  // real number of edges of a vertex
    std::string get_smiles() const;
    Cycles find_cycles_connected_graph();
    Cycles find_fused_rings();
    Cycles find_rings();
    VertexRingMap vertex_rings();
    Cliques max_weight_clique(const int);
    template <class Vertex2>
    Matches match(const Graph<Vertex2>&) const;
    bool isomorphic(Graph& g) {
        Matches m = match(g);
        if (m.size() > 0 && m[0].first.size() == g.size() &&
            this->size() == g.size())
            return true;
        return false;
    }
    template <class P>
    friend std::ostream& operator<<(std::ostream& stream, const Graph<P>& g);
};

template <typename T>
T intersection(T p1, T p2) {
    T v;
    set_intersection(p1.begin(), p1.end(), p2.begin(), p2.end(),
                     inserter(v, v.begin()));
    return v;
}

template <class Vertex>
template <class T>  // vertices is any container of unique_ptr<Vertex> or
                    // Vertex*
                    void Graph<Vertex>::__init(const T& vertices,
                                               const bool ict) {
    std::map<const Vertex*, node_id> idx;
    for (auto& v : vertices) {
        idx[&*v] = this->size();
        this->add(&*v);
#ifdef STATCHEM_DEBUG_MESSAGES
        for (auto& adj_v : *v) {
            //~ dbgmsg("pv1 = " << &*v << " (" << v << ") " << " pv2 = " <<
            //&adj_v);
            dbgmsg("pv1 = " << &*v << " pv2 = " << &adj_v);
            //~ dbgmsg("v1 = " << v->atom_number() << " v2 = " <<
            //adj_v.atom_number());
        }
#endif
    }
    if (ict) {  // init connection table if requested
        init_conn();
        __num_edges.resize(this->size());
        for (auto& v : *this) {
            const node_id i1 = idx[&v];
            for (auto& adj_v : v) {
                if (idx.count(&adj_v)) {
                    const node_id i2 = idx.at(&adj_v);
                    dbgmsg("new2");
                    dbgmsg("__init : set_conn " << i1 << " " << i2 << " v1 "
                                                << v << " v2 " << adj_v);
                    set_conn(i1, i2);
                    __num_edges[i1]++;  // this counts the real number of edges
                                        // if this is subgraph
                }
            }
        }
    }
    dbgmsg("exiting __init");
}

template <class Vertex>
typename Graph<Vertex>::Path reconstruct_path(
    Vertex& goal, const std::map<Vertex*, Vertex*>& came_from) {
    typename Graph<Vertex>::Path path;
    path.push_back(&goal);
    dbgmsg("path = ");
    dbgmsg(*path.back());
    while (came_from.find(path.back()) != came_from.end()) {
        path.push_back(came_from.at(path.back()));
        dbgmsg(*path.back());
    }
    dbgmsg("----------------------");
    return path;
}

template <class Vertex>
typename Graph<Vertex>::Path find_path(Vertex& start, Vertex& goal) {
    std::queue<Vertex*> openset;
    openset.push(&start);
    typename Graph<Vertex>::VertexSet closedset;
    std::map<Vertex*, Vertex*> came_from;
    while (!openset.empty()) {
        Vertex& curr = *openset.front();
        openset.pop();
        dbgmsg(curr);
        closedset.insert(&curr);
        if (&curr == &goal) {
            dbgmsg("finished successfully " << curr << " == " << goal);
            return reconstruct_path(curr, came_from);
        }
        for (auto& adj_v : curr) {
            if (closedset.find(&adj_v) == closedset.end()) {
                came_from[&adj_v] = &curr;
                openset.push(&adj_v);
            }
        }
    }
    return typename Graph<Vertex>::Path();  // failure - return empty path
}

template <class Vertex>
template <class Vertex2>
bool Graph<Vertex>::__match(
    std::vector<node_id>& c1, std::vector<node_id>& c2, Matches& m,
    UllSubState<Graph<Vertex>, Graph<Vertex2>>* s) const {
    if (s->IsGoal()) {
        int n = s->CoreLen();
        s->GetCoreSet(c1, c2);
#ifdef STATCHEM_DEBUG_MESSAGES
        for (auto& v : c1) dbgmsg(v);
#endif
        m.push_back(std::pair<std::vector<node_id>, std::vector<node_id>>(
            std::vector<node_id>(c1.begin(), c1.begin() + n),
            std::vector<node_id>(c2.begin(), c2.begin() + n)));
        //~ MatchedVertices mv;
        //~ for (int i = 0; i < n; ++i) {
        //~ node_id nid1 = c1[i];
        //~ node_id nid2 = c2[i];
        //~ mv[c1[i]] = c2[i];
        //~ }
        //~ m.push_back(mv);
        return false;
    }
    if (s->IsDead()) return false;
    node_id n1 = NULL_NODE, n2 = NULL_NODE;
    while (s->NextPair(n1, n2, n1, n2)) {
        dbgmsg("__match : n1 = " << n1 << " n2 = " << n2);
        if (s->IsFeasiblePair(n1, n2)) {
            UllSubState<Graph<Vertex>, Graph<Vertex2>>* s1 = s->Clone();
            s1->AddPair(n1, n2);
            dbgmsg("__match : n1 = "
                   << n1 << " n2 = " << n2 << " s1->IsGoal() = "
                   << s1->IsGoal() << " c1.size() = " << c1.size()
                   << " c2.size() = " << c2.size());
            if (__match(c1, c2, m, s1)) {
                s1->BackTrack();
                delete s1;
                return true;
            } else {
                s1->BackTrack();
                delete s1;
            }
        }
    }
    return false;
}

template <class Vertex>
template <class Vertex2>
typename Graph<Vertex>::Matches Graph<Vertex>::match(
    const Graph<Vertex2>& other) const {
    UllSubState<Graph<Vertex>, Graph<Vertex2>> s0(*this, other);
    const Graph<Vertex>& g1 = s0.GetGraph1();
    const Graph<Vertex2>& g2 = s0.GetGraph();
    /* Choose a conservative dimension for the arrays */
    int n;
    if (g1.size() < g2.size())
        n = g2.size();
    else
        n = g1.size();
    dbgmsg(n);
    dbgmsg(g1);
    dbgmsg(g2);
    std::vector<node_id> c1(n), c2(n);
    Matches m;
    __match(c1, c2, m, &s0);
    return m;
}

template <class Vertex>
typename Graph<Vertex>::Cliques Graph<Vertex>::max_weight_clique(
    const int iter) {
    std::unique_ptr<int[]> weight(new int[this->size()]);
    for (int i = 0; i < this->size(); ++i)
        weight[i] = this->element(i).weight();
    std::vector<std::vector<int>> qmax;
    MNTS m(qmax, __conn, weight.get(), iter);
    Cliques clique;
    for (auto& rows : qmax) {
        dbgmsg("found max weight clique of " << std::to_string(rows.size())
                                             << " vertices");
        clique.push_back(Path());
        for (auto& vnum : rows) {
            clique.back().push_back(&this->operator[](vnum));
            dbgmsg("clique push vertex = " << vnum);
        }
    }
    return clique;
}

//~ template<class Vertex>
//~ typename Graph<Vertex>::Cliques Graph<Vertex>::max_clique() {
//~ Benchmark::reset();
//~ vector<vector<int>> qmax;
//~ Maxclique m(__conn, this->size());
//~ int qmax_sz;
//~ m.mcq(qmax, qmax_sz);
//~ Cliques clique;
//~ for (auto &rows : qmax) {
//~ dbgmsg("found max weight clique of " << to_string(rows.size())
//~ << " vertices");
//~ clique.push_back(Path());
//~ for (auto &vnum : rows) {
//~ clique.back().push_back(&this->operator[](vnum));
//~ dbgmsg("clique push vertex = " << vnum);
//~ }
//~ }
//~ cout << "time to find max.weight clique " << Benchmark::seconds_from_start()
//<< " wallclock seconds" << endl;
//~ return clique;
//~ }

template <class Vertex>
typename Graph<Vertex>::Cycles Graph<Vertex>::find_fused_rings() {
    Cycles fused;
    Cycles cycles = find_cycles_connected_graph();
    bool mergeable = true;
    // merge cycles until no more can be merged
    while (mergeable) {
        mergeable = false;
        for (typename Cycles::iterator it = cycles.begin(); it != cycles.end();
             ++it) {
            const VertexSet& first = *it;
            VertexSet current(first.begin(), first.end());
            typename Cycles::iterator it2 = it;
            for (++it2; it2 != cycles.end();) {
                const VertexSet& second = *it2;
                VertexSet inter;
                set_intersection(first.begin(), first.end(), second.begin(),
                                 second.end(), inserter(inter, inter.begin()));
                // merge first with second if > 1 vertices in common
                if (inter.size() > 1) {
                    mergeable = true;
                    current.insert(second.begin(), second.end());
                    it2 = cycles.erase(it2);
                } else {
                    ++it2;
                }
            }
            fused.insert(current);
        }
        if (mergeable) {
            cycles = fused;
            fused.clear();
        }
    }
#ifdef STATCHEM_DEBUG_MESSAGES
    dbgmsg("+++++++++++++FUSED RINGS++++++++++++++++++");
    for (auto& cycle : fused) {
        dbgmsg("new fused ring : ");
        for (auto& pv : cycle)
            //~ dbgmsg(pv->print());
            //~ dbgmsg(*pv);
            dbgmsg(pv->get_label());
    }
#endif
    return fused;
}

template <class Vertex>
typename Graph<Vertex>::Cycles Graph<Vertex>::find_rings() {
    Cycles rings;
    Cycles cycles = find_cycles_connected_graph();
    std::vector<VertexSet> v(cycles.begin(), cycles.end());
    sort(v.begin(), v.end(), [](const VertexSet& i, const VertexSet& j) {
        return i.size() < j.size();
    });
#ifndef NDEBUG
    for (typename std::vector<VertexSet>::iterator it = v.begin(); it != v.end();
         it++) {
        dbgmsg("size = " << it->size());
        for (auto it2 = it->begin(); it2 != it->end(); it2++)
            //~ dbgmsg((*it2)->print() << " ");
            //~ dbgmsg(**it2);
            dbgmsg((*it2)->get_label() << " ");
    }
#endif
    VertexSet mapped;
    // go over cycles sorted by increasing size
    for (typename std::vector<VertexSet>::iterator it = v.begin();
         it != v.end(); it++) {
        VertexSet& cycle = *it;
        VertexSet inter;
        // find cycles that are not made out of smaller cycles
        set_intersection(cycle.begin(), cycle.end(), mapped.begin(),
                         mapped.end(), inserter(inter, inter.begin()));
        if (inter.size() == cycle.size()) continue;
        mapped.insert(cycle.begin(), cycle.end());
        rings.insert(cycle);
    }
    return rings;
}

template <class Vertex>
typename Graph<Vertex>::VertexRingMap Graph<Vertex>::vertex_rings() {
    VertexRingMap ring_map;
    Cycles rings = find_rings();

    for (auto& vert : *this) {
        ring_map[&vert] = std::vector<size_t>();
        for (const auto& ring : rings) {
            if (ring.count(&vert)) ring_map[&vert].push_back(ring.size());
        }
    }

    return ring_map;
}

template <class Vertex>
void Graph<Vertex>::__expand(Vertex& v, Path& path, Cycles& cycles,
                             VertexSet& visited) {
    path.push_back(&v);
    visited.insert(&v);
    dbgmsg("new vertex");
    for (auto& adj_v : v) {
        if (!visited.count(&adj_v)) {
            __expand(adj_v, path, cycles, visited);
        } else {
            typename Path::iterator it = find(path.begin(), path.end(), &adj_v);
            VertexSet cycle(it, path.end());
            if (cycle.size() > 2) {
                cycles.insert(cycle);
            }
        }
    }
    visited.erase(&v);
    path.pop_back();
}

template <class Vertex>
typename Graph<Vertex>::Cycles Graph<Vertex>::find_cycles_connected_graph() {
    Path path;
    Cycles cycles;
    VertexSet visited;
    dbgmsg("expanding the vertices");
    __expand(this->first(), path, cycles, visited);
    dbgmsg("after expanding the vertices");
    dbgmsg("CYCLES FOUND : ");
#ifdef STATCHEM_DEBUG_MESSAGES
    for (auto& cycle : cycles) {
        dbgmsg("CYCLE : ");
        for (auto& pvertex : cycle) {
            dbgmsg(*pvertex);
        }
    }
#endif
    //~ dbgmsg("CYCLES FOUND : " << cycles);
    return cycles;
}

template <class Vertex>
std::string Graph<Vertex>::get_smiles() const {
    std::stringstream ss;
    std::map<Vertex*, int> idx;
    for (auto& v : *this) ss << v.get_label() << " ";
    return ss.str();
}

template <class P>
std::ostream& operator<<(std::ostream& stream, const Graph<P>& g) {
    std::map<const P*, int> idx;
    int i = 0;
    for (auto& vertex : g) idx[&vertex] = i++;
    for (auto& vertex : g) {
        stream << "vertex label = " << vertex.get_label()
               << " index = " << idx.at(&vertex) << " edges = [";
        for (auto& adj_v : vertex) {
            if (idx.count(&adj_v)) {
                stream << "{label = " << adj_v.get_label()
                       << " index = " << idx.at(&adj_v) << "} ";
            }
        }
        stream << "]" << std::endl;
    }
    if (!g.__conn.empty()) {
        int edge_size = 0;
        for (auto& v1 : g) {
            for (auto& v2 : g) {
                //~ if (g.__conn->size() > 0 && (*g.__conn)[idx[&v1]][idx[&v2]]
                //== true) edge_size++;
                if (g.__conn[idx[&v1]][idx[&v2]] == true) edge_size++;
            }
        }
        stream << "p " << g.size() << " " << edge_size << std::endl;
        for (size_t i = 0; i < g.size(); i++) {
            for (size_t j = i + 1; j < g.size(); j++) {
                //~ if (g.__conn->size() > 0 && (*g.__conn)[i][j] == true)
                //stream << "e " << i + 1 << " " << j + 1 << endl;
                if (g.__conn[i][j] == true)
                    stream << "e " << i + 1 << " " << j + 1 << std::endl;
            }
        }
    }
    for (auto& vertex : g) {
        stream << "w " << idx[&vertex] << " " << vertex.weight() << std::endl;
    }
    return stream;
}
}
}

#endif
