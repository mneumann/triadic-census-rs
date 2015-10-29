#![feature(test)]
extern crate petgraph;

// XXX: Self-loops?

use petgraph::{Graph, Directed};
use petgraph::graph::{NodeIndex, IndexType};
use std::collections::HashSet;
use std::mem;
use std::collections::BTreeMap;
use std::cmp;
use std::convert::From;

#[cfg(test)]
extern crate test;

#[repr(u8)]
#[derive(Debug)]
pub enum TriadType {
    T003 = 0,
    T012,
    T102,
    T021D,

    T021U,
    T021C,
    T111D,
    T111U,

    T030T,
    T030C,
    T201,
    T120D,

    T120U,
    T120C,
    T210,
    T300,
}

impl TriadType {
    #[inline(always)]
    fn from_u8(i: u8) -> TriadType {
        assert!(i < 16);
        unsafe { mem::transmute(i) }
    }
}

pub type NodeIdx = usize;

/// Interface required for Triadic Census.
/// Nodes are externally represented by `NodeIdx`.
pub trait DirectedGraph {
    fn node_count(&self) -> NodeIdx;
    fn has_edge(&self, src: NodeIdx, dst: NodeIdx) -> bool;
    fn each_undirected_neighbor<F: FnMut(NodeIdx)>(&self, node: NodeIdx, callback: F);

    // excluding u and v
    fn set_of_undirected_neighbors(&self, u: NodeIdx, v: NodeIdx) -> HashSet<NodeIdx> {
        let mut s = HashSet::new(); // XXX: Use bitset?
        self.each_undirected_neighbor(u, |nu| { s.insert(nu); });
        self.each_undirected_neighbor(v, |nv| { s.insert(nv); });
        let _ = s.remove(&u);
        let _ = s.remove(&v);
        return s;
    }
}

fn tricode<G: DirectedGraph>(graph: &G, v: NodeIdx, u: NodeIdx, w: NodeIdx) -> usize {
    let mut tricode = 0;
    if graph.has_edge(v, u) {
        tricode |= 1 << 0;
    }
    if graph.has_edge(u, v) {
        tricode |= 1 << 1;
    }
    if graph.has_edge(v, w) {
        tricode |= 1 << 2;
    }
    if graph.has_edge(w, v) {
        tricode |= 1 << 3;
    }
    if graph.has_edge(u, w) {
        tricode |= 1 << 4;
    }
    if graph.has_edge(w, u) {
        tricode |= 1 << 5;
    }
    return tricode;
}

const TRITYPES: [u8; 64] = [0, 1, 1, 2, 1, 3, 5, 7, 1, 5, 4, 6, 2, 7, 6, 10, 1, 5, 3, 7, 4, 8, 8,
                            12, 5, 9, 8, 13, 6, 13, 11, 14, 1, 4, 5, 6, 5, 8, 9, 13, 3, 8, 8, 11,
                            7, 12, 13, 14, 2, 6, 7, 10, 6, 11, 13, 14, 7, 13, 12, 14, 10, 14, 14,
                            15];


#[derive(Debug)]
pub struct TriadicCensus {
    census: [u64; 16],
}

impl TriadicCensus {
    /// Calculates the `distance` between two TriadicCensus structs, treated as
    /// points in a 16-dimensional space.
    pub fn distance(a: &TriadicCensus, b: &TriadicCensus) -> f64 {
        let mut sum = 0.0f64;
        for i in 0..16 {
            let va = a.census[i];
            let vb = b.census[i];
            let diff = (if va > vb {
                va - vb
            } else {
                vb - va
            }) as f64;
            sum += diff * diff;
        }
        sum.sqrt()
    }

    #[inline(always)]
    fn new() -> TriadicCensus {
        TriadicCensus { census: [0; 16] }
    }

    #[inline(always)]
    fn add(&mut self, triad_type: TriadType, cnt: u64) {
        self.census[triad_type as usize] += cnt;
    }

    #[inline(always)]
    fn set(&mut self, triad_type: TriadType, cnt: u64) {
        self.census[triad_type as usize] = cnt;
    }

    #[inline(always)]
    pub fn get(&self, triad_type: TriadType) -> u64 {
        self.census[triad_type as usize]
    }

    #[inline(always)]
    pub fn as_slice<'a>(&'a self) -> &'a [u64] {
        &self.census[..]
    }

}

impl<'a, G: DirectedGraph> From<&'a G> for TriadicCensus {
    /// XXX: Optimize: find_edge() is O(k) with k the number of edges.
    fn from(graph: &G) -> Self {
        let n = graph.node_count();

        let mut census = TriadicCensus::new();
        let mut neighbors_v = HashSet::new();

        for v in 0..n {
            neighbors_v.clear();
            graph.each_undirected_neighbor(v, |u| {
                if u > v {
                    neighbors_v.insert(u);
                }
            });

            for &u in neighbors_v.iter() {
                let tri_type = if graph.has_edge(v, u) && graph.has_edge(u, v) {
                    TriadType::T102
                } else {
                    TriadType::T012
                };

                let mut s_cnt = 0;
                let s = graph.set_of_undirected_neighbors(u, v);
                for &w in s.iter() {
                    s_cnt += 1;
                    if u < w ||
                       (v < w && w < u && !(graph.has_edge(v, w) || graph.has_edge(w, v))) {
                        let tri_type = TriadType::from_u8(TRITYPES[tricode(graph, v, u, w)]);
                        census.add(tri_type, 1);
                    }
                }

                let cnt = n as i64 - s_cnt as i64 - 2;
                assert!(cnt >= 0);
                census.add(tri_type, cnt as u64);
            }
        }

        let mut sum = 0;
        for &cnt in &census.as_slice()[1..] {
            sum += cnt;
        }
        let n = n as u64;
        // Integer division below is guaranteed to produce an integral result
        census.set(TriadType::T003, (n * (n - 1) * (n - 2)) / 6 - sum);
        return census;
    }
}

pub struct SimpleDigraph {
    g: Graph<(), (), Directed>,
}

impl SimpleDigraph {
    pub fn new() -> SimpleDigraph {
        SimpleDigraph { g: Graph::new() }
    }

    pub fn add_node(&mut self) -> NodeIdx {
        self.g.add_node(()).index()
    }

    pub fn add_edge(&mut self, src: NodeIdx, dst: NodeIdx) {
        self.g.add_edge(NodeIndex::new(src), NodeIndex::new(dst), ());
    }

    pub fn from(num_nodes: usize, edge_list: &[(NodeIdx, NodeIdx)]) -> SimpleDigraph {
        let mut g = SimpleDigraph::new();
        for _ in 0..num_nodes {
            let _ = g.add_node();
        }
        for &(src, dst) in edge_list {
            g.add_edge(src, dst);
        }
        g
    }
}

impl DirectedGraph for SimpleDigraph {
    fn node_count(&self) -> NodeIdx {
        self.g.node_count()
    }

    #[inline]
    fn has_edge(&self, src: NodeIdx, dst: NodeIdx) -> bool {
        self.g.find_edge(NodeIndex::new(src), NodeIndex::new(dst)).is_some()
    }

    #[inline]
    fn each_undirected_neighbor<F: FnMut(NodeIdx)>(&self, node: NodeIdx, mut callback: F) {
        for n in self.g.neighbors_undirected(NodeIndex::new(node)) {
            callback(n.index())
        }
    }
}

#[inline]
fn calc_index(n: usize, src: NodeIdx, dst: NodeIdx) -> (usize, u64) {
    assert!(src < n && dst < n);
    let idx = src * n + dst;
    let idx_64 = idx / 64;
    let bit_idx = idx % 64;
    (idx_64, 1u64 << bit_idx as u64)
}

pub struct OptSparseDigraph {
    g: SimpleDigraph,
    edges: BTreeMap<usize, u64>,
}

impl OptSparseDigraph {
    pub fn from(graph: SimpleDigraph) -> OptSparseDigraph {
        let n = graph.node_count();

        let mut edges: BTreeMap<usize, u64> = BTreeMap::new();

        for edge in graph.g.raw_edges().iter() {
            let src = edge.source().index();
            let dst = edge.target().index();
        
            // treat as adjacency matrix (with `src` using rows and `dst` columns)
            let (idx, bit_pattern) = calc_index(n, src, dst);

            let entry = edges.entry(idx);
            *entry.or_insert(0) |= bit_pattern;
        }

        OptSparseDigraph {g: graph, edges: edges}
    }
}

impl DirectedGraph for OptSparseDigraph {
    fn node_count(&self) -> NodeIdx {
        self.g.node_count()
    }

    #[inline]
    fn has_edge(&self, src: NodeIdx, dst: NodeIdx) -> bool {
        let (idx, bit_pattern) = calc_index(self.g.node_count(), src, dst);

        match self.edges.get(&idx) {
            Some(&v) => {
                (v & bit_pattern) != 0
            }
            _ => false
        }
    }

    #[inline]
    fn each_undirected_neighbor<F: FnMut(NodeIdx)>(&self, node: NodeIdx, callback: F) {
        self.g.each_undirected_neighbor(node, callback);
    }
}

pub struct OptDenseDigraph {
    g: SimpleDigraph,
    n: usize,
    matrix: Box<[u64]>,
}

impl OptDenseDigraph {
    pub fn from(graph: SimpleDigraph) -> OptDenseDigraph {
        let n = graph.node_count();

        let vec_len = (n * n) / 64 + cmp::min(1, ((n*n) % 64));
        let mut matrix: Vec<u64> = (0..vec_len).map(|_| 0).collect();

        for edge in graph.g.raw_edges().iter() {
            let src = edge.source().index();
            let dst = edge.target().index();
        
            // treat as adjacency matrix (with `src` using rows and `dst` columns)
            let (idx, bit_pattern) = calc_index(n, src, dst);

            matrix[idx] |= bit_pattern;
        }

        OptDenseDigraph {g: graph, n: n, matrix: matrix.into_boxed_slice()}
    }
}

impl DirectedGraph for OptDenseDigraph {
    fn node_count(&self) -> NodeIdx {
        self.n
    }

    #[inline]
    fn has_edge(&self, src: NodeIdx, dst: NodeIdx) -> bool {
        let (idx, bit_pattern) = calc_index(self.n, src, dst);
        (self.matrix[idx] & bit_pattern) != 0
    }

    #[inline]
    fn each_undirected_neighbor<F: FnMut(NodeIdx)>(&self, node: NodeIdx, callback: F) {
        self.g.each_undirected_neighbor(node, callback);
    }
}


#[test]
fn test_simple() {
    let mut g = SimpleDigraph::new();

    let n1 = g.add_node();
    let n2 = g.add_node();
    let n3 = g.add_node();

    let c1 = TriadicCensus::from(&g);
    assert_eq!(&[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c1.as_slice());

    g.add_edge(n1, n3);
    let c2 = TriadicCensus::from(&g);
    assert_eq!(&[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c2.as_slice());

    let d = TriadicCensus::distance(&c1, &c2); // should be sqrt(2)
    assert!(d > 1.414 && d < 1.415);

    g.add_edge(n3, n1);
    let c3 = TriadicCensus::from(&g);
    assert_eq!(&[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c3.as_slice());

    g.add_edge(n2, n3);
    let c4 = TriadicCensus::from(&g);
    assert_eq!(&[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c4.as_slice());
}

#[cfg(test)]
fn line_graph(n: u32) -> SimpleDigraph {
    let mut g = SimpleDigraph::new();

    let mut prev = g.add_node();
    for _ in 0..n - 1 {
        let cur = g.add_node();
        g.add_edge(prev, cur);
        prev = cur;
    }
    return g;
}

#[cfg(test)]
fn circular_graph(n: u32) -> SimpleDigraph {
    let mut g = SimpleDigraph::new();

    let first = g.add_node();
    let mut prev = first;
    for _ in 0..n - 1 {
        let cur = g.add_node();
        g.add_edge(prev, cur);
        prev = cur;
    }
    g.add_edge(prev, first);
    return g;
}

#[test]
fn test_line_graph() {
    // Result compared with igraph.triad_census().
    let c = TriadicCensus::from(&line_graph(20));
    assert_eq!(&[816, 306, 0, 0, 0, 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c.as_slice());

    let c = TriadicCensus::from(&line_graph(40));
    assert_eq!(&[8436, 1406, 0, 0, 0, 38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c.as_slice());
}

#[test]
fn test_circular_graph() {
    // Result compared with igraph.triad_census().
    let c = TriadicCensus::from(&circular_graph(20));
    assert_eq!(&[800, 320, 0, 0, 0, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c.as_slice());

    let c = TriadicCensus::from(&circular_graph(40));
    assert_eq!(&[8400, 1440, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c.as_slice());

}

#[cfg(test)]
mod example_graphs;

#[bench]
fn bench_erdos_renyi_graph(b: &mut test::Bencher) {
    let g = SimpleDigraph::from(example_graphs::GRAPH1_NODES, &example_graphs::GRAPH1_EDGES);
    b.iter(|| {
        let c = TriadicCensus::from(&g);
        assert_eq!(&[2484, 14525, 7954, 7156, 7237, 14346, 15426, 15413, 14041, 4778, 8454, 7492, 7614, 15161, 16641, 2978], c.as_slice());
    });
}

#[bench]
fn bench_erdos_renyi_graph_opt_sparse(b: &mut test::Bencher) {
    let g = SimpleDigraph::from(example_graphs::GRAPH1_NODES, &example_graphs::GRAPH1_EDGES);
    let g = OptSparseDigraph::from(g);
    b.iter(|| {
        let c = TriadicCensus::from(&g);
        assert_eq!(&[2484, 14525, 7954, 7156, 7237, 14346, 15426, 15413, 14041, 4778, 8454, 7492, 7614, 15161, 16641, 2978], c.as_slice());
    });
}

#[bench]
fn bench_erdos_renyi_graph_opt_dense(b: &mut test::Bencher) {
    let g = SimpleDigraph::from(example_graphs::GRAPH1_NODES, &example_graphs::GRAPH1_EDGES);
    let g = OptDenseDigraph::from(g);
    b.iter(|| {
        let c = TriadicCensus::from(&g);
        assert_eq!(&[2484, 14525, 7954, 7156, 7237, 14346, 15426, 15413, 14041, 4778, 8454, 7492, 7614, 15161, 16641, 2978], c.as_slice());
    });
}
