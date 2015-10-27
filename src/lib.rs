#![feature(test)]
extern crate petgraph;

// XXX: Self-loops?

use petgraph::{Graph, Directed};
use petgraph::graph::{NodeIndex, IndexType};
use std::collections::HashSet;
use std::mem;

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

    /// XXX: Optimize: find_edge() is O(k) with k the number of edges.
    pub fn from_graph<G: DirectedGraph>(graph: &G) -> TriadicCensus {
        let n = graph.node_count();

        let mut census = TriadicCensus::new();

        for v in 0..n {
            let mut neighbors_v = HashSet::new();
            graph.each_undirected_neighbor(v, |u| {
                if u > v {
                    neighbors_v.insert(u);
                }
            });

            for u in neighbors_v {
                let mut s = HashSet::new(); // XXX: Use bitset? // reuse
                graph.each_undirected_neighbor(u, |nu| {
                    s.insert(nu);
                });
                graph.each_undirected_neighbor(v, |nv| {
                    s.insert(nv);
                });
                let _ = s.remove(&u);
                let _ = s.remove(&v);

                let tri_type = if graph.has_edge(v, u) && graph.has_edge(u, v) {
                    TriadType::T102
                } else {
                    TriadType::T012
                };

                let cnt = n as i64 - s.len() as i64 - 2;
                assert!(cnt >= 0);
                census.add(tri_type, cnt as u64);
                for w in s {
                    if u < w ||
                       (v < w && w < u && !(graph.has_edge(v, w) || graph.has_edge(w, v))) {
                        let tri_type = TriadType::from_u8(TRITYPES[tricode(graph, v, u, w)]);
                        census.add(tri_type, 1);
                    }
                }
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

pub struct SimpleDiGraph {
    g: Graph<(), (), Directed>,
}

impl SimpleDiGraph {
    pub fn new() -> SimpleDiGraph {
        SimpleDiGraph { g: Graph::new() }
    }

    pub fn add_node(&mut self) -> NodeIdx {
        self.g.add_node(()).index()
    }

    pub fn add_edge(&mut self, src: NodeIdx, dst: NodeIdx) {
        self.g.add_edge(NodeIndex::new(src), NodeIndex::new(dst), ());
    }

    pub fn from(num_nodes: usize, edge_list: &[(NodeIdx, NodeIdx)]) -> SimpleDiGraph {
        let mut g = SimpleDiGraph::new();
        for _ in 0..num_nodes {
            let _ = g.add_node();
        }
        for &(src, dst) in edge_list {
            g.add_edge(src, dst);
        }
        g
    }
}

impl DirectedGraph for SimpleDiGraph {
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

#[test]
fn test_simple() {
    let mut g = SimpleDiGraph::new();

    let n1 = g.add_node();
    let n2 = g.add_node();
    let n3 = g.add_node();

    let c1 = TriadicCensus::from_graph(&g);
    assert_eq!(&[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c1.as_slice());

    g.add_edge(n1, n3);
    let c2 = TriadicCensus::from_graph(&g);
    assert_eq!(&[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c2.as_slice());

    let d = TriadicCensus::distance(&c1, &c2); // should be sqrt(2)
    assert!(d > 1.414 && d < 1.415);

    g.add_edge(n3, n1);
    let c3 = TriadicCensus::from_graph(&g);
    assert_eq!(&[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c3.as_slice());

    g.add_edge(n2, n3);
    let c4 = TriadicCensus::from_graph(&g);
    assert_eq!(&[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c4.as_slice());
}

#[cfg(test)]
fn line_graph(n: u32) -> SimpleDiGraph {
    let mut g = SimpleDiGraph::new();

    let mut prev = g.add_node();
    for _ in 0..n - 1 {
        let cur = g.add_node();
        g.add_edge(prev, cur);
        prev = cur;
    }
    return g;
}

#[cfg(test)]
fn circular_graph(n: u32) -> SimpleDiGraph {
    let mut g = SimpleDiGraph::new();

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
    let c = TriadicCensus::from_graph(&line_graph(20));
    assert_eq!(&[816, 306, 0, 0, 0, 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c.as_slice());

    let c = TriadicCensus::from_graph(&line_graph(40));
    assert_eq!(&[8436, 1406, 0, 0, 0, 38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c.as_slice());
}

#[test]
fn test_circular_graph() {
    // Result compared with igraph.triad_census().
    let c = TriadicCensus::from_graph(&circular_graph(20));
    assert_eq!(&[800, 320, 0, 0, 0, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c.as_slice());

    let c = TriadicCensus::from_graph(&circular_graph(40));
    assert_eq!(&[8400, 1440, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               c.as_slice());

}

#[cfg(test)]
mod example_graphs;

#[bench]
fn bench_erdos_renyi_graph(b: &mut test::Bencher) {
    let g = SimpleDiGraph::from(example_graphs::GRAPH1_NODES, &example_graphs::GRAPH1_EDGES);
    b.iter(|| {
        let _c = TriadicCensus::from_graph(&g);
    });
}
