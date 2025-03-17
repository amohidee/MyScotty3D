#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>

#include <iostream>

#pragma optimize( "", off )
/*
 * triangulate: split all non-boundary faces into triangles.
 *
 * Works on all valid meshes.
 */
void Halfedge_Mesh::triangulate() {
	//A2G1: triangulation
	std::cout << "num faces " << faces.size() << std::endl;
	for (FaceRef f = faces.begin(); f != faces.end(); ++f) {
		std::cout << "f id= " << f->id << std::endl;
		std::cout << "f degrees= " << f->degree() << std::endl;
		std::cout << "f is boundary= " << f->boundary << std::endl;
		if (f->degree() <= 3 || f->boundary) {
			continue;
		}
		HalfedgeRef h_start = f->halfedge;
		std::cout << "h start position:: " << h_start->vertex->position << std::endl;
		HalfedgeRef h = h_start;
		FaceRef curr_face = f;
		
		VertexRef base_vertex = h_start->vertex;
		
		while (h->next->next->next != h_start) {
			// create new 
			std::cout << "h position:: " << h->vertex->position << std::endl;
			EdgeRef e_new = emplace_edge();
			HalfedgeRef h1 = emplace_halfedge();
			HalfedgeRef h2 = emplace_halfedge();
			FaceRef new_face = emplace_face();
			
			e_new->halfedge = h1;

			// connect new halfedges
			
			h1->set_tnvef(h2, h->next->next, base_vertex, e_new, f);
			h2->set_tnvef(h1, h, h->next->next->vertex, e_new, new_face);
			
			// connect geometry for new face
			h->face = new_face;
			h->next->face = new_face;
			h->next->next = h2;
			new_face->halfedge = h;
			f->halfedge = h1;
			
			// set new h
			h = h1;
			std::cout << "new face: id = " << new_face->id << ", degrees = " << new_face->degree() << std::endl;

		}
		h->next->next->next = h;
	
	}
	
}
#pragma optimize( "", on )

/*
 * linear_subdivide: split faces into quads without moving anything.
 *
 * Works on all valid meshes.
 *
 * (NOTE: uses catmark_subdivide_helper for subdivision)
 */
void Halfedge_Mesh::linear_subdivide() {
	std::unordered_map< VertexCRef, Vec3 > vertex_positions;
	std::unordered_map< EdgeCRef, Vec3 > edge_vertex_positions;
	std::unordered_map< FaceCRef, Vec3 > face_vertex_positions;

	//A2G2: linear subdivision

	// For every vertex, assign its current position to vertex_positions[v]:
	for (VertexRef v = vertices.begin(); v != vertices.end(); ++v) {
        vertex_positions[v] = v->position;
    }

	//(TODO)

    // For every edge, assign the midpoint of its adjacent vertices to edge_vertex_positions[e]:
	// (you may wish to investigate the helper functions of Halfedge_Mesh::Edge)

	for (EdgeRef e = edges.begin(); e != edges.end(); ++e) {
        edge_vertex_positions[e] = e->center();
    }

	//(TODO)

    // For every *non-boundary* face, assign the centroid (i.e., arithmetic mean) to face_vertex_positions[f]:
	// (you may wish to investigate the helper functions of Halfedge_Mesh::Face)

	for (FaceRef f = faces.begin(); f != faces.end(); ++f) {
        if (!f->boundary) {
            face_vertex_positions[f] = f->center();
        }
    }

	//(TODO)


	//use the helper function to actually perform the subdivision:
	catmark_subdivide_helper(vertex_positions, edge_vertex_positions, face_vertex_positions);
}

/*
 * catmark_subdivide: split faces into quads with positions calculated by
 *   the Catmull-Clark ruleset.
 *
 * Works on all valid meshes.
 *
 * (NOTE: uses catmark_subdivide_helper for subdivision)
 */
void Halfedge_Mesh::catmark_subdivide() {
	std::unordered_map< VertexCRef, Vec3 > vertex_positions;
	std::unordered_map< EdgeCRef, Vec3 > edge_vertex_positions;
	std::unordered_map< FaceCRef, Vec3 > face_vertex_positions;

	//A2G3: Catmull-Clark Subdivision

	// This routine should end up looking a lot like linear_subdivide
	// above, with the exception that the positions are a bit trickier
	// to compute.

	//Overview of the rules:
	// https://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface

	// Faces
	for (FaceRef f = faces.begin(); f != faces.end(); ++f) {
        if (!f->boundary) {
            face_vertex_positions[f] = f->center();
        }
    }

	// Edges
	for (EdgeRef e = edges.begin(); e != edges.end(); ++e) {
        Vec3 midpoint = e->center();
        if (!e->on_boundary()) {
            Vec3 face_avg = (face_vertex_positions[e->halfedge->face] + face_vertex_positions[e->halfedge->twin->face]) * 0.5f;
            edge_vertex_positions[e] = (midpoint + face_avg) * 0.5f;
        } else {
            edge_vertex_positions[e] = midpoint;
        }
    }

	// Vertices
	for (VertexRef v = vertices.begin(); v != vertices.end(); ++v) {
        if (v->on_boundary()) {
			// boundary case: weighting of 3/4 on the original vertex and 1/8 on the neighboring boundary vertices.
            Vec3 boundary_sum(0.0f, 0.0f, 0.0f);
            float boundary_count = 0.0f;
            HalfedgeRef h = v->halfedge;
            do {
                if (h->edge->on_boundary()) {
                    boundary_sum += h->twin->vertex->position;
                    boundary_count++;
                }
                h = h->twin->next;
            } while (h != v->halfedge);
            vertex_positions[v] = (v->position * 0.75f) + (boundary_sum * (0.25f / (boundary_count)));
        } else {
            Vec3 Q(0.0f, 0.0f, 0.0f); // face center average
            Vec3 R(0.0f, 0.0f, 0.0f); // edge center average
            float n = 0.0f;

            HalfedgeRef h = v->halfedge;
            do {
                Q += h->face->center();
                R += h->edge->center();
                n++;
                h = h->twin->next;
            } while (h != v->halfedge);

            Q /= n;
            R /= n;
            vertex_positions[v] = (Q + (2.0f * R) + (n - 3.0f) * v->position) / n;
        }
    }
	
	//Now, use the provided helper function to actually perform the subdivision:
	catmark_subdivide_helper(vertex_positions, edge_vertex_positions, face_vertex_positions);

}

/*
 * loop_subdivide: sub-divide non-boundary faces with the Loop subdivision rule
 * 
 * If all non-boundary faces are triangles:
 *   subdivides mesh using the Loop subdivision rule
 *   returns true
 * Otherwise:
 *   does not change mesh
 *   returns false
 *
 * Do note that this requires a working implementation of edge split and edge flip
 */
bool Halfedge_Mesh::loop_subdivide() {

	//preamble: check for any non-triangular non-boundary faces:
	for (FaceCRef f = faces.begin(); f != faces.end(); ++f) {
		if (f->boundary) continue; //ignore boundary faces for this check
		if (f->halfedge->next->next->next != f->halfedge) {
			//found a non-triangular face!
			return false;
		}
	}

	//if execution reaches this point, all non-boundary faces are triangular, so proceed to subdivide:

	// A2Go1: Loop subdivision.
	// Optional! Only one of {A2Go1, A2Go2, A2Go3} is required!

	// Each vertex and edge of the original mesh can be associated with a
	// vertex in the new (subdivided) mesh. Therefore, our strategy for
	// computing the subdivided vertex locations is to *first* compute the
	// new positions using the connectivity of the original (coarse) mesh.
	// Navigating this mesh will be much easier than navigating the new
	// subdivided (fine) mesh, which has more elements to traverse.  We
	// will then assign vertex positions in the new mesh based on the
	// values we computed for the original mesh.
    
	// Compute new positions for all the vertices in the input mesh using
	// the Loop subdivision rule and store them in vertex_new_pos.
	[[maybe_unused]]
	std::unordered_map< VertexRef, Vec3 > vertex_new_pos;
	for (VertexRef v = vertices.begin(); v != vertices.end(); ++v) {
		if (v->on_boundary()) {
			// boundary case: weighting of 3/4 on the original vertex and 1/8 on the neighboring boundary vertices.
            Vec3 boundary_sum(0.0f, 0.0f, 0.0f);
            float boundary_count = 0.0f;
            HalfedgeRef h = v->halfedge;
            do {
                if (h->edge->on_boundary()) {
                    boundary_sum += h->twin->vertex->position;
                    boundary_count++;
                }
                h = h->twin->next;
            } while (h != v->halfedge);
            vertex_new_pos[v] = (v->position * 0.75f) + (boundary_sum * (0.25f / (boundary_count)));
        } else {
			int n = 0;
			Vec3 neighbor_sum(0.0f, 0.0f, 0.0f);
			HalfedgeRef h = v->halfedge;
			do {
				neighbor_sum += h->twin->vertex->position;
				n++;
				h = h->twin->next;
			} while (h != v->halfedge);
			
			float u = (n == 3) ? (3.0f / 16.0f) : (3.0f / (8.0f * n));
			vertex_new_pos[v] = (1 - n * u) * v->position + u * neighbor_sum;
		}
    }

	// std::cout << "vertex_new_pos map contents:\n";
    // for (const auto& pair : vertex_new_pos) {
    //     VertexRef v = pair.first;  // The vertex reference
    //     Vec3 pos = pair.second;    // The new position associated with the vertex
    //     std::cout << "Vertex ID: " << v->id << " -> Position: (" 
    //               << pos.x << ", " << pos.y << ", " << pos.z << ")\n";
    // }

	    
	// Next, compute the subdivided vertex positions associated with edges, and
	// store them in edge_new_pos:
	[[maybe_unused]]
	std::unordered_map< EdgeRef, Vec3 > edge_new_pos;
	for (EdgeRef e = edges.begin(); e != edges.end(); ++e) {
		if (e->on_boundary()) {
			VertexRef A = e->halfedge->vertex;
			VertexRef B = e->halfedge->twin->vertex;
			edge_new_pos[e] = (0.5f) * (A->position + B->position);
		} else {
			VertexRef A = e->halfedge->vertex;
			VertexRef B = e->halfedge->twin->vertex;
			VertexRef C = e->halfedge->next->next->vertex;
			VertexRef D = e->halfedge->twin->next->next->vertex;
			
			edge_new_pos[e] = (3.0f / 8.0f) * (A->position + B->position) + (1.0f / 8.0f) * (C->position + D->position);
		}
        
    }

	// std::cout << "edge_new_pos map contents:\n";
    // for (const auto& pair : edge_new_pos) {
    //     EdgeRef e = pair.first;  // The edge reference
    //     Vec3 pos = pair.second;  // The new position associated with the edge
    //     std::cout << "Edge ID: " << e->id << " -> Position: (" 
    //               << pos.x << ", " << pos.y << ", " << pos.z << ")\n";
    // }
    
	// Next, we're going to split every edge in the mesh, in any order, placing
	// the split vertex at the recorded edge_new_pos.
	//
	// We'll later need to distinguish edges that align with old edges to new
	// edges added by splitting. So store references to the new edges:
	[[maybe_unused]]
	std::vector< EdgeRef > new_edges;
	size_t num_edges = edges.size();
    EdgeRef E = edges.begin();
	// std::cout << std::endl << "edges: ";
	// for (EdgeRef edge = edges.begin(); edge != edges.end(); ++edge) {
	// 	VertexRef A = edge->halfedge->vertex;
	// 	VertexRef B = edge->halfedge->twin->vertex;
	// 	std::cout << "(" << A->position << ", " << B->position << "), ";
	// }
    for (size_t i = 0; i < num_edges; i++) {
        EdgeRef next_edge = E;
        ++next_edge;
        VertexRef old_v1 = E->halfedge->vertex;
        VertexRef old_v2 = E->halfedge->twin->vertex;
        auto v_opt = split_edge(E);
        if (v_opt.has_value()) {
			VertexRef v_new = v_opt.value();
            v_new->position = edge_new_pos[E];
			// HalfedgeRef last_elem = --halfedges.end();
            // new_edges.push_back(last_elem->edge);
			HalfedgeRef h = v_new->halfedge;
			do {
				EdgeRef new_edge = h->edge;
				VertexRef u = E->halfedge->vertex;
				VertexRef v = E->halfedge->twin->vertex;
				if (new_edge->halfedge->vertex != old_v1 && new_edge->halfedge->twin->vertex != old_v1 &&
                    new_edge->halfedge->vertex != old_v2 && new_edge->halfedge->twin->vertex != old_v2) {
                    new_edges.push_back(new_edge);
                }
				h = h->twin->next;
			} while (h != v_new->halfedge);
		}
        E = next_edge;
		// std::cout << std::endl << "edges: ";
		// for (EdgeRef edge = edges.begin(); edge != edges.end(); ++edge) {
		// 	VertexRef A = edge->halfedge->vertex;
		// 	VertexRef B = edge->halfedge->twin->vertex;
		// 	std::cout << "(" << A->position << ", " << B->position << "), ";
		// }
    }

	// std::cout << std::endl << "new edges: ";
	// for (EdgeRef edge : new_edges) {
	// 	VertexRef A = edge->halfedge->vertex;
	// 	VertexRef B = edge->halfedge->twin->vertex;
    //     std::cout << "(" << A->position << ", " << B->position << "), ";
    // }

	// Also note that in this loop, we only want to iterate over edges of the
	// original mesh. Otherwise, we'll end up splitting edges that we just split
	// (and the loop will never end!)

	// Now flip any new edge that connects a new and old vertex.
	// To check if a vertex is new, you can use a simple helper that
	// checks if has an entry in vertex_new_pos:
	[[maybe_unused]]
	auto is_new = [&vertex_new_pos](VertexRef v) -> bool {
		return !vertex_new_pos.count(v);
	};

    // Now flip any new edge that connects an old and new vertex.
    for (EdgeRef e : new_edges) {
        if (is_new(e->halfedge->vertex) != is_new(e->halfedge->twin->vertex)) {
            flip_edge(e);
			VertexRef A = e->halfedge->vertex;
			VertexRef B = e->halfedge->twin->vertex;
			// std::cout << "flipped: (" << A->position << ", " << B->position << "), " << std::endl;
        }
    }
    // Finally, copy new vertex positions into the Vertex::position.
	for (VertexRef v = vertices.begin(); v != vertices.end(); ++v) {
        if (vertex_new_pos.count(v)) {
            v->position = vertex_new_pos[v];
        }
    }

	// std::cout << std::endl << "vertices (" << vertices.size() << "): ";
	// for (VertexRef v = vertices.begin(); v != vertices.end(); ++v) {
    //     std::cout << v->position << std::endl;
    // }

	// std::cout << std::endl << "edges: ";
	// for (EdgeRef edge = edges.begin(); edge != edges.end(); ++edge) {
	// 	VertexRef A = edge->halfedge->vertex;
	// 	VertexRef B = edge->halfedge->twin->vertex;
    //     std::cout << "(" << A->position << ", " << B->position << "), ";
    // }

	

	return true;
}

//isotropic_remesh: improves mesh quality through local operations.
// Do note that this requires a working implementation of EdgeSplit, EdgeFlip, and EdgeCollapse
void Halfedge_Mesh::isotropic_remesh(Isotropic_Remesh_Parameters const &params) {

	//A2Go2: Isotropic Remeshing
	// Optional! Only one of {A2Go1, A2Go2, A2Go3} is required!

	// Compute the mean edge length. This will be the "target length".

    // Repeat the four main steps for `outer_iterations` iterations:

    // -> Split edges much longer than the target length.
	//     ("much longer" means > target length * params.longer_factor)

    // -> Collapse edges much shorter than the target length.
	//     ("much shorter" means < target length * params.shorter_factor)

    // -> Flip each edge if it improves vertex degree.

    // -> Finally, apply some tangential smoothing to the vertex positions.
	//     This means move every vertex in the plane of its normal,
	//     toward the centroid of its neighbors, by params.smoothing_step of
	//     the total distance (so, smoothing_step of 1 would move all the way,
	//     smoothing_step of 0 would not move).
	// -> Repeat the tangential smoothing part params.smoothing_iterations times.

	//NOTE: many of the steps in this function will be modifying the element
	//      lists they are looping over. Take care to avoid use-after-free
	//      or infinite-loop problems.

}

struct Edge_Record {
	Edge_Record() {
	}
	Edge_Record(std::unordered_map<uint32_t, Mat4>& VQ, Halfedge_Mesh::EdgeRef e) : edge(e) {
		
		// Compute the combined quadric from the edge endpoints.
        // -> Build the 3x3 linear system whose solution minimizes the quadric error
        //    associated with these two endpoints.
        // -> Use this system to solve for the optimal position, and store it in
        //    Edge_Record::optimal.
        // -> Also store the cost associated with collapsing this edge in
        //    Edge_Record::cost.
	}
	Halfedge_Mesh::EdgeRef edge;
	Vec3 optimal;
	float score;
};

bool operator<(const Edge_Record& r1, const Edge_Record& r2) {
	if (r1.score != r2.score) {
		return (r1.score < r2.score);
	}
	Halfedge_Mesh::EdgeRef e1 = r1.edge;
	Halfedge_Mesh::EdgeRef e2 = r2.edge;
	return &*e1 < &*e2;
}

template<class T> struct MutablePriorityQueue {
	void insert(const T& item) {
		queue.insert(item);
	}
	void remove(const T& item) {
		if (queue.find(item) != queue.end()) {
			queue.erase(item);
		}
	}
	const T& top() const {
		return *(queue.begin());
	}
	void pop() {
		queue.erase(queue.begin());
	}
	size_t size() {
		return queue.size();
	}

	std::set<T> queue;
};

/*
 * simplify: reduce edge count through collapses
 *  ratio: proportion of original faces to retain
 *
 * you may choose to have your implementation work only on triangle meshes,
 *  in which case it may return 'false' if there are non-triangular
 *  non-boundary faces
 *
 * returns false if it ran out of edges to collapse
 * returns true otherwise
 * 
 * Do note that this requires a working implementation of EdgeCollapse
 */
bool Halfedge_Mesh::simplify(float ratio) {

	//A2Go3: simplification
	// Optional! Only one of {A2Go1, A2Go2, A2Go3} is required!

	std::unordered_map<uint32_t, Mat4> face_quadrics;
	std::unordered_map<uint32_t, Mat4> vertex_quadrics;
	std::unordered_map<uint32_t, Edge_Record> edge_records;
	MutablePriorityQueue<Edge_Record> queue;

	// Compute initial quadrics for each face by writing the plane equation for
    // the face in homogeneous coordinates. These quadrics should be stored in
    // face_quadrics
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in vertex_quadrics
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an Edge_Record for each edge and sticking it in the
    //    queue. You may want to use the above MutablePriorityQueue<Edge_Record> for this.
    // -> Until reaching the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

    return false;
}

/*
 * catmark_subdivide_helper: add vertex in every edge and non-boundary face, set positions from parameters
 *
 * Works on all valid meshes.
 */
void Halfedge_Mesh::catmark_subdivide_helper(
	std::unordered_map< VertexCRef, Vec3 > const &vertex_positions, //positions for vertices after subdivision
	std::unordered_map< EdgeCRef, Vec3 > const &edge_vertex_positions, //positions for new vertices added in each edge
	std::unordered_map< FaceCRef, Vec3 > const &face_vertex_positions //positions for new vertices added in each face
	) {

	//check that positions were supplied for every vertex:
	for (VertexCRef v = vertices.begin(); v != vertices.end(); ++v) {
		if (!vertex_positions.count(v)) {
			throw std::runtime_error("No vertex position supplied for vertex with id " + std::to_string(v->id) + ".");
		}
	}

	//check that positions were supplied for every edge:
	for (EdgeCRef e = edges.begin(); e != edges.end(); ++e) {
		if (!edge_vertex_positions.count(e)) {
			throw std::runtime_error("No edge vertex position supplied for edge with id " + std::to_string(e->id) + ".");
		}
	}

	//check that positions were supplied for every (non-boundary) face:
	for (FaceCRef f = faces.begin(); f != faces.end(); ++f) {
		if (f->boundary) {
			if (face_vertex_positions.count(f)) {
				throw std::runtime_error("Extraneous vertex position was supplied for boundary face with id " + std::to_string(f->id) + ".");
			}
		} else {
			if (!face_vertex_positions.count(f)) {
				throw std::runtime_error("No vertex position supplied for face with id " + std::to_string(f->id) + ".");
			}
		}
	}


	{ //check that mesh is in a valid state to start with:
		auto error = validate();
		if (error) {
			throw std::runtime_error("catmark_subdivide_helper called on invalid mesh: " + error.value().second);
		}
	}

	if (vertices.empty() || edges.empty() || faces.empty()) {
		//empty mesh must be empty:
		assert(vertices.empty() && edges.empty() && faces.empty());
		return;
	}

	//store the old last vertex, face, and edge to allow iterating over only the old elements later:
	//(this works because the emplace_* functions add to the end of the element lists)
	VertexRef last_old_vertex = std::prev(vertices.end());
	EdgeRef last_old_edge = std::prev(edges.end());
	FaceRef last_old_face = std::prev(faces.end());
	
	//(can't store .end() iterators because emplace_back puts things "before the end")

	//------------------------
	//split every edge:
	//old halfedges stay connected to their vertices
	//old edge stays connected to e->halfedge->vertex

	//before:
	//     -----h---->
	//  v1 -----e----- v2
	//     <----t-----
	//after:
	//     --h->    --h2->
	//  v1 --e-- vm --e2-- v2
	//     <-t2-    <--t--

	for (EdgeRef e = edges.begin(); e != std::next(last_old_edge); ++e) {
		HalfedgeRef h = e->halfedge;
		HalfedgeRef t = h->twin; assert(t->edge == e);
		VertexRef v1 = h->vertex;
		VertexRef v2 = t->vertex;

		//new elements:
		VertexRef vm = emplace_vertex();
		HalfedgeRef h2 = emplace_halfedge();
		HalfedgeRef t2 = emplace_halfedge();
		EdgeRef e2 = emplace_edge(e->sharp);

		//middle vertex:
		vm->halfedge = h2; //could also use t2
		vm->position = edge_vertex_positions.at(e);
		interpolate_data({v1, v2}, vm);

		//second edge:
		e2->halfedge = h2;

		//second halfedge:
		h2->next = h->next;
		h2->twin = t;
		h2->vertex = vm;
		h2->edge = e2;
		h2->face = h->face;
		interpolate_data({h, h->next}, h2);

		//second twin halfedge:
		t2->next = t->next;
		t2->twin = h;
		t2->vertex = vm;
		t2->edge = e;
		t2->face = t->face;
		interpolate_data({t, t->next}, t2);

		//fix up pointers for existing halfedges:
		h->next = h2;
		h->twin = t2;

		t->next = t2;
		t->twin = h2;
		t->edge = e2;
	}


	//---------------------------
	//split (non-boundary) faces:

	//before:
	//
	//  v0 <-h7- v7 <-h6- v6
	//  |                 ^
	//  h0                h5
	//  v                 |
	//  v1       f        v5
	//  |                 ^
	//  h1                h4
	//  v                 |
	//  v2 -h2-> v3 -h3-> v4
	//
	//after:
	//  v0 <-h7- v7 <-h6- v6
	//  |        |        ^
	//  h0   f   e3   f3  h5
	//  v  --c-> |        |
	//  v1 --e0- vm --e2- v5 
	//  |  <-t-- |        ^
	//  h1   f1  e1  f2   h4
	//  v        |        |
	//  v2 -h2-> v3 -h3-> v4
	//
	// (each new eN has new halfedges as you'd expect,
	//  with eN->halfedge being directed toward the central vertex.)

	for (FaceRef f = faces.begin(); f != std::next(last_old_face); ++f) {
		if (f->boundary) continue; //skip boundary faces

		//get face halfedges:
		std::vector< HalfedgeRef > face_halfedges;
		{
			HalfedgeRef h = f->halfedge;
			do {
				face_halfedges.emplace_back(h);
				h = h->next;
			} while (h != f->halfedge);
			assert(face_halfedges.size() % 2 == 0); //should always be pairs of halfedges along subdivided edges!
		}

		//get face vertices and corners to interpolate data from:
		// (skip the odd vertices/halfedges -- they were just added)
		std::vector< HalfedgeCRef > from_corners;
		std::vector< VertexCRef > from_vertices;
		for (uint32_t i = 0; i < face_halfedges.size(); i += 2) {
			from_corners.emplace_back(face_halfedges[i]);
			from_vertices.emplace_back(face_halfedges[i]->vertex);
		}

		//add central vertex:
		VertexRef vm = emplace_vertex();
		vm->position = face_vertex_positions.at(f);
		interpolate_data(from_vertices, vm);

		//add halfedges and edges around the central vertex:
		std::vector< EdgeRef > inner_edges;
		for (uint32_t i = 0; i + 1 < face_halfedges.size(); i += 2) {
			EdgeRef e = emplace_edge(false);
			HalfedgeRef c = emplace_halfedge();
			HalfedgeRef t = emplace_halfedge();

			e->halfedge = c;

			//halfedge coming from the side:
			c->twin = t;
			//c->next will be set later
			c->vertex = face_halfedges[i+1]->vertex;
			c->edge = e;
			//c->face will be set later
			interpolate_data({face_halfedges[i+1]}, c); //just copy the data

			//halfedge coming from the center:
			t->twin = c;
			//t->next will be set later
			t->vertex = vm;
			t->edge = e;
			//t->face will be set later
			interpolate_data(from_corners, t);

			if (i == 0) vm->halfedge = t;

			//save edge for later connection:
			inner_edges.emplace_back(e);
		}

		//hook up pointers for all the quads:
		for (uint32_t i = 0; i + 1 < face_halfedges.size(); i += 2) {
			HalfedgeRef h0 = face_halfedges[i];
			HalfedgeRef h1 = inner_edges.at(i/2)->halfedge;
			HalfedgeRef h2 = inner_edges.at((i/2 == 0 ? inner_edges.size()-1 : i/2-1))->halfedge->twin;
			HalfedgeRef h3 = face_halfedges[(i == 0 ? face_halfedges.size()-1 : i-1)];

			//connect halfedges around the face:
			h0->next = h1;
			h1->next = h2;
			h2->next = h3;
			assert(h3->next == h0); //already connected and part of the face

			//connect halfedges to the face:
			if (i == 0) {
				//first face re-uses f:
				assert(f->halfedge == h0);
				assert(h0->face == f);
				h1->face = f;
				h2->face = f;
				assert(h3->face == f);
			} else {
				//other faces made fresh:
				FaceRef n = emplace_face(false);
				n->halfedge = h0;
				h0->face = n;
				h1->face = n;
				h2->face = n;
				h3->face = n;
			}
		}

	}

	//--------------------------
	//update positions for vertices
	for (VertexRef v = vertices.begin(); v != std::next(last_old_vertex); ++v) {
		v->position = vertex_positions.at(v);
	}

	{ //PARANOIA: sanity check:
		auto ret = validate();
		if (ret) {
			warn("After subdivide, validate says:\n  %s", ret.value().second.c_str());
		}
		assert(!ret && "subdivide helper should never break topology");
	}
}

/*
 * flip_orientation: flip direction of all halfedges
 *
 * works on all valid meshes.
 */
void Halfedge_Mesh::flip_orientation() {

	//store new h->vertex and v->halfedge pointers:
	std::unordered_map<Halfedge const *, VertexRef> he_to_v;
	std::unordered_map<Vertex const *, HalfedgeRef> v_to_he;
	for (auto &he : halfedges) {
		he_to_v[&he] = he.twin->vertex;
	}
	for (auto &v : vertices) {
		v_to_he[&v] = v.halfedge->twin;
	}

	//reverse all face loops:
	for (auto &face : faces) {
		//read off halfedges around face:
		std::vector<HalfedgeRef> hs;
		std::vector<Vec2> uvs;
		std::vector<Vec3> normals;

		HalfedgeRef h = face.halfedge;
		do {
			hs.emplace_back(h);
			uvs.emplace_back(h->corner_uv);
			normals.emplace_back(h->corner_normal);
			h = h->next;
		} while (h != face.halfedge);

		//reverse face ordering:
		for (uint32_t i = 0; i < hs.size(); ++i) {
			hs[(i+1)%hs.size()]->next = hs[i];
			hs[i]->corner_uv = uvs[(i+1)%hs.size()];
			hs[i]->corner_normal = normals[(i+1)%hs.size()];
		}
	}

	//update h->vertex and v->halfedge pointers:
	for (auto &he : halfedges) {
		he.vertex = he_to_v.at(&he);
	}
	for (auto &v : vertices) {
		v.halfedge = v_to_he.at(&v);
	}
}

/*
 * set_corner_normals: compute face-corner normals based on `sharp` flag and smoothing threshold
 *
 * works on all valid meshes.
 */
void Halfedge_Mesh::set_corner_normals(float threshold) {
	//first, figure out which edges to consider sharp for this operation:
	std::unordered_set< Edge const * > sharp_edges;
	sharp_edges.reserve(edges.size());

	//all edges between boundary and non-boundary get marked sharp regardless of mode:
	for (auto const &edge : edges) {
		if (edge.halfedge->face->boundary != edge.halfedge->twin->face->boundary) {
			sharp_edges.emplace(&edge);
		}
	}

	if (threshold >= 180.0f) {
		//"smooth mode" -- all other edges are considered smooth
	} else {
		//"flat mode" / "auto mode" -- any edges which are marked sharp or have face angle <= threshold get marked sharp:
		float cos_threshold = std::cos( Radians( std::clamp(threshold, 0.0f, 180.0f) ) );
		if (threshold <= 0.0f) cos_threshold = 2.0f; //make sure everything is sharp
		for (auto const &edge : edges) {
			//get adjacent halfedges:
			HalfedgeRef h1 = edge.halfedge;
			HalfedgeRef h2 = h1->twin;
			if (h1->face->boundary || h2->face->boundary) {
				//don't care about edges boundary-boundary, and inside-boundary already marked.
				//thus: nothing to do here
			} else if (edge.sharp) {
				//flagged as sharp, so mark it sharp:
				sharp_edges.emplace(&edge);
			} else {
				//inside-inside edge, non-marked, check angle:
				Vec3 n1 = h1->face->normal();
				Vec3 n2 = h2->face->normal();
				float cos = dot(n1,n2);
				if (cos <= cos_threshold) {
					//treat as sharp:
					sharp_edges.emplace(&edge);
				}
			}
		}
	}

	//clear current corner normals:
	for (auto h = halfedges.begin(); h != halfedges.end(); ++h) {
		h->corner_normal = Vec3{0.0f, 0.0f, 0.0f};
	}

	//now circulate all vertices to set normals:
	for (auto const &v : vertices) {
		//get halfedge leaving this vertex:
		HalfedgeRef begin = v.halfedge;
		assert(&*begin->vertex == &v);

		//circulate begin until it is at a sharp edge (thus, the next corner starts a smoothing group):
		do {
			if (sharp_edges.count(&*begin->edge)) break;
			begin = begin->twin->next;
		} while (begin != v.halfedge); //could be all one big happy smoothing group

		//store all corners around the vertex:
		struct Corner {
			HalfedgeRef in; //halfedge pointing to v
			HalfedgeRef out; //halfedge pointing away from v
			Vec3 weighted_normal; //face normal at corner, weighted... somehow (see below)
		};

		std::vector< std::vector< Corner > > groups;
		HalfedgeRef h = begin;
		do {
			//start a new smoothing group on sharp edges (or at the very first edge):
			if (h == begin || sharp_edges.count(&*h->edge)) {
				groups.emplace_back();
			}
			//add corner after h to current smoothing group:
			Corner corner;
			corner.in = h->twin;
			corner.out = h->twin->next;
			assert(corner.in->face == corner.out->face); //PARANOIA
			{ //compute some sort of weighted normal:
				assert(&*corner.in->vertex != &v);
				assert(&*corner.in->twin->vertex == &v);
				Vec3 from = corner.in->vertex->position - v.position;

				assert(&*corner.out->vertex == &v);
				assert(&*corner.out->twin->vertex != &v);
				Vec3 to = corner.out->twin->vertex->position - v.position;
				/*
				//basic area weighting (weird with non-flat faces and reflex vertices):
				corner.weighted_normal = cross(to - v.position, from - v.position);
				*/
				/*//sort sort of angle weighting thing -- this never works as well as one would hope:
				//...still needs work for reflex angles also
				float angle = std::atan2(cross(from,to).norm(), dot(from, to));
				corner.weighted_normal = angle * corner.in->face->normal();
				*/
				//some other sort of slightly fancy area weighting:
				corner.weighted_normal = cross(to - v.position, from - v.position).norm() * corner.in->face->normal();
			}
			groups.back().emplace_back(corner);

			//advance h:
			h = h->twin->next;
		} while (h != begin);

		//compute weighted normals per-corner:
		for (auto const &group : groups) {
			assert(!group.empty());
			if (group[0].in->face->boundary) {
				//boundary group.
				//PARANOIA:
				for (auto const &corner : group) {
					assert(corner.in->face->boundary);
				}
				//no need for normals on boundary corners
				continue;
			}
			//compute weighted normal:
			Vec3 sum = Vec3{0.0f, 0.0f, 0.0f};
			for (auto const &corner : group) {
				sum += corner.weighted_normal;
			}
			//normalize:
			sum = sum.unit();
			//assign to all corners in group:
			for (auto const &corner : group) {
				assert(&*corner.out->vertex == &v);
				corner.out->corner_normal = sum;
			}
		}
	}

	//normals computed!
}

/*
 * set_corner_uvs_per_face: set uv coordinates to map texture per-face
 */
void Halfedge_Mesh::set_corner_uvs_per_face() {
	//clear existing UVs:
	for (auto &halfedge : halfedges) {
		halfedge.corner_uv = Vec2(0.0f, 0.0f);
	}
	
	//set UVs per-face:
	for (auto const &face : faces) {
		if (face.boundary) continue;

		//come up with a plane perpendicular-ish to the face:
		Vec3 n = face.normal();
		Vec3 p1;
		if (std::abs(n.x) < std::abs(n.y) && std::abs(n.x) < std::abs(n.z)) {
			p1 = Vec3(1.0f, 0.0f, 0.0f);
		} else if (std::abs(n.y) < std::abs(n.z)) {
			p1 = Vec3(0.0f, 1.0f, 0.0f);
		} else {
			p1 = Vec3(0.0f, 0.0f, 1.0f);
		}
		p1 = (p1 - dot(p1, n) * n).unit();
		Vec3 p2 = cross(n, p1);

		//find bounds of face on plane:
		Vec2 min = Vec2(std::numeric_limits< float >::infinity(), std::numeric_limits< float >::infinity());
		Vec2 max = Vec2(-std::numeric_limits< float >::infinity(), -std::numeric_limits< float >::infinity());
		HalfedgeRef v = face.halfedge;
		do {
			Vec2 pt = Vec2(dot(p1, v->vertex->position), dot(p2, v->vertex->position));
			min = hmin(min, pt);
			max = hmax(max, pt);
			v = v->next;
		} while (v != face.halfedge);

		//set corner uvs based on position within bounds:
		do {
			Vec2 pt = Vec2(dot(p1, v->vertex->position), dot(p2, v->vertex->position));
			v->corner_uv = Vec2(
				(pt.x - min.x) / (max.x - min.x),
				(pt.y - min.y) / (max.y - min.y)
			);
			v = v->next;
		} while (v != face.halfedge);
	}
}

/*
 * set_corner_uvs_project: set uv coordinates to map texture by projection to a plane
 */
void Halfedge_Mesh::set_corner_uvs_project(Vec3 origin, Vec3 u_axis, Vec3 v_axis) {

	u_axis /= u_axis.norm_squared();
	v_axis /= v_axis.norm_squared();

	for (auto &halfedge : halfedges) {
		if (halfedge.face->boundary) {
			halfedge.corner_uv = Vec2(0.0f, 0.0f);
		} else {
			halfedge.corner_uv = Vec2(
				dot(halfedge.vertex->position - origin, u_axis),
				dot(halfedge.vertex->position - origin, v_axis)
			);
		}
	}
}
