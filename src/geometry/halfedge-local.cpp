
#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it cannot perform an operation (i.e., because
    the resulting mesh does not have a valid representation).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
 * add_face: add a standalone face to the mesh
 *  sides: number of sides
 *  radius: distance from vertices to origin
 *
 * We provide this method as an example of how to make new halfedge mesh geometry.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::add_face(uint32_t sides, float radius) {
	//faces with fewer than three sides are invalid, so abort the operation:
	if (sides < 3) return std::nullopt;


	std::vector< VertexRef > face_vertices;
	//In order to make the first edge point in the +x direction, first vertex should
	// be at -90.0f - 0.5f * 360.0f / float(sides) degrees, so:
	float const start_angle = (-0.25f - 0.5f / float(sides)) * 2.0f * PI_F;
	for (uint32_t s = 0; s < sides; ++s) {
		float angle = float(s) / float(sides) * 2.0f * PI_F + start_angle;
		VertexRef v = emplace_vertex();
		v->position = radius * Vec3(std::cos(angle), std::sin(angle), 0.0f);
		face_vertices.emplace_back(v);
	}

	assert(face_vertices.size() == sides);

	//assemble the rest of the mesh parts:
	FaceRef face = emplace_face(false); //the face to return
	FaceRef boundary = emplace_face(true); //the boundary loop around the face

	std::vector< HalfedgeRef > face_halfedges; //will use later to set ->next pointers

	for (uint32_t s = 0; s < sides; ++s) {
		//will create elements for edge from a->b:
		VertexRef a = face_vertices[s];
		VertexRef b = face_vertices[(s+1)%sides];

		//h is the edge on face:
		HalfedgeRef h = emplace_halfedge();
		//t is the twin, lies on boundary:
		HalfedgeRef t = emplace_halfedge();
		//e is the edge corresponding to h,t:
		EdgeRef e = emplace_edge(false); //false: non-sharp

		//set element data to something reasonable:
		//(most ops will do this with interpolate_data(), but no data to interpolate here)
		h->corner_uv = a->position.xy() / (2.0f * radius) + 0.5f;
		h->corner_normal = Vec3(0.0f, 0.0f, 1.0f);
		t->corner_uv = b->position.xy() / (2.0f * radius) + 0.5f;
		t->corner_normal = Vec3(0.0f, 0.0f,-1.0f);

		//thing -> halfedge pointers:
		e->halfedge = h;
		a->halfedge = h;
		if (s == 0) face->halfedge = h;
		if (s + 1 == sides) boundary->halfedge = t;

		//halfedge -> thing pointers (except 'next' -- will set that later)
		h->twin = t;
		h->vertex = a;
		h->edge = e;
		h->face = face;

		t->twin = h;
		t->vertex = b;
		t->edge = e;
		t->face = boundary;

		face_halfedges.emplace_back(h);
	}

	assert(face_halfedges.size() == sides);

	for (uint32_t s = 0; s < sides; ++s) {
		face_halfedges[s]->next = face_halfedges[(s+1)%sides];
		face_halfedges[(s+1)%sides]->twin->next = face_halfedges[s]->twin;
	}

	return face;
}


/*
 * bisect_edge: split an edge without splitting the adjacent faces
 *  e: edge to split
 *
 * returns: added vertex
 *
 * We provide this as an example for how to implement local operations.
 * (and as a useful subroutine!)
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {
	// Phase 0: draw a picture
	//
	// before:
	//    ----h--->
	// v1 ----e--- v2
	//   <----t---
	//
	// after:
	//    --h->    --h2->
	// v1 --e-- vm --e2-- v2
	//    <-t2-    <--t--
	//

	// Phase 1: collect existing elements
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->vertex;
	VertexRef v2 = t->vertex;

	// Phase 2: Allocate new elements, set data
	VertexRef vm = emplace_vertex();
	vm->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, vm); //set bone_weights

	EdgeRef e2 = emplace_edge();
	e2->sharp = e->sharp; //copy sharpness flag

	HalfedgeRef h2 = emplace_halfedge();
	interpolate_data({h, h->next}, h2); //set corner_uv, corner_normal

	HalfedgeRef t2 = emplace_halfedge();
	interpolate_data({t, t->next}, t2); //set corner_uv, corner_normal

	// The following elements aren't necessary for the bisect_edge, but they are here to demonstrate phase 4
    FaceRef f_not_used = emplace_face();
    HalfedgeRef h_not_used = emplace_halfedge();

	// Phase 3: Reassign connectivity (careful about ordering so you don't overwrite values you may need later!)

	vm->halfedge = h2;

	e2->halfedge = h2;

	assert(e->halfedge == h); //unchanged

	//n.b. h remains on the same face so even if h->face->halfedge == h, no fixup needed (t, similarly)

	h2->twin = t;
	h2->next = h->next;
	h2->vertex = vm;
	h2->edge = e2;
	h2->face = h->face;

	t2->twin = h;
	t2->next = t->next;
	t2->vertex = vm;
	t2->edge = e;
	t2->face = t->face;
	
	h->twin = t2;
	h->next = h2;
	assert(h->vertex == v1); // unchanged
	assert(h->edge == e); // unchanged
	//h->face unchanged

	t->twin = h2;
	t->next = t2;
	assert(t->vertex == v2); // unchanged
	t->edge = e2;
	//t->face unchanged


	// Phase 4: Delete unused elements
    erase_face(f_not_used);
    erase_halfedge(h_not_used);

	// Phase 5: Return the correct iterator
	return vm;
}


/*
 * split_edge: split an edge and adjacent (non-boundary) faces
 *  e: edge to split
 *
 * returns: added vertex. vertex->halfedge should lie along e
 *
 * Note that when splitting the adjacent faces, the new edge
 * should connect to the vertex ccw from the ccw-most end of e
 * within the face.
 *
 * Do not split adjacent boundary faces.
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(EdgeRef e) {
	// A2L2 (REQUIRED): split_edge

    // collect
    HalfedgeRef h1 = e->halfedge;
    HalfedgeRef h2 = h1->twin;

	if (h1->face->boundary && h2->face->boundary) {
		return std::nullopt;
	}
    

    VertexRef v1 = h1->vertex;
    VertexRef v2 = h2->vertex;
    VertexRef v3 = h1->next->next->vertex;
    VertexRef v4 = h2->next->next->vertex;
    
    FaceRef f1 = h1->face;
    FaceRef f2 = h2->face;

    // create new stuff
    VertexRef v_new = *bisect_edge(e);

	HalfedgeRef h2_new = v_new->halfedge;
	HalfedgeRef h1_new = h2_new->twin;
	f1 = h1_new->face;
	f2 = h2_new->face; 
	h1 = h1_new->next;
	h2 = h1->twin;
	//EdgeRef e1 = h2_new->edge;

	v1 = h1_new->vertex;
    v2 = h2->vertex;
    v3 = h1->next->next->vertex;
    v4 = h2_new->next->next->vertex;

	if (!h1->face->boundary) {
		//std::cout << "boundary 1";
		EdgeRef e2 = emplace_edge();

		FaceRef f3 = emplace_face();

		HalfedgeRef h3 = emplace_halfedge();
    	HalfedgeRef h3_twin = emplace_halfedge();
		
		// Set halfedge connectivity

		

		h3->set_tnvef(h3_twin, h1->next->next, v_new, e2, f1);
    	h3_twin->set_tnvef(h3, h1, v3, e2, f3);

		h1->face = f3;
		h1_new->face = f1;
		h1_new->next = h3;
		h1->next->next = h3_twin;
		h1->next->face = f3;

		// set new edge connectivity

		e2->halfedge = h3;

		// set new face connectivity

		f1->halfedge = h1_new;
		f3->halfedge = h1;

		// std::cout << "f1_id: " << f1->id << ", f3_id: " << f3->id << std::endl;

		// std::cout << "h1_id: " << h1->id << ", h1_new_id: " << h1_new->id << "h3_id: " << h3->id << ", h3_twin_id: " << h3_twin->id << std::endl;

		HalfedgeRef h = f1->halfedge;
		// std::cout << "f1 ids = ";
		// do {
		// 	h = h->next;
		// 	std::cout << h->id << ", ";
		// }
		// while (h != f1->halfedge);
		// std::cout << std::endl;

		

	}

	if(!h2->face->boundary) {
		//std::cout << "boundary 2";
		EdgeRef e3 = emplace_edge();

		FaceRef f4 = emplace_face();

		
		HalfedgeRef h4 = emplace_halfedge();
		HalfedgeRef h4_twin = emplace_halfedge();

		// Set halfedge connectivity
		
		
		
		h4->set_tnvef(h4_twin, h2_new->next->next, v_new, e3, f2);
		h4_twin->set_tnvef(h4, h2_new, v4, e3, f4);

		h2->face = f2;
		h2_new->face = f4;
		h2->next = h4;
		h2_new->next->next = h4_twin;
		h2_new->next->face = f4;

		// set new edge connectivity
		
		e3->halfedge = h4;

		// set new face connectivity
		
		f2->halfedge = h2;
		f4->halfedge = h2_new;

		// std::cout << "f2_id: " << f2->id << ", f4_id: " << f4->id << std::endl;

		// std::cout << "h2_id: " << h2->id << ", h2_new_id: " << h2_new->id << "h4_id: " << h4->id << ", h4_twin_id: " << h4_twin->id << std::endl;


		// HalfedgeRef h = f2->halfedge;
		// std::cout << "f2 ids = ";
		// do {
		// 	h = h->next;
		// 	std::cout << h->id << ", ";
		// }
		// while (h != f2->halfedge);
		// std::cout << std::endl;
	}

    
    

    return v_new;
}



/*
 * inset_vertex: divide a face into triangles by placing a vertex at f->center()
 *  f: the face to add the vertex to
 *
 * returns:
 *  std::nullopt if insetting a vertex would make mesh invalid
 *  the inset vertex otherwise
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
	// A2Lx4 (OPTIONAL): inset vertex
	
	(void)f;
    return std::nullopt;
}


/* [BEVEL NOTE] Note on the beveling process:

	Each of the bevel_vertex, bevel_edge, and extrude_face functions do not represent
	a full bevel/extrude operation. Instead, they should update the _connectivity_ of
	the mesh, _not_ the positions of newly created vertices. In fact, you should set
	the positions of new vertices to be exactly the same as wherever they "started from."

	When you click on a mesh element while in bevel mode, one of those three functions
	is called. But, because you may then adjust the distance/offset of the newly
	beveled face, we need another method of updating the positions of the new vertices.

	This is where bevel_positions and extrude_positions come in: these functions are
	called repeatedly as you move your mouse, the position of which determines the
	amount / shrink parameters. These functions are also passed an array of the original
	vertex positions, stored just after the bevel/extrude call, in order starting at
	face->halfedge->vertex, and the original element normal, computed just *before* the
	bevel/extrude call.

	Finally, note that the amount, extrude, and/or shrink parameters are not relative
	values -- you should compute a particular new position from them, not a delta to
	apply.
*/

/*
 * bevel_vertex: creates a face in place of a vertex
 *  v: the vertex to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(VertexRef v) {
	//A2Lx5 (OPTIONAL): Bevel Vertex
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in bevel_vertex_helper (A2Lx5h)

	(void)v;
    return std::nullopt;
}

/*
 * bevel_edge: creates a face in place of an edge
 *  e: the edge to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(EdgeRef e) {
	//A2Lx6 (OPTIONAL): Bevel Edge
	// Reminder: This function does not update the vertex positions.
	// remember to also fill in bevel_edge_helper (A2Lx6h)

	(void)e;
    return std::nullopt;
}

/*
 * extrude_face: creates a face inset into a face
 *  f: the face to inset
 *
 * returns: reference to the inner face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::extrude_face(FaceRef f) {
	//A2L4: Extrude Face
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in extrude_helper (A2L4h)

	if (f->boundary) {
        return std::nullopt;
    }

    // collect
    std::vector<VertexRef> original_vertices;
    std::vector<HalfedgeRef> original_halfedges;
    HalfedgeRef h = f->halfedge;
    do {
        original_vertices.emplace_back(h->vertex);
        original_halfedges.emplace_back(h);
        h = h->next;
    } while (h != f->halfedge);

    size_t num_vertices = original_vertices.size();

	assert(num_vertices == f->degree());

	// Print IDs of original vertices
    std::cout << "Original Vertices IDs: ";
    for (const auto& v : original_vertices) {
        std::cout << v->id << " ";
    }
    std::cout << std::endl;

    // Print IDs of original halfedges
    std::cout << "Original Halfedges IDs: ";
    for (const auto& he : original_halfedges) {
        std::cout << he->id << " ";
    }
    std::cout << std::endl;
    
    // create
    std::vector<VertexRef> new_vertices(num_vertices);
    std::vector<EdgeRef> side_edges(num_vertices);
	std::vector<EdgeRef> face_edges(num_vertices);
    std::vector<FaceRef> side_faces(num_vertices);
    std::vector<HalfedgeRef> new_face_halfedges_inner(num_vertices);
	std::vector<HalfedgeRef> new_face_halfedges_outer(num_vertices);
	std::vector<HalfedgeRef> new_quads_halfedges(num_vertices);
	std::vector<HalfedgeRef> new_quads_halfedges_twins(num_vertices);

    for (size_t i = 0; i < num_vertices; i++) {
        new_vertices[i] = emplace_vertex();
        new_vertices[i]->position = original_vertices[i]->position;
        
        side_edges[i] = emplace_edge();
		face_edges[i] = emplace_edge();

        side_faces[i] = emplace_face();

        new_face_halfedges_inner[i] = emplace_halfedge();
		new_face_halfedges_outer[i] = emplace_halfedge();
		new_quads_halfedges[i] = emplace_halfedge();
		new_quads_halfedges_twins[i] = emplace_halfedge();
    }

    // FaceRef inner_face = emplace_face();
    
    // iterate, collect, and connect
    for (size_t i = 0; i < num_vertices; i++) {
		
        size_t next = (i + 1) % num_vertices;
		size_t prev = (i == 0) ? num_vertices - 1 : i - 1;

        
        // Retrieve original halfedges
        HalfedgeRef h_old = original_halfedges[i];
        
        // collect side quad halfedges
        HalfedgeRef h_new = new_quads_halfedges[i];
        HalfedgeRef h_new_twin = new_quads_halfedges_twins[i];

        
        // connect side quad halfedges
        h_new->set_tnvef(h_new_twin, h_old, new_vertices[i], side_edges[i], side_faces[i]);
        h_new_twin->set_tnvef(h_new, new_face_halfedges_outer[prev], original_vertices[i], side_edges[i], side_faces[prev]);
        
        
        side_faces[i]->halfedge = h_new;
        side_edges[i]->halfedge = h_new;
		
		new_vertices[i]->halfedge = h_new;

        h_old->next = new_quads_halfedges_twins[next];
        h_old->face = side_faces[i];

		// collect face halfedges

		HalfedgeRef face_halfedge = new_face_halfedges_outer[i];
        HalfedgeRef face_halfedge_twin = new_face_halfedges_inner[i];

		// connect

		face_halfedge->set_tnvef(face_halfedge_twin, h_new, new_vertices[next], face_edges[i], side_faces[i]);
        face_halfedge_twin->set_tnvef(face_halfedge, new_face_halfedges_inner[next], new_vertices[i], face_edges[i], f);

		face_edges[i]->halfedge = face_halfedge;
		f->halfedge = face_halfedge_twin;

		std::cout << "i = " << i << ", h_new = " << h_new->id << ", h_new_twin = " << h_new_twin->id << ", face halfend = " << face_halfedge->id << ", face halfedge twin = " << face_halfedge_twin->id << ", h_old = " << h_old->id << ", h_old next: " << h_old->next->id << std::endl;
		std::cout << "side face = " << side_faces[i]->id << ", side edge = " << side_edges[i]->id << ", new vertex = " << new_vertices[i]->id << ", old vertex = " << original_vertices[i]->id << ", old vertex halfedge= " << original_vertices[i]->halfedge->id<< ", face edge = " << face_edges[i]->id << std::endl;
    }

	std::cout << "inner face = " << f->id << std::endl;

    return f;
}

/*
 * flip_edge: rotate non-boundary edge ccw inside its containing faces
 *  e: edge to flip
 *
 * if e is a boundary edge, does nothing and returns std::nullopt
 * if flipping e would create an invalid mesh, does nothing and returns std::nullopt
 *
 * otherwise returns the edge, post-rotation
 *
 * does not create or destroy mesh elements.
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(EdgeRef e) {
	//A2L1: Flip Edge

	if (e->on_boundary()) {
		return std::nullopt;
	}

	// collect
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->next->vertex;
	VertexRef v2 = t->next->vertex;
	VertexRef v3 = h->next->next->vertex;
	VertexRef v4 = t->next->next->vertex;
	FaceRef f1 = h->face;
	FaceRef f2 = t->face;


	HalfedgeRef h_prev = h;
	do {
		if (h_prev->next == h) {
			break;
		}
		h_prev = h_prev->next;
	} while (h_prev != h);
	assert(h_prev != h);

	HalfedgeRef t_prev = t;
	do {
		if (t_prev->next == t) {
			break;
		}
		t_prev = t_prev->next;
	} while (t_prev != t);
	assert(t_prev != t);

	HalfedgeRef h_next = h->next;
	HalfedgeRef t_next = t->next;


	// disconnect
	h_prev->next = t_next;
	t_prev->next = h_next;
	
	v1->halfedge = h_next;
	v2->halfedge = t_next;
	f1->halfedge = h;
	f2->halfedge = t;

	t->next->face = f1;
	h->next->face = f2;

	


	// connect
	h->next = h_next->next;
	t->next = t_next->next;
	t_next->next = h;
	h_next->next = t;
	t->vertex = v3;
	h->vertex = v4;
	

	return std::make_optional(e);
}


/*
 * make_boundary: add non-boundary face to boundary
 *  face: the face to make part of the boundary
 *
 * if face ends up adjacent to other boundary faces, merge them into face
 *
 * if resulting mesh would be invalid, does nothing and returns std::nullopt
 * otherwise returns face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::make_boundary(FaceRef face) {
	//A2Lx7: (OPTIONAL) make_boundary

	return std::nullopt; //TODO: actually write this code!
}

/*
 * dissolve_vertex: merge non-boundary faces adjacent to vertex, removing vertex
 *  v: vertex to merge around
 *
 * if merging would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_vertex(VertexRef v) {
	// A2Lx1 (OPTIONAL): Dissolve Vertex

    return std::nullopt;
}

/*
 * dissolve_edge: merge the two faces on either side of an edge
 *  e: the edge to dissolve
 *
 * merging a boundary and non-boundary face produces a boundary face.
 *
 * if the result of the merge would be an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_edge(EdgeRef e) {
	// A2Lx2 (OPTIONAL): dissolve_edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data
	
    return std::nullopt;
}

/* collapse_edge: collapse edge to a vertex at its middle
 *  e: the edge to collapse
 *
 * if collapsing the edge would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(EdgeRef e) {
	//A2L3: Collapse Edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

	// if (e->on_boundary()) {
    //     return std::nullopt;
    // }

    // collect
    HalfedgeRef h1 = e->halfedge;
    HalfedgeRef h2 = h1->twin;
	HalfedgeRef h_new = h1->next;
    
    VertexRef v1 = h1->vertex;
    VertexRef v2 = h2->vertex;

	FaceRef f1 = h1->face;
	FaceRef f2 = h2->face;
	uint32_t f1_degree = f1->degree();
	uint32_t f2_degree = f2->degree();

	std::cout << "v1_id: " << v1->id << ", v2_id: " << v2->id << std::endl;

	std::cout << "h1_id: " << h1->id << ", h2_id: " << h2->id << std::endl;

	std::cout << "f1_id: " << f1->id << " ,f1_degree: " << f1->degree() << ", f2_id: " << f2->id << " ,f2_degree: " << f2->degree() << std::endl;

    if (f1->degree() <= 2 || f2->degree() <= 2) {
        return std::nullopt;
    }

	if (f1->degree() == 3) {
		bool all_boundary = true;
		HalfedgeRef h = f1->halfedge;
		do {
			if (!h->edge->on_boundary()) {
				all_boundary = false; 
			}
			h = h->next;
		} while (h != f1->halfedge);
		if (all_boundary) {
			return std::nullopt;
		}
    }

	if (f2->degree() == 3) {
		bool all_boundary = true;
		HalfedgeRef h = f2->halfedge;
		do {
			if (!h->edge->on_boundary()) {
				all_boundary = false; 
			}
			h = h->next;
		} while (h != f2->halfedge);
		if (all_boundary) {
			return std::nullopt;
		}
    }

    // Create
    VertexRef v_new = emplace_vertex();
    v_new->position = (v1->position + v2->position) * 0.5f;
    interpolate_data({v1, v2}, v_new);

	std::cout << "v_new_id: " << v_new->id << std::endl;


	// helper for triangle faces
	auto merge_edges = [&](FaceRef f) {
        //if (f->degree() == 3) {
			std::cout << "merging " << f->id << " with degree " << f->degree() << std::endl;

			// collect
			HalfedgeRef h = e->halfedge;
			if (h->face != f) {
				h = h->twin;
			}
            HalfedgeRef h_next = h->next;
			HalfedgeRef h_next_twin = h_next->twin;
            HalfedgeRef h_prev = h_next->next;
			HalfedgeRef h_prev_twin = h_prev->twin;
			std::cout << "h_next: " << h_next->id << ", h_prev: " << h_prev->id << std::endl;
			std::cout << "h_next_twin: " << h_next_twin->id << ", h_prev_twin: " << h_prev_twin->id << std::endl;

			VertexRef v3 = h_prev->vertex;
            
            EdgeRef e_merge = h_next->edge;
            EdgeRef e_remove = h_prev->edge;

			// REMOVE: h_next, h_prev
			// KEEP: h_next_twin, h_prev_twin
            
            // interpolate
            interpolate_data({h_next_twin, h_prev}, h_next_twin);
            interpolate_data({h_next, h_prev_twin}, h_prev_twin);
            
            // connect

			// vertices new halfedge
			v_new->halfedge =h_prev_twin;
			v3->halfedge = h_next_twin;

			
			e_merge->halfedge = h_next_twin;
			
			// h prev twin new values
			h_prev_twin->edge = e_merge;
			h_prev_twin->vertex = v_new;

			// set new twins
			h_next_twin->twin = h_prev_twin;
			h_prev_twin->twin = h_next_twin;
            
            // delete
            erase_halfedge(h_next);
            erase_halfedge(h_prev);
			erase_edge(e_remove);
            erase_face(f);

			std::cout << "finish merge: " << std::endl;
        //}
    };

	
    

    // connect halfedges to new vertex
	std::cout << "connecting v1" << std::endl;
	std::cout << "ids: ";
    HalfedgeRef h = v1->halfedge;
    do {
		std::cout << h->id << ", ";
        h->vertex = v_new;
        h = h->twin->next;
    } while (h != v1->halfedge);

	std::cout << "\nconnecting v2" << std::endl;
    h = v2->halfedge;
    do {
        h->vertex = v_new;
        h = h->twin->next;
    } while (h != v2->halfedge);


	// if triangle, merge edges
	if (f1_degree == 3) {
		merge_edges(f1);
	} else {
		f1->halfedge = h1->next;
		std::cout << "connecting f1" << std::endl;
		std::cout << "ids: ";
		h = h1;
		do {
			std::cout << h->id << ", ";
			if (h->next == h1) {
				h->next = h->next->next;
				break;
			}
			h = h->next;
		} while (h != h1);
	}
	if (f2_degree == 3) {
		merge_edges(f2);
	} else {
		f2->halfedge = h2->next;
		h = h2;
		do {
			if (h->next == h2) {
				h->next = h->next->next;
				break;
			}
			h = h->next;
		} while (h != h2);
	}
	std::cout << "v_new halfedge id: " << v_new->halfedge->id << std::endl; 
	// make sure halfedges poiting to old halfedges are redirected
	

	std::cout << "\nconnecting f2" << std::endl;
	

	// v_new halfedge
	if (f1_degree != 3 && f2_degree != 3) {
		v_new->halfedge = h_new;
	}

	

    // delete
    erase_edge(e);
	erase_halfedge(h1);
	erase_halfedge(h2);
    erase_vertex(v1);
    erase_vertex(v2);

	std::cout << "v_new halfedge id: " << v_new->halfedge->id << std::endl; 

    return std::make_optional(v_new);
	
    return std::nullopt;
}

/*
 * collapse_face: collapse a face to a single vertex at its center
 *  f: the face to collapse
 *
 * if collapsing the face would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(FaceRef f) {
	//A2Lx3 (OPTIONAL): Collapse Face

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

    return std::nullopt;
}

/*
 * weld_edges: glue two boundary edges together to make one non-boundary edge
 *  e, e2: the edges to weld
 *
 * if welding the edges would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns e, updated to represent the newly-welded edge
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::weld_edges(EdgeRef e, EdgeRef e2) {
	//A2Lx8: Weld Edges

	//Reminder: use interpolate_data() to merge bone_weights data on vertices!

    return std::nullopt;
}



/*
 * bevel_positions: compute new positions for the vertices of a beveled vertex/edge
 *  face: the face that was created by the bevel operation
 *  start_positions: the starting positions of the vertices
 *     start_positions[i] is the starting position of face->halfedge(->next)^i
 *  direction: direction to bevel in (unit vector)
 *  distance: how far to bevel
 *
 * push each vertex from its starting position along its outgoing edge until it has
 *  moved distance `distance` in direction `direction`. If it runs out of edge to
 *  move along, you may choose to extrapolate, clamp the distance, or do something
 *  else reasonable.
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after bevel_vertex or bevel_edge.
 * (So you can assume the local topology is set up however your bevel_* functions do it.)
 *
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::bevel_positions(FaceRef face, std::vector<Vec3> const &start_positions, Vec3 direction, float distance) {
	//A2Lx5h / A2Lx6h (OPTIONAL): Bevel Positions Helper
	
	// The basic strategy here is to loop over the list of outgoing halfedges,
	// and use the preceding and next vertex position from the original mesh
	// (in the start_positions array) to compute an new vertex position.
	
}

/*
 * extrude_positions: compute new positions for the vertices of an extruded face
 *  face: the face that was created by the extrude operation
 *  move: how much to translate the face
 *  shrink: amount to linearly interpolate vertices in the face toward the face's centroid
 *    shrink of zero leaves the face where it is
 *    positive shrink makes the face smaller (at shrink of 1, face is a point)
 *    negative shrink makes the face larger
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after extrude_face.
 * (So you can assume the local topology is set up however your extrude_face function does it.)
 *
 * Using extrude face in the GUI will assume a shrink of 0 to only extrude the selected face
 * Using bevel face in the GUI will allow you to shrink and increase the size of the selected face
 * 
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::extrude_positions(FaceRef face, Vec3 move, float shrink) {
	//A2L4h: Extrude Positions Helper

	//General strategy:
	// use mesh navigation to get starting positions from the surrounding faces,
	// compute the centroid from these positions + use to shrink,
	// offset by move

	// collect
	std::vector<VertexRef> face_vertices;
	HalfedgeRef h = face->halfedge;
	do
	{
		face_vertices.emplace_back(h->vertex);
		h = h->next;
	} while (h != face->halfedge);

	Vec3 center_pt = face->center();

	// shrink
	for (auto v : face_vertices)
	{
		v->position = shrink * center_pt + (1 - shrink) * v->position;
		v->position += move;
	}
	
}

