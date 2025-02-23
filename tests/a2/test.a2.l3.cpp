#include "test.h"
#include "geometry/halfedge.h"

static void expect_collapse(Halfedge_Mesh &mesh, Halfedge_Mesh::EdgeRef edge, Halfedge_Mesh const &after) {
	if (auto ret = mesh.collapse_edge(edge)) {
		if (auto msg = mesh.validate()) {
			throw Test::error("Invalid mesh: " + msg.value().second);
		}
		// check mesh shape:
		if (auto difference = Test::differs(mesh, after, Test::CheckAllBits)) {
			throw Test::error("Resulting mesh did not match expected: " + *difference);
		}
	} else {
		throw Test::error("collapse_edge rejected operation!");
	}
}

/*
BASIC CASE

Initial mesh:
0--1\
|  | \
2--3--4
|  | /
5--6/

Collapse Edge on Edge: 2-3

After mesh:
0-----1\
 \   /  \
  \ /    \
   2------3
  / \    /
 /   \  /
4-----5/
*/
Test test_a2_l3_collapse_edge_basic_simple("a2.l3.collapse_edge.basic.simple", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		  Vec3(-1.0f, 1.0f, 0.0f), 	Vec3(1.1f, 1.0f, 0.0f),
		 Vec3(-1.2f, 0.0f, 0.0f),   	 Vec3(1.2f, 0.0f, 0.0f),  Vec3(2.3f, 0.0f, 0.0f),
		Vec3(-1.4f,-1.0f, 0.0f), 		Vec3(1.5f, -1.0f, 0.0f)
	}, {
		{0, 2, 3, 1}, 
		{2, 5, 6, 3}, 
		{1, 3, 4}, 
		{3, 6, 4}
	});

	Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		  Vec3(-1.0f, 1.0f, 0.0f), 	Vec3(1.1f, 1.0f, 0.0f),
		 			Vec3(0.0f, 0.0f, 0.0f),  			Vec3(2.3f, 0.0f, 0.0f),
		Vec3(-1.4f,-1.0f, 0.0f), 		Vec3(1.5f, -1.0f, 0.0f)
	}, {
		{0, 2, 1}, 
		{2, 4, 5}, 
		{1, 2, 3}, 
		{2, 5, 3}
	});

	expect_collapse(mesh, edge, after);
});

/*
EDGE CASE

Initial mesh:
0--1\
|\ | \
| \2--3
|  | /
4--5/

Collapse Edge on Edge: 0-1

After mesh:
    0--\
   / \  \
  /   \  \
 /     1--2
/      | /
3------4/
*/
Test test_a2_l3_collapse_edge_edge_boundary("a2.l3.collapse_edge.edge.boundary", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                         Vec3(1.2f, 0.0f, 0.0f),  Vec3(2.3f, 0.0f, 0.0f),
		Vec3(-1.4f,-0.7f, 0.0f), Vec3(1.5f, -1.0f, 0.0f)
	}, {
		{0, 2, 1}, 
		{0, 4, 5, 2}, 
		{1, 2, 3}, 
		{2, 5, 3}
	});

	Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		       Vec3(0.05f, 1.05f, 0.0f), 
		                         Vec3(1.2f, 0.0f, 0.0f),  Vec3(2.3f, 0.0f, 0.0f),
		Vec3(-1.4f,-0.7f, 0.0f), Vec3(1.5f, -1.0f, 0.0f)
	}, {
		{0, 1, 2}, 
		{0, 3, 4, 1}, 
		{1, 4, 2}
	});

	expect_collapse(mesh, edge, after);
});


/*
0---1---2
|   |   |
|   |   |
3---4---5
|   |   |
|   |   |
6---7---8

collapse 3-4

0---1---2
|   |   |
 \ /    |
  3-----4
 / \    |
|   |   |
5---6---7

*/
Test test_a2_l3_collapse_foursquares_internal("a2.l3.collapse_edge.foursquares.internal", []() {
    Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
        Vec3(-1.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),
        Vec3(-1.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(-1.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f),
    }, {
        {0, 3, 4, 1},
        {3, 6, 7, 4},
        {1, 4, 5, 2},
        {4, 7, 8, 5}
    });

    Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

    Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
        Vec3(-1.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),
        Vec3(-0.5f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(-1.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f),     
    }, {
        {0, 3, 1},
        {3, 5, 6},
        {1, 3, 4, 2},
        {3, 6, 7, 4}
    });
    expect_collapse(mesh, edge, after);
});

/*
0---1---2
|   |   |
|   |   |
3---4---5
|   |   |
|   |   |
6---7---8

collapse 0-1

  0-----1
 / \    |
|   |   |
2---3---4
|   |   |
|   |   |
5---6---7

*/
Test test_a2_l3_collapse_foursquares_topleft_full("a2.l3.collapse_edge.foursquares.topleft.full", []() {
    Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
        Vec3(-1.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),
        Vec3(-1.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(-1.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f),
    }, {
        {0, 3, 4, 1},
        {3, 6, 7, 4},
        {1, 4, 5, 2},
        {4, 7, 8, 5}
    });

    Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->next->edge;

    Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
        Vec3(-0.5f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),
        Vec3(-1.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(-1.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f),     
    }, {
        {0, 2, 3},
        {0, 3, 4, 1},
        {2, 5, 6, 3},
        {3, 6, 7, 4}
    });
    expect_collapse(mesh, edge, after);
});

/*
    0---1
  / |   |
 /  |   |
2---3---4
|   |   |
|   |   |
5---6---7

collapse 2-3

    0---1
    |   |
   /    |
  2-----3
 / \    |
/   \   |
4---5---6

*/
Test test_a2_l3_collapse_foursquares_topleft_triangle("a2.l3.collapse_edge.foursquares.topleft.triangle", []() {
    Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
        Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),
        Vec3(-1.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(-1.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f),
    }, {
        {0, 2, 3},
        {0, 3, 4, 1},
        {2, 5, 6, 3},
        {3, 6, 7, 4}
    });

    Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

    Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
        Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),
        Vec3(-0.5f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(-1.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f),     
    }, {
        {2, 4, 5},
        {0, 2, 3, 1},
        {2, 5, 6, 3}
    });
    expect_collapse(mesh, edge, after);
});

/*
0---1---2
|   |   |
|   |   |
3---4---5
 \  |   |
   \|   |
    6---7

collapse 3-4

0---1---2
|   |   |
 \ /    |
  3-----4
  \     |
   \    |
    5---6
*/

Test test_a2_l3_collapse_foursquares_botleft_triangle("a2.l3.collapse_edge.foursquares.botleft.triangle", []() {
    Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
        Vec3(-1.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),
        Vec3(-1.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f),
    }, {
        {0, 3, 4, 1},
        {1, 4, 5, 2},
        {3, 6, 4},
        {4, 6, 7, 5}
    });

    Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

    Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
        Vec3(-1.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),
        Vec3(-0.5f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f),
    }, {
        {0, 3, 1},
        {1, 3, 4, 2},
        {3, 5, 6, 4}
    });

    expect_collapse(mesh, edge, after);
});

/*
    0---1
  / |   |
 /  |   |
2---3---4
 \  |   |
  \ |   |
    5---6

collapse 2-3

    0---1
    /   |
   /    |
  2-----3
   \    |
    \   |
    4---5

*/

Test test_a2_l3_collapse_foursquares_triangles_internal("a2.l3.collapse_edge.foursquares.triangles.internal", []() {
    Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
        Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),
        Vec3(-1.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f),
    }, {
        {0, 2, 3},
        {2, 5, 3},
        {0, 3, 4, 1},
        {3, 5, 6, 4}
    });

    Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

    Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
        Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),
        Vec3(-0.5f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f),       
    }, {
        {0, 2, 3, 1},
        {2, 4, 5, 3}
    });
    expect_collapse(mesh, edge, after);
});

/*
   1
 / | \
0  |  3
 \ | /
   2

collapse 0 - 2

   1
   | \
  /_ 3
 0
*/

Test test_a2_l3_collapse_twotriangles_boundary("a2.l3.collapse_edge.twotriangles.boundary", []() {
    Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
        Vec3(-1.0f, 0.0f, 0.0f), 
        Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f),
        Vec3(1.0f, 0.0f, 0.0f)
    }, {
        {0, 2, 1},
        {1, 2, 3}
    });
    Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->edge;

    Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
        Vec3(-0.5f, -0.5f, 0.0f), 
        Vec3(0.0f, 1.0f, 0.0f),
        Vec3(1.0f, 0.0f, 0.0f)
    }, {
        {0, 2, 1}
    });
    expect_collapse(mesh, edge, after);
});

/*
  4
 / \
2---3
| \ |
|  \|
0---1
collapse 0 - 1

  3
 / \
1---2
 \ /
  0
*/

Test test_a2_l3_collapse_pentagon_boundary("a2.l3.collapse_edge.pentagon.boundary", []() {
    Halfedge_Mesh pentagon = Halfedge_Mesh::from_indexed_faces({
        Vec3(0.0f, 0.0f, 0.0f), Vec3(2.0f, 0.0f, 0.0f),
        Vec3(0.0f, 1.0f, 0.0f), Vec3(2.0f, 1.0f, 0.0f),
        Vec3(1.0f, 2.0f, 0.0f)

    }, {
        {0, 1, 2},
        {1, 3, 2},
        {2, 3, 4}
    });

    Halfedge_Mesh::EdgeRef edge = pentagon.halfedges.begin()->edge;
    
    Halfedge_Mesh doubletriangle = Halfedge_Mesh::from_indexed_faces({
        Vec3(1.0f, 0.0f, 0.0f),
        Vec3(0.0f, 1.0f, 0.0f), Vec3(2.0f, 1.0f, 0.0f),
        Vec3(1.0f, 2.0f, 0.0f)

    }, {
        {0, 2, 1},
        {1, 2, 3}
    }); 
    expect_collapse(pentagon, edge, doubletriangle);
});

Test test_a2_l3_collapse_triangle("a2.l3.collapse_edge.triangle", []() {
    Halfedge_Mesh triangle = Halfedge_Mesh::from_indexed_faces({
        Vec3(0.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        Vec3(0.0f, 1.0f, 0.0f), 
    }, {
        {0, 1, 2},
    });

    Halfedge_Mesh::EdgeRef edge = triangle.halfedges.begin()->edge;

	Halfedge_Mesh doubletriangle = Halfedge_Mesh::from_indexed_faces({
        Vec3(1.0f, 0.0f, 0.0f),
        Vec3(0.0f, 1.0f, 0.0f), Vec3(2.0f, 1.0f, 0.0f),
        Vec3(1.0f, 2.0f, 0.0f)

    }, {
        {0, 2, 1},
        {1, 2, 3}
    }); 
    
    expect_collapse(triangle, edge, doubletriangle);
});