#include "test.h"
#include "rasterizer/pipeline.h"
#include "rasterizer/programs.h"

#include <limits>
#include <iomanip>
#include <algorithm>
#include <unordered_set>

using TestPipeline = Pipeline< PrimitiveType::Lines, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat >;

namespace std {
	template< >
	struct hash< Vec2 > {
		size_t operator()(const Vec2 &v) const {
			static hash< float > hf;
			size_t x = hf(v.x);
			size_t y = hf(v.y);
			return x ^ (y << 16) ^ (y >> (sizeof(y)*8-16));
		}
	};
}

//check that line produces exactly the listed fragments:
void check_line_covers(std::string const &desc, std::vector< Vec2 > const &line_strip, std::unordered_set< Vec2 > const &expected) {

	std::unordered_set< Vec2 > got;
	for (uint32_t i = 0; i + 1 < line_strip.size(); ++i) {
		TestPipeline::ClippedVertex a,b;
		a.fb_position = Vec3(line_strip[i].x, line_strip[i].y, 0.25f);
		a.inv_w = 1.0f;
		a.attributes.fill(1.0f);
		b.fb_position = Vec3(line_strip[i+1].x, line_strip[i+1].y, 0.75f);
		b.inv_w = 1.0f;
		b.attributes.fill(2.0f);
		TestPipeline::rasterize_line(a, b, [&](TestPipeline::Fragment const &frag){
			got.emplace(frag.fb_position.x, frag.fb_position.y);
		});
	}

	std::vector< std::string > raster;

	raster.emplace_back(".");

	uint32_t out_of_raster = 0;

	auto draw = [&raster,&out_of_raster](Vec2 const &px, char c) {
		int32_t x = int32_t(std::floor(px.x));
		int32_t y = int32_t(std::floor(px.y));

		if (x < 0 || y < 0 || x > 10 || y > 10) {
			++out_of_raster;
			return;
		}

		if (uint32_t(y) >= raster.size()) {
			raster.resize(y+1, "");
		}
		if (uint32_t(x) >= raster[y].size()) {
			raster[y].resize(x+1, '.');
		}
		raster[y][x] = c;
	};

	uint32_t matched = 0;
	uint32_t missed = 0;
	uint32_t extra = 0;

	for (auto const &f : got) {
		if ((f.x - std::floor(f.x) != 0.5f) || (f.y - std::floor(f.y) != 0.5f)) {
			throw Test::error("Rasterizing '" + desc + "', got fragment at (" + std::to_string(f.x) + ", " + std::to_string(f.y) + "), which isn't at a pixel center.");
		}
		if (expected.count(f)) {
			draw(f, '#');
			++matched;
		} else {
			draw(f, '!');
			++extra;
		}
	}
	for (auto const &f : expected) {
		if (!got.count(f)) {
			draw(f, '?');
			++missed;
		}
	}

	if (extra > 0 || missed > 0) {
		//failed!
		std::string info = "Example '" + desc + "' missed " + std::to_string(missed) + " ('?'); had " + std::to_string(extra) + " extra ('!'); and matched " + std::to_string(matched) + " ('#') fragments:";

		//square up the raster:
		size_t width = 0;
		for (auto const &row : raster) {
			width = std::max(width, row.size());
		}
		for (auto &row : raster) {
			row.resize(width, '.');
		}

		for (uint32_t y = static_cast<uint32_t>(raster.size()) - 1; y < static_cast<uint32_t>(raster.size()); --y) {
			info += "\n    " + raster[y];
		}

		if (out_of_raster) info += "\n    (" + std::to_string(out_of_raster) + " out-of-range fragments not plotted.)";

		puts(""); //because "test..."
		info("%s", info.c_str());

		throw Test::error("Example '" + desc + "' didn't match expected.");
	}

	//if nothing extra and nothing missed, success!
	assert(matched == expected.size());
}

//check that line produces exactly the fragments drawn in a fancy picture:
void check_line_covers(std::string const &desc, std::initializer_list< Vec2 > const &line_strip, std::initializer_list< std::string > const &raster_) {
	//convert raster to set of points ( with lower-left being (0,0) ):
	std::vector< std::string > raster(raster_);
	std::unordered_set< Vec2 > expected;
	for (uint32_t y = 0; y < raster.size(); ++y) {
		std::string const &row = raster[raster.size()-1-y];
		for (uint32_t x = 0; x < row.size(); ++x) {
			if (row[x] != '.') {
				expected.emplace(x + 0.5f, y + 0.5f);
			}
		}
	}
	//use list-of-points version:
	check_line_covers(desc, line_strip, expected);
}

//--------------------------------------------------
//entering/exiting diamond at (1,1):
// only lines that *exit* the diamond should produce a fragment.


Test test_a1_task2_diamond_inside("a1.task2.diamond.inside", []() {
	check_line_covers(
		"line inside diamond (1,1)",
		{ Vec2(1.5f, 1.25f), Vec2(1.25f, 1.5f) },
		{"...",
		 "...",
		 "..."}
	);
});


Test test_a1_task2_diamond_outside("a1.task2.diamond.outside", []() {
	check_line_covers(
		"line outside diamond (1,1)",
		{ Vec2(1.125f, 1.25f), Vec2(1.25f, 1.125f) },
		{"...",
		 "...",
		 "..."}
	);
});


//----------------------------
//simple horizontal and vertical lines (set up so that no enter/exit logic needed):

Test test_a1_task2_simple_horizontal("a1.task2.simple.horizontal", []() {
	check_line_covers(
		"horizontal line from (1.125, 1.125) to (4.875, 1.125)",
		{ Vec2(1.125f, 1.125f), Vec2(4.875f, 1.125f) },
		{"......",
		 ".####.",
		 "......"}
	);
});


Test test_a1_task2_simple_vertical("a1.task2.simple.vertical", []() {
	check_line_covers(
		"vertical line from (1.125, 1.125) to (1.125, 4.875)",
		{ Vec2(1.125f, 1.125f), Vec2(1.125f, 4.875f) },
		{"...",
		 ".#.",
		 ".#.",
		 ".#.",
		 ".#.",
		 "..."}
	);
});

Test test_a1_task2_diamond_edge("a1.task2.diamond.edge", []() {
    check_line_covers(
        "line beginning and ending on diamond edge from (0.75, 0.75) to (1.25, 1.25)",
        { Vec2(0.75f, 0.75f), Vec2(1.25f, 1.25f) },
        {"...",
         "...",
         "#.."}
    );
});

Test test_a1_task2_diamond_edge_same("a1.task2.diamond.edge.same", []() {
    check_line_covers(
        "line beginning and ending on edges of the _same_ diamond from (0.25, 0.25) to (0.75, 0.75)",
        { Vec2(0.25f, 0.25f), Vec2(0.75f, 0.75f) },
        {"...",
         "...",
         "..."}
    );
});

Test test_a1_task2_diamond_edge_vertex("a1.task2.diamond.edge.vertex", []() {
    check_line_covers(
        "line beginning on vertex edge and ending on right vertex of the same diamond (0.25, 0.25) to (1.0, 0.5)",
        { Vec2(0.25f, 0.25f), Vec2(1.0f, 0.5f) },
        {"...",
         "...",
         "#.."}
    );
});

Test test_a1_task2_diamond_vertex_top("a1.task2.diamond.vertex.top", []() {
    check_line_covers(
        "line tangent to top diamond vertex from (0.0, 1.0) to (2.8, 1.0)",
        { Vec2(0.0f, 1.0f), Vec2(2.8f, 1.0f) },
        {"...",
         "###",
         "..."}
    );
});

Test test_a1_task2_diamond_vertex_left("a1.task2.diamond.vertex.left", []() {
    check_line_covers(
        "line tangent to left diamond vertex from (1.0, .125) to (1.0, 2.875)",
        { Vec2(1.0f, .125f), Vec2(1.0f, 2.875f) },
        {".#.",
         ".#.",
         ".#."}
    );
});


Test test_a1_task2_single_pixel_1("a1.task2.single.pixel.1", []() {
    check_line_covers(
        "line within one pixel from (2.5, 2.5) to (2.975, 2.125)",
        { Vec2(2.5f, 2.5f), Vec2(2.975f, 2.125f) },
        {"....",
         "..#.",
         "....",
         "...."}
    );
});

Test test_a1_task2_single_pixel_2("a1.task2.single.pixel.2", []() {
    check_line_covers(
        "line within one pixel from (1.125, 1.25) to (1.875, 1.25)",
        { Vec2(1.125f, 1.25f), Vec2(1.875f, 1.25f) },
        {"...",
         ".#.",
         "..."}
    );
});

Test test_a1_task2_single_pixel_3("a1.task2.single.pixel.3", []() {
    check_line_covers(
        "just one point inside the diamond",
        { Vec2(1.5f, 1.5f), Vec2(1.5f, 1.5f) },
        {"...",
         "...",
         "..."}
    );
});

Test test_a1_task2_single_pixel_4("a1.task2.single.pixel.4", []() {
    check_line_covers(
        "just one point outside the diamond",
        { Vec2(1.125f, 1.25f), Vec2(1.125f, 1.25f) },
        {"...",
         "...",
         "..."}
    );
});

Test test_a1_task2_single_pixel_5("a1.task2.single.pixel.5", []() {
    check_line_covers(
        "line within one pixel from (1.125, 1.25) to (1.5, 1.25)",
        { Vec2(1.125f, 1.25f), Vec2(1.5f, 1.25f) },
        {"...",
         "...",
         "..."}
    );
});

Test test_a1_task2_single_pixel_6("a1.task2.single.pixel.6", []() {
    check_line_covers(
        "line within one pixel from (1.5, 1.25) to (1.875, 1.3)",
        { Vec2(1.5f, 1.25f), Vec2(1.875f, 1.3f) },
        {"...",
         ".#.",
         "..."}
    );
});


