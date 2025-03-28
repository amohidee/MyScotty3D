// clang-format off
#include "pipeline.h"

#include <iostream>

#include "../lib/log.h"
#include "../lib/mathlib.h"
#include "framebuffer.h"
#include "sample_pattern.h"
template<PrimitiveType primitive_type, class Program, uint32_t flags>
void Pipeline<primitive_type, Program, flags>::run(std::vector<Vertex> const& vertices,
                                                   typename Program::Parameters const& parameters,
                                                   Framebuffer* framebuffer_) {
	// Framebuffer must be non-null:
	assert(framebuffer_);
	auto& framebuffer = *framebuffer_;

	// A1T7: sample loop
	// TODO: update this function to rasterize to *all* sample locations in the framebuffer.
	//  	 This will probably involve inserting a loop of the form:
	// 		 	std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	//      	for (uint32_t s = 0; s < samples.size(); ++s) { ... }
	//   	 around some subset of the code.
	// 		 You will also need to transform the input and output of the rasterize_* functions to
	// 	     account for the fact they deal with pixels centered at (0.5,0.5).

	std::vector<ShadedVertex> shaded_vertices;
	shaded_vertices.reserve(vertices.size());

	//--------------------------
	// shade vertices:
	for (auto const& v : vertices) {
		ShadedVertex sv;
		Program::shade_vertex(parameters, v.attributes, &sv.clip_position, &sv.attributes);
		shaded_vertices.emplace_back(sv);
	}

	//--------------------------
	// assemble + clip + homogeneous divide vertices:
	std::vector<ClippedVertex> clipped_vertices;

	// reserve some space to avoid reallocations later:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		// clipping lines can never produce more than one vertex per input vertex:
		clipped_vertices.reserve(shaded_vertices.size());
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		// clipping triangles can produce up to 8 vertices per input vertex:
		clipped_vertices.reserve(shaded_vertices.size() * 8);
	}
	// clang-format off

	//coefficients to map from clip coordinates to framebuffer (i.e., "viewport") coordinates:
	//x: [-1,1] -> [0,width]
	//y: [-1,1] -> [0,height]
	//z: [-1,1] -> [0,1] (OpenGL-style depth range)
	Vec3 const clip_to_fb_scale = Vec3{
		framebuffer.width / 2.0f,
		framebuffer.height / 2.0f,
		0.5f
	};
	Vec3 const clip_to_fb_offset = Vec3{
		0.5f * framebuffer.width,
		0.5f * framebuffer.height,
		0.5f
	};

	// helper used to put output of clipping functions into clipped_vertices:
	auto emit_vertex = [&](ShadedVertex const& sv) {
		ClippedVertex cv;
		float inv_w = 1.0f / sv.clip_position.w;
		cv.fb_position = clip_to_fb_scale * inv_w * sv.clip_position.xyz() + clip_to_fb_offset;
		cv.inv_w = inv_w;
		cv.attributes = sv.attributes;
		clipped_vertices.emplace_back(cv);
	};

	// actually do clipping:
	if constexpr (primitive_type == PrimitiveType::Lines) {
		for (uint32_t i = 0; i + 1 < shaded_vertices.size(); i += 2) {
			clip_line(shaded_vertices[i], shaded_vertices[i + 1], emit_vertex);
		}
	} else if constexpr (primitive_type == PrimitiveType::Triangles) {
		for (uint32_t i = 0; i + 2 < shaded_vertices.size(); i += 3) {
			clip_triangle(shaded_vertices[i], shaded_vertices[i + 1], shaded_vertices[i + 2], emit_vertex);
		}
	} else {
		static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
	}

	//--------------------------
	// rasterize primitives:

	std::vector< Vec3 > const &samples = framebuffer.sample_pattern.centers_and_weights;
	// actually do rasterization:
	// if constexpr (primitive_type == PrimitiveType::Lines) {
	// 	for (uint32_t i = 0; i + 1 < clipped_vertices.size(); i += 2) {
	// 		rasterize_line(clipped_vertices[i], clipped_vertices[i + 1], emit_fragment);
	// 	}
	// } else if constexpr (primitive_type == PrimitiveType::Triangles) {
	// 	for (uint32_t i = 0; i + 2 < clipped_vertices.size(); i += 3) {
	// 		rasterize_triangle(clipped_vertices[i], clipped_vertices[i + 1], clipped_vertices[i + 2], emit_fragment);
	// 	}
	// } else {
	// 	static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
	// }

	// rasterization sampling
	// for each sample, shift the clipped vertex array
	for (uint32_t s = 0; s < samples.size(); ++s) {
		std::vector<Fragment> fragments;

		

		// helper used to put output of rasterization functions into fragments:
		auto emit_fragment = [&](Fragment const& f) { fragments.emplace_back(f); };
        Vec2 sample_offset = Vec2(samples[s].x, samples[s].y);

        std::vector<ClippedVertex> adjusted_vertices = clipped_vertices;
        for (auto& v : adjusted_vertices) {
			// offset vertex
            v.fb_position.x += (sample_offset.x);
            v.fb_position.y += (sample_offset.y);
        }

        if constexpr (primitive_type == PrimitiveType::Lines) {
            for (uint32_t i = 0; i + 1 < adjusted_vertices.size(); i += 2) {
                rasterize_line(adjusted_vertices[i], adjusted_vertices[i + 1], emit_fragment);
            }
        } else if constexpr (primitive_type == PrimitiveType::Triangles) {
            for (uint32_t i = 0; i + 2 < adjusted_vertices.size(); i += 3) {
                rasterize_triangle(adjusted_vertices[i], adjusted_vertices[i + 1], adjusted_vertices[i + 2], emit_fragment);
            }
        } else {
			static_assert(primitive_type == PrimitiveType::Lines, "Unsupported primitive type.");
		}

		//--------------------------
		// depth test + shade + blend fragments:
		uint32_t out_of_range = 0; // check if rasterization produced fragments outside framebuffer 
								// (indicates something is wrong with clipping)

		std::cout << "num frags: " << fragments.size();
		for (auto const& f : fragments) {

			// fragment location (in pixels):
			int32_t x = (int32_t)std::floor(f.fb_position.x);
			int32_t y = (int32_t)std::floor(f.fb_position.y);

			// if clipping is working properly, this condition shouldn't be needed;
			// however, it prevents crashes while you are working on your clipping functions,
			// so we suggest leaving it in place:
			if (x < 0 || (uint32_t)x >= framebuffer.width || 
				y < 0 || (uint32_t)y >= framebuffer.height) {
				//std::cout << "x: " << x << ",  " << "y: " << y << "\n";
				++out_of_range;
				continue;
			}


			// local names that refer to destination sample in framebuffer:
			float& fb_depth = framebuffer.depth_at(x, y, s);
			Spectrum& fb_color = framebuffer.color_at(x, y, s);


			// depth test:
			if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Always) {
				// "Always" means the depth test always passes.
			} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Never) {
				// "Never" means the depth test never passes.
				continue; //discard this fragment
			} else if constexpr ((flags & PipelineMask_Depth) == Pipeline_Depth_Less) {
				// "Less" means the depth test passes when the new fragment has depth less than the stored depth.
				// A1T4: Depth_Less
				// TODO: implement depth test! We want to only emit fragments that have a depth less than the stored depth, hence "Depth_Less".
				if (f.fb_position.z >= fb_depth) {
					continue;
				}
			} else {
				static_assert((flags & PipelineMask_Depth) <= Pipeline_Depth_Always, "Unknown depth test flag.");
			}

			// if depth test passes, and depth writes aren't disabled, write depth to depth buffer:
			if constexpr (!(flags & Pipeline_DepthWriteDisableBit)) {
				fb_depth = f.fb_position.z;
			}

			// shade fragment:
			ShadedFragment sf;
			sf.fb_position = f.fb_position;
			Program::shade_fragment(parameters, f.attributes, f.derivatives, &sf.color, &sf.opacity);

			// write color to framebuffer if color writes aren't disabled:
			if constexpr (!(flags & Pipeline_ColorWriteDisableBit)) {
				// blend fragment:
				if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Replace) {
					fb_color = sf.color;
				} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Add) {
					// A1T4: Blend_Add
					// TODO: framebuffer color should have fragment color multiplied by fragment opacity added to it.
					fb_color += sf.color * sf.opacity; //<-- replace this line
				} else if constexpr ((flags & PipelineMask_Blend) == Pipeline_Blend_Over) {
					// A1T4: Blend_Over
					// TODO: set framebuffer color to the result of "over" blending (also called "alpha blending") the fragment color over the framebuffer color, using the fragment's opacity
					// 		 You may assume that the framebuffer color has its alpha premultiplied already, and you just want to compute the resulting composite color
					fb_color = sf.color + (1.0f - sf.opacity) * fb_color; //<-- replace this line
				} else {
					static_assert((flags & PipelineMask_Blend) <= Pipeline_Blend_Over, "Unknown blending flag.");
				}
			}
		}
		if (out_of_range > 0) {
			if constexpr (primitive_type == PrimitiveType::Lines) {
				warn("Produced %d fragments outside framebuffer; this indicates something is likely "
					"wrong with the clip_line function.",
					out_of_range);
			} else if constexpr (primitive_type == PrimitiveType::Triangles) {
				std::cout << "fb width: " << framebuffer.width << ", " << "fb height: " << framebuffer.height << "\n";
				warn("Produced %d fragments outside framebuffer; this indicates something is likely "
					"wrong with the clip_triangle function.",
					out_of_range);
			}
		}
	}
}

// -------------------------------------------------------------------------
// clipping functions

// helper to interpolate between vertices:
template<PrimitiveType p, class P, uint32_t F>
auto Pipeline<p, P, F>::lerp(ShadedVertex const& a, ShadedVertex const& b, float t) -> ShadedVertex {
	ShadedVertex ret;
	ret.clip_position = (b.clip_position - a.clip_position) * t + a.clip_position;
	for (uint32_t i = 0; i < ret.attributes.size(); ++i) {
		ret.attributes[i] = (b.attributes[i] - a.attributes[i]) * t + a.attributes[i];
	}
	return ret;
}

/*
 * clip_line - clip line to portion with -w <= x,y,z <= w, emit vertices of clipped line (if non-empty)
 *  	va, vb: endpoints of line
 *  	emit_vertex: call to produce truncated line
 *
 * If clipping shortens the line, attributes of the shortened line should respect the pipeline's interpolation mode.
 * 
 * If no portion of the line remains after clipping, emit_vertex will not be called.
 *
 * The clipped line should have the same direction as the full line.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::clip_line(ShadedVertex const& va, ShadedVertex const& vb,
                                      std::function<void(ShadedVertex const&)> const& emit_vertex) {
	// Determine portion of line over which:
	// 		pt = (b-a) * t + a
	//  	-pt.w <= pt.x <= pt.w
	//  	-pt.w <= pt.y <= pt.w
	//  	-pt.w <= pt.z <= pt.w
	// ... as a range [min_t, max_t]:

	float min_t = 0.0f;
	float max_t = 1.0f;

	// want to set range of t for a bunch of equations like:
	//    a.x + t * ba.x <= a.w + t * ba.w
	// so here's a helper:
	auto clip_range = [&min_t, &max_t](float l, float dl, float r, float dr) {
		// restrict range such that:
		// l + t * dl <= r + t * dr
		// re-arranging:
		//  l - r <= t * (dr - dl)
		if (dr == dl) {
			// want: l - r <= 0
			if (l - r > 0.0f) {
				// works for none of range, so make range empty:
				min_t = 1.0f;
				max_t = 0.0f;
			}
		} else if (dr > dl) {
			// since dr - dl is positive:
			// want: (l - r) / (dr - dl) <= t
			min_t = std::max(min_t, (l - r) / (dr - dl));
		} else { // dr < dl
			// since dr - dl is negative:
			// want: (l - r) / (dr - dl) >= t
			max_t = std::min(max_t, (l - r) / (dr - dl));
		}
	};

	// local names for clip positions and their difference:
	Vec4 const& a = va.clip_position;
	Vec4 const& b = vb.clip_position;
	Vec4 const ba = b - a;

	// -a.w - t * ba.w <= a.x + t * ba.x <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.x, ba.x);
	clip_range(a.x, ba.x, a.w, ba.w);
	// -a.w - t * ba.w <= a.y + t * ba.y <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.y, ba.y);
	clip_range(a.y, ba.y, a.w, ba.w);
	// -a.w - t * ba.w <= a.z + t * ba.z <= a.w + t * ba.w
	clip_range(-a.w, -ba.w, a.z, ba.z);
	clip_range(a.z, ba.z, a.w, ba.w);

	if (min_t < max_t) {
		if (min_t == 0.0f) {
			emit_vertex(va);
		} else {
			ShadedVertex out = lerp(va, vb, min_t);
			// don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
				out.attributes = va.attributes;
			}
			emit_vertex(out);
		}
		if (max_t == 1.0f) {
			emit_vertex(vb);
		} else {
			ShadedVertex out = lerp(va, vb, max_t);
			// don't interpolate attributes if in flat shading mode:
			if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
				out.attributes = va.attributes;
			}
			emit_vertex(out);
		}
	}
}

/*
 * clip_triangle - clip triangle to portion with -w <= x,y,z <= w, emit resulting shape as triangles (if non-empty)
 *  	va, vb, vc: vertices of triangle
 *  	emit_vertex: call to produce clipped triangles (three calls per triangle)
 *
 * If clipping truncates the triangle, attributes of the new vertices should respect the pipeline's interpolation mode.
 * 
 * If no portion of the triangle remains after clipping, emit_vertex will not be called.
 *
 * The clipped triangle(s) should have the same winding order as the full triangle.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::clip_triangle(
	ShadedVertex const& va, ShadedVertex const& vb, ShadedVertex const& vc,
	std::function<void(ShadedVertex const&)> const& emit_vertex) {
	// A1EC: clip_triangle
	// TODO: correct code!
	emit_vertex(va);
	emit_vertex(vb);
	emit_vertex(vc);
}

// -------------------------------------------------------------------------
// rasterization functions

/*
 * rasterize_line:
 * calls emit_fragment( frag ) for every pixel "covered" by the line (va.fb_position.xy, vb.fb_position.xy).
 *
 *    a pixel (x,y) is "covered" by the line if it exits the inscribed diamond:
 * 
 *        (x+0.5,y+1)
 *        /        \
 *    (x,y+0.5)  (x+1,y+0.5)
 *        \        /
 *         (x+0.5,y)
 *
 *    to avoid ambiguity, we consider diamonds to contain their left and bottom points
 *    but not their top and right points. 
 * 
 * 	  since 45 degree lines breaks this rule, our rule in general is to rasterize the line as if its
 *    endpoints va and vb were at va + (e, e^2) and vb + (e, e^2) where no smaller nonzero e produces 
 *    a different rasterization result. 
 *    We will not explicitly check for 45 degree lines along the diamond edges (this will be extra credit),
 *    but you should be able to handle 45 degree lines in every other case (such as starting from pixel centers)
 *
 * for each such diamond, pass Fragment frag to emit_fragment, with:
 *  - frag.fb_position.xy set to the center (x+0.5,y+0.5)
 *  - frag.fb_position.z interpolated linearly between va.fb_position.z and vb.fb_position.z
 *  - frag.attributes set to va.attributes (line will only be used in Interp_Flat mode)
 *  - frag.derivatives set to all (0,0)
 *
 * when interpolating the depth (z) for the fragments, you may use any depth the line takes within the pixel
 * (i.e., you don't need to interpolate to, say, the closest point to the pixel center)
 *
 * If you wish to work in fixed point, check framebuffer.h for useful information about the framebuffer's dimensions.
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::rasterize_line(
	ClippedVertex const& va, ClippedVertex const& vb,
	std::function<void(Fragment const&)> const& emit_fragment) {
	if constexpr ((flags & PipelineMask_Interp) != Pipeline_Interp_Flat) {
		assert(0 && "rasterize_line should only be invoked in flat interpolation mode.");
	}
	// A1T2: rasterize_line

	// TODO: Check out the block comment above this function for more information on how to fill in
	// this function!
	// The OpenGL specification section 3.5 may also come in handy.

	// Extract integer coordinates
    float x1 = va.fb_position.x;
    float y1 = va.fb_position.y;
    float x2 = vb.fb_position.x;
    float y2 = vb.fb_position.y;

    float z1 = va.fb_position.z;
    float z2 = vb.fb_position.z;

	// Compute differences
    float dx = x2 - x1;
    float dy = y2 - y1;
	float abs_dx = std::abs(dx);
    float abs_dy = std::abs(dy);
    //int sx = (dx > 0) ? 1 : -1;
    //int sy = (dy > 0) ? 1 : -1;

	if (abs_dx > abs_dy) {
		// X MAJOR
		// i = X axis
		// left-to-right order
		if (x1 > x2) {
			std::swap(x1, x2);
			std::swap(y1, y2);
			std::swap(z1, z2);
		}
		float t1 = std::floor(x1);
		float t2 = std::floor(x2);
		// calculate beginning and end diamond exit
		float px1 = std::floor(x1) + 0.5f;
		float py1 = std::floor(y1) + 0.5f;
		float px2 = std::floor(x2) + 0.5f;
		float py2 = std::floor(y2) + 0.5f;
		//same x square edge case
		if (px1 == px2) {
			if ((x1 - px1) + std::abs(y1 - py1) < 0.5f && std::abs(x2 - px2) + std::abs(y2 - py2) >= 0.5f) {
				Vec3 pixel_center(px1, py1, z1);


				Fragment frag;
				frag.fb_position = pixel_center;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				emit_fragment(frag);
			}
			return;
		}
		if ((x1 - px1) + std::abs(y1 - py1) < 0.5f) {
			Vec3 pixel_center(px1, py1, z1);


			Fragment frag;
			frag.fb_position = pixel_center;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			emit_fragment(frag);
		}
		if (-(x2 - px2) + std::abs(y2 - py2) < 0.5f) {
			Vec3 pixel_center(px2, py2, z1);


			Fragment frag;
			frag.fb_position = pixel_center;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			emit_fragment(frag);
		}
		for (float u = t1; u < t2; u++) {
			float w = ((u + 0.5f) - x1) / (x2 -x1);
			float v = w * (y2-y1) + y1;
			float z = z1 * (1 - w) + z2 * w;

			Vec3 pixel_center(std::floor(u) + 0.5f, std::floor(v) + 0.5f, z);


			Fragment frag;
			frag.fb_position = pixel_center;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			emit_fragment(frag);


		}
	} else {
		// Y MAJOR
		// i = Y axis
		if (y1 > y2) {
			std::swap(x1, x2);
			std::swap(y1, y2);
			std::swap(z1, z2);
		}
		float t1 = std::floor(y1);
		float t2 = std::floor(y2);
		// calculate beginning and end diamond exit
		float px1 = std::floor(x1) + 0.5f;
		float py1 = std::floor(y1) + 0.5f;
		float px2 = std::floor(x2) + 0.5f;
		float py2 = std::floor(y2) + 0.5f;
		// same y square edge case
		if (py1 == py2) {
			if (std::abs(x1 - px1) + std::abs(y1 - py1) < 0.5f && std::abs(x2 - px2) + std::abs(y2 - py2) >= 0.5f) {
				Vec3 pixel_center(px1, py1, z1);


				Fragment frag;
				frag.fb_position = pixel_center;
				frag.attributes = va.attributes;
				frag.derivatives.fill(Vec2(0.0f, 0.0f));
				emit_fragment(frag);
			}
			return;
		}
		if (std::abs(x1 - px1) + (y1 - py1) < 0.5f) {
			Vec3 pixel_center(px1, py1, z1);


			Fragment frag;
			frag.fb_position = pixel_center;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			emit_fragment(frag);
		}
		if (std::abs(x2 - px2) + -(y2 - py2) < 0.5f) {
			Vec3 pixel_center(px2, py2, z1);


			Fragment frag;
			frag.fb_position = pixel_center;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			emit_fragment(frag);
		}
		for (float u = t1; u < t2; u++) {
			float w = ((u + 0.5f) - y1) / (y2 -y1);
			float v = w * (x2-x1) + x1;
			float z = z1 * (1 - w) + z2 * w;

			Vec3 pixel_center(std::floor(v) + 0.5f, std::floor(u) + 0.5f, z);

			Fragment frag;
			frag.fb_position = pixel_center;
			frag.attributes = va.attributes;
			frag.derivatives.fill(Vec2(0.0f, 0.0f));
			emit_fragment(frag);
		}
	}

}

/*
 * rasterize_triangle(a,b,c,emit) calls 'emit(frag)' at every location
 *  	(x+0.5,y+0.5) (where x,y are integers) covered by triangle (a,b,c).
 *
 * The emitted fragment should have:
 * - frag.fb_position.xy = (x+0.5, y+0.5)
 * - frag.fb_position.z = linearly interpolated fb_position.z from a,b,c (NOTE: does not depend on Interp mode!)
 * - frag.attributes = depends on Interp_* flag in flags:
 *   - if Interp_Flat: copy from va.attributes
 *   - if Interp_Smooth: interpolate as if (a,b,c) is a 2D triangle flat on the screen
 *   - if Interp_Correct: use perspective-correct interpolation
 * - frag.derivatives = derivatives w.r.t. fb_position.x and fb_position.y of the first frag.derivatives.size() attributes.
 *
 * Notes on derivatives:
 * 	The derivatives are partial derivatives w.r.t. screen locations. That is:
 *    derivatives[i].x = d/d(fb_position.x) attributes[i]
 *    derivatives[i].y = d/d(fb_position.y) attributes[i]
 *  You may compute these derivatives analytically or numerically.
 *
 *  See section 8.12.1 "Derivative Functions" of the GLSL 4.20 specification for some inspiration. (*HOWEVER*, the spec is solving a harder problem, and also nothing in the spec is binding on your implementation)
 *
 *  One approach is to rasterize blocks of four fragments and use forward and backward differences to compute derivatives.
 *  To assist you in this approach, keep in mind that the framebuffer size is *guaranteed* to be even. (see framebuffer.h)
 *
 * Notes on coverage:
 *  If two triangles are on opposite sides of the same edge, and a
 *  fragment center lies on that edge, rasterize_triangle should
 *  make sure that exactly one of the triangles emits that fragment.
 *  (Otherwise, speckles or cracks can appear in the final render.)
 * 
 *  For degenerate (co-linear) triangles, you may consider them to not be on any side of an edge.
 * 	Thus, even if two degnerate triangles share an edge that contains a fragment center, you don't need to emit it.
 *  You will not lose points for doing something reasonable when handling this case
 *
 *  This is pretty tricky to get exactly right!
 *
 */
template<PrimitiveType p, class P, uint32_t flags>
void Pipeline<p, P, flags>::rasterize_triangle(
	ClippedVertex const& va, ClippedVertex const& vb, ClippedVertex const& vc,
	std::function<void(Fragment const&)> const& emit_fragment) {
	// NOTE: it is okay to restructure this function to allow these tasks to use the
	//  same code paths. Be aware, however, that all of them need to remain working!
	//  (e.g., if you break Flat while implementing Correct, you won't get points
	//   for Flat.)
	if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Flat) {
		// A1T3: flat triangles
		// TODO: rasterize triangle (see block comment above this function).

        float ax = va.fb_position.x, ay = va.fb_position.y; 
        float bx = vb.fb_position.x, by = vb.fb_position.y; 
        float cx = vc.fb_position.x, cy = vc.fb_position.y; 

		float az = va.fb_position.z;
		float bz = vb.fb_position.z;
		float cz = vc.fb_position.z;

		 // bounding box
        float minX = std::floor(std::min({ax, bx, cx}));
        float minY = std::floor(std::min({ay, by, cy}));
        float maxX = std::ceil(std::max({ax, bx, cx}));
        float maxY = std::ceil(std::max({ay, by, cy}));

		
		auto cross_product = [](Vec2 v1, Vec2 v2) -> float {
        	return v1.x * v2.y - v1.y * v2.x;
    	};

		auto point_in_triangle = [&](float x , float y) -> bool {
            Vec2 ac (cx - ax, cy - ay);
            Vec2 ab (bx - ax, by - ay);
            Vec2 cb (bx - cx, by - cy);
            Vec2 ca (ax - cx, ay - cy);
            Vec2 ba (ax - bx, ay - by);
            Vec2 bc (cx - bx, cy - by);

            Vec2 aq (x - ax, y - ay);
            Vec2 cq (x - cx, y - cy);
            Vec2 bq (x - bx, y - by);

            float w1 = cross_product(ac, ab) * cross_product(ac, aq);
            float w2 = cross_product(cb, ca) * cross_product(cb, cq);
            float w3 = cross_product(ba, bc) * cross_product(ba, bq);

            return (w1 > 0 && w2 > 0 && w3 > 0) || (w1 < 0 && w2 < 0 && w3 < 0);
        };

		// winding > 0 = CCW, winding < 0 = CW
		Vec2 ab (bx - ax, by - ay);
		Vec2 ac (cx - ax, cy - ay);
		float winding = cross_product(ab, ac);

		if (winding == 0) {
			return;
		}

		auto is_top_left_edge = [](float v1x, float v1y, float v2x, float v2y, float v3x, float v3y, float winding) -> bool {
			if (winding < 0) {
				// top edge
				if (v1y == v2y && v1y > v3y) return true;
				// left edge
				if (v1y < v2y) return true;
			} 
			if (winding > 0) {
				// top edge
				if (v1y == v2y && v1y > v3y) return true;
				// left edge
				if (v1y > v2y) return true;
			}
            
            return false;
        };

		bool top_left_ab = is_top_left_edge(ax, ay, bx, by, cx, cy, winding);
        bool top_left_bc = is_top_left_edge(bx, by, cx, cy, ax, ay, winding);
        bool top_left_ca = is_top_left_edge(cx, cy, ax, ay, bx, by, winding);

		// from stackoverflow
		// check if c is between a and b
		auto isBetween = [](Vec2 a, Vec2 b, Vec2 c) {
			float crossproduct = (c.y - a.y) * (b.x - a.x) - (c.x - a.x) * (b.y - a.y);

			// compare versus epsilon for floating point values, or != 0 if using integers
			if (std::abs(crossproduct) > 0) {
				return false;
			}
				

			float dotproduct = (c.x - a.x) * (b.x - a.x) + (c.y - a.y)*(b.y - a.y);
			if (dotproduct < 0) {
				return false;
			}
				

			float squaredlengthba = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y);
			if (dotproduct > squaredlengthba) {
				return false;
			}
				

			return true;

		};
			

        float area = 0.5f * std::abs(cross_product(ab, ac));
        if (area == 0.0f) return;

        // Iterate over pixels in bounding box
        for (float y = minY; y <= maxY; ++y) {
            for (float x = minX; x <= maxX; ++x) {
				// pixel
				float px = x + 0.5f;
				float py = y + 0.5f;
				//std::cout << "px: " << px << ", py: " << py << ", "; 

				Vec2 pix (px, py);
				Vec2 a (ax, ay);
				Vec2 b (bx, by);
				Vec2 c (cx, cy);

				bool on_ab = isBetween(a, b, pix); // (py - ay) * (bx - ax) == (px - ax) * (by - ay);
				bool on_bc = isBetween(b, c, pix); //(py - by) * (cx - bx) == (px - bx) * (cy - by);
				bool on_ca = isBetween(c, a, pix); //(py - cy) * (ax - cx) == (px - cx) * (ay - cy);
                
				bool on_top_left = (on_ab && top_left_ab) || (on_bc && top_left_bc) || (on_ca && top_left_ca);
				bool on_non_top_left = (on_ab && !top_left_ab) || (on_bc && !top_left_bc) || (on_ca && !top_left_ca);

				if (on_non_top_left) {
					on_top_left = false;
				}

				//std::cout << std::boolalpha;
				//std::cout << "point on ab: " << on_ab << ", ";
				//std::cout << "point on bc: " << on_bc << ", ";
				//std::cout << "point on ca: " << on_ca << ", ";
				//std::cout << "point in triangle: " << point_in_triangle(px, py) << ", ";
				//std::cout << "point on top left: " << on_top_left << "\n";
                if (point_in_triangle(px, py) || on_top_left) {
                    // TODO: interpolate z
					Vec2 pb (bx - px, by - py);
					Vec2 pc (cx - px, cy - py);
					Vec2 pa (ax - px, ay - py);
					float phi_a = 0.5f * std::abs(cross_product(pb, pc)) / area;
					float phi_b = 0.5f * std::abs(cross_product(pa, pc)) / area;
					float phi_c = 0.5f * std::abs(cross_product(pa, pb)) / area;

                    float pz = phi_a * az + phi_b * bz + phi_c * cz;
					// float fakez = (va.fb_position.z + vb.fb_position.z + vc.fb_position.z) / 3.0f;

                    Fragment frag;
                    frag.fb_position = Vec3(px, py, pz);
                    frag.attributes = va.attributes;
                    frag.derivatives.fill(Vec2(0.0f, 0.0f));
                    emit_fragment(frag);
                }
            }
        }

	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Smooth) {
		// A1T5: screen-space smooth triangles
		// TODO: rasterize triangle (see block comment above this function).

		float ax = va.fb_position.x, ay = va.fb_position.y; 
        float bx = vb.fb_position.x, by = vb.fb_position.y; 
        float cx = vc.fb_position.x, cy = vc.fb_position.y; 

		float az = va.fb_position.z;
		float bz = vb.fb_position.z;
		float cz = vc.fb_position.z;

		 // bounding box
        float minX = std::floor(std::min({ax, bx, cx}));
        float minY = std::floor(std::min({ay, by, cy}));
        float maxX = std::ceil(std::max({ax, bx, cx}));
        float maxY = std::ceil(std::max({ay, by, cy}));

		auto cross_product = [](Vec2 v1, Vec2 v2) -> float {
        	return v1.x * v2.y - v1.y * v2.x;
    	};

		auto point_in_triangle = [&](float x , float y) -> bool {
            Vec2 ac (cx - ax, cy - ay);
            Vec2 ab (bx - ax, by - ay);
            Vec2 cb (bx - cx, by - cy);
            Vec2 ca (ax - cx, ay - cy);
            Vec2 ba (ax - bx, ay - by);
            Vec2 bc (cx - bx, cy - by);

            Vec2 aq (x - ax, y - ay);
            Vec2 cq (x - cx, y - cy);
            Vec2 bq (x - bx, y - by);

            float w1 = cross_product(ac, ab) * cross_product(ac, aq);
            float w2 = cross_product(cb, ca) * cross_product(cb, cq);
            float w3 = cross_product(ba, bc) * cross_product(ba, bq);

            return (w1 > 0 && w2 > 0 && w3 > 0) || (w1 < 0 && w2 < 0 && w3 < 0);
        };

		// winding > 0 = CCW, winding < 0 = CW
		Vec2 ab (bx - ax, by - ay);
		Vec2 ac (cx - ax, cy - ay);
		float winding = cross_product(ab, ac);

		if (winding == 0) {
			return;
		}

		auto is_top_left_edge = [](float v1x, float v1y, float v2x, float v2y, float v3x, float v3y, float winding) -> bool {
			if (winding < 0) {
				// top edge
				if (v1y == v2y && v1y > v3y) return true;
				// left edge
				if (v1y < v2y) return true;
			} 
			if (winding > 0) {
				// top edge
				if (v1y == v2y && v1y > v3y) return true;
				// left edge
				if (v1y > v2y) return true;
			}
            
            return false;
        };

		bool top_left_ab = is_top_left_edge(ax, ay, bx, by, cx, cy, winding);
        bool top_left_bc = is_top_left_edge(bx, by, cx, cy, ax, ay, winding);
        bool top_left_ca = is_top_left_edge(cx, cy, ax, ay, bx, by, winding);

		// from stackoverflow
		// check if c is between a and b
		auto isBetween = [](Vec2 a, Vec2 b, Vec2 c) {
			float crossproduct = (c.y - a.y) * (b.x - a.x) - (c.x - a.x) * (b.y - a.y);

			// compare versus epsilon for floating point values, or != 0 if using integers
			if (std::abs(crossproduct) > 0) {
				return false;
			}
				

			float dotproduct = (c.x - a.x) * (b.x - a.x) + (c.y - a.y)*(b.y - a.y);
			if (dotproduct < 0) {
				return false;
			}
				

			float squaredlengthba = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y);
			if (dotproduct > squaredlengthba) {
				return false;
			}
				

			return true;

		};
			

        float area = 0.5f * std::abs(cross_product(ab, ac));
        if (area == 0.0f) return;

        // iterate over pixels in bounding box
        for (float y = minY; y <= maxY; y++) {
            for (float x = minX; x <= maxX; x++) {

				// pixel
				float px = x + 0.5f;
				float py = y + 0.5f;
				//std::cout << "px: " << px << ", py: " << py << ", "; 

				Vec2 pix (px, py);
				Vec2 a (ax, ay);
				Vec2 b (bx, by);
				Vec2 c (cx, cy);

				bool on_ab = isBetween(a, b, pix); // (py - ay) * (bx - ax) == (px - ax) * (by - ay);
				bool on_bc = isBetween(b, c, pix); //(py - by) * (cx - bx) == (px - bx) * (cy - by);
				bool on_ca = isBetween(c, a, pix); //(py - cy) * (ax - cx) == (px - cx) * (ay - cy);
				
				bool on_top_left = (on_ab && top_left_ab) || (on_bc && top_left_bc) || (on_ca && top_left_ca);
				bool on_non_top_left = (on_ab && !top_left_ab) || (on_bc && !top_left_bc) || (on_ca && !top_left_ca);

				if (on_non_top_left) {
					on_top_left = false;
				}

				

				if (point_in_triangle(px, py) || on_top_left) {
					// if point in triangle
					Vec2 pb (bx - px, by - py);
					Vec2 pc (cx - px, cy - py);
					Vec2 pa (ax - px, ay - py);
					float phi_a = 0.5f * std::abs(cross_product(pb, pc)) / area;
					float phi_b = 0.5f * std::abs(cross_product(pa, pc)) / area;
					float phi_c = 0.5f * std::abs(cross_product(pa, pb)) / area;

					float pz = phi_a * az + phi_b * bz + phi_c * cz;

					// calculate attributes
					std::array<float, FA> frag_attributes;
					std::array<Vec2, FD> frag_derivatives;
					for (size_t i = 0; i < FA; i++) {
						frag_attributes[i] = phi_a * va.attributes[i] + phi_b * vb.attributes[i] + phi_c * vc.attributes[i];

					}

					for (size_t i = 0; i < FD; i++) {
						// derivative implementation attempt 2 (take derivative of barycentric formula)
						float dL0_dx = (by - cy) / (2 * area);
						float dL0_dy = (cx - bx) / (2 * area);

						float dL1_dx = (cy - ay) / (2 * area);
						float dL1_dy = (ax - cx) / (2 * area);

						float dL2_dx = (ay - by) / (2 * area);
						float dL2_dy = (bx - ax) / (2 * area);
						float A0 = va.attributes[i];
						float A1 = vb.attributes[i];
						float A2 = vc.attributes[i];

						float dA_dx = A0 * dL0_dx + A1 * dL1_dx + A2 * dL2_dx;
						float dA_dy = A0 * dL0_dy + A1 * dL1_dy + A2 * dL2_dy;
						frag_derivatives[i] = Vec2(dA_dx, dA_dy);

					}
					
					Fragment frag;
					frag.fb_position = Vec3(px, py, pz);
					frag.derivatives = frag_derivatives;
					frag.attributes = frag_attributes;
					emit_fragment(frag);

					// std::cout << "Fragment at (" 
					// << frag.fb_position.x << ", " 
					// << frag.fb_position.y << ", " 
					// << frag.fb_position.z << ")";

				}
				
            }
        }

		
	} else if constexpr ((flags & PipelineMask_Interp) == Pipeline_Interp_Correct) {
		// A1T5: perspective correct triangles
		// TODO: rasterize triangle (block comment above this function).

		float ax = va.fb_position.x, ay = va.fb_position.y; 
        float bx = vb.fb_position.x, by = vb.fb_position.y; 
        float cx = vc.fb_position.x, cy = vc.fb_position.y; 

		float az = va.fb_position.z;
		float bz = vb.fb_position.z;
		float cz = vc.fb_position.z;

		 // bounding box
        float minX = std::floor(std::min({ax, bx, cx}));
        float minY = std::floor(std::min({ay, by, cy}));
        float maxX = std::ceil(std::max({ax, bx, cx}));
        float maxY = std::ceil(std::max({ay, by, cy}));

		auto cross_product = [](Vec2 v1, Vec2 v2) -> float {
        	return v1.x * v2.y - v1.y * v2.x;
    	};

		auto point_in_triangle = [&](float x , float y) -> bool {
            Vec2 ac (cx - ax, cy - ay);
            Vec2 ab (bx - ax, by - ay);
            Vec2 cb (bx - cx, by - cy);
            Vec2 ca (ax - cx, ay - cy);
            Vec2 ba (ax - bx, ay - by);
            Vec2 bc (cx - bx, cy - by);

            Vec2 aq (x - ax, y - ay);
            Vec2 cq (x - cx, y - cy);
            Vec2 bq (x - bx, y - by);

            float w1 = cross_product(ac, ab) * cross_product(ac, aq);
            float w2 = cross_product(cb, ca) * cross_product(cb, cq);
            float w3 = cross_product(ba, bc) * cross_product(ba, bq);

            return (w1 > 0 && w2 > 0 && w3 > 0) || (w1 < 0 && w2 < 0 && w3 < 0);
        };

		// winding > 0 = CCW, winding < 0 = CW
		Vec2 ab (bx - ax, by - ay);
		Vec2 ac (cx - ax, cy - ay);
		float winding = cross_product(ab, ac);

		if (winding == 0) {
			return;
		}

		auto is_top_left_edge = [](float v1x, float v1y, float v2x, float v2y, float v3x, float v3y, float winding) -> bool {
			if (winding < 0) {
				// top edge
				if (v1y == v2y && v1y > v3y) return true;
				// left edge
				if (v1y < v2y) return true;
			} 
			if (winding > 0) {
				// top edge
				if (v1y == v2y && v1y > v3y) return true;
				// left edge
				if (v1y > v2y) return true;
			}
            
            return false;
        };

		bool top_left_ab = is_top_left_edge(ax, ay, bx, by, cx, cy, winding);
        bool top_left_bc = is_top_left_edge(bx, by, cx, cy, ax, ay, winding);
        bool top_left_ca = is_top_left_edge(cx, cy, ax, ay, bx, by, winding);

		// from stackoverflow
		// check if c is between a and b
		auto isBetween = [](Vec2 a, Vec2 b, Vec2 c) {
			float crossproduct = (c.y - a.y) * (b.x - a.x) - (c.x - a.x) * (b.y - a.y);

			// compare versus epsilon for floating point values, or != 0 if using integers
			if (std::abs(crossproduct) > 0) {
				return false;
			}
				

			float dotproduct = (c.x - a.x) * (b.x - a.x) + (c.y - a.y)*(b.y - a.y);
			if (dotproduct < 0) {
				return false;
			}
				

			float squaredlengthba = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y);
			if (dotproduct > squaredlengthba) {
				return false;
			}
				

			return true;

		};
			

        float area = 0.5f * std::abs(cross_product(ab, ac));
        if (area == 0.0f) return;

        // iterate over pixels in bounding box
        for (float y = minY; y <= maxY; y++) {
            for (float x = minX; x <= maxX; x++) {

				// pixel
				float px = x + 0.5f;
				float py = y + 0.5f;
				//std::cout << "px: " << px << ", py: " << py << ", "; 

				Vec2 pix (px, py);
				Vec2 a (ax, ay);
				Vec2 b (bx, by);
				Vec2 c (cx, cy);

				bool on_ab = isBetween(a, b, pix); // (py - ay) * (bx - ax) == (px - ax) * (by - ay);
				bool on_bc = isBetween(b, c, pix); //(py - by) * (cx - bx) == (px - bx) * (cy - by);
				bool on_ca = isBetween(c, a, pix); //(py - cy) * (ax - cx) == (px - cx) * (ay - cy);
				
				bool on_top_left = (on_ab && top_left_ab) || (on_bc && top_left_bc) || (on_ca && top_left_ca);
				bool on_non_top_left = (on_ab && !top_left_ab) || (on_bc && !top_left_bc) || (on_ca && !top_left_ca);

				if (on_non_top_left) {
					on_top_left = false;
				}

				

				if (point_in_triangle(px, py) || on_top_left) {
					// if point in triangle
					Vec2 pb (bx - px, by - py);
					Vec2 pc (cx - px, cy - py);
					Vec2 pa (ax - px, ay - py);
					float phi_a = 0.5f * std::abs(cross_product(pb, pc)) / area;
					float phi_b = 0.5f * std::abs(cross_product(pa, pc)) / area;
					float phi_c = 0.5f * std::abs(cross_product(pa, pb)) / area;

					float pz = phi_a * az + phi_b * bz + phi_c * cz;

					// interpolate inverse 4th component
					float w_a = va.inv_w;
					float w_b = vb.inv_w;
					float w_c = vc.inv_w;

					float inter_inv_w = phi_a * w_a + phi_b * w_b + phi_c * w_c;

					// calculate attributes
					std::array<float, FA> frag_attributes;
					std::array<Vec2, FD> frag_derivatives;
					for (size_t i = 0; i < FA; i++) {
						// interpolate phi / w
						float inter_phi_w = phi_a * va.attributes[i] * w_a + phi_b * vb.attributes[i] * w_b + phi_c * vc.attributes[i] * w_c;
						frag_attributes[i] = inter_phi_w / inter_inv_w;

					}

					for (size_t i = 0; i < FD; i++) {
						// derivative implementation attempt 2 (take derivative of barycentric formula)
						float dL0_dx = (by - cy) / (2 * area);
						float dL0_dy = (cx - bx) / (2 * area);

						float dL1_dx = (cy - ay) / (2 * area);
						float dL1_dy = (ax - cx) / (2 * area);

						float dL2_dx = (ay - by) / (2 * area);
						float dL2_dy = (bx - ax) / (2 * area);
						float A0 = va.attributes[i] * w_a;
						float A1 = vb.attributes[i] * w_b;
						float A2 = vc.attributes[i] * w_c;

						float dA_dx = (A0 * dL0_dx + A1 * dL1_dx + A2 * dL2_dx) / inter_inv_w;
						float dA_dy = (A0 * dL0_dy + A1 * dL1_dy + A2 * dL2_dy) / inter_inv_w;
						frag_derivatives[i] = Vec2(dA_dx, dA_dy);

					}
					
					Fragment frag;
					frag.fb_position = Vec3(px, py, pz);	
					frag.derivatives = frag_derivatives;
					frag.attributes = frag_attributes;
					emit_fragment(frag);
				}
				
            }
        }
	}
}

//-------------------------------------------------------------------------
// compile instantiations for all programs and blending and testing types:

#include "programs.h"

template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Smooth>;
template struct Pipeline<PrimitiveType::Triangles, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Correct>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Add | Pipeline_Depth_Less | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Always | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Never | Pipeline_Interp_Flat>;
template struct Pipeline<PrimitiveType::Lines, Programs::Lambertian,
                         Pipeline_Blend_Over | Pipeline_Depth_Less | Pipeline_Interp_Flat>;