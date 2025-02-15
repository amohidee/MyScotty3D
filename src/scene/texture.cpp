
#include "texture.h"

#include <iostream>

namespace Textures {


Spectrum sample_nearest(HDR_Image const &image, Vec2 uv) {
	//clamp texture coordinates, convert to [0,w]x[0,h] pixel space:
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

	//the pixel with the nearest center is the pixel that contains (x,y):
	int32_t ix = int32_t(std::floor(x));
	int32_t iy = int32_t(std::floor(y));

	//texture coordinates of (1,1) map to (w,h), and need to be reduced:
	ix = std::min(ix, int32_t(image.w) - 1);
	iy = std::min(iy, int32_t(image.h) - 1);

	return image.at(ix, iy);
}

Spectrum sample_bilinear(HDR_Image const &image, Vec2 uv) {
	// A1T6: sample_bilinear
	//TODO: implement bilinear sampling strategy on texture 'image'

	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
    float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

    int x0 = std::clamp(int32_t(std::floor(x - 0.5)), 0, int32_t(image.w) - 1);
    int y0 = std::clamp(int32_t(std::floor(y - 0.5)), 0, int32_t(image.h) - 1);
    int x1 = std::clamp(x0 + 1, 0, int32_t(image.w) - 1);
    int y1 = std::clamp(y0 + 1, 0, int32_t(image.h) - 1);

    float dx = std::max(x - (x0 + 0.5F), 0.0f);
    float dy = std::max(y - (y0 + 0.5F), 0.0f);


	//std::cout << "(" << x0 << ", " << y0 << ", " << x1 << ", " << y1 << ")\n";
	if (x1 > 1999 || y1 > 500 || x0 > 1999 || y0 > 500) {
		std::cout << "(" << x0 << ", " << y0 << ", " << x1 << ", " << y1 << ")\n";
	}

    // 4 nearest texels
    Spectrum t00 = image.at(x0, y0); 
    Spectrum t10 = image.at(x1, y0); 
    Spectrum t01 = image.at(x0, y1);
    Spectrum t11 = image.at(x1, y1);

    // interpolation
    Spectrum tx = (1.0f - dx) * t00 + dx * t10; 
    Spectrum ty = (1.0f - dx) * t01 + dx * t11;
    Spectrum t  = (1.0f - dy) * tx + dy * ty;

    return t;
}


Spectrum sample_trilinear(HDR_Image const &base, std::vector< HDR_Image > const &levels, Vec2 uv, float lod) {
	// A1T6: sample_trilinear
	//TODO: implement trilinear sampling strategy on using mip-map 'levels'

	size_t max_level = levels.size();
	
    float clamped_lod = std::clamp(lod, 0.0f, float(max_level));

    int d = int(std::floor(clamped_lod));

    int d1 = std::min(d + 1, int(max_level));

    float dd = clamped_lod - d;


    Spectrum td = (d == 0) ? sample_bilinear(base, uv) : sample_bilinear(levels[d - 1], uv);
    Spectrum td1 = (d1 == 0) ? sample_bilinear(base, uv): sample_bilinear(levels[d1 - 1], uv);

	auto output = (1.0f - dd) * td + dd * td1;

    return output;
}

/*
 * generate_mipmap- generate mipmap levels from a base image.
 *  base: the base image
 *  levels: pointer to vector of levels to fill (must not be null)
 *
 * generates a stack of levels [1,n] of sizes w_i, h_i, where:
 *   w_i = max(1, floor(w_{i-1})/2)
 *   h_i = max(1, floor(h_{i-1})/2)
 *  with:
 *   w_0 = base.w
 *   h_0 = base.h
 *  and n is the smalles n such that w_n = h_n = 1
 *
 * each level should be calculated by downsampling a blurred version
 * of the previous level to remove high-frequency detail.
 *
 */
void generate_mipmap(HDR_Image const &base, std::vector< HDR_Image > *levels_) {
	assert(levels_);
	auto &levels = *levels_;


	{ // allocate sublevels sufficient to scale base image all the way to 1x1:
		int32_t num_levels = static_cast<int32_t>(std::log2(std::max(base.w, base.h)));
		assert(num_levels >= 0);

		levels.clear();
		levels.reserve(num_levels);

		uint32_t width = base.w;
		uint32_t height = base.h;
		for (int32_t i = 0; i < num_levels; ++i) {
			assert(!(width == 1 && height == 1)); //would have stopped before this if num_levels was computed correctly

			width = std::max(1u, width / 2u);
			height = std::max(1u, height / 2u);

			levels.emplace_back(width, height);
		}
		assert(width == 1 && height == 1);
		assert(levels.size() == uint32_t(num_levels));
	}

	//now fill in the levels using a helper:
	//downsample:
	// fill in dst to represent the low-frequency component of src
	auto downsample = [](HDR_Image const &src, HDR_Image &dst) {
		//dst is half the size of src in each dimension:
		assert(std::max(1u, src.w / 2u) == dst.w);
		assert(std::max(1u, src.h / 2u) == dst.h);

		// A1T6: generate
		//TODO: Write code to fill the levels of the mipmap hierarchy by downsampling

		//Be aware that the alignment of the samples in dst and src will be different depending on whether the image is even or odd.
		for (uint32_t y = 0; y < dst.h; ++y) {
            for (uint32_t x = 0; x < dst.w; ++x) {
                // Compute corresponding source coordinates
                uint32_t src_x = x * 2;
                uint32_t src_y = y * 2;

				if (src_x + 1 > 1999 || src_y + 1 > 999 || src_x > 1999 || src_y > 999) {
					std::cout << "(" << src_x << ", " << src_y << ")\n";
				}

                // Gather pixel samples (handling odd sizes)
                Spectrum t00 = src.at(src_x, src_y);
                Spectrum t10 = (src_x + 1 < src.w) ? src.at(src_x + 1, src_y) : t00;
                Spectrum t01 = (src_y + 1 < src.h) ? src.at(src_x, src_y + 1) : t00;
                Spectrum t11 = (src_x + 1 < src.w && src_y + 1 < src.h) ? src.at(src_x + 1, src_y + 1) : t00;

                // Compute downsampled pixel (box filter average)
                dst.at(x, y) = (t00 + t10 + t01 + t11) / 4.0f;
            }
        }

	};

	std::cout << "Regenerating mipmap (" << levels.size() << " levels): [" << base.w << "x" << base.h << "]";
	std::cout.flush();
	for (uint32_t i = 0; i < levels.size(); ++i) {
		HDR_Image const &src = (i == 0 ? base : levels[i-1]);
		HDR_Image &dst = levels[i];
		std::cout << " -> [" << dst.w << "x" << dst.h << "]"; std::cout.flush();

		downsample(src, dst);
	}
	std::cout << std::endl;
	
}

Image::Image(Sampler sampler_, HDR_Image const &image_) {
	sampler = sampler_;
	image = image_.copy();
	update_mipmap();
}

Spectrum Image::evaluate(Vec2 uv, float lod) const {
	if (image.w == 0 && image.h == 0) return Spectrum();
	if (sampler == Sampler::nearest) {
		return sample_nearest(image, uv);
	} else if (sampler == Sampler::bilinear) {
		return sample_bilinear(image, uv);
	} else {
		return sample_trilinear(image, levels, uv, lod);
	}
}

void Image::update_mipmap() {
	if (sampler == Sampler::trilinear) {
		generate_mipmap(image, &levels);
	} else {
		levels.clear();
	}
}

GL::Tex2D Image::to_gl() const {
	return image.to_gl(1.0f);
}

void Image::make_valid() {
	update_mipmap();
}

Spectrum Constant::evaluate(Vec2 uv, float lod) const {
	return color * scale;
}

} // namespace Textures
bool operator!=(const Textures::Constant& a, const Textures::Constant& b) {
	return a.color != b.color || a.scale != b.scale;
}

bool operator!=(const Textures::Image& a, const Textures::Image& b) {
	return a.image != b.image;
}

bool operator!=(const Texture& a, const Texture& b) {
	if (a.texture.index() != b.texture.index()) return false;
	return std::visit(
		[&](const auto& data) { return data != std::get<std::decay_t<decltype(data)>>(b.texture); },
		a.texture);
}
