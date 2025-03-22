
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <stack>
#include <iostream>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims; ///< number of primitives in the bucket
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh

	// Keep these
    nodes.clear();
    primitives = std::move(prims);

    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.
    if (primitives.empty()) return;
	//TODO
    root_idx = build_recursive(0, primitives.size(), max_leaf_size);
}

template<typename Primitive>
size_t BVH<Primitive>::build_recursive(size_t start, size_t end, size_t max_leaf_size) {
    size_t count = end - start;
    // std::cout << "(start, end, count) = (" << start << ", " << end << ", " << count << ")" << std::endl;

    // compute sub bbox for recursive call
    BBox bounds;
    BBox centroid_bounds;

    for (size_t i = start; i < end; i++) {
        bounds.enclose(primitives[i].bbox());
    }
    for (size_t i = start; i < end; i++) {
        centroid_bounds.enclose(primitives[i].bbox().center());
    }

    

    
    struct Bucket {
        BBox bbox;
        size_t count = 0;
    };

    const int num_buckets = 12;
    float best_cost = std::numeric_limits<float>::infinity();
    int best_axis = -1;
    int best_split = -1;


    // base case: create leaf
    if (count <= max_leaf_size) {
        return new_node(bounds, start, count, 0, 0);
    }

    // else recursive case
    for (int axis = 0; axis < 3; axis++) {
        Bucket buckets[num_buckets];

        if (centroid_bounds.max[axis] - centroid_bounds.min[axis] < 1e-8f) {
            continue;
        }

        // compute_bucket(p.centroid)
        for (size_t i = start; i < end; i++) {
            Vec3 centroid = primitives[i].bbox().center();
            float relative_pos = (centroid[axis] - centroid_bounds.min[axis]) / (centroid_bounds.max[axis] - centroid_bounds.min[axis]);
            int bucket_num = std::min(int(relative_pos * num_buckets), num_buckets - 1);
            buckets[bucket_num].bbox.enclose(primitives[i].bbox());
            buckets[bucket_num].count++;
        }

        // for each |B| - 1 possible partitions
        for (int split = 1; split < num_buckets; split++) {
            BBox left_bbox;
			BBox right_bbox;
            size_t left_count = 0;
			size_t right_count = 0;

            // split into left and right bbox

            for (int i = 0; i < split; i++) {
                left_bbox.enclose(buckets[i].bbox);
                left_count += buckets[i].count;
            }

            for (int i = split; i < num_buckets; i++) {
                right_bbox.enclose(buckets[i].bbox);
                right_count += buckets[i].count;
            }

            if (left_count == 0 || right_count == 0) continue;

            // calculate surface area using SAH

            float SA_left = left_bbox.surface_area();
            float SA_right = right_bbox.surface_area();
            float SA_total = bounds.surface_area();

            float cost = 0.125f + (SA_left / SA_total) * left_count + 
                         (SA_right / SA_total) * right_count;

            if (cost < best_cost) {
                best_cost = cost;
                best_axis = axis;
                best_split = split;
            }
        }
    }

    // recurse on lowest cost partition found (or make node leaf)
    if (best_axis == -1) {
        // std::cout << "best axis = -1" << std::endl;
        return new_node(bounds, start, count, 0, 0);
    }

    float relative_split = (centroid_bounds.max[best_axis] - centroid_bounds.min[best_axis]) * (float(best_split) / num_buckets);

    float split_val = centroid_bounds.min[best_axis] + relative_split;

    // std::partition call
    auto mid = std::partition(primitives.begin() + start, primitives.begin() + end,
        [best_axis, split_val](const Primitive& p) { return p.bbox().center()[best_axis] < split_val; });

    size_t mid_idx = mid - primitives.begin();

    // edge case: mid_idx == start or end???

    size_t left = build_recursive(start, mid_idx, max_leaf_size);
    size_t right = build_recursive(mid_idx, end, max_leaf_size);

    return new_node(bounds, start, count, left, right);
}

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

	//TODO: replace this code with a more efficient traversal:
    Trace ret;
    if (nodes.empty()) {
        return ret;
    }
    find_closest_hit(ray, &nodes[root_idx], &ret);
    return ret;
}

template<typename Primitive>
void BVH<Primitive>::find_closest_hit(const Ray& ray, const Node* node, Trace* closest) const {

    // first check if ray even hits bounding box
    Vec2 t_bounds = ray.dist_bounds;
    if (!node->bbox.hit(ray, t_bounds)) return;

    if (node->is_leaf()) {
        // leaf case
        for (size_t i = node->start; i < node->start + node->size; i++) {
            Trace hit = primitives[i].hit(ray);
            *closest = Trace::min(*closest, hit);
        }
    } else {
        // internal node
        const Node* left = &nodes[node->l];
        const Node* right = &nodes[node->r];

        Vec2 t_left = ray.dist_bounds;
        Vec2 t_right = ray.dist_bounds;

        bool hit_left = left->bbox.hit(ray, t_left);
        bool hit_right = right->bbox.hit(ray, t_right);

        if (hit_left && hit_right) {
            bool left_first = t_left.x <= t_right.x;

            const Node* first = left_first ? left : right;
            const Node* second = left_first ? right : left;
            float t_second = left_first ? t_right.x : t_left.x;

            find_closest_hit(ray, first, closest);

            if (!closest->hit || t_second < closest->distance) {
                find_closest_hit(ray, second, closest);
            }
        } else if (hit_left) {
            find_closest_hit(ray, left, closest);
        } else if (hit_right) {
            find_closest_hit(ray, right, closest);
        }
    }
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
