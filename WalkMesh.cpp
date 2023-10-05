#include "WalkMesh.hpp"

#include "read_write_chunk.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

WalkMesh::WalkMesh(std::vector< glm::vec3 > const &vertices_, std::vector< glm::vec3 > const &normals_, std::vector< glm::uvec3 > const &triangles_)
	: vertices(vertices_), normals(normals_), triangles(triangles_) {

	//construct next_vertex map (maps each edge to the next vertex in the triangle):
	next_vertex.reserve(triangles.size()*3);
	auto do_next = [this](uint32_t a, uint32_t b, uint32_t c) {
		auto ret = next_vertex.insert(std::make_pair(glm::uvec2(a,b), c));
		assert(ret.second);
	};
	for (auto const &tri : triangles) {
		do_next(tri.x, tri.y, tri.z);
		do_next(tri.y, tri.z, tri.x);
		do_next(tri.z, tri.x, tri.y);
	}
}

// project pt to the plane of triangle a,b,c and return the barycentric weights of the projected point:
// Credit goes to graphics class from last year for inspiration
glm::vec3 barycentric_weights(glm::vec3 const &a, glm::vec3 const &b, glm::vec3 const &c, glm::vec3 const &pt) {
    // make a the origin
    glm::vec3 pa = pt - a;
    glm::vec3 ba = b - a;
    glm::vec3 ca = c - a;

    // cross product is area - get fraction of cross products to get fraction of area
    glm::vec3 areaVec = glm::cross(ba, ca);
    glm::vec3 areaNorm = glm::normalize(areaVec);
    float area = glm::length(areaVec);
    glm::vec3 pb = glm::cross(ba, pa) / area;
    glm::vec3 pc = glm::cross(pa, ca) / area;

    // projection to 'divide' vectors
    float x = glm::dot(pb, areaNorm);
    float y = glm::dot(pc, areaNorm);

    return glm::vec3(1.0f - (x + y), y, x);
}

WalkPoint WalkMesh::nearest_walk_point(glm::vec3 const &world_point) const {
	assert(!triangles.empty() && "Cannot start on an empty walkmesh");

	WalkPoint closest;
	float closest_dis2 = std::numeric_limits< float >::infinity();

	for (auto const &tri : triangles) {
		//find closest point on triangle:

		glm::vec3 const &a = vertices[tri.x];
		glm::vec3 const &b = vertices[tri.y];
		glm::vec3 const &c = vertices[tri.z];

		//get barycentric coordinates of closest point in the plane of (a,b,c):
		glm::vec3 coords = barycentric_weights(a,b,c, world_point);

		//is that point inside the triangle?
		if (coords.x >= 0.0f && coords.y >= 0.0f && coords.z >= 0.0f) {
			//yes, point is inside triangle.
			float dis2 = glm::length2(world_point - to_world_point(WalkPoint(tri, coords)));
			if (dis2 < closest_dis2) {
				closest_dis2 = dis2;
				closest.indices = tri;
				closest.weights = coords;
			}
		} else {
			//check triangle vertices and edges:
			auto check_edge = [&world_point, &closest, &closest_dis2, this](uint32_t ai, uint32_t bi, uint32_t ci) {
				glm::vec3 const &a = vertices[ai];
				glm::vec3 const &b = vertices[bi];

				//find closest point on line segment ab:
				float along = glm::dot(world_point-a, b-a);
				float max = glm::dot(b-a, b-a);
				glm::vec3 pt;
				glm::vec3 coords;
				if (along < 0.0f) {
					pt = a;
					coords = glm::vec3(1.0f, 0.0f, 0.0f);
				} else if (along > max) {
					pt = b;
					coords = glm::vec3(0.0f, 1.0f, 0.0f);
				} else {
					float amt = along / max;
					pt = glm::mix(a, b, amt);
					coords = glm::vec3(1.0f - amt, amt, 0.0f);
				}

				float dis2 = glm::length2(world_point - pt);
				if (dis2 < closest_dis2) {
					closest_dis2 = dis2;
					closest.indices = glm::uvec3(ai, bi, ci);
					closest.weights = coords;
				}
			};
			check_edge(tri.x, tri.y, tri.z);
			check_edge(tri.y, tri.z, tri.x);
			check_edge(tri.z, tri.x, tri.y);
		}
	}
	assert(closest.indices.x < vertices.size());
	assert(closest.indices.y < vertices.size());
	assert(closest.indices.z < vertices.size());
	return closest;
}

// Adapted from @stroucki's code as shown in class
void WalkMesh::walk_in_triangle(WalkPoint const &start, glm::vec3 const &step, WalkPoint *end_, float *time_) const {
	assert(end_);
	auto &end = *end_;

	assert(time_);
	auto &time = *time_;

	glm::vec3 step_coords;
	{
        glm::vec3 start_world = to_world_point(start);
        glm::vec3 end_world = start_world + step;
        //project 'step' into a barycentric-coordinates direction:
        step_coords = barycentric_weights(vertices[start.indices[0]], vertices[start.indices[1]], vertices[start.indices[2]], end_world) - start.weights;
    }

    //figure out which edge (if any) is crossed first.
    // set time and end appropriately.
    glm::vec3 intersection = -start.weights / step_coords;

    if (step_coords.x >= 0.0f)
        intersection.x = std::numeric_limits<float>::infinity();
    if (step_coords.y >= 0.0f)
        intersection.y = std::numeric_limits<float>::infinity();
    if (step_coords.z >= 0.0f)
        intersection.z = std::numeric_limits<float>::infinity();

    int min;
    if (intersection.x <= intersection.y && intersection.x <= intersection.z) {
        min = 0;
        time = intersection.x;
    } else if (intersection.y <= intersection.z && intersection.y <= intersection.x) {
        min = 1;
        time = intersection.y;
    } else {
        min = 2;
        time = intersection.z;
    }

    if (time >= 1.0f) {
        time = 1.0f;
        min = -1;
    }

    //printf("TIME = %f\n", time);

    glm::vec3 end_weights = start.weights + step_coords * time;

	//Remember: our convention is that when a WalkPoint is on an edge,
	// then wp.weights.z == 0.0f (so will likely need to re-order the indices)
    if (min == 0) {
        end.indices[0] = start.indices[1];
        end.indices[1] = start.indices[2];
        end.indices[2] = start.indices[0];
        end.weights[0] = end_weights[1];
        end.weights[1] = end_weights[2];
        end.weights[2] = 0.0f;
    } else if (min == 1) {
        end.indices[0] = start.indices[2];
        end.indices[1] = start.indices[0];
        end.indices[2] = start.indices[1];
        end.weights[0] = end_weights[2];
        end.weights[1] = end_weights[0];
        end.weights[2] = 0.0f;
    } else if (min == 2) {
        end.indices = start.indices;
        end.weights = end_weights;
        end.weights[2] = 0.0f;
    } else {
        end.indices = start.indices;
        end.weights = end_weights;
    }
}

bool WalkMesh::cross_edge(WalkPoint const &start, WalkPoint *end_, glm::quat *rotation_) const {
	assert(end_);
	auto &end = *end_;

	assert(rotation_);
	auto &rotation = *rotation_;

	assert(start.weights.z == 0.0f); //*must* be on an edge.
	glm::uvec2 edge = glm::uvec2(start.indices);

    auto twin = next_vertex.find(glm::uvec2(edge.y, edge.x));
	//check if 'edge' is a non-boundary edge:
	if (twin != next_vertex.end()) {
		//it is!

        unsigned int third_index = twin->second;
        glm::vec3 third = vertices[third_index];
		//make 'end' represent the same (world) point, but on triangle (edge.y, edge.x, [other point]):
        end.indices.x = edge.y;
        end.indices.y = edge.x;
        end.indices.z = third_index;
        end.weights = barycentric_weights(vertices[edge.y], vertices[edge.x], third, to_world_point(start));

		//make 'rotation' the rotation that takes (start.indices)'s normal to (end.indices)'s normal:
        glm::vec3 x1 = vertices[start.indices.x];
        glm::vec3 y1 = vertices[start.indices.y];
        glm::vec3 z1 = vertices[start.indices.z];
        glm::vec3 n1 = glm::normalize(glm::cross(z1 - x1, y1 - x1));

        glm::vec3 x2 = vertices[end.indices.x];
        glm::vec3 y2 = vertices[end.indices.y];
        glm::vec3 z2 = vertices[end.indices.z];
        glm::vec3 n2 = glm::normalize(glm::cross(z2 - x2, y2 - x2));

        glm::vec3 c = glm::cross(n1, n2);
        if (glm::length(c) > 0.00001f)
            rotation = glm::normalize(glm::angleAxis(asin(glm::length(c)), glm::normalize(c)));
        else
            rotation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);

		return true;
	} else {
		end = start;
		rotation = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
		return false;
	}
}


WalkMeshes::WalkMeshes(std::string const &filename) {
	std::ifstream file(filename, std::ios::binary);

	std::vector< glm::vec3 > vertices;
	read_chunk(file, "p...", &vertices);

	std::vector< glm::vec3 > normals;
	read_chunk(file, "n...", &normals);

	std::vector< glm::uvec3 > triangles;
	read_chunk(file, "tri0", &triangles);

	std::vector< char > names;
	read_chunk(file, "str0", &names);

	struct IndexEntry {
		uint32_t name_begin, name_end;
		uint32_t vertex_begin, vertex_end;
		uint32_t triangle_begin, triangle_end;
	};

	std::vector< IndexEntry > index;
	read_chunk(file, "idxA", &index);

	if (file.peek() != EOF) {
		std::cerr << "WARNING: trailing data in walkmesh file '" << filename << "'" << std::endl;
	}

	//-----------------

	if (vertices.size() != normals.size()) {
		throw std::runtime_error("Mis-matched position and normal sizes in '" + filename + "'");
	}

	for (auto const &e : index) {
		if (!(e.name_begin <= e.name_end && e.name_end <= names.size())) {
			throw std::runtime_error("Invalid name indices in index of '" + filename + "'");
		}
		if (!(e.vertex_begin <= e.vertex_end && e.vertex_end <= vertices.size())) {
			throw std::runtime_error("Invalid vertex indices in index of '" + filename + "'");
		}
		if (!(e.triangle_begin <= e.triangle_end && e.triangle_end <= triangles.size())) {
			throw std::runtime_error("Invalid triangle indices in index of '" + filename + "'");
		}

		//copy vertices/normals:
		std::vector< glm::vec3 > wm_vertices(vertices.begin() + e.vertex_begin, vertices.begin() + e.vertex_end);
		std::vector< glm::vec3 > wm_normals(normals.begin() + e.vertex_begin, normals.begin() + e.vertex_end);

		//remap triangles:
		std::vector< glm::uvec3 > wm_triangles; wm_triangles.reserve(e.triangle_end - e.triangle_begin);
		for (uint32_t ti = e.triangle_begin; ti != e.triangle_end; ++ti) {
			if (!( (e.vertex_begin <= triangles[ti].x && triangles[ti].x < e.vertex_end)
			    && (e.vertex_begin <= triangles[ti].y && triangles[ti].y < e.vertex_end)
			    && (e.vertex_begin <= triangles[ti].z && triangles[ti].z < e.vertex_end) )) {
				throw std::runtime_error("Invalid triangle in '" + filename + "'");
			}
			wm_triangles.emplace_back(
				triangles[ti].x - e.vertex_begin,
				triangles[ti].y - e.vertex_begin,
				triangles[ti].z - e.vertex_begin
			);
		}
		
		std::string name(names.begin() + e.name_begin, names.begin() + e.name_end);

		auto ret = meshes.emplace(name, WalkMesh(wm_vertices, wm_normals, wm_triangles));
		if (!ret.second) {
			throw std::runtime_error("WalkMesh with duplicated name '" + name + "' in '" + filename + "'");
		}

	}
}

WalkMesh const &WalkMeshes::lookup(std::string const &name) const {
	auto f = meshes.find(name);
	if (f == meshes.end()) {
		throw std::runtime_error("WalkMesh with name '" + name + "' not found.");
	}
	return f->second;
}
