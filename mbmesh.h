#ifndef MBMSH_H
#define MBMSH_H

#include <cstring>
#include <cstdint>
#include <array>
#include <vector>
#include <stdexcept>
#include <istream>
#include <ostream>
#include <iomanip>
#include <sstream>

namespace mbmesh {

const unsigned int NUM_BONES = 64;
const uint64_t MESH_HEADER_MAGIC = 0x30306873656d626d;
const uint64_t ANIMATION_HEADER_MAGIC = 0x30306d696e61626d;

class Mesh {
public:
    struct Vertex {
        struct Position {
            float x,y,z;
        } pos;
        struct Normal {
            float x,y,z;
        } norm;
        struct TexCoord {
            float u,v;
        } tex;
        typedef float Weight;
        typedef std::array<Weight, 4> WeightList;
        typedef uint32_t BoneIndex;
        typedef std::array<BoneIndex, 4> BoneList;
        WeightList weights;
        BoneList bones;
        bool operator<(const Vertex &rhs) const {
            return memcmp(this, &rhs, sizeof(*this)) < 0;
        }
        void dump_text(std::ostream& out) const;
    };

    struct Triangle {
        typedef uint32_t ElementIndex;
        ElementIndex a,b,c;
    };

    struct Bone {
        Vertex::BoneIndex parent;
        std::array<float, 16> mat;
    };

    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
    std::array<Bone, NUM_BONES> bones;

    struct Header {
        uint64_t magic;
        uint32_t num_vertices;
        uint32_t num_triangles;
    };

    Mesh();
    Mesh(const std::string &filename);
    void write(const std::string &filename) const;
    void dump_text(std::ostream &) const;
};

struct Orientation {
    struct Location {
        float x,y,z;
    } loc;
    struct Rotation {
        float w,x,y,z;
    } rot;
    struct Scale {
        float x,y,z;
    } scale;
    Orientation() :
        loc {0.0f, 0.0f, 0.0f},
        rot {1.0f,0.0f,0.0f,0.0f},
    scale {1.0f,1.0f, 1.0f} {
    }
};

struct RigOrientation {
    std::array<Orientation, NUM_BONES> bone_orientations;
    Orientation operator[](int i) const {
        return bone_orientations[i];
    }
};

class Animation {
public:
    Animation();
    Animation(const std::string& filename);
    RigOrientation at(const float time) const;
    void write(const std::string& filename) const;
    void dump_text(std::ostream &out) const;

private:
    struct Curve {
        uint32_t start_index;
        uint32_t end_index;
    };
    struct BoneAnim {
        struct Location {
            Curve x;
            Curve y;
            Curve z;
        } loc;
        struct Rotation {
            Curve w;
            Curve x;
            Curve y;
            Curve z;
        } rot;
        struct Scale {
            Curve x;
            Curve y;
            Curve z;
        } scale;
    };

    struct Header {
        uint64_t magic;
        float duration;
        uint32_t data_count;
    };

    float duration;
    std::array<BoneAnim, NUM_BONES> bone_anims;
    std::vector<float> data;

};

class Pose {
public:
    std::array< std::array<float, 16>, NUM_BONES > matrices;
    Pose() {
        for (auto &m : matrices) {
            m = {{1.f, 0, 0, 0, 0, 1.f, 0, 0, 0, 0, 1.f, 0, 0, 0, 0, 1.f }};
        }
    }
    Pose(const Mesh &mesh, const RigOrientation &rig);
    operator float*() {
        return matrices[0].data();
    }

};

}

#endif /* MBMSH_H */
