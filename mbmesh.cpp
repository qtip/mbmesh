#include "mbmesh.h"
#include <cmath>
#include <stdio.h>
#include <ostream>
#include <vector>
using std::vector;
#include <algorithm>
using std::rotate;
using std::upper_bound;
#include <array>
using std::array;
using std::get;
#include <iterator>
using std::advance;
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <memory>
using std::unique_ptr;
#include <physfs.h>

using namespace mbmesh;

/*===== Mesh ======*/

void Mesh::Vertex::dump_text(std::ostream &out) const {
    out<<std::fixed;
    out<<"  pos:     ["<<pos.x<<","<<pos.y<<","<<pos.z<<"]\n";
    out<<"  norm:    ["<<norm.x<<","<<norm.y<<","<<norm.z<<"]\n";
    out<<"  tex:     ["<<tex.u<<","<<tex.v<<"]\n";
    out<<"  weights: ["<<weights[0]<<","<<weights[1]<<","<<weights[2]<<","<<weights[3]<<"]\n";
    out<<"  bones:   ["<<bones[0]<<","<<bones[1]<<","<<bones[2]<<","<<bones[3]<<"]\n";
}

Mesh::Mesh() {}

Mesh::Mesh(const string &filename) {
    unique_ptr<PHYSFS_File, int (&)(PHYSFS_File*)> physfp(PHYSFS_openRead(filename.c_str()), PHYSFS_close);
    unique_ptr<FILE, int (&)(FILE*)> fp(fopen(filename.c_str(), "rb"), fclose);

    if (!fp.get() && !physfp.get()) {
        std::stringstream what;
        what<<"couldn't load file "<<filename;
        throw std::runtime_error(what.str());
    }

    Header header;
    if (physfp.get()) {
        PHYSFS_read(physfp.get(), reinterpret_cast<char*>(&header), 1, sizeof(header));
    } else {
        fread(reinterpret_cast<char*>(&header), sizeof(header), 1, fp.get());
    }
    if (header.magic != MESH_HEADER_MAGIC) {
        std::stringstream what;
        what<<"invalid format, got 0x"<<std::hex<<header.magic<<", expected 0x"<<MESH_HEADER_MAGIC;
        what<<" (";
        what.write(reinterpret_cast<const char*>(&MESH_HEADER_MAGIC), sizeof MESH_HEADER_MAGIC);
        what<<")";
        throw std::runtime_error(what.str());
    }
    vertices.resize(header.num_vertices);
    triangles.resize(header.num_triangles);
    if (physfp.get()) {
        PHYSFS_read(physfp.get(), reinterpret_cast<char*>(vertices.data()), sizeof(Vertex), header.num_vertices);
        PHYSFS_read(physfp.get(), reinterpret_cast<char*>(triangles.data()), sizeof(Triangle), header.num_triangles);
        PHYSFS_read(physfp.get(), reinterpret_cast<char*>(bones.data()), sizeof(Bone), NUM_BONES);
    } else {
        fread(reinterpret_cast<char*>(vertices.data()), header.num_vertices, sizeof(Vertex), fp.get());
        fread(reinterpret_cast<char*>(triangles.data()), header.num_triangles, sizeof(Triangle), fp.get());
        fread(reinterpret_cast<char*>(bones.data()), NUM_BONES, sizeof(Bone), fp.get());
    }
}

void Mesh::write(const std::string &filename) const {
    Header header;
    ofstream stream(filename, ofstream::binary);
    header.magic = MESH_HEADER_MAGIC;
    header.num_vertices = vertices.size();
    header.num_triangles = triangles.size();
    stream.write(reinterpret_cast<const char*>(&header), sizeof header);
    stream.write(reinterpret_cast<const char*>(vertices.data()), sizeof(Vertex) * header.num_vertices);
    stream.write(reinterpret_cast<const char*>(triangles.data()), sizeof(Triangle) * header.num_triangles);
    stream.write(reinterpret_cast<const char*>(bones.data()), sizeof(Bone) * NUM_BONES);
}

void Mesh::dump_text(std::ostream &out) const {
    out<<std::fixed;
    out<<"# ";
    out.write(reinterpret_cast<const char*>(&MESH_HEADER_MAGIC), sizeof MESH_HEADER_MAGIC);
    out<<"\n---\n";
    out<<"\n# num_vertices: "<<vertices.size()<<"\n";
    out<<"vertices:\n";
    for (const Vertex &vtx : vertices) {
        out<<"-\n";
        vtx.dump_text(out);
    }
    if (vertices.size() == 0) {
        out<<"  []\n";
    }
    out<<"\n# num_triangles: "<<triangles.size()<<"\n";
    out<<"triangles:\n";
    for (const Triangle &tri : triangles) {
        out<<"- ["<<tri.a<<","<<tri.b<<","<<tri.c<<"]\n";
    }
    if (triangles.size() == 0) {
        out<<"  []\n";
    }
    out<<"\n# num_bones: "<<NUM_BONES<<"\n";
    out<<"bones:\n";
    for (const Bone &bone : bones) {
        out<<"-\n";
        out<<"  parent: "<<bone.parent<<"\n";
        out<<"  mat:    ["<<bone.mat[0]<<","<<bone.mat[1]<<","<<bone.mat[2]<<","<<bone.mat[3]<<",\n";
        out<<"           "<<bone.mat[4]<<","<<bone.mat[5]<<","<<bone.mat[6]<<","<<bone.mat[7]<<",\n";
        out<<"           "<<bone.mat[8]<<","<<bone.mat[9]<<","<<bone.mat[10]<<","<<bone.mat[11]<<",\n";
        out<<"           "<<bone.mat[12]<<","<<bone.mat[13]<<","<<bone.mat[14]<<","<<bone.mat[15]<<"]\n";
    }
    if (bones.size() == 0) {
        out<<"  []\n";
    }
}

/*======= Animation =======*/

Animation::Animation() { }

Animation::Animation(const string& filename) {
    unique_ptr<PHYSFS_File, int (&)(PHYSFS_File*)> physfp(PHYSFS_openRead(filename.c_str()), PHYSFS_close);
    unique_ptr<FILE, int (&)(FILE*)> fp(fopen(filename.c_str(), "rb"), fclose);

    if (!fp.get() && !physfp.get()) {
        std::stringstream what;
        what<<"couldn't load file "<<filename;
        throw std::runtime_error(what.str());
    }

    Header header;
    if(physfp.get()) {
        PHYSFS_read(physfp.get(), reinterpret_cast<char*>(&header), 1, sizeof(header));
    } else {
        std::fread(reinterpret_cast<char*>(&header), sizeof(header), 1, fp.get());
    }

    if (header.magic != ANIMATION_HEADER_MAGIC) {
        std::stringstream what;
        what<<"invalid format, got 0x"<<std::hex<<header.magic<<", expected 0x"<<ANIMATION_HEADER_MAGIC;
        what<<" (";
        what.write(reinterpret_cast<const char*>(&ANIMATION_HEADER_MAGIC), sizeof ANIMATION_HEADER_MAGIC);
        what<<")";
        throw std::runtime_error(what.str());
    }
    duration = header.duration;
    data.resize(header.data_count);

    if(physfp.get()) {
        PHYSFS_read(physfp.get(), reinterpret_cast<char*>(bone_anims.data()), NUM_BONES, sizeof(BoneAnim));
        PHYSFS_read(physfp.get(), reinterpret_cast<char*>(data.data()), header.data_count, sizeof(float));
    } else {
        std::fread(reinterpret_cast<char*>(bone_anims.data()), sizeof(BoneAnim), NUM_BONES, fp.get());
        std::fread(reinterpret_cast<char*>(data.data()), sizeof(float), header.data_count, fp.get());
    }
}

/**
 *  An iterator that advances the wrapped iterator N times for every 1 advancement.
 */
template <class Parent, int N>
class skip_iterator { 

public:
    Parent iter;
    const Parent end;
    typedef typename Parent::difference_type difference_type;
    typedef typename Parent::value_type value_type;
    typedef typename Parent::reference reference;
    typedef typename Parent::pointer pointer;
    typedef std::forward_iterator_tag iterator_category;

    skip_iterator(){} // default
    skip_iterator(const skip_iterator& other) : // copy
        iter(other.iter), end(other.end) {
    }
    skip_iterator(Parent iter,  const Parent end) :
        iter(iter), end(end) {
    }

    skip_iterator& operator=(const skip_iterator& other) { // assignment
        iter = other.iter;
        return *this;
    }

    bool operator==(const skip_iterator& other) const { return !(*this != other); }
    bool operator!=(const skip_iterator& other) const { return iter != other.iter && iter < end; }

    skip_iterator& operator++(){std::advance(iter, N); return *this;}
    skip_iterator& operator--(){std::advance(iter, -N); return *this;}

    reference operator*() const { return *iter; }
    pointer operator->() const { return &*iter; }
};

float bezerp(const float a, const float b, const float c, const float d, const float x) {
    return x*x*x*d - (x*x*x - 3.0f*x*x + 3.0f*x - 1.0f)*a + 3.0f*(x*x*x - 2.0f*x*x + x)*b - 3.0f*(x*x*x - x*x)*c;
}

float bezerp(const std::vector<float> &collection, const uint32_t start, const uint32_t end, const float time) {
    using std::lower_bound;
    if (time == collection[start]) {
        // time matches the first frame exactly--no need to interpolate
        return collection[start];
    } else {
        // find the sub-curve (cubic) that `time` is inbetween
        std::vector<float>::difference_type diff_start = static_cast<std::vector<float>::difference_type>(start);
        std::vector<float>::difference_type diff_end = static_cast<std::vector<float>::difference_type>(end);
        typedef skip_iterator<std::vector<float>::const_iterator, 6> knot_iterator;
        knot_iterator begin(collection.begin()+diff_start, collection.begin()+diff_end);
        knot_iterator end(collection.begin()+diff_end, collection.begin()+diff_end);
        knot_iterator right_knot_iter = lower_bound(begin, end, time);
        knot_iterator left_knot_iter(right_knot_iter);
        --left_knot_iter;
        float left_knot = *left_knot_iter;
        float left_handle = *(left_knot_iter.iter + 2);
        float right_handle = *(left_knot_iter.iter + 4);
        float right_knot = *right_knot_iter;
        // bisect a few times to find a parametric `t` such that bezerp(sub_cubic_spline_times, t) == time
        float left = 0.0f;
        float right = 1.0f;
        float t = 0.5f;
        for (int n = 0; n < 15; n++) {
            float v1 = bezerp(left_knot, left_handle, right_handle, right_knot, left);
            float v2 = bezerp(left_knot, left_handle, right_handle, right_knot, t);
            right = v1 < time != v2 < time ? t : right;
            left = v1 < time != v2 < time ? left : t;
            t = (right - left) / 2.0f + left;
        }
        return bezerp(*(left_knot_iter.iter + 1), *(left_knot_iter.iter+3), *(left_knot_iter.iter+5), *(right_knot_iter.iter+1), t);
    }
}

#include <iostream>
RigOrientation Animation::at(const float time) const {
    RigOrientation r;
    const float looped_time = fmod(time, duration);
    for (unsigned int bone = 0; bone < NUM_BONES; ++bone) {
        const BoneAnim &bone_anim = bone_anims[bone];
        if (bone_anim.loc.x.start_index != bone_anim.loc.x.end_index) {
            r.bone_orientations[bone].loc.x = bezerp(data, bone_anim.loc.x.start_index, bone_anim.loc.x.end_index, looped_time);
        }
        if (bone_anim.loc.y.start_index != bone_anim.loc.y.end_index) {
            r.bone_orientations[bone].loc.y = bezerp(data, bone_anim.loc.y.start_index, bone_anim.loc.y.end_index, looped_time);
        }
        if (bone_anim.loc.z.start_index != bone_anim.loc.z.end_index) {
            r.bone_orientations[bone].loc.z = bezerp(data, bone_anim.loc.z.start_index, bone_anim.loc.z.end_index, looped_time);
        }
        if (bone_anim.rot.w.start_index != bone_anim.rot.w.end_index) {
            r.bone_orientations[bone].rot.w = bezerp(data, bone_anim.rot.w.start_index, bone_anim.rot.w.end_index, looped_time);
        }
        if (bone_anim.rot.x.start_index != bone_anim.rot.x.end_index) {
            r.bone_orientations[bone].rot.x = bezerp(data, bone_anim.rot.x.start_index, bone_anim.rot.x.end_index, looped_time);
        }
        if (bone_anim.rot.y.start_index != bone_anim.rot.y.end_index) {
            r.bone_orientations[bone].rot.y = bezerp(data, bone_anim.rot.y.start_index, bone_anim.rot.y.end_index, looped_time);
        }
        if (bone_anim.rot.z.start_index != bone_anim.rot.z.end_index) {
            r.bone_orientations[bone].rot.z = bezerp(data, bone_anim.rot.z.start_index, bone_anim.rot.z.end_index, looped_time);
        }
        if (bone_anim.scale.x.start_index != bone_anim.scale.x.end_index) {
            r.bone_orientations[bone].scale.x = bezerp(data, bone_anim.scale.x.start_index, bone_anim.scale.x.end_index, looped_time);
        }
        if (bone_anim.scale.y.start_index != bone_anim.scale.y.end_index) {
            r.bone_orientations[bone].scale.y = bezerp(data, bone_anim.scale.y.start_index, bone_anim.scale.y.end_index, looped_time);
        }
        if (bone_anim.scale.z.start_index != bone_anim.scale.z.end_index) {
            r.bone_orientations[bone].scale.z = bezerp(data, bone_anim.scale.z.start_index, bone_anim.scale.z.end_index, looped_time);
        }
        // normalize rotation quaternion
        // TODO: only do this when necessary
        float mag = sqrt(r.bone_orientations[bone].rot.w*r.bone_orientations[bone].rot.w +
            r.bone_orientations[bone].rot.x*r.bone_orientations[bone].rot.x +
            r.bone_orientations[bone].rot.y*r.bone_orientations[bone].rot.y +
            r.bone_orientations[bone].rot.z*r.bone_orientations[bone].rot.z);
        r.bone_orientations[bone].rot.w /= mag;
        r.bone_orientations[bone].rot.x /= mag;
        r.bone_orientations[bone].rot.y /= mag;
        r.bone_orientations[bone].rot.z /= mag;


    }
    return r;
}

void Animation::write(const string& filename) const {
    // no real need to write to physfs
    ofstream stream(filename, ofstream::binary);
    Header header;
    header.magic = ANIMATION_HEADER_MAGIC;
    header.data_count = data.size();
    header.duration = duration;
    stream.write(reinterpret_cast<const char*>(&header), sizeof header);
    stream.write(reinterpret_cast<const char*>(bone_anims.data()), NUM_BONES*sizeof(Animation::BoneAnim));
    stream.write(reinterpret_cast<const char*>(data.data()), header.data_count*sizeof(float));
}

void Animation::dump_text(std::ostream &out) const {
    out<<std::fixed;
    out<<"# ";
    out.write(reinterpret_cast<const char*>(&ANIMATION_HEADER_MAGIC), sizeof ANIMATION_HEADER_MAGIC);
    out<<"\n---\n\n";
    out<<"duration: "<<duration<<"\n";
    out<<"bones:\n";
    for (unsigned int b = 0; b < NUM_BONES; ++b) {
        out<<"- # bone "<<b<<"\n";
        out<<"  loc:\n";
        out<<"    x: [ ";
        for (unsigned int d = this->bone_anims[b].loc.x.start_index; d < this->bone_anims[b].loc.x.end_index; ++d) {
            out<<this->data[d]<<(d < bone_anims[b].loc.x.end_index-1 ? ", " : "");
        }
        out<<" ]\n";

        out<<"    y: [ ";
        for (unsigned int d = this->bone_anims[b].loc.y.start_index; d < this->bone_anims[b].loc.y.end_index; ++d) {
            out<<this->data[d]<<(d < bone_anims[b].loc.y.end_index-1 ? ", " : "");
        }
        out<<" ]\n";

        out<<"    z: [ ";
        for (unsigned int d = this->bone_anims[b].loc.z.start_index; d < this->bone_anims[b].loc.z.end_index; ++d) {
            out<<this->data[d]<<(d < bone_anims[b].loc.z.end_index-1 ? ", " : "");
        }
        out<<" ]\n";
        out<<"  rot:\n";
        out<<"    w: [ ";
        for (unsigned int d = this->bone_anims[b].rot.w.start_index; d < this->bone_anims[b].rot.w.end_index; ++d) {
            out<<this->data[d]<<(d < bone_anims[b].rot.w.end_index-1 ? ", " : "");
        }
        out<<" ]\n";

        out<<"    x: [ ";
        for (unsigned int d = this->bone_anims[b].rot.x.start_index; d < this->bone_anims[b].rot.x.end_index; ++d) {
            out<<this->data[d]<<(d < bone_anims[b].rot.x.end_index-1 ? ", " : "");
        }
        out<<" ]\n";

        out<<"    y: [ ";
        for (unsigned int d = this->bone_anims[b].rot.y.start_index; d < this->bone_anims[b].rot.y.end_index; ++d) {
            out<<this->data[d]<<(d < bone_anims[b].rot.y.end_index-1 ? ", " : "");
        }
        out<<" ]\n";

        out<<"    z: [ ";
        for (unsigned int d = this->bone_anims[b].rot.z.start_index; d < this->bone_anims[b].rot.z.end_index; ++d) {
            out<<this->data[d]<<(d < bone_anims[b].rot.z.end_index-1 ? ", " : "");
        }
        out<<" ]\n";
        out<<"  scale:\n";
        out<<"    x: [ ";
        for (unsigned int d = this->bone_anims[b].scale.x.start_index; d < this->bone_anims[b].scale.x.end_index; ++d) {
            out<<this->data[d]<<(d < bone_anims[b].scale.x.end_index-1 ? ", " : "");
        }
        out<<" ]\n";

        out<<"    y: [ ";
        for (unsigned int d = this->bone_anims[b].scale.y.start_index; d < this->bone_anims[b].scale.y.end_index; ++d) {
            out<<this->data[d]<<(d < bone_anims[b].scale.y.end_index-1 ? ", " : "");
        }
        out<<" ]\n";

        out<<"    z: [ ";
        for (unsigned int d = this->bone_anims[b].scale.z.start_index; d < this->bone_anims[b].scale.z.end_index; ++d) {
            out<<this->data[d]<<(d < bone_anims[b].scale.z.end_index-1 ? ", " : "");
        }
        out<<" ]\n";
    }
}

/*************** Pose *******************/

inline void mat_multiply(const array<float, 16> &lhs, array<float, 16> &rhs){
    array<float, 16> tmp;
    get<0>(tmp) = get<0>(lhs) * get<0>(rhs) + get<4>(lhs) * get<1>(rhs) + get<8>(lhs) * get<2>(rhs) + get<12>(lhs) * get<3>(rhs);
    get<1>(tmp) = get<1>(lhs) * get<0>(rhs) + get<5>(lhs) * get<1>(rhs) + get<9>(lhs) * get<2>(rhs) + get<13>(lhs) * get<3>(rhs);
    get<2>(tmp) = get<2>(lhs) * get<0>(rhs) + get<6>(lhs) * get<1>(rhs) + get<10>(lhs) * get<2>(rhs) + get<14>(lhs) * get<3>(rhs);
    get<3>(tmp) = get<3>(lhs) * get<0>(rhs) + get<7>(lhs) * get<1>(rhs) + get<11>(lhs) * get<2>(rhs) + get<15>(lhs) * get<3>(rhs);
    get<4>(tmp) = get<0>(lhs) * get<4>(rhs) + get<4>(lhs) * get<5>(rhs) + get<8>(lhs) * get<6>(rhs) + get<12>(lhs) * get<7>(rhs);
    get<5>(tmp) = get<1>(lhs) * get<4>(rhs) + get<5>(lhs) * get<5>(rhs) + get<9>(lhs) * get<6>(rhs) + get<13>(lhs) * get<7>(rhs);
    get<6>(tmp) = get<2>(lhs) * get<4>(rhs) + get<6>(lhs) * get<5>(rhs) + get<10>(lhs) * get<6>(rhs) + get<14>(lhs) * get<7>(rhs);
    get<7>(tmp) = get<3>(lhs) * get<4>(rhs) + get<7>(lhs) * get<5>(rhs) + get<11>(lhs) * get<6>(rhs) + get<15>(lhs) * get<7>(rhs);
    get<8>(tmp) = get<0>(lhs) * get<8>(rhs) + get<4>(lhs) * get<9>(rhs) + get<8>(lhs) * get<10>(rhs) + get<12>(lhs) * get<11>(rhs);
    get<9>(tmp) = get<1>(lhs) * get<8>(rhs) + get<5>(lhs) * get<9>(rhs) + get<9>(lhs) * get<10>(rhs) + get<13>(lhs) * get<11>(rhs);
    get<10>(tmp) = get<2>(lhs) * get<8>(rhs) + get<6>(lhs) * get<9>(rhs) + get<10>(lhs) * get<10>(rhs) + get<14>(lhs) * get<11>(rhs);
    get<11>(tmp) = get<3>(lhs) * get<8>(rhs) + get<7>(lhs) * get<9>(rhs) + get<11>(lhs) * get<10>(rhs) + get<15>(lhs) * get<11>(rhs);
    get<12>(tmp) = get<0>(lhs) * get<12>(rhs) + get<4>(lhs) * get<13>(rhs) + get<8>(lhs) * get<14>(rhs) + get<12>(lhs) * get<15>(rhs);
    get<13>(tmp) = get<1>(lhs) * get<12>(rhs) + get<5>(lhs) * get<13>(rhs) + get<9>(lhs) * get<14>(rhs) + get<13>(lhs) * get<15>(rhs);
    get<14>(tmp) = get<2>(lhs) * get<12>(rhs) + get<6>(lhs) * get<13>(rhs) + get<10>(lhs) * get<14>(rhs) + get<14>(lhs) * get<15>(rhs);
    get<15>(tmp) = get<3>(lhs) * get<12>(rhs) + get<7>(lhs) * get<13>(rhs) + get<11>(lhs) * get<14>(rhs) + get<15>(lhs) * get<15>(rhs);
    rhs = tmp;
}

inline bool mat_inverse(const array<float, 16> &m, array<float, 16> &out) {
    array<float, 16> inv;
    float det;

    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
             m[4]  * m[11] * m[14] +
             m[8]  * m[6]  * m[15] -
             m[8]  * m[7]  * m[14] -
             m[12] * m[6]  * m[11] +
             m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
              m[4]  * m[10] * m[13] +
              m[8]  * m[5] * m[14] -
              m[8]  * m[6] * m[13] -
              m[12] * m[5] * m[10] +
              m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
             m[1]  * m[11] * m[14] +
             m[9]  * m[2] * m[15] -
             m[9]  * m[3] * m[14] -
             m[13] * m[2] * m[11] +
             m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
             m[0]  * m[11] * m[13] +
             m[8]  * m[1] * m[15] -
             m[8]  * m[3] * m[13] -
             m[12] * m[1] * m[11] +
             m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
             m[0]  * m[7] * m[14] +
             m[4]  * m[2] * m[15] -
             m[4]  * m[3] * m[14] -
             m[12] * m[2] * m[7] +
             m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
              m[0]  * m[6] * m[13] +
              m[4]  * m[1] * m[14] -
              m[4]  * m[2] * m[13] -
              m[12] * m[1] * m[6] +
              m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
             m[1] * m[7] * m[10] +
             m[5] * m[2] * m[11] -
             m[5] * m[3] * m[10] -
             m[9] * m[2] * m[7] +
             m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
              m[0] * m[7] * m[9] +
              m[4] * m[1] * m[11] -
              m[4] * m[3] * m[9] -
              m[8] * m[1] * m[7] +
              m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    memcpy(out.data(), inv.data(), sizeof(inv));

    return true;

}


Pose::Pose(const Mesh &mesh, const RigOrientation &rig) {
    typedef array<float, 16> Mat16;
    array< array<float, 16>, NUM_BONES> tmp;
    // TODO: Fix this to require fewer matrix multiplications by building matrices in the correct order
    // turn loc rot scale into matrices for this mesh
    unsigned int bone = 0;
    for(Mesh::Vertex::BoneIndex j = 0; j < NUM_BONES; ++j) {
        tmp[j] = {{1.f, 0, 0, 0, 0, 1.f, 0, 0, 0, 0, 1.f, 0, 0, 0, 0, 1.f }};
        // TODO: precompute this per mesh!
        mat_inverse(mesh.bones[j].mat, tmp[j]);

        if (rig[j].scale.x != 1.0f || rig[j].scale.y != 1.0f || rig[j].scale.z != 1.0f) {
            Mat16 scale = {{ rig[j].scale.x, 0, 0, 0, 0, rig[j].scale.y, 0, 0, 0, 0, rig[j].scale.z, 0, 0, 0, 0, 1.0f }};
            mat_multiply(scale, tmp[j]);
        }
        if (rig[j].rot.w != 1.0f || rig[j].rot.x != 0.0f || rig[j].rot.y != 0.0f || rig[j].rot.z != 0.0f) {

            // calculate coefficients
            float x2 = rig[j].rot.x + rig[j].rot.x;
            float y2 = rig[j].rot.y + rig[j].rot.y; 
            float z2 = rig[j].rot.z + rig[j].rot.z;
            float xx = rig[j].rot.x * x2;
            float xy = rig[j].rot.x * y2;
            float xz = rig[j].rot.x * z2;
            float yy = rig[j].rot.y * y2;
            float yz = rig[j].rot.y * z2;
            float zz = rig[j].rot.z * z2;
            float wx = rig[j].rot.w * x2;
            float wy = rig[j].rot.w * y2;
            float wz = rig[j].rot.w * z2;


            Mat16 rotate = {{
                1.0f - (yy + zz), 
                xy + wz,
                xz - wy,
                0.0f,

                xy - wz,
                1.0f - (xx + zz),
                yz + wx,
                0.0f,

                xz + wy,
                yz - wx,
                1.0f - (xx + yy),
                0.0f,

                0.0f,
                0.0f,
                0.0f,
                1.0f
            }};
            mat_multiply(rotate, tmp[j]);
        }
        if (rig[j].loc.x != 1.0f || rig[j].loc.y != 1.0f || rig[j].loc.z != 1.0f) {
            Mat16 translate = {{1.f, 0, 0, 0, 0, 1.f, 0, 0, 0, 0, 1.f, 0, rig[j].loc.x, rig[j].loc.y, rig[j].loc.z, 1.f }};
            mat_multiply(translate, tmp[j]);
        }
        // store this transformation matrix in the temporary array
        mat_multiply(mesh.bones[j].mat, tmp[j]);
    }
    // build tranform chains
    for(Mesh::Vertex::BoneIndex j = 0; j < NUM_BONES; ++j) {
        matrices[j] = {{1.f, 0, 0, 0, 0, 1.f, 0, 0, 0, 0, 1.f, 0, 0, 0, 0, 1.f }};
        auto current = j;
        mat_multiply(tmp[current], matrices[j]);
        while (mesh.bones[current].parent != current) {
            current = mesh.bones[current].parent;
            mat_multiply(tmp[current], matrices[j]);
        }
    }
}
