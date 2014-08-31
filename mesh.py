import sys
from collections import namedtuple
from struct import Struct

class Mesh(object):
    MAGIC = 0x30306873656d626d
    NUM_BONES = 64
    Vertex = namedtuple('Vertex', ['pos', 'norm', 'tex', 'weights', 'bones'])
    Position = namedtuple('Position', ['x', 'y', 'z'])
    Normal = namedtuple('Normal', ['x', 'y', 'z'])
    TexCoord = namedtuple('TexCoord', ['u', 'v'])
    Triangle = namedtuple('Triangle', ['a', 'b', 'c'])
    WeightList = namedtuple('WeightList', ['w0', 'w1', 'w2', 'w3'])
    BoneList = namedtuple('BoneList', ['b0', 'b1', 'b2', 'b3'])
    Bone = namedtuple('Bone', ['parent', 'mat'])

    HEADER_STRUCT = Struct('<QLL')
    VERTEX_STRUCT = Struct('<ffffffffffffLLLL')
    TRIANGLE_STRUCT = Struct('<LLL')
    BONE_STRUCT = Struct('<Lffffffffffffffff')

    def __init__(self):
        self.vertices = []
        self.triangles = []
        self.bones = []

    @classmethod
    def load(cls, f):
        mesh = Mesh()
        magic, num_vertices, num_triangles = cls.HEADER_STRUCT.unpack(f.read(cls.HEADER_STRUCT.size))
        if magic != cls.MAGIC:
            raise RuntimeError("unsupported mbmesh version")
        for v in range(num_vertices):
            x,y,z,nx,ny,nz,u,v,w0,w1,w2,w3,b0,b1,b2,b3 = cls.VERTEX_STRUCT.unpack(f.read(cls.VERTEX_STRUCT.size))
            vertex = Mesh.Vertex(
                    Mesh.Position(x,y,z),
                    Mesh.Normal(nx,ny,nz),
                    Mesh.TexCoord(u,v),
                    Mesh.WeightList(w0,w1,w2,w3),
                    Mesh.BoneList(b0,b1,b2,b3))
            mesh.vertices.append(vertex)
        for t in range(num_triangles):
            a,b,c = cls.TRIANGLE_STRUCT.unpack(f.read(cls.TRIANGLE_STRUCT.size))
            triangle = cls.Triangle(a,b,c)
            mesh.triangles.append(triangle)
        for b in range(cls.NUM_BONES):
            parent, *mat = cls.BONE_STRUCT.unpack(f.read(cls.BONE_STRUCT.size))
            bone = cls.Bone(parent,mat)
            mesh.bones.append(bone)
        return mesh

    def write(self, out):
        header = self.HEADER_STRUCT.pack(self.MAGIC, len(self.vertices), len(self.triangles))
        out.write(header)
        for v in self.vertices:
            v_data = []
            v_data.extend(v.pos)
            v_data.extend(v.norm)
            v_data.extend(v.tex)
            v_data.extend(v.weights)
            v_data.extend(v.bones)
            out.write(self.VERTEX_STRUCT.pack(*v_data))
        for t in self.triangles:
            out.write(self.TRIANGLE_STRUCT.pack(*t))
        for j in self.bones:
            j_data = []
            j_data.append(j.parent)
            j_data.extend(j.mat)
            out.write(self.BONE_STRUCT.pack(*j_data))
        # fill in missing bones
        for j in range(len(self.bones), self.NUM_BONES):
            out.write(self.BONE_STRUCT.pack(j, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))


    def dump_text(self, out=sys.stdout):
        header = self.HEADER_STRUCT.pack(self.MAGIC, len(self.vertices), len(self.triangles))
        out.write("# ")
        out.write(header[:8].decode('utf-8'))
        out.write("\n---\n")
        out.write("\n# num_vertices: {}\n".format(len(self.vertices)))
        out.write("vertices:\n")
        for v in self.vertices:
            out.write("-\n")
            out.write("  pos:     [{:f},{:f},{:f}]\n".format(*v.pos))
            out.write("  norm:    [{:f},{:f},{:f}]\n".format(*v.norm))
            out.write("  tex:     [{:f},{:f}]\n".format(*v.tex))
            out.write("  weights: [{:f},{:f},{:f},{:f}]\n".format(*v.weights))
            out.write("  bones:   [{},{},{},{}]\n".format(*v.bones))
        out.write("\n# num_triangles: {}\n".format(len(self.triangles)))
        out.write("triangles:\n")
        for t in self.triangles:
            out.write("- [{},{},{}]\n".format(*t))
        out.write("\n# num_bones: {}\n".format(self.NUM_BONES))
        out.write("bones:\n")
        for j in self.bones:
            out.write("-\n")
            out.write("  parent: {}\n".format(j.parent))
            out.write("  mat:    [{:f},{:f},{:f},{:f},\n".format(*j.mat[:4]))
            out.write("           {:f},{:f},{:f},{:f},\n".format(*j.mat[4:8]))
            out.write("           {:f},{:f},{:f},{:f},\n".format(*j.mat[8:12]))
            out.write("           {:f},{:f},{:f},{:f}]\n".format(*j.mat[12:16]))

if __name__ == "__main__":
    m = Mesh()
    m.dump_text(sys.stdout)
