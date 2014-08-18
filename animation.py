import sys

from collections import namedtuple
from struct import Struct, pack, unpack

class Animation(object):
    MAGIC = 0x30306d696e61626d
    NUM_BONES = 64
    Curve = namedtuple('Curve', ['start_index', 'end_index'])
    class Rotation(object):
        def __init__(self, w, x, y, z):
            self.w = w
            self.x = x
            self.y = y
            self.z = z
    class Scale(object):
        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z
    class Location(object):
        def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z

    class BoneAnim(object):
        def __init__(self, loc, rot, scale):
            self.loc = loc
            self.rot = rot
            self.scale = scale
    
    BONE_ANIM_STRUCT = Struct('<iiiiiiiiiiiiiiiiiiii')

    HEADER_STRUCT = Struct('<QfL')

    def __init__(self):
        self.data = []
        self.duration = 0.0
        self.bone_anims = []
        for i in range(self.NUM_BONES):
            self.bone_anims.append(self.BoneAnim(
                self.Location(
                    self.Curve(0,0),
                    self.Curve(0,0),
                    self.Curve(0,0)),
                self.Rotation(
                    self.Curve(0,0),
                    self.Curve(0,0),
                    self.Curve(0,0),
                    self.Curve(0,0)),
                self.Scale(
                    self.Curve(0,0),
                    self.Curve(0,0),
                    self.Curve(0,0))))

    @classmethod
    def load(cls, f):
        animation = cls()
        magic, duration, data_count = cls.HEADER_STRUCT.unpack(f.read(cls.HEADER_STRUCT.size))
        animation.duration = duration
        if magic != cls.MAGIC:
            raise RuntimeError("unsupported mbanim version")
        for b in range(cls.NUM_BONES):
            lx0, lx1, ly0, ly1, lz0, lz1, rw0, rw1, rx0, rx1, ry0, ry1, rz0, rz1, sx0, sx1, sy0, sy1, sz0, sz1 = cls.BONE_ANIM_STRUCT.unpack(f.read(cls.BONE_ANIM_STRUCT.size))
            animation.bone_anims[b] = cls.BoneAnim(
                cls.Location(
                    cls.Curve(lx0,lx1),
                    cls.Curve(ly0,ly1),
                    cls.Curve(lz0,lz1)),
                cls.Rotation(
                    cls.Curve(rw0,rw1),
                    cls.Curve(rx0,rx1),
                    cls.Curve(ry0,ry1),
                    cls.Curve(rz0,rz1)),
                cls.Scale(
                    cls.Curve(sx0,sx1),
                    cls.Curve(sy0,sy1),
                    cls.Curve(sz0,sz1)))
        data = unpack('f'*data_count, f.read(4*data_count))
        animation.data = list(data)
        return animation

    def write(self, out):
        assert len(self.bone_anims) == self.NUM_BONES
        header = self.HEADER_STRUCT.pack(self.MAGIC, self.duration, len(self.data))
        out.write(header)
        for b in self.bone_anims:
            out.write(self.BONE_ANIM_STRUCT.pack(
                b.loc.x.start_index, b.loc.x.end_index,
                b.loc.y.start_index, b.loc.y.end_index,
                b.loc.z.start_index, b.loc.z.end_index,
                b.rot.w.start_index, b.rot.w.end_index,
                b.rot.x.start_index, b.rot.x.end_index,
                b.rot.y.start_index, b.rot.y.end_index,
                b.rot.z.start_index, b.rot.z.end_index,
                b.scale.x.start_index, b.scale.x.end_index,
                b.scale.y.start_index, b.scale.y.end_index,
                b.scale.z.start_index, b.scale.z.end_index,
            ))
        out.write(pack('f'*len(self.data), *self.data))

    def dump_text(self, out=sys.stdout):
        header = self.HEADER_STRUCT.pack(self.MAGIC, self.duration, len(self.data))
        out.write("# ")
        out.write(header[:8].decode('utf-8'))
        out.write("\n---\n\n")
        out.write("duration: {self.duration:f}\n".format(self=self))
        out.write("bones:\n")
        for b in range(self.NUM_BONES):
            out.write("- # bone {:d}\n".format(b))
            out.write("  loc:\n")
            out.write("    x: [ {} ]\n".format(', '.join('{:f}'.format(d) for d in self.data[self.bone_anims[b].loc.x.start_index:self.bone_anims[b].loc.x.end_index])))
            out.write("    y: [ {} ]\n".format(', '.join('{:f}'.format(d) for d in self.data[self.bone_anims[b].loc.y.start_index:self.bone_anims[b].loc.y.end_index])))
            out.write("    z: [ {} ]\n".format(', '.join('{:f}'.format(d) for d in self.data[self.bone_anims[b].loc.z.start_index:self.bone_anims[b].loc.z.end_index])))
            out.write("  rot:\n")
            out.write("    w: [ {} ]\n".format(', '.join('{:f}'.format(d) for d in self.data[self.bone_anims[b].rot.w.start_index:self.bone_anims[b].rot.w.end_index])))
            out.write("    x: [ {} ]\n".format(', '.join('{:f}'.format(d) for d in self.data[self.bone_anims[b].rot.x.start_index:self.bone_anims[b].rot.x.end_index])))
            out.write("    y: [ {} ]\n".format(', '.join('{:f}'.format(d) for d in self.data[self.bone_anims[b].rot.y.start_index:self.bone_anims[b].rot.y.end_index])))
            out.write("    z: [ {} ]\n".format(', '.join('{:f}'.format(d) for d in self.data[self.bone_anims[b].rot.z.start_index:self.bone_anims[b].rot.z.end_index])))
            out.write("  scale:\n")
            out.write("    x: [ {} ]\n".format(', '.join('{:f}'.format(d) for d in self.data[self.bone_anims[b].scale.x.start_index:self.bone_anims[b].scale.x.end_index])))
            out.write("    y: [ {} ]\n".format(', '.join('{:f}'.format(d) for d in self.data[self.bone_anims[b].scale.y.start_index:self.bone_anims[b].scale.y.end_index])))
            out.write("    z: [ {} ]\n".format(', '.join('{:f}'.format(d) for d in self.data[self.bone_anims[b].scale.z.start_index:self.bone_anims[b].scale.z.end_index])))
