import sys
import itertools
from contextlib import contextmanager

import bpy
from bpy_extras.io_utils import ExportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty
from bpy.types import Operator

if 'mesh' in locals() or 'animation' in locals():
    import importlib
    importlib.reload(mesh)
    importlib.reload(animation)

from . import mesh
from . import animation

Mesh = mesh.Mesh
Animation = animation.Animation

@contextmanager
def temp_mesh(obj, *args, **kwargs):
    """
    Generate a temporary mesh inside a 'with' context.
    The parameters are the same as
    http://www.blender.org/documentation/blender_python_api_2_71_0/bpy.types.Object.html#bpy.types.Object.to_mesh
    
    >>> with temp_mesh(obj=myobject, scene=bpy.context.scene, apply_modifiers=True, settings='PREVIEW') as m:
    >>>     do_whatever_you_want_with(m) # won't affect the scene
    """
    mesh = obj.to_mesh(*args, **kwargs)
    try:
        yield mesh
    finally:
        bpy.data.meshes.remove(mesh)

def armature_for_obj(obj):
    armatures = [m.object for m in obj.modifiers if m.type == "ARMATURE"]
    if not armatures:
        return None
    elif len(armatures) == 1:
        return armatures[0]
    else:
        raise RuntimeError("expected 0 or 1 armature, got {}".format(len(armatures)))
    return armatures[0]

def pairs(iterable):
    i,j = itertools.tee(iterable)
    next(j)
    while True:
        yield next(i), next(j)

def flat_fcurve(fcurve):
    for kf_a, kf_b in pairs(fcurve.keyframe_points):
        yield kf_a.co.x
        yield kf_a.co.y
        yield kf_a.handle_right.x
        yield kf_a.handle_right.y
        yield kf_b.handle_left.x
        yield kf_b.handle_left.y
    yield fcurve.keyframe_points[-1].co.x
    yield fcurve.keyframe_points[-1].co.y

def flat_matrix(matrix):
    for vec in matrix.col:
        for val in vec:
            yield val

class BoneResolver(object):
    def __init__(self, obj):
        self.obj = obj
        self.armature = armature_for_obj(obj)
        self.name_lookup = dict()
        if self.armature:
            self.bone_names = [b.name for b in self.armature.data.bones]
            self.bone_names.sort()
            for i, bone_name in enumerate(self.bone_names):
                self.name_lookup[bone_name] = i
        else:
            self.bone_names = []
    def __getitem__(self,index):
        return self.name_lookup.__getitem__(self.obj.vertex_groups[index.group].name)
    def __contains__(self, item):
        return self.name_lookup.__contains__(self.obj.vertex_groups[item.group].name)
    @property
    def bones(self):
        return (self.armature.data.bones[bone_name] for bone_name in self.bone_names)
    @property
    def matrices(self):
        return (tuple(flat_matrix(self.armature.data.bones[bone_name].matrix_local)) for bone_name in self.bone_names)
    @property
    def parents(self):
        for i, bone_name in enumerate(self.bone_names):
            parent_bone = self.armature.data.bones[bone_name].parent
            if parent_bone != None:
                yield self.name_lookup[parent_bone.name]
            else:
                yield i

def make_vert(obj, mesh, face, i):
    """
    Given an object, its mesh (or a temporary one), its tessface,
    and face vert index ( [0-2] for tris, [0-3] for quads,
    return a Mesh.Vertex object
    """
    vert = mesh.vertices[face.vertices[i]]
    co = obj.matrix_world * vert.co
    pos = Mesh.Position(co.x, co.y, co.z)
    norm = Mesh.Normal(vert.normal.x, vert.normal.y, vert.normal.z)
    if mesh.uv_textures:
        uv = mesh.tessface_uv_textures[0].data[face.index].uv[i]
        tex = Mesh.TexCoord(uv[0], uv[1])
    else:
        tex = Mesh.TexCoord( [0.0, 1.0][i % 2], [0.0, 1.0][(i % 4) // 2])

    resolver = BoneResolver(obj)
    if resolver.armature:
        bone_weights = []
        for g in vert.groups:
            if g in resolver:
                bone_weights.append((g.weight, resolver[g]))
        bone_weights.sort()
        bone_weights.reverse()
        while len(bone_weights) < 4:
            bone_weights.append((0.0,0))
        bone_weights = bone_weights[:4]
        weights = Mesh.WeightList(*(bw[0] for bw in bone_weights))
        bones = Mesh.BoneList(*(bw[1] for bw in bone_weights))
    else:
        weights = Mesh.WeightList(0.0, 0.0, 0.0, 0.0)
        bones = Mesh.BoneList(0, 0, 0, 0)
    return Mesh.Vertex(pos, norm, tex, weights, bones)

def make_mesh(scene, obj):
    """
    Given a scene, a blender MESH object
    return a Mesh, triangulating all quads
    """
    mesh = Mesh()

    with temp_mesh(obj, scene, True, 'PREVIEW') as tmesh:
        verts = []
        vert_lookup = dict() # key: vert, value: vert index
        def get_vert_id(vert):
            if vert not in vert_lookup:
                vert_lookup[vert] = len(verts)
                verts.append(vert)
            return vert_lookup[vert]
        for face in tmesh.tessfaces:
            a = get_vert_id(make_vert(obj, tmesh, face, 0))
            b = get_vert_id(make_vert(obj, tmesh, face, 1))
            c = get_vert_id(make_vert(obj, tmesh, face, 2))
            mesh.triangles.append(Mesh.Triangle(a,b,c))

            # triangulate quads
            if len(face.vertices) == 4:
                b = c
                c = get_vert_id(make_vert(obj, tmesh, face, 3))
                mesh.triangles.append(Mesh.Triangle(a,b,c))
        mesh.vertices = verts
    resolver = BoneResolver(obj)
    if resolver.armature:
        for i, (parent, mat) in enumerate(zip(resolver.parents, resolver.matrices)):
            mesh.bones.append(Mesh.Bone(parent, mat))
    return mesh

def make_anim(scene, obj):
    """
    Given a scene, a blender ARMATURE object
    return an Animation
    """
    anim = Animation()

    armature = obj.data
    action = obj.animation_data.action

    bone_names = [b.name for b in armature.bones]
    bone_names.sort()
    get_bone_index = {name: num for num, name in enumerate(bone_names)}

    # lookup tables for FCurve's array_index property
    rot_dimension = "wxyz"
    loc_dimension = "xyz"
    scale_dimension = "xyz"

    for group in action.groups:
        bone_index = get_bone_index[group.name]

        time_min, time_max = float('inf'), float('-inf')
        for fcurve in group.channels:
            start_index = len(anim.data)
            anim.data.extend(flat_fcurve(fcurve))
            end_index = len(anim.data)
            time_min = min(time_min, min(anim.data[start_index:end_index:6]))
            time_max = max(time_max, max(anim.data[start_index:end_index:6]))

            curve = Animation.Curve(start_index, end_index)

            channel_name = fcurve.data_path[fcurve.data_path.rfind('.')+1:]
            if channel_name == 'rotation_quaternion':
                setattr(anim.bone_anims[bone_index].rot, rot_dimension[fcurve.array_index], curve)
            elif channel_name == 'location':
                setattr(anim.bone_anims[bone_index].loc, loc_dimension[fcurve.array_index], curve)
            elif channel_name == 'scale':
                setattr(anim.bone_anims[bone_index].scale, scale_dimension[fcurve.array_index], curve)
            else:
                raise RuntimeError("Unknown dimension {}".format(repr(channel_name)))

    # make sure the animation starts exactly at 0
    for i in range(0, len(anim.data), 2):
        anim.data[i] -= time_min
    # scale the animation time to seconds
    for i in range(0, len(anim.data), 2):
        anim.data[i] /= scene.render.fps
    # set the duration to the largest time value
    anim.duration = (time_max-time_min)/scene.render.fps

    return anim

class MBMeshExport(Operator, ExportHelper):
    """This appears in the tooltip of the operator and in the generated docs"""
    bl_idname = "export_mesh.mbmesh"  # important since its how bpy.ops.import_test.some_data is constructed
    bl_label = "Export .mbmesh"

    # ExportHelper mixin class uses this
    filename_ext = ".mbmesh"

    filter_glob = StringProperty(
            default="*.mbmesh",
            options={'HIDDEN'},
            )

    # List of operator properties, the attributes will be assigned
    # to the class instance from the operator settings before calling.
    use_setting = BoolProperty(
            name="Example Boolean",
            description="Example Tooltip",
            default=True,
            )

    type = EnumProperty(
            name="Example Enum",
            description="Choose between two items",
            items=(('OPT_A', "First Option", "Description one"),
                   ('OPT_B', "Second Option", "Description two")),
            default='OPT_A',
            )

    @classmethod
    def poll(cls, context):
        return context.active_object.type == 'MESH'

    def execute(self, context):
        print("running MBMeshExport.execute...")
        if context.active_object.type != 'MESH':
            raise RuntimeError("No mesh object selected")
        mesh = make_mesh(context.scene, context.active_object)

        with open(self.filepath, 'wb') as f:
            mesh.write(f)
        print("MBMeshExport.execute finished")

        return {'FINISHED'}

class MBAnimationExport(Operator, ExportHelper):
    """This appears in the tooltip of the operator and in the generated docs"""
    bl_idname = "export_anim.mbanim"  # important since its how bpy.ops.import_test.some_data is constructed
    bl_label = "Export .mbanim"

    # ExportHelper mixin class uses this
    filename_ext = ".mbanim"

    filter_glob = StringProperty(
            default="*.mbanim",
            options={'HIDDEN'},
            )

    # List of operator properties, the attributes will be assigned
    # to the class instance from the operator settings before calling.
    use_setting = BoolProperty(
            name="Example Boolean",
            description="Example Tooltip",
            default=True,
            )

    type = EnumProperty(
            name="Example Enum",
            description="Choose between two items",
            items=(('OPT_A', "First Option", "Description one"),
                   ('OPT_B', "Second Option", "Description two")),
            default='OPT_A',
            )

    @classmethod
    def poll(cls, context):
        return context.active_object.type == 'ARMATURE' and context.active_object.animation_data.action

    def execute(self, context):
        print("running {self.__class__.__name__}.execute...".format(self=self))
        if context.active_object.type != 'ARMATURE':
            raise NameError("No armature selected")
        obj = context.active_object
        anim = make_anim(context.scene, context.active_object)

        with open(self.filepath, 'wb') as f:
            anim.write(f)
        print("{self.__class__.__name__}.execute finished".format(self=self))

        return {'FINISHED'}
