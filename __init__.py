bl_info = {
    "name": "Export MBMesh",
    "author": "Daniel Snider",
    "location": "File > Export > MBMesh (.mbmesh)",
    "version": (0,0,0),
    "blender": (2, 71, 0),
    "description": "Export MBMesh format",
    "category": "Import-Export"}

if "bpy" in locals():
    import importlib
    if 'blenderexport' in locals():
        importlib.reload(blenderexport)

try:
    import bpy
except ImportError:
    pass
if "bpy" in locals():
    from . import blenderexport

    def menu_func_export(self, context):
        self.layout.operator(blenderexport.MBMeshExport.bl_idname, text="MBMesh (.mbmesh)")
        self.layout.operator(blenderexport.MBAnimationExport.bl_idname, text="MBAnim (.mbanim)")

    def register():
        bpy.utils.register_module(__name__)
        bpy.types.INFO_MT_file_export.append(menu_func_export)

    def unregister():
        bpy.utils.unregister_module(__name__)
        bpy.types.INFO_MT_file_export.remove(menu_func_export)

    if __name__ == "__main__":
        register()
else:
    from .mesh import Mesh
    from .animation import Animation
