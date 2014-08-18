# mbmesh

A very small mesh/animation format, loader, and exporter designed to fit my
exact needs to render skeletally animated meshes in OpenGL.

### Install

TODO: this section

### Blender addon

Install into your blender addon directory:

```sh
cd $HOME/.config/blender/$(blender --version | head -n1 | awk '{ print $2 }')/scripts/addons && git clone https://github.com/qtip/mbmesh.git
```

Activate the addon in the addon tab under user preferences, then export by choosing `File > Export > MBMesh (.mbmesh)`

## Use

```c++
#include <vector>
using std::vector;
#include "mbmesh.h"
using mbmesh::Mesh;
using mbmesh::Animation;
using mbmesh::Pose;

struct Assets {
    vector<Mesh> meshes;
    vector<Animation> animations;
} assets;

int main(int argv, char **argv) {
    const int HERO = assets.meshes.size();
    assets.meshes.emplace_back("./hero.mbmesh");
    const int WALK = assets.meshes.size();
    assets.animations.emplace_back("./walk.mbanim");

    // ...

    glBindBuffer(GL_ARRAY_BUFFER, array_buffer);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(Mesh::Vertex) * assets.meshes[HERO].vertices.size(),
        assets.meshes[HERO].vertices.data(),
        GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_array_buffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
        sizeof(Mesh::Triangle) * assets.meshes[HERO].triangles.size(),
        assets.meshes[HERO].triangles.data(),
        GL_STATIC_DRAW);

    // ...

    Pose pose(assets.meshes[HERO], assets.animations[WALK].at(time));
    glUniformMatrix4fv(bone_matrix_uniform_location, 64, GL_FALSE, pose);
    glDrawElements(GL_TRIANGLES, assets.meshes[HERO].triangles.size()*3, GL_UNSIGNED_INT, (void*)0);

    // ...

    return 0;
}
```
