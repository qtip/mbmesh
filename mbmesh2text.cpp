#include <stdexcept>
using std::exception;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <physfs.h>

#include "mbmesh.h"
using mbmesh::Mesh;

int main(int argc, char *argv[]) {
    PHYSFS_init(argv[0]);
    try {
        if (argc == 3) {
            Mesh mesh(argv[1]);
            ofstream out(argv[2]);
            mesh.dump_text(out);
        } else {
            cout<<"usage: "<<argv[0]<<" INFILE OUTFILE\n";
            PHYSFS_deinit();
            return 1;
        }
    } catch(exception& e) {
        cerr<<"fatal error: "<<e.what()<<endl;
        PHYSFS_deinit();
        return 3;
    } catch(...) {
        cerr<<"fatal error"<<endl;
        PHYSFS_deinit();
        return 4;
    }
    return 0;
}
