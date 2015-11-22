#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>
#include <openbabel/stereo/stereo.h>
//#include <openbabel/canon.h>
#include <openbabel/ring.h>
#include <emscripten/bind.h>
#include <string>
#include <vector>

using namespace emscripten;
using namespace OpenBabel;
using namespace std;

class Mol
{
public:
    Mol(OBMol obmol);
    ~Mol();

    static Mol* fromSmiles(string smiles);

    // return the canonicalindex of the atoms
    std::vector<unsigned int> canonicalindex();
    std::vector<unsigned int> Mol_NMR_FP();
    std::vector<unsigned int> Atom_NMR_FP();
    std::vector<unsigned int> BFS();
    //void getRingcode(OBAtom *root);

private:
    OBMol* obmol;
    OBMol getMol() const {return *obmol;}; // access to the RDkit mol object in private function not to share with javascript!

};


Mol* passThrough(Mol* ptr) { return ptr; }

// Binding code
EMSCRIPTEN_BINDINGS(obmol) {
    // register the vectors
    register_vector<string>("VectorString");
    register_vector<unsigned int>("VectorUint");

    class_<Mol>("Mol")
    .function("canonicalindex",&Mol::canonicalindex, allow_raw_pointers())
    .function("Mol_NMR_FP",&Mol::Mol_NMR_FP, allow_raw_pointers())
    .function("Atom_NMR_FP",&Mol::Atom_NMR_FP, allow_raw_pointers())
    .function("BFS",&Mol::BFS, allow_raw_pointers())
    //.function("getRingcode",&Mol::getRingcode, allow_raw_pointers())
    .class_function("fromSmiles", &Mol::fromSmiles, allow_raw_pointers());
    
   
}








