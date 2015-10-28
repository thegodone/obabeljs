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

class Molecule
{
public:
    Molecule(OBMol *obmol);
    
    ~Molecule();

    static Molecule* fromSmiles(string smiles);



    // return the canonicalindex of the atoms
    std::vector<unsigned int> canonicalindex();
    std::vector<unsigned int> Mol_NMR_FP();
    std::vector<unsigned int> Atom_NMR_FP();
    std::vector<unsigned int> BFS();
    void getRingcode(OBAtom *root, OBMol mol);


private:
    OBMol* obmol;
    OBMol getMol() const {return *obmol;}; // access to the RDkit mol object in private function not to share with javascript!

};


Molecule* passThrough(Molecule* ptr) { return ptr; }

// Binding code
EMSCRIPTEN_BINDINGS(obmol) {
    // register the vectors
    register_vector<string>("VectorString");
    register_vector<unsigned int>("VectorUint");

    class_<Molecule>("Molecule")
    .function("canonicalindex",&Molecule::canonicalindex, allow_raw_pointers())
    .function("Mol_NMR_FP",&Molecule::Mol_NMR_FP, allow_raw_pointers())
    .function("Atom_NMR_FP",&Molecule::Atom_NMR_FP, allow_raw_pointers())
    .function("BFS",&Molecule::BFS, allow_raw_pointers())
    .function("getRingcode",&Molecule::getRingcode, allow_raw_pointers())
    .class_function("fromSmiles", &Molecule::fromSmiles, allow_raw_pointers());
    
   
}








