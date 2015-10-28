#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/canon.h>
#include <openbabel/ring.h>
#include "obabeljs.h"

using namespace OpenBabel;
using namespace std;
    
Molecule::Molecule(OBMol* mol): obmol(mol) {}
/**
 * @brief [delete]
 * @details [action on rdmol during delete process]
 * @return [description]
 */
Molecule::~Molecule() {
  if(obmol != 0) {
    //printf("Destroy called\n");
    delete obmol;
    obmol =0;
  }
}

Molecule* Molecule::fromSmiles(string smiles) {
   
          OBMol obmol;
          OBConversion conv;
          conv.SetInFormat("smi");
          conv.SetOutFormat("smi"); // for canonical use can output format!
          bool ok = conv.ReadString(&obmol, smiles); // not sure it will work!
          return new Molecule(&obmol);
}

// return the canonicalindex of the atoms
std::vector<unsigned int> Molecule::canonicalindex()
{

  OBGraphSym gs = OBGraphSym(obmol);
  std::vector<unsigned int> symclasses;
  gs.GetSymmetry(symclasses);

  std::vector<unsigned int> canonlabels;
  CanonicalLabels(obmol, symclasses, canonlabels); // need to add canon.h header to see this function properly.

  return canonlabels;
}


std::vector<unsigned int> Molecule::Mol_NMR_FP()
{
          OBSmartsPattern sp;

          std::string array[] = {"[CX4H0]","[CX4H1]","[CX4H2]","[CX4H3]","[$([CX3H0]=*)]","[$([CX3H2]=*),$([CX3H1]=*)]","[$([CX2]#*),$(*=[CX2]=*)]","[$([cX3H0](a)(a)A),$([cX3H0](a)(a)a)]","[$([cX3H1](a)a)]","[NX3;H0]","[NX3;H1]", "[NX3;H2]","[$([NX2H1]=*),$([NX2]=*)]","[$([NX1]#*)]","[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]","[$([nX3H1+0]),$([nX3H0+0]),$([nX2,nX3+])]","[NX1]~[NX2]~[NX2,NX1]","[OX2H0]","[OX2H1]","[$([OX1]=*)]","[PX3;$([H2]),$([H1]),$([H0])]", "[PX4;$([H2](=[OX1])),$([H1](=[OX1])),$([H0](=[OX1]))]","[$([SX2H1]),$([SX2H0])]","[$([SX1]=*)]","[$([SX3](=[OX1])),$([SX3+]([OX1-]))]","[$([SX4](=[OX1])(=[OX1])),$([SX4+2]([OX1-])([OX1-]))]","[FX1]","[ClX1]","[BrX1]","[IX1]"};

          std::vector<std::string> mypat(std::begin(array), std::end(array));

          unsigned int t = 30;

          std::vector<unsigned int> fpNMR;
          for (int i =0; i<mypat.size();i++){
            sp.Init(mypat[i]);
            // complete Matching 
            bool ok = sp.Match(*obmol, false);
            // Substructure mapping
            std::vector<std::vector<int> > mapListU;
            fpNMR[i]=sp.GetUMapList().size(); //  return unique matches
          }

          return fpNMR;
}



std::vector<unsigned int> Molecule::Atom_NMR_FP()
{
          OBSmartsPattern sp;
          std::string array[] = {"[CX4H0]","[CX4H1]","[CX4H2]","[CX4H3]","[$([CX3H0]=*)]","[$([CX3H2]=*),$([CX3H1]=*)]",
          "[$([CX2]#*),$(*=[CX2]=*)]","[$([cX3H0](a)(a)A),$([cX3H0](a)(a)a)]","[$([cX3H1](a)a)]","[NX3;H0]","[NX3;H1]", 
          "[NX3;H2]","[$([NX2H1]=*),$([NX2]=*)]","[$([NX1]#*)]","[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",
          "[$([nX3H1+0]),$([nX3H0+0]),$([nX2,nX3+])]","[NX1]~[NX2]~[NX2,NX1]","[OX2H0]","[OX2H1]","[$([OX1]=*)]",
          "[PX3;$([H2]),$([H1]),$([H0])]", "[PX4;$([H2](=[OX1])),$([H1](=[OX1])),$([H0](=[OX1]))]","[$([SX2H1]),$([SX2H0])]",
          "[$([SX1]=*)]","[$([SX3](=[OX1])),$([SX3+]([OX1-]))]","[$([SX4](=[OX1])(=[OX1])),$([SX4+2]([OX1-])([OX1-]))]",
          "[FX1]","[ClX1]","[BrX1]","[IX1]"};
          
          std::vector<std::string> mypat(std::begin(array), std::end(array));

          OBAtom* atom;

          std::vector<unsigned int> fpatomNMR;

          for (int i =0; i<mypat.size();i++){
            sp.Init(mypat[i]);
            // complete Matching 
            bool ok = sp.Match(*obmol, false);
            // Substructure mapping
            std::vector<std::vector<int> > mapListU;
            mapListU = sp.GetUMapList(); // Unique matches
            for (int m(0); m < mapListU.size(); ++m){
              for (int a(0); a < mapListU[m].size(); ++a){
                atom = obmol->GetAtom(mapListU[m][a]);
                fpatomNMR[i]=atom->GetIdx();
              }
            }
          }
          return fpatomNMR;
}



std::vector<unsigned int> Molecule::BFS()
{
        OBElementTable etab;

        std::vector<unsigned int> bfs;



        FOR_BFS_OF_MOL(a, obmol)
        {
          // hydrogenandring info of an atom:
          std::cout << a->GetHyb();
          std::cout << a->ImplicitHydrogenCount();
          std::cout << a->HasAlphaBetaUnsat();
          std::cout << a->MemberOfRingCount();
          std::cout << a->MemberOfRingSize();
          std::cout << a->CountRingBonds();
          std::cout << a->IsAromatic();                              
          printf("\n");

          // atom definition vector
          int atomicnum = a->GetAtomicNum();
          std::cout << atomicnum;
          std::cout << a->GetIdx();
          std::cout << a->GetHvyValence();
          std::cout << a->GetHeteroValence();
          std::cout << a->GetValence();
          std::cout << etab.GetVdwRad(atomicnum) ;
          std::cout << etab.GetElectroNeg(atomicnum);
          std::cout << etab.GetIonization(atomicnum);
          std::cout << etab.GetElectronAffinity(atomicnum);
          printf("\n");


          FOR_BONDS_OF_ATOM(b, *a)
          {
               // The variable b behaves like OBBond* when used with -> and * but
               // but needs to be explicitly converted when appearing as a parameter
               // in a function call - use &*b
           std::cout << a->GetIdx() << "   - bondtype:" << b->GetBO();
         }
         printf("\n");
        
         std::cout << "Depth:" << a.CurrentDepth() << "\n";

         std::cout << "Atomidx:"  << a->GetIdx() << "AtomType:" << a->GetType() << "\n";
       }

       printf("\n");

       return bfs;

}


void Molecule::getRingcode(OBAtom *root, OBMol mol)
        {
          int count;
          char buffer[10000];
          std::ostringstream ofs;
          OBAtom *atom;
          vector<OBRing*> vr;

          vr = mol.GetSSSR(); 
          vector<OBNodeBase*>::iterator i;

          vector<OBRing*>::iterator k;
          vector<int>::iterator j;

          int c=0;
          for (k = vr.begin();k != vr.end();++k) {
            c+=1;
              if ((*k)->_pathset[root->GetIdx()]){
                  sprintf(buffer,"%3d,%3d",c,(*k)->Size());
                  ofs << buffer;
              }
          }

        cout << ofs.str().c_str();
        return;

}










