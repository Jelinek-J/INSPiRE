#pragma once

#include <string>
#include <map>
#include <list>

namespace inspire {
  namespace backend {
    // Characteristic of atom that are specific for the given ATOM/ HETATM records,
    // i.e. x, y, z, occupancy, tempFactor and charge
    struct Characteristic {
      // Atom  serial number (columns 7-11 in PDB files)
      int SERIAL_NUMBER;
      // Orthogonal coordinates for X in Angstroms (columns 31-38 in PDB files)
      double X;
      // Orthogonal coordinates for Y in Angstroms (columns 39-46 in PDB files)
      double Y;
      // Orthogonal coordinates for Z in Angstroms (columns 47-54 in PDB files)
      double Z;
      // Occupancy (columns 55-60 in PDB files)
      float OCCUPANCY;
      // Temperature factor (columns 61-66 in PDB files)
      float TEMPERATURE;
      // Charge on the atom (columns 79-80 in PDB files)
      signed char CHARGE;
    };

    // Atom ~ a group of ATOM/ HETATM records with the same name (belonging to the same aminoacid, chain etc.)
    struct Atom {
      // Element symbol (columns 77-78 in PDB files)
      std::string ELEMENT;

      // Alternate location indicator (column 17 in PDB files) => characteristic
      std::map<char, Characteristic> ALTERNATIVE_LOCATIONS;
    };

    // Aminoacid ~ a basic aminoacid -
    // - e.g. group of ATOM records with the same combination of resSeq and iCode (and belonging to the same chain, model etc.)
    struct Aminoacid {
      // Residue name (columns 18-20 in PDB files for ATOM lines)
      std::string RESIDUE_NAME;
      // Atom name (columns 13-16 in PDB files) => atom
      std::map<std::string, Atom> ATOMS;

      Aminoacid(std::string residue_name) : RESIDUE_NAME(residue_name) { }
    };

    // Modified aminoacid ~ a group of HETATM records within a chain block prior 'TER' line
    struct ModifiedAminoacid : Aminoacid {
      // NOTE: Modification-parent aminoacid determination moved to Feature section as there are multiple level to what consider as related.
      // E.g. '_chem_comp.mon_nstd_parent_comp_id'; '_chem_comp.one_letter_code'; '_chem_comp.name' or '_chem_comp.pdbx_synonyms' similarity.
      //
      // Residue name (columns 18-20 in PDB files)
      // Note: for parents Aminoacid's residue_name, there is used code of original unmodified Aminoacids, or XAA if not known
      //std::string MODIFICATION;

      ModifiedAminoacid(std::string residue_name) : Aminoacid(residue_name) { }
    };

    // Chain ~ ATOM / HETATM records with the same chainID
    struct Chain {
      // [Residue sequence number (columns 23-26 in PDB files); Code for insertion of residues (column 27 in PDB files)] => aminoacid
      std::map<std::pair<int, char>, Aminoacid> AMINOACIDS;
    };

    // Model ~ block between 'MODEL' and 'EMDMDL' lines, resp. all chains for a single-model file
    struct Model {
      // Chain identifier (column 22 in PDB files) => chain
      std::map<char, Chain> CHAINS;
    };

    // Coordinates in 3D space, i.e. [x;y;z]
    using Coordinate = std::tuple<double, double, double>;
    // Coeficients for linear equaiton [x;y;z;1], i.e. transformation matrix 3D => 1D
    using TransformationVector = std::tuple<double, double, double, double>;
    // Transformation matrix for 3D space
    using TransformationMatrix = std::tuple<TransformationVector, TransformationVector, TransformationVector>;

    // Biomolecule ~ item from Remark 350 in PDB record,
    // i.e. description and a list of transformations together with affected chains
    struct Biomolecule {
      // Description of the assembly
      std::string COMMENT;
      // Transformation id => [transformation; affected chains]
      std::map<int, std::pair<TransformationMatrix, std::string>> TRANSFORMATIONS;

      Biomolecule() { };
      Biomolecule(std::string comment) : COMMENT(comment) { };
    };

    // Protein ~ single PDB record
    struct Protein {
      // Unique identifier (not checked) of the protein
      std::string ID_CODE;

      // dictionary: model_serial_number => model
      std::map<int, Model> MODELS;

      // Crystallographic symmetries (remark 290 in PDB files): operator number => transformation
      std::map<int, TransformationMatrix> CRYSTALLOGRAPHIC_SYMMETRIES;

      // Biological assemblies (remark 350 in pdb files): id => biomolecule
      std::map<int, Biomolecule> BIOMOLECULES;

      static double transform(TransformationVector vector, double x, double y, double z) {
        return std::get<0>(vector)*x+std::get<1>(vector)*y+std::get<2>(vector)*z+std::get<3>(vector);
      }
    };
  }
}