// Symbol table
#include "kinematics_scrape_flag.h"
// SYMBOL_TABLE_ENTRY       KinPkg    0 buildMoleculeUsingOligomerEvent BUILD-MOLECULE-USING-OLIGOMER-EVENT private ; cName=_sym_buildMoleculeUsingOligomerEvent lispName=BUILD-MOLECULE-USING-OLIGOMER-EVENT
// SYMBOL_TABLE_ENTRY       KinPkg    1 ChiDihedral_O                  CHI-DIHEDRAL                   export ; cName=_sym_ChiDihedral_O lispName=CHI-DIHEDRAL
// SYMBOL_TABLE_ENTRY   KeywordPkg    2 moleculeId                     MOLECULE-ID                    export ; cName=_sym_moleculeId lispName=MOLECULE-ID
// SYMBOL_TABLE_ENTRY       KinPkg    3 BackboneDependentRotamerLibrary_O BACKBONE-DEPENDENT-ROTAMER-LIBRARY export ; cName=_sym_BackboneDependentRotamerLibrary_O lispName=BACKBONE-DEPENDENT-ROTAMER-LIBRARY
// SYMBOL_TABLE_ENTRY       KinPkg    4 BackboneDependentRotamerSet_O  BACKBONE-DEPENDENT-ROTAMER-SET export ; cName=_sym_BackboneDependentRotamerSet_O lispName=BACKBONE-DEPENDENT-ROTAMER-SET
// SYMBOL_TABLE_ENTRY       KinPkg    5 origin                         ORIGIN                         private ; cName=_sym_origin lispName=ORIGIN
// SYMBOL_TABLE_ENTRY       KinPkg    6 AtomTemplate_O                 ATOM-TEMPLATE                  export ; cName=_sym_AtomTemplate_O lispName=ATOM-TEMPLATE
// SYMBOL_TABLE_ENTRY       KinPkg    7 resizeMoleculesEvent           RESIZE-MOLECULES-EVENT         private ; cName=_sym_resizeMoleculesEvent lispName=RESIZE-MOLECULES-EVENT
// SYMBOL_TABLE_ENTRY       KinPkg    8 Atom_O                         ATOM                           export ; cName=_sym_Atom_O lispName=ATOM
// SYMBOL_TABLE_ENTRY       KinPkg    9 ChainNode_O                    CHAIN-NODE                     export ; cName=_sym_ChainNode_O lispName=CHAIN-NODE
// SYMBOL_TABLE_ENTRY       KinPkg   10 RotamerSetBase_O               ROTAMER-SET-BASE               export ; cName=_sym_RotamerSetBase_O lispName=ROTAMER-SET-BASE
// SYMBOL_TABLE_ENTRY   KeywordPkg   11 residueId                      RESIDUE-ID                     export ; cName=_sym_residueId lispName=RESIDUE-ID
// SYMBOL_TABLE_ENTRY       KinPkg   12 CheckpointAtom_O               CHECKPOINT-ATOM                export ; cName=_sym_CheckpointAtom_O lispName=CHECKPOINT-ATOM
// SYMBOL_TABLE_ENTRY       KinPkg   13 sigmas                         SIGMAS                         private ; cName=_sym_sigmas lispName=SIGMAS
// SYMBOL_TABLE_ENTRY       KinPkg   14 atom                           ATOM                           private ; cName=_sym_atom lispName=ATOM
// SYMBOL_TABLE_ENTRY       KinPkg   15 MonomerNode_O                  MONOMER-NODE                   export ; cName=_sym_MonomerNode_O lispName=MONOMER-NODE
// SYMBOL_TABLE_ENTRY       KinPkg   16 FoldTree_O                     FOLD-TREE                      export ; cName=_sym_FoldTree_O lispName=FOLD-TREE
// SYMBOL_TABLE_ENTRY       KinPkg   17 count                          COUNT                          private ; cName=_sym_count lispName=COUNT
// SYMBOL_TABLE_ENTRY       KinPkg   18 RotamerSet_O                   ROTAMER-SET                    export ; cName=_sym_RotamerSet_O lispName=ROTAMER-SET
// SYMBOL_TABLE_ENTRY       KinPkg   19 AtomTree_O                     ATOM-TREE                      export ; cName=_sym_AtomTree_O lispName=ATOM-TREE
// SYMBOL_TABLE_ENTRY       KinPkg   20 CheckpointOutPlugAtom_O        CHECKPOINT-OUT-PLUG-ATOM       export ; cName=_sym_CheckpointOutPlugAtom_O lispName=CHECKPOINT-OUT-PLUG-ATOM
// SYMBOL_TABLE_ENTRY       KinPkg   21 MonomerBaseNode_O              MONOMER-BASE-NODE              export ; cName=_sym_MonomerBaseNode_O lispName=MONOMER-BASE-NODE
// SYMBOL_TABLE_ENTRY       KinPkg   22 AggregateNode_O                AGGREGATE-NODE                 export ; cName=_sym_AggregateNode_O lispName=AGGREGATE-NODE
// SYMBOL_TABLE_ENTRY       KinPkg   23 delayed                        DELAYED                        private ; cName=_sym_delayed lispName=DELAYED
// SYMBOL_TABLE_ENTRY       KinPkg   24 jump                           JUMP                           private ; cName=_sym_jump lispName=JUMP
// SYMBOL_TABLE_ENTRY   KeywordPkg   25 oldValue                       OLD-VALUE                      export ; cName=_sym_oldValue lispName=OLD-VALUE
// SYMBOL_TABLE_ENTRY       KinPkg   26 BondedAtomTemplate_O           BONDED-ATOM-TEMPLATE           export ; cName=_sym_BondedAtomTemplate_O lispName=BONDED-ATOM-TEMPLATE
// SYMBOL_TABLE_ENTRY       KinPkg   27 probability                    PROBABILITY                    private ; cName=_sym_probability lispName=PROBABILITY
// SYMBOL_TABLE_ENTRY       KinPkg   28 rootBonded                     ROOT-BONDED                    private ; cName=_sym_rootBonded lispName=ROOT-BONDED
// SYMBOL_TABLE_ENTRY       KinPkg   29 RingClosingMonomerNode_O       RING-CLOSING-MONOMER-NODE      export ; cName=_sym_RingClosingMonomerNode_O lispName=RING-CLOSING-MONOMER-NODE
// SYMBOL_TABLE_ENTRY       KinPkg   30 Rotamer_O                      ROTAMER                        export ; cName=_sym_Rotamer_O lispName=ROTAMER
// SYMBOL_TABLE_ENTRY       KinPkg   31 BondId_O                       BOND-ID                        export ; cName=_sym_BondId_O lispName=BOND-ID
// SYMBOL_TABLE_ENTRY       KinPkg   32 Checkpoint_O                   CHECKPOINT                     export ; cName=_sym_Checkpoint_O lispName=CHECKPOINT
// SYMBOL_TABLE_ENTRY       KinPkg   33 phi                            PHI                            private ; cName=_sym_phi lispName=PHI
// SYMBOL_TABLE_ENTRY       KinPkg   34 psi                            PSI                            private ; cName=_sym_psi lispName=PSI
// SYMBOL_TABLE_ENTRY       KinPkg   35 indices                        INDICES                        private ; cName=_sym_indices lispName=INDICES
// SYMBOL_TABLE_ENTRY       KinPkg   36 ChiList_O                      CHI-LIST                       export ; cName=_sym_ChiList_O lispName=CHI-LIST
// SYMBOL_TABLE_ENTRY       KinPkg   37 dihedrals                      DIHEDRALS                      private ; cName=_sym_dihedrals lispName=DIHEDRALS
// SYMBOL_TABLE_ENTRY   KeywordPkg   38 newValue                       NEW-VALUE                      export ; cName=_sym_newValue lispName=NEW-VALUE
// SYMBOL_TABLE_ENTRY       KinPkg   39 DelayedBondedAtomTemplate_O    DELAYED-BONDED-ATOM-TEMPLATE   export ; cName=_sym_DelayedBondedAtomTemplate_O lispName=DELAYED-BONDED-ATOM-TEMPLATE
// SYMBOL_TABLE_ENTRY       KinPkg   40 Conformation_O                 CONFORMATION                   export ; cName=_sym_Conformation_O lispName=CONFORMATION
// SYMBOL_TABLE_ENTRY       KinPkg   41 RootBondedAtomTemplate_O       ROOT-BONDED-ATOM-TEMPLATE      export ; cName=_sym_RootBondedAtomTemplate_O lispName=ROOT-BONDED-ATOM-TEMPLATE
// SYMBOL_TABLE_ENTRY       KinPkg   42 MonomerId_O                    MONOMER-ID                     export ; cName=_sym_MonomerId_O lispName=MONOMER-ID
// SYMBOL_TABLE_ENTRY       KinPkg   43 bonded                         BONDED                         private ; cName=_sym_bonded lispName=BONDED
// SYMBOL_TABLE_ENTRY       KinPkg   44 rotamer                        ROTAMER                        private ; cName=_sym_rotamer lispName=ROTAMER
#ifdef KeywordPkg_SYMBOLS
DO_SYMBOL(_sym_moleculeId,2,KeywordPkg,"MOLECULE-ID",true);
DO_SYMBOL(_sym_residueId,11,KeywordPkg,"RESIDUE-ID",true);
DO_SYMBOL(_sym_oldValue,25,KeywordPkg,"OLD-VALUE",true);
DO_SYMBOL(_sym_newValue,38,KeywordPkg,"NEW-VALUE",true);
#endif
#ifdef KinPkg_SYMBOLS
DO_SYMBOL(_sym_buildMoleculeUsingOligomerEvent,0,KinPkg,"BUILD-MOLECULE-USING-OLIGOMER-EVENT",false);
DO_SYMBOL(_sym_ChiDihedral_O,1,KinPkg,"CHI-DIHEDRAL",true);
DO_SYMBOL(_sym_BackboneDependentRotamerLibrary_O,3,KinPkg,"BACKBONE-DEPENDENT-ROTAMER-LIBRARY",true);
DO_SYMBOL(_sym_BackboneDependentRotamerSet_O,4,KinPkg,"BACKBONE-DEPENDENT-ROTAMER-SET",true);
DO_SYMBOL(_sym_origin,5,KinPkg,"ORIGIN",false);
DO_SYMBOL(_sym_AtomTemplate_O,6,KinPkg,"ATOM-TEMPLATE",true);
DO_SYMBOL(_sym_resizeMoleculesEvent,7,KinPkg,"RESIZE-MOLECULES-EVENT",false);
DO_SYMBOL(_sym_Atom_O,8,KinPkg,"ATOM",true);
DO_SYMBOL(_sym_ChainNode_O,9,KinPkg,"CHAIN-NODE",true);
DO_SYMBOL(_sym_RotamerSetBase_O,10,KinPkg,"ROTAMER-SET-BASE",true);
DO_SYMBOL(_sym_CheckpointAtom_O,12,KinPkg,"CHECKPOINT-ATOM",true);
DO_SYMBOL(_sym_sigmas,13,KinPkg,"SIGMAS",false);
DO_SYMBOL(_sym_atom,14,KinPkg,"ATOM",false);
DO_SYMBOL(_sym_MonomerNode_O,15,KinPkg,"MONOMER-NODE",true);
DO_SYMBOL(_sym_FoldTree_O,16,KinPkg,"FOLD-TREE",true);
DO_SYMBOL(_sym_count,17,KinPkg,"COUNT",false);
DO_SYMBOL(_sym_RotamerSet_O,18,KinPkg,"ROTAMER-SET",true);
DO_SYMBOL(_sym_AtomTree_O,19,KinPkg,"ATOM-TREE",true);
DO_SYMBOL(_sym_CheckpointOutPlugAtom_O,20,KinPkg,"CHECKPOINT-OUT-PLUG-ATOM",true);
DO_SYMBOL(_sym_MonomerBaseNode_O,21,KinPkg,"MONOMER-BASE-NODE",true);
DO_SYMBOL(_sym_AggregateNode_O,22,KinPkg,"AGGREGATE-NODE",true);
DO_SYMBOL(_sym_delayed,23,KinPkg,"DELAYED",false);
DO_SYMBOL(_sym_jump,24,KinPkg,"JUMP",false);
DO_SYMBOL(_sym_BondedAtomTemplate_O,26,KinPkg,"BONDED-ATOM-TEMPLATE",true);
DO_SYMBOL(_sym_probability,27,KinPkg,"PROBABILITY",false);
DO_SYMBOL(_sym_rootBonded,28,KinPkg,"ROOT-BONDED",false);
DO_SYMBOL(_sym_RingClosingMonomerNode_O,29,KinPkg,"RING-CLOSING-MONOMER-NODE",true);
DO_SYMBOL(_sym_Rotamer_O,30,KinPkg,"ROTAMER",true);
DO_SYMBOL(_sym_BondId_O,31,KinPkg,"BOND-ID",true);
DO_SYMBOL(_sym_Checkpoint_O,32,KinPkg,"CHECKPOINT",true);
DO_SYMBOL(_sym_phi,33,KinPkg,"PHI",false);
DO_SYMBOL(_sym_psi,34,KinPkg,"PSI",false);
DO_SYMBOL(_sym_indices,35,KinPkg,"INDICES",false);
DO_SYMBOL(_sym_ChiList_O,36,KinPkg,"CHI-LIST",true);
DO_SYMBOL(_sym_dihedrals,37,KinPkg,"DIHEDRALS",false);
DO_SYMBOL(_sym_DelayedBondedAtomTemplate_O,39,KinPkg,"DELAYED-BONDED-ATOM-TEMPLATE",true);
DO_SYMBOL(_sym_Conformation_O,40,KinPkg,"CONFORMATION",true);
DO_SYMBOL(_sym_RootBondedAtomTemplate_O,41,KinPkg,"ROOT-BONDED-ATOM-TEMPLATE",true);
DO_SYMBOL(_sym_MonomerId_O,42,KinPkg,"MONOMER-ID",true);
DO_SYMBOL(_sym_bonded,43,KinPkg,"BONDED",false);
DO_SYMBOL(_sym_rotamer,44,KinPkg,"ROTAMER",false);
#endif
