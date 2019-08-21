
import numpy
import MDAnalysis


############ beg Classes

class LipidType():
    def __init__(self, defList, syst):
        self.ldef = defList
        self.syst = syst

    def setCountUpper(self, count):
        self.ldef[2] = count

    def setCountLower(self, count):
        self.ldef[3] = count

    def __str__(self):
        return self.ldef[0] + " - upper count " + str(self.ldef[2]) + " - lower count " + str(self.ldef[2])

    def getLeafletSelection(self, syst):
        if self.ldef[1] == 0:
            #return syst.select_atoms("resname %s and name %s" % (self.ldef[0], self.getHeadBeads()))
            return self.syst.select_atoms("resname %s and name %s" % (self.ldef[0], self.ldef[9]))
        else:
            return MDAnalysis.core.groups.AtomGroup([], syst)

    def getNoneLeafletSelection(self, syst):
        if self.ldef[1] == 1:
            return self.syst.select_atoms("resname %s and name %s" % (self.ldef[0], self.ldef[9]))
        else:
            return MDAnalysis.core.groups.AtomGroup([], syst)

    def getTailBeads(self):
        # ldef[8] is the bonded list
        # 1) Get all bonds for tails (Linker should be included so exclude first ldef[6] - 1
        cList = (self.ldef[8].split())[self.ldef[6] - 1:]
        # 2) Split the bonds into beads and flatten into 1D array
        #cList = numpy.array(map(lambda x:x.split("-"), cList)).flatten()
        cList = [item for sublist in map(lambda x:x.split("-"), cList) for item in sublist]
        # 3) make unique and return as a space seperated sring
        return " ".join(numpy.unique(cList))

    def getHeadBeads(self):
        # ldef[8] is the bonded list
        # 1) Get all bonds for head (Linker should be included so include ldef[6]
        cList = (self.ldef[8].split())[:self.ldef[6]]
        # 2) Split the bonds into beads and flatten into 1D array
        #cList = numpy.array(map(lambda x:x.split("-"), cList)).flatten()
        cList = [item for sublist in map(lambda x:x.split("-"), cList) for item in sublist]
        # 3) make unique and return as a space seperated sring
        return " ".join(numpy.unique(cList))

############ end Classes


def lipid_masterlist():

    # Define all lipids
    lipidTypeList = []
    # Structure of lipid type definition list
    #  0 Martini name (3-4 char)
    #  1 Flippes -- 0 if none of these lipid flip-flop / 1 lipid flip-flop  (Note then we need to update the count each frame and can not use for leaflet finding)
    #  2 Lipid count upper/outer leaflet (in current frame) --  -1 is not set
    #  3 Lipid count lower/inner leaflet (in current frame) --  -1 is not set
    #  4 Lipid charge
    #  5 Headgroup group ID (see lipidGroup list below)
    #  6 Headgroup bond count (this is how many of the bond list should excluded from tail order and to find headgroup vs tail bead names
    #  7 Lipid tails (on letter code)
    #  8 Bond list (used for order param calculations and to find
    #  9 Representative headgroup bead - use for lipid tilt and to define leaflets


    # WARNING lipid definitions have changed over time so a few lipids have diffrent versions - that can be selected here
    #version = 0 # This is default current version included notmally by always or in the else
    version = 1 # This is the PM inital version e.g. used in the normal (not new) complex PM simulations

    # GroupID         0     1     2     3     4        5     6     7       8      9      10     11
    lipidGroup = ["CHOL", "PC", "PE", "SM", "PS", "Glyco", "PI", "PA", "PIPs", "CER", "Lyso", "DAG"]

    # PC
    if version == 1:
        lipidTypeList.append(["POPX", 0,  -1, -1,  0, 1, 3, "CCDC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
        lipidTypeList.append(["POPC", 0,  -1, -1,  0, 1, 3, "CCDC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
        lipidTypeList.append(["XOPC", 0,  -1, -1,  0, 1, 3, "CCDC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
        lipidTypeList.append(["DOPC", 0,  -1, -1,  0, 1, 3, "CCDC   CCDC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-D3B D3B-C4B", "NC3"])
        lipidTypeList.append(["PUPC", 0,  -1, -1,  0, 1, 3, "DDDDDC CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-C6A C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    else:
        lipidTypeList.append(["POPX", 0,  -1, -1,  0, 1, 3, "CDCC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
        lipidTypeList.append(["POPC", 0,  -1, -1,  0, 1, 3, "CDCC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
        lipidTypeList.append(["XOPC", 0,  -1, -1,  0, 1, 3, "CDCC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
        lipidTypeList.append(["DOPC", 0,  -1, -1,  0, 1, 3, "CDCC   CDCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-D2B D2B-C3B C3B-C4B", "NC3"])
        lipidTypeList.append(["PUPC", 0,  -1, -1,  0, 1, 3, "DDDDD  CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["DIPC", 0,  -1, -1,  0, 1, 3, "CDDC   CDDC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-D2B D2B-D3B D3B-C4B", "NC3"])
    lipidTypeList.append(    ["DIPX", 0,  -1, -1,  0, 1, 3, "CDDC   CDDC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-D2B D2B-D3B D3B-C4B", "NC3"])
    lipidTypeList.append(    ["PIPX", 0,  -1, -1,  0, 1, 3, "CDDC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["PIPC", 0,  -1, -1,  0, 1, 3, "CDDC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["XIPC", 0,  -1, -1,  0, 1, 3, "CDDC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["PEPC", 0,  -1, -1,  0, 1, 3, "CDDCC  CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A C4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["PAPC", 0,  -1, -1,  0, 1, 3, "DDDDC  CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["DAPC", 0,  -1, -1,  0, 1, 3, "DDDDC  DDDDC ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-D1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         D1B-D2B D2B-D3B D3B-D4B D4B-C5B", "NC3"])
    lipidTypeList.append(    ["DPPX", 0,  -1, -1,  0, 1, 3, "CCCC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["DPPC", 0,  -1, -1,  0, 1, 3, "CCCC   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["PFPC", 0,  -1, -1,  0, 1, 3, "CDDD   CCCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-D4A                 C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["OIPC", 0,  -1, -1,  0, 1, 3, "CDDC   CDCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-D2B D2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["OIPC", 0,  -1, -1,  0, 1, 3, "CDDC   CDCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-D2B D2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["OUPC", 0,  -1, -1,  0, 1, 3, "DDDDD  CDCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         C1B-D2B D2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["OUPC", 0,  -1, -1,  0, 1, 3, "DDDDD  CDCC  ", "NC3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         C1B-D2B D2B-C3B C3B-C4B", "NC3"])

    # PE
    if version == 1:
        lipidTypeList.append(["POPE", 0,  -1, -1,  0, 2, 3, "CCDC   CCCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NH3"])
        lipidTypeList.append(["XOPE", 0,  -1, -1,  0, 2, 3, "CCDC   CCCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NH3"])
        lipidTypeList.append(["DOPE", 0,  -1, -1,  0, 2, 3, "CCDC   CCDC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-D3B D3B-C4B", "NH3"])
        lipidTypeList.append(["PUPE", 0,  -1, -1,  0, 2, 3, "DDDDDC CCCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-C6A C1B-C2B C2B-C3B C3B-C4B", "NH3"])
        lipidTypeList.append(["DUPE", 0,  -1, -1,  0, 2, 3, "DDDDDC DDDDDC", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-D1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-C6A D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-C6B", "NH3"])
    else:
        lipidTypeList.append(["POPE", 0,  -1, -1,  0, 2, 3, "CDCC   CCCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NH3"])
        lipidTypeList.append(["XOPE", 0,  -1, -1,  0, 2, 3, "CDCC   CCCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NH3"])
        lipidTypeList.append(["DOPE", 0,  -1, -1,  0, 2, 3, "CDCC   CDCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-D2B D2B-C3B C3B-C4B", "NH3"])
        lipidTypeList.append(["PUPE", 0,  -1, -1,  0, 2, 3, "DDDDD  CCCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         C1B-C2B C2B-C3B C3B-C4B", "NH3"])
        lipidTypeList.append(["DUPE", 0,  -1, -1,  0, 2, 3, "DDDDD  DDDDD ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-D1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         D1B-D2B D2B-D3B D3B-D4B D4B-D5B", "NH3"])
    lipidTypeList.append(    ["PIPE", 0,  -1, -1,  0, 2, 3, "CDDC   CCCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "NH3"])
    lipidTypeList.append(    ["PQPE", 0,  -1, -1,  0, 2, 3, "CDDDC  CCCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "NH3"])
    lipidTypeList.append(    ["PAPE", 0,  -1, -1,  0, 2, 3, "DDDDC  CCCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "NH3"])
    lipidTypeList.append(    ["XAPE", 0,  -1, -1,  0, 2, 3, "DDDDC  CCCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "NH3"])
    lipidTypeList.append(    ["DAPE", 0,  -1, -1,  0, 2, 3, "DDDDC  DDDDC ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-D1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         D1B-D2B D2B-D3B D3B-D4B D4B-C5B", "NH3"])
    lipidTypeList.append(    ["OIPE", 0,  -1, -1,  0, 2, 3, "CDDC   CDCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-D2B D2B-C3B C3B-C4B", "NH3"])
    lipidTypeList.append(    ["OIPE", 0,  -1, -1,  0, 2, 3, "CDDC   CDCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-D2B D2B-C3B C3B-C4B", "NH3"])
    lipidTypeList.append(    ["OAPE", 0,  -1, -1,  0, 2, 3, "CDDC   CDCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-D2B D2B-C3B C3B-C4B", "NH3"])
    lipidTypeList.append(    ["OAPE", 0,  -1, -1,  0, 2, 3, "CDDC   CDCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-D2B D2B-C3B C3B-C4B", "NH3"])
    lipidTypeList.append(    ["OUPE", 0,  -1, -1,  0, 2, 3, "DDDDD  CDCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         C1B-D2B D2B-C3B C3B-C4B", "NH3"])
    lipidTypeList.append(    ["OUPE", 0,  -1, -1,  0, 2, 3, "DDDDD  CDCC  ", "NH3-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         C1B-D2B D2B-C3B C3B-C4B", "NH3"])

    # PS
    if version == 1:
        lipidTypeList.append(["POPS", 0,  -1, -1, -1, 4, 3, "CCDC   CCCC  ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "CNO"])
        lipidTypeList.append(["PUPS", 0,  -1, -1, -1, 4, 3, "DDDDDC CCCC  ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-C6A C1B-C2B C2B-C3B C3B-C4B", "CNO"])
        lipidTypeList.append(["DUPS", 0,  -1, -1, -1, 4, 3, "DDDDDC DDDDDC", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-D1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-C6A D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-C6B", "CNO"])
    else:
        lipidTypeList.append(["POPS", 0,  -1, -1, -1, 4, 3, "CDCC   CCCC  ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "CNO"])
        lipidTypeList.append(["PUPS", 0,  -1, -1, -1, 4, 3, "DDDDD  CCCC  ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         C1B-C2B C2B-C3B C3B-C4B", "CNO"])
        lipidTypeList.append(["DUPS", 0,  -1, -1, -1, 4, 3, "DDDDD  DDDDD ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-D1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         D1B-D2B D2B-D3B D3B-D4B D4B-D5B", "CNO"])
    lipidTypeList.append(    ["PIPS", 0,  -1, -1, -1, 4, 3, "CDDC   CCCC  ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "CNO"])
    lipidTypeList.append(    ["PQPS", 0,  -1, -1, -1, 4, 3, "CDDDC  CCCC  ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "CNO"])
    lipidTypeList.append(    ["PAPS", 0,  -1, -1, -1, 4, 3, "DDDDC  CCCC  ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "CNO"])
    lipidTypeList.append(    ["XAPS", 0,  -1, -1, -1, 4, 3, "DDDDC  CCCC  ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "CNO"])
    lipidTypeList.append(    ["DAPS", 0,  -1, -1, -1, 4, 3, "DDDDC  DDDDC ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-D1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         D1B-D2B D2B-D3B D3B-D4B D4B-C5B", "CNO"])
    lipidTypeList.append(    ["OUPS", 0,  -1, -1, -1, 4, 3, "DDDDD  CDCC  ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         C1B-D2B D2B-C3B C3B-C4B", "CNO"])
    lipidTypeList.append(    ["DPPS", 0,  -1, -1, -1, 4, 3, "CCCC   CCCC  ", "CNO-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "CNO"])

    # PA
    if version == 1:
        lipidTypeList.append(["POPA", 0,  -1, -1, -2, 7, 2, "CCDC   CCCC  ", "        PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "PO4"])
        lipidTypeList.append(["PUPA", 0,  -1, -1, -2, 7, 2, "DDDDDC CCCC  ", "        PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-C6A C1B-C2B C2B-C3B C3B-C4B", "PO4"])
    else:
        lipidTypeList.append(["POPA", 0,  -1, -1, -2, 7, 2, "CDCC   CCCC  ", "        PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "PO4"])
        lipidTypeList.append(["PUPA", 0,  -1, -1, -2, 7, 2, "DDDDD  CCCC  ", "        PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         C1B-C2B C2B-C3B C3B-C4B", "PO4"])
    lipidTypeList.append(    ["PIPA", 0,  -1, -1, -2, 7, 2, "CDDC   CCCC  ", "        PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "PO4"])
    lipidTypeList.append(    ["PAPA", 0,  -1, -1, -2, 7, 2, "DDDDC  CCCC  ", "        PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "PO4"])

    # DG Calculate inner/outer together
    if version == 1:
        lipidTypeList.append(["PODG", 1,  -1, -1,  0,11, 1, "CCDC   CCCC  ", "                GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "GL1"])
        lipidTypeList.append(["PUDG", 1,  -1, -1,  0,11, 1, "DDDDDC CCCC  ", "                GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-C6A C1B-C2B C2B-C3B C3B-C4B", "GL1"])
    else:
        lipidTypeList.append(["PODG", 1,  -1, -1,  0,11, 1, "CDCC   CCCC  ", "                GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "GL1"])
        lipidTypeList.append(["PUDG", 1,  -1, -1,  0,11, 1, "DDDDD  CCCC  ", "                GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A         C1B-C2B C2B-C3B C3B-C4B", "GL1"])
    lipidTypeList.append(    ["PIDG", 1,  -1, -1,  0,11, 1, "CDDC   CCCC  ", "                GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "GL1"])
    lipidTypeList.append(    ["PADG", 1,  -1, -1,  0,11, 1, "DDDDC  CCCC  ", "                GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "GL1"])

    # LPC
    if version == 1:
        lipidTypeList.append(["PPC",  0,  -1, -1,  0,10, 2, "CCCC         ", "NC3-PO4 PO4-GL1         GL1-C1A          C1A-C2A C2A-C3A C3A-C4A                                        ", "NC3"])
        lipidTypeList.append(["OPC",  0,  -1, -1,  0,10, 2, "CCDC         ", "NC3-PO4 PO4-GL1         GL1-C1A          C1A-C2A C2A-D3A D3A-C4A                                        ", "NC3"])
        lipidTypeList.append(["IPC",  0,  -1, -1,  0,10, 2, "CDDC         ", "NC3-PO4 PO4-GL1         GL1-C1A          C1A-D2A D2A-D3A D3A-C4A                                        ", "NC3"])
        lipidTypeList.append(["APC",  0,  -1, -1,  0,10, 2, "DDDDC        ", "NC3-PO4 PO4-GL1         GL1-D1A          D1A-D2A D2A-D3A D3A-D4A D4A-C5A                                ", "NC3"])
        lipidTypeList.append(["UPC",  0,  -1, -1,  0,10, 2, "DDDDDC       ", "NC3-PO4 PO4-GL1         GL1-D1A          D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-C6A                        ", "NC3"])
    else:
        lipidTypeList.append(["PPC",  0,  -1, -1,  0,10, 3, "CCCC         ", "NC3-PO4 PO4-GL1 GL1-GL2 GL2-C1B          C1B-C2B C2B-C3B C3B-C4B                                        ", "NC3"])
        lipidTypeList.append(["OPC",  0,  -1, -1,  0,10, 3, "CDCC         ", "NC3-PO4 PO4-GL1 GL1-GL2 GL2-C1B          C1B-D2B D2B-C3B C3B-C4B                                        ", "NC3"])
        lipidTypeList.append(["IPC",  0,  -1, -1,  0,10, 3, "CDDC         ", "NC3-PO4 PO4-GL1 GL1-GL2 GL2-C1B          C1B-D2B D2B-D3B D3B-C4B                                        ", "NC3"])
        lipidTypeList.append(["APC",  0,  -1, -1,  0,10, 3, "DDDDC        ", "NC3-PO4 PO4-GL1 GL1-GL2 GL2-D1B          D1B-D2B D2B-D3B D3B-D4B D4B-C5B                                ", "NC3"])
        lipidTypeList.append(["UPC",  0,  -1, -1,  0,10, 3, "DDDDD        ", "NC3-PO4 PO4-GL1 GL1-GL2 GL2-D1B          D1B-D2B D2B-D3B D3B-D4B D4B-D5B                                ", "NC3"])
    lipidTypeList.append(    ["PPE",  0,  -1, -1,  0,10, 3, "CCCC         ", "NH3-PO4 PO4-GL1 GL1-GL2 GL2-C1B          C1B-C2B C2B-C3B C3B-C4B                                        ", "NH3"])
    lipidTypeList.append(    ["IPE",  0,  -1, -1,  0,10, 3, "CDDC         ", "NH3-PO4 PO4-GL1 GL1-GL2 GL2-C1B          C1B-D2B D2B-D3B D3B-C4B                                        ", "NH3"])

    # SM
    if version == 1:
        lipidTypeList.append(["POSM", 0,  -1, -1,  0, 3, 3, "TCC    CCDC  ", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-D3B D3B-C4B", "NC3"])
    else:
        lipidTypeList.append(["POSM", 0,  -1, -1,  0, 3, 3, "TCC    CDCC  ", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-D2B D2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["DPSM", 0,  -1, -1,  0, 3, 3, "TCC    CCCC  ", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["XPSM", 0,  -1, -1,  0, 3, 3, "TCC    CCCC  ", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-C4B", "NC3"])
    lipidTypeList.append(    ["DBSM", 0,  -1, -1,  0, 3, 3, "TCCC   CCCCC ", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B C4B-C5B", "NC3"])
    lipidTypeList.append(    ["DXSM", 0,  -1, -1,  0, 3, 3, "TCCCC  CCCCCC", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A C3A-C4A C4A-C5A         C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B", "NC3"])
    lipidTypeList.append(    ["PBSM", 0,  -1, -1,  0, 3, 3, "TCC    CCCCC ", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-C4B C4B-C5B", "NC3"])
    lipidTypeList.append(    ["PGSM", 0,  -1, -1,  0, 3, 3, "TCC    CCDCC ", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-D3B D3B-C4B C4B-C5B", "NC3"])
    lipidTypeList.append(    ["PNSM", 0,  -1, -1,  0, 3, 3, "TCC    CCCDCC", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B", "NC3"])
    lipidTypeList.append(    ["BNSM", 0,  -1, -1,  0, 3, 3, "TCCC   CCCDCC", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B", "NC3"])
    lipidTypeList.append(    ["XNSM", 0,  -1, -1,  0, 3, 3, "TCCCC  CCCDCC", "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A C3A-C4A C4A-C5A         C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B", "NC3"])

    # CE - Calculate inner/outer together
    lipidTypeList.append(    ["DPCE", 1,  -1, -1,  0, 9, 1, "TCC    CCCC  ", "                AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-C4B", "AM1"])
    lipidTypeList.append(    ["DBCE", 1,  -1, -1,  0, 9, 1, "TCCC   CCCCC ", "                AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-C4B C4B-C5B", "AM1"])
    lipidTypeList.append(    ["DXCE", 1,  -1, -1,  0, 9, 1, "TCCCC  CCCCCC", "                AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A C3A-C4A C4A-C5A         C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B", "AM1"])
    lipidTypeList.append(    ["POCE", 1,  -1, -1,  0, 9, 1, "TCC    CDCC  ", "                AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-D2B D2B-C3B C3B-C4B", "AM1"])
    lipidTypeList.append(    ["PNCE", 1,  -1, -1,  0, 9, 1, "TCC    CCCDCC", "                AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B", "AM1"])
    lipidTypeList.append(    ["XNCE", 1,  -1, -1,  0, 9, 1, "TCCCC  CCCDCC", "                AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A C3A-C4A C4A-C5A         C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B", "AM1"])

    # PI
    if version == 1:
        lipidTypeList.append(["POPI", 0,  -1, -1, -1, 6, 6, "CCDC   CCCC  ", "C1-C2 C2-C3 C3-C1 C1-CP   CP-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "C1"])
        lipidTypeList.append(["PIPI", 0,  -1, -1, -1, 6, 6, "CDDC   CCCC  ", "C1-C2 C2-C3 C3-C1 C1-CP   CP-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "C1"])
        lipidTypeList.append(["PAPI", 0,  -1, -1, -1, 6, 6, "DDDDC  CCCC  ", "C1-C2 C2-C3 C3-C1 C1-CP   CP-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "C1"])
        lipidTypeList.append(["PUPI", 0,  -1, -1, -1, 6, 6, "DDDDDC CCCC  ", "C1-C2 C2-C3 C3-C1 C1-CP   CP-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-C6A C1B-C2B C2B-C3B C3B-C4B", "C1"])
    else:
        lipidTypeList.append(["POPI", 0,  -1, -1, -1, 6, 6, "CDCC   CCCC  ", "C1-C2 C2-C3 C3-C1 C1-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A           C1B-C2B C2B-C3B C3B-C4B", "C1"])
        lipidTypeList.append(["PIPI", 0,  -1, -1, -1, 6, 6, "CDDC   CCCC  ", "C1-C2 C2-C3 C3-C1 C1-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-D3A D3A-C4A           C1B-C2B C2B-C3B C3B-C4B", "C1"])
        lipidTypeList.append(["PAPI", 0,  -1, -1, -1, 6, 6, "DDDDC  CCCC  ", "C1-C2 C2-C3 C3-C1 C1-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A   C1B-C2B C2B-C3B C3B-C4B", "C1"])
        lipidTypeList.append(["PUPI", 0,  -1, -1, -1, 6, 6, "DDDDD  CCCC  ", "C1-C2 C2-C3 C3-C1 C1-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-D5A   C1B-C2B C2B-C3B C3B-C4B", "C1"])

    # PIPs
    if version == 1:
        lipidTypeList.append(["POP1", 0,  -1, -1, -3, 8, 7, "CCDC   CCCC  ", "            P1-C1 C1-C2 C2-C3 C3-C1 C1-CP CP-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "C1"])
        lipidTypeList.append(["POP2", 0,  -1, -1, -5, 8, 8, "CCDC   CCCC  ", "      P2-C1 P1-C1 C1-C2 C2-C3 C3-C1 C1-CP CP-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "C1"])
        lipidTypeList.append(["POP3", 0,  -1, -1, -7, 8, 9, "CCDC   CCCC  ", "P3-C1 P2-C1 P1-C1 C1-C2 C2-C3 C3-C1 C1-CP CP-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-C2A C2A-D3A D3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "C1"])
    else:
        lipidTypeList.append(["POP1", 0,  -1, -1, -3, 8, 7, "CDCC   CCCC  ", "            P1-C1 C1-C2 C2-C3 C3-C1 C1-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "C1"])
        lipidTypeList.append(["POP2", 0,  -1, -1, -5, 8, 8, "CDCC   CCCC  ", "      P2-C1 P1-C1 C1-C2 C2-C3 C3-C1 C1-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "C1"])
        lipidTypeList.append(["POP3", 0,  -1, -1, -7, 8, 9, "CDCC   CCCC  ", "P3-C1 P2-C1 P1-C1 C1-C2 C2-C3 C3-C1 C1-PO4 PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B  C1A-D2A D2A-C3A C3A-C4A                 C1B-C2B C2B-C3B C3B-C4B", "C1"])
    lipidTypeList.append(    ["PAP1", 0,  -1, -1, -3, 8, 7, "DDDDC  CCCC  ", "            P1-C1 C1-C2 C2-C3 C3-C1 C1-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "C1"])
    lipidTypeList.append(    ["PAP2", 0,  -1, -1, -5, 8, 8, "DDDDC  CCCC  ", "      P2-C1 P1-C1 C1-C2 C2-C3 C3-C1 C1-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "C1"])
    lipidTypeList.append(    ["PAP3", 0,  -1, -1, -7, 8, 9, "DDDDC  CCCC  ", "P3-C1 P2-C1 P1-C1 C1-C2 C2-C3 C3-C1 C1-PO4 PO4-GL1 GL1-GL2 GL1-D1A GL2-C1B  D1A-D2A D2A-D3A D3A-D4A D4A-C5A         C1B-C2B C2B-C3B C3B-C4B", "C1"])

    # Glygo - NOTE the headgroup bead connectivity is not correct, insted the beads are just listed
    lipidTypeList.append(    ["DPG1", 0,  -1, -1, -1, 5, 18, "TCC    CCCC  ", "GM1 GM2 GM3 GM4 GM5 GM6 GM7 GM8 GM9 GM10 GM11 GM12 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                  C1B-C2B C2B-C3B C3B-C4B", "GM1"])
    lipidTypeList.append(    ["DBG1", 0,  -1, -1, -1, 5, 18, "TCCC   CCCCC ", "GM1 GM2 GM3 GM4 GM5 GM6 GM7 GM8 GM9 GM10 GM11 GM12 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                  C1B-C2B C2B-C3B C3B-C4B C4B-C5B", "GM1"])
    lipidTypeList.append(    ["DXG1", 0,  -1, -1, -1, 5, 18, "TCCCC  CCCCCC", "GM1 GM2 GM3 GM4 GM5 GM6 GM7 GM8 GM9 GM10 GM11 GM12 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A C3A-C4A C4A-C5A  C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B", "GM1"])
    lipidTypeList.append(    ["POG1", 0,  -1, -1, -1, 5, 18, "TCC    CDCC  ", "GM1 GM2 GM3 GM4 GM5 GM6 GM7 GM8 GM9 GM10 GM11 GM12 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                  C1B-D2B D2B-C3B C3B-C4B", "GM1"])
    lipidTypeList.append(    ["PNG1", 0,  -1, -1, -1, 5, 18, "TCC    CCCDCC", "GM1 GM2 GM3 GM4 GM5 GM6 GM7 GM8 GM9 GM10 GM11 GM12 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                  C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B", "GM1"])
    lipidTypeList.append(    ["XNG1", 0,  -1, -1, -1, 5, 18, "TCCCC  CCCDCC", "GM1 GM2 GM3 GM4 GM5 GM6 GM7 GM8 GM9 GM10 GM11 GM12 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A C3A-C4A C4A-C5A  C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B", "GM1"])

    lipidTypeList.append(    ["DPG3", 0,  -1, -1, -1, 5, 12, "TCC    CCCC  ", "GM1 GM2 GM3 GM4 GM5 GM6 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-C4B", "GM1"])
    lipidTypeList.append(    ["DBG3", 0,  -1, -1, -1, 5, 12, "TCCC   CCCCC ", "GM1 GM2 GM3 GM4 GM5 GM6 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-C4B C4B-C5B", "GM1"])
    lipidTypeList.append(    ["DXG3", 0,  -1, -1, -1, 5, 12, "TCCCC  CCCCCC", "GM1 GM2 GM3 GM4 GM5 GM6 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A C3A-C4A C4A-C5A         C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B", "GM1"])
    lipidTypeList.append(    ["POG3", 0,  -1, -1, -1, 5, 12, "TCC    CDCC  ", "GM1 GM2 GM3 GM4 GM5 GM6 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-D2B D2B-C3B C3B-C4B", "GM1"])
    lipidTypeList.append(    ["PNG3", 0,  -1, -1, -1, 5, 12, "TCC    CCCDCC", "GM1 GM2 GM3 GM4 GM5 GM6 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B", "GM1"])
    lipidTypeList.append(    ["XNG3", 0,  -1, -1, -1, 5, 12, "TCCCC  CCCDCC", "GM1 GM2 GM3 GM4 GM5 GM6 GM13 GM14 GM15 GM16 GM17  AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A C3A-C4A C4A-C5A         C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B", "GM1"])

    lipidTypeList.append(    ["DPGS", 0,  -1, -1,  0, 5, 5, "TCC    CCCC  ", "C1-C2 C1-C3 C2-C3 C3-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-C4B", "C1"])
    lipidTypeList.append(    ["PNGS", 0,  -1, -1,  0, 5, 5, "TCC    CCCDCC", "C1-C2 C1-C3 C2-C3 C3-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B", "C1"])
    lipidTypeList.append(    ["POGS", 0,  -1, -1,  0, 5, 5, "TCC    CDCC  ", "C1-C2 C1-C3 C2-C3 C3-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-D2B D2B-C3B C3B-C4B", "C1"])
    lipidTypeList.append(    ["DBGS", 0,  -1, -1,  0, 5, 5, "TCCC   CCCCC ", "C1-C2 C1-C3 C2-C3 C3-AM1 AM1-AM2 AM1-T1A AM2-C1B  T1A-C2A C2A-C3A                         C1B-C2B C2B-C3B C3B-C4B C4B-C5B", "C1"])

    # Sterols - NOTE the headgroup bead connectivity is not correct, insted the beads are just listed
    lipidTypeList.append(    ["CHOL", 1,  -1, -1,  0, 0, 1, ""             , "ROH R1 R2 R3 R4 R5 C1 C2", "ROH"])
    lipidTypeList.append(    ["XHOL", 1,  -1, -1,  0, 0, 1, ""             , "ROH R1 R2 R3 R4 R5 C1 C2", "ROH"])


    return lipidTypeList
