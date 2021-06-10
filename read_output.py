import re
import numpy as np

from qcelemental.models import Molecule
from qcelemental.models.results import AtomicResultProperties
from qcelemental.molparse import regex

from qcengine.programs.util.pdict import PreservingDict

test = "static"
response_file = test + "/response.out"
txt = open(response_file, "r")

data = txt.read()


def harvest_outfile_pass(outtext):

    psivar = PreservingDict()

    NUMBER = r"(?x:" + regex.NUMBER + ")"  # NUMBER
    NUMSPACE = NUMBER + r"\s*"  # NUMBER + SPACE

    OPTIONS = [
        r"Number of Response States:",
        r"Number of Ground States:",
        r"k =",
    ]
    PSIVAR = ["NUM STATES", "NUM ORBITALS", "K"]
    optDict = dict(zip(OPTIONS, PSIVAR))
    txt.close()

    for var, VAR in optDict.items():
        mobj = re.search(r"^\s*" + var + r"\s*" + NUMBER + r"\s*$", data, re.MULTILINE)
        print(mobj)
        if mobj:
            psivar[VAR] = mobj.group(1)
    # Grab the Orbital Energies  There are NUM ORBITALS
    num_states = int(psivar["NUM STATES"])
    num_orbitals = int(psivar["NUM ORBITALS"])

    print(num_states)
    print(num_orbitals)
    # print(NUMSPACE)
    NUMSPACEORB = str()
    for i in range(num_orbitals):
        NUMSPACEORB += NUMSPACE
    # print(NUMSPACEORB)

    var = r"Orbital Energies: \[\*\]"
    VAR = "ORBITAL ENERGIES"
    mobj = re.search(
        r"^\s*" + var + r"\s*" + NUMSPACEORB + r"$",
        data,
        re.MULTILINE,
    )

    # print(mobj)

    if mobj:
        oe_list = []
        for i in range(num_orbitals):
            oe_list.append(mobj.group(i + 1))

        psivar[VAR] = np.array(oe_list, dtype=float)

    def grab_tensor(var, VAR, row, col, psivar):
        first_line = r"^\s*" + var + r"\s+"
        NUMBER = r"(?x:" + regex.NUMBER + ")"  # NUMBER
        NUMSPACE = NUMBER + r"\s*"  # NUMBER + SPACE
        print(first_line)

        CAPTURE_LINE = str()
        for j in range(col):
            CAPTURE_LINE += NUMSPACE
        total = first_line
        for i in range(row):
            front = r"^\[" + str(i) + r",\*\]\s*"
            line = front + CAPTURE_LINE
            total += line
            print(line)

        mobj = re.search(
            total,
            data,
            re.MULTILINE,
        )
        print(mobj)
        if mobj:
            oe_list = []
            for i in range(row):
                for j in range(col):
                    oe_list.append(mobj.group(i + 1))
            tensor = np.array(oe_list)
            psivar[VAR] = tensor.reshape((row, col))
        return psivar

    psivar = grab_tensor(r"Ground state overlap:", "OVERLAP", 5, 5, psivar)
    psivar = grab_tensor(r"Ground state hamiltonian:", "HAMILTONIAN", 5, 5, psivar)
    return psivar


# Orbital Energies: [*] -32.77244263  -1.93039104  -0.85040982  -0.85040981  -0.85040978

# Cheat Sheet for Regex
# \s One whitespace
# * zero or more
# ^ Start of string
# $ End of string
