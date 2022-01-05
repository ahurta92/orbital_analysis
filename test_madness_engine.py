import qcengine as qcng
import qcelemental as qcel

mol = qcel.models.Molecule.from_data("""
     O  0.0  0.000  -0.129
     H  0.0 -1.494  1.027
     H  0.0  1.494  1.027
 """)
program ="madness"
basis=None
keywords={"dft__k": 7, "dft__aobasis": "sto-3g", "dft__econv": 1.0000e-05}

inp = qcel.models.AtomicInput(
    molecule=mol,
    driver="energy",
    model={"method": "hf", "basis": None},
    keywords=keywords
)
qcng.get_program("madness")

ret=qcng.compute(inp,"madness")
