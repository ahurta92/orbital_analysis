import re

test = "static"
response_file = test + "/response.out"
txt = open(response_file, "r")

data = txt.read()

reNUMBER = r"""(
    (?:[-+]?\d*\.\d+(?:[DdEe][-+]?\d+)?) |   # .num with optional sign, exponent, wholenum
    (?:[-+]?\d+\.\d*(?:[DdEe][-+]?\d+)?) |   # num. with optional sign, exponent, decimals
    (?:[-+]?\d+(?:[DdEe][-+]?\d+)?)          # num with optional sign, exponent
         )"""
#
NUMBER = r"(?x:" + reNUMBER + ")"
OPTIONS = [r"Orbital Energies: [*]"]
PSIVAR = ["ORBITAL ENERGIES"]
optDict = dict(zip(OPTIONS, PSIVAR))

mobj = re.search(r"^\s+" + var + r"\s*" + NUMBER + r"s*$", outtext, re.MULTILINE)
re.search
txt.close()
# Orbital Energies: [*] -32.77244263  -1.93039104  -0.85040982  -0.85040981  -0.85040978

# Cheat Sheet for Regex
# \s One whitespace
# * zero or more
# ^ Start of string
# $ End of string
