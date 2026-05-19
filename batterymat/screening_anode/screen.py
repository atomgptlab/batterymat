from jarvis.core.atoms import Atoms
from jarvis.analysis.thermodynamics.energetics import unary_energy
from jarvis.db.figshare import data
import pandas as pd
from alignn.ff.ff import (
    AlignnAtomwiseCalculator,
    default_path,
    wt01_path,
    ForceField,
)
from jarvis.core.specie import Specie
from jarvis.db.figshare import data
import pandas as pd
from jarvis.db.jsonutils import loadjson, dumpjson

model_path = wt01_path()
calculator = AlignnAtomwiseCalculator(path=model_path, stress_wt=0.3)


def atom_to_energy(atoms):
    num_atoms = atoms.num_atoms
    atoms = atoms.ase_converter()
    atoms.calc = calculator
    forces = atoms.get_forces()
    energy = atoms.get_potential_energy()
    stress = atoms.get_stress()
    return energy  # ,forces,stress

# JARVIS-DFT 3D dataset
dft3d = data("dft_3d")
df = pd.DataFrame(dft3d)

# Simple count
ions = ["Li", "Na", "K", "Ca", "Mg", "Al"]
for i in ions:
    dff = df[df["dimensionality"] == "intercalated ion"]
    dff[i] = dff["atoms"].apply(
        lambda x: True if i in x["elements"] else False
    )
    print(i, len(dff[dff[i] == True]))

# Ref:
# https://onlinelibrary.wiley.com/doi/full/10.1002/adma.201601357
# https://www.sciencedirect.com/science/article/pii/S0378775321011162 
# https://www.sciencedirect.com/science/article/pii/S0378775319307980?via%3Dihub
def get_voltage_profile(atoms=[], element="Li", charge=None):
    if charge is None:
        charge = abs(Specie(element).element_property("min_oxid_s"))
    chem_pot = charge * unary_energy(element)

    # print('charge',element,charge)
    elements = atoms.elements
    coords = atoms.frac_coords
    lattice_mat = atoms.lattice_mat

    new_elements = []
    new_coords = []
    options = []
    for ii, i in enumerate(elements):
        if i == element:
            options.append(ii)
    #print('options',options)
    comb = []
    for i in range(len(options)):
        comb.append(options[0 : i + 1])
    #print('comb',comb)

    mem = []
    info = {}
    info["atoms"] = atoms
    info["comb"] = []
    perf_en = atom_to_energy(atoms)
    info["en"] = perf_en

    mem.append(info)
    for ii, i in enumerate(comb):
        new_elements = []
        new_coords = []
        for jj, j in enumerate(elements):
            if jj not in i:
                new_elements.append(j)
                new_coords.append(coords[jj])
        new_atoms = Atoms(
            coords=new_coords,
            lattice_mat=lattice_mat,
            elements=new_elements,
            cartesian=False,
        )
        en = atom_to_energy(new_atoms)
        info = {}
        info["en"] = en
        info["atoms"] = new_atoms  # .to_dict()
        info["comb"] = i
        mem.append(info)
    max_en = max([i["en"] for i in mem])
    voltages = []
    for ii, i in enumerate(mem):
        if ii < len(mem) - 1:
            temp = (
                mem[ii + 1]["en"] - mem[ii]["en"] - unary_energy(element) #0.925
            )  # ,len(i['comb']))
            voltages.append(temp)
    return voltages


# voltages=get_voltage_profile(atoms=atoms_nmc)

df = pd.DataFrame(dft3d)
# Make jsin/csv files for each ion
ions = ["Li", "Na", "K", "Ca", "Mg", "Al"]
for ion in ions:
    dff = df[df["dimensionality"] == "intercalated ion"]
    dff[ion] = dff["atoms"].apply(
        lambda x: True if ion in x["elements"] else False
    )
    dffion = dff[dff[ion] == True]
    mem = []
    for i, ii in dffion.iterrows():
        try:
            atoms = Atoms.from_dict(ii["atoms"])
            voltages = get_voltage_profile(atoms=atoms,element=ion)
            print(
                ion,
                ii["jid"],
                ii["formula"],
                ii["ehull"],
                max(voltages),
                voltages,
            )
            info = {}
            info["jid"] = ii["jid"]
            info["formula"] = ii["formula"]
            info["voltages"] = voltages
            info["max_voltages"] = max(voltages)
            mem.append(info)
        except:
            pass
    name = ion + ".json"
    dumpjson(data=mem, filename=name)
