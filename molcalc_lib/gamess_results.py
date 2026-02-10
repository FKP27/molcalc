from typing import List, Tuple

import numpy as np

from ppqm import chembridge
from ppqm import units
from ppqm.utils import files as misc


def _assign_orbital_occupancies(
    energies: np.ndarray, n_electrons: int, tol_hartree: float = 1e-4
) -> np.ndarray:
    """
    Assign orbital occupancies using Aufbau + Hund's rule for degenerate sets.
    Returns an array of 0/1/2 occupancies with the same length as energies.
    """
    n_orb = len(energies)
    occ = np.zeros(n_orb, dtype=int)
    if n_orb == 0 or n_electrons <= 0:
        return occ

    order = np.argsort(energies)
    remaining = int(n_electrons)
    i = 0

    while i < n_orb and remaining > 0:
        group = [order[i]]
        e0 = energies[order[i]]
        i += 1

        while i < n_orb and abs(energies[order[i]] - e0) <= tol_hartree:
            group.append(order[i])
            i += 1

        # Fill singly first (Hund's rule)
        for idx in group:
            if remaining <= 0:
                break
            occ[idx] += 1
            remaining -= 1

        # Then pair
        for idx in group:
            if remaining <= 0:
                break
            if occ[idx] < 2:
                occ[idx] += 1
                remaining -= 1

    return occ


def _get_homo_lumo_indices(
    energies: np.ndarray, occ: np.ndarray, tol_hartree: float = 1e-4
) -> Tuple[List[int], List[int]]:
    """Return lists of HOMO/LUMO indices, grouping degenerate levels."""
    homo_indices: List[int] = []
    lumo_indices: List[int] = []

    if len(energies) == 0 or len(occ) == 0:
        return homo_indices, lumo_indices

    filled = np.where(occ > 0)[0]
    if filled.size == 0:
        return homo_indices, lumo_indices

    homo_idx = int(filled.max())
    homo_energy = energies[homo_idx]
    homo_indices = [
        i for i, e in enumerate(energies) if abs(e - homo_energy) <= tol_hartree
    ]

    # LUMO: first unoccupied above HOMO
    for i in range(homo_idx + 1, len(energies)):
        if occ[i] == 0:
            lumo_energy = energies[i]
            lumo_indices = [
                j
                for j, e in enumerate(energies)
                if abs(e - lumo_energy) <= tol_hartree
            ]
            break

    return homo_indices, lumo_indices


def view_gamess_calculation(calculation):
    """

    arg:
        calculation - SQLAlchemy GamessCalculation

    return:
        data - dict of information needed for the view template

    """

    # enthalpy = Column(Float)
    # charges = Column(String)
    #
    # islinear = Column(String)
    # vibjsmol = Column(String)
    # vibfreq = Column(String)
    # vibintens = Column(String)
    #
    # thermo = Column(String)
    #
    # orbitals = Column(String)
    # orbitalstxt = Column(String)
    #
    # soltotal = Column(Float)
    # solpolar = Column(Float)
    # solnonpolar = Column(Float)
    # solsurface = Column(Float)
    # soldipole = Column(String)
    # soldipoletotal = Column(Float)

    # Convert model to dictionary
    data = calculation.__dict__

    if calculation.name is None:
        data["name"] = calculation.smiles
    else:
        data["name"] = calculation.name

    # Total electron count from stored SDF (sum atomic numbers - formal charge)
    data["electron_count"] = None
    try:
        sdfstr = data.get("sdf")
        if sdfstr:
            molobj = chembridge.sdfstr_to_molobj(sdfstr)
            if molobj is not None:
                atoms = chembridge.get_atoms(molobj, type=int)
                charge = chembridge.get_charge(molobj)
                data["electron_count"] = int(np.sum(atoms) - charge)
    except Exception:
        data["electron_count"] = None

    fmt = "{:.2f}"

    data["enthalpy"] = fmt.format(data["enthalpy"] * units.calories_to_joule)

    #               E         H         G         CV        CP        S
    #            KJ/MOL    KJ/MOL    KJ/MOL   J/MOL-K   J/MOL-K   J/MOL-K
    #  ELEC.      0.000     0.000     0.000     0.000     0.000     0.000
    #  TRANS.     3.718     6.197   -36.542    12.472    20.786   143.348
    #  ROT.       3.718     3.718   -15.045    12.472    12.472    62.932
    #  VIB.     119.279   119.279   119.164     2.252     2.252     0.385
    #  TOTAL    126.716   129.194    67.577    27.195    35.509   206.665
    #  VIB. THERMAL CORRECTION E(T)-E(0) = H(T)-H(0) =        99.870 J/MOL

    thermotable = calculation.thermo
    thermotable = misc.load_array(thermotable)

    data["h_elect"] = fmt.format(thermotable[0, 1])
    data["h_trans"] = fmt.format(thermotable[1, 1])
    data["h_rotat"] = fmt.format(thermotable[2, 1])
    data["h_vibra"] = fmt.format(thermotable[3, 1])
    data["h_total"] = fmt.format(thermotable[4, 1])

    data["cp_elect"] = fmt.format(thermotable[0, 4])
    data["cp_trans"] = fmt.format(thermotable[1, 4])
    data["cp_rotat"] = fmt.format(thermotable[2, 4])
    data["cp_vibra"] = fmt.format(thermotable[3, 4])
    data["cp_total"] = fmt.format(thermotable[4, 4])

    data["s_elect"] = fmt.format(thermotable[0, 5])
    data["s_trans"] = fmt.format(thermotable[1, 5])
    data["s_rotat"] = fmt.format(thermotable[2, 5])
    data["s_vibra"] = fmt.format(thermotable[3, 5])
    data["s_total"] = fmt.format(thermotable[4, 5])

    # Molecular orbitals format + occupancy
    orbitals = misc.load_array(data["orbitals"])
    if data["electron_count"] is None:
        data["orbital_occupancies"] = []
        data["homo_indices"] = []
        data["lumo_indices"] = []
    else:
        occ = _assign_orbital_occupancies(orbitals, data["electron_count"])
        data["orbital_occupancies"] = occ.tolist()
        homo_indices, lumo_indices = _get_homo_lumo_indices(orbitals, occ)
        data["homo_indices"] = homo_indices
        data["lumo_indices"] = lumo_indices
    data["orbitals"] = [fmt.format(x) for x in orbitals * units.hartree_to_ev]

    # Vibrational Frequencies format
    data["vibfreq"] = misc.load_array(data["vibfreq"])
    islinear = int(data["islinear"]) == int(1)
    offset = 5
    if not islinear:
        offset = 6
    data["vibfreq"] = data["vibfreq"][offset:]
    data["vibfreq"] = [fmt.format(x) for x in data["vibfreq"]]
    data["viboffset"] = offset

    # Solvation calculations
    if data["charges"] is None:
        data["has_solvation"] = False

    else:

        data["has_solvation"] = True

        dipoles = misc.load_array(data["soldipole"])

        data["dipolex"] = dipoles[0]
        data["dipoley"] = dipoles[1]
        data["dipolez"] = dipoles[2]

        data["soltotal"] = fmt.format(
            data["soltotal"] * units.calories_to_joule
        )
        data["solpolar"] = fmt.format(
            data["solpolar"] * units.calories_to_joule
        )
        data["solnonpolar"] = fmt.format(
            data["solnonpolar"] * units.calories_to_joule
        )
        data["solsurface"] = fmt.format(data["solsurface"])
        data["soldipoletotal"] = fmt.format(data["soldipoletotal"])

        charges = misc.load_array(data["charges"])
        charges = np.array(charges)
        charge = np.sum(charges)
        charge = np.round(charge, decimals=0)

        data["charge"] = int(charge)

    return data
