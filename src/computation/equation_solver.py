from fractions import Fraction
from math import gcd, lcm

from computation.chem_eqn import extract
import numpy as np

import re

def stoch_matrix(equation:str)-> np.ndarray:
    reactants, products = extract(equation)

    elements = list(next(iter(reactants.values())).keys())
    compound_entries = [
        (compound, counts, 1) for compound, counts in reactants.items()
    ] + [
        (compound, counts, -1) for compound, counts in products.items()
    ]

    matrix = np.array(
        [
            [sign * counts[element] for _, counts, sign in compound_entries]
            for element in elements
        ],
        dtype=int,
    )
    return matrix


def solver(matrix) -> list[int]:
    _, _, vh = np.linalg.svd(matrix, full_matrices=True)
    rank = np.linalg.matrix_rank(matrix)
    nullspace = vh[rank:]

    if len(nullspace) != 1:
        raise ValueError(
            "Expected one independent balancing solution; "
            f"the matrix has {len(nullspace)}."
        )

    solution = nullspace[0]
    if np.all(solution < 0):
        solution = -solution
    if np.any(solution <= 0):
        raise ValueError("The nullspace solution does not give positive coefficients.")

    solution = solution / np.min(solution)
    rational_solution = [
        Fraction(float(value)).limit_denominator(10000) for value in solution
    ]
    common_denominator = lcm(*(value.denominator for value in rational_solution))
    coefficients = [
        value.numerator * (common_denominator // value.denominator)
        for value in rational_solution
    ]
    common_factor = gcd(*coefficients)
    coefficients = [value // common_factor for value in coefficients]

    if np.any(matrix @ np.array(coefficients) != 0):
        raise ValueError("Could not convert the numerical solution to exact coefficients.")

    return coefficients

def decorate(equation:str, compounds:list[str], coefficients:list[int]):
    coefficient_by_compound: dict[str, list[int]] = {}
    for compound, coefficient in zip(compounds, coefficients):
        match = re.fullmatch(r"(\d*)(.+)", compound)
        if match is None:
            raise ValueError(f"Invalid compound: {compound}")

        original_coefficient = int(match.group(1)) if match.group(1) else 1
        coefficient_by_compound.setdefault(compound, []).append(
            coefficient * original_coefficient
        )

    occurrence_by_compound: dict[str, int] = {}
    pattern = "|".join(
        re.escape(compound)
        for compound in sorted(coefficient_by_compound, key=len, reverse=True)
    )
    pattern = rf"(?<![A-Za-z0-9])(?:{pattern})(?![A-Za-z0-9])"

    def replace_compound(match: re.Match[str]) -> str:
        compound = match.group(0)
        occurrence = occurrence_by_compound.get(compound, 0)
        occurrence_by_compound[compound] = occurrence + 1
        coefficient = coefficient_by_compound[compound][occurrence]

        original = re.fullmatch(r"\d*(.+)", compound)
        if original is None:
            raise ValueError(f"Invalid compound: {compound}")

        prefix = str(coefficient) if coefficient != 1 else ""
        return f"{prefix}{original.group(1)}"

    return re.sub(pattern, replace_compound, equation)

if __name__ == "__main__":
    "Na3PO4 + Zn(NO3)2 -> Zn3(PO4)2 + NaNO3"
    equation = input("enter your unbalanced equation: \n")

    matrix = stoch_matrix(equation)
    reactants, products = extract(equation)
    compounds = list(reactants) + list(products)
    elements = sorted(
        {
            element
            for counts in list(reactants.values()) + list(products.values())
            for element in counts
        }
    )
    coefficients = solver(matrix)

    print(decorate(equation, compounds, coefficients))
    