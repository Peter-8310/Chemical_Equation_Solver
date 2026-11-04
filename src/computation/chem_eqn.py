from collections import defaultdict

def parse_formula(formula: str):
    counts = defaultdict(int)
    i = 0

    while i < len(formula):
        if formula[i] == '(':
            j = formula.find(')', i)
            if j == -1:
                raise ValueError(f"Unmatched parentheses in formula: {formula}")

            inner_formula = formula[i + 1:j]
            multiplier_start = j + 1
            multiplier = ""
            while multiplier_start < len(formula) and formula[multiplier_start].isdigit():
                multiplier += formula[multiplier_start]
                multiplier_start += 1

            multiplier_value = int(multiplier) if multiplier else 1
            for element, count in parse_formula(inner_formula).items():
                counts[element] += count * multiplier_value

            i = multiplier_start
            continue

        if formula[i].isupper():
            element = formula[i]
            i += 1
            if i < len(formula) and formula[i].islower():
                element += formula[i]
                i += 1

            number = ""
            while i < len(formula) and formula[i].isdigit():
                number += formula[i]
                i += 1

            counts[element] += int(number) if number else 1
            continue

        i += 1

    return dict(counts)


def get_all_elements(reaction: str):
    all_elements = set()
    for side in reaction.split("->"):
        for term in side.split("+"):
            term = term.strip()
            if not term:
                continue

            while term and term[0].isdigit():
                term = term[1:]

            all_elements.update(parse_formula(term).keys())

    return sorted(all_elements)


def extract(
    reaction: str,
) -> tuple[dict[str, dict[str, int]], dict[str, dict[str, int]]]:
    left_side, right_side = [side.strip() for side in reaction.split("->")]
    all_elements = get_all_elements(reaction)

    def count_side(side: str):
        compounds = {}

        for term in side.split("+"):
            term = term.strip()
            if not term:
                continue

            coefficient = ""
            index = 0
            while index < len(term) and term[index].isdigit():
                coefficient += term[index]
                index += 1

            multiplier = int(coefficient) if coefficient else 1
            formula = term[index:]
            counts = {element: 0 for element in all_elements}

            for element, count in parse_formula(formula).items():
                counts[element] = count * multiplier

            compounds[term] = counts

        return compounds

    return count_side(left_side), count_side(right_side)

if __name__ == "__main__":
    equation = "Na3PO4 + Zn(NO3)2 -> Zn3(PO4)2 + 6NaNO3"
    reactants, products = extract(equation)
    print("Reactants:")
    for compound, counts in reactants.items():
        print(f"  {compound}: {counts}")

    print("\nProducts:")
    for compound, counts in products.items():
        print(f"  {compound}: {counts}")
