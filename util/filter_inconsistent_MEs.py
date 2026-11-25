import sys


def extraer_subexps(exp):
    stack = []
    matching_braces = {"{": "}", "{@": "}@"}

    for n, tok in enumerate(exp):
        if tok.startswith('{'):
            stack.append((tok, n+1))
        elif tok.startswith('}'):
            if len(stack) > 0:
                open_tok, pos_start = stack.pop()
                if matching_braces.get(open_tok) != tok:
                    return "-1 " + tok + " " + str(n)
            else:
                return "-1 " + exp[n+1] + " " + str(n+1) if (n + 1) < len(exp) else "-1 end"

    return str(len(stack)) + " " + str(stack) if len(stack) > 0 else 0


def main(path):
    consistent = []
    removed = 0
    with open(path, encoding="utf-8") as index:
        for line in index:
            fields = line.split()
            tokens = fields[1].replace("¬", " ").split() if len(fields) >= 2 else []
            if tokens and extraer_subexps(tokens) == 0:
                consistent.append(line)
            else:
                removed += 1

    with open(path, "w", encoding="utf-8") as index:
        index.writelines(consistent)
    print("Removed:", removed, "inconsistent spots...")


if __name__ == "__main__":
    main(sys.argv[1])
