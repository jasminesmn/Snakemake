try:
    import sys
    print("Script naam: " + sys.argv[0] + "\nParameter: " + sys.argv[1] +
          "\nDit script voert de analyse uit voor de mutaties van het DNA van de "
          "geinfecteerde sluipwesp tegenover het DNA van de niet geinfecteerde sluipwesp. "
          "\n\nDit script kan alleen aangeroepen worden met snakemake met de volgende commandline: "
          "snakemake --snakefile Snakemake_s1113710 rule analyse\n"
          "\nWanneer het script wordt aangeroepen met python3 volgt de volgende foutmelding: \n")
except IndexError:
    pass


def variant(line, deleties, inserties, mutaties):
    """
    Deze functie checkt of de regel daadwerkelijk een variant is.
    Daarna checkt deze functie of het een insertie/deletie is en telt dan
    de insertie/deletie op bij een totaal aan inserties/deleties.
    Als het geen indel is dan wordt de functie mutatie_checken(mutaties, line)
    aangeroepen.
    :param line: Line uit de vcf file
    :param deleties: Integer van het totaal aantal deleties
    :param inserties: Integer van het totaal aantal inserties
    :param mutaties: Dictionary van alle mutaties met hun totaal per mutatie.
    :return: Updated aantal inserties en deleties, en updated mutaties
    dictionary
    """
    if not line.startswith("#"):
        line = line.split()
        if len(line[3]) > len(line[4]):
            inserties += 1
        elif len(line[3]) < len(line[4]):
            deleties += 1
        else:
            mutaties = mutatie_checken(mutaties, line)
    return deleties, inserties, mutaties


def mutatie_checken(mutaties, variant):
    """
    Deze functie checkt per variant die geen indel is wat voor mutatie het is.
    Het soort mutatie wordt in de dictionary opgeteld.
    :param mutaties: Dictionary van alle mutaties met hun totaal per mutatie.
    :param variant: Line uit de vcf file die egen indel is
    :return: Updated mutaties dictionary
    """
    if variant[3] == "A":
        if variant[4] == "T":
            mutaties["A->T"] += 1
        elif variant[4] == "C":
            mutaties["A->C"] += 1
        elif variant[4] == "G":
            mutaties["A->G"] += 1
    elif variant[3] == "T":
        if variant[4] == "A":
            mutaties["T->A"] += 1
        elif variant[4] == "C":
            mutaties["T->C"] += 1
        elif variant[4] == "G":
            mutaties["T->G"] += 1
    elif variant[3] == "C":
        if variant[4] == "A":
            mutaties["C->A"] += 1
        elif variant[4] == "G":
            mutaties["C->G"] += 1
        elif variant[4] == "T":
            mutaties["C->T"] += 1
    elif variant[3] == "G":
        if variant[4] == "A":
            mutaties["G->A"] += 1
        elif variant[4] == "C":
            mutaties["G->C"] += 1
        elif variant[4] == "T":
            mutaties["G->T"] += 1
    return mutaties


def gcd(deleties, inserties):
    """
    Deze functie berekend de Greatest Common Deviser om uiteindelijk de
    verhouding uit te rekenen.
    :param deleties: Aantal mutaties
    :param inserties: Aantal inserties
    :return: de GCD
    """
    if inserties == 0:
        return deleties
    else:
        return gcd(inserties, deleties % inserties)


def functies(infile, outfile):
    """
    Deze functie is een soort main(). Er wordt door de VCF file gelooped om de
    varianten te analyseren. De informatie wordt weggeschreven naar een bestand.
    :param infile: de naam van de VCF file
    :param outfile: de naam van het output bestand
    :return: NVT
    """
    deleties, inserties = 0, 0
    mutaties = {"A->T": 0, "A->C": 0, "A->G": 0, "T->A": 0, "T->C": 0,
                "T->G": 0, "C->A": 0, "C->G": 0, "C->T": 0, "G->A": 0,
                "G->C": 0, "G->T": 0}
    with open(infile, "r") as vcf_file:
        while True:
            line = vcf_file.readline()
            if not line:
                break
            deleties, inserties, mutaties = \
                variant(line.strip("\n"), deleties, inserties,
                        mutaties)
            GCD = gcd(deleties, inserties)
            with open(outfile, "w") as analyse_file:
                try:
                    analyse_file.write(
                        "Aantal deleties ten opzichte van het referentiegenoom: " +
                        str(deleties) +
                        "\nAantal inserties ten opzichte van het "
                        "referentiegenoom: " + str(inserties) +
                        "\nVerhouding deleties/inserties = "
                        + str(int(deleties / GCD)) + ":" +
                        str(int(inserties / GCD)) + "\nMutatie rate: \n")
                    for mutatie in mutaties:
                        analyse_file.write(mutatie + " : " +
                                           str(mutaties[mutatie]) + "\n")
                except ZeroDivisionError:
                    continue

            analyse_file.close()

if __name__ == '__main__':
    functies(snakemake.input.vcf, snakemake.output.analyse_file)

