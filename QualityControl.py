try:
    import sys
    print("Script naam: " + sys.argv[0] + "\nParameter: " + sys.argv[1] +
          "\nDit script genereert per invoerbestand de volgende output dat in een "
          "tekstbestand wordt weggeschreven: \n1. Het aantal reads"
          "\n2. De gemiddelde lengte van de reads \n3. De minimumlengte (oftewel: de kortste) van de reads"
          "\n4. Maximumlengte (oftewel: de langste) van de reads "
          "\n5. De GC content (oftewel: het percentage G+C)"
          "\n6. De GC content per positie in de read "
          "\n\nDit script kan alleen aangeroepen worden met snakemake met de volgende commandlines: "
          "snakemake --snakefile Snakemake_s1113710 rule QC !!OF!! snakemake --snakefile Snakemake_s1113710 rule QC_trimmed\n"
          "\nWanneer het script wordt aangeroepen met python3 volgt de volgende foutmelding: \n")
except IndexError:
    pass

def lengte(sequentie, totale_lengte):
    """
    Deze functie telt per read de lengte van de read bij de totale lengte van
    alle reads.
    :param sequentie: De sequentie van de read
    :param totale_lengte: Integer van de totale lengte
    :return: Updated totale_lengte
    """
    totale_lengte += len(sequentie)
    return totale_lengte


def informatie(sequentie, kortste_read, langste_read):
    """
    Deze functie houdt de langste en de kortste read bij. Wanneer de volgende
    read langer/korter is dan de kortste/langste read die al bekend is,
    Veranderd de kortste_read/langste_read in de huidige read in de loop.
    :param sequentie: De sequentie van de read
    :param kortste_read: Kortste read
    :param langste_read: Langste read
    :return: Kortste read, langste read
    """
    if len(sequentie) < kortste_read:
        kortste_read = len(sequentie)
    if len(sequentie) > langste_read:
        langste_read = len(sequentie)
    return kortste_read, langste_read


def GCaantal(sequentie, totaalGC):
    """
    Deze functie telt alle Gtjes en Ctjes per read op bij het totaal aantal
    Gtjes en Ctjes in alle reads.
    :param sequentie: De sequentie van de read
    :param totaalGC: Integer van het totaal aantal Gtjes en Ctjes
    :return: Updated totaalGC
    """
    countGC = 0
    for nucleotide in sequentie:
        if nucleotide == "G" or nucleotide == "C":
            countGC += 1
    totaalGC += countGC
    return totaalGC


def GCposities(sequentie, GC_per_positie, langste_read):
    """
    Deze functie loop door de sequentie van de read heen. Er is een dictionary
    die per positie laat zien hoeveel Gtjes en Ctjes geteld zijn. In de loop
    wordt dit aantal namelijk gewijzigd in de dictionary. De lengte van de
    langste read wordt gebruikt om de lengte van deze dictionary te bepalen.
    :param sequentie: De sequentie van de read
    :param GC_per_positie: Dictionairy met per positie het GC content
    :param langste_read: Lengte van de langste read.
    :return: Updated GC_per_positie dictionary
    """
    posities = 0
    if len(GC_per_positie) < langste_read:
        posities = langste_read - len(GC_per_positie)
    for positie in range(posities):
        GC_per_positie.append([0, 0])
    for position in range(len(sequentie)):
        try:
            GC_per_positie[position][1] += 1
            if sequentie[position] == "G" or sequentie[position] == "C":
                GC_per_positie[position][0] += 1
        except IndexError:
            continue
    return GC_per_positie

def aanroep(sequentie, kortste_read, langste_read, totale_lengte, totaalGC,
            GC_per_positie):
    """
    Deze functie roept alle informatie functies aan.
    :param sequentie: De sequentie van de read
    :param kortste_read: Kortste read
    :param langste_read: Langste read
    :param totale_lengte: Integer van de totale lengte
    :param totaalGC: Integer van het totaal aantal Gtjes en Ctjes
    :param GC_per_positie: Dictionairy met per positie het GC content
    :return: Alle parameters worden geupdate gereturned
    """
    totale_lengte = lengte(sequentie, totale_lengte)
    kortste_read, langste_read = informatie(sequentie, kortste_read,
                                            langste_read)
    totaalGC = GCaantal(sequentie, totaalGC)
    GC_per_positie = GCposities(sequentie, GC_per_positie, langste_read)
    return totale_lengte, kortste_read, langste_read, totaalGC, GC_per_positie

def wegschrijven(totale_lengte, reads, kortste_read, langste_read, totaalGC,
                 GC_per_positie, input_bestand, output_bestand):
    """
    Deze functie schrijft alle informatie weg naar een bestand.
    :param totale_lengte: Integer van de totale lengte
    :param reads: Aantal reads
    :param kortste_read: Kortste read
    :param langste_read: Langste read
    :param totaalGC: Integer van het totaal aantal Gtjes en Ctjes
    :param GC_per_positie: Dictionairy met per positie het GC content
    :param input_bestand: Naam input bestand
    :param output_bestand: Naam output bestand
    :return:
    """
    gem_lengte = round((totale_lengte / reads), 2)
    gem_GC = round((totaalGC / totale_lengte * 100), 2)

    fastqc = open(output_bestand, "w")
    fastqc.write("Quality rapport voor: " + input_bestand + "\n\n" +
                 "Lengte kortste read: " + str(kortste_read) + "\n" +
                 "Lengte langste read: " + str(langste_read) + "\n" +
                 "Gemiddelde lengte van de reads: " + str(gem_lengte) + "\n" +
                 "Aantal reads: " + str(reads) + "\n" +
                 "Gemiddelde GC percentage: " + str(gem_GC) + "%" + "\n\n")

    with open(output_bestand, "a") as fastqc:
        for positie, percentage in enumerate(GC_per_positie):
            fastqc.write("Gemiddelde GC percentage op positie " +
                         str(positie + 1) + ": " +
                         str(round((percentage[0] / percentage[1] * 100),
                                   2)) + "%" + "\n")

def functies(infile, outfile):
    """
    Deze functie is een soort main(). Hier wordt er door de reads file gelooped
    om alle informatie te verkrijgen van de file en per read.
    :param infile: Naam input bestand
    :param outfile: Naam output bestand
    :return: NVT
    """
    kortste_read, langste_read = 1000, 0
    totale_lengte, totaalGC, reads = 0, 0, 0
    GC_per_positie = []

    with open(infile, "r") as file:
        while True:
            header = file.readline()
            sequentie = file.readline()
            tab = file.readline()
            phredscore = file.readline()
            if not sequentie:
                break
            reads += 1
            totale_lengte, kortste_read, langste_read, totaalGC, \
            GC_per_positie = aanroep(sequentie.strip("\n"), kortste_read,
                                     langste_read, totale_lengte, totaalGC,
                                     GC_per_positie)
    wegschrijven(totale_lengte, reads, kortste_read, langste_read, totaalGC,
                 GC_per_positie, infile, outfile)

if __name__ == '__main__':
    functies(snakemake.input.reads, snakemake.output.qc)
