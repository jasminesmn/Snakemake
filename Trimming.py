try:
    import sys
    print("Script naam: " + sys.argv[0] + "\nParameter: " + sys.argv[1] +
          "\nDit script trimt de reads door middel van de methode: sliding window."
          "\nHiervoor is een cutoff van 28 gehanteerd, ook worden reads met een "
          "gemiddelde phredscore lager dan 28 niet meegenomen in de getrimde reads."
          "\nGetrimde reads die korter zijn dan 20 worden ook niet meegenomen in de"
          "getrimde reads file. "
          "\n\nDit script kan alleen aangeroepen worden met snakemake met de volgende commandline: "
          "snakemake --snakefile Snakemake_s1113710 rule trimming\n"
          "\nWanneer het script wordt aangeroepen met python3 volgt de volgende foutmelding: \n")
except IndexError:
    pass


def trimming(sequentie, phredscore):
    """
    De read wordt hier eerst forward getrimd, daarna backwards.
    :param sequentie: De nucleotide sequentie van de read
    :param phredscore: De phredscore per nucleotide in de read (in ASCII)
    :return: Geheel getrimde sequentie, phredscores in ASCII, en phredscore in
    integers
    """
    sequentie, phredscore, phredscores = trimforward(sequentie, phredscore)
    sequentie, phredscore, phredscores = trimbackwards(sequentie, phredscore,
                                                       phredscores)
    return sequentie, phredscore, phredscores


def trimforward(sequentie, phredscore):
    """
    Er wordt een cuttof gebruikt van 28. Dit omdat de gemiddelde phredscore van
    de reads ongeveer 28 is. Gesequencede nucleotiden met een phredscore lager
    dan 28 zullen niet goed genoeg zijn. Er is gewerkt met een window van 3 om
    met precieze gemiddelden te werken.
    :param sequentie: De nucleotide sequentie van de read
    :param phredscore: De phredscore per nucleotide in de read (in ASCII)
    :return: Getrimde sequentie, phredscores in ASCII en phredscore in integers
    """
    phredscores, sequentielijst, phredscorelijst = [], list(sequentie), \
                                                   list(phredscore)
    cutoff, window, position, trim = 28, 3, 0, 0
    for score in phredscore:
        phredscores.append(ord(score) - 33)
    for scores in range(len(phredscores)):
        if (sum(phredscores[position:scores + window]) / window) >= cutoff:
            trim = position
            break
        position += 1
    if trim != 0:
        del sequentielijst[:trim]
        del phredscorelijst[:trim]
        del phredscores[:trim]
    return "".join(sequentielijst), "".join(phredscorelijst), phredscores


def trimbackwards(sequentie, phredscore, phredscores):
    """
    Er wordt een cuttof gebruikt van 28. Dit omdat de gemiddelde phredscore van
    de reads ongeveer 28 is. Gesequencede nucleotiden met een phredscore lager
    dan 28 zullen niet goed genoeg zijn. Er is gewerkt met een window van 3 om
    met precieze gemiddelden te werken.
    :param sequentie: De nucleotide sequentie van de read
    :param phredscore: De phredscore per nucleotide in de read (in ASCII)
    :param phredscores: De phredscore per nucleotide in de read
    :return: Getrimde sequentie, phredscores in ASCII en phredscore in integers
    """
    sequentielijst, phredscorelijst = list(sequentie), list(phredscore)
    cutoff, window, position, trim = 28, 3, len(phredscores), 0
    for scores in range(len(phredscores), 0, -1):
        if (sum(phredscores[scores - window:position]) / window) >= cutoff:
            break
        trim += 1
        position -= 1
    if trim != 0:
        del sequentielijst[-trim:]
        del phredscorelijst[-trim:]
        del phredscores[-trim:]
    return "".join(sequentielijst), "".join(phredscorelijst), phredscores


def functies(infile_fw, infile_rv, outfile_fw, outfile_rv):
    """
    Deze functie is een soort main(). Hier wordt door de fastq files gelooped
    om per read te trimmen.
    :param infile_fw: Naam forward fastq file
    :param infile_rv: Naam reverse fastq file
    :param outfile_fw: Naam forward fastq file getrimd
    :param outfile_rv: Naam reverse fastq file getrimd
    :return: NVT
    """
    trimmed_fw_file = open(outfile_fw, "w")
    trimmed_rv_file = open(outfile_rv, "w")
    file_fw, file_rv = open(infile_fw, "r"), open(infile_rv, "r")

    while True:
        header_fw, header_rv = file_fw.readline(), file_rv.readline()
        sequentie_fw, sequentie_rv = file_fw.readline(), file_rv.readline()
        tab_fw, tab_rv = file_fw.readline(), file_rv.readline()
        phredscore_fw, phredscore_rv = file_fw.readline(), file_rv.readline()

        if not sequentie_fw and not sequentie_rv:
            break
        sequentie_fw, phredscore_fw, phredscores_fw = \
            trimming(sequentie_fw.strip("\n"), phredscore_fw.strip("\n"))
        sequentie_rv, phredscore_rv, phredscores_rv = \
            trimming(sequentie_rv.strip("\n"), phredscore_rv.strip("\n"))
        try:
            gem_phred_fw = sum(phredscores_fw) / len(phredscores_fw)
            gem_phred_rv = sum(phredscores_rv) / len(phredscores_rv)
        except ZeroDivisionError:
            continue
        """Het is belangrijk dat er genoeg variatie kan zitten in een 
        sequentie. Wanneer de read een sequentie heeft van langer dan 20 
        nucleotiden zijn er 4^20 mogelijke variaties in de sequentie. Dit staat
        gelijk aan 1099511627776 mogelijke variaties. In verhouding met de 
        grote van het referentie genoom zijn dit ruim voldoende variaties. 
        Hierom is gekozen voor een drempelwaarde van 20. 
        Omdat er niet alleen bij het begin en eind van de reads slechte phred-
        scores zijn is er gekozen om ook reads met een gemiddelde phredscore
        lager dan 28 uit te sluiten. Omdat er over het algemeen over een
        gemiddelde phredscore van 28 gesproken kan worden, is er voor deze 
        drempelwaarde gekozen. Het staat dus gelijk aan de cutoff van de 
        sliding window trimming methode. """
        if len(sequentie_fw) > 20 and len(sequentie_rv) > 20\
                and gem_phred_rv > 28 and gem_phred_fw > 28:
            trimmed_fw_file.write(header_fw + sequentie_fw + "\n" + tab_fw +
                                  phredscore_fw + "\n")
            trimmed_rv_file.write(header_rv + sequentie_rv + "\n" + tab_rv +
                                  phredscore_rv + "\n")


if __name__ == '__main__':
    functies(snakemake.input.fw, snakemake.input.rv, snakemake.output.fw,
             snakemake.output.rv)


