Deze pipeline laat een aantal python en shell scripts achter elkaar runnen om uiteindelijk een consensus sequentie van de geïnfecteerde sluipwesp te genereren. Er worden zelfgeschreven scripts gebruikt, maar ook bestaande applicaties. 
Voor de analyse van de varianten in het DNA van de geïnfecteerde sluipwesp is een VCF file gebruikt. Deze VCF file is eerder in een rule aangemaakt.

De pipeline wordt aangeroepen met de volgende commandline:
snakemake --snakefile Snakefile_s1113710

Als de gebruiker een specifieke uitvoerdirectory wilt voor de uitvoerbestanden kan de pipeline worden aangeroepen met de volgende commandline:
snakemake -d runjasmine --snakefile Snakefile_1113710
Alle uitvoerbestanden worden dan in de directory runjasmine/ geplaatst. 


QC rapports: bngsa_nietinfected_1.QC, bngsa_nietinfected_2.QC, bngsa_nietinfected_1_trimmed.QC, bngsa_nietinfected_2_trimmed.QC
Trimmed reads file: bngsa_nietinfected_1_trimmed.fastq, bngsa_nietinfected_2_trimmed.fastq
SAM bestand: alignment.sam
VCF bestanden: alignment_varianten.vcf.gz, alignment_varianten.vcf
Consensus sequentie: sluipwesp_consensus.fa
Analyse rapport: infected_sluipwest_analyse.txt

