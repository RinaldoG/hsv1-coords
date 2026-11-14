#!/usr/bin/env python3
"""Build data/derived/aliases_literature.tsv, the literature half of the HSV-1 gene
name dictionary (milestone M1).

Every alias carries a citation. A name with no citation does not appear here: it is
listed instead under "Unverified candidates" in data/derived/alias_notes.md.

Three classes of source are used, and the citation column says which:

  * database records held locally and checksummed in data/raw/MANIFEST.tsv
    (NC_001806.2, BK012101.1, JN555585.1, the Zenodo GTF),
  * UniProtKB/Swiss-Prot entries for Human herpesvirus 1 strain 17 syn+ (taxon 10299),
    whose protein "alternative names" are curated synonyms,
  * primary and review literature, each with a PMID or DOI.

Run: python3 src/build_aliases.py
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "data", "derived", "aliases_literature.tsv")

# --------------------------------------------------------------------------- #
# Citation registry. Key -> (citation text, citation identifier).
# --------------------------------------------------------------------------- #

CIT = {
    # ---- database records (local copies checksummed in data/raw/MANIFEST.tsv) ----
    "REF": (
        "NCBI RefSeq NC_001806.2, Human alphaherpesvirus 1 strain 17, complete genome, "
        "provisional RefSeq record, last modified 13 August 2018; /gene, /product and "
        "/note qualifiers of the feature table",
        "NC_001806.2",
    ),
    "TPA": (
        "GenBank BK012101.1 (TPA_asm), Human alphaherpesvirus 1 strain 17, complete "
        "genome, third-party annotation deposited 7 November 2020 by Whisnant AW, "
        "Juerges CS, Erhard F and Doelken L; /gene, /product and /note qualifiers of "
        "the CDS feature table",
        "BK012101.1",
    ),
    "GTF": (
        "Whisnant AW, Juerges CS, Hennig T, Erhard F, Doelken L et al. (2020) "
        "HSV1_Annotation.gtf, dataset accompanying Nat Commun 11:2038, Zenodo; "
        "gene_id and transcript_id attributes, seqname JN555585",
        "10.5281/zenodo.3724023",
    ),
    "W20": (
        "Whisnant AW, Jurges CS, Hennig T, Wyler E, Prusty B, Rutkowski AJ et al. "
        "(2020) Integrative functional genomics decodes herpes simplex virus 1. "
        "Nat Commun 11:2038",
        "10.1038/s41467-020-15992-5",
    ),
    "DAVISON2011": (
        "Davison AJ (2011) Evolution of sexually transmitted and sexually "
        "transmissible human herpesviruses. Ann N Y Acad Sci 1230:E37-E49 "
        "(the submitting reference of JN555585.1 and NC_001806.2, whose COMMENT "
        "states that a standard protein nomenclature has been applied so that "
        "orthologs have the same name in all herpesviruses)",
        "PMID:22417106",
    ),
    # ---- UniProtKB/Swiss-Prot entries (taxon 10299), read via the REST API ----
    "UP:P06492": ("UniProtKB/Swiss-Prot P06492 (VP16_HHV11), protein alternative names", "P06492"),
    "UP:P08393": ("UniProtKB/Swiss-Prot P08393 (ICP0_HHV11), protein alternative names", "P08393"),
    "UP:P08392": ("UniProtKB/Swiss-Prot P08392 (ICP4_HHV11), protein alternative names", "P08392"),
    "UP:P10238": ("UniProtKB/Swiss-Prot P10238 (ICP27_HHV11), protein alternative names", "P10238"),
    "UP:P04485": ("UniProtKB/Swiss-Prot P04485 (ICP22_HHV11), protein alternative names", "P04485"),
    "UP:P03170": ("UniProtKB/Swiss-Prot P03170 (ICP47_HHV11), protein alternative names", "P03170"),
    "UP:P04487": ("UniProtKB/Swiss-Prot P04487 (RNB_HHV11), protein alternative names", "P04487"),
    "UP:P04290": ("UniProtKB/Swiss-Prot P04290 (UL13_HHV11), protein alternative names", "P04290"),
    "UP:P10231": ("UniProtKB/Swiss-Prot P10231 (TEG5_HHV11), protein alternative names", "P10231"),
    "UP:P10230": ("UniProtKB/Swiss-Prot P10230 (TEG1_HHV11), protein alternative names", "P10230"),
    "UP:P10210": ("UniProtKB/Swiss-Prot P10210 (SCAF_HHV11), protein alternative names", "P10210"),
    "UP:P36313": ("UniProtKB/Swiss-Prot P36313 (ICP34_HHV11), protein alternative names", "P36313"),
    "UP:P08543": ("UniProtKB/Swiss-Prot P08543 (RIR1_HHV11), gene and protein names", "P08543"),
    "UP:P04296": ("UniProtKB/Swiss-Prot P04296 (DNBI_HHV11), gene names", "P04296"),
    "UP:P10225": ("UniProtKB/Swiss-Prot P10225 (SHUT_HHV11), protein alternative names", "P10225"),
    "UP:P10226": ("UniProtKB/Swiss-Prot P10226 (PAP_HHV11), protein alternative names", "P10226"),
    "UP:P10193": ("UniProtKB/Swiss-Prot P10193 (OBP_HHV11), protein alternative names", "P10193"),
    "UP:P10192": ("UniProtKB/Swiss-Prot P10192 (HEPA_HHV11), protein alternative names", "P10192"),
    "UP:P10186": ("UniProtKB/Swiss-Prot P10186 (UNG_HHV11), protein alternative names", "P10186"),
    "UP:P68331": ("UniProtKB/Swiss-Prot P68331 (GK_HHV11), protein alternative names", "P68331"),
    "UP:P04295": ("UniProtKB/Swiss-Prot P04295 (TRM3_HHV11), protein alternative names", "P04295"),
    "UP:P0DTH5": ("UniProtKB/Swiss-Prot P0DTH5 (KITH_HHV11), gene names", "P0DTH5"),
    "UP:O09800": ("UniProtKB/Swiss-Prot O09800 (GN_HHV11), gene names", "O09800"),
    "UP:O09802": ("UniProtKB/Swiss-Prot O09802 (US8A_HHV11), gene and protein names", "O09802"),
    "UP:P10234": ("UniProtKB/Swiss-Prot P10234 (DUT_HHV11), protein alternative names", "P10234"),
    "UP:P10229": ("UniProtKB/Swiss-Prot P10229 (EV45_HHV11), protein alternative names", "P10229"),
    "UP:P06481": ("UniProtKB/Swiss-Prot P06481 (US9_HHV11), protein alternative names", "P06481"),
    "UP:P04488": ("UniProtKB/Swiss-Prot P04488 (GE_HHV11), protein alternative names", "P04488"),
    "UP:P06484": ("UniProtKB/Swiss-Prot P06484 (GG_HHV11), protein alternative names", "P06484"),
    "UP:P06491": ("UniProtKB/Swiss-Prot P06491 (MCP_HHV11), protein name and short name", "P06491"),
    "UP:P10185": ("UniProtKB/Swiss-Prot P10185 (GL_HHV11), protein short name", "P10185"),
    "UP:P04294": ("UniProtKB/Swiss-Prot P04294 (AN_HHV11), protein recommended name", "P04294"),
    "UP:P10236": ("UniProtKB/Swiss-Prot P10236 (PRIM_HHV11), protein recommended name", "P10236"),
    # ---- primary literature ----
    "EVERETT1984": (
        "Everett RD (1984) Trans activation of transcription by herpes virus products: "
        "requirement for two HSV-1 immediate-early polypeptides for maximum activity. "
        "EMBO J 3:3135-3141",
        "PMID:6098466",
    ),
    "PATERSON1990": (
        "Paterson T, Everett RD (1990) A prominent serine-rich region in Vmw175, the "
        "major transcriptional regulator protein of herpes simplex virus type 1, is "
        "not essential for virus growth in tissue culture. J Gen Virol 71:1775-1783",
        "PMID:2167931",
    ),
    "COUSENS1989": (
        "Cousens DJ, Greaves R, Goding CR, O'Hare P (1989) The C-terminal 79 amino "
        "acids of the herpes simplex virus regulatory protein, Vmw65, efficiently "
        "activate transcription in yeast and mammalian cells in chimaeric DNA-binding "
        "proteins. EMBO J 8:2337-2342",
        "PMID:2676518",
    ),
    "ACE1989": (
        "Ace CI, McKee TA, Ryan JM, Cameron JM, Preston CM (1989) Construction and "
        "characterization of a herpes simplex virus type 1 mutant unable to transinduce "
        "immediate-early gene expression. J Virol 63:2260-2269",
        "PMID:2539517",
    ),
    "EVERETTMAUL1994": (
        "Everett RD, Maul GG (1994) HSV-1 IE protein Vmw110 causes redistribution of "
        "PML. EMBO J 13:5062-5069",
        "PMID:7957072",
    ),
    "PRESTON1992": (
        "Preston VG, Rixon FJ, McDougall IM, McGregor M, al Kobaisi MF (1992) "
        "Processing of the herpes simplex virus assembly protein ICP35 near its "
        "carboxy terminal end requires the product of the whole of the UL26 reading "
        "frame. Virology 186:87-98",
        "PMID:1309284",
    ),
    "LIUROIZMAN1991": (
        "Liu FY, Roizman B (1991) The herpes simplex virus 1 gene encoding a protease "
        "also contains within its coding domain the gene encoding the more abundant "
        "substrate. J Virol 65:5149-5156",
        "PMID:1654435",
    ),
    "CARTER1996": (
        "Carter KL, Roizman B (1996) The promoter and transcriptional unit of a novel "
        "herpes simplex virus 1 alpha gene are contained in, and encode a protein in "
        "frame with, the open reading frame of the alpha 22 gene. J Virol 70:172-178",
        "PMID:8523523",
    ),
    "JOVASEVIC2010": (
        "Jovasevic V, Roizman B (2010) The novel HSV-1 US5-1 RNA is transcribed off a "
        "domain encoding US5, US4, US3, US2 and alpha22. Virol J 7:103",
        "PMID:20492679",
    ),
    "HUBENTHAL1987": (
        "Hubenthal-Voss J, Starr L, Roizman B (1987) The herpes simplex virus origins "
        "of DNA synthesis in the S component are each contained in a transcribed open "
        "reading frame. J Virol 61:3349-3355",
        "PMID:2822945",
    ),
    # ---- authoritative reference works that adjudicate rather than merely use ----
    "ICTV2021": (
        "Gatherer D, Depledge DP, Hartley CA, Szpara ML, Vaz PK, Benko M, Brandt CR, "
        "Bryant NA, Dastjerdi A, Doszpoly A, Gompels UA, Inoue N, Jarosinski KW, Kaul R, "
        "Lacoste V, Norberg P, Origgi FC, Orton RJ, Pellett PE, Schmid DS, Spatz SJ, "
        "Stewart JP, Trimpert J, Waltzek TB, Davison AJ, ICTV Report Consortium. ICTV "
        "Report chapter: Orthoherpesviridae (posted August 2021, updated May 2023), "
        "Table 2, shared genes in members of the family Orthoherpesviridae. Summary "
        "paper: ICTV Virus Taxonomy Profile: Herpesviridae 2021, J Gen Virol 102:001673",
        "10.1099/jgv.0.001673",
    ),
    "DAVISON2010": (
        "Davison AJ (2010) Herpesvirus systematics. Vet Microbiol 143:52-69. The source "
        "of the standardised-protein-name policy that NC_001806.2 and JN555585.1 "
        "implement: 'My tack has been to retain original gene names so that workers can "
        "find their way round an accession, and to apply standard names to the encoded "
        "proteins'",
        "PMID:20346601",
    ),
    "FIELDS6CH59": (
        "Pellett PE, Roizman B (2013) Herpesviridae. In: Knipe DM, Howley PM, eds. Fields "
        "Virology, 6th edition, chapter 59, from p. 1802. Philadelphia: Lippincott "
        "Williams and Wilkins. Set ISBN 9781451105636",
        "ISBN:9781451105636",
    ),
    "FIELDS6CH60": (
        "Roizman B, Knipe DM, Whitley RJ (2013) Herpes Simplex Viruses. In: Knipe DM, "
        "Howley PM, eds. Fields Virology, 6th edition, chapter 60, pp. 1823-1897. "
        "Philadelphia: Lippincott Williams and Wilkins. Set ISBN 9781451105636. The "
        "classical gene inventory is a prose list under the heading 'Open Reading "
        "Frames', not a numbered table",
        "ISBN:9781451105636",
    ),
    # ---- convention-origin papers ----
    "HONESS1973": (
        "Honess RW, Roizman B (1973) Proteins specified by herpes simplex virus. XI. "
        "Identification and relative molar rates of synthesis of structural and "
        "nonstructural herpes virus polypeptides in the infected cell. J Virol "
        "12:1347-1365 (the origin of the ICP numbering, by SDS-PAGE mobility; ICP stands "
        "for infected cell polypeptide)",
        "PMID:4357511",
    ),
    "SPEAR1972": (
        "Spear PG, Roizman B (1972) Proteins specified by herpes simplex virus. V. "
        "Purification and structural proteins of the herpesvirion. J Virol 9:143-159 "
        "(the origin of the virion-protein band numbering)",
        "PMID:4110104",
    ),
    "GIBSON1972": (
        "Gibson W, Roizman B (1972) Proteins specified by herpes simplex virus. 8. "
        "Characterization and composition of multiple capsid forms of subtypes 1 and 2. "
        "J Virol 10:1044-1052 (first assignment of the capsid band numbers 5, 19, 21, "
        "22a, 23 and 24)",
        "PMID:4344252",
    ),
    "GIBSON1974": (
        "Gibson W, Roizman B (1974) Proteins specified by herpes simplex virus. Staining "
        "and radiolabeling properties of B capsid and virion proteins in polyacrylamide "
        "gels. J Virol 13:155-165",
        "PMID:4129836",
    ),
    "PRESTON1978": (
        "Preston VG, Davison AJ, Marsden HS, Timbury MC, Subak-Sharpe JH, Wilkie NM "
        "(1978) Recombinants between herpes simplex virus types 1 and 2: analyses of "
        "genome structures and expression of immediate early polypeptides. J Virol "
        "28:499-517 (earliest located use of the Glasgow Vmw convention)",
        "PMID:214575",
    ),
    "MCGEOCH1985": (
        "McGeoch DJ, Dolan A, Donald S, Rixon FJ (1985) Sequence determination and "
        "genetic content of the short unique region in the genome of herpes simplex virus "
        "type 1. J Mol Biol 181:1-13 (the paper that introduced the US1 to US12 numbering)",
        "PMID:2984429",
    ),
    "MCGEOCH1988": (
        "McGeoch DJ, Dalrymple MA, Davison AJ, Dolan A, Frame MC, McNab D, Perry LJ, "
        "Scott JE, Taylor P (1988) The complete DNA sequence of the long unique region in "
        "the genome of herpes simplex virus type 1. J Gen Virol 69:1531-1574 (the paper "
        "that introduced the UL1 to UL56 numbering)",
        "PMID:2839594",
    ),
    "MCGEOCH1986SR": (
        "McGeoch DJ, Dolan A, Donald S, Brauer DH (1986) Complete DNA sequence of the "
        "short repeat region in the genome of herpes simplex virus type 1. Nucleic Acids "
        "Res 14:1727-1745",
        "PMID:3005980",
    ),
    "PERRY1988": (
        "Perry LJ, McGeoch DJ (1988) The DNA sequences of the long repeat region and "
        "adjoining parts of the long unique region in the genome of herpes simplex virus "
        "type 1. J Gen Virol 69:2831-2846",
        "PMID:2846760",
    ),
    # ---- per-name primary papers ----
    "CAMPBELL1984": (
        "Campbell ME, Palfreyman JW, Preston CM (1984) Identification of herpes simplex "
        "virus DNA sequences which encode a trans-acting polypeptide responsible for "
        "stimulation of immediate early transcription. J Mol Biol 180:1-19",
        "PMID:6096556",
    ),
    "MCKNIGHT1987": (
        "McKnight JL, Pellett PE, Jenkins FJ, Roizman B (1987) Characterization and "
        "nucleotide sequence of two herpes simplex virus 1 genes whose products modulate "
        "alpha-trans-inducing factor-dependent activation of alpha genes. J Virol "
        "61:992-1001",
        "PMID:3029433",
    ),
    "PELLETT1985": (
        "Pellett PE, McKnight JL, Jenkins FJ, Roizman B (1985) Nucleotide sequence and "
        "predicted amino acid sequence of a protein encoded in a small herpes simplex "
        "virus DNA fragment capable of trans-inducing alpha genes. Proc Natl Acad Sci USA "
        "82:5870-5874",
        "PMID:2994050",
    ),
    "WAGNER1988": (
        "Wagner EK, Devi-Rao G, Feldman LT, Dobson AT, Zhang YF, Flanagan WM, Stevens JG "
        "(1988) Physical characterization of the herpes simplex virus latency-associated "
        "transcript in neurons. J Virol 62:1194-1202",
        "PMID:2831380",
    ),
    "STEVENS1987": (
        "Stevens JG, Wagner EK, Devi-Rao GB, Cook ML, Feldman LT (1987) RNA complementary "
        "to a herpesvirus alpha gene mRNA is prominent in latently infected neurons. "
        "Science 235:1056-1059",
        "PMID:2434993",
    ),
    "READ1983": (
        "Read GS, Frenkel N (1983) Herpes simplex virus mutants defective in the "
        "virion-associated shutoff of host polypeptide synthesis and exhibiting abnormal "
        "synthesis of alpha (immediate early) viral polypeptides. J Virol 46:498-512",
        "PMID:6302315",
    ),
    "READ1993": (
        "Read GS, Karr BM, Knight K (1993) Isolation of a herpes simplex virus type 1 "
        "mutant with a deletion in the virion host shutoff gene and identification of "
        "multiple forms of the vhs (UL41) polypeptide. J Virol 67:7149-7160",
        "PMID:8230437",
    ),
    "HUTCHINSONGK1992": (
        "Hutchinson L, Goldsmith K, Snoddy D, Ghosh H, Graham FL, Johnson DC (1992) "
        "Identification and characterization of a novel herpes simplex virus glycoprotein, "
        "gK, involved in cell fusion. J Virol 66:5603-5609",
        "PMID:1323714",
    ),
    "HUTCHINSONGL1992": (
        "Hutchinson L, Browne H, Wargent V, Davis-Poynter N, Primorac S, Goldsmith K, "
        "Minson AC, Johnson DC (1992) A novel herpes simplex virus glycoprotein, gL, forms "
        "a complex with glycoprotein H (gH) and affects normal folding and surface "
        "expression of gH. J Virol 66:2240-2250",
        "PMID:1312629",
    ),
    "BAUCKE1979": (
        "Baucke RB, Spear PG (1979) Membrane proteins specified by herpes simplex viruses. "
        "V. Identification of an Fc-binding glycoprotein. J Virol 32:779-789",
        "PMID:229267",
    ),
    "FRAME1986": (
        "Frame MC, Marsden HS, McGeoch DJ (1986) Novel herpes simplex virus type 1 "
        "glycoproteins identified by antiserum against a synthetic oligopeptide from the "
        "predicted product of gene US4. J Gen Virol 67:745-751",
        "PMID:3007661",
    ),
    "BARKER1992": (
        "Barker DE, Roizman B (1992) The unique sequence of the herpes simplex virus 1 L "
        "component contains an additional translated open reading frame designated "
        "UL49.5. J Virol 66:562-566",
        "PMID:1309263",
    ),
    "PRESTONVG1984": (
        "Preston VG, Palfreyman JW, Dutia BM (1984) Identification of a herpes simplex "
        "virus type 1 polypeptide which is a component of the virus-induced ribonucleotide "
        "reductase. J Gen Virol 65:1457-1466",
        "PMID:6088678",
    ),
    "GOLDSTEIN1988": (
        "Goldstein DJ, Weller SK (1988) Herpes simplex virus type 1-induced ribonucleotide "
        "reductase activity is dispensable for virus growth and DNA synthesis: isolation "
        "and characterization of an ICP6 lacZ insertion mutant. J Virol 62:196-205",
        "PMID:2824847",
    ),
    "CHUNG1989": (
        "Chung TD, Wymer JP, Smith CC, Kulka M, Aurelian L (1989) Protein kinase activity "
        "associated with the large subunit of herpes simplex virus type 2 ribonucleotide "
        "reductase (ICP10). J Virol 63:3389-3398",
        "PMID:2545912",
    ),
    "KNOPF1987": (
        "Knopf CW (1987) The nucleotide sequence of the herpes simplex virus type 1 late "
        "gene ICP18.5 of strain Angelotti. Nucleic Acids Res 15:8109-8110",
        "PMID:2823222",
    ),
    "LIUROIZMAN1991B": (
        "Liu FY, Roizman B (1991) The promoter, transcriptional unit, and coding sequence "
        "of herpes simplex virus 1 family 35 proteins are contained within and in frame "
        "with the UL26 open reading frame. J Virol 65:206-212",
        "PMID:1845885",
    ),
    "MCGREGOR1996": (
        "McGregor F, Phelan A, Dunlop J, Clements JB (1996) Regulation of herpes simplex "
        "virus poly (A) site usage and the action of immediate-early protein IE63 in the "
        "early-late switch. J Virol 70:1931-1940",
        "PMID:8627719",
    ),
    "JACKERS1992": (
        "Jackers P, Defechereux P, Baudoux L, Lambert C, Massaer M, Merville-Louis MP, "
        "Rentier B, Piette J (1992) Characterization of regulatory functions of the "
        "varicella-zoster virus gene 63-encoded protein. J Virol 66:3899-3903",
        "PMID:1316489",
    ),
    "ELLIOTT1992": (
        "Elliott GD, Meredith DM (1992) The herpes simplex virus type 1 tegument protein "
        "VP22 is encoded by gene UL49. J Gen Virol 73:723-726",
        "PMID:1312128",
    ),
    "LAGUNOFF1994": (
        "Lagunoff M, Roizman B (1994) Expression of a herpes simplex virus 1 open reading "
        "frame antisense to the gamma(1)34.5 gene and transcribed by an RNA 3' coterminal "
        "with the unspliced latency-associated transcript. J Virol 68:6021-6028",
        "PMID:8057477",
    ),
    "RANDALL1997": (
        "Randall G, Lagunoff M, Roizman B (1997) The product of ORF O located within the "
        "domain of herpes simplex virus 1 genome transcribed during latent infection binds "
        "to and inhibits in vitro binding of infected cell protein 4 to its cognate DNA "
        "site. Proc Natl Acad Sci USA 94:10379-10384",
        "PMID:9294219",
    ),
    "OHARE1985": (
        "O'Hare P, Hayward GS (1985) Three trans-acting regulatory proteins of herpes "
        "simplex virus modulate immediate-early gene expression in a pathway involving "
        "positive and negative feedback regulation. J Virol 56:723-733",
        "PMID:2999428",
    ),
    "ZHANG1991": (
        "Zhang Y, Sirko DA, McKnight JL (1991) Role of herpes simplex virus type 1 UL46 "
        "and UL47 in alpha TIF-mediated transcriptional induction: characterization of "
        "three viral deletion mutants. J Virol 65:829-841",
        "PMID:1846201",
    ),
    "HAGMANN1997": (
        "Hagmann M, Georgiev O, Schaffner W (1997) The VP16 paradox: herpes simplex "
        "virus VP16 contains a long-range activation domain but within the natural "
        "multiprotein complex activates only from promoter-proximal positions. "
        "J Virol 71:5952-5962",
        "PMID:9223485",
    ),
    "MINAGAWA1994": (
        "Minagawa H, Tanaka S, Toh Y, Mori R (1994) Detection of herpes simplex virus "
        "type 1-encoded RNA by polymerase chain reaction: different pattern of viral RNA "
        "detection in latently infected murine trigeminal ganglia following in vitro or "
        "in vivo reactivation. J Gen Virol 75:647-650",
        "PMID:8126462",
    ),
    "MULLEN1995": (
        "Mullen MA, Gerstberger S, Ciufo DM, Mosca JD, Hayward GS (1995) Evaluation of "
        "colocalization interactions between the IE110, IE175, and IE63 transactivator "
        "proteins of herpes simplex virus within subcellular punctate structures. "
        "J Virol 69:476-491",
        "PMID:7983744",
    ),
    "ROBERTS1988": (
        "Roberts MS, Boundy A, O'Hare P, Pizzorno MC, Ciufo DM, Hayward GS (1988) Direct "
        "correlation between a negative autoregulatory response element at the cap site "
        "of the herpes simplex virus type 1 IE175 (alpha 4) promoter and a specific "
        "binding site for the IE175 (ICP4) protein. J Virol 62:4307-4320",
        "PMID:2845144",
    ),
    "SPECTOR1991": (
        "Spector D, Purves F, Roizman B (1991) Role of alpha-transinducing factor (VP16) "
        "in the induction of alpha genes within the context of viral genomes. "
        "J Virol 65:3504-3513",
        "PMID:1645782",
    ),
    "POON1995": (
        "Poon AP, Roizman B (1995) The phenotype in vitro and in infected cells of herpes "
        "simplex virus 1 alpha trans-inducing factor (VP16) carrying temperature-"
        "sensitive mutations introduced by substitution of cysteines. J Virol 69:7658-7667",
        "PMID:7494274",
    ),
    "WEINHEIMER1992": (
        "Weinheimer SP, Boyd BA, Durham SK, Resnick JL, O'Boyle DR (1992) Deletion of the "
        "VP16 open reading frame of herpes simplex virus type 1. J Virol 66:258-269",
        "PMID:1309245",
    ),
    "LIU1999": (
        "Liu Y, Gong W, Huang CC, Herr W, Cheng X (1999) Crystal structure of the "
        "conserved core of the herpes simplex virus transcriptional regulatory protein "
        "VP16. Genes Dev 13:1692-1703",
        "PMID:10398682",
    ),
    "EVERETT1995": (
        "Everett R, O'Hare P, O'Rourke D, Barlow P, Orr A (1995) Point mutations in the "
        "herpes simplex virus type 1 Vmw110 RING finger helix affect activation of gene "
        "expression, viral growth, and interaction with PML-containing nuclear "
        "structures. J Virol 69:7339-7344",
        "PMID:7474166",
    ),
    "ALLEN1997": (
        "Allen KE, Everett RD (1997) Mutations which alter the DNA binding properties of "
        "the herpes simplex virus type 1 transactivating protein Vmw175 also affect its "
        "ability to support virus replication. J Gen Virol 78:2913-2922",
        "PMID:9367379",
    ),
    "ZHANGMCKNIGHT1993": (
        "Zhang Y, McKnight JL (1993) Herpes simplex virus type 1 UL46 and UL47 deletion "
        "mutants lack VP11 and VP12 or VP13 and VP14, respectively, and exhibit altered "
        "viral thymidine kinase expression. J Virol 67:1482-1492",
        "PMID:8382306",
    ),
    "MCKNIGHT1994": (
        "McKnight JL, Doerr M, Zhang Y (1994) An 85-kilodalton herpes simplex virus type "
        "1 alpha trans-induction factor (VP16)-VP13/14 fusion protein retains the "
        "transactivation and structural properties of the wild-type molecule during virus "
        "infection. J Virol 68:1750-1757",
        "PMID:8107236",
    ),
    "SINGH1996": (
        "Singh M, Fraefel C, Bello LJ, Lawrence WC, Schwyzer M (1996) Identification and "
        "characterization of BICP27, an early protein of bovine herpesvirus 1 which may "
        "stimulate mRNA 3' processing. J Gen Virol 77:615-625",
        "PMID:8627249",
    ),
    "VILLANUEVA2021": (
        "Villanueva-Valencia JR, Tsimtsirakis E, Evilevitch A (2021) Role of HSV-1 Capsid "
        "Vertex-Specific Component (CVSC) and Viral Terminal DNA in Capsid Docking at the "
        "Nuclear Pore. Viruses 13:2515",
        "10.3390/v13122515",
    ),
    "JIN2014": (
        "Jin F, Li S, Zheng K, Zhuo C, Ma K, Chen M, Wang Q, Zhang P, Fan J, Ren Z, Wang "
        "Y (2014) Silencing herpes simplex virus type 1 capsid protein encoding genes by "
        "siRNA: a promising antiviral therapeutic approach. PLoS One 9:e96623",
        "10.1371/journal.pone.0096623",
    ),
    "JAMBUNATHAN2021": (
        "Jambunathan N, Clark CM, Musarrat F, Chouljenko VN, Rudd J, Kousoulas KG (2021) "
        "Two Sides to Every Story: Herpes Simplex Type-1 Viral Glycoproteins gB, gD, "
        "gH/gL, gK, and Cellular Receptors Function as Key Players in Membrane Fusion. "
        "Viruses 13:1849",
        "10.3390/v13091849",
    ),
    "IBRAHIM2026": (
        "Ibrahim R, Fidouh N, Bunel V, Bouscarat F, Issaurat P, Fernandez J, Joly V, "
        "Rioux C, Boutolleau D, Ghosn J (2026) Clinical and genotypic characteristics of "
        "antiviral-resistant herpes simplex virus 1 and 2 infections: a case series. "
        "BMC Infect Dis 26:1072",
        "10.1186/s12879-026-13310-w",
    ),
    "DUAN2019": (
        "Duan Y, Zeng J, Fan S, Liao Y, Feng M, Wang L, Zhang Y, Li Q (2019) Herpes "
        "Simplex Virus Type 1-Encoded miR-H2-3p Manipulates Cytosolic DNA-Stimulated "
        "Antiviral Innate Immune Response by Targeting DDX41. Viruses 11:756",
        "10.3390/v11080756",
    ),
    "BRAUN1984": (
        "Braun DK, Roizman B, Pereira L (1984) Characterization of post-translational "
        "products of herpes simplex virus gene 35 proteins binding to the surfaces of "
        "full capsids but not empty capsids. J Virol 49:142-153",
        "PMID:6317887",
    ),
    "UMBACH2008": (
        "Umbach JL, Kramer MF, Jurak I, Karnowski HW, Coen DM, Cullen BR (2008) "
        "MicroRNAs expressed by herpes simplex virus 1 during latent infection regulate "
        "viral mRNAs. Nature 454:780-783",
        "PMID:18596690",
    ),
    "CUI2006": (
        "Cui C, Griffiths A, Li G, Silva LM, Kramer MF, Gaasterland T, Wang XJ, Coen DM "
        "(2006) Prediction and identification of herpes simplex virus 1-encoded "
        "microRNAs. J Virol 80:5499-5508",
        "PMID:16699030",
    ),
    "SCHWYZER1994": (
        "Schwyzer M, Wirth UV, Vogt B, Fraefel C (1994) BICP22 of bovine herpesvirus 1 "
        "is encoded by a spliced 1.7 kb RNA which exhibits immediate early and late "
        "transcription kinetics. J Gen Virol 75:1703-1711",
        "PMID:8021599",
    ),
    "PERNG2002": (
        "Perng GC, Maguen B, Jin L, Mott KR, Kurylo J, BenMohamed L et al. (2002) A "
        "novel herpes simplex virus type 1 transcript (AL-RNA) antisense to the 5' end "
        "of the latency-associated transcript produces a protein in infected rabbits. "
        "J Virol 76:8003-8010",
        "PMID:12134005",
    ),
}

# --------------------------------------------------------------------------- #
# The dictionary. canonical -> list of (alias, alias_type, cit_key, inventory,
# confidence, notes).
#
# alias_type is one of: systematic, protein, functional, molecular_weight, other.
# inventory is one of: classical, whisnant2020, both.
# confidence is one of: established, probable, ambiguous.
# --------------------------------------------------------------------------- #

S, P, F, M, O = "systematic", "protein", "functional", "molecular_weight", "other"
EST, PROB, AMB = "established", "probable", "ambiguous"

ROWS = []


DECLINE = "decline"
RESOLVE = "resolve"
"""What the resolver does with this row.

Round 10: the `confidence` column was loaded into the dictionary and never read by any
code path, so four rows flagged `ambiguous`, with notes saying in terms that they are
not resolved here, resolved silently and confidently. `ICP10` is the HSV-2 name for the
positional homologue and returned UL39 with no signal; `AP` is attributable to either
UL26 or UL26.5 and returned UL26.

Blanket refusal would be wrong, because `ambiguous` covers two different situations: a
string that cannot be attributed (refuse) and a string that is correctly attributed but
sits one letter from another gene's name (answer, and say so). So each flagged row states
which it is, rather than the resolver inferring it from the prose of the note.
"""


def add(canonical, alias, atype, cit, inventory, conf, notes="", resolution=None):
    """Record one gene/alias pair.

    An `ambiguous` row must state its resolution explicitly. Round 11: `resolution`
    defaulted to `RESOLVE`, so the guard here could never fire and
    `test_every_ambiguous_row_states_what_to_do_with_it` could never fail. Marking a new
    row `ambiguous` made it resolve silently, which is the behaviour round 10 added this
    column to stop: six ambiguous rows still resolved, and the next one added would have
    too. The default is None now, so the choice has to be made at the call site.
    """
    if conf == AMB and resolution is None:
        raise ValueError(
            f"{alias} ({canonical}) is marked ambiguous and states no resolution. "
            f"Pass DECLINE if the string cannot be attributed to one gene, or RESOLVE "
            f"if it can and the ambiguity is about something else, such as a near-miss "
            f"with another gene's name.")
    ROWS.append((canonical, alias, atype, cit, inventory, conf, notes,
                 resolution or RESOLVE))


# --- repeat-region genes -------------------------------------------------- #
add("RL1", "RL1", S, "REF", "classical", EST, "systematic name in the repeat long region")
add("RL1", "ICP34.5", P, "REF", "both", EST, "RefSeq /product 'neurovirulence protein ICP34.5'")
add("RL1", "ICP34.5", P, "TPA", "whisnant2020", EST, "CDS /note 'icp34.5' on both TRL1 and IRL1")
add("RL1", "gamma(1)34.5", F, "UP:P36313", "classical", EST,
    "UniProt alternative name 'protein gamma(1)34.5'; Roizman kinetic-class naming")
add("RL1", "TRL1", S, "TPA", "whisnant2020", EST,
    "the 2020 re-annotation names the terminal-repeat copy separately; classical RL1 is "
    "annotated once, the re-annotation annotates TRL1 and IRL1")
add("RL1", "IRL1", S, "TPA", "whisnant2020", EST, "internal-repeat copy of RL1")

add("RL2", "RL2", S, "REF", "classical", EST, "")
add("RL2", "ICP0", P, "REF", "both", EST, "RefSeq /product 'ubiquitin E3 ligase ICP0'")
add("RL2", "ICP0", P, "TPA", "whisnant2020", EST, "CDS /note 'ICP0' on both TRL2 and IRL2")
add("RL2", "IE110", M, "UP:P08393", "classical", EST, "UniProt 'Immediate-early protein IE110'")
add("RL2", "IE110", M, "TPA", "whisnant2020", EST, "CDS /note 'IE110'")
add("RL2", "Vmw110", M, "UP:P08393", "classical", EST, "UniProt alternative name 'VMW110'")
add("RL2", "Vmw110", M, "EVERETTMAUL1994", "classical", EST,
    "title uses 'Vmw110' for the HSV-1 IE protein; Glasgow molecular-weight convention")
add("RL2", "Vmw110", M, "EVERETT1984", "classical", EST,
    "abstract glosses it directly: 'Vmw110 (ICP 0)'")
add("RL2", "alpha-0", F, "UP:P08393", "classical", EST,
    "UniProt 'Alpha-0 protein'; Roizman kinetic-class naming")
add("RL2", "TRL2", S, "TPA", "whisnant2020", EST, "terminal-repeat copy")
add("RL2", "IRL2", S, "TPA", "whisnant2020", EST, "internal-repeat copy")

add("RS1", "RS1", S, "REF", "classical", EST, "")
add("RS1", "ICP4", P, "REF", "both", EST, "RefSeq /product 'transcriptional regulator ICP4'")
add("RS1", "ICP4", P, "TPA", "whisnant2020", EST, "CDS /note 'icp4' on both IRS1 and TRS1")
add("RS1", "IE175", M, "UP:P08392", "classical", EST, "UniProt 'Transcriptional activator IE175'")
add("RS1", "IE175", M, "TPA", "whisnant2020", EST, "CDS /note 'ie175'")
add("RS1", "Vmw175", M, "PATERSON1990", "classical", EST,
    "title: 'Vmw175, the major transcriptional regulator protein of herpes simplex virus "
    "type 1'. Not carried by UniProt P08392, so the citation is the primary paper")
add("RS1", "Vmw175", M, "EVERETT1984", "classical", EST,
    "abstract glosses it directly: 'Vmw175 (ICP 4)'")
add("RS1", "alpha-4", F, "UP:P08392", "classical", EST, "UniProt 'Alpha-4 protein'")
add("RS1", "IRS1", S, "TPA", "whisnant2020", EST, "internal-repeat copy")
add("RS1", "TRS1", S, "TPA", "whisnant2020", EST, "terminal-repeat copy")

add("LAT", "LAT", O, "REF", "both", EST,
    "RefSeq /gene='LAT' and ncRNA /product='LAT'; annotated twice, once truncated as "
    "complement(<1..7569)")
add("LAT", "LAT_gene", O, "GTF", "whisnant2020", EST,
    "gene_id in the Zenodo GTF; the underscore form is a file-level identifier, not a "
    "literature name")
add("LAT", "LAT-ORF-O", O, "TPA", "whisnant2020", EST,
    "BK012101.1 /gene for the compound LAT/ORF-O/ORF-P locus")
add("LAT", "latency-associated transcript", F, "CUI2006", "classical", EST,
    "abstract gloss, hyphenated: 'the transcription start site of the "
    "latency-associated transcript (LAT)'")
add("LAT", "latency associated transcript", F, "UMBACH2008", "classical", EST,
    "abstract gloss, unhyphenated: 'a non-coding RNA, the latency associated transcript "
    "(LAT)'. Both spellings are in use and a dictionary that normalises on the hyphen "
    "will miss one of them")
add("IRL2.5", "pri-miRNA-H1", O, "TPA", "whisnant2020", EST,
    "BK012101.1 /note on the IRL2.5 feature. The LAT-region miRNAs themselves are not "
    "given canonical entries in this pass: see alias_notes.md for the scope decision")
add("IRL2.5", "LAT 0.7kb", O, "TPA", "whisnant2020", EST, "BK012101.1 /note on IRL2.5")
add("TRL2.5", "pri-miRNA-H1", O, "TPA", "whisnant2020", EST,
    "BK012101.1 /note on the TRL2.5 feature, the terminal-repeat copy of IRL2.5")
add("TRL2.5", "LAT 0.7kb", O, "TPA", "whisnant2020", EST, "BK012101.1 /note on TRL2.5")

add("ORF-O", "ORF-O", O, "TPA", "whisnant2020", EST,
    "BK012101.1 CDS, note 'orphan'; antisense to the 5' end of LAT")
add("ORF-P", "ORF-P", O, "TPA", "whisnant2020", EST, "BK012101.1 CDS, note 'orphan'")
add("AL-RNA", "AL-RNA", O, "PERNG2002", "classical", PROB,
    "named in the primary paper as a transcript antisense to the 5' end of LAT")
add("AL-RNA", "AL-RNA", O, "TPA", "whisnant2020", AMB,
    "BK012101.1 /note cites 'Perng et al., JVI 2002, not observed': the re-annotation "
    "records the name but did not detect the feature",
    # identity mapping: the name resolves to itself, so answering adds nothing and withholding helps nobody
    RESOLVE)
add("Ori-S", "Ori-S", O, "HUBENTHAL1987", "classical", EST,
    "cited by name in the BK012101.1 /note 'Hubenthal-Voss,Starr,Roizman JVI 1987'")
add("Ori-S", "Ori-S-TRS1", S, "TPA", "whisnant2020", EST, "BK012101.1 /gene for the compound locus")

# --- unique long region --------------------------------------------------- #
UL_SIMPLE = {
    # canonical: (RefSeq product/alias string, alias, alias_type, extra cit, note)
    "UL1": [("gL", P, "REF", "RefSeq /product 'envelope glycoprotein L'"),
            ("gL", P, "UP:P10185", "UniProt short name 'gL'"),
            ("glycoprotein L", P, "REF", "RefSeq /product, written out")],
    "UL2": [("UDG", P, "UP:P10186", "UniProt short name"),
            ("UNG", P, "UP:P10186", "UniProt alternative short name"),
            ("uracil-DNA glycosylase", F, "REF", "RefSeq /product")],
    "UL5": [("helicase-primase helicase subunit", F, "REF", "RefSeq /product")],
    "UL6": [("portal protein", P, "REF", "RefSeq /product 'capsid portal protein'")],
    "UL8": [("HEPA", P, "UP:P10192", "UniProt short name for the helicase-primase associated protein"),
            ("helicase-primase subunit", F, "REF", "RefSeq /product")],
    "UL9": [("OBP", P, "UP:P10193", "UniProt short name"),
            ("OriBP", P, "UP:P10193", "UniProt alternative short name"),
            ("origin-binding protein", F, "REF",
             "RefSeq /product 'DNA replication origin-binding helicase'")],
    "UL10": [("gM", P, "REF", "RefSeq /product 'envelope glycoprotein M'"),
             ("gM", P, "TPA", "CDS /note 'glycoprotein m'")],
    "UL12": [("alkaline nuclease", F, "UP:P04294", "UniProt recommended name"),
             ("DNase", F, "REF", "RefSeq /product 'deoxyribonuclease'")],
    "UL13": [("Vmw57", M, "UP:P04290", "UniProt alternative name 'Virion protein VMW57'")],
    "UL15": [("terminase large subunit", F, "UP:P04295", "UniProt alternative name"),
             ("DNA packaging terminase subunit 1", F, "REF", "RefSeq /product")],
    "UL18": [("VP23", P, "REF", "RefSeq /note 'VP23'"),
             ("VP23", P, "TPA", "CDS /note 'vp23'"),
             ("capsid triplex subunit 2", F, "REF", "RefSeq /product")],
    "UL19": [("VP5", P, "REF", "RefSeq /note 'VP5'"),
             ("VP5", P, "TPA", "CDS /note 'vp5'"),
             ("MCP", P, "UP:P06491", "UniProt short name for the major capsid protein"),
             ("major capsid protein", F, "REF", "RefSeq /product")],
    "UL22": [("gH", P, "REF", "RefSeq /product 'envelope glycoprotein H'"),
             ("gH", P, "TPA", "CDS /note 'glycoprotein h'")],
    "UL23": [("TK", F, "UP:P0DTH5", "UniProt primary gene name 'TK'"),
             ("thymidine kinase", F, "REF", "RefSeq /product")],
    "UL27": [("gB", P, "REF", "RefSeq /product 'envelope glycoprotein B'"),
             ("gB", P, "TPA", "CDS /note 'glycoprotein b'")],
    "UL28": [("ICP18.5", P, "TPA", "CDS /note 'icp18.5'"),
             ("DNA packaging terminase subunit 2", F, "REF", "RefSeq /product")],
    "UL29": [("ICP8", P, "REF", "RefSeq /note 'ICP8'"),
             ("ICP8", P, "UP:P04296", "UniProt gene synonym 'ICP8'"),
             ("DBP", P, "UP:P04296", "UniProt primary gene name"),
             ("major DNA-binding protein", F, "UP:P04296", "UniProt recommended name"),
             ("single-stranded DNA-binding protein", F, "REF", "RefSeq /product")],
    "UL30": [("DNA polymerase catalytic subunit", F, "REF", "RefSeq /product"),
             ("Pol", F, "TPA", "CDS /note 'DNA polymerase'")],
    "UL31": [("nuclear egress lamina protein", F, "REF", "RefSeq /product")],
    "UL34": [("nuclear egress membrane protein", F, "REF", "RefSeq /product")],
    "UL35": [("VP26", P, "REF", "RefSeq /note 'VP26'"),
             ("VP26", P, "TPA", "CDS /note 'vp26'"),
             ("small capsid protein", F, "REF", "RefSeq /product")],
    "UL36": [("VP1/2", P, "REF", "RefSeq /note 'VP1/2'"),
             ("large tegument protein", F, "REF", "RefSeq /product")],
    "UL38": [("VP19C", P, "REF", "RefSeq /note 'VP19C'"),
             ("VP19C", P, "TPA", "CDS /note 'vp19c'"),
             ("capsid triplex subunit 1", F, "REF", "RefSeq /product")],
    "UL39": [("ICP6", P, "UP:P08543", "UniProt gene synonym 'ICP6'"),
             ("ICP6", P, "TPA", "CDS /note 'icp6'"),
             ("RR1", F, "TPA", "CDS /note 'rr-1'"),
             ("R1", F, "UP:P08543", "UniProt short name 'R1'"),
             ("ribonucleotide reductase subunit 1", F, "REF", "RefSeq /product"),
             ("ribonucleotide reductase large subunit", F, "REF",
              "conventional name in the literature for the UL39 product. NOT the "
              "RefSeq /product string, which is 'ribonucleotide reductase subunit "
              "1'; this row previously claimed it was")],
    "UL40": [("RR2", F, "TPA", "CDS /note 'rr-2'"),
             ("ribonucleotide reductase subunit 2", F, "REF", "RefSeq /product"),
             ("ribonucleotide reductase small subunit", F, "REF",
              "conventional name in the literature for the UL40 product. NOT the "
              "RefSeq /product string, which is 'ribonucleotide reductase subunit "
              "2'; this row previously claimed it was")],
    "UL41": [("vhs", F, "UP:P10225", "UniProt alternative short name 'Vhs'"),
             ("vhs", F, "TPA", "CDS /note 'vhs'"),
             ("virion host shutoff protein", F, "UP:P10225", "UniProt recommended name"),
             ("tegument host shutoff protein", F, "REF", "RefSeq /product")],
    "UL42": [("PAP", P, "UP:P10226", "UniProt alternative short name"),
             ("DNA polymerase processivity subunit", F, "REF", "RefSeq /product"),
             ("DNA polymerase processivity factor", F, "REF",
              "conventional name in the literature for the UL42 product. NOT the "
              "RefSeq /product string, which is 'DNA polymerase processivity "
              "subunit'; this row previously claimed it was")],
    "UL44": [("gC", P, "REF", "RefSeq /product 'envelope glycoprotein C'"),
             ("gC", P, "TPA", "CDS /note 'glycoprotein c'")],
    "UL45": [("18 kDa protein", M, "UP:P10229", "UniProt alternative name")],
    "UL46": [("VP11/12", P, "REF", "RefSeq /product 'tegument protein VP11/12'"),
             ("VP11/12", P, "UP:P10230", "UniProt alternative name"),
             ("VP11", P, "TPA", "CDS /note 'vp11'"),
             ("VP12", P, "TPA", "CDS /note 'vp12'")],
    "UL47": [("VP13/14", P, "REF", "RefSeq /product 'tegument protein VP13/14'"),
             ("VP13/14", P, "UP:P10231", "UniProt alternative name"),
             ("VP13", P, "TPA", "CDS /note 'vp13'"),
             ("VP14", P, "TPA", "CDS /note 'vp14'"),
             ("Vmw82/81", M, "UP:P10231", "UniProt alternative name 'VMW82/81'")],
    "UL49": [("VP22", P, "REF", "RefSeq /product 'tegument protein VP22'"),
             ("VP22", P, "TPA", "CDS /note 'vp22'")],
    "UL50": [("dUTPase", F, "UP:P10234", "UniProt short name"),
             ("deoxyuridine triphosphatase", F, "REF", "RefSeq /product")],
    "UL52": [("DNA primase", F, "UP:P10236", "UniProt recommended name"),
             ("helicase-primase primase subunit", F, "REF", "RefSeq /product")],
    "UL53": [("gK", P, "REF", "RefSeq /product 'envelope glycoprotein K'"),
             ("gK", P, "TPA", "CDS /note 'glycoprotein k'"),
             ("syncytial protein", F, "UP:P68331", "UniProt alternative name")],
}

# The genes whose only non-systematic name is the reference record's descriptive
# /product string. Included so that every gene in the classical inventory carries at
# least one cited non-systematic name, rather than the set quietly stopping at the
# genes that happen to be famous.
for _canon, _prod in {
    "UL3": "nuclear protein UL3",
    "UL4": "nuclear protein UL4",
    "UL20": "envelope protein UL20",
    "UL21": "tegument protein UL21",
    "UL43": "envelope protein UL43",
    "UL45": "membrane protein UL45",
    "UL55": "nuclear protein UL55",
    "UL56": "membrane protein UL56",
    "US2": "virion protein US2",
    "US8A": "membrane protein US8A",
    "US9": "membrane protein US9",
    "US11": "tegument protein US11",
}.items():
    add(_canon, _prod, F, "REF", "", EST,
        "the reference record's /product qualifier, which for this gene is the only "
        "non-systematic name on offer: no ICP, VP, Vmw or glycoprotein name has been "
        "assigned to it")

for canon, items in UL_SIMPLE.items():
    for alias, atype, cit, note in items:
        inv = "both" if cit == "TPA" else "classical"
        add(canon, alias, atype, cit, inv, EST, note)

# UL26 / UL26.5: the protease and its substrate, the hardest VP mapping
add("UL26", "VP24", P, "REF", "classical", EST, "RefSeq /note 'VP24' (the protease, N-terminal region)")
add("UL26", "VP24", P, "TPA", "whisnant2020", EST, "CDS /note 'vp24'")
add("UL26", "VP21", P, "REF", "classical", EST,
    "RefSeq /note 'VP21' (the minor scaffold protein, the remainder of the protein)")
add("UL26", "P40", M, "UP:P10210", "classical", EST, "UniProt alternative name 'Capsid protein P40'")
add("UL26", "P40", M, "TPA", "whisnant2020", EST, "CDS /note 'p40'")
add("UL26", "assemblin", P, "UP:P10210", "classical", EST,
    "UniProt names the cleaved protease chain 'Assemblin'")
add("UL26", "pPR", P, "UP:P10210", "classical", EST, "UniProt alternative short name for the protease precursor")
add("UL26", "capsid maturation protease", F, "REF", "classical", EST, "RefSeq /product")
add("UL26", "capsid scaffolding protein", F, "UP:P10210", "classical", AMB,
    "UniProt recommended name for UL26. One letter from the RefSeq /product for UL26.5, "
    "'capsid scaffold protein'. Two different genes, two near-identical strings: a string "
    "match on either will silently pick the wrong gene. See alias_notes.md",
    # UniProt's recommended name for UL26; the ambiguity is the one letter between it and UL26.5's product string, not the attribution
    RESOLVE)
add("UL26", "AP", P, "UP:P10210", "classical", AMB,
    "UniProt names the cleaved 'Assembly protein (AP)' chain of UL26. That chain covers "
    "the same polypeptide region as the separate UL26.5 product (VP22a, ICP35), so AP is "
    "attributable to either gene depending on whether the full-length precursor or the "
    "nested gene is meant. Not resolved here", DECLINE)
add("UL26.5", "UL26.5", S, "REF", "both", EST,
    "a decimal name that predates the 2020 re-annotation and is in the classical inventory")
add("UL26.5", "VP22a", P, "REF", "classical", EST,
    "RefSeq /note 'VP22a' on UL26.5. Do not merge with VP22 (UL49): different genes, "
    "different proteins")
add("UL26.5", "ICP35", P, "PRESTON1992", "classical", EST,
    "title and abstract: 'the herpes simplex virus (HSV) type 1 assembly protein ICP35', "
    "with the abstract placing it against UL26 ('the UL26 gene specified a protein "
    "distinct from ICP35')")
add("UL26.5", "ICP35", P, "LIUROIZMAN1991", "classical", EST,
    "abstract: 'the UL26.5 gene, which specifies a capsid protein designated infected "
    "cell protein 35'")
add("UL26.5", "capsid scaffold protein", F, "REF", "classical", AMB,
    "RefSeq /product for UL26.5. One letter from UniProt's recommended name for UL26, "
    "'capsid scaffolding protein'. See alias_notes.md",
    # the RefSeq product for UL26.5; same near-miss, same reasoning
    RESOLVE)
add("UL26.5", "infected cell protein 35", P, "LIUROIZMAN1991", "classical", EST,
    "abstract, written out in full: 'the UL26.5 gene, which specifies a capsid protein "
    "designated infected cell protein 35'")

# UL48, the worked example of the whole problem
add("UL48", "UL48", S, "REF", "both", EST, "")
add("UL48", "VP16", P, "REF", "both", EST, "RefSeq /product 'transactivating tegument protein VP16'")
add("UL48", "VP16", P, "TPA", "whisnant2020", EST, "CDS /note 'vp16'")
add("UL48", "VP16", P, "UP:P06492", "classical", EST, "UniProt recommended name 'Tegument protein VP16'")
add("UL48", "alpha-TIF", F, "UP:P06492", "classical", EST, "UniProt alternative short name 'Alpha-TIF'")
add("UL48", "alpha-TIF", F, "TPA", "whisnant2020", EST, "CDS /note 'alpha-tif'")
add("UL48", "alpha trans-inducing protein", F, "UP:P06492", "classical", EST,
    "UniProt alternative name, written out")
add("UL48", "Vmw65", M, "UP:P06492", "classical", EST, "UniProt alternative name 'Vmw65'")
add("UL48", "Vmw65", M, "COUSENS1989", "classical", EST,
    "title: 'the herpes simplex virus regulatory protein, Vmw65'")
add("UL48", "Vmw65", M, "ACE1989", "classical", EST,
    "abstract: 'the gene encoding the transinducing factor Vmw65'")
add("UL48", "ICP25", P, "UP:P06492", "classical", PROB,
    "UniProt alternative name 'ICP25'. No primary paper for this string was located, so "
    "the only citation is the UniProt record")

# UL49A / UL49.5, where the two inventories disagree on the name
add("UL49A", "UL49A", S, "REF", "classical", EST, "the classical record's name for this gene")
add("UL49A", "UL49.5", S, "TPA", "whisnant2020", EST,
    "the 2020 re-annotation names it UL49.5; its CDS /note reads 'ul49a; glycoprotein n', "
    "which is the record stating the equivalence itself")
add("UL49A", "UL49.5", S, "UP:O09800", "both", EST, "UniProt lists gene names 'gN UL49.5 UL49A'")
add("UL49A", "gN", P, "REF", "both", EST, "RefSeq /product 'envelope glycoprotein N'")
add("UL49A", "gN", P, "TPA", "whisnant2020", EST, "CDS /note 'glycoprotein n'")

# UL54, and the IE63 collision
add("UL54", "UL54", S, "REF", "both", EST, "")
add("UL54", "ICP27", P, "REF", "both", EST, "RefSeq /note 'ICP27'")
add("UL54", "ICP27", P, "TPA", "whisnant2020", EST, "CDS /note 'icp27'")
add("UL54", "ICP27", P, "UP:P10238", "classical", EST, "UniProt alternative name 'Infected cell protein 27'")
add("UL54", "IE63", M, "UP:P10238", "classical", EST, "UniProt 'Immediate-early protein IE63'")
add("UL54", "IE63", M, "TPA", "whisnant2020", AMB,
    "CDS /note 'icp27; ie63'. The same record also puts 'ie63' on US1/ICP22, so the "
    "string IE63 is not unique within one authoritative record. See alias_notes.md",
    # two canonicals exist for IE63, so resolve() raises AmbiguousName by the multi-target path; declining here would hide which genes collide
    RESOLVE)
add("UL54", "Vmw63", M, "UP:P10238", "classical", EST, "UniProt alternative name 'VMW63'")
add("UL54", "mRNA export factor", F, "UP:P10238", "classical", EST, "UniProt recommended name")

# --------------------------------------------------------------------------- #
# Literature attestations. Added after the database layer because several of them
# upgrade a row whose only source was a curated record, and because the spelling
# variants below are exactly what a dictionary has to absorb from real text.
# Access level for each is recorded in alias_testcases.md.
# --------------------------------------------------------------------------- #

# UL48, the worked example. Four conventions in one parenthesis, three times over.
add("UL48", "ICP25", P, "ZHANG1991", "", EST,
    "abstract: 'the alpha trans-induction factor (alpha TIF, ICP25, VP16, Vmw65)'. This "
    "is the primary-literature citation that moves ICP25 off the UniProt record alone")
add("UL48", "ICP25", P, "HAGMANN1997", "", EST,
    "abstract: 'the viral transactivator VP16 (Vmw65, alphaTIF, ICP25)'")
add("UL48", "ICP25", P, "MINAGAWA1994", "", EST, "abstract: 'VP16 (ICP25; Vmw65)'")
add("UL48", "VP16", P, "WEINHEIMER1992", "", EST,
    "abstract: 'VP16 (also called Vmw65 and alpha TIF)'. The gloss is prose, not a comma "
    "list, so a parenthesis-only extractor misses it")
add("UL48", "Vmw65", M, "LIU1999", "", EST, "abstract: 'virion protein VP16 (Vmw65; alphaTIF)'")
add("UL48", "alpha TIF", F, "ZHANG1991", "", EST,
    "spaced spelling. PMC renders the Greek letter as the word 'alpha'; the printed paper "
    "will show the Greek character")
add("UL48", "alphaTIF", F, "HAGMANN1997", "", EST, "closed-up spelling")
add("UL48", "alpha trans-induction factor", F, "ZHANG1991", "", EST, "written out, 'induction'")
add("UL48", "alpha trans-inducing factor", F, "POON1995", "", EST,
    "written out, 'inducing'. Title and abstract: 'alpha trans-inducing factor (alpha "
    "TIF, VP16, Vmw65)'")
add("UL48", "alpha-transinducing factor", F, "SPECTOR1991", "", EST,
    "written out, closed up and hyphenated. Title: 'Role of alpha-transinducing factor "
    "(VP16)'. Three spellings of one functional name across three papers")

# RL2 / ICP0
add("RL2", "alpha 0", F, "ROBERTS1988", "", EST,
    "abstract: 'the herpes simplex virus type 1 IE110 (alpha 0)'. Roizman's alpha-gene "
    "kinetic numbering, a fifth convention distinct from both ICP and IE numbering")
add("RL2", "Vmw110", M, "MINAGAWA1994", "", EST, "abstract: 'infected cell protein 0 (ICP0) (Vmw110)'")
add("RL2", "Vmw110", M, "EVERETT1995", "", EST, "abstract: 'Vmw110 (also known as ICP0)'")
add("RL2", "IE110", M, "MULLEN1995", "", EST,
    "title and abstract: 'IE175 (ICP4), IE110 (ICP0), and IE63 (ICP27)'")
add("RL2", "infected cell protein 0", P, "MINAGAWA1994", "", EST, "written out")
add("RL2", "infected-cell polypeptide 0", P, "DUAN2019", "", EST,
    "written out differently: 'infected-cell polypeptide 0 (ICP0)'. Polypeptide, not "
    "protein, and hyphenated")
add("RL2", "ICP0", P, "W20", "", EST,
    "the 2020 paper glosses the classical name against the systematic one: 'The RL2 locus "
    "encodes the major viral immediate early protein ICP0'")

# RS1 / ICP4
add("RS1", "alpha 4", F, "ROBERTS1988", "", EST, "title: 'the ... IE175 (alpha 4) promoter'")
add("RS1", "alpha 4", F, "SPECTOR1991", "", EST, "abstract: 'alpha 4 gene'")
add("RS1", "Vmw175", M, "MINAGAWA1994", "", EST, "abstract: 'ICP4 (Vmw175)'")
add("RS1", "Vmw175", M, "ALLEN1997", "", EST, "abstract: 'immediate early protein Vmw175 (ICP4)'")
add("RS1", "IE175", M, "MULLEN1995", "", EST, "abstract: 'IE175 (ICP4)'")
add("RS1", "IE175", M, "ROBERTS1988", "", EST, "title, glossed twice, to alpha 4 and to ICP4")

# RL1 / ICP34.5
add("RL1", "gamma1 34.5", F, "DUAN2019", "", EST,
    "the paper writes 'STING by ICP27, gamma1 34.5 and UL46' having used ICP34.5 "
    "elsewhere: the same gene named two ways in neighbouring sentences. PMC renders the "
    "Greek letter as the word")
add("RL1", "ICP34.5", P, "W20", "", EST,
    "the 2020 paper glosses it against the systematic name: 'Within the ICP34.5 (RL1) "
    "locus'")

# UL54 / ICP27
add("UL54", "alpha 27", F, "SPECTOR1991", "", EST, "abstract: 'alpha 27 promoter'")
add("UL54", "IE63", M, "MULLEN1995", "", EST,
    "title and abstract: 'IE63 (ICP27)'. An HSV-1 paper assigning IE63 to ICP27, which is "
    "the weight of evidence against the BK012101.1 note that also puts ie63 on US1")
add("UL54", "Vmw63", M, "SINGH1996", "", PROB,
    "abstract: 'homology to ICP27 (Vmw63) of herpes simplex virus 1'. The gloss is "
    "explicit but it sits in a bovine herpesvirus 1 paper describing a homologue, not an "
    "HSV-1 primary paper. UniProt P10238 carries the same equivalence")
add("UL54", "ICP27", P, "W20", "", EST,
    "the 2020 paper glosses it: 'the major immediate early protein ICP27 (UL54)'")

# The capsid VP series, where the gel bands do not map cleanly onto genes
add("UL19", "ICP5", P, "BRAUN1984", "", EST,
    "abstract: 'Surface iodination labeled the major capsid protein (ICP5) of empty "
    "capsids' and 'ICP5, the major capsid protein'. The only attestation of this string "
    "located in an accessible source; it is in neither UniProt nor the two feature tables")
add("UL19", "VP5", P, "VILLANUEVA2021", "", EST,
    "'The HSV-1 capsid is comprised mainly of the major capsid protein VP5 (UL19)'")
add("UL19", "VP5", P, "JIN2014", "", EST, "glossed in both directions: 'VP5 (UL19)' and 'UL19(VP5)'")
add("UL19", "VP5", P, "W20", "", EST, "'the major capsid protein (VP5, UL19)'")
add("UL18", "VP23", P, "VILLANUEVA2021", "", EST,
    "'composed of two VP23 (UL18) proteins and one VP19C (UL38) protein'")
add("UL18", "VP23", P, "JIN2014", "", EST, "glossed in both directions: 'VP23 (UL18)' and 'UL18(VP23)'")
add("UL38", "VP19C", P, "VILLANUEVA2021", "", EST, "same sentence as UL18/VP23")
add("UL38", "VP19C", P, "JIN2014", "", EST, "'a null mutant of UL38 (VP19C)'")
add("UL35", "VP26", P, "JIN2014", "", EST, "table row 'UL35(VP26)'")
add("UL35", "VP26", P, "W20", "", EST,
    "'a fusion protein of mCherry on the N-terminus of the UL35 gene product (VP26)'")
add("UL26", "VP24", P, "JIN2014", "", EST, "'the two products encoded by UL26, VP24 and VP21'")
add("UL26", "VP21", P, "JIN2014", "", EST, "same sentence")
add("UL26.5", "ICP35", P, "JIN2014", "", EST, "table row 'UL26.5(ICP35)'")
add("UL26.5", "ICP35", P, "BRAUN1984", "", EST,
    "abstract: 'a genetically and immunologically related family of structural (gamma) "
    "polypeptides ... designated as infected cell polypeptides (ICP) 35', resolved into "
    "bands ICP35a to ICP35f")
add("UL26.5", "virion protein 22a", P, "BRAUN1984", "", EST,
    "abstract: 'One of these corresponds to virion protein 22a, a component of capsids "
    "containing DNA'. This is the VP22a name written out, and the sentence is what ties "
    "the ICP35 family to the VP22a band")
add("UL6", "pUL6", O, "JIN2014", "", EST,
    "'a dodecamer of pUL6'. The pUL prefix is a general protein-form convention applied "
    "to any UL gene, not a name specific to UL6; it is recorded here because a dictionary "
    "that does not strip the p will fail on it")

# The VP11/12 and VP13/14 slash forms
add("UL46", "VP11", P, "ZHANGMCKNIGHT1993", "", EST,
    "title and abstract: 'UL46 encoded the viron tegument phosphoproteins VP11 and VP12 "
    "(VP11/12)'")
add("UL46", "VP12", P, "ZHANGMCKNIGHT1993", "", EST, "same sentence")
add("UL46", "VP11/12", P, "ZHANGMCKNIGHT1993", "", EST,
    "the slash form is a distinct string shape: a dictionary must either hold it whole or "
    "decompose it")
add("UL47", "VP13", P, "ZHANGMCKNIGHT1993", "", EST,
    "'UL47 also encoded two virion tegument phosphoproteins, VP13 and VP14 (VP13/14)'")
add("UL47", "VP14", P, "ZHANGMCKNIGHT1993", "", EST, "same sentence")
add("UL47", "VP13/14", P, "ZHANGMCKNIGHT1993", "", EST, "slash form")
add("UL47", "VP13/14", P, "MCKNIGHT1994", "", EST, "title uses the slash form")

# Glycoproteins, glossed to UL numbers
add("UL1", "gL", P, "JAMBUNATHAN2021", "", EST,
    "'Glycoprotein gL, a 224aa glycoprotein encoded by the UL1 gene' and 'UL1 (gL)'")
add("UL22", "gH", P, "JAMBUNATHAN2021", "", EST, "'UL22 (gH)'")
add("UL53", "gK", P, "JAMBUNATHAN2021", "", EST, "'UL53 (gK)'")

# Antiviral-resistance usage: TK and UL23 and UL30 in one paper
add("UL23", "TK", F, "IBRAHIM2026", "", EST,
    "glossed in both directions: 'those in UL23 (TK)' and 'the viral thymidine kinase "
    "(TK), which is encoded by the UL23 gene'")
add("UL23", "thymidine kinase", F, "IBRAHIM2026", "", EST, "'UL23 (thymidine kinase)'")
add("UL23", "thymidine kinase", F, "MINAGAWA1994", "", EST, "named alongside ICP0, ICP4, ICP27 and VP16")
add("UL23", "TK", F, "DUAN2019", "", EST, "'HSV-1 thymidine kinase (TK)'")
add("UL30", "DNA polymerase", F, "IBRAHIM2026", "", EST, "'UL30 (DNA polymerase)'")
add("UL30", "DNA pol", F, "IBRAHIM2026", "", EST, "'HSV DNA pol (UL30)', the abbreviated form")

# UL41 / VHS
add("UL41", "VHS", F, "DUAN2019", "", EST,
    "'HSV-1 virion host shutoff (VHS) protein encoded by UL41 gene'. The paper uses both "
    "VHS and vhs, so case folding matters")
add("UL41", "virion host shutoff", F, "DUAN2019", "", EST, "written out in the same sentence")

# LAT variants
add("LAT", "LATs", O, "MINAGAWA1994", "", EST,
    "abstract uses the plural 'Latency-associated transcripts' and 'LATs'. A dictionary "
    "keyed on the singular will miss it")
add("LAT", "latency-associated transcript", F, "DUAN2019", "", EST,
    "'regulation of the latency-associated transcript (LAT)'")

# UL49.5 in the paper, not only in the deposition
add("UL49A", "UL49.5", S, "W20", "", EST,
    "named in the main text: 'Six of the previously identified CDS (UL11, UL49.5, US5, "
    "US9, US12, and RL2 iso1) are under 100 aa'")

# --------------------------------------------------------------------------- #
# Name-origin papers, where the paper names the gene rather than the convention.
# Papers that founded a whole convention without naming an individual gene are
# recorded in alias_notes.md, not here, because they attest no specific string.
# --------------------------------------------------------------------------- #

# UL48: the naming papers for alpha-TIF
add("UL48", "Vmw65", M, "CAMPBELL1984", "", EST,
    "abstract: 'Vmw65 is the virion component which trans-activates HSV IE "
    "transcription'. The paper that identified the protein behind the function")
add("UL48", "alpha-TIF", F, "MCKNIGHT1987", "", EST,
    "the paper that defines the abbreviation. Abstract: 'This protein, designated as the "
    "alpha-trans-inducing factor (alpha-TIF), maps within the BamHI F fragment'")
add("UL48", "alpha-trans-inducing factor", F, "MCKNIGHT1987", "", EST,
    "fully hyphenated spelling, from the same sentence and the paper's title")
add("UL48", "alpha trans-inducing factor", F, "PELLETT1985", "", EST,
    "earliest located use of the spelled-out term: 'the alpha trans-inducing factor was a "
    "structural component of the virion'")
add("UL48", "Vmw65", M, "OHARE1985", "", EST,
    "abstract uses 'IE175K protein', 'IE110K protein' and 'Vmw65' together")

# LAT: the discovery paper and the naming paper are different papers
add("LAT", "LAT", O, "WAGNER1988", "", EST,
    "the paper that coined the abbreviation. Abstract: 'Most of this latency-associated "
    "transcript (LAT) was isolated in the poly(A)- fraction from ganglia'. Earliest "
    "MEDLINE-indexed use of the phrase")
add("LAT", "latency-associated transcript", F, "WAGNER1988", "", EST, "same sentence")
add("LAT", "ICP-0-related transcripts", O, "STEVENS1987", "", PROB,
    "the discovery paper, which does NOT use LAT or latency-associated transcript: its "
    "abstract calls them 'ICP-0-related transcripts' and 'this RNA species'. Recorded so "
    "the dictionary does not credit the name to the wrong paper. Weak as an alias")

# UL41: note the original expansion differs from the modern one
add("UL41", "vhs", F, "READ1983", "", EST,
    "the naming paper. Abstract: 'They were consequently designated as "
    "virion-associated host shutoff (vhs) mutants'")
add("UL41", "virion-associated host shutoff", F, "READ1983", "", EST,
    "the ORIGINAL expansion. The now-standard 'virion host shutoff' drops "
    "'associated', so the two expansions are not the same string")
add("UL41", "vhs", F, "READ1993", "", EST,
    "the paper that ties the name to the gene; title: 'the vhs (UL41) polypeptide'")

# Glycoproteins: the letter-assignment papers that could be verified
add("UL53", "gK", P, "HUTCHINSONGK1992", "", EST,
    "the naming paper. Abstract: 'A single 40-kDa protein, designated gK (the ninth "
    "HSV-1 glycoprotein to be described)'")
add("UL1", "gL", P, "HUTCHINSONGL1992", "", EST,
    "the naming paper. Abstract: 'This novel glycoprotein is the 10th HSV-1 glycoprotein "
    "to be described and was named glycoprotein L (gL)'")
add("UL1", "glycoprotein L", P, "HUTCHINSONGL1992", "", EST, "written out, same sentence")
add("US4", "gG", P, "FRAME1986", "", EST,
    "abstract: 'the product of gene US4 of HSV-1 is the equivalent of gG'")
add("US8", "gE1", P, "BAUCKE1979", "", PROB,
    "abstract: 'one polypeptide, designated pE, is processed to yield gE1, which is in "
    "turn processed to yield gE2'. Note the string is gE1, not gE and not gE-1: three "
    "spellings, and this paper attests only the first")

# UL49.5 is a CLASSICAL name, not a 2020 one. This row is what makes the
# inventory column come out as "both" rather than "whisnant2020".
add("UL49A", "UL49.5", S, "BARKER1992", "", EST,
    "the gene-designation paper, 28 years before the re-annotation. Title: 'an "
    "additional translated open reading frame designated UL49.5'. So UL49.5 belongs to "
    "the classical inventory as well, and treating it as a 2020 name would be wrong")

# UL39 / ribonucleotide reductase, including the HSV-2 naming difference
add("UL39", "Vmw136", M, "PRESTONVG1984", "", EST,
    "abstract: 'a polypeptide of 140 000 mol. wt. (Vmw136, ICP6)'. A clean Vmw to ICP "
    "cross-reference, and the only located attestation of this string")
add("UL39", "ICP6", P, "PRESTONVG1984", "", EST, "same parenthesis")
add("UL39", "ICP6", P, "GOLDSTEIN1988", "", EST,
    "title: 'an ICP6 lacZ insertion mutant'; confirms ICP6 as the reductase large "
    "subunit gene")
add("UL39", "RR1", F, "CHUNG1989", "", EST, "abstract: 'ribonucleotide reductase (RR1)'")
add("UL39", "ICP10", P, "CHUNG1989", "", AMB,
    "ICP10 is the HSV-2 name for the positional homologue, not an HSV-1 name. The "
    "abstract uses both in one sentence: 'The large subunit of the herpes simplex virus "
    "type 2 (HSV-2) ribonucleotide reductase (RR1) ... The HSV-1 RR1 was shown to be "
    "negative for PK activity'. The HSV-1 and HSV-2 ICP series were numbered "
    "independently by gel mobility, so ICP6 and ICP10 are homologues with different "
    "numbers. Held as a cross-species ambiguity, not merged", DECLINE)

add("UL28", "ICP18.5", P, "KNOPF1987", "", EST,
    "title: 'the herpes simplex virus type 1 late gene ICP18.5 of strain Angelotti'")

# UL26.5 / the family 35 proteins
add("UL26.5", "family 35 proteins", O, "LIUROIZMAN1991B", "", EST,
    "the paper that ties the name to the gene. Abstract: 'The second transcriptional "
    "unit, designated UL26.5 ... encodes the family 35 proteins'")
add("UL26.5", "protein 22a", P, "GIBSON1974", "", EST,
    "abstract: 'protein 22a of B capsids, which does not have an electrophoretically "
    "identical counterpart in virions, shares a relatively unique set of staining and "
    "radiolabeling properties with virion protein 22, which has a slightly more rapid "
    "electrophoretic mobility'. This is the sentence that establishes 22a and 22 as "
    "different species, in 1974, before either was assigned to a gene")
add("UL49", "VP22", P, "ELLIOTT1992", "", EST,
    "the gene assignment. Abstract: 'The results clearly indicate that the protein "
    "encoded by UL49 is VP22'. With GIBSON1974 on UL26.5 this is what keeps VP22 and "
    "VP22a apart: different genes, not a spelling variant")

# UL54 / IE63, and the cross-species origin of the collision
add("UL54", "IE63", M, "MCGREGOR1996", "", EST,
    "abstract opens: 'The essential herpes simplex virus type 1 (HSV-1) immediate-early "
    "IE63 (ICP27) is pleiotropic in function'. An unambiguous HSV-1 assignment")
add("UL54", "IE63", M, "PERRY1988", "", EST,
    "abstract places the gene in UL, not US: 'the region of UL adjacent to IRL also "
    "contains three genes, one encoding the immediate early protein IE63'")

# RS1 and RL2: the IE weight numbering, and its K-suffixed variant forms
add("RS1", "IE175", M, "MCGEOCH1986SR", "", EST,
    "abstract: 'This contains immediate early gene 3, encoding the IE175 protein, an "
    "important transcriptional activator of later virus genes'")
add("RL2", "IE gene 1", O, "EVERETT1984", "", EST,
    "abstract: 'the product of IE gene 1, Vmw110 (ICP 0)'. The IE GENE numbering (1 to 5) "
    "is a separate scheme from the IE molecular-weight numbering (IE110, IE175, IE63): "
    "RL2 is IE gene 1 and also IE110, and the two numbers have nothing to do with each "
    "other")
add("RS1", "IE gene 3", O, "MCGEOCH1986SR", "", EST,
    "abstract: 'This contains immediate early gene 3, encoding the IE175 protein'. So RS1 "
    "is IE gene 3 and IE175 at once")
add("US12", "alpha47", F, "FIELDS6CH60", "", PROB,
    "the chapter's prose gene list writes the entry as 'alpha47/US12 (88 codons, alpha, "
    "nonessential)'. Access: a scanned copy read in part, so reported rather than quoted "
    "from the printed page")
add("RS1", "IE175K", M, "OHARE1985", "", EST,
    "abstract writes 'IE175K protein', with a K suffix. A dictionary keyed on IE175 "
    "alone will miss it")
add("RL2", "IE110K", M, "OHARE1985", "", EST, "abstract writes 'IE110K protein'")
add("RL2", "Vmw IE 110", M, "PRESTON1978", "", EST,
    "earliest located use of the Vmw convention, and it is spaced: 'the immediate early "
    "polypeptide Vmw IE 110 of HSV-1 and its functionally equivalent polypeptide, Vmw IE "
    "118, of HSV-2'. Note Vmw IE 118 is the HSV-2 counterpart, not an HSV-1 name")

# ORF O and ORF P: the Roizman-lab naming papers
add("ORF-P", "ORF P", O, "LAGUNOFF1994", "", EST,
    "the naming paper, and it is spaced not hyphenated. Abstract: 'One ORF (ORF P), "
    "coincident with but antisense to the gamma(1)34.5 gene, was expressed but only "
    "under conditions in which ICP4 was not functional'")
add("ORF-O", "ORF O", O, "RANDALL1997", "", EST,
    "the naming paper, spaced. Title: 'The product of ORF O located within the domain of "
    "herpes simplex virus 1 genome transcribed during latent infection'")

# --- unique short region -------------------------------------------------- #
add("US1", "US1", S, "REF", "both", EST, "")
add("US1", "IE63", M, "JACKERS1992", "", AMB,
    "NOT an HSV-1 attestation. This varicella-zoster virus paper is the likely source of "
    "the collision: 'Varicella-zoster virus (VZV) gene 63 encodes a protein (IE63) ... "
    "which has amino acid similarities to the immediate-early (IE) protein 22 (ICP22) of "
    "herpes simplex virus type 1'. VZV IE63 is the homologue of HSV-1 ICP22, which is "
    "how ie63 plausibly reached US1 in BK012101.1. The HSV-1 reading of IE63 is UL54. "
    "Recorded as ambiguity with its explanation, not resolved away",
    # as above, the IE63 collision is reported by the multi-target path
    RESOLVE)
add("US1", "ICP22", P, "REF", "both", EST, "RefSeq /product 'regulatory protein ICP22'")
add("US1", "ICP22", P, "TPA", "whisnant2020", EST, "CDS /note 'icp22'")
add("US1", "ICP22", P, "UP:P04485", "classical", EST, "UniProt 'Transcriptional regulator ICP22'")
add("US1", "IE68", M, "UP:P04485", "classical", EST, "UniProt 'Immediate-early protein IE68'")
add("US1", "IE63", M, "TPA", "whisnant2020", AMB,
    "CDS /note 'icp22; ie63'. This conflicts with UniProt P10238, which assigns IE63 to "
    "UL54/ICP27, and with UniProt P04485, which gives IE68 for ICP22. Recorded as an "
    "ambiguity, not resolved. See alias_notes.md",
    # as above
    RESOLVE)
add("US1", "Vmw68", M, "SCHWYZER1994", "classical", EST,
    "abstract gloss: 'designated BICP22 because of its homology to ICP22 (Vmw68) of "
    "herpes simplex virus type 1'. The gloss is in a bovine herpesvirus 1 paper making a "
    "comparative statement, not an HSV-1 primary paper, but it is explicit and it is the "
    "only attestation of the string located. UniProt P04485 gives IE68, not Vmw68")
add("US1", "alpha22", F, "JOVASEVIC2010", "classical", EST,
    "title: 'a domain encoding US5, US4, US3, US2 and alpha22'")
add("US1", "alpha 22", F, "CARTER1996", "classical", EST,
    "title: 'the open reading frame of the alpha 22 gene'")
add("US1.5", "US1.5", S, "TPA", "whisnant2020", EST,
    "BK012101.1 CDS /note 'truncated isoform of us1'")
add("US1.5", "US1.5", S, "CARTER1996", "classical", EST,
    "the primary paper for the shorter in-frame alpha gene inside alpha22; it predates "
    "the 2020 re-annotation, so this name belongs to both inventories")

US_SIMPLE = {
    "US2": [],
    "US3": [("protein kinase US3", F, "REF", "RefSeq /product 'serine/threonine protein kinase US3'")],
    "US4": [("gG", P, "REF", "RefSeq /product 'envelope glycoprotein G'"),
            ("gG", P, "TPA", "CDS /note 'glycoprotein g'"),
            ("gG-1", P, "UP:P06484", "UniProt alternative name; the -1 marks the HSV-1 form")],
    "US5": [("gJ", P, "REF", "RefSeq /product 'envelope glycoprotein J'"),
            ("gJ", P, "TPA", "CDS /note 'glycoprotein j'")],
    "US6": [("gD", P, "REF", "RefSeq /product 'envelope glycoprotein D'"),
            ("gD", P, "TPA", "CDS /note 'glycoprotein d'")],
    "US7": [("gI", P, "REF", "RefSeq /product 'envelope glycoprotein I'"),
            ("gI", P, "TPA", "CDS /note 'glycoprotein i'")],
    "US8": [("gE", P, "REF", "RefSeq /product 'envelope glycoprotein E'"),
            ("gE", P, "TPA", "CDS /note 'glycoprotein e'"),
            ("gE-1", P, "UP:P04488", "UniProt alternative name; the -1 marks the HSV-1 form")],
    "US9": [("10 kDa protein", M, "UP:P06481", "UniProt alternative name")],
    "US10": [("virion protein US10", P, "REF", "RefSeq /product")],
    "US11": [("Vmw21", M, "UP:P04487", "UniProt alternative name 'Vmw21'"),
             ("Vmw21", M, "TPA", "CDS /note 'vmw21'")],
}
for canon, items in US_SIMPLE.items():
    for alias, atype, cit, note in items:
        inv = "both" if cit == "TPA" else "classical"
        add(canon, alias, atype, cit, inv, EST, note)

add("US8A", "US8A", S, "REF", "classical", EST, "the classical record's name")
add("US8A", "US8.5", S, "UP:O09802", "both", EST,
    "UniProt gene name 'US8.5' and protein name 'Protein US8.5' for the same gene")
add("US8A", "US8.5", S, "TPA", "whisnant2020", EST,
    "the 2020 re-annotation uses US8.5; BK012101.1 /note reads 'encodes US8A; ... us8a'")

add("US12", "US12", S, "REF", "both", EST, "")
add("US12", "ICP47", P, "REF", "both", EST, "RefSeq /product 'TAP transporter inhibitor ICP47'")
add("US12", "ICP47", P, "TPA", "whisnant2020", EST, "CDS /note 'icp47'")
add("US12", "IE12", M, "UP:P03170", "classical", EST, "UniProt 'Immediate-early protein IE12'")
add("US12", "IE12", M, "TPA", "whisnant2020", EST, "CDS /note 'ie12'")
add("US12", "Vmw12", M, "UP:P03170", "classical", EST, "UniProt alternative name 'Vmw12'")
add("US12", "IE-5", M, "UP:P03170", "classical", PROB,
    "UniProt alternative name 'Immediate-early-5'. No primary paper located for the "
    "string, so the UniProt record is the only citation")

add("US5.1", "US5.1", S, "GTF", "whisnant2020", EST, "gene_id in the Zenodo GTF")
add("US5.1", "US5-1", O, "JOVASEVIC2010", "classical", EST,
    "the primary paper writes it US5-1; the re-annotation writes it US5.1. BK012101.1 "
    "/note cites 'Jovasevic & Roizman Virol J 2010. not observed in TiSS or PacBio data'")

# --- plain systematic names, cited to the reference record ----------------- #
PLAIN = (
    ["UL%d" % n for n in range(1, 57)]
    + ["US%d" % n for n in range(1, 13)]
    + ["US8A", "UL26.5", "UL49A"]
)
# --------------------------------------------------------------------------- #
# The ICTV standardised protein names, from Table 2 of the ICTV Orthoherpesviridae
# report chapter: the 43 genes detectably conserved across the family, with the
# protein name ICTV applies so that orthologues carry one name in every
# herpesvirus. These are the same strings the reference record's /product
# qualifiers carry, which is not a coincidence: Davison 2010 states the policy and
# annotated the record. Double-citing them here is deliberate, because it is the
# difference between a name a database happens to use and a name a taxonomy body
# has settled.
#
# Note what this table does NOT contain, verified by regex over the whole chapter:
# no ICP number, no VP number, no Vmw string, no alpha-TIF, and no US gene at all.
# It is an authority for the standardised protein layer only.
# --------------------------------------------------------------------------- #

ICTV_CORE = {
    "UL30": "DNA polymerase catalytic subunit",
    "UL42": "DNA polymerase processivity subunit",
    "UL9": "DNA replication origin-binding helicase",
    "UL5": "helicase-primase helicase subunit",
    "UL8": "helicase-primase subunit",
    "UL52": "helicase-primase primase subunit",
    "UL29": "single-stranded DNA-binding protein",
    "UL23": "thymidine kinase",
    "UL39": "ribonucleotide reductase subunit 1",
    "UL40": "ribonucleotide reductase subunit 2",
    "UL50": "deoxyuridine triphosphatase",
    "UL2": "uracil-DNA glycosylase",
    "UL12": "deoxyribonuclease",
    "UL15": "DNA packaging terminase subunit 1",
    "UL28": "DNA packaging terminase subunit 2",
    "UL6": "capsid portal protein",
    "UL25": "DNA packaging tegument protein UL25",
    "UL32": "DNA packaging protein UL32",
    "UL33": "DNA packaging protein UL33",
    "UL17": "DNA packaging tegument protein UL17",
    "UL31": "nuclear egress lamina protein",
    "UL34": "nuclear egress membrane protein",
    "UL19": "major capsid protein",
    "UL18": "capsid triplex subunit 2",
    "UL38": "capsid triplex subunit 1",
    "UL35": "small capsid protein",
    "UL26": "capsid maturation protease",
    "UL26.5": "capsid scaffold protein",
    "UL7": "tegument protein UL7",
    "UL11": "myristylated tegument protein",
    "UL14": "tegument protein UL14",
    "UL16": "tegument protein UL16",
    "UL36": "large tegument protein",
    "UL37": "tegument protein UL37",
    "UL51": "tegument protein UL51",
    "UL27": "envelope glycoprotein B",
    "UL1": "envelope glycoprotein L",
    "UL22": "envelope glycoprotein H",
    "UL10": "envelope glycoprotein M",
    "UL49A": "envelope glycoprotein N",
    "UL13": "tegument serine/threonine protein kinase",
    "UL54": "multifunctional expression regulator",
    "UL24": "nuclear protein UL24",
}
assert len(ICTV_CORE) == 43, len(ICTV_CORE)

for canon, pname in ICTV_CORE.items():
    add(canon, pname, F, "ICTV2021", "", EST,
        "ICTV Table 2 standardised protein name, in the HSV-1 column of the 43 "
        "conserved core genes. The chapter's body states 'A subset of 43 genes is "
        "detectably conserved among family members'; Davison 2010 says 44, so the count "
        "itself is contested between two ICTV-authored sources")
    add(canon, pname, F, "REF", "", EST,
        "the same string as the reference record's /product qualifier. Davison 2010 "
        "states the policy behind it: original gene names retained, standard names "
        "applied to the encoded proteins")

# Davison's own review writes two genes with decimal names where the ICTV chapter
# does not. Both are ICTV-authored and Davison is on both, so this is a live conflict
# between two authorities rather than an error in either.
add("UL49A", "UL49.5", S, "DAVISON2010", "", EST,
    "Davison 2010 writes UL49.5 in its body text while the ICTV report chapter's Table 2 "
    "writes UL49A for the same gene. Two ICTV-authored sources, Davison on both, "
    "disagreeing on the name. Recorded as a dual attestation, not resolved")
add("UL49A", "UL49A", S, "ICTV2021", "", EST,
    "the ICTV Table 2 spelling, against Davison 2010's UL49.5")
add("US8A", "US8.5", S, "DAVISON2010", "", EST, "Davison 2010 body text uses US8.5")

# Fields Virology, the reference chapter that adjudicates the synonymy directly.
add("UL48", "VP16", P, "FIELDS6CH59", "", PROB,
    "the chapter's discussion of HSV naming synonymy sets out VP16 = ICP25 = vmw65 = "
    "alpha-TIF = UL48 as one equivalence. Access: a scanned copy, read in part, with no "
    "verbatim quotation captured, so the gloss is reported rather than quoted. Confirming "
    "it needs the printed page. This is the strongest adjudicating source for the UL48 "
    "set and the weakest-verified citation in the file")
add("RS1", "ICP4", P, "FIELDS6CH59", "", PROB,
    "the same discussion sets out ICP4 = vmw175 = alpha4. Same access caveat as the UL48 "
    "row above")
add("UL48", "UL48", S, "FIELDS6CH60", "", PROB,
    "the classical inventory this dictionary calls 'classical' is the gene list in this "
    "chapter, which the 2020 re-annotation names as its own baseline: its supplementary "
    "rule 1 reads 'All viral ORFs and transcript mentioned in the 6th Edition of Fields "
    "of Virology were included'. Access: table of contents verified, body read in part "
    "from a scan")

# Read the re-annotation's gene_id set from the GTF so neither list can drift from the
# file that defines it.
GTF_PATH = os.path.join(REPO, "data", "raw", "HSV1_Annotation.gtf")
ids = set()
if os.path.exists(GTF_PATH):
    import re as _re

    for line in open(GTF_PATH):
        m = _re.search(r'gene_id "([^"]*)"', line)
        if m:
            ids.add(m.group(1))
else:
    sys.stderr.write("warning: %s absent, run 'make inputs' first\n" % GTF_PATH)

have = {c for (c, a, *_rest) in ROWS if a == c}
for name in PLAIN:
    if name not in have:
        add(name, name, S, "REF", "", EST,
            "systematic positional name, as carried in the reference record's feature table")
    if name in ids:
        add(name, name, S, "GTF", "", EST,
            "the 2020 re-annotation keeps this classical systematic name unchanged: it is "
            "a gene_id in HSV1_Annotation.gtf")
    # The paper that introduced the numbering this name belongs to. Confidence is
    # probable, not established, and the reason is in the note: the numbering is
    # verified from the abstract, the individual string is not.
    if _re.fullmatch(r"UL\d+", name):
        add(name, name, S, "MCGEOCH1988", "", PROB,
            "the paper that introduced the UL numbering. Its abstract states 'Fifty-six "
            "genes were identified' but names no individual gene, and the full text is "
            "paywalled with no PMC copy, so what is verified is the numbering scheme and "
            "its 56-gene extent, not the printed string '" + name + "'")
    elif _re.fullmatch(r"US\d+", name):
        add(name, name, S, "MCGEOCH1985", "", PROB,
            "the paper that introduced the US numbering. Its abstract states 'We consider "
            "that this region contains 12 genes' but names no individual gene, so what is "
            "verified is the numbering scheme and its 12-gene extent, not the printed "
            "string '" + name + "'")

# --- the re-annotation's own novel names ---------------------------------- #
novel = []
if ids:
    classical_names = set(PLAIN) | {
        "RL1", "RL2", "RS1", "LAT", "LAT_gene", "ORF-O", "ORF-P", "ORF-O/P",
        "AL-RNA", "Ori-S", "TRL1", "TRL2", "IRL1", "IRL2", "IRS1", "TRS1",
        "US1.5", "US5.1", "US8.5", "UL49.5",
    }
    novel = sorted(i for i in ids if i not in classical_names)
    for name in novel:
        note = ("open reading frame or transcript resolved by the 2020 re-annotation and "
                "absent from the classical inventory; the decimal suffix places it "
                "between the two flanking classical genes")
        if name.startswith("pri-miRNA"):
            note = "primary miRNA transcript named in the 2020 re-annotation"
        elif name in ("IRL1.5", "TRL1.5", "IRL2.5", "TRL2.5", "IRL3", "TRL3"):
            note = ("repeat-region feature resolved by the 2020 re-annotation; the two "
                    "copies are named separately")
        add(name, name, S, "GTF", "", EST, note)

# --------------------------------------------------------------------------- #

HEADER = [
    "canonical_systematic_name",
    "alias",
    "alias_type",
    "citation",
    "citation_identifier",
    "inventory",
    "confidence",
    "notes",
    "resolution",
    "verbatim_in_citation",
]

# Which deposited record each citation identifier is, for the verbatim check below.
LOCAL_RECORDS = {
    "NC_001806.2": "ref_NC_001806.2.gb",
    "NC_001806.1": "ref_NC_001806.1.gb",
    "BK012101.1": "alt_BK012101.1.gb",
    "JN555585.1": "base_JN555585.1.gb",
}


def verbatim_in(ident: str, alias: str, raw: str) -> str:
    """Whether this alias string actually occurs in the record it cites.

    `yes`, `no`, or `` when the citation is a paper rather than a file here.

    Round 10: 11 aliases cited a deposited record that does not contain them. Eight are
    conventional abbreviations whose notes disclose the record's own wording, which is
    honest but is not what README.md promises when it says every alias carries a
    citation. Three claimed to be the RefSeq /product string and were not: the record
    says "ribonucleotide reductase subunit 1", not "large subunit". Those three are
    corrected; the column makes the distinction explicit for the rest, and the suite
    checks it rather than leaving it to the next reader to notice.
    """
    filename = LOCAL_RECORDS.get(ident)
    if not filename:
        return ""
    path = os.path.join(raw, filename)
    if not os.path.exists(path):
        return ""
    with open(path) as fh:
        text = fh.read().lower()
    # A word-boundary match, not a substring. Round 11: `gL`, `gM`, `gH`, `gN`, `gK`,
    # `gG`, `gD`, `gE`, `gI` and `Pol` were all marked attested because the record
    # contains "glycoprotein ..." and "polymerase", so 15 of the 220 yes rows were wrong
    # and the column contradicted itself within one class of alias: gB said no.
    pattern = r"(?<![0-9a-z])" + re.escape(alias.lower()) + r"(?![0-9a-z])"
    return "yes" if re.search(pattern, text) else "no"

# Which gene inventory each source speaks for. The inventory column is DERIVED from
# this rather than asserted per row, so a gene/alias pair reads "both" only when it is
# genuinely attested on both sides. The two inventories are never silently merged.
#
# The 2020 re-annotation is represented by BK012101.1, the Zenodo GTF and the paper.
# Everything else (the reference record, UniProt curation keyed to the classical gene
# set, and every primary paper here, all of which predate 2020) speaks for the
# classical inventory.
WHISNANT_SOURCES = {"TPA", "GTF", "W20"}


def inventory_for(citkeys):
    w = any(k in WHISNANT_SOURCES for k in citkeys)
    c = any(k not in WHISNANT_SOURCES for k in citkeys)
    if w and c:
        return "both"
    return "whisnant2020" if w else "classical"


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    pair_sources = {}
    for canon, alias, _atype, citkey, _inv, _conf, _notes, _res in ROWS:
        pair_sources.setdefault((canon, alias), set()).add(citkey)

    seen = set()
    out = []
    for canon, alias, atype, citkey, _inv, conf, notes, resolution in ROWS:
        if citkey not in CIT:
            sys.exit("unknown citation key: %s (row %s/%s)" % (citkey, canon, alias))
        inv = inventory_for(pair_sources[(canon, alias)])
        citation, ident = CIT[citkey]
        key = (canon, alias, citkey)
        if key in seen:
            continue
        seen.add(key)
        row = [canon, alias, atype, citation, ident, inv, conf, notes,
               resolution,
               verbatim_in(ident, alias, os.path.join(REPO, "data", "raw"))]
        for cell in row:
            if "\t" in cell or "\n" in cell:
                sys.exit("tab or newline in cell: %r" % cell)
            # Em dash and en dash, written as escapes so that a global dash sweep over
            # this file cannot turn the check into one that rejects ordinary hyphens.
            if "—" in cell or "–" in cell:
                sys.exit("em dash or en dash not allowed: %r" % cell)
        out.append(row)

    out.sort(key=lambda r: (r[0], r[1], r[4]))
    with open(OUT, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for row in out:
            fh.write("\t".join(row) + "\n")

    genes = {r[0] for r in out}
    aliases = {(r[0], r[1]) for r in out}
    print("wrote %s" % OUT)
    print("  %d rows, %d canonical genes, %d distinct gene/alias pairs"
          % (len(out), len(genes), len(aliases)))
    for inv in ("classical", "whisnant2020", "both"):
        print("  inventory=%-12s %d rows" % (inv, sum(1 for r in out if r[5] == inv)))
    for conf in ("established", "probable", "ambiguous"):
        print("  confidence=%-12s %d rows" % (conf, sum(1 for r in out if r[6] == conf)))
    uncited = [r for r in out if not r[3] or not r[4]]
    print("  rows with no citation: %d" % len(uncited))
    print("  novel re-annotation names carried: %d" % len(novel))


if __name__ == "__main__":
    main()
