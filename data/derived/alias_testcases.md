# M1 acceptance test cases: papers that name the same genes differently

Thirteen published papers, each using a different naming convention or a different
spelling for genes the dictionary also holds. These are the M1 gate's "reproduces the
synonymy used in at least three published papers" cases, and they are written as tests
rather than as a reading list: each one says what the dictionary must resolve to pass.

Every name string below was read in the text named under **Access**. Strings that appear
only in a paper's reference list are not counted as that paper's own usage. Where a
string is quoted, it is quoted as the accessed text prints it.

**One typography caveat that matters for exact-string testing.** MEDLINE and PMC render
Greek letters as spelled-out words, so where a case below records `alpha TIF`,
`alpha 4` or `gamma1 34.5`, the printed paper almost certainly prints the Greek
character. The dictionary is tested against the PMC rendering, which is what a text
pipeline will actually be handed. The typography of the printed page was not checked.

Access levels used below:

- **full text read** - the body of the paper was downloaded and searched.
- **abstract only** - PMC or Europe PMC served the abstract and reference list; the body
  is a scanned PDF or paywalled. Only names in the title and abstract are recorded.

---

## Test case 1: all four conventions in one parenthesis

Zhang Y, Sirko DA, McKnight JL (1991) Role of herpes simplex virus type 1 UL46 and UL47
in alpha TIF-mediated transcriptional induction: characterization of three viral deletion
mutants. *J Virol* 65:829-841. PMID 1846201. DOI 10.1128/jvi.65.2.829-841.1991.

**Access: abstract only.**

The abstract prints: **"the alpha trans-induction factor (alpha TIF, ICP25, VP16, Vmw65)"**

| Gene | Name strings this paper uses |
|---|---|
| UL48 | `alpha trans-induction factor`, `alpha TIF`, `ICP25`, `VP16`, `Vmw65` |
| UL46 | `UL46` |
| UL47 | `UL47` |
| UL23 | `thymidine kinase` |

Conventions: systematic, protein (twice, `ICP25` and `VP16`), functional, molecular
weight, all for one gene, in one 42-character parenthesis.

**Pass condition:** all five UL48 strings resolve to canonical `UL48`, and `ICP25` is not
treated as a misspelling of `ICP27`. `ICP27` is a different gene (UL54).

---

## Test case 2: five genes, three conventions, one abstract

Minagawa H, Tanaka S, Toh Y, Mori R (1994) Detection of herpes simplex virus type
1-encoded RNA by polymerase chain reaction: different pattern of viral RNA detection in
latently infected murine trigeminal ganglia following in vitro or in vivo reactivation.
*J Gen Virol* 75:647-650. PMID 8126462. DOI 10.1099/0022-1317-75-3-647.

**Access: abstract only.**

The abstract prints: **"we detected transcripts for infected cell protein 0 (ICP0)
(Vmw110) 24 h after explantation, but RNAs encoding ICP4 (Vmw175), ICP27, thymidine
kinase and VP16 (ICP25; Vmw65) remained undetectable"**

| Gene | Name strings this paper uses |
|---|---|
| RL2 | `infected cell protein 0`, `ICP0`, `Vmw110` |
| RS1 | `ICP4`, `Vmw175` |
| UL54 | `ICP27` |
| UL23 | `thymidine kinase` |
| UL48 | `VP16`, `ICP25`, `Vmw65` |
| LAT | `Latency-associated transcripts`, `LATs` |

**Pass condition:** five genes resolve from a single sentence. The plural `LATs` resolves
to `LAT`. The semicolon separator inside `(ICP25; Vmw65)` is handled, as is the
double-parenthesis form `ICP0) (Vmw110)`.

This is the only accessed source in the set that glosses `Vmw110` to `ICP0` and `Vmw175`
to `ICP4` in the same sentence.

---

## Test case 3: the IE series glossed to the ICP series

Mullen MA, Gerstberger S, Ciufo DM, Mosca JD, Hayward GS (1995) Evaluation of
colocalization interactions between the IE110, IE175, and IE63 transactivator proteins of
herpes simplex virus within subcellular punctate structures. *J Virol* 69:476-491.
PMID 7983744. DOI 10.1128/jvi.69.1.476-491.1995.

**Access: abstract only.**

The abstract prints: **"three herpes simplex virus-encoded nuclear transactivator
proteins, IE175 (ICP4), IE110 (ICP0), and IE63 (ICP27)"**

| Gene | Name strings this paper uses |
|---|---|
| RS1 | `IE175`, `ICP4` |
| RL2 | `IE110`, `ICP0`, and the typo `1E110` (digit one for letter I) |
| UL54 | `IE63`, `ICP27` |

**Pass condition:** the three IE numbers resolve to three different genes. This case is
also the load-bearing evidence on the `IE63` ambiguity: an HSV-1 paper assigns `IE63` to
`ICP27`, which is UL54, against the `BK012101.1` note that also puts `ie63` on US1. The
dictionary must return UL54 as the established reading and surface the US1 usage as an
ambiguity rather than dropping either.

The `1E110` typo is in the published abstract. It is kept as a dirty-input case: the
dictionary is not required to resolve it, but a pipeline built on the dictionary should
not crash on it.

---

## Test case 4: IE numbering bridged to Roizman's alpha-gene numbering

Roberts MS, Boundy A, O'Hare P, Pizzorno MC, Ciufo DM, Hayward GS (1988) Direct
correlation between a negative autoregulatory response element at the cap site of the
herpes simplex virus type 1 IE175 (alpha 4) promoter and a specific binding site for the
IE175 (ICP4) protein. *J Virol* 62:4307-4320. PMID 2845144.
DOI 10.1128/jvi.62.11.4307-4320.1988.

**Access: abstract only** (title and abstract).

The title alone glosses one gene two ways: **"IE175 (alpha 4) promoter"** and
**"IE175 (ICP4) protein"**. The abstract adds: "its own protein product (IE175 or ICP4)"
and "the herpes simplex virus type 1 IE110 (alpha 0)".

| Gene | Name strings this paper uses |
|---|---|
| RS1 | `IE175`, `alpha 4`, `ICP4` |
| RL2 | `IE110`, `alpha 0` |

Conventions: molecular weight (IE), protein (ICP), and the kinetic-class alpha numbering.

**Pass condition:** `alpha 4` resolves to RS1 and `alpha 0` to RL2. Note that the alpha
digits do not track the ICP digits: `alpha 4` is `ICP4` but `alpha 0` is `ICP0` and
`alpha 27` is `ICP27`, while `alpha 22` is `ICP22`, so the coincidence is partial and must
not be implemented as a rule. Each pairing is a dictionary entry with its own citation.

This convention was not in the M1 brief's list of four. It is attested here, in test case
5 and in test case 11, and the dictionary holds it.

---

## Test case 5: the Roizman lab writing alpha first

Spector D, Purves F, Roizman B (1991) Role of alpha-transinducing factor (VP16) in the
induction of alpha genes within the context of viral genomes. *J Virol* 65:3504-3513.
PMID 1645782. DOI 10.1128/jvi.65.7.3504-3513.1991.

**Access: abstract only.**

The abstract prints: **"the five alpha genes are induced by alpha-transinducing factor
(alpha TIF; VP16), a virion protein"**

| Gene | Name strings this paper uses |
|---|---|
| UL48 | `alpha-transinducing factor`, `alpha TIF`, `VP16` |
| UL54 | `alpha 27 promoter`, `alpha 27` |
| RS1 | `alpha 4 gene`, `alpha 4` |

**Pass condition:** `alpha 27` resolves to UL54 and `alpha 4` to RS1. The spelling
`alpha-transinducing` (closed up, hyphenated) must resolve as well as
`alpha trans-induction` in test case 1 and `alpha trans-inducing` in
Poon AP, Roizman B (1995) *J Virol* 69:7658-7667, PMID 7494274, whose abstract prints
"alpha trans-inducing factor (alpha TIF, VP16, Vmw65)". Three spellings of one functional
name across three papers from overlapping groups.

---

## Test case 6: the 2020 re-annotation and its own naming layer

Whisnant AW, Jurges CS, Hennig T, Wyler E, Prusty B, Rutkowski AJ, L'hernault A,
Djakovic L, Gobel M, Doring K, Menegatti J, Antrobus R, Matheson NJ, Kunzig FWH,
Mastrobuoni G, Bielow C, Kempa S, Liang C, Dandekar T, Zimmer R, Landthaler M, Grasser F,
Lehner PJ, Friedel CC, Erhard F, Dolken L (2020) Integrative functional genomics decodes
herpes simplex virus 1. *Nat Commun* 11:2038. PMID 32341360. PMCID PMC7184758.
DOI 10.1038/s41467-020-15992-5.

**Access: full text read** (JATS XML from Europe PMC). Supplementary Data 1 to 9 were not
retrieved: they are separate files not carried in the XML.

This paper is both a test case and the origin of the second inventory, so two things are
tested against it.

**First, its stated relationship to the classical names, quoted verbatim:**

> "we fully maintained the current nomenclature for all ORFs in the reference annotation
> and attributed each ORF to the most highly expressed transcript in its vicinity"

So the classical UL/US names are kept, not renumbered. The new naming is additive.

**Second, its own naming rules for novel features, quoted verbatim:**

> "Alternative transcript isoforms initiating within less than 500 nucleotides were
> labeled with "*" (extended) or "#" (truncated), e.g., UL13 mRNA *1. Finally,
> alternatively spliced transcripts were labeled with iso 1 and iso 2. Short ORFs
> (<100 aa), were named upstream ORF (uORF), upstream overlapping ORF (uoORF), internal
> ORF (iORF) and downstream ORF (dORF) in relation to the next neighboring large ORF. Any
> ORF for which no transcript could be identified to be responsible for its translation
> was labeled as orphan."

> "In accordance, any transcript, which was not found to encode an ORF or sORF within its
> first 500 nt was also labeled as orphan."

The figure 6d legend adds the ORF-type vocabulary: "CDS coding sequence, NTE N-terminal
extended ORFs, NTT N-terminal truncated ORFs, uORF upstream ORF, uoORF upstream
overlapping ORF, iORF internal ORF, dORF downstream ORF, sORF short ORF".

**Third, its cross-convention glosses:**

| Gene | Name strings this paper uses | Gloss as printed |
|---|---|---|
| RL2 | `RL2`, `ICP0`, `RL2A`, `RL2 iso1` | "The RL2 locus encodes the major viral immediate early protein ICP0" |
| RL1 | `RL1`, `ICP34.5`, `RL1A` | "Within the ICP34.5 (RL1) locus" |
| UL54 | `UL54`, `ICP27` | "the major immediate early protein ICP27 (UL54)" |
| UL19 | `UL19`, `VP5`, `MCP` | "the major capsid protein (VP5, UL19)" |
| UL35 | `UL35`, `VP26` | "mCherry on the N-terminus of the UL35 gene product (VP26)" |
| RS1 | `ICP4` only | `RS1` does not appear in the main text |
| US1, US12 | `ICP22`, `ICP47` | "the spliced immediate early (IE) genes ICP0, ICP22 and ICP47" |
| UL1 | `glycoprotein L` | "fusion proteins (e.g., between ICP0 and glycoprotein L)" |
| UL49A | `UL49.5` | "previously identified CDS (UL11, UL49.5, US5, US9, US12, and RL2 iso1)" |
| sequence | `JN555585` | "HSV-1 ( JN555585 )" |

**Pass condition:** every classical name in the table resolves unchanged, `UL49.5`
resolves to the same gene as `UL49A`, and the `inventory` column distinguishes names this
paper added from names it inherited. `VP16`, `alpha-TIF` and every `Vmw` string are
absent from this paper's main text: the dictionary must not claim it as a source for them.

`UL15.4` and `UL20.6` do not occur in the main text either. They are in the deposited
annotation (`BK012101.1`, the Zenodo GTF) and are cited to those, not to the paper.

---

## Test case 7: systematic to VP, glossed in both directions

Jin F, Li S, Zheng K, Zhuo C, Ma K, Chen M, Wang Q, Zhang P, Fan J, Ren Z, Wang Y (2014)
Silencing herpes simplex virus type 1 capsid protein encoding genes by siRNA: a promising
antiviral therapeutic approach. *PLoS One* 9:e96623. PMID 24794394.
DOI 10.1371/journal.pone.0096623.

**Access: full text read.**

| Gene | Name strings this paper uses | Gloss as printed |
|---|---|---|
| UL18 | `UL18`, `VP23` | both directions: "capsid proteins VP23 (UL18)" and "UL18 (VP23) and UL19 (VP5) can be potential drug targets"; tables print `UL18(VP23)` |
| UL19 | `UL19`, `VP5` | "VP5 (UL19)"; "a major capsid protein (VP5)"; tables print `UL19(VP5)` |
| UL35 | `UL35`, `VP26` | tables print `UL35(VP26)` |
| UL38 | `UL38`, `VP19C` | "a null mutant of UL38 (VP19C)"; tables print `UL38(VP19C)` |
| UL26.5 | `UL26.5`, `ICP35` | tables print `UL26.5(ICP35)`; "the scaffolding protein (product of the UL26.5 gene)" |
| UL26 | `UL26`, `VP24`, `VP21` | "the two products encoded by UL26, VP24 and VP21" |
| UL6 | `pUL6` | "a dodecamer of pUL6" |
| UL23 | `TK`, `Thymidine Kinase` | "Thymidine Kinase (TK)" |
| others | `ICP4`, `ICP6`, `UL5`, `UL8`, `UL9`, `UL12`, `UL29`, `UL30`, `UL39`, `UL42`, `UL52` | |

**Pass condition:** the no-space table form `UL18(VP23)` resolves as readily as the
spaced form, the `pUL6` protein prefix resolves to UL6, and `UL26` keeps `VP24`/`VP21`
while `UL26.5` keeps `ICP35`. A dictionary that merges UL26 and UL26.5 fails this case.

---

## Test case 8: glycoprotein letters, with asymmetric glossing

Jambunathan N, Clark CM, Musarrat F, Chouljenko VN, Rudd J, Kousoulas KG (2021) Two Sides
to Every Story: Herpes Simplex Type-1 Viral Glycoproteins gB, gD, gH/gL, gK, and Cellular
Receptors Function as Key Players in Membrane Fusion. *Viruses* 13:1849. PMID 34578430.
DOI 10.3390/v13091849.

**Access: full text read.**

| Gene | Name strings this paper uses | Gloss as printed |
|---|---|---|
| UL1 | `gL`, `UL1` | "Glycoprotein gL, a 224aa glycoprotein encoded by the UL1 gene"; "UL1 (gL)" |
| UL22 | `gH`, `UL22` | "UL22 (gH)"; "Glycoprotein H (gH)" |
| UL53 | `gK`, `UL53` | "The UL20 and UL53 (gK) genes are among the most conserved"; "Glycoprotein K (gK)" |
| UL27 | `gB` only | 97 occurrences of `gB`, and `UL27` never appears |
| US6 | `gD` only | 51 occurrences of `gD`, and no `US` name appears anywhere in the paper |
| UL44, US8, US7, UL10 | `gC`, `gE`, `gI`, `gM` | "including gE, gI, gM, gK, and the UL20 and UL45 proteins" |
| UL20, UL45, UL11, UL16, UL21 | systematic only | "the tegument proteins UL11, UL16, and UL21 bind to the cytoplasmic tail of gE" |

**Pass condition:** the glycoprotein letters resolve without needing the paper to have
glossed them. This is the realistic case: the same paper that spells out `UL1 (gL)` never
gives `gB` a UL number, so a dictionary that only learns pairings from in-paper glosses
would resolve gL and fail on gB. The dictionary must carry `gB` to UL27 from a database
record instead.

---

## Test case 9: a modern clinical paper, systematic and functional together

Ibrahim R, Fidouh N, Bunel V, Bouscarat F, Issaurat P, Fernandez J, Joly V, Rioux C,
Boutolleau D, Ghosn J (2026) Clinical and genotypic characteristics of antiviral-resistant
herpes simplex virus 1 and 2 infections: a case series. *BMC Infect Dis* 26:1072.
PMID 42015059. DOI 10.1186/s12879-026-13310-w.

**Access: full text read.**

| Gene | Name strings this paper uses | Gloss as printed |
|---|---|---|
| UL23 | `UL23` (40 times), `TK` (23), `thymidine kinase` (7) | both directions: "sequencing the UL23 (thymidine kinase) and UL30 (DNA polymerase) genes", "those in UL23 (TK)", and "the viral thymidine kinase (TK), which is encoded by the UL23 gene" |
| UL30 | `UL30` (36 times), `DNA polymerase`, `DNA pol` | "mutations in the gene encoding HSV DNA pol (UL30)" |

**Pass condition:** three distinct string forms for one gene (`UL23`, `TK`,
`thymidine kinase`) all resolve, and the abbreviated `DNA pol` resolves to UL30. This is
the case that matters most for downstream use, because resistance genotyping is where a
coordinate and name layer gets consumed in practice.

---

## Test case 10: all four conventions plus the kinetic class and the miRNA layer

Duan Y, Zeng J, Fan S, Liao Y, Feng M, Wang L, Zhang Y, Li Q (2019) Herpes Simplex Virus
Type 1-Encoded miR-H2-3p Manipulates Cytosolic DNA-Stimulated Antiviral Innate Immune
Response by Targeting DDX41. *Viruses* 11:756. PMID 31443275. DOI 10.3390/v11080756.

**Access: full text read.**

| Gene | Name strings this paper uses | Gloss as printed |
|---|---|---|
| RL2 | `ICP0` | "HSV-1-encoded infected-cell polypeptide 0 (ICP0) is an immediate early (IE, alpha) protein" |
| UL41 | `VHS`, `vhs`, `UL41` | "HSV-1 virion host shutoff (VHS) protein encoded by UL41 gene" |
| UL23 | `TK`, `thymidine kinase` | "HSV-1 thymidine kinase (TK), glycoprotein C (gC), and ICP0 expression" |
| LAT | `LAT` | "regulation of the latency-associated transcript (LAT)" |
| RL1 | `ICP34.5`, `gamma1 34.5` | the same gene named both ways in neighbouring sentences: "STING by ICP27, gamma1 34.5 and UL46" |
| miRNAs | `miR-H2`, `miR-H1`, `miR-H4`, `miR-H6`, `miR-H8`, `miR-H27`, `miR-H28`, `miR-H29`, `miR-H2-3p`, `miR-H4-5p`, `miR-H6-3p` | "miR-H2-3p, miR-H4-5p, and miR-H6-3p target viral ICP0, ICP34.5, and ICP4 mRNA" |
| others | `UL30`, `UL37`, `UL46`, `US3`, `US11`, `UL36USP`, `ICP4`, `ICP27`, `ICP47`, `VP22`, `gC` | "cGAS by UL37, UL41, and VP22" |

**Pass condition:** `VHS` and `vhs` both resolve to UL41, so matching is case-insensitive,
and `gamma1 34.5` resolves to RL1 alongside `ICP34.5`. The `miR-H` names are recorded in
the dictionary only where a source attests them (see `alias_notes.md` on the miRNA scope
decision); this case is where that gap is visible.

---

## Test case 11: the earliest accessible Vmw-to-ICP gloss

Everett RD (1984) Trans activation of transcription by herpes virus products: requirement
for two HSV-1 immediate-early polypeptides for maximum activity. *EMBO J* 3:3135-3141.
PMID 6098466. PMCID PMC557829.

**Access: abstract read in full; the PMC body is a scanned PDF and was searched only
through the PMC article page, which exposes abstract and references.**

The abstract prints: **"the HSV-1 IE polypeptide Vmw175 (ICP 4) is essential"** and
**"Vmw175 plus the product of IE gene 1, Vmw110 (ICP 0), is a much better activator"**

| Gene | Name strings this paper uses |
|---|---|
| RS1 | `Vmw175`, `ICP 4` |
| RL2 | `Vmw110`, `ICP 0`, `IE gene 1` |
| IE transcripts | `IE RNA 1`, `IE RNA 2`, `IE RNA 3`, `IE RNA 4` |

**Pass condition:** the spaced forms `ICP 4` and `ICP 0` resolve as readily as `ICP4` and
`ICP0`. Note the IE *gene* numbering (`IE gene 1` to `IE gene 5`) is a sixth numbering,
separate from the IE *molecular weight* numbering (`IE110`, `IE175`, `IE63`). This paper
uses the gene numbering; test cases 3 and 4 use the weight numbering. Only three of the
five Vmw numbers appear anywhere in the paper, so it is not a source for `Vmw63`,
`Vmw68` or `Vmw12`.

---

## Test case 12: the ICP35 family and the VP22a band, in one abstract

Braun DK, Roizman B, Pereira L (1984) Characterization of post-translational products of
herpes simplex virus gene 35 proteins binding to the surfaces of full capsids but not
empty capsids. *J Virol* 49:142-153. PMID 6317887. DOI 10.1128/jvi.49.1.142-153.1984.

**Access: abstract only.**

The abstract prints: **"a genetically and immunologically related family of structural
(gamma) polypeptides of herpes simplex virus 1 designated as infected cell polypeptides
(ICP) 35"**, then **"One of these corresponds to virion protein 22a, a component of
capsids containing DNA"**, and **"Surface iodination labeled the major capsid protein
(ICP5) of empty capsids"** plus **"ICP5, the major capsid protein"**.

| Gene | Name strings this paper uses |
|---|---|
| UL26.5 | `infected cell polypeptides (ICP) 35`, `ICP35`, `ICP35a` to `ICP35f`, `virion protein 22a` |
| UL19 | `ICP5`, `major capsid protein` |

**Pass condition:** `ICP5` resolves to UL19. This is the single accessible attestation of
`ICP5` located anywhere: the string is in neither UniProt, nor `NC_001806.2`, nor
`BK012101.1`, and it has no title hits in Europe PMC. Without this paper, `ICP5` would be
an unverified candidate.

The same abstract is also why `VP22a` sits on UL26.5 rather than UL26: the ICP35 family it
describes is the UL26.5 product, and one of its members is the `virion protein 22a` band.
It demonstrates the gel-band problem directly, resolving one ICP number into six bands
(`ICP35a` to `ICP35f`) and about twenty two-dimensional spots.

---

## Test case 13: VP to systematic, with the paper stating its own naming policy

Villanueva-Valencia JR, Tsimtsirakis E, Evilevitch A (2021) Role of HSV-1 Capsid
Vertex-Specific Component (CVSC) and Viral Terminal DNA in Capsid Docking at the Nuclear
Pore. *Viruses* 13:2515. PMID 34960783. DOI 10.3390/v13122515.

**Access: full text read.**

This paper declares its convention in a parenthesis, which is the rarest and most useful
thing a paper can do for an alias dictionary:

> "by tether tegument protein UL36, which binds Nup358 localized on the cytoplasmic NPC
> filaments **(gene names are used to refer to their respective proteins)**"

| Gene | Name strings this paper uses | Gloss as printed |
|---|---|---|
| UL19 | `VP5`, `UL19` | "the major capsid protein VP5 (UL19)" |
| UL18 | `VP23`, `UL18` | "two VP23 (UL18) proteins and one VP19C (UL38) protein" |
| UL38 | `VP19C`, `UL38` | same sentence |
| UL6 | `UL6` | "a dodecameric ring of UL6 protein"; "portal forming proteins (UL6)" |
| UL25, UL17, UL36 | `UL25`, `UL17`, `UL36` | "a singlet of UL17 and duplexes of both UL25 and UL36" |
| UL35 | `VP26` | present, with no UL35 gloss |
| UL48 | `VP16` | one occurrence |

**Pass condition:** the systematic names resolve when used as protein names, which is
what the declared policy means in practice.

---

## Also attested, held as supporting rather than as separate test cases

These were read and are cited in `aliases_literature.tsv`, but they duplicate a
convention already covered above.

| Paper | Access | What it adds |
|---|---|---|
| Hagmann M, Georgiev O, Schaffner W (1997) *J Virol* 71:5952-5962, PMID 9223485 | abstract only | "VP16 (Vmw65, alphaTIF, ICP25)" - the closed-up `alphaTIF` spelling |
| Weinheimer SP, Boyd BA, Durham SK, Resnick JL, O'Boyle DR (1992) *J Virol* 66:258-269, PMID 1309245 | abstract only | "VP16 (also called Vmw65 and alpha TIF)" - a prose gloss marker, which a regex tuned to `(A, B, C)` misses |
| Liu Y, Gong W, Huang CC, Herr W, Cheng X (1999) *Genes Dev* 13:1692-1703, PMID 10398682 | abstract only | "VP16 (Vmw65; alphaTIF)" - semicolon separator |
| McKnight JL, Doerr M, Zhang Y (1994) *J Virol* 68:1750-1757, PMID 8107236 | abstract only | "alpha trans-induction factor (alpha TIF or VP16)" - an `or` separator, and the `VP13/14` slash form in the title |
| Everett R, O'Hare P, O'Rourke D, Barlow P, Orr A (1995) *J Virol* 69:7339-7344, PMID 7474166 | abstract only | "Vmw110 (also known as ICP0)" |
| Allen KE, Everett RD (1997) *J Gen Virol* 78:2913-2922, PMID 9367379 | abstract only | "immediate early protein Vmw175 (ICP4)" |
| Zhang Y, McKnight JL (1993) *J Virol* 67:1482-1492, PMID 8382306 | abstract only | `VP11`, `VP12`, `VP11/12` for UL46 and `VP13`, `VP14`, `VP13/14` for UL47: the slash forms, which are a distinct string shape |
| Ace CI, McKee TA, Ryan JM, Cameron JM, Preston CM (1989) *J Virol* 63:2260-2269, PMID 2539517 | abstract only | "the transinducing factor Vmw65", and the `IE RNA 1` to `IE RNA 4` numbering |
| Schwyzer M, Wirth UV, Vogt B, Fraefel C (1994) *J Gen Virol* 75:1703-1711, PMID 8021599 | abstract only | "homology to ICP22 (Vmw68) of herpes simplex virus type 1": the only located attestation of `Vmw68`, in a bovine herpesvirus 1 paper |
| Singh M, Fraefel C, Bello LJ, Lawrence WC, Schwyzer M (1996) *J Gen Virol* 77:615-625, PMID 8627249 | abstract only | "homology to ICP27 (Vmw63) of herpes simplex virus 1": likewise a bovine herpesvirus 1 paper |
| Preston VG, Rixon FJ, McDougall IM, McGregor M, al Kobaisi MF (1992) *Virology* 186:87-98, PMID 1309284 | abstract only | `ICP35`, `UL26` and `Vmw175` in one abstract: three conventions, and it states "the UL26 gene specified a protein distinct from ICP35" |
| Liu FY, Roizman B (1991) *J Virol* 65:5149-5156, PMID 1654435 | abstract only | "the UL26.5 gene, which specifies a capsid protein designated infected cell protein 35" |
| Umbach JL, Kramer MF, Jurak I, Karnowski HW, Coen DM, Cullen BR (2008) *Nature* 454:780-783, PMID 18596690 | abstract only | "the latency associated transcript (LAT)" unhyphenated, with `miR-H2-3p`, `miR-H6`, `ICP0`, `ICP4` |
| Cui C, Griffiths A, Li G, Silva LM, Kramer MF, Gaasterland T, Wang XJ, Coen DM (2006) *J Virol* 80:5499-5508, PMID 16699030 | abstract only | "the latency-associated transcript (LAT)" hyphenated. With Umbach 2008, this is the pair that proves both spellings are in print |
| Jovasevic V, Roizman B (2010) *Virol J* 7:103, PMID 20492679 | title and abstract | "a domain encoding US5, US4, US3, US2 and alpha22": the alpha numbering applied to US1, and the `US5-1` hyphenated form against the re-annotation's `US5.1` |
| Carter KL, Roizman B (1996) *J Virol* 70:172-178, PMID 8523523 | title and abstract | "the open reading frame of the alpha 22 gene": the primary source for US1.5 |
| Hubenthal-Voss J, Starr L, Roizman B (1987) *J Virol* 61:3349-3355, PMID 2822945 | title and abstract | Ori-S as a transcribed open reading frame; cited by name inside the `BK012101.1` note |
| Perng GC, Maguen B, Jin L, Mott KR, Kurylo J, BenMohamed L et al. (2002) *J Virol* 76:8003-8010, PMID 12134005 | title and abstract | `AL-RNA` antisense to the 5' end of LAT; the `BK012101.1` note cites it and records the feature as not observed |

---

## What these cases prove, and what they do not

**They prove** the dictionary resolves real printed strings, including: four conventions
inside one parenthesis; glosses in both directions; four separator styles (comma,
semicolon, `or`, "also called"); three spellings of `alpha TIF`; the no-space table form
`UL18(VP23)`; the `pUL` protein prefix; the slash forms `VP11/12` and `VP13/14`; the
plural `LATs`; spaced ICP numbers; case variation between `VHS` and `vhs`; both
hyphenations of latency-associated transcript; and two numbering systems that share
digits with neither the ICP series nor each other.

**They do not prove** anything about papers from 1972 to 1988 that established these
conventions in the first place. Every one of the foundational papers is either paywalled
or available only as a page scan, so the test set is weighted towards what is accessible
rather than what is authoritative. `alias_notes.md` records which sources could not be
reached and what that costs.
