# Published coordinates resolved by hand

Derived by `src/build_published_coordinates.py` from `data/published_coordinates.tsv` and
`data/published_coordinate_checks.tsv`. Do not edit: `make published-coordinates` rewrites it.

This is the record of M3's one criterion a test run cannot settle by itself.
The criterion, written in the build plan before any code:

> Correctly resolves a published coordinate from at least five papers across at least three strains to the canonical reference, checked by hand against each paper.

**Met:** 6 papers against the 5 required, across 4 strains against the 3 required (17, F, KOS, SC16).
24 papers were opened in all. Every resolution is checked by 131 mechanical checks over the
deposited records, the genome FASTA and the built artefacts, of which 131 pass.

**What is transcription and what is checked.** The published coordinate, the quoted sentence and
the fetch route in each row are **hand transcription** from the paper's full text, flagged in the
`transcribed` column the same way `data/conventions.tsv` flags its `paper_transcribed` rows.
Nothing in the suite re-reads a paper: continuous integration has no network. Everything
downstream of the transcription is checked, and no expected value in the checks table was
produced by calling the function it checks.

**Quotations are literal slices** of the fetched bytes, cut by string search rather than retyped,
so a subscript in the publisher's markup can come through as a space: Lee et al.'s "both TR L
repeats" is how that sentence reads once the tags are stripped.

## The counted cases

| Paper | Strain | Accession as the paper states it | How the frame is known | Published | Resolved in NC_001806.2 | Checks |
|---|---|---|---|---|---|---|
| Amelio AL, McAnany PK, Bloom DC. A chromatin insulator-like element in the herpes simplex virus type 1 latency-associated transcript region binds CCCTC-binding factor and displays enhancer-blocking and silencing activities. Journal of virology 2006;80:2358-2368 | 17 | NC_001806 | determined from paper claim | 120503-120635 (CTCCC cluster); 120461-120665 (ChIP amplicon); 117067-117156 and 96331-96517 (two further primer pairs) | 120505-120637, not well defined: the counterpart interval is 5707-5911 on the opposite strand. Amplicon 120463-120667 | 24 |
| Creech CC, Neumann DM. Changes to euchromatin on LAT and ICP4 following reactivation are more prevalent in an efficiently reactivating strain of HSV-1. PloS one 2010;5:e15416 | 17 | X14112.1 | stated with version | 119326-119397 (LAT 5' exon amplicon); 120797-120820 (LAT I-1 primer, printed as a 25-mer) | 119327-119398; the primer start resolves to 120799, not well defined, counterpart 5575 on the opposite strand | 9 |
| Pinnoji RC, Bedadala GR, George B, Holland TC, Hill JM, Hsia SC. Repressor element-1 silencing transcription factor/neuronal restrictive silencer factor (REST/NRSF) can regulate HSV-1 immediate-early transcription via histone modification. Virology journal 2007;4:56 | 17 | X14112 | stated accession held | 132082-132103 (RE-1/NRSE) | 132085-132106, not well defined: the counterpart interval is 146091-146112 on the opposite strand, where the same six bases are annotated as the US12 TATA box | 9 |
| Singh N, Zachariah S, Phillips AT, Tscharke D. Lytic promoter activity during herpes simplex virus latency is dependent on genome location. Journal of virology 2024;98:e0125824 | KOS | JQ673480 | stated accession held | 55985-56282 (gB promoter); arm boundaries 11610/11611, 53024/53025, 116144/116145 | 56044-56341; 11662/11663; 53083/53084; 116224/116225 | 26 |
| Tongmuang N, Krishnan M, Connor V, Crump C, Jensen LE. UL56 Is Essential for Herpes Simplex Virus-1 Virulence In Vivo but Is Dispensable for Induction of Host-Protective Immunity. Vaccines 2024;12:837 | SC16 | MN159383.1 | stated with version | 115941-116645 (UL56 CDS); 114423-118177 (targeting arm region) | 116222-116926; the arm's outer end 118177 resolves to 118446 but is declined, counterpart 7928 | 7 |
| Kitabatake T, Ishioka K, Suzuki T, Miyazaki N, Suzutani T. Attenuated neurovirulence of dUTPase-deficient herpes simplex virus type 1 in suckling mice but not in newborn mice. Fukushima journal of medical science 2026;72:131-140 | F | GU734771.1 | stated with version | 106888-106947 and 106908-106967 (two primer windows); the substitution 106918G>A, called the 4th codon of UL50 (W4Ter) | 106992-107051 and 107012-107071; the substitution at 107022. Well defined, with a uniform +104 offset across the window | 14 |

### Why each counted case's frame is known

- **amelio2006_ctrl2** (determined_from_paper_claim): Six primers the paper prints each occur exactly once and start or end on their printed coordinate in the 152,261 bp record, at three loci, and shift by +1 or +2 in NC_001806.2; so the paper's bare NC_001806 is the .1 version
- **creech2010_lat5exon** (stated_with_version): The figure legends name X14112.1 with its version suffix
- **pinnoji2007_nrse** (stated_accession_held): The corpus holds exactly one version of X14112, 152,261 bp and identical to NC_001806.1 base for base; the overlap with the ICP22 (US1) TATA box that the paper states holds base-exact in that frame, where the TATA box is the last six bases of the published interval, and does not hold in NC_001806.2
- **singh2024_kos_arms** (stated_accession_held): The corpus holds exactly one version of JQ673480; four guide-RNA protospacers the paper prints each occur exactly once in it and once in the reference, and in every case the reference occurrence is where the library lifts the deposit occurrence to
- **tongmuang2024_ul56** (stated_with_version): The paper writes the accession with its version suffix in the same clause as the coordinates
- **kitabatake2026_ul50** (stated_with_version): The accession is written with its version suffix in the same sentence as the coordinates, and again in Methods. The printed revertant primer reproduces all 60 bases of GU734771.1 at the paper's own start coordinate with no differences, and the printed mutagenic primer differs from the record at exactly one base, which is the position the paper names

## Cases deliberately not counted

Recorded because a clean run of five is not what happened, and because a criterion is worth
nothing if the cases that failed it are not visible.

| Paper | Strain | Status | Why it is not counted |
|---|---|---|---|
| Lassalle F, Beale MA, Bharucha T, Williams CA, Williams RJ, Cudini J, et al | 17 | supporting | The published coordinates are already in the canonical frame, so the crosswalk does no work and this cannot count toward a criterion about resolving a coordinate; kept because it confirms this repository's own transcription of the paper and because all five boundaries turn out to be gene or CDS edges rather than the repeat extents the paper's own caption describes |
| Jiao X, Sui H, Lyons C, Tran B, Sherman BT, Imamichi T | McKrae | supporting | The paper does not state which record its number is in, and the two candidates agree here, so the answer is right either way but the frame is not stated and the lift does no work; counting it would count luck |
| Lee JS, Raja P, Pan D, Pesola JM, Coen DM, Knipe DM | KOS | supporting | The library declines, which is the correct answer and is what the fourth M3 criterion requires, but a decline is not a resolution so it cannot count toward a criterion about resolving a coordinate; the paper itself says the deletion is in both repeat copies and the two reference copies are identical across the window |
| Singh P, Collins MF, Johns RN, Manuel KA, Ye ZA, Bloom DC, et al | 17 | not_counted | The frame cannot be settled from the papers, and the published endpoints are internally inconsistent: 120,500 to 120,635 inclusive is 136 bases against the stated 135, and the same group's 2024 figure legend prints 120,503 for the same element |
| Chau VQ, Kolb AW, Miller DL, Yannuzzi NA, Brandt CR | BP-K15 | not_counted | Frame not stated, and the paper's own two numbers disagree: 126,983 to 132,514 inclusive is 5,532 bp, described in the same sentence as an 8,511 bp segment. BP-K15 is also not one of the 183 corpus genomes |
| Colgrove RC, Liu X, Griffiths A, Raja P, Deluca NA, Newman RM, et al | KOS | not_counted | Frame not stated. OriS is also genuinely present in two copies, so a single number would be incomplete even with the frame known |
| Szpara ML, Parsons L, Enquist LW | F, H129 | not_counted | No nucleotide coordinate is printed anywhere in the article text, so there is nothing to resolve. Worth recording because this is the source paper for two of the corpus genomes |
| Dremel SE, Koparde VN, Arbuckle JH, Hogan CH, Kristie TM, Krug LT, et al | KOS | supporting | The library declines, correctly: both endpoints are inside the long inverted repeat and the two reference copies are identical across the window. A decline is not a resolution, so it cannot count, but it is the second independent case confirming that lifting an interval's endpoints does not preserve its length: 678 bases in the deposit against 677 between the two correctly lifted ends, because a homopolymer differs |
| Lee CK, Ng SY, Chai CN, Lim YF, Hu TJ, Lee OF, et al | clinical isolate | supporting | The published coordinates are already in the canonical frame, so the crosswalk does no work and this cannot count toward a criterion about resolving a coordinate. Kept because it is the clearest example in the set of a paper citing its reference version correctly, and because the window checks out as a real US4 interval |
| Morse AM, Calabro KR, Fear JM, Bloom DC, McIntyre LM | 17 | not_counted | No usable coordinate in the article text |
| Chang W, Hao M, Qiu J, Sherman BT, Imamichi T | 17 | not_counted | No usable coordinate in the article text |
| Bellizzi A, Çakır S, Donadoni M, Sariyer R, Liao S, Liu H, et al | 17 | not_counted | No usable coordinate in the article text |
| Whisnant AW, Jürges CS, Hennig T, Wyler E, Prusty B, Rutkowski AJ, et al | 17 | not_counted | No usable coordinate in the article text |
| Boldogkői Z, Szűcs A, Balázs Z, Sharon D, Snyder M, Tombácz D | 17 | not_counted | No usable coordinate in the article text |
| Jurak I, Kramer MF, Mellor JC, van Lint AL, Roth FP, Knipe DM, et al | 17 | not_counted | Not retrieved. Named in this repository's own round 18 note as the miRNA case that motivated the WrongReferenceVersion refusal, so its coordinates would have been worth having |
| Umbach JL, Kramer MF, Jurak I, Karnowski HW, Coen DM, Cullen BR | 17 | not_counted | Not retrieved. Named alongside Jurak 2010 in the same round 18 note |
| Washington SD, Singh P, Johns RN, Edwards TG, Mariani M, Frietze S, et al | 17 | not_counted | Not retrieved. Singh 2022 attributes its delta-CTRL2 virus to this paper, so it is the link that would settle whether the 120,500 figure originates here and in which frame |
| Washington SD, Musarrat F, Ertel MK, Backes GL, Neumann DM | 17 | not_counted | Not retrieved. The companion paper to Washington 2019, for the same reason |

## Every check, and what decides it

### amelio2006_ctrl2

| Check | Kind | Decided against | Result |
|---|---|---|---|
| convert1 | version_convert | convert(120503, NC_001806.1, NC_001806.2) = 120505 | pass |
| convert2 | version_convert | convert(120635, NC_001806.1, NC_001806.2) = 120637 | pass |
| convert3 | version_convert | convert(120461, NC_001806.1, NC_001806.2) = 120463 | pass |
| convert4 | version_convert | convert(120665, NC_001806.1, NC_001806.2) = 120667 | pass |
| convert5 | version_convert | convert(117067, NC_001806.1, NC_001806.2) = 117068 | pass |
| convert6 | version_convert | convert(96331, NC_001806.1, NC_001806.2) = 96332 | pass |
| primer_ctrl2_f_old | sequence_unique_at | X14112.1:120461 = CTCTGTGGTTAACACCAGAG | pass |
| primer_ctrl2_f_new | sequence_unique_at | NC_001806.2:120463 = CTCTGTGGTTAACACCAGAG | pass |
| primer_ctrl2_r_old | sequence_unique_at | X14112.1:120646 = GCGATACATCCAACACAGAC | pass |
| primer_ctrl2_r_new | sequence_unique_at | NC_001806.2:120648 = GCGATACATCCAACACAGAC | pass |
| primer_lat_f_old | sequence_unique_at | X14112.1:117067 = GCATGCGTCGCCCAAC | pass |
| primer_lat_f_new | sequence_unique_at | NC_001806.2:117068 = GCATGCGTCGCCCAAC | pass |
| primer_lat_r_old | sequence_unique_at | X14112.1:117137 = GATCACATGCAATCTAACTG | pass |
| primer_lat_r_new | sequence_unique_at | NC_001806.2:117138 = GATCACATGCAATCTAACTG | pass |
| primer_gc_f_old | sequence_unique_at | X14112.1:96331 = CCTTGCCGTGGTCCTGTGGA | pass |
| primer_gc_f_new | sequence_unique_at | NC_001806.2:96332 = CCTTGCCGTGGTCCTGTGGA | pass |
| primer_gc_r_old | sequence_unique_at | X14112.1:96498 = CATCGACCCCAAACCCCAAC | pass |
| primer_gc_r_new | sequence_unique_at | NC_001806.2:96499 = CATCGACCCCAAACCCCAAC | pass |
| amplicon_counterpart_start | counterpart_pair | NC_001806.2:120463=C against 5911=G; complementary True; repeat map alternatives (5911,) | pass |
| amplicon_counterpart_end | counterpart_pair | NC_001806.2:120667=C against 5707=G; complementary True; repeat map alternatives (5707,) | pass |
| declines | lift_declines | declined; detail 120463\|5911, wanted 120463\|5911 | pass |
| reiteration_set_4 | feature_location | ref_NC_001806.2.gb repeat_region/note=RL reiteration set 4 -> 5732..5878\|120496..120642 | pass |
| lat_gene | feature_location | ref_NC_001806.2.gb gene/gene=LAT -> complement(<1..7569)\|118777..127151 | pass |
| lat_ncrna | feature_location | ref_NC_001806.2.gb ncRNA/gene=LAT -> complement(join(<1..4953,6910..7569))\|join(118805..119464,121421..127151) | pass |

### creech2010_lat5exon

| Check | Kind | Decided against | Result |
|---|---|---|---|
| convert1 | version_convert | convert(119326, NC_001806.1, NC_001806.2) = 119327 | pass |
| convert2 | version_convert | convert(119397, NC_001806.1, NC_001806.2) = 119398 | pass |
| convert3 | version_convert | convert(120797, NC_001806.1, NC_001806.2) = 120799 | pass |
| primer_old | sequence_unique_at | X14112.1:120797 = GACACGGATTGGCTGGTGTAGTGGG | pass |
| primer_new | sequence_unique_at | NC_001806.2:120799 = GACACGGATTGGCTGGTGTAGTGGG | pass |
| counterpart | counterpart_pair | NC_001806.2:120799=G against 5575=C; complementary True; repeat map alternatives (5575,) | pass |
| declines | lift_declines | declined; detail 120799\|5575, wanted 120799\|5575 | pass |
| rl2_mrna | feature_location | ref_NC_001806.2.gb mRNA/gene=RL2 -> join(2114..2318,3084..3750,3887..5699)\|complement(join(120675..122487,122624..123290,124056..124287)) | pass |
| no_lat_in_old_record | feature_absent | ref_NC_001806.1.gb holds 0 ncRNA feature(s) | pass |

### pinnoji2007_nrse

| Check | Kind | Decided against | Result |
|---|---|---|---|
| convert1 | version_convert | convert(132082, NC_001806.1, NC_001806.2) = 132085 | pass |
| convert2 | version_convert | convert(132103, NC_001806.1, NC_001806.2) = 132106 | pass |
| element_old | sequence_unique_at | X14112.1:132082 = GGGTCTCTCCGGCGCACATAAA | pass |
| element_new | sequence_unique_at | NC_001806.2:132085 = GGGTCTCTCCGGCGCACATAAA | pass |
| element_counterpart_copy | sequence_unique_at | NC_001806.2:146091 = TTTATGTGCGCCGGAGAGACCC | pass |
| counterpart | counterpart_pair | NC_001806.2:132085=G against 146112=C; complementary True; repeat map alternatives (146112,) | pass |
| tata_old | feature_location | ref_NC_001806.1.gb regulatory/gene=US1 -> 132098..132103\|133941..133946 | pass |
| tata_new | feature_location | ref_NC_001806.2.gb regulatory/gene=US1 -> 132101..132106\|133944..133949 | pass |
| tata_other_copy | feature_location | ref_NC_001806.2.gb regulatory/gene=US12 -> complement(146091..146096) | pass |

### singh2024_kos_arms

| Check | Kind | Decided against | Result |
|---|---|---|---|
| lift1 | lift | lift(55985, JQ673480.1, NC_001806.2) = 56044 | pass |
| lift2 | lift | lift(56282, JQ673480.1, NC_001806.2) = 56341 | pass |
| lift3 | lift | lift(11610, JQ673480.1, NC_001806.2) = 11662 | pass |
| lift4 | lift | lift(11611, JQ673480.1, NC_001806.2) = 11663 | pass |
| lift5 | lift | lift(53024, JQ673480.1, NC_001806.2) = 53083 | pass |
| lift6 | lift | lift(53025, JQ673480.1, NC_001806.2) = 53084 | pass |
| lift7 | lift | lift(116144, JQ673480.1, NC_001806.2) = 116224 | pass |
| lift8 | lift | lift(116145, JQ673480.1, NC_001806.2) = 116225 | pass |
| spacer_ul26_27_kos | sequence_unique_at | JQ673480.1:52856 = TTTGTCACGGGAAAGGAAAG | pass |
| spacer_ul26_27_ref | sequence_unique_at | NC_001806.2:52919 = TTTGTCACGGGAAAGGAAAG | pass |
| spacer_ul26_27_lift | lift | lift(52856, JQ673480.1, NC_001806.2) = 52919 | pass |
| spacer_ul27_kos | sequence_unique_at | JQ673480.1:53022 = TCACAGGTCGTCCTCGTCGGC | pass |
| spacer_ul27_ref | sequence_unique_at | NC_001806.2:53081 = TCACAGGTCGTCCTCGTCGGC | pass |
| spacer_ul27_lift | lift | lift(53022, JQ673480.1, NC_001806.2) = 53081 | pass |
| spacer_ul55_56_kos | sequence_unique_at | JQ673480.1:116044 = CCAGGCGTGGTGTGAGTTTG | pass |
| spacer_ul55_56_ref | sequence_unique_at | NC_001806.2:116124 = CCAGGCGTGGTGTGAGTTTG | pass |
| spacer_ul55_56_lift | lift | lift(116044, JQ673480.1, NC_001806.2) = 116124 | pass |
| spacer_ul56_kos | sequence_unique_at | JQ673480.1:116132 = ACAGGGGCGCTTACCGCCAC | pass |
| spacer_ul56_ref | sequence_unique_at | NC_001806.2:116212 = ACAGGGGCGCTTACCGCCAC | pass |
| spacer_ul56_lift | lift | lift(116132, JQ673480.1, NC_001806.2) = 116212 | pass |
| cds_ul3 | feature_location | ref_NC_001806.2.gb CDS/gene=UL3 -> 10991..11665 | pass |
| cds_ul27 | feature_location | ref_NC_001806.2.gb CDS/gene=UL27 -> complement(53081..55795) | pass |
| cds_ul56 | feature_location | ref_NC_001806.2.gb CDS/gene=UL56 -> complement(116222..116926) | pass |
| stop_ul3 | sequence_at | NC_001806.2:11663 = TAA | pass |
| stop_ul27 | sequence_at | NC_001806.2:53081 = TCA | pass |
| stop_ul56 | sequence_at | NC_001806.2:116222 = TTA | pass |

### tongmuang2024_ul56

| Check | Kind | Decided against | Result |
|---|---|---|---|
| lift1 | lift | lift(115941, MN159383.1, NC_001806.2) = 116222 | pass |
| lift2 | lift | lift(116645, MN159383.1, NC_001806.2) = 116926 | pass |
| lift3 | lift | lift(114423, MN159383.1, NC_001806.2) = 114704 | pass |
| arm_outer_declines | lift_declines | declined; detail 118446\|7928, wanted 118446\|7928 | pass |
| cds_ul56 | feature_location | ref_NC_001806.2.gb CDS/gene=UL56 -> complement(116222..116926) | pass |
| orf_translates | orf_on_minus_strand | MN159383.1:115941-116645 on the minus strand gives 234 residues, MASE...LWR* | pass |
| orf_matches_reference | orf_matches_translation | 234 residues against the reference's 234, 231 identical | pass |

### lassalle2020_excluded

| Check | Kind | Decided against | Result |
|---|---|---|---|
| cds_ul1 | feature_location | ref_NC_001806.2.gb CDS/gene=UL1 -> 9338..10012 | pass |
| cds_ul56 | feature_location | ref_NC_001806.2.gb CDS/gene=UL56 -> complement(116222..116926) | pass |
| cds_us1 | feature_location | ref_NC_001806.2.gb CDS/gene=US1 -> 132647..133909 | pass |
| cds_us12 | feature_location | ref_NC_001806.2.gb CDS/gene=US12 -> complement(145316..145582) | pass |
| cds_rs1 | feature_location | ref_NC_001806.2.gb CDS/gene=RS1 -> complement(127235..131131)\|147066..150962 | pass |
| boundaries_are_cds_edges | positions_are_feature_edges | 5 of 5 positions are edges of a CDS feature in ref_NC_001806.2.gb | pass |
| boundaries_are_not_repeat_edges | positions_are_feature_edges | 0 of 5 positions are edges of a repeat_region feature in ref_NC_001806.2.gb | pass |
| trs_counterpart | counterpart_pair | NC_001806.2:150962=A against 127235=T; complementary True; repeat map alternatives (127235,) | pass |
| irs_window | window_revcomp_equals | NC_001806.2:127235 TTACAGCACCCCGTCCCCCT... against NC_001806.2:150962 TTACAGCACCCCGTCCCCCT... | pass |

### jiao2019_ul53

| Check | Kind | Decided against | Result |
|---|---|---|---|
| lift_identity | lift | lift(112361, MN136524.1, NC_001806.2) = 112361 | pass |
| lift_identity_back | lift | lift(112361, NC_001806.2, MN136524.1) = 112361 | pass |
| cds_ul53_ref | feature_location | ref_NC_001806.2.gb CDS/gene=UL53 -> 112180..113196 | pass |
| cds_ul53_jn | feature_location | base_JN555585.1.gb CDS/gene=UL53 -> 112180..113196 | pass |
| codon_ref | sequence_at | NC_001806.2:112360 = CTA | pass |
| codon_mckrae | sequence_at | MN136524.1:112360 = CCA | pass |
| old_frame_base | sequence_at | X14112.1:112361 = A | pass |
| codon_api | codon_positions | codon(UL53, 61).positions = (112360, 112361, 112362) | pass |

### lee2018_dctrl2

| Check | Kind | Decided against | Result |
|---|---|---|---|
| declines_start | lift_declines | declined; detail 120302\|6072, wanted 120302\|6072 | pass |
| declines_end | lift_declines | declined; detail 120672\|5702, wanted 120672\|5702 | pass |
| copies_identical_start | window_revcomp_equals | NC_001806.2:120302 GTTAACCCCCCCCGCCCTTTGTGG... against NC_001806.2:6072 GTTAACCCCCCCCGCCCTTTGTGG... | pass |
| copies_identical_end | window_revcomp_equals | NC_001806.2:120672 AAAAGATACAAAAGTAAACCTTTA... against NC_001806.2:5702 AAAAGATACAAAAGTAAACCTTTA... | pass |
| published_start_in_the_named_deposit | sequence_unique_at | KT899744.1:120136 = GTTAACCCCCCCCCCGCCCTTTGTGGCGGCCAGAGGCCAG | pass |
| counterpart_start | counterpart_pair | NC_001806.2:120302=G against 6072=C; complementary True; repeat map alternatives (6072,) | pass |
| counterpart_end | counterpart_pair | NC_001806.2:120672=A against 5702=T; complementary True; repeat map alternatives (5702,) | pass |

### singh2022_dctrl2

| Check | Kind | Decided against | Result |
|---|---|---|---|
| frames_read_alike | window_equals | X14112.1:120503 CTCCCACGCACCCCCACTCCCACG... against NC_001806.2:120505 CTCCCACGCACCCCCACTCCCACG... | pass |
| convert_if_old | version_convert | convert(120500, NC_001806.1, NC_001806.2) = 120502 | pass |
| stated_length_is_wrong | interval_length | 120500 to 120635 inclusive is 136 bases | pass |

### chau2023_segment

| Check | Kind | Decided against | Result |
|---|---|---|---|
| version_matters_here | version_convert | convert(126983, NC_001806.1, NC_001806.2) = 126986 | pass |
| start_names_two_places | counterpart_pair | NC_001806.2:126983=C against 151214=G; complementary True; repeat map alternatives (151214,) | pass |
| end_names_two_places | counterpart_pair | NC_001806.2:132514=G against 145683=C; complementary True; repeat map alternatives (145683,) | pass |
| stated_length_disagrees | interval_length | 126983 to 132514 inclusive is 5532 bases | pass |

### colgrove2016_oris

| Check | Kind | Decided against | Result |
|---|---|---|---|
| start_names_two_places | counterpart_pair | NC_001806.2:131870=C against 146327=G; complementary True; repeat map alternatives (146327,) | pass |
| end_names_two_places | counterpart_pair | NC_001806.2:132250=A against 145947=T; complementary True; repeat map alternatives (145947,) | pass |

### kitabatake2026_ul50

| Check | Kind | Decided against | Result |
|---|---|---|---|
| lift1 | lift | lift(106888, GU734771.1, NC_001806.2) = 106992 | pass |
| lift2 | lift | lift(106918, GU734771.1, NC_001806.2) = 107022 | pass |
| lift3 | lift | lift(106947, GU734771.1, NC_001806.2) = 107051 | pass |
| lift4 | lift | lift(106908, GU734771.1, NC_001806.2) = 107012 | pass |
| lift5 | lift | lift(106967, GU734771.1, NC_001806.2) = 107071 | pass |
| revertant_primer_in_the_deposit | sequence_unique_at | GU734771.1:106888 = CATCGTCCCTAACAGGAAGATGAGTCAGTGGGGATCCGGGGCGATCCTTGTCCAGCCGGA | pass |
| revertant_primer_in_the_reference | sequence_unique_at | NC_001806.2:106992 = CATCGTCCCTAACAGGAAGATGAGTCAGTGGGGATCCGGGGCGATCCTTGTCCAGCCGGA | pass |
| second_window_in_the_deposit | sequence_unique_at | GU734771.1:106908 = TGAGTCAGTGGGGATCCGGGGCGATCCTTGTCCAGCCGGACAGCTTGGGTCGGGGGTACG | pass |
| second_window_in_the_reference | sequence_unique_at | NC_001806.2:107012 = TGAGTCAGTGGGGATCCGGGGCGATCCTTGTCCAGCCGGACAGCTTGGGTCGGGGGTACG | pass |
| cds_ul50 | feature_location | ref_NC_001806.2.gb CDS/gene=UL50 -> 107011..108126 | pass |
| codon4_reference | sequence_at | NC_001806.2:107020 = TGG | pass |
| codon4_deposit | sequence_at | GU734771.1:106916 = TGG | pass |
| codon4_api | codon_positions | codon(UL50, 4).positions = (107020, 107021, 107022) | pass |
| window_is_60_bases | interval_length | 106888 to 106947 inclusive is 60 bases | pass |

### dremel2025_circlat

| Check | Kind | Decided against | Result |
|---|---|---|---|
| declines_start | lift_declines | declined; detail 6151\|120223, wanted 6151\|120223 | pass |
| declines_end | lift_declines | declined; detail 6827\|119547, wanted 6827\|119547 | pass |
| copies_identical | window_revcomp_equals | NC_001806.2:6151 ACGTGGGTGTTGGGGGGTGGGGCG... against NC_001806.2:120223 ACGTGGGTGTTGGGGGGTGGGGCG... | pass |
| counterpart_start | counterpart_pair | NC_001806.2:6151=A against 120223=T; complementary True; repeat map alternatives (120223,) | pass |
| counterpart_end | counterpart_pair | NC_001806.2:6827=G against 119547=C; complementary True; repeat map alternatives (119547,) | pass |
| length_in_the_deposit | interval_length | 6137 to 6814 inclusive is 678 bases | pass |

### lee2023_us4_contig

| Check | Kind | Decided against | Result |
|---|---|---|---|
| cds_us4 | feature_location | ref_NC_001806.2.gb CDS/gene=US4 -> 136747..137463 | pass |
| us4_starts_with_met | sequence_at | NC_001806.2:136747 = ATG | pass |
| contig_is_193_bases | interval_length | 136823 to 137015 inclusive is 193 bases | pass |
