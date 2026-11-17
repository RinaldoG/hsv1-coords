"""The gate on M3's by-hand criterion: published coordinates resolved from the papers.

The criterion is quoted in `build_published_coordinates.CRITERION` and in
`ROADMAP.md`. It is the one M3 criterion a test run cannot settle on its own, because
it needs papers read, and `PITFALLS.md` root cause 4 is this project marking exactly
that kind of criterion met twice without checking. So the record of the reading is two
committed tables, and this is what holds them to account.

What this file does and does not claim.

- It does **not** assert that any paper says what a row claims. The published
  coordinate, the quote and the fetch route are hand transcription, flagged as such in
  the `transcribed` column, and the suite has no network. Asserting otherwise would be
  a gate that cannot fail dressed up as provenance.
- It **does** assert everything downstream: that each resolution the table records is
  the resolution the library produces, that the sequence evidence is where the table
  says it is in the records themselves, that a counted case carries both a resolution
  and independent evidence, that every value a flag can take carries its own burden,
  and that the criterion's two thresholds are met by counts derived from the table
  rather than by numbers typed beside them.
- It **does** assert the triangle that makes a case evidence rather than assertion: a
  sequence the paper prints sits at the published coordinate in the record the paper
  names, the same sequence sits at the resolved coordinate in the reference, and the
  library maps the first to the second. All three must hold together.

`TheGateHasPower` is the part worth reading first. Standing rule 28: a gate is accepted
only after a mutation has made it fail, and the question is not whether it asserts
something but what reachable state falsifies it. The sweep there perturbs **every** one
of the check rows in turn and requires each perturbation to be caught, and reports the
residue as a number, because one mutation target measures one mutation target
(rule 57).
"""
from __future__ import annotations

import copy
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import build_published_coordinates as pc  # noqa: E402

RAW = ROOT / "data" / "raw"
NEEDED = ("genomes.fasta", "ref_NC_001806.2.gb", "ref_NC_001806.1.gb",
          "base_JN555585.1.gb")


def _artefacts_present() -> list[str]:
    return [n for n in NEEDED if not (RAW / n).exists()]


# --------------------------------------------------------------- validation helpers
#
# Each returns a list of problems rather than asserting, so that `TheGateHasPower`
# can feed it a deliberately broken table and require a problem back. A helper that
# only ever ran against the real table would be a check whose power nothing measured.

def correspondence_problems(cases: list[dict], checks: list[dict]) -> list[str]:
    """Every check names a case, every case that should have checks has them.

    Both directions. Rule 62: a row a gate cannot match is a finding, not a row to
    skip, and a one-way membership test is blind to half of what can go wrong.
    """
    out = []
    names = [c["case"] for c in cases]
    for check in checks:
        if check["case"] not in names:
            out.append(f"check {check['case']}/{check['check']} names no case")
    seen = set()
    for check in checks:
        key = (check["case"], check["check"])
        if key in seen:
            out.append(f"duplicate check key {key}")
        seen.add(key)
    for case in cases:
        mine = [k for k in checks if k["case"] == case["case"]]
        if case["status"] in ("counted", "supporting") and not mine:
            out.append(f"case {case['case']} is {case['status']} with no checks, so "
                       f"it is prose")
        if case["status"] == "not_counted" and case["resolution_kind"] != "none" \
                and not mine:
            out.append(f"case {case['case']} records a resolution_kind of "
                       f"{case['resolution_kind']} and carries no check")
    if len(names) != len(set(names)):
        out.append("the cases table repeats a case key")
    return out


def vocabulary_problems(cases: list[dict], checks: list[dict]) -> list[str]:
    out = []
    for case in cases:
        for column, allowed in (("status", pc.STATUSES), ("frame_basis", pc.BASES),
                                ("resolution_kind", pc.RESOLUTION_KINDS),
                                ("transcribed", pc.TRANSCRIBED)):
            if case[column] not in allowed:
                out.append(f"{case['case']}.{column} is {case[column]!r}, which is "
                           f"not one of {allowed}")
        for column in ("case", "strain", "frame_basis", "status", "status_reason",
                       "citation", "pmcid", "source", "transcribed"):
            if not case[column]:
                out.append(f"{case['case']}.{column} is empty")
    for check in checks:
        if check["kind"] not in pc.CHECK_KINDS:
            out.append(f"{check['case']}/{check['check']} has kind "
                       f"{check['kind']!r}, which the runner does not decide")
    return out


def burden_problems(cases: list[dict], checks: list[dict],
                    genomes: dict[str, str]) -> list[str]:
    """What each value of `status` and `frame_basis` owes.

    Standing rule 49: every value a flag can take carries its own burden, and the
    burden is what the word means. A two-valued field with an obligation on only one
    value is answered by writing the other, and that is how four earlier rounds of
    hardening were walked past in the conventions table.
    """
    out = []
    for case in cases:
        mine = [k for k in checks if k["case"] == case["case"]]
        kinds = {k["kind"] for k in mine}
        basis, status = case["frame_basis"], case["status"]

        if status == "counted":
            if basis in pc.BASES_THAT_CANNOT_COUNT:
                out.append(f"{case['case']} is counted on a frame_basis of {basis}, "
                           f"which cannot support a counted case: a paper that does "
                           f"not state which record its numbers are in has not "
                           f"published a resolvable coordinate")
            if case["resolution_kind"] not in pc.KINDS_THAT_RESOLVE:
                out.append(f"{case['case']} is counted with resolution_kind "
                           f"{case['resolution_kind']}, which is not a resolution: "
                           f"the criterion is about resolving a coordinate, and an "
                           f"identity, a decline and an unresolvable case each leave "
                           f"the crosswalk untested")
            if not kinds & set(pc.RESOLVING_CHECKS):
                out.append(f"{case['case']} is counted with no check tying its "
                           f"published coordinate to a resolved one")
            if not kinds & set(pc.EVIDENCE_CHECKS):
                out.append(f"{case['case']} is counted with no evidence from outside "
                           f"the library, so its resolution rests on the library "
                           f"agreeing with itself")

        if basis == "stated_with_version":
            if "." not in case["frame_accession"].split(",")[0]:
                out.append(f"{case['case']} claims stated_with_version and its "
                           f"frame_accession {case['frame_accession']!r} carries no "
                           f"version suffix")
        elif basis == "stated_accession_held":
            accession = case["frame_accession"]
            if "." in accession:
                out.append(f"{case['case']} claims stated_accession_held, which is "
                           f"for an accession written WITHOUT a version suffix, and "
                           f"{accession!r} carries one")
            held = sorted(k for k in genomes if k.rsplit(".", 1)[0] == accession)
            if len(held) != 1:
                out.append(f"{case['case']} claims stated_accession_held and the "
                           f"corpus holds {len(held)} versions of {accession}: "
                           f"{held}. The claim is that there is exactly one, so "
                           f"the absent version suffix cannot matter")
            if not kinds & {"sequence_unique_at", "feature_location"}:
                out.append(f"{case['case']} claims stated_accession_held with nothing "
                           f"confirming the paper's coordinate against that record")
        elif basis == "determined_from_paper_claim":
            # The burden is the strongest of the three, because the frame is inferred
            # rather than read: the inference has to be shown, in the old frame, more
            # than once, or it is one coincidence.
            old_frame = [k for k in mine
                         if k["kind"] == "sequence_unique_at"
                         and k["record"] in ("X14112.1", "NC_001806.1")]
            if len(old_frame) < 3:
                out.append(f"{case['case']} infers its frame from the paper's own "
                           f"claim and shows {len(old_frame)} unique sequence "
                           f"placements in the old frame; three is the floor, because "
                           f"one or two could be coincidence")
            if not case["frame_evidence"]:
                out.append(f"{case['case']} infers its frame and states no evidence")
        elif basis in pc.BASES_THAT_CANNOT_COUNT:
            if not case["frame_evidence"] and basis != "none_printed":
                out.append(f"{case['case']} has basis {basis} and says nothing about "
                           f"what was tried")

        if case["transcribed"] == "by_hand_fetch_failed" and status == "counted":
            out.append(f"{case['case']} is counted on a paper that was not retrieved")
        if not case["status_reason"]:
            out.append(f"{case['case']} states no reason for its status of {status}")
    return out


def criterion_problems(cases: list[dict]) -> list[str]:
    """The criterion's own two thresholds, against counts derived from the table.

    Rule 47: a count is compared against the thing it counts, never against a typed
    number. The two thresholds ARE typed, because they are the criterion's words and
    not a measurement; what must not be typed is how many papers and strains there are.
    """
    out = []
    kept = pc.counted(cases)
    papers = len({c["pmcid"] for c in kept})
    if papers != len(kept):
        out.append(f"{len(kept)} counted rows cover only {papers} distinct papers, so "
                   f"one paper is being counted twice")
    if papers < pc.PAPERS_REQUIRED:
        out.append(f"the criterion asks for at least {pc.PAPERS_REQUIRED} papers and "
                   f"{papers} are counted")
    strains = pc.strains_counted(cases)
    if len(strains) < pc.STRAINS_REQUIRED:
        out.append(f"the criterion asks for at least {pc.STRAINS_REQUIRED} strains and "
                   f"the counted cases cover {len(strains)}: {strains}")
    return out


def triangle_problems(cases: list[dict], checks: list[dict]) -> list[str]:
    """For every motif placed in two records, the library must map one onto the other.

    This is the assertion that makes a counted case evidence. Neither half alone says
    anything: a motif can sit at the published coordinate in one record and at some
    other coordinate in the reference, and the lift can be internally consistent while
    pointing somewhere else. Only the two together pin the resolution.

    **The first version of this skipped any motif placed in fewer than two records**,
    which is root cause 30 again: a containment check cannot see an omission. Deleting
    the reference-side placement of a motif left the suite green, because the motif
    then had one record and was skipped, and the case still closed its other triangle.
    So the rule is stated the other way round now. Every motif placed in a record that
    is NOT the reference owes a reference placement, and the number of triangles is
    compared with the number of rows that owe one. A reference-only motif owes nothing,
    which is right: the counterpart copy of a repeat has no deposit-side partner.
    """
    out = []
    from hsv1coords import lift
    from hsv1coords.versions import VersionCrosswalk
    version = VersionCrosswalk.load()
    for case in pc.counted(cases):
        mine = [k for k in checks if k["case"] == case["case"]]
        placed: dict[str, list[dict]] = {}
        for check in mine:
            if check["kind"] == "sequence_unique_at":
                placed.setdefault(check["motif"], []).append(check)
        triangles = 0
        owed = sum(1 for rows in placed.values() for r in rows
                   if r["record"] != "NC_001806.2")
        for motif, rows in placed.items():
            reference = [r for r in rows if r["record"] == "NC_001806.2"]
            others = [r for r in rows if r["record"] != "NC_001806.2"]
            if not others:
                # A reference-only motif, such as the other copy of a repeat. It owes
                # no triangle, and it is not counted in `owed` either.
                continue
            if not reference:
                out.append(f"{case['case']}: motif {motif} is placed in "
                           f"{sorted(r['record'] for r in others)} with no placement "
                           f"in NC_001806.2, so nothing closes its triangle")
                continue
            want = int(reference[0]["position"])
            for row in others:
                source, position = row["record"], int(row["position"])
                if source in ("X14112.1", "NC_001806.1"):
                    got = version.convert(position, "NC_001806.1", "NC_001806.2")
                    how = "VersionCrosswalk.convert"
                else:
                    got = lift(position, source, "NC_001806.2")
                    how = "lift"
                if got != want:
                    out.append(f"{case['case']}: {motif} sits at {source}:{position} "
                               f"and at NC_001806.2:{want}, but {how} gives {got}")
                triangles += 1
        if triangles != owed:
            out.append(f"{case['case']}: {owed} sequence placements are in a record "
                       f"other than the reference and {triangles} triangles closed, "
                       f"so a placement has lost its reference partner")
        if not triangles and case["frame_basis"] != "stated_with_version":
            # Pinnoji prints no usable sequence (its figure is an image), so its
            # triangle is closed by the two records' own TATA_box features instead,
            # which `burden_problems` requires. Named rather than exempted silently:
            # rule 46, an exemption inside a gate is named and counted.
            if case["case"] not in TRIANGLE_EXEMPT:
                out.append(f"{case['case']} closes no sequence triangle and is not "
                           f"in TRIANGLE_EXEMPT")
    return out


TRIANGLE_EXEMPT = {
    "pinnoji2007_nrse": "The paper's own sequence figure is an image, so no printed "
                        "motif is available. Its frame and its resolution are closed "
                        "instead by the US1 TATA_box feature of both reference "
                        "records, which sits at the last six bases of the published "
                        "interval in each frame, three apart, which is the offset the "
                        "version crosswalk applies there.",
}
"""Counted cases that close no sequence triangle, each with its reason.

Rule 59: an exemption is a named, counted list, preferably empty. This one has a
single entry and `test_the_exemption_list_is_not_padded` holds it to that.
"""


# -------------------------------------------------------------------------- the tests

class Tables(unittest.TestCase):
    """The two tables are committed inputs, so an absent one is a failure, not a skip."""

    def test_both_tables_are_present_and_not_empty(self):
        for path in (pc.CASES, pc.CHECKS):
            self.assertTrue(path.exists(), f"{path} is missing; it is a committed "
                                           f"input, not a built artefact")
        self.assertGreater(len(pc.load_cases()), 0)
        self.assertGreater(len(pc.load_checks()), 0)

    def test_the_tables_correspond_in_both_directions(self):
        self.assertEqual(
            correspondence_problems(pc.load_cases(), pc.load_checks()), [])

    def test_every_controlled_value_is_in_its_vocabulary(self):
        self.assertEqual(vocabulary_problems(pc.load_cases(), pc.load_checks()), [])


class ArtefactBackedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        absent = _artefacts_present()
        if absent:
            raise unittest.SkipTest(
                f"{absent} are fetched, not committed. Run `make inputs genomes`. "
                f"This skip is registered in test_m5_conventions.py, where "
                f"test_a_checkout_holds_every_artefact_the_gates_skip_without turns "
                f"it into a failure in a checkout.")
        cls.art = pc.Artefacts()
        cls.cases = pc.load_cases()
        cls.checks = pc.load_checks()
        cls.results = pc.run_all(cls.art)


class EveryCheckHolds(ArtefactBackedTests):
    def test_every_check_holds(self):
        failed = [(row["case"], row["check"], why)
                  for row, ok, why in self.results if not ok]
        self.assertEqual(failed, [],
                         f"{len(failed)} of {len(self.results)} checks fail")

    def test_the_checks_cover_every_kind_the_runner_decides(self):
        """A kind with no row is a branch of the runner nothing exercises."""
        used = {row["kind"] for row in self.checks}
        self.assertEqual(sorted(set(pc.CHECK_KINDS) - used), [],
                         "the runner decides a kind no row uses, so that branch is "
                         "untested; delete the branch or add the row")

    def test_the_count_of_checks_is_reported_against_what_it_counts(self):
        """Rule 47 and rule 16: derived, printed, and compared with the table."""
        per_case = {}
        for row in self.checks:
            per_case[row["case"]] = per_case.get(row["case"], 0) + 1
        print(f"\n{len(self.checks)} checks over {len(per_case)} of "
              f"{len(self.cases)} cases, "
              f"{sum(1 for _, ok, _ in self.results if ok)} passing; per case "
              f"{dict(sorted(per_case.items()))}")
        self.assertEqual(sum(per_case.values()), len(self.checks))
        for case in pc.counted(self.cases):
            self.assertGreaterEqual(
                per_case.get(case["case"], 0), 5,
                f"{case['case']} is counted on fewer than five checks")


class EveryStatusCarriesItsBurden(ArtefactBackedTests):
    def test_every_burden_is_discharged(self):
        self.assertEqual(
            burden_problems(self.cases, self.checks, self.art.genomes), [])

    def test_the_exemption_list_is_not_padded(self):
        named = set(TRIANGLE_EXEMPT)
        cases = {c["case"] for c in pc.counted(self.cases)}
        self.assertEqual(named - cases, set(),
                         "TRIANGLE_EXEMPT names a case that is not counted")
        self.assertLessEqual(len(named), 1,
                             f"TRIANGLE_EXEMPT has grown to {len(named)} entries; an "
                             f"exemption list that grows is the gate being walked past")
        for case, reason in TRIANGLE_EXEMPT.items():
            self.assertGreater(len(reason), 80, f"{case}'s exemption states no reason")


class TheTrianglesClose(ArtefactBackedTests):

    TRIANGLES_PER_CASE = {
        "amelio2006_ctrl2": 6,
        "creech2010_lat5exon": 1,
        "pinnoji2007_nrse": 1,
        "singh2024_kos_arms": 4,
        "tongmuang2024_ul56": 0,
        "kitabatake2026_ul50": 2,
    }
    """How many closed triangles each counted case carries, pinned per case.

    Rule 60 and rule 47: a per-item measurement is pinned per item, and a count is
    compared against the thing it counts rather than against a floor. A floor would
    not have caught what this did.

    Why it exists. `triangle_problems` is written as an obligation on the non-reference
    side: a motif placed in a deposit owes a placement in the reference. That catches a
    deleted reference-side row. It cannot catch a deleted DEPOSIT-side row, because the
    motif then becomes reference-only, owes nothing, and the count of what is owed
    falls with it, so both sides of that comparison move together and it stays green.
    That is root cause 41's shape inside my own fix for root cause 30, found by
    deleting the row and watching the suite pass.

    A number here that falls is a deleted placement and needs a reason in the commit
    message. A number that rises is new evidence, which is ordinary, and raising it is
    the commit that earned it.
    """

    def test_every_printed_motif_pins_its_own_resolution(self):
        self.assertEqual(triangle_problems(self.cases, self.checks), [])

    def test_no_case_carries_fewer_triangles_than_it_did(self):
        counted = {}
        for case in pc.counted(self.cases):
            placed: dict[str, list[dict]] = {}
            for check in self.checks:
                if (check["case"] == case["case"]
                        and check["kind"] == "sequence_unique_at"):
                    placed.setdefault(check["motif"], []).append(check)
            counted[case["case"]] = sum(
                1 for rows in placed.values()
                if any(r["record"] != "NC_001806.2" for r in rows)
                and any(r["record"] == "NC_001806.2" for r in rows))
        self.assertEqual(
            counted, self.TRIANGLES_PER_CASE,
            "the number of closed sequence triangles has changed. A fall is a deleted "
            "placement: the evidence tying a paper's printed sequence to its own "
            "coordinate in the record it names, or tying that to the reference. "
            "Neither half is optional.")


class TheCriterionIsMet(ArtefactBackedTests):
    """M3's fourth criterion, named so that this test is what marks it met."""

    def test_at_least_five_papers_across_at_least_three_strains(self):
        self.assertEqual(criterion_problems(self.cases), [],
                         f"the criterion is: {pc.CRITERION}")

    def test_the_criterion_text_is_quoted_in_the_roadmap(self):
        """Rule 11: a correction is not applied until something executable changes.

        If the roadmap's wording of the criterion and this module's drift apart, the
        thing being marked met stops being the thing that was written down.
        """
        roadmap = (ROOT / "ROADMAP.md").read_text()
        flat = " ".join(roadmap.split())
        self.assertIn(" ".join(pc.CRITERION.split()), flat,
                      "ROADMAP.md no longer quotes the criterion this module gates")

    def test_the_criterion_survives_dropping_every_inferred_frame(self):
        """The margin, asserted rather than described.

        One counted case has a reference version this repository inferred from the
        paper's own printed sequence rather than read from the paper's words. When
        that was the only route to the fifth paper, the criterion rested on an
        inference, which is thin for something a specialist reads. It no longer does,
        and this is what holds that open: drop every `determined_from_paper_claim`
        case and the criterion must still be met.

        If this test is what fails after a future edit, the honest repair is to find
        another paper, not to relax the test. Deleting it would put the dependence
        back without saying so.
        """
        cases = [c for c in self.cases
                 if c["frame_basis"] != "determined_from_paper_claim"]
        inferred = len(self.cases) - len(cases)
        self.assertGreater(inferred, 0,
                           "no case is marked determined_from_paper_claim, so this "
                           "test is vacuous; it must be deleted or the label restored")
        self.assertEqual(
            criterion_problems(cases), [],
            f"the criterion is met only because of the {inferred} case(s) whose "
            f"coordinate frame was inferred rather than stated. That was true before "
            f"17 November 2026 and was recorded as the margin; it must not become "
            f"true again silently.")

    def test_the_strain_count_is_over_normalised_strains(self):
        """Three papers write 17, 17syn+ and 17Syn+ for one strain.

        Counting those as three would clear the strain threshold on spelling. The
        `strain` column is normalised and `strain_as_stated` keeps each paper's words,
        so this asserts the two columns are not the same thing.
        """
        kept = pc.counted(self.cases)
        stated = {c["strain_as_stated"] for c in kept if c["strain_as_stated"]}
        self.assertGreater(len(stated), len(pc.strains_counted(self.cases)),
                           "every counted case states its strain exactly as the "
                           "normalised column does, so this gate cannot see the "
                           "spelling problem it exists for")


class TheGateHasPower(ArtefactBackedTests):
    """Four controls. Each must turn a named assertion above red.

    Standing rule 28: "does it assert something" is the wrong question, "what
    reachable state falsifies it" is the right one. Rule 81: when a fix is written, run
    its test without the fix.
    """

    @staticmethod
    def _perturb(row: dict) -> dict | None:
        """A meaningful single-field change for one check row, or None if there is
        no obvious one for its kind."""
        bad = copy.deepcopy(row)
        kind = row["kind"]
        if kind in ("sequence_at", "sequence_unique_at"):
            motif = row["motif"]
            flipped = ("A" if motif[0] != "A" else "C") + motif[1:]
            bad["motif"] = flipped
        elif kind == "feature_location":
            bad["expect"] = row["expect"] + "|1..2"
        elif kind == "lift_declines":
            primary, _, alternative = row["expect"].partition("|")
            bad["expect"] = f"{int(primary) + 1}|{alternative}"
        elif kind in ("window_equals", "window_revcomp_equals"):
            record, _, position = row["record"].partition(":")
            bad["record"] = f"{record}:{int(position) + 7}"
        elif kind == "codon_positions":
            bad["expect"] = ",".join(str(int(p) + 1)
                                     for p in row["expect"].split(","))
        elif kind in ("orf_on_minus_strand", "orf_matches_translation"):
            bad["position"] = str(int(row["position"]) + 1)
        elif row["expect"].isdigit():
            bad["expect"] = str(int(row["expect"]) + 1)
        elif kind == "positions_are_feature_edges":
            bad["expect"] = str(int(row["expect"]) + 1)
        else:
            return None
        return bad

    def test_every_check_row_can_be_made_to_fail(self):
        """The residue, as a number. Rule 57: one mutation target measures one
        mutation target, so this sweeps all of them and reports what survived."""
        survived, unperturbable = [], []
        for row in self.checks:
            bad = self._perturb(row)
            if bad is None:
                unperturbable.append(f"{row['case']}/{row['check']} ({row['kind']})")
                continue
            ok, _ = pc.run_check(bad, self.art)
            if ok:
                survived.append(f"{row['case']}/{row['check']} ({row['kind']})")
        print(f"\nmutation sweep: {len(self.checks) - len(survived) - len(unperturbable)}"
              f" of {len(self.checks)} check rows caught, {len(survived)} survived, "
              f"{len(unperturbable)} had no perturbation")
        self.assertEqual(survived, [], "a perturbed check row still passes, so that "
                                       "row asserts less than it appears to")
        self.assertEqual(unperturbable, [], "a check row has no perturbation, so "
                                            "nothing here measures whether it can "
                                            "fail")

    def test_flipping_a_status_to_counted_is_caught(self):
        """A case whose frame the paper does not state cannot be promoted quietly."""
        cases = copy.deepcopy(self.cases)
        victim = next(c for c in cases if c["frame_basis"] == "not_stated")
        victim["status"] = "counted"
        problems = burden_problems(cases, self.checks, self.art.genomes)
        self.assertNotEqual(problems, [],
                            f"{victim['case']} was relabelled counted on a "
                            f"frame_basis of not_stated and nothing complained")

    def test_dropping_a_counted_case_s_resolving_checks_is_caught(self):
        cases = copy.deepcopy(self.cases)
        victim = pc.counted(cases)[0]["case"]
        thinned = [k for k in copy.deepcopy(self.checks)
                   if not (k["case"] == victim
                           and k["kind"] in pc.RESOLVING_CHECKS)]
        self.assertNotEqual(
            burden_problems(cases, thinned, self.art.genomes), [],
            f"every check tying {victim}'s published coordinate to a resolved one was "
            f"deleted and nothing complained")

    def test_relabelling_a_strain_cannot_manufacture_the_strain_count(self):
        """The three strain-17 papers are where this would be easiest."""
        cases = copy.deepcopy(self.cases)
        kept = [c for c in cases if c["status"] == "counted"]
        for c in kept:
            c["strain"] = "17"
        self.assertNotEqual(
            criterion_problems(cases), [],
            "every counted case was relabelled as one strain and the criterion still "
            "reported as met, so the strain threshold is not being checked")


if __name__ == "__main__":
    unittest.main()
