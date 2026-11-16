.PHONY: all inputs genomes names aliases annotation crosswalk versions repeatmap gff3 agreement verify-crosswalk package-data testcases check test validate-gff3 clean

# Order matters twice over. `check` validates the inputs, so it runs before anything is
# built from them. `verify-crosswalk` checks the COMMITTED block files against a fresh
# alignment, so it runs BEFORE `crosswalk` overwrites them: round 12 measured a moved
# block failing the verifier standalone and passing `make all`, because the rebuild
# silently restored the file first. CI has always had this order; this target did not.
all: inputs genomes check verify-crosswalk names aliases annotation crosswalk versions \
     repeatmap gff3 agreement package-data testcases test

inputs:        ## fetch and checksum every input into data/raw
	python3 -I src/fetch_inputs.py

genomes:       ## fetch the complete genome set and build the strain inventory
	python3 -I src/fetch_genomes.py

names:         ## extract and reconcile gene names across the deposited records
	python3 -I src/extract_names.py

aliases:       ## rebuild the cited gene alias dictionary
	python3 -I src/build_aliases.py

annotation:    ## build the classed, sourced, versioned annotation table
	python3 -I src/build_annotation.py

crosswalk:     ## align every genome to the reference and build the coordinate blocks
	python3 -I src/build_crosswalk.py

versions:      ## map coordinates between NC_001806.1 and NC_001806.2
	python3 -I src/build_version_crosswalk.py

repeatmap:     ## build the repeat and diploid-locus map
	python3 -I src/build_repeat_map.py

gff3:          ## render the annotation as GFF3 (one file per source) and the repeat map as BED
	python3 -I src/build_gff3.py

agreement:     ## pin every genome's measured base agreement with the reference
	python3 -I src/build_base_agreement.py

verify-crosswalk:  ## check every COMMITTED block file against a fresh alignment
	python3 -I src/verify_crosswalk.py

package-data:  ## stage the built artefacts inside the package, for a wheel
	python3 -I src/stage_package_data.py

testcases:     ## regenerate the machine-readable synonymy test cases from the prose doc
	python3 -I src/build_testcases.py

check:         ## verify the inputs match the manifest and the assumptions still hold
	python3 -I src/check_inputs.py

test:          ## run the milestone gates
	python3 -I -m unittest discover -s tests -v

# Not in `all`: it needs two tools the library does not, genometools (`gt`) and the
# `gff3` validator from PyPI. See src/validate_gff3.py for how to install them without
# touching the system, and for what it accepts as a pass.
validate-gff3: ## check every GFF3 with two independent external parsers
	python3 -I src/validate_gff3.py

clean:
	# MANIFEST.tsv is tracked and is the only committed pin of the inputs: deleting it
	# meant `make inputs` re-pinned silently to whatever upstream served, with the CI
	# diff as the only guard. It is kept, and so is the directory marker.
	find data/raw -mindepth 1 ! -name .gitkeep ! -name MANIFEST.tsv -delete
	rm -f data/derived/names_by_source.tsv data/derived/name_reconciliation.tsv \
	      data/derived/mechanical_report.md
