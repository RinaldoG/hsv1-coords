.PHONY: all inputs names testcases check test clean

all: inputs names testcases check test

inputs:        ## fetch and checksum every input into data/raw
	python3 -I src/fetch_inputs.py

names:         ## extract and reconcile gene names across the deposited records
	python3 -I src/extract_names.py

testcases:     ## regenerate the machine-readable synonymy test cases from the prose doc
	python3 -I src/build_testcases.py

check:         ## verify the inputs match the manifest and the assumptions still hold
	python3 -I src/check_inputs.py

test:          ## run the milestone gates
	python3 -I -m unittest discover -s tests -v

clean:
	rm -rf data/raw/* data/derived/names_by_source.tsv data/derived/name_reconciliation.tsv \
	       data/derived/mechanical_report.md
	touch data/raw/.gitkeep
