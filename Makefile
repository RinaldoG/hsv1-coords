.PHONY: inputs check test clean

inputs:        ## fetch and checksum every input into data/raw
	python3 -I src/fetch_inputs.py

check:         ## verify the inputs match the manifest
	python3 -I src/check_inputs.py

test:
	python3 -I -m unittest discover -s tests -v

clean:
	rm -rf data/raw/* data/derived/*
	touch data/raw/.gitkeep data/derived/.gitkeep
