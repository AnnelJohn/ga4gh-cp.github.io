CSV to Family Phenopacket converter
=

This is a demo of how to use the Python `phenopackets` library to write a custom converter to extract `Family`
Phenopackets from CSV files.

Running
---
To run the demo, you will either need to have the minimum python version and dependencies, as listed in the
`pyproject.toml` installed on your system, use a virtual environment, or best of all use `uv`. `uv` can be installed
from https://docs.astral.sh/uv/. To run using `uv` from the project root:

```shell
uv run csv_to_phenopackets.py families.csv output
```

Validation
---

Validating the output file requires the
`phenopacket-tools-cli`(https://github.com/phenopackets/phenopacket-tools/releases/latest)
to be installed. This will require java 17+ to be installed on your system (`brew install java` on MacOS, or download
from https://adoptium.net/en-GB/temurin/releases or https://www.oracle.com/java/technologies/downloads/).

```shell
wget https://github.com/phenopackets/phenopacket-tools/releases/download/v1.0.0-RC3/phenopacket-tools-cli-1.0.0-RC3-distribution.zip
unzip phenopacket-tools-cli-1.0.0-RC3-distribution.zip
rm phenopacket-tools-cli-1.0.0-RC3-distribution.zip
```

To validate the output files, run:

```shell
java -jar phenopacket-tools-cli-1.0.0-RC3/phenopacket-tools-cli-1.0.0-RC3.jar validate output/*.json
```

