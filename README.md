# PolymerSolv-Studio
Educational Python application for evidence-based polymer–solvent selection using HSP/RED screening and a multi-criteria evidence framework.
The application combines Hansen solubility parameter (HSP) screening
with literature evidence, processing conditions, environmental,
health and safety requirements, and evidence completeness.

## Important limitation

HSP distance and relative energy difference (RED) are screening
indicators. They do not prove dissolution, processing suitability,
safety, environmental acceptability or regulatory compliance.

The software does not automatically select a “best” solvent. Users
must evaluate the available evidence and justify the final decision.

## Teaching cases

The application includes three teaching cases:

- Nafion processing
- Cellulose acetate membrane fabrication
- Polystyrene and expanded-polystyrene dissolution and recovery

Users may also formulate their own polymer–solvent selection problem.

## Learning workflow

PolymerSolv Studio follows six stages:

1. **DEFINE:** State the application, alternatives, criteria and constraints.
2. **INVESTIGATE:** Record relevant literature and experimental conditions.
3. **MODEL:** Calculate and interpret HSP distance and RED.
4. **COMPARE:** Compare model results with literature evidence.
5. **DECIDE:** Justify a recommendation using the available evidence.
6. **REFLECT:** Identify uncertainty and propose experimental verification.

## System requirements

- Python 3.x
- A graphical desktop environment
- Required packages listed in `requirements.txt`

The application is intended for desktop use. Python installations on
Windows normally include Tkinter. Some Linux distributions require
Tkinter to be installed separately.

## Downloading the software

Select the green **Code** button on the GitHub repository page and
choose **Download ZIP**.

Extract the ZIP file before running the application.

## Installation

Open a terminal or command prompt inside the extracted folder.

Create a virtual environment:

### Windows

```bash
python -m venv .venv
.venv\Scripts\activate
```

### Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the required packages:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

If Tkinter is unavailable on Ubuntu or Debian Linux, it can usually be
installed using:

```bash
sudo apt install python3-tk
```

## Running the application

Run:

```bash
python PolymerSolv-Studio-GUI-Application-Python.py
```

On systems where Python 3 is called `python3`, use:

```bash
python3 PolymerSolv-Studio-GUI-Application-Python.py
```

## Outputs

Depending on the selected workflow and export options, the application
can produce:

- Material-selection briefs
- Evidence-audit files
- Literature-evidence CSV files
- PDF reports
- HSP projection figures
- Domain-specific RED figures
- Binary-mixture sensitivity figures

Representative outputs may be added to the `example-output` folder.

## Evidence responsibility

Literature records are entered by users and are not retrieved or
verified automatically. Users must verify:

- Bibliographic information
- Polymer identity and grade
- Solvent identity and composition
- Concentration, temperature and processing conditions
- Safety and environmental information
- Regulatory and experimental claims

Missing information must not be treated as evidence of suitability.

## Educational status

PolymerSolv Studio is an educational decision-support application. It
is not a replacement for laboratory testing, professional engineering
judgement, safety assessment or regulatory review.

The proposed learning outcomes and assessment rubric require classroom
evaluation before claims about educational effectiveness can be made.

## Citation

If this software supports published work, teaching or course material,
please cite the associated journal article.

The complete citation and DOI will be added after publication.

## Software version

This repository contains PolymerSolv Studio version 1.0.0.

## Licence

The source code is distributed under the MIT License. See the
`LICENSE` file for details.

## Contact

For questions, corrections or suggestions, use the **Issues** section
of this GitHub repository.
