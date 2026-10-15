"""
PolymerSolv Studio — Evidence-Based Polymer–Solvent Material Selection
Modified teaching version following the second-review action list.

Core principle
--------------
HSP/RED output is evidence for an engineering decision; it is not the decision.

The program implements a guided six-stage workflow:
DEFINE → INVESTIGATE → MODEL → COMPARE → DECIDE → REFLECT

Dependencies
------------
Required: numpy, matplotlib
Optional for GUI: tkinter (normally included with standard Windows Python)

Run:
    python PolymerSolv_Studio_Modified.py
    python PolymerSolv_Studio_Modified.py --test
    python PolymerSolv_Studio_Modified.py --export-all

Important data note
-------------------
The supplied case data are teaching records, not a certified engineering database.
Each record identifies whether an item is adopted, source-derived, qualitative, or
expert-derived. Instructors must replace provisional records with source-verified
values appropriate to the polymer grade, temperature, process, and jurisdiction.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.patches import Ellipse

try:
    import tkinter as tk
    from tkinter import messagebox, ttk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    TK_AVAILABLE = True
except Exception:
    tk = None
    ttk = None
    messagebox = None
    FigureCanvasTkAgg = None
    TK_AVAILABLE = False

APP_NAME = "PolymerSolv Studio"
APP_VERSION = "2.0 — evidence-audited teaching version"
HSP_UNITS = "MPa½"
DEFAULT_WEIGHTS = [30.0, 25.0, 25.0, 10.0, 10.0]
WEIGHT_LABELS = [
    "Thermodynamic screening",
    "Process evidence",
    "EHS / regulatory evidence",
    "Sustainability indicator",
    "Cost / availability evidence",
]

DARK = {
    "figure": "#0B0F19",
    "axes": "#151E2E",
    "text": "#F8FAFC",
    "subtext": "#94A3B8",
    "grid": "#243247",
    "border": "#334155",
    "blue": "#38BDF8",
    "green": "#10B981",
    "amber": "#F59E0B",
    "red": "#EF4444",
    "purple": "#A78BFA",
}


# -----------------------------------------------------------------------------
# Evidence and case database
# -----------------------------------------------------------------------------

def source_record(citation: str, value_type: str, temperature: str = "Not reported",
                  method: str = "Adopted teaching record", uncertainty: str = "Not reported") -> dict:
    """Create transparent provenance metadata for a database entry."""
    return {
        "citation": citation,
        "value_type": value_type,
        "temperature": temperature,
        "method": method,
        "uncertainty": uncertainty,
    }


def candidate(hsp: Tuple[float, float, float], bp_c: Optional[float],
              process_evidence: str, ehs_evidence: str,
              sustainability_evidence: str, cost_evidence: str,
              literature_evidence: str, hard_gate: dict,
              source: dict, color: str) -> dict:
    """Build a candidate record whose evidence is inspectable by students."""
    return {
        "hsp": np.array(hsp, dtype=float),
        "bp_c": bp_c,
        "process_evidence": process_evidence,
        "ehs_evidence": ehs_evidence,
        "sustainability_evidence": sustainability_evidence,
        "cost_evidence": cost_evidence,
        "literature_evidence": literature_evidence,
        "hard_gate": hard_gate,
        "source": source,
        "color": color,
    }


# Scores are deliberately marked as qualitative instructional indicators.
# They are not measured physical properties and can be inspected in the audit tab.
CASES_DATABASE: Dict[str, dict] = {
    "nafion": {
        "id": "nafion",
        "name": "Case 1: Nafion Processing",
        "discipline": "Energy and electrochemical engineering",
        "problem_statement": (
            "Evaluate solvent and solvent-mixture choices for Nafion processing. "
            "Use HSP results together with reported performance, process evidence, "
            "and EHS information. Investigate any difference between model predictions "
            "and reported behavior."
        ),
        "application": "Preparation of Nafion-based coatings or membranes for electrochemical devices.",
        "student_guidance": "High: use the supplied evidence, then explain the model–literature discrepancy.",
        "default_blend": ("Isopropanol", "Water", 80.0),
        "criteria_prompt": (
            "Define what ‘suitable’ means for the intended process. Consider compatibility, "
            "drying or removal, safety controls, reported morphology, and evidence quality."
        ),
        "polymer": {
            "name": "Nafion 117 perfluorosulfonic-acid ionomer",
            "type": "multi-domain teaching representation",
            "material_specification": "Nafion 117; equivalent weight and pretreatment must be recorded before use in a real study.",
            "source": source_record(
                "Kusoglu and Weber (2017), Chemical Reviews 117, 987–1104; adopted HSP teaching parameters.",
                "Adopted literature-derived teaching values",
                "Not uniformly reported", "Domain representation for instructional screening", "Parameter applicability must be checked"
            ),
            "domains": {
                "Backbone": {
                    "hsp": np.array([16.4, 10.5, 8.9]),
                    "R0": 8.8,
                    "description": "Hydrophobic fluorocarbon-rich region",
                    "color": "#38BDF8",
                    "source": source_record("Adopted from the prior validated teaching dataset; verify against the primary parameter source.", "Adopted", "Not reported", "HSP representation", "Not reported"),
                },
                "Side-chain": {
                    "hsp": np.array([15.2, 11.7, 15.9]),
                    "R0": 7.7,
                    "description": "Hydrophilic sulfonic-acid-containing region",
                    "color": "#F43F5E",
                    "source": source_record("Adopted from the prior validated teaching dataset; verify against the primary parameter source.", "Adopted", "Not reported", "HSP representation", "Not reported"),
                },
            },
        },
        "solvents": {
            "DMAc": candidate(
                (16.8, 11.5, 10.2), 165.0,
                "Boiling point is high relative to volatile alcohols; solvent removal may require time or vacuum.",
                "Hazard classification and local handling controls must be checked against the current SDS and local regulations.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Reported as a useful high-temperature processing solvent in polymer-membrane practice; conditions and grade matter.",
                {"status": "WARN", "rule": "No universal hard gate assigned", "reason": "Use only with course-defined controls and a current SDS review.", "source": "Instructor must specify local EHS policy and SDS date."},
                source_record("Hansen (2007) adopted solvent coordinates; physical-property source must be verified before deployment.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#34D399",
            ),
            "DMF": candidate(
                (17.4, 13.7, 11.3), 153.0,
                "Moderate boiling point; removal and residual-solvent controls should be evaluated.",
                "Hazard classification and legal restrictions must be checked for the institution and jurisdiction.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Historically reported for polymer processing; suitability is process-specific.",
                {"status": "WARN", "rule": "No universal hard gate assigned", "reason": "Current local regulatory and SDS review required.", "source": "Instructor must specify local EHS policy and SDS date."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property and classification records.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#60A5FA",
            ),
            "DMSO": candidate(
                (18.4, 16.4, 10.2), 189.0,
                "High boiling point may make removal difficult; evaporation behavior is not predicted by HSP.",
                "Check current SDS and process-specific exposure assessment.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Reported solvent behavior must be checked at the proposed temperature and concentration.",
                {"status": "WARN", "rule": "Process evidence needed", "reason": "High boiling point can create a recovery/drying trade-off.", "source": "Boiling-point record; verify with reliable property source."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#F87171",
            ),
            "Ethanol": candidate(
                (15.8, 8.8, 19.4), 78.0,
                "Volatile solvent; film formation depends on evaporation history and formulation.",
                "Flammability controls must be assessed from the current SDS.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Commonly used in ionomer dispersions and inks; behavior is formulation-dependent.",
                {"status": "WARN", "rule": "Flammability-control review", "reason": "A flammable solvent requires process-specific controls.", "source": "Current SDS and institutional laboratory rules."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#A78BFA",
            ),
            "Isopropanol": candidate(
                (15.8, 6.1, 16.4), 82.0,
                "Volatile solvent; evaporation and mixture composition influence coating formation.",
                "Flammability controls must be assessed from the current SDS.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Reported as a component in Nafion dispersions; standalone HSP screening does not determine final morphology.",
                {"status": "WARN", "rule": "Flammability-control review", "reason": "A flammable solvent requires process-specific controls.", "source": "Current SDS and institutional laboratory rules."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#FBBF24",
            ),
            "Water": candidate(
                (15.5, 16.0, 42.3), 100.0,
                "Acts differently from an organic solvent; can affect hydration, aggregation, and phase evolution.",
                "No special chemical hazard classification is assumed here; review full process hazards separately.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Useful as a component of ionomer dispersion systems; not a molecular-solubility guarantee for the backbone.",
                {"status": "PASS", "rule": "No case-specific mandatory exclusion", "reason": "Status does not imply successful membrane processing.", "source": "Teaching-case assumption; inspect all processing evidence."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#38BDF8",
            ),
        },
        "conflict_case": {
            "system": "Isopropanol/water mixture at 80:20 volume ratio",
            "model_prediction": "The adopted HSP model produces separate domain RED values and may classify the composite screening result as borderline or unfavorable.",
            "reported_observation": "Reported Nafion dispersions can be processable in alcohol/water media; the precise outcome depends on ionomer form, concentration, temperature, mixing history, and casting conditions.",
            "evidence_source": "Welch et al. (2012), ACS Macro Letters 1, 1403–1407; Safronova et al. (2021), Membranes 11, 782.",
            "student_prompt": "Compare the calculated domain results with the reported observation. What can the HSP model explain, and what information is still needed before recommending this mixture?",
            "limitations": "HSP screening does not establish particle morphology, phase evolution, evaporation, processing history, or kinetic stability.",
        },
    },
    "cellulose_acetate": {
        "id": "cellulose_acetate",
        "name": "Case 2: Cellulose Acetate Membrane Fabrication",
        "discipline": "Separation processes and membrane engineering",
        "problem_statement": (
            "Evaluate solvent candidates for cellulose acetate membrane fabrication considering compatibility, "
            "processability, literature evidence, and applicable regulatory/EHS constraints."
        ),
        "application": "Phase-inversion or solution-casting studies for cellulose acetate membranes.",
        "student_guidance": "Moderate: compare candidates and distinguish mandatory constraints from preferences.",
        "default_blend": ("Acetone", "Ethyl lactate", 60.0),
        "criteria_prompt": (
            "Define the processing target (for example, solution casting or phase inversion) and separate mandatory "
            "requirements from desirable but negotiable criteria."
        ),
        "polymer": {
            "name": "Cellulose acetate, membrane-grade teaching record",
            "type": "single-domain teaching representation",
            "material_specification": "Degree of substitution, molar mass, and supplier grade must be verified before using values for an experimental decision.",
            "source": source_record(
                "Hansen (2007) and adopted teaching values; verify against a source matching the selected cellulose acetate grade.",
                "Adopted literature-derived teaching values", "Not reported", "Single-domain screening representation", "Grade-dependent"
            ),
            "domains": {
                "Polymer": {
                    "hsp": np.array([18.6, 12.7, 11.0]),
                    "R0": 7.6,
                    "description": "Cellulose acetate teaching representation",
                    "color": "#10B981",
                    "source": source_record("Adopted teaching record; verify degree of substitution and source before use.", "Adopted", "Not reported", "HSP representation", "Grade-dependent"),
                },
            },
        },
        "solvents": {
            "NMP": candidate(
                (18.0, 12.3, 7.2), 202.0,
                "High boiling point may influence solvent removal and membrane morphology; measure or source process data.",
                "Applicable classifications and restrictions are jurisdiction-, date-, and use-specific. Consult current SDS and regulations.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Often reported as a strong cellulose-acetate solvent; thermodynamic screening alone is not a final decision.",
                {"status": "WARN", "rule": "Instructor-defined regulatory review", "reason": "This tool does not claim blanket prohibition. A course may define a hard gate only with a documented jurisdiction, date, and scope.", "source": "Current regulation and SDS must be cited by the instructor/student."},
                source_record("Hansen (2007) adopted solvent coordinates; classification and regulation require current source checks.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#EF4444",
            ),
            "Acetone": candidate(
                (15.5, 10.4, 7.0), 56.0,
                "High volatility can affect evaporation rate and skin formation; process evidence is required.",
                "Flammability controls must be assessed from the current SDS.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Reported as a solvent or co-solvent in some cellulose-acetate formulations; assess casting conditions.",
                {"status": "WARN", "rule": "Flammability and volatility review", "reason": "Volatility is a process trade-off, not an automatic rejection.", "source": "Current SDS and process-specific evidence."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#F59E0B",
            ),
            "Ethyl lactate": candidate(
                (16.0, 7.6, 12.5), 154.0,
                "Intermediate-to-high boiling point may moderate evaporation; membrane formation must be demonstrated experimentally.",
                "Review current SDS and concentration-specific handling requirements.",
                "Bio-derived origin is not treated as proof of overall environmental superiority; life-cycle and recovery evidence are needed.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Reported for polymer processing in selected studies; check the actual polymer grade and formulation.",
                {"status": "WARN", "rule": "Evidence-completeness review", "reason": "Process and recovery evidence are needed before a sustainability conclusion.", "source": "Literature/SDS/process data selected by instructor or student."},
                source_record("Hansen (2007) adopted solvent coordinates; verify solvent and process literature.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#10B981",
            ),
            "2-MeTHF": candidate(
                (16.9, 5.1, 4.4), 80.0,
                "Volatility can support removal but may affect casting and process safety.",
                "Flammability and peroxide-forming considerations must be evaluated from current sources.",
                "Bio-derived origin is not treated as proof of overall environmental superiority.",
                "Relative cost/availability is not quantified in this teaching database.",
                "May act as a co-solvent or swelling agent depending on composition; verify with experimental evidence.",
                {"status": "WARN", "rule": "Flammability/control review", "reason": "Use depends on the selected process and control measures.", "source": "Current SDS and process-specific evidence."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#3B82F6",
            ),
            "Dichloromethane": candidate(
                (18.2, 6.3, 6.1), 40.0,
                "Very volatile; may produce rapid phase evolution and requires robust containment.",
                "Hazard classification and local restrictions must be checked against current authoritative sources.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "May be reported in historic membrane recipes; process use requires specific evidence and controls.",
                {"status": "WARN", "rule": "Jurisdiction-specific review", "reason": "The software does not claim a universal ban.", "source": "Current local regulation and SDS."},
                source_record("Hansen (2007) adopted solvent coordinates; classification requires current source verification.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#94A3B8",
            ),
            "Water": candidate(
                (15.5, 16.0, 42.3), 100.0,
                "Often relevant as a non-solvent or coagulation-bath component in phase inversion; its role is process-specific.",
                "No special chemical hazard classification is assumed here; assess full process context.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Commonly used as a non-solvent in phase inversion; high RED can be useful for a non-solvent role.",
                {"status": "PASS", "rule": "No case-specific mandatory exclusion", "reason": "Pass does not mean it is a good molecular solvent.", "source": "Teaching-case role definition."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#38BDF8",
            ),
        },
        "conflict_case": {
            "system": "A low-RED candidate compared with a blend selected for process and EHS reasons",
            "model_prediction": "The lowest RED candidate has the strongest modelled thermodynamic affinity under the adopted HSP parameters.",
            "reported_observation": "Membrane fabrication also depends on solvent removal, nonsolvent exchange, polymer concentration, additives, and applicable EHS requirements.",
            "evidence_source": "Students must identify case-specific primary literature and current regulatory/SDS sources.",
            "student_prompt": "Which criteria are hard constraints in your stated context? Which are soft preferences? Explain why a low RED alone cannot determine the selected formulation.",
            "limitations": "The tool does not predict membrane pore structure, flux, mechanical properties, or compliance with any jurisdiction.",
        },
    },
    "polystyrene": {
        "id": "polystyrene",
        "name": "Case 3: Polystyrene/EPS Dissolution and Recovery",
        "discipline": "Circular materials and process design",
        "problem_statement": (
            "Evaluate whether alternative solvents can provide technically and environmentally credible options for dissolution "
            "of polystyrene or expanded-polystyrene (EPS) feedstock, including recovery of polymer and solvent."
        ),
        "application": "Dissolution, volume reduction, and recovery of polystyrene/EPS feedstock.",
        "student_guidance": "Lower: define the operational meaning of circular recovery and justify trade-offs with evidence.",
        "default_blend": ("d-Limonene", "Ethyl acetate", 80.0),
        "criteria_prompt": (
            "Define the feedstock, required recovery route, solvent handling boundary, and the measurable basis for environmental and process comparisons."
        ),
        "polymer": {
            "name": "General-purpose atactic polystyrene teaching record",
            "type": "single-domain teaching representation",
            "material_specification": "This case can represent general PS or EPS only after the feedstock, additives, contamination, and molecular-weight distribution are specified.",
            "source": source_record(
                "Hansen (2007), adopted teaching values; verify against a source matching the selected PS/EPS grade.",
                "Adopted literature-derived teaching values", "Not reported", "Single-domain screening representation", "Grade- and temperature-dependent"
            ),
            "domains": {
                "Polymer": {
                    "hsp": np.array([18.5, 4.5, 2.9]),
                    "R0": 5.3,
                    "description": "Atactic polystyrene teaching representation",
                    "color": "#8B5CF6",
                    "source": source_record("Adopted teaching record; verify grade and source before use.", "Adopted", "Not reported", "HSP representation", "Grade-dependent"),
                },
            },
        },
        "solvents": {
            "Toluene": candidate(
                (18.0, 1.4, 2.0), 111.0,
                "Moderate boiling point can support thermal recovery; actual recovery duty and emissions must be quantified.",
                "Hazard classifications and control measures require a current SDS and local regulatory review.",
                "No universal environmental conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Frequently reported as a strong polystyrene solvent; results depend on feedstock and operating conditions.",
                {"status": "WARN", "rule": "Exposure and VOC-control review", "reason": "A candidate can be technically strong while requiring substantial controls.", "source": "Current SDS and local process-safety requirements."},
                source_record("Hansen (2007) adopted solvent coordinates; verify current SDS and property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#EF4444",
            ),
            "d-Limonene": candidate(
                (17.2, 1.8, 4.3), 176.0,
                "High boiling point may complicate thermal solvent recovery; alternative separation routes should be evaluated.",
                "Current SDS and oxidation/product-quality concerns must be evaluated.",
                "Bio-based origin is an input to assessment, not a universal environmental rating.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Reported as a polystyrene/EPS solvent in selected studies; recovery conditions remain central to the decision.",
                {"status": "WARN", "rule": "Recovery-route evidence required", "reason": "High boiling point creates an explicit recovery trade-off.", "source": "Process-specific recovery data and current SDS."},
                source_record("Hansen (2007) adopted solvent coordinates; verify source and recovery literature.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#10B981",
            ),
            "Ethyl acetate": candidate(
                (15.8, 5.3, 7.2), 77.0,
                "Volatility may aid recovery but affects flammability and mixture behavior.",
                "Flammability controls must be assessed from a current SDS.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "May be considered as a co-solvent or process component; verify actual dissolution and recovery evidence.",
                {"status": "WARN", "rule": "Flammability-control review", "reason": "Process suitability depends on handling and recovery design.", "source": "Current SDS and process evidence."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#3B82F6",
            ),
            "MEK": candidate(
                (16.0, 9.0, 5.1), 80.0,
                "Volatility may aid removal but can alter process safety and solvent loss.",
                "Flammability and exposure controls must be assessed from a current SDS.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Reported behavior is feedstock- and condition-dependent.",
                {"status": "WARN", "rule": "Flammability-control review", "reason": "Requires process-specific handling controls.", "source": "Current SDS."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#F59E0B",
            ),
            "Cyclohexane": candidate(
                (16.8, 0.0, 0.2), 81.0,
                "Boiling point alone does not establish dissolution kinetics, recovery, or suitability.",
                "Flammability and exposure controls must be assessed from a current SDS.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "Reported solvent quality can depend strongly on temperature and polymer molecular weight.",
                {"status": "WARN", "rule": "Temperature/process evidence required", "reason": "Model result must be checked under selected conditions.", "source": "Current SDS and primary process literature."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#A78BFA",
            ),
            "Ethanol": candidate(
                (15.8, 8.8, 19.4), 78.0,
                "Can be considered as a precipitation/non-solvent component; role depends on the recovery route.",
                "Flammability controls must be assessed from a current SDS.",
                "No universal sustainability conclusion is assigned by this tool.",
                "Relative cost/availability is not quantified in this teaching database.",
                "High RED may be useful for precipitation rather than dissolution.",
                {"status": "WARN", "rule": "Flammability-control review", "reason": "Pass/fail depends on its defined process role.", "source": "Current SDS and recovery design evidence."},
                source_record("Hansen (2007) adopted solvent coordinates; verify physical-property source.", "Adopted", "Reference conditions", "Database value", "Not reported"),
                "#38BDF8",
            ),
        },
        "conflict_case": {
            "system": "Conventional aromatic solvent versus alternative solvent with a different recovery burden",
            "model_prediction": "More than one candidate may have RED below 1; this does not establish that they have equal process or environmental performance.",
            "reported_observation": "Dissolution, polymer recovery, residual solvent, emissions, and energy demand depend on the operating route and feedstock.",
            "evidence_source": "Students must identify process-specific primary literature and credible EHS/environmental sources.",
            "student_prompt": "Define ‘circular recovery’ operationally. Which evidence would allow you to compare solvent recovery, residual solvent, emissions, and polymer quality?",
            "limitations": "HSP does not calculate life-cycle impact, solvent recovery energy, occupational exposure, or recycled polymer properties.",
        },
    },
}


# -----------------------------------------------------------------------------
# HSP calculations, transparent scoring, statuses, and audit information
# -----------------------------------------------------------------------------

def calc_ra(point_a: Iterable[float], point_b: Iterable[float]) -> float:
    """Hansen distance: Ra = sqrt(4ΔD² + ΔP² + ΔH²)."""
    a = np.asarray(point_a, dtype=float)
    b = np.asarray(point_b, dtype=float)
    d = a - b
    return float(math.sqrt(4.0 * d[0] ** 2 + d[1] ** 2 + d[2] ** 2))


def calc_red(point: Iterable[float], target: Iterable[float], radius: float) -> float:
    """Relative Energy Difference, RED = Ra/R0."""
    if radius <= 0:
        raise ValueError("R0 must be positive.")
    return calc_ra(point, target) / float(radius)


def calc_mixed_hsp(solvents: Dict[str, dict], names: List[str], ratios: List[float]) -> np.ndarray:
    """Linear volume-fraction mixing of δD, δP and δH; fractions are normalized."""
    if len(names) != len(ratios) or not names:
        raise ValueError("Supply the same non-zero number of solvent names and ratios.")
    fractions = np.asarray(ratios, dtype=float)
    if np.any(fractions < 0) or fractions.sum() <= 0:
        raise ValueError("Ratios must be non-negative and sum to a value greater than zero.")
    fractions = fractions / fractions.sum()
    mixed = np.zeros(3, dtype=float)
    for name, fraction in zip(names, fractions):
        mixed += fraction * solvents[name]["hsp"]
    return mixed


def red_status(red: float) -> Tuple[str, str]:
    """A screening interpretation, not an experimental success classification."""
    if red < 1.0:
        return "Inside adopted HSP sphere", DARK["green"]
    if red <= 1.2:
        return "Boundary screening region", DARK["amber"]
    return "Outside adopted HSP sphere", DARK["red"]


def evaluate_domains(polymer: dict, hsp_point: Iterable[float]) -> Dict[str, dict]:
    """Compute domain-specific Ra and RED without concealing multi-domain values."""
    result = {}
    for name, domain in polymer["domains"].items():
        ra = calc_ra(hsp_point, domain["hsp"])
        red = ra / domain["R0"]
        label, color = red_status(red)
        result[name] = {
            "ra": ra,
            "red": red,
            "R0": domain["R0"],
            "screening": label,
            "color": color,
        }
    return result


def composite_red(domain_results: Dict[str, dict]) -> Tuple[float, str]:
    """Return arithmetic mean only as a labelled teaching indicator, never a final decision."""
    values = [row["red"] for row in domain_results.values()]
    return float(np.mean(values)), "Unweighted arithmetic mean of domain RED values; not a physical prediction or decision rule."


def evidence_status(candidate_record: dict) -> str:
    """Return transparent evidence status based on documented fields, not presumed performance."""
    gate = candidate_record.get("hard_gate", {}).get("status", "WARN")
    source = candidate_record.get("source", {})
    if gate == "FAIL":
        return "HARD FAIL"
    if source.get("value_type", "").lower().startswith("expert"):
        return "WARN: expert-derived data"
    if not source.get("citation") or source.get("citation") == "Not reported":
        return "WARN: source incomplete"
    return "PASS: evidence record present" if gate == "PASS" else "WARN: review required"


def hard_gate_status(candidate_record: dict) -> Tuple[str, str]:
    """Hard gate has precedence over all soft scoring."""
    gate = candidate_record["hard_gate"]
    status = gate.get("status", "WARN").upper()
    if status not in {"PASS", "WARN", "FAIL"}:
        status = "WARN"
    return status, gate.get("reason", "No reason documented.")


def map_red_to_indicator(avg_red: float) -> Tuple[float, str]:
    """Explicit pedagogical normalization from RED to 0–100 screening indicator."""
    if avg_red <= 0.6:
        return 100.0, "RED ≤ 0.60 mapped to 100"
    if avg_red <= 1.0:
        return 100.0 - (avg_red - 0.6) * 75.0, "Linear map between RED 0.60 and 1.00"
    if avg_red <= 1.5:
        return 70.0 - (avg_red - 1.0) * 100.0, "Linear map between RED 1.00 and 1.50"
    return max(0.0, 20.0 - (avg_red - 1.5) * 20.0), "Decreasing map above RED 1.50"


def qualitative_indicator(record: dict, field: str) -> Tuple[Optional[float], str, str]:
    """
    Convert documented teaching evidence into an explicitly qualitative indicator.
    Returns score, status, and explanation. None means that no score is claimed.
    """
    text = record.get(field, "").strip()
    if not text:
        return None, "MISSING", "No evidence text was supplied."
    # Indicators are intentionally not numerical claims about the physical property.
    # Their sole purpose is to make the assumptions inspectable and adjustable.
    if field == "process_evidence":
        bp = record.get("bp_c")
        if bp is None:
            return None, "MISSING", "No boiling-point/process record is available."
        if 80 <= bp <= 150:
            return 70.0, "QUALITATIVE", "Teaching indicator: moderate boiling-point range; process evidence still required."
        if 55 <= bp < 80 or 150 < bp <= 185:
            return 55.0, "QUALITATIVE", "Teaching indicator: volatility/recovery trade-off; not a measured processability score."
        return 40.0, "QUALITATIVE", "Teaching indicator: substantial volatility or recovery trade-off; not a measured processability score."
    if field == "ehs_evidence":
        return 50.0, "QUALITATIVE", "Neutral placeholder indicator. Inspect the cited SDS/regulatory evidence; do not interpret as a safety rating."
    if field == "sustainability_evidence":
        return 50.0, "QUALITATIVE", "Neutral placeholder indicator. No universal environmental score is claimed."
    if field == "cost_evidence":
        return None, "MISSING", "Cost/availability is not quantified in the teaching database."
    return None, "MISSING", "No mapping is defined."


def compute_evidence_profile(candidate_record: dict, domain_results: Dict[str, dict], weights: List[float]) -> dict:
    """
    Calculate a transparent conditional comparison index.

    This output is a sensitivity aid, not a recommendation. Missing values are excluded
    from the normalized denominator and reported explicitly. A hard FAIL cannot obtain
    a PASS overall status.
    """
    raw_weights = np.asarray(weights, dtype=float)
    if len(raw_weights) != 5 or np.any(raw_weights < 0) or raw_weights.sum() <= 0:
        raw_weights = np.asarray(DEFAULT_WEIGHTS, dtype=float)
    normalized_requested = raw_weights / raw_weights.sum()

    avg_red, composite_note = composite_red(domain_results)
    thermo, thermo_note = map_red_to_indicator(avg_red)
    process, process_status, process_note = qualitative_indicator(candidate_record, "process_evidence")
    ehs, ehs_status, ehs_note = qualitative_indicator(candidate_record, "ehs_evidence")
    sustain, sustain_status, sustain_note = qualitative_indicator(candidate_record, "sustainability_evidence")
    cost, cost_status, cost_note = qualitative_indicator(candidate_record, "cost_evidence")

    labels = WEIGHT_LABELS
    scores = [thermo, process, ehs, sustain, cost]
    statuses = ["MODEL-DERIVED", process_status, ehs_status, sustain_status, cost_status]
    notes = [thermo_note + "; " + composite_note, process_note, ehs_note, sustain_note, cost_note]

    usable = [(score, weight) for score, weight in zip(scores, normalized_requested) if score is not None]
    if usable:
        denominator = sum(weight for _, weight in usable)
        conditional_index = sum(score * weight for score, weight in usable) / denominator
    else:
        denominator = 0.0
        conditional_index = None

    gate_status, gate_reason = hard_gate_status(candidate_record)
    missing = [label for label, score in zip(labels, scores) if score is None]
    if gate_status == "FAIL":
        overall_status = "HARD FAIL — mandatory constraint violated"
    elif missing:
        overall_status = "WARN — comparison incomplete"
    elif gate_status == "WARN" or any(s in {"QUALITATIVE", "MISSING"} for s in statuses[1:]):
        overall_status = "WARN — evidence/trade-off review required"
    else:
        overall_status = "PASS — available evidence satisfies stated gate"

    return {
        "conditional_index": conditional_index,
        "overall_status": overall_status,
        "hard_gate": gate_status,
        "hard_gate_reason": gate_reason,
        "average_red": avg_red,
        "composite_note": composite_note,
        "labels": labels,
        "scores": scores,
        "statuses": statuses,
        "notes": notes,
        "requested_weights": normalized_requested.tolist(),
        "available_weight_fraction": denominator,
        "missing_criteria": missing,
    }


def weight_sensitivity(case_data: dict, weights: List[float], variation: float = 0.20) -> dict:
    """Check whether the top conditional comparison index changes when one weight changes ±variation."""
    solvents = case_data["solvents"]
    polymer = case_data["polymer"]
    baseline = {}
    for name, record in solvents.items():
        profile = compute_evidence_profile(record, evaluate_domains(polymer, record["hsp"]), weights)
        if profile["conditional_index"] is not None and profile["hard_gate"] != "FAIL":
            baseline[name] = profile["conditional_index"]
    baseline_top = max(baseline, key=baseline.get) if baseline else None

    variations = []
    base = np.asarray(weights, dtype=float)
    for i, label in enumerate(WEIGHT_LABELS):
        for direction in (-1, 1):
            trial = base.copy()
            trial[i] = max(0.0, trial[i] * (1 + direction * variation))
            ranked = {}
            for name, record in solvents.items():
                profile = compute_evidence_profile(record, evaluate_domains(polymer, record["hsp"]), trial.tolist())
                if profile["conditional_index"] is not None and profile["hard_gate"] != "FAIL":
                    ranked[name] = profile["conditional_index"]
            top = max(ranked, key=ranked.get) if ranked else None
            variations.append({"criterion": label, "direction": direction, "top": top})

    changed = any(item["top"] != baseline_top for item in variations if baseline_top is not None)
    return {
        "baseline_top": baseline_top,
        "changed": changed,
        "message": (
            "The conditional comparison leader changes under a ±20% single-weight test; treat the decision as weight-sensitive."
            if changed else
            "The conditional comparison leader does not change under the implemented ±20% single-weight test. This does not prove physical suitability."
        ),
        "variations": variations,
    }


def audit_rows(case_data: dict) -> List[dict]:
    """Generate the requested decision-engine audit table."""
    rows = []
    polymer = case_data["polymer"]
    rows.append({
        "component": "HSP values", "what": "δD, δP, δH polymer/solvent coordinates", "units": HSP_UNITS,
        "source": polymer["source"]["citation"], "type": polymer["source"]["value_type"],
        "normalization": "None; used in Ra equation", "status_logic": "Evidence source must be inspected", "missing_data": "Warn; do not calculate with undocumented values.",
    })
    rows.append({
        "component": "R0", "what": "HSP interaction radius, domain-specific where used", "units": HSP_UNITS,
        "source": "Domain-specific source record in the case database", "type": "Adopted teaching value",
        "normalization": "RED = Ra/R0", "status_logic": "R0 must be > 0", "missing_data": "Stop RED calculation and flag missing radius.",
    })
    rows.append({
        "component": "RED", "what": "HSP distance normalized by R0", "units": "dimensionless",
        "source": "Equation implemented in calc_ra and calc_red", "type": "Model-derived",
        "normalization": "Ra = sqrt(4ΔD²+ΔP²+ΔH²); RED = Ra/R0", "status_logic": "<1 inside adopted sphere; 1–1.2 boundary; >1.2 outside. Screening only.", "missing_data": "No RED displayed if HSP/R0 is missing.",
    })
    rows.append({
        "component": "Process evidence", "what": "Boiling point plus qualitative process note", "units": "°C plus narrative",
        "source": "Candidate source record", "type": "Qualitative teaching indicator",
        "normalization": "Transparent 0–100 teaching indicator based on boiling-point bands; not measured processability", "status_logic": "WARN where process data are incomplete or trade-offs exist", "missing_data": "Excluded from conditional index; warning displayed.",
    })
    rows.append({
        "component": "EHS/regulatory evidence", "what": "Narrative hazard and regulatory evidence", "units": "categorical/narrative",
        "source": "Current SDS and jurisdiction-specific regulation required", "type": "Qualitative / source-dependent",
        "normalization": "Neutral placeholder 50 only when a narrative exists; not a safety score", "status_logic": "HARD FAIL only if instructor documents a true mandatory rule", "missing_data": "Excluded from conditional index; warning displayed.",
    })
    rows.append({
        "component": "Sustainability indicator", "what": "Narrative evidence only", "units": "categorical/narrative",
        "source": "Student/instructor must cite metric or source", "type": "Qualitative / source-dependent",
        "normalization": "Neutral placeholder 50; no universal green rating claimed", "status_logic": "Warn unless evidence scope is documented", "missing_data": "Excluded from conditional index; warning displayed.",
    })
    rows.append({
        "component": "Cost/availability evidence", "what": "Price, availability, or procurement data", "units": "currency/date or narrative",
        "source": "Not supplied in this teaching database", "type": "Missing by default",
        "normalization": "None until a dated, comparable source is supplied", "status_logic": "Warn", "missing_data": "Excluded from conditional index; warning displayed.",
    })
    rows.append({
        "component": "Weights", "what": "Student-adjustable criterion priorities", "units": "%",
        "source": "Default teaching baseline: 30/25/25/10/10", "type": "Instructional assumption",
        "normalization": "Weights normalized across available indicators", "status_logic": "Sensitivity test reports ranking changes", "missing_data": "Invalid weights reset to baseline.",
    })
    rows.append({
        "component": "Hard constraints", "what": "Documented mandatory rule", "units": "PASS/WARN/FAIL",
        "source": "Case-specific gate record", "type": "Policy/process-context dependent",
        "normalization": "No soft-score compensation permitted for HARD FAIL", "status_logic": "FAIL overrides conditional index", "missing_data": "Warn; no hard gate is assumed without documentation.",
    })
    rows.append({
        "component": "Conditional comparison index", "what": "Weighted mean of available model/teaching indicators", "units": "0–100",
        "source": "compute_evidence_profile", "type": "Derived teaching aid",
        "normalization": "Weights renormalized after missing indicators are excluded", "status_logic": "Never labelled ‘best’ or experimental success", "missing_data": "Available-weight fraction is reported.",
    })
    return rows


# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------

def style_axis(ax):
    ax.set_facecolor(DARK["axes"])
    for spine in ax.spines.values():
        spine.set_color(DARK["border"])
    ax.tick_params(colors=DARK["subtext"], labelsize=8)
    ax.xaxis.label.set_color(DARK["text"])
    ax.yaxis.label.set_color(DARK["text"])
    ax.grid(True, color=DARK["grid"], linestyle="--", linewidth=0.6, alpha=0.7)


def plot_hsp_projections(case_data: dict, mixture_hsp: np.ndarray, mixture_label: str):
    polymer = case_data["polymer"]
    solvents = case_data["solvents"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), facecolor=DARK["figure"])
    fig.suptitle(f"HSP evidence view — {case_data['name']}", color="white", fontsize=12, fontweight="bold")
    planes = [(0, 1, "δD (MPa½)", "δP (MPa½)"), (0, 2, "δD (MPa½)", "δH (MPa½)"), (1, 2, "δP (MPa½)", "δH (MPa½)")]
    for ax, (xidx, yidx, xlabel, ylabel) in zip(axes, planes):
        style_axis(ax)
        for dname, domain in polymer["domains"].items():
            hsp = domain["hsp"]
            r0 = domain["R0"]
            width = 2 * (r0 / 2 if xidx == 0 else r0)
            height = 2 * (r0 / 2 if yidx == 0 else r0)
            ell = Ellipse((hsp[xidx], hsp[yidx]), width=width, height=height,
                          edgecolor=domain["color"], facecolor=domain["color"], alpha=0.13, linestyle="--", linewidth=1.6)
            ax.add_patch(ell)
            ax.scatter(hsp[xidx], hsp[yidx], s=62, color=domain["color"], edgecolors="white", linewidths=0.7, zorder=5)
            ax.text(hsp[xidx] + 0.25, hsp[yidx] + 0.35, dname, color=domain["color"], fontsize=7.5, fontweight="bold")
        for sname, sinfo in solvents.items():
            h = sinfo["hsp"]
            ax.scatter(h[xidx], h[yidx], s=42, color=sinfo["color"], edgecolors="white", linewidths=0.5, zorder=4)
            ax.text(h[xidx] + 0.2, h[yidx] + 0.25, sname, color=sinfo["color"], fontsize=6.8)
        ax.scatter(mixture_hsp[xidx], mixture_hsp[yidx], s=150, marker="*", color=DARK["amber"], edgecolors="white", linewidths=0.8, zorder=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title("HSP projection", color="white", fontsize=9)
    fig.text(0.5, 0.02, f"Mixture: {mixture_label}. Ellipses represent adopted HSP spheres; they do not predict processing success.",
             ha="center", color=DARK["subtext"], fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))
    return fig


def plot_red_evidence(case_data: dict, mixture_hsp: np.ndarray, mixture_label: str):
    polymer = case_data["polymer"]
    solvents = case_data["solvents"]
    domains = list(polymer["domains"].items())
    fig, axes = plt.subplots(1, len(domains), figsize=(6.4 * len(domains), 4.8), facecolor=DARK["figure"], squeeze=False)
    fig.suptitle("Domain-specific RED screening evidence", color="white", fontsize=12, fontweight="bold")
    points = [(name, info["hsp"], info["color"]) for name, info in solvents.items()]
    points.append((f"Mixture\n{mixture_label}", mixture_hsp, DARK["amber"]))
    for ax, (dname, domain) in zip(axes[0], domains):
        style_axis(ax)
        values = [calc_red(point, domain["hsp"], domain["R0"]) for _, point, _ in points]
        names = [name for name, _, _ in points]
        colors = [color for _, _, color in points]
        bars = ax.bar(names, values, color=colors, edgecolor=DARK["border"], width=0.62)
        ax.axhline(1.0, color="white", linestyle="--", linewidth=1.1, label="RED = 1 boundary")
        ax.axhline(1.2, color=DARK["amber"], linestyle=":", linewidth=1.1, label="RED = 1.2 boundary")
        for bar, value in zip(bars, values):
            label, color = red_status(value)
            ax.text(bar.get_x() + bar.get_width() / 2, value + 0.03, f"{value:.2f}", ha="center", va="bottom", fontsize=7.4, color=color)
        ax.set_ylim(0, max(1.4, max(values) * 1.2))
        ax.set_ylabel("RED = Ra/R0")
        ax.set_title(f"{dname} domain (R0 = {domain['R0']:.1f})", color="white", fontsize=10)
        ax.tick_params(axis="x", rotation=25)
        ax.legend(fontsize=7, facecolor=DARK["axes"], labelcolor="white")
    fig.text(0.5, 0.02, "RED is a thermodynamic screening indicator. It does not establish dissolution rate, morphology, safety, or compliance.",
             ha="center", color=DARK["subtext"], fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))
    return fig


def plot_mixture_sensitivity(case_data: dict, solvent_a: str, solvent_b: str, selected_ratio: float):
    solvents = case_data["solvents"]
    polymer = case_data["polymer"]
    x = np.linspace(0, 100, 151)
    fig, ax = plt.subplots(figsize=(8.2, 4.8), facecolor=DARK["figure"])
    style_axis(ax)
    all_domain = {name: [] for name in polymer["domains"]}
    average = []
    for fraction in x:
        mixed = calc_mixed_hsp(solvents, [solvent_a, solvent_b], [fraction, 100 - fraction])
        result = evaluate_domains(polymer, mixed)
        values = []
        for name, row in result.items():
            all_domain[name].append(row["red"])
            values.append(row["red"])
        average.append(np.mean(values))
    for name, curve in all_domain.items():
        color = polymer["domains"][name]["color"]
        ax.plot(x, curve, color=color, linewidth=1.8, label=f"{name} RED")
    if len(all_domain) > 1:
        ax.plot(x, average, color=DARK["amber"], linewidth=2.0, linestyle="--", label="Arithmetic mean (teaching indicator)")
    ax.axhline(1.0, color="white", linestyle="--", linewidth=1.1, label="RED = 1 boundary")
    ax.axhline(1.2, color=DARK["amber"], linestyle=":", linewidth=1.1, label="RED = 1.2 boundary")
    index = int(round(selected_ratio))
    index = max(0, min(index, len(x) - 1))
    ax.axvline(selected_ratio, color=DARK["amber"], linewidth=1.2)
    ax.set_xlabel(f"Volume fraction of {solvent_a} (%) — remainder is {solvent_b}")
    ax.set_ylabel("RED")
    ax.set_title(f"Binary mixture sensitivity: {solvent_a}/{solvent_b}", color="white", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7.5, facecolor=DARK["axes"], labelcolor="white")
    fig.text(0.5, 0.02, "Linear volume-fraction mixing is assumed for δD, δP and δH. Explore sensitivity; do not select the minimum RED automatically.",
             ha="center", color=DARK["subtext"], fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return fig


def plot_evidence_profile(case_data: dict, weights: List[float]):
    solvents = case_data["solvents"]
    polymer = case_data["polymer"]
    selected = list(solvents.keys())[:5]
    categories = ["Thermo", "Process", "EHS", "Sustain.", "Cost"]
    fig, ax = plt.subplots(figsize=(9.2, 5.0), facecolor=DARK["figure"])
    style_axis(ax)
    x = np.arange(len(categories))
    width = 0.14
    for i, name in enumerate(selected):
        profile = compute_evidence_profile(solvents[name], evaluate_domains(polymer, solvents[name]["hsp"]), weights)
        values = [v if v is not None else np.nan for v in profile["scores"]]
        ax.bar(x + (i - (len(selected)-1)/2) * width, values, width=width, label=name, color=solvents[name]["color"], alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, color=DARK["text"])
    ax.set_ylim(0, 105)
    ax.set_ylabel("0–100 conditional teaching indicator")
    ax.set_title("Inspectable evidence indicators (not universal property scores)", color="white", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7.5, facecolor=DARK["axes"], labelcolor="white")
    fig.text(0.5, 0.02, "Process/EHS/sustainability indicators are qualitative teaching mappings. Missing cost data are shown as absent rather than invented.",
             ha="center", color=DARK["subtext"], fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return fig


# -----------------------------------------------------------------------------
# Export functions
# -----------------------------------------------------------------------------

def _as_text(value: Optional[str]) -> str:
    return value.strip() if value and value.strip() else "Not provided by student"


def export_audit_table(case_data: dict, filename: str) -> str:
    rows = audit_rows(case_data)
    with open(filename, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return filename


def export_candidate_evidence_csv(case_data: dict, filename: str) -> str:
    fields = [
        "candidate", "deltaD_MPa_half", "deltaP_MPa_half", "deltaH_MPa_half", "boiling_point_C",
        "process_evidence", "ehs_evidence", "sustainability_evidence", "cost_evidence", "literature_evidence",
        "hard_gate_status", "hard_gate_rule", "hard_gate_reason", "hard_gate_source",
        "source_citation", "source_value_type", "source_temperature", "source_method", "source_uncertainty",
    ]
    with open(filename, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for name, record in case_data["solvents"].items():
            h = record["hsp"]
            src = record["source"]
            gate = record["hard_gate"]
            writer.writerow({
                "candidate": name, "deltaD_MPa_half": f"{h[0]:.2f}", "deltaP_MPa_half": f"{h[1]:.2f}", "deltaH_MPa_half": f"{h[2]:.2f}",
                "boiling_point_C": "" if record["bp_c"] is None else f"{record['bp_c']:.1f}",
                "process_evidence": record["process_evidence"], "ehs_evidence": record["ehs_evidence"],
                "sustainability_evidence": record["sustainability_evidence"], "cost_evidence": record["cost_evidence"],
                "literature_evidence": record["literature_evidence"],
                "hard_gate_status": gate["status"], "hard_gate_rule": gate["rule"], "hard_gate_reason": gate["reason"], "hard_gate_source": gate["source"],
                "source_citation": src["citation"], "source_value_type": src["value_type"], "source_temperature": src["temperature"],
                "source_method": src["method"], "source_uncertainty": src["uncertainty"],
            })
    return filename


def export_material_selection_brief(case_data: dict, blend_info: Tuple[str, str, float, float],
                                    weights: List[float], student_fields: dict,
                                    filename: str) -> str:
    """Export the complete student-facing Material Selection Brief."""
    polymer = case_data["polymer"]
    solvents = case_data["solvents"]
    sol_a, sol_b, ratio_a, ratio_b = blend_info
    mixture = calc_mixed_hsp(solvents, [sol_a, sol_b], [ratio_a, ratio_b])
    mixture_domains = evaluate_domains(polymer, mixture)
    sensitivity = weight_sensitivity(case_data, weights)

    lines = []
    lines.append("# Material Selection Brief")
    lines.append("")
    lines.append(f"**Case:** {case_data['name']}  ")
    lines.append(f"**Student:** {_as_text(student_fields.get('student_name'))}  ")
    lines.append(f"**Date:** {date.today().isoformat()}  ")
    lines.append("")
    lines.append("## 1. Engineering problem and application")
    lines.append(case_data["problem_statement"])
    lines.append("")
    lines.append(f"**Application:** {case_data['application']}")
    lines.append(f"**Student problem definition:** {_as_text(student_fields.get('problem_definition'))}")
    lines.append(f"**Selected polymer:** {_as_text(student_fields.get('define_polymer'))}")
    lines.append(f"**Polarity and solvent class:** {_as_text(student_fields.get('define_polarity'))}; {_as_text(student_fields.get('define_class'))}")
    lines.append(f"**Catalog suggestion:** {_as_text(student_fields.get('define_suggested'))} (not validated for dissolution)")
    lines.append(f"**Definition of suitable / decision criteria:** {_as_text(student_fields.get('criteria_definition'))}")
    lines.append("")
    lines.append("## 2. Candidate evidence and sources")
    lines.append("| Candidate | HSP (δD, δP, δH) | Process evidence | EHS/regulatory evidence | Sustainability evidence | Literature/empirical evidence | Gate |")
    lines.append("|---|---|---|---|---|---|---|")
    for name, record in solvents.items():
        h = record["hsp"]
        gate, reason = hard_gate_status(record)
        lines.append(
            f"| {name} | ({h[0]:.1f}, {h[1]:.1f}, {h[2]:.1f}) {HSP_UNITS} | {record['process_evidence']} | "
            f"{record['ehs_evidence']} | {record['sustainability_evidence']} | {record['literature_evidence']} | {gate}: {reason} |"
        )
    lines.append("")
    lines.append("**Evidence audit note:** HSP entries are adopted teaching values. Use the exported evidence CSV/audit table to inspect citations, methods, temperatures, and uncertainty notes.")
    lines.append("")
    lines.append("## Student literature investigation")
    if student_fields.get("no_literature"):
        lines.append("No suitable literature data were found or the student could not identify relevant literature data.")
    elif not student_fields.get("literature_records"):
        lines.append("No literature records entered; search status not established.")
    else:
        for i, rec in enumerate(student_fields["literature_records"], 1):
            lines.append(f"### Student-entered paper {i}")
            for key, label in (("doi", "DOI / URL"),("title", "Title"),("authors", "Authors"),("year", "Year"),
                    ("polymer", "Polymer"),("solvent_a", "Solvent A"),("solvent_b", "Solvent B"),
                    ("ratio", "Ratio and basis"),("concentration", "Concentration"),("temperature", "Temperature"),
                    ("time", "Dissolution time"),("observation", "Reported observation"),("notes", "Notes")):
                lines.append(f"- {label}: {_as_text(rec.get(key))}")
    lines.append("These paper records are student-entered and have not been verified by the software.")
    lines.append("")
    lines.append("## Literature vs Model comparison")
    lines.append("| Paper | Paper polymer | Paper solvents | Paper ratio | Model solvents | Model ratio | Model RED | Comparison |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for rec in student_fields.get("literature_comparison", []):
        lines.append("| " + " | ".join(str(rec[key]).replace("|", "\\|").replace("\n", " ") for key in
            ("paper", "polymer", "paper_system", "paper_ratio", "model_system", "model_ratio", "red", "comparison")) + " |")
    if not student_fields.get("literature_comparison"):
        lines.append("| N/A | N/A | N/A | N/A | N/A | N/A | N/A | No direct comparison available |")
    lines.append("")
    lines.append("## 3. HSP/RED model evidence")
    lines.append("The program uses: `Ra = sqrt(4ΔδD² + ΔδP² + ΔδH²)` and `RED = Ra/R0`.")
    lines.append("For mixtures, δD, δP, and δH are mixed linearly by normalized volume fraction. This is an assumption, not a validated prediction for every mixture.")
    lines.append(f"**Selected mixture:** {ratio_a:.0f}% {sol_a} / {ratio_b:.0f}% {sol_b}")
    lines.append(f"**Mixture HSP:** ({mixture[0]:.2f}, {mixture[1]:.2f}, {mixture[2]:.2f}) {HSP_UNITS}")
    lines.append("| Domain | Ra | R0 | RED | Screening interpretation |")
    lines.append("|---|---:|---:|---:|---|")
    for name, row in mixture_domains.items():
        lines.append(f"| {name} | {row['ra']:.3f} | {row['R0']:.3f} | {row['red']:.3f} | {row['screening']} |")
    mean_red, mean_note = composite_red(mixture_domains)
    lines.append(f"\n**Composite teaching indicator:** {mean_red:.3f}. {mean_note}")
    lines.append("")
    lines.append("## 4. Comparison of credible alternatives")
    lines.append("| Candidate | Conditional comparison index | Status | Missing criteria | Gate reason |")
    lines.append("|---|---:|---|---|---|")
    for name, record in solvents.items():
        profile = compute_evidence_profile(record, evaluate_domains(polymer, record["hsp"]), weights)
        index = "Not calculated" if profile["conditional_index"] is None else f"{profile['conditional_index']:.1f}"
        missing = ", ".join(profile["missing_criteria"]) or "None"
        lines.append(f"| {name} | {index} | {profile['overall_status']} | {missing} | {profile['hard_gate_reason']} |")
    lines.append("")
    lines.append("**Weight-sensitivity observation:** " + sensitivity["message"])
    lines.append(f"**Student observation after changing weights:** {_as_text(student_fields.get('weight_observation'))}")
    lines.append("")
    lines.append("## 5. Model–evidence comparison")
    conflict = case_data["conflict_case"]
    lines.append(f"**Focus system:** {conflict['system']}")
    lines.append(f"**Model evidence:** {conflict['model_prediction']}")
    lines.append(f"**Reported/literature evidence:** {conflict['reported_observation']}")
    lines.append(f"**Evidence source:** {conflict['evidence_source']}")
    lines.append(f"**Student agreement/conflict analysis:** {_as_text(student_fields.get('conflict_analysis'))}")
    lines.append("")
    lines.append("## 6. Recommendation and supporting evidence")
    lines.append(f"**Recommendation / shortlist:** {_as_text(student_fields.get('recommendation'))}")
    lines.append(f"**Computational evidence:** {_as_text(student_fields.get('computational_evidence'))}")
    lines.append(f"**Literature/practical evidence:** {_as_text(student_fields.get('practical_evidence'))}")
    lines.append(f"**Main trade-off:** {_as_text(student_fields.get('main_tradeoff'))}")
    lines.append("")
    lines.append("## 7. Reflection and verification")
    lines.append(f"**Confidence level (evidence-confidence, not probability):** {_as_text(student_fields.get('confidence'))}")
    lines.append(f"**Main uncertainty or evidence gap:** {_as_text(student_fields.get('main_uncertainty'))}")
    lines.append(f"**What would change my decision?** {_as_text(student_fields.get('change_decision'))}")
    lines.append(f"**Proposed verification experiment and measurable endpoint:** {_as_text(student_fields.get('verification'))}")
    lines.append("")
    lines.append("## 8. Required student sentence")
    lines.append(
        "I recommend ___ because the computational evidence indicates ___, while the literature/practical evidence indicates ___. "
        "The main uncertainty is ___. I would change this decision if ___. The next verification step is ___."
    )
    lines.append("")
    lines.append("---")
    lines.append("Generated by PolymerSolv Studio. The report records available evidence and assumptions; it does not certify engineering suitability or regulatory compliance.")

    Path(filename).write_text("\n".join(lines), encoding="utf-8")
    return filename


def export_case_bundle(case_data: dict, output_dir: str = "PolymerSolv_exports") -> List[str]:
    """Headless export of figures, audit table, evidence table, and blank student brief."""
    directory = Path(output_dir) / case_data["id"]
    directory.mkdir(parents=True, exist_ok=True)
    sol_a, sol_b, ratio_a = case_data["default_blend"]
    ratio_b = 100.0 - ratio_a
    mix = calc_mixed_hsp(case_data["solvents"], [sol_a, sol_b], [ratio_a, ratio_b])
    mix_label = f"{ratio_a:.0f}% {sol_a} / {ratio_b:.0f}% {sol_b}"
    weights = DEFAULT_WEIGHTS[:]
    outputs = []
    figures = [
        ("01_HSP_projections.png", plot_hsp_projections(case_data, mix, mix_label)),
        ("02_RED_evidence.png", plot_red_evidence(case_data, mix, mix_label)),
        ("03_mixture_sensitivity.png", plot_mixture_sensitivity(case_data, sol_a, sol_b, ratio_a)),
        ("04_evidence_indicators.png", plot_evidence_profile(case_data, weights)),
    ]
    for name, figure in figures:
        path = directory / name
        figure.savefig(path, dpi=300, bbox_inches="tight", facecolor=figure.get_facecolor())
        plt.close(figure)
        outputs.append(str(path))
    outputs.append(export_audit_table(case_data, str(directory / "decision_engine_audit.csv")))
    outputs.append(export_candidate_evidence_csv(case_data, str(directory / "candidate_evidence.csv")))
    fields = {
        "student_name": "", "problem_definition": "", "criteria_definition": "", "weight_observation": "",
        "conflict_analysis": "", "recommendation": "", "computational_evidence": "", "practical_evidence": "",
        "main_tradeoff": "", "confidence": "", "main_uncertainty": "", "change_decision": "", "verification": "",
    }
    outputs.append(export_material_selection_brief(case_data, (sol_a, sol_b, ratio_a, ratio_b), weights, fields, str(directory / "Material_Selection_Brief_Template.md")))
    return outputs


# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------

class PolymerSolvGUI:
    """A direct, guided GUI that preserves student judgment and makes assumptions visible."""

    def __init__(self, root):
        self.root = root
        self.root.title(f"{APP_NAME} — Evidence-Based Material Selection")
        self.root.geometry("1400x900")
        self.root.minsize(1120, 720)
        self.root.configure(bg=DARK["figure"])
        self.current_case_id = "nafion"
        self.case_data = CASES_DATABASE[self.current_case_id]
        self.weights = DEFAULT_WEIGHTS[:]
        self.solvent_a, self.solvent_b, self.ratio_a = self.case_data["default_blend"]
        self.student_widgets = {}
        self.setup_styles()
        self.build_header()
        self.build_tabs()
        self.refresh_all()

    def setup_styles(self):
        self.style = ttk.Style()
        try:
            self.style.theme_use("clam")
        except Exception:
            pass
        self.style.configure(".", background=DARK["figure"], foreground=DARK["text"], font=("Segoe UI", 9))
        self.style.configure("TNotebook", background=DARK["figure"], borderwidth=0)
        self.style.configure("TNotebook.Tab", background="#1E293B", foreground=DARK["subtext"], padding=(12, 7))
        self.style.map("TNotebook.Tab", background=[("selected", "#0369A1")], foreground=[("selected", "white")])
        self.style.configure("Treeview", background=DARK["axes"], foreground=DARK["text"], fieldbackground=DARK["axes"], rowheight=24)
        self.style.configure("Treeview.Heading", background="#1E293B", foreground=DARK["blue"], font=("Segoe UI", 8, "bold"))
        self.style.map("Treeview", background=[("selected", "#0284C7")], foreground=[("selected", "white")])
        self.style.configure("Primary.TButton", background="#0369A1", foreground="white", padding=(9, 5))
        self.style.map("Primary.TButton", background=[("active", "#0284C7")])
        self.style.configure("Good.TLabel", background="#064E3B", foreground="#D1FAE5", padding=(7, 4))
        self.style.configure("Warn.TLabel", background="#78350F", foreground="#FEF3C7", padding=(7, 4))
        self.style.configure("Fail.TLabel", background="#7F1D1D", foreground="#FEE2E2", padding=(7, 4))
        self.style.configure("Info.TLabel", background="#0C4A6E", foreground="#E0F2FE", padding=(7, 4))

    def build_header(self):
        frame = tk.Frame(self.root, bg="#0F172A", height=64)
        frame.pack(fill="x")
        text = tk.Frame(frame, bg="#0F172A")
        text.pack(side="left", padx=16, pady=9)
        tk.Label(text, text=APP_NAME, bg="#0F172A", fg=DARK["blue"], font=("Segoe UI", 15, "bold")).pack(anchor="w")
        tk.Label(text, text="Guided evidence-based material-selection workflow | HSP/RED is evidence, not the decision", bg="#0F172A", fg=DARK["subtext"], font=("Segoe UI", 8)).pack(anchor="w")
        actions = tk.Frame(frame, bg="#0F172A")
        actions.pack(side="right", padx=15, pady=12)
        ttk.Button(actions, text="Export Material Selection Brief", style="Primary.TButton", command=self.export_report).pack(side="left", padx=4)
        ttk.Button(actions, text="Export Audit + Evidence CSV", style="Primary.TButton", command=self.export_tables).pack(side="left", padx=4)

    def build_tabs(self):
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=10, pady=10)
        self.tab_define = ttk.Frame(self.notebook)
        self.tab_investigate = ttk.Frame(self.notebook)
        self.tab_model = ttk.Frame(self.notebook)
        self.tab_compare = ttk.Frame(self.notebook)
        self.tab_decide = ttk.Frame(self.notebook)
        self.tab_reflect = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_define, text="1. DEFINE")
        self.notebook.add(self.tab_investigate, text="2. INVESTIGATE")
        self.notebook.add(self.tab_model, text="3. MODEL")
        self.notebook.add(self.tab_compare, text="4. COMPARE")
        self.notebook.add(self.tab_decide, text="5. DECIDE")
        self.notebook.add(self.tab_reflect, text="6. REFLECT")
        self.build_define_tab()
        self.build_investigate_tab()
        self.build_model_tab()
        self.build_compare_tab()
        self.build_decide_tab()
        self.build_reflect_tab()

    def directive(self, parent, text):
        label = ttk.Label(parent, text=text, style="Info.TLabel", wraplength=1280, justify="left")
        label.pack(fill="x", padx=10, pady=(10, 6))
        return label

    def build_define_tab(self):
        self.directive(self.tab_define, "DEFINE — Select an open engineering problem. Write what ‘suitable’ means before reviewing a ranking.")
        cards = tk.Frame(self.tab_define, bg=DARK["figure"])
        cards.pack(fill="x", padx=10, pady=4)
        self.case_var = tk.StringVar(value="nafion")
        metadata = [
            ("nafion", "Case 1: Nafion processing", "Evaluate solvent and mixture choices and investigate model–literature differences.", DARK["blue"]),
            ("cellulose_acetate", "Case 2: Cellulose acetate membranes", "Evaluate candidates using compatibility, processability, and applicable regulatory constraints.", DARK["green"]),
            ("polystyrene", "Case 3: Polystyrene/EPS recovery", "Evaluate technically and environmentally credible dissolution and recovery routes.", DARK["purple"]),
        ]
        for cid, title, body, color in metadata:
            card = tk.Frame(cards, bg=DARK["axes"], bd=1, relief="solid")
            card.pack(side="left", fill="both", expand=True, padx=5, pady=5)
            tk.Radiobutton(card, text=title, variable=self.case_var, value=cid, command=self.change_case,
                           bg=DARK["axes"], fg="white", selectcolor="#0F172A", activebackground=DARK["axes"], activeforeground="white",
                           font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=10, pady=(10, 4))
            tk.Label(card, text=body, bg=DARK["axes"], fg=DARK["subtext"], wraplength=350, justify="left", font=("Segoe UI", 8)).pack(anchor="w", padx=10, pady=(0, 10))
            tk.Frame(card, bg=color, height=3).pack(fill="x")
        body = tk.Frame(self.tab_define, bg=DARK["figure"])
        body.pack(fill="both", expand=True, padx=10, pady=6)
        left = tk.Frame(body, bg=DARK["axes"], bd=1, relief="solid")
        left.pack(side="left", fill="both", expand=True, padx=(0, 5))
        right = tk.Frame(body, bg=DARK["axes"], bd=1, relief="solid")
        right.pack(side="right", fill="both", expand=True, padx=(5, 0))
        tk.Label(left, text="Case context", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=(10, 5))
        self.define_context = tk.Text(left, height=12, bg="#0F172A", fg="white", wrap="word", relief="flat", font=("Segoe UI", 9))
        self.define_context.pack(fill="both", expand=True, padx=12, pady=(0, 12))
        tk.Label(right, text="Your problem definition", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=(10, 5))
        tk.Label(right, text="Write one paragraph: application, target outcome, and boundaries.", bg=DARK["axes"], fg=DARK["subtext"], font=("Segoe UI", 8)).pack(anchor="w", padx=12)
        self.student_widgets["problem_definition"] = tk.Text(right, height=6, bg="#0F172A", fg="white", wrap="word", relief="flat", font=("Segoe UI", 9))
        self.student_widgets["problem_definition"].pack(fill="x", padx=12, pady=6)
        tk.Label(right, text="Define ‘suitable’ and state criteria / hard constraints.", bg=DARK["axes"], fg=DARK["subtext"], font=("Segoe UI", 8)).pack(anchor="w", padx=12)
        self.student_widgets["criteria_definition"] = tk.Text(right, height=6, bg="#0F172A", fg="white", wrap="word", relief="flat", font=("Segoe UI", 9))
        self.student_widgets["criteria_definition"].pack(fill="x", padx=12, pady=(6, 12))

    def build_investigate_tab(self):
        self.directive(self.tab_investigate, "INVESTIGATE — Inspect source-linked candidate evidence. Missing evidence is a finding, not a reason to invent a score.")
        panel = tk.Frame(self.tab_investigate, bg=DARK["figure"])
        panel.pack(fill="both", expand=True, padx=10, pady=5)
        columns = ("Candidate", "HSP source", "Type", "T / method", "Gate", "Evidence status")
        self.investigate_tree = ttk.Treeview(panel, columns=columns, show="headings", height=14)
        widths = [150, 380, 150, 220, 120, 200]
        for col, width in zip(columns, widths):
            self.investigate_tree.heading(col, text=col)
            self.investigate_tree.column(col, width=width, anchor="w")
        self.investigate_tree.pack(fill="both", expand=True, side="left")
        scroll = ttk.Scrollbar(panel, orient="vertical", command=self.investigate_tree.yview)
        scroll.pack(side="right", fill="y")
        self.investigate_tree.configure(yscrollcommand=scroll.set)
        detail = tk.Frame(self.tab_investigate, bg=DARK["axes"], bd=1, relief="solid")
        detail.pack(fill="x", padx=10, pady=(5, 10))
        tk.Label(detail, text="Selected candidate evidence", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=(8, 3))
        self.investigate_detail = tk.Text(detail, height=8, bg="#0F172A", fg="white", wrap="word", relief="flat", font=("Segoe UI", 8))
        self.investigate_detail.pack(fill="x", padx=12, pady=(0, 10))
        self.investigate_tree.bind("<<TreeviewSelect>>", self.show_candidate_detail)

    def build_model_tab(self):
        self.directive(self.tab_model, "MODEL — Explore HSP space, domain-specific RED, and composition sensitivity. The model does not predict morphology, kinetics, evaporation, or compliance.")
        controls = tk.Frame(self.tab_model, bg=DARK["axes"], bd=1, relief="solid")
        controls.pack(fill="x", padx=10, pady=(0, 6))
        tk.Label(controls, text="Solvent A", bg=DARK["axes"], fg="white").pack(side="left", padx=(12, 4), pady=8)
        self.solvent_a_combo = ttk.Combobox(controls, state="readonly", width=18)
        self.solvent_a_combo.pack(side="left", padx=4)
        self.solvent_a_combo.bind("<<ComboboxSelected>>", lambda e: self.update_model())
        tk.Label(controls, text="Solvent B", bg=DARK["axes"], fg="white").pack(side="left", padx=(12, 4))
        self.solvent_b_combo = ttk.Combobox(controls, state="readonly", width=18)
        self.solvent_b_combo.pack(side="left", padx=4)
        self.solvent_b_combo.bind("<<ComboboxSelected>>", lambda e: self.update_model())
        tk.Label(controls, text="Volume % A", bg=DARK["axes"], fg="white").pack(side="left", padx=(12, 4))
        self.ratio_scale = tk.Scale(controls, from_=0, to=100, orient="horizontal", command=lambda v: self.update_model(),
                                    bg=DARK["axes"], fg="white", troughcolor="#0F172A", highlightthickness=0, length=210)
        self.ratio_scale.pack(side="left", padx=4)
        self.model_summary = ttk.Label(controls, text="", style="Info.TLabel")
        self.model_summary.pack(side="left", padx=12, pady=6, fill="x", expand=True)
        plot_frame = tk.Frame(self.tab_model, bg=DARK["figure"])
        plot_frame.pack(fill="both", expand=True, padx=10, pady=4)
        self.model_plot_host = plot_frame
        self.model_view_var = tk.StringVar(value="RED evidence")
        selector = tk.Frame(self.tab_model, bg=DARK["figure"])
        selector.pack(fill="x", padx=10, pady=(0, 10))
        for label in ("RED evidence", "HSP projections", "Mixture sensitivity"):
            tk.Radiobutton(selector, text=label, variable=self.model_view_var, value=label, command=self.render_model_plot,
                           bg=DARK["figure"], fg="white", selectcolor="#0F172A", activebackground=DARK["figure"], activeforeground="white").pack(side="left", padx=12)

    def build_compare_tab(self):
        self.directive(self.tab_compare, "COMPARE — Use transparent conditional indicators and hard gates. A high index is not proof of suitability; a HARD FAIL cannot be compensated by soft criteria.")
        weights = tk.Frame(self.tab_compare, bg=DARK["axes"], bd=1, relief="solid")
        weights.pack(fill="x", padx=10, pady=(0, 6))
        tk.Label(weights, text="Default/reproducible weights", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(side="left", padx=12, pady=8)
        self.weight_entries = []
        for label, value in zip(WEIGHT_LABELS, DEFAULT_WEIGHTS):
            box = tk.Frame(weights, bg=DARK["axes"])
            box.pack(side="left", padx=6)
            tk.Label(box, text=label.replace(" ", "\n"), bg=DARK["axes"], fg=DARK["subtext"], justify="center", font=("Segoe UI", 7)).pack()
            entry = tk.Spinbox(box, from_=0, to=100, width=5, bg="#0F172A", fg="white", insertbackground="white", command=self.apply_weights)
            entry.delete(0, "end")
            entry.insert(0, str(int(value)))
            entry.pack()
            self.weight_entries.append(entry)
        ttk.Button(weights, text="Apply weights", style="Primary.TButton", command=self.apply_weights).pack(side="left", padx=10)
        ttk.Button(weights, text="Reset baseline", style="Primary.TButton", command=self.reset_weights).pack(side="left", padx=4)
        self.sensitivity_label = ttk.Label(self.tab_compare, text="", style="Warn.TLabel", wraplength=1280)
        self.sensitivity_label.pack(fill="x", padx=10, pady=(0, 6))
        body = tk.Frame(self.tab_compare, bg=DARK["figure"])
        body.pack(fill="both", expand=True, padx=10, pady=4)
        cols = ("Candidate", "Avg RED", "Thermo", "Process", "EHS", "Sustain.", "Cost", "Index", "Gate", "Status")
        self.compare_tree = ttk.Treeview(body, columns=cols, show="headings", height=13)
        widths = [125, 75, 75, 75, 75, 75, 65, 75, 80, 255]
        for col, width in zip(cols, widths):
            self.compare_tree.heading(col, text=col)
            self.compare_tree.column(col, width=width, anchor="center")
        self.compare_tree.pack(fill="both", expand=True, side="left")
        scroll = ttk.Scrollbar(body, orient="vertical", command=self.compare_tree.yview)
        scroll.pack(side="right", fill="y")
        self.compare_tree.configure(yscrollcommand=scroll.set)
        note = tk.Label(self.tab_compare, text="Indicator note: Process, EHS, and sustainability values are explicit qualitative teaching mappings. Cost is missing by default and excluded from the conditional index rather than guessed.",
                        bg=DARK["figure"], fg=DARK["subtext"], wraplength=1280, justify="left", font=("Segoe UI", 8))
        note.pack(fill="x", padx=12, pady=(0, 10))

    def build_decide_tab(self):
        self.directive(self.tab_decide, "DECIDE — Select a candidate or shortlist, then make a defensible evidence-based recommendation. The software will warn about gates; it will not choose for you.")
        outer = tk.Frame(self.tab_decide, bg=DARK["figure"])
        outer.pack(fill="both", expand=True, padx=10, pady=4)
        left = tk.Frame(outer, bg=DARK["axes"], bd=1, relief="solid")
        left.pack(side="left", fill="both", expand=True, padx=(0, 5))
        right = tk.Frame(outer, bg=DARK["axes"], bd=1, relief="solid")
        right.pack(side="right", fill="both", expand=True, padx=(5, 0))
        tk.Label(left, text="Model–literature comparison prompt", bg=DARK["axes"], fg=DARK["amber"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=(10, 4))
        self.conflict_text = tk.Text(left, height=23, bg="#0F172A", fg="white", wrap="word", relief="flat", font=("Segoe UI", 9))
        self.conflict_text.pack(fill="both", expand=True, padx=12, pady=(0, 12))
        tk.Label(right, text="Defensible recommendation", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=(10, 4))
        tk.Label(right, text="Recommendation / shortlist", bg=DARK["axes"], fg=DARK["subtext"], font=("Segoe UI", 8)).pack(anchor="w", padx=12)
        self.recommendation_combo = ttk.Combobox(right, state="readonly")
        self.recommendation_combo.pack(fill="x", padx=12, pady=4)
        self.recommendation_combo.bind("<<ComboboxSelected>>", self.check_recommendation_gate)
        self.recommendation_status = ttk.Label(right, text="Select an option to inspect its documented gate.", style="Info.TLabel", wraplength=570)
        self.recommendation_status.pack(fill="x", padx=12, pady=4)
        for key, label, height in [
            ("computational_evidence", "Computational evidence (RED/domain/sensitivity)", 3),
            ("practical_evidence", "Literature/practical/EHS evidence", 3),
            ("main_tradeoff", "Main trade-off", 3),
        ]:
            tk.Label(right, text=label, bg=DARK["axes"], fg=DARK["subtext"], font=("Segoe UI", 8)).pack(anchor="w", padx=12, pady=(6, 1))
            widget = tk.Text(right, height=height, bg="#0F172A", fg="white", wrap="word", relief="flat", font=("Segoe UI", 8))
            widget.pack(fill="x", padx=12)
            self.student_widgets[key] = widget

    def build_reflect_tab(self):
        self.directive(self.tab_reflect, "REFLECT — State uncertainty, what would change your decision, and a verification experiment with a measurable endpoint.")
        panel = tk.Frame(self.tab_reflect, bg=DARK["axes"], bd=1, relief="solid")
        panel.pack(fill="both", expand=True, padx=10, pady=4)
        fields = [
            ("student_name", "Student name", 1),
            ("weight_observation", "What happened when you changed weights?", 3),
            ("conflict_analysis", "Where do model and literature/practical evidence agree or conflict?", 4),
            ("confidence", "Evidence-confidence level (Low / Medium / High) and why — not a probability", 2),
            ("main_uncertainty", "Main uncertainty or missing evidence", 3),
            ("change_decision", "What would change my decision?", 3),
            ("verification", "Verification experiment and measurable endpoint", 4),
        ]
        for key, label, height in fields:
            tk.Label(panel, text=label, bg=DARK["axes"], fg=DARK["subtext"], font=("Segoe UI", 8, "bold")).pack(anchor="w", padx=14, pady=(8, 1))
            if key == "student_name":
                widget = tk.Entry(panel, bg="#0F172A", fg="white", insertbackground="white", relief="flat", font=("Segoe UI", 9))
                widget.pack(fill="x", padx=14)
            else:
                widget = tk.Text(panel, height=height, bg="#0F172A", fg="white", wrap="word", relief="flat", font=("Segoe UI", 8))
                widget.pack(fill="x", padx=14)
            self.student_widgets[key] = widget
        tk.Label(panel, text="Required student-facing sentence: I recommend ___ because the computational evidence indicates ___, while the literature/practical evidence indicates ___. The main uncertainty is ___. I would change this decision if ___. The next verification step is ___.",
                 bg=DARK["axes"], fg=DARK["amber"], wraplength=1240, justify="left", font=("Segoe UI", 9, "bold")).pack(fill="x", padx=14, pady=(10, 12))

    def get_student_field(self, key: str) -> str:
        widget = self.student_widgets.get(key)
        if widget is None:
            return ""
        if isinstance(widget, tk.Entry):
            return widget.get().strip()
        return widget.get("1.0", "end").strip()

    def change_case(self):
        self.current_case_id = self.case_var.get()
        self.case_data = CASES_DATABASE[self.current_case_id]
        self.solvent_a, self.solvent_b, self.ratio_a = self.case_data["default_blend"]
        self.refresh_all()

    def refresh_all(self):
        case = self.case_data
        self.define_context.delete("1.0", "end")
        self.define_context.insert("end", f"Engineering problem:\n{case['problem_statement']}\n\nApplication:\n{case['application']}\n\nGuidance level:\n{case['student_guidance']}\n\nStudent task:\n{case['criteria_prompt']}\n\nPolymer specification:\n{case['polymer']['material_specification']}")
        names = list(case["solvents"].keys())
        self.solvent_a_combo["values"] = names
        self.solvent_b_combo["values"] = names
        self.solvent_a_combo.set(self.solvent_a)
        self.solvent_b_combo.set(self.solvent_b)
        self.ratio_scale.set(self.ratio_a)
        self.recommendation_combo["values"] = names + [f"Custom blend: {self.solvent_a}/{self.solvent_b}", "Shortlist / insufficient evidence"]
        self.recommendation_combo.set("Shortlist / insufficient evidence")
        self.populate_investigate()
        self.populate_compare()
        self.populate_conflict()
        self.update_model()
        self.check_recommendation_gate()

    def populate_investigate(self):
        self.investigate_tree.delete(*self.investigate_tree.get_children())
        for name, record in self.case_data["solvents"].items():
            src = record["source"]
            gate, _ = hard_gate_status(record)
            evidence = evidence_status(record)
            self.investigate_tree.insert("", "end", iid=name, values=(name, src["citation"], src["value_type"], f"{src['temperature']}; {src['method']}", gate, evidence))
        self.investigate_detail.delete("1.0", "end")
        self.investigate_detail.insert("end", "Select a candidate row to inspect its source, evidence fields, and gate rationale.")

    def show_candidate_detail(self, event=None):
        selected = self.investigate_tree.selection()
        if not selected:
            return
        name = selected[0]
        record = self.case_data["solvents"][name]
        src = record["source"]
        gate = record["hard_gate"]
        h = record["hsp"]
        content = (
            f"CANDIDATE: {name}\n\n"
            f"HSP: δD={h[0]:.2f}, δP={h[1]:.2f}, δH={h[2]:.2f} {HSP_UNITS}\n"
            f"Boiling point record: {record['bp_c']} °C\n"
            f"HSP/source citation: {src['citation']}\n"
            f"Value type: {src['value_type']}\nTemperature: {src['temperature']}\nMethod: {src['method']}\nUncertainty: {src['uncertainty']}\n\n"
            f"PROCESS EVIDENCE\n{record['process_evidence']}\n\n"
            f"EHS / REGULATORY EVIDENCE\n{record['ehs_evidence']}\n\n"
            f"SUSTAINABILITY INDICATOR EVIDENCE\n{record['sustainability_evidence']}\n\n"
            f"COST / AVAILABILITY EVIDENCE\n{record['cost_evidence']}\n\n"
            f"LITERATURE / EMPIRICAL EVIDENCE\n{record['literature_evidence']}\n\n"
            f"HARD GATE: {gate['status']}\nRule: {gate['rule']}\nReason: {gate['reason']}\nSource: {gate['source']}"
        )
        self.investigate_detail.delete("1.0", "end")
        self.investigate_detail.insert("end", content)

    def update_model(self):
        if not hasattr(self, "solvent_a_combo"):
            return
        self.solvent_a = self.solvent_a_combo.get() or self.solvent_a
        self.solvent_b = self.solvent_b_combo.get() or self.solvent_b
        self.ratio_a = float(self.ratio_scale.get())
        ratio_b = 100.0 - self.ratio_a
        if self.solvent_a == self.solvent_b:
            self.model_summary.config(text="Choose two different solvents to make a binary-mixture sensitivity plot.", style="Warn.TLabel")
            return
        mix = calc_mixed_hsp(self.case_data["solvents"], [self.solvent_a, self.solvent_b], [self.ratio_a, ratio_b])
        domains = evaluate_domains(self.case_data["polymer"], mix)
        text = "; ".join(f"{name} RED={row['red']:.3f}" for name, row in domains.items())
        self.model_summary.config(text=f"Mixture: {self.ratio_a:.0f}% {self.solvent_a} / {ratio_b:.0f}% {self.solvent_b} | {text}", style="Info.TLabel")
        self.render_model_plot()

    def render_model_plot(self):
        if not hasattr(self, "model_plot_host") or self.solvent_a == self.solvent_b:
            return
        for widget in self.model_plot_host.winfo_children():
            widget.destroy()
        ratio_b = 100.0 - self.ratio_a
        mix = calc_mixed_hsp(self.case_data["solvents"], [self.solvent_a, self.solvent_b], [self.ratio_a, ratio_b])
        label = f"{self.ratio_a:.0f}% {self.solvent_a} / {ratio_b:.0f}% {self.solvent_b}"
        view = self.model_view_var.get()
        if view == "HSP projections":
            figure = plot_hsp_projections(self.case_data, mix, label)
        elif view == "Mixture sensitivity":
            figure = plot_mixture_sensitivity(self.case_data, self.solvent_a, self.solvent_b, self.ratio_a)
        else:
            figure = plot_red_evidence(self.case_data, mix, label)
        canvas = FigureCanvasTkAgg(figure, master=self.model_plot_host)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        plt.close(figure)

    def read_weights(self):
        try:
            trial = [float(entry.get()) for entry in self.weight_entries]
            if any(x < 0 for x in trial) or sum(trial) <= 0:
                raise ValueError
            return trial
        except Exception:
            return DEFAULT_WEIGHTS[:]

    def apply_weights(self):
        self.weights = self.read_weights()
        self.populate_compare()

    def reset_weights(self):
        for entry, value in zip(self.weight_entries, DEFAULT_WEIGHTS):
            entry.delete(0, "end")
            entry.insert(0, str(int(value)))
        self.weights = DEFAULT_WEIGHTS[:]
        self.populate_compare()

    def populate_compare(self):
        self.compare_tree.delete(*self.compare_tree.get_children())
        polymer = self.case_data["polymer"]
        for name, record in self.case_data["solvents"].items():
            profile = compute_evidence_profile(record, evaluate_domains(polymer, record["hsp"]), self.weights)
            score_text = ["—" if value is None else f"{value:.0f}" for value in profile["scores"]]
            index = "—" if profile["conditional_index"] is None else f"{profile['conditional_index']:.1f}"
            self.compare_tree.insert("", "end", values=(name, f"{profile['average_red']:.2f}", *score_text, index, profile["hard_gate"], profile["overall_status"]))
        sensitivity = weight_sensitivity(self.case_data, self.weights)
        style = "Warn.TLabel" if sensitivity["changed"] else "Info.TLabel"
        self.sensitivity_label.config(text="Weight sensitivity: " + sensitivity["message"], style=style)

    def populate_conflict(self):
        conflict = self.case_data["conflict_case"]
        content = (
            f"FOCUS SYSTEM\n{conflict['system']}\n\n"
            f"MODEL EVIDENCE\n{conflict['model_prediction']}\n\n"
            f"REPORTED / LITERATURE EVIDENCE\n{conflict['reported_observation']}\n\n"
            f"SOURCE(S)\n{conflict['evidence_source']}\n\n"
            f"STUDENT TASK\n{conflict['student_prompt']}\n\n"
            f"WHAT THE MODEL DOES NOT ESTABLISH\n{conflict['limitations']}"
        )
        self.conflict_text.delete("1.0", "end")
        self.conflict_text.insert("end", content)

    def check_recommendation_gate(self, event=None):
        choice = self.recommendation_combo.get() if hasattr(self, "recommendation_combo") else ""
        if choice in self.case_data["solvents"]:
            gate, reason = hard_gate_status(self.case_data["solvents"][choice])
            if gate == "FAIL":
                self.recommendation_status.config(text=f"HARD FAIL: {choice} violates a documented mandatory gate. A high conditional index cannot override this. Reason: {reason}", style="Fail.TLabel")
            elif gate == "WARN":
                self.recommendation_status.config(text=f"WARN: {choice} requires evidence/trade-off review. {reason}", style="Warn.TLabel")
            else:
                self.recommendation_status.config(text=f"PASS on the available documented gate only. This does not prove experimental success. {reason}", style="Good.TLabel")
        else:
            self.recommendation_status.config(text="No candidate is selected as automatically best. Explain why a shortlist or custom mixture is justified.", style="Info.TLabel")

    def collect_student_fields(self):
        keys = [
            "student_name", "problem_definition", "criteria_definition", "weight_observation", "conflict_analysis",
            "computational_evidence", "practical_evidence", "main_tradeoff", "confidence", "main_uncertainty",
            "change_decision", "verification",
        ]
        values = {key: self.get_student_field(key) for key in keys}
        values["recommendation"] = self.recommendation_combo.get()
        return values

    def export_report(self):
        ratio_b = 100.0 - self.ratio_a
        filename = f"Material_Selection_Brief_{self.current_case_id}.md"
        export_material_selection_brief(self.case_data, (self.solvent_a, self.solvent_b, self.ratio_a, ratio_b), self.weights, self.collect_student_fields(), filename)
        if messagebox:
            messagebox.showinfo("Export complete", f"Saved: {os.path.abspath(filename)}")

    def export_tables(self):
        audit_name = f"Decision_Engine_Audit_{self.current_case_id}.csv"
        evidence_name = f"Candidate_Evidence_{self.current_case_id}.csv"
        export_audit_table(self.case_data, audit_name)
        export_candidate_evidence_csv(self.case_data, evidence_name)
        if messagebox:
            messagebox.showinfo("Export complete", f"Saved:\n{os.path.abspath(audit_name)}\n{os.path.abspath(evidence_name)}")


# -----------------------------------------------------------------------------
# Verification and command-line entry point
# -----------------------------------------------------------------------------

def run_self_tests():
    print("=" * 72)
    print("POLYMERSOLV STUDIO v2 — SELF-TESTS")
    print("=" * 72)
    naf = CASES_DATABASE["nafion"]
    mixture = calc_mixed_hsp(naf["solvents"], ["Isopropanol", "Water"], [80, 20])
    red_values = evaluate_domains(naf["polymer"], mixture)
    assert abs(red_values["Backbone"]["red"] - 1.475) < 0.05
    assert abs(red_values["Side-chain"]["red"] - 0.886) < 0.05
    print("PASS 1: Nafion domain-specific mixture RED calculations reproduce reference teaching values.")
    ca = CASES_DATABASE["cellulose_acetate"]
    red_nmp = calc_red(ca["solvents"]["NMP"]["hsp"], ca["polymer"]["domains"]["Polymer"]["hsp"], ca["polymer"]["domains"]["Polymer"]["R0"])
    assert abs(red_nmp - 0.527) < 0.05
    print("PASS 2: Cellulose-acetate/NMP RED reference calculation is internally consistent.")
    ps = CASES_DATABASE["polystyrene"]
    profile = compute_evidence_profile(ps["solvents"]["Toluene"], evaluate_domains(ps["polymer"], ps["solvents"]["Toluene"]["hsp"]), DEFAULT_WEIGHTS)
    assert profile["hard_gate"] == "WARN"
    assert profile["conditional_index"] is not None
    print("PASS 3: Conditional index reports warning-level evidence rather than an automatic recommendation.")
    for cid, case in CASES_DATABASE.items():
        audit = audit_rows(case)
        assert len(audit) >= 9
        assert all("component" in row for row in audit)
    print("PASS 4: Decision-engine audit records are available for every case.")
    print("All self-tests passed.")


def main():
    parser = argparse.ArgumentParser(description="PolymerSolv Studio — evidence-based polymer–solvent material selection")
    parser.add_argument("--test", action="store_true", help="Run calculation and audit self-tests without a GUI")
    parser.add_argument("--export-all", action="store_true", help="Export figures, evidence tables, audit tables, and brief templates for all cases")
    args = parser.parse_args()
    if args.test:
        run_self_tests()
        return
    if args.export_all:
        for case in CASES_DATABASE.values():
            paths = export_case_bundle(case)
            print(f"Exported {case['id']}:\n  " + "\n  ".join(paths))
        return
    if not TK_AVAILABLE:
        print("Tkinter is not available in this Python installation. Run with --test or --export-all, or install/use a Python build that includes Tkinter.")
        return
    root = tk.Tk()
    PolymerSolvPages(root)
    root.mainloop()



# -----------------------------------------------------------------------------
# Six-stage page update. All original Model functions and Model GUI methods above
# are inherited, not altered. The entry fields below are student-editable.
# -----------------------------------------------------------------------------
import copy
import re

# USER-EDITABLE: extend polymer and solvent choices here. The solvent list is a
# polarity catalog, NOT evidence that a solvent dissolves any polymer.
POLYMER_CHOICES = ("Nafion", "Cellulose acetate", "Polystyrene (PS)", "PVDF", "PVA", "PES", "PAN", "PMMA", "Other / Custom")
SOLVENT_POLARITY_CATALOG = {
    "DMAc": ("Polar", "Aprotic"), "DMF": ("Polar", "Aprotic"), "DMSO": ("Polar", "Aprotic"),
    "NMP": ("Polar", "Aprotic"), "Acetone": ("Polar", "Aprotic"), "MEK": ("Polar", "Aprotic"),
    "Ethyl_Lactate": ("Polar", "Protic"), "Ethanol": ("Polar", "Protic"),
    "Isopropanol": ("Polar", "Protic"), "Water": ("Polar", "Protic"),
    "Ethyl_Acetate": ("Polar", "Aprotic"), "2_MeTHF": ("Polar", "Aprotic"),
    "DCM": ("Polar", "Aprotic"), "Toluene": ("Nonpolar", "Aprotic"),
    "d_Limonene": ("Nonpolar", "Aprotic"), "Cyclohexane": ("Nonpolar", "Aprotic"),
}
LITERATURE_FIELDS = (
    ("doi", "DOI or source URL"), ("title", "Title of paper"),
    ("authors", "Authors"), ("year", "Publication year"),
    ("polymer", "Polymer in paper"), ("solvent_a", "Solvent A"),
    ("solvent_b", "Solvent B (optional)"), ("ratio", "Reported ratio, with basis, e.g. 70/30 v/v"),
    ("concentration", "Polymer concentration with units"),
    ("temperature", "Temperature with units"), ("time", "Dissolution time with units"),
    ("observation", "Paper's reported observation"), ("notes", "Student notes / caveats"),
)


def _canonical_solvent(name):
    aliases = {"ipa": "Isopropanol", "isopropyl alcohol": "Isopropanol",
        "ethyl acetate": "Ethyl_Acetate", "ethyl lactate": "Ethyl_Lactate",
        "d-limonene": "d_Limonene", "h2o": "Water"}
    key = name.strip().lower().replace("_", " ")
    return aliases.get(key, next((x for x in SOLVENT_POLARITY_CATALOG if x.lower().replace("_", " ") == key), None))


def _volume_ratio(text):
    """Only explicitly identified volume/volume percentages may be subtracted."""
    if not re.search(r"(?:v/v|vol/vol|vol%)", text, re.I):
        return None
    match = re.search(r"(\d+(?:\.\d+)?)\s*[/ :]\s*(\d+(?:\.\d+)?)", text)
    if not match:
        return None
    a, b = float(match.group(1)), float(match.group(2))
    return (a, b) if 0 <= a <= 100 and 0 <= b <= 100 and abs(a + b - 100) < 1e-6 else None


class PolymerSolvPages(PolymerSolvGUI):
    """Retains the original model, its widgets, calculations, plots, and exports."""
    def __init__(self, root):
        self.literature_records = []
        self.no_literature = False
        self.response_state = {}
        super().__init__(root)

    def build_tabs(self):
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=10, pady=10)
        self.tab_define = ttk.Frame(self.notebook)
        self.tab_investigate = ttk.Frame(self.notebook)
        self.tab_model = ttk.Frame(self.notebook)
        self.tab_compare = ttk.Frame(self.notebook)
        self.tab_decide = ttk.Frame(self.notebook)
        self.tab_reflect = ttk.Frame(self.notebook)
        for tab, title in ((self.tab_define, "1. DEFINE"), (self.tab_investigate, "2. INVESTIGATE"),
                           (self.tab_model, "3. MODEL"), (self.tab_compare, "4. COMPARE"),
                           (self.tab_decide, "5. DECIDE"), (self.tab_reflect, "6. REFLECT")):
            self.notebook.add(tab, text=title)
        self.build_define_tab()
        self.build_investigate_tab()
        self.build_model_tab()  # original protected Model page, unchanged
        self.build_compare_tab()
        self.build_decide_tab()
        self.build_reflect_tab()
        self._add_navigation()
        self.notebook.bind("<<NotebookTabChanged>>", self._tab_changed)

    def build_define_tab(self):
        self.directive(self.tab_define, "DEFINE — State your own problem; catalog suggestions are NOT literature evidence.")
        container = tk.Frame(self.tab_define, bg=DARK["figure"])
        container.pack(fill="both", expand=True, padx=10, pady=5)
        left = tk.Frame(container, bg=DARK["axes"], bd=1, relief="solid")
        left.pack(side="left", fill="both", expand=True, padx=(0, 5))
        right = tk.Frame(container, bg=DARK["axes"], bd=1, relief="solid")
        right.pack(side="right", fill="both", expand=True, padx=(5, 0))
        tk.Label(left, text="Teaching case for the existing Model", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=8)
        self.case_var = tk.StringVar(value="nafion")
        for cid, case in CASES_DATABASE.items():
            tk.Radiobutton(left, text=case["name"], variable=self.case_var, value=cid,
                command=self.change_case, bg=DARK["axes"], fg="white", selectcolor="#0F172A",
                activebackground=DARK["axes"], activeforeground="white").pack(anchor="w", padx=10)
        self.define_context = tk.Text(left, bg="#0F172A", fg="white", wrap="word", relief="flat", font=("Segoe UI", 9))
        self.define_context.pack(fill="both", expand=True, padx=12, pady=10)
        tk.Label(right, text="Your problem definition", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=10)
        self.student_widgets["problem_definition"] = tk.Text(right, height=7, bg="#0F172A", fg="white", insertbackground="white", wrap="word", relief="flat", font=("Segoe UI", 10))
        self.student_widgets["problem_definition"].pack(fill="x", padx=12, pady=(0, 12))
        tk.Label(right, text="Polymer type (USER-EDITABLE catalog above)", bg=DARK["axes"], fg=DARK["subtext"]).pack(anchor="w", padx=12)
        self.polymer_var = tk.StringVar(value="Nafion")
        widget = ttk.Combobox(right, textvariable=self.polymer_var, values=POLYMER_CHOICES, state="readonly")
        widget.pack(fill="x", padx=12, pady=5)
        widget.bind("<<ComboboxSelected>>", self._polymer_changed)
        self.custom_var = tk.StringVar()
        self.custom_line = tk.Frame(right, bg=DARK["axes"])
        tk.Label(self.custom_line, text="Custom polymer name", bg=DARK["axes"], fg="white").pack(anchor="w")
        tk.Entry(self.custom_line, textvariable=self.custom_var, bg="#0F172A", fg="white", insertbackground="white").pack(fill="x")
        tk.Label(right, text="Solvent polarity", bg=DARK["axes"], fg=DARK["subtext"]).pack(anchor="w", padx=12)
        self.polarity_var = tk.StringVar(value="Any / Not specified")
        widget = ttk.Combobox(right, textvariable=self.polarity_var, values=("Polar", "Nonpolar", "Any / Not specified"), state="readonly")
        widget.pack(fill="x", padx=12, pady=5); widget.bind("<<ComboboxSelected>>", self._update_suggestions)
        tk.Label(right, text="Solvent class", bg=DARK["axes"], fg=DARK["subtext"]).pack(anchor="w", padx=12)
        self.class_var = tk.StringVar(value="Any / Not specified")
        widget = ttk.Combobox(right, textvariable=self.class_var, values=("Protic", "Aprotic", "Any / Not specified"), state="readonly")
        widget.pack(fill="x", padx=12, pady=5); widget.bind("<<ComboboxSelected>>", self._update_suggestions)
        tk.Label(right, text="Suggested solvents (polarity matches, NOT validated dissolving solvents)", bg=DARK["axes"], fg=DARK["subtext"], wraplength=500).pack(anchor="w", padx=12)
        self.suggest_var = tk.StringVar()
        self.suggest_combo = ttk.Combobox(right, textvariable=self.suggest_var, state="readonly")
        self.suggest_combo.pack(fill="x", padx=12, pady=6)
        self.define_note = tk.Label(right, text="", bg=DARK["axes"], fg=DARK["amber"], wraplength=500, justify="left")
        self.define_note.pack(fill="x", padx=12, pady=6)
        self._update_suggestions()

    def _polymer_changed(self, event=None):
        if self.polymer_var.get() == "Other / Custom":
            self.custom_line.pack(fill="x", padx=12, pady=6)
        else:
            self.custom_line.pack_forget()
        self._update_suggestions()

    def _update_suggestions(self, event=None):
        values = [name for name, (polarity, kind) in SOLVENT_POLARITY_CATALOG.items()
                  if (self.polarity_var.get() == "Any / Not specified" or polarity == self.polarity_var.get())
                  and (self.class_var.get() == "Any / Not specified" or kind == self.class_var.get())]
        self.suggest_combo["values"] = values
        if self.suggest_var.get() not in values:
            self.suggest_var.set(values[0] if values else "")
        self.define_note.config(text=f"{len(values)} catalog match(es). The Model still uses {self.case_data['polymer']['name']}; if you selected another polymer, its RED is not a prediction for your polymer.")

    def build_investigate_tab(self):
        self.directive(self.tab_investigate, "INVESTIGATE — Search scientific papers and enter the paper's actual information; no automated retrieval or invented citations.")
        shell = tk.Frame(self.tab_investigate, bg=DARK["figure"])
        shell.pack(fill="both", expand=True, padx=10, pady=5)
        left = tk.Frame(shell, bg=DARK["axes"]); left.pack(side="left", fill="both", expand=True, padx=(0, 5))
        right = tk.Frame(shell, bg=DARK["axes"]); right.pack(side="right", fill="both", expand=True, padx=(5, 0))
        self.no_lit_var = tk.BooleanVar(value=False)
        tk.Checkbutton(left, text="No suitable literature data found / could not find relevant data",
            variable=self.no_lit_var, command=self._toggle_no_lit, bg=DARK["axes"], fg=DARK["amber"],
            selectcolor="#0F172A", activebackground=DARK["axes"], activeforeground=DARK["amber"]).pack(anchor="w", padx=10, pady=5)
        form_canvas=tk.Canvas(left,bg=DARK["axes"],highlightthickness=0)
        scroll=ttk.Scrollbar(left,orient="vertical",command=form_canvas.yview)
        form_canvas.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right",fill="y");form_canvas.pack(side="left",fill="both",expand=True)
        form=tk.Frame(form_canvas,bg=DARK["axes"])
        win=form_canvas.create_window((0,0),window=form,anchor="nw")
        form.bind("<Configure>",lambda e:form_canvas.configure(scrollregion=form_canvas.bbox("all")))
        form_canvas.bind("<Configure>",lambda e:form_canvas.itemconfigure(win,width=e.width))
        self.lit_fields={}
        solvent_values=tuple(sorted({n for c in CASES_DATABASE.values() for n in c["solvents"]}))
        for key,label in LITERATURE_FIELDS:
            tk.Label(form,text=label,bg=DARK["axes"],fg=DARK["subtext"]).pack(anchor="w",padx=9,pady=(5,1))
            if key=="polymer": w=ttk.Combobox(form,values=POLYMER_CHOICES,state="normal")
            elif key in ("solvent_a","solvent_b"): w=ttk.Combobox(form,values=solvent_values,state="normal")
            elif key in ("observation","notes"): w=tk.Text(form,height=3,bg="#0F172A",fg="white",insertbackground="white",wrap="word")
            else: w=tk.Entry(form,bg="#0F172A",fg="white",insertbackground="white")
            w.pack(fill="x",padx=9,pady=2);self.lit_fields[key]=w
        bar=tk.Frame(form,bg=DARK["axes"]);bar.pack(fill="x",padx=9,pady=8)
        ttk.Button(bar,text="Add / update paper",command=self._save_paper).pack(side="left",padx=3)
        ttk.Button(bar,text="Clear form",command=self._clear_paper_form).pack(side="left",padx=3)
        self.lit_tree=ttk.Treeview(right,columns=("Title","Polymer","Solvents","Ratio","DOI"),show="headings",height=11)
        for col in ("Title","Polymer","Solvents","Ratio","DOI"): self.lit_tree.heading(col,text=col);self.lit_tree.column(col,width=130)
        self.lit_tree.pack(fill="both",expand=True,padx=9,pady=10)
        self.lit_tree.bind("<<TreeviewSelect>>",self._load_paper)
        self.lit_status=tk.Label(right,bg=DARK["axes"],fg=DARK["amber"],text="No student-entered papers.")
        self.lit_status.pack(anchor="w",padx=9)
        ttk.Button(right,text="Delete selected paper",command=self._delete_paper).pack(anchor="w",padx=9,pady=7)
        self.edit_paper_index=None
        self._refresh_papers()

    def _field(self,key):
        w=self.lit_fields[key]
        return (w.get("1.0","end-1c") if isinstance(w,tk.Text) else w.get()).strip()

    def _set_field(self,key,value):
        w=self.lit_fields[key]
        if isinstance(w,tk.Text):w.delete("1.0","end");w.insert("1.0",value)
        else:w.delete(0,"end");w.insert(0,value)

    def _clear_paper_form(self):
        self.edit_paper_index=None
        for key in self.lit_fields:self._set_field(key,"")
        self.lit_tree.selection_remove(self.lit_tree.selection())

    def _toggle_no_lit(self):
        if self.no_lit_var.get() and self.literature_records:
            self.no_lit_var.set(False)
            messagebox.showwarning("Existing papers", "Delete recorded papers before marking literature unavailable.")
        self.no_literature=self.no_lit_var.get()
        self._refresh_papers()

    def _save_paper(self):
        r={key:self._field(key) for key,_ in LITERATURE_FIELDS}
        if self.no_lit_var.get():messagebox.showwarning("No data", "Uncheck the no-data option before entering a paper.");return
        if not r["title"] or not r["polymer"] or not r["solvent_a"] or not (r["doi"] or r["authors"]):
            messagebox.showwarning("Incomplete", "Enter paper title, polymer, solvent A, and DOI/URL or authors.");return
        if r["year"] and not re.fullmatch(r"\d{4}",r["year"]):
            messagebox.showwarning("Year", "Enter a four-digit year or leave blank.");return
        if r["ratio"] and not r["solvent_b"]:
            messagebox.showwarning("Ratio", "Specify solvent B for a mixture ratio or leave ratio blank.");return
        if r["ratio"] and not re.search(r"(?:v/v|w/w|mol/mol|vol%)",r["ratio"],re.I):
            messagebox.showwarning("Basis", "Write a ratio basis, e.g. 70/30 v/v; otherwise leave blank and describe in notes.");return
        if self.edit_paper_index is None:self.literature_records.append(r)
        else:self.literature_records[self.edit_paper_index]=r
        self._clear_paper_form();self._refresh_papers()

    def _load_paper(self,event=None):
        selected=self.lit_tree.selection()
        if not selected:return
        self.edit_paper_index=int(selected[0])
        for key,value in self.literature_records[self.edit_paper_index].items():self._set_field(key,value)

    def _delete_paper(self):
        selected=self.lit_tree.selection()
        if selected:del self.literature_records[int(selected[0])]
        self._clear_paper_form();self._refresh_papers()

    def _refresh_papers(self):
        self.lit_tree.delete(*self.lit_tree.get_children())
        for i,r in enumerate(self.literature_records):
            system=r["solvent_a"]+(" / "+r["solvent_b"] if r["solvent_b"] else "")
            self.lit_tree.insert("","end",iid=str(i),values=(r["title"],r["polymer"],system,r["ratio"] or "N/A",r["doi"] or "N/A"))
        self.lit_status.config(text="No literature found (student selected)." if self.no_lit_var.get() else f"{len(self.literature_records)} student-entered paper(s).")

    def populate_investigate(self):
        self._refresh_papers()  # Don't delete records when switching case.

    def build_compare_tab(self):
        self.directive(self.tab_compare, "COMPARE — Compare student-entered papers with current Model inputs and outputs; non-comparable values are N/A.")
        wrapper=tk.Frame(self.tab_compare,bg=DARK["axes"]);wrapper.pack(fill="x",padx=10,pady=6)
        cols=("Paper","Paper polymer","Paper solvents","Paper ratio","Model solvents","Model ratio","Model RED","Comparison")
        self.lit_compare_tree=ttk.Treeview(wrapper,columns=cols,show="headings",height=7)
        for col,width in zip(cols,(140,110,160,115,150,100,165,350)):
            self.lit_compare_tree.heading(col,text=col);self.lit_compare_tree.column(col,width=width)
        self.lit_compare_tree.pack(fill="x",padx=8,pady=8)
        self.comparison_status=tk.Label(wrapper,text="",bg=DARK["axes"],fg=DARK["amber"],wraplength=1250,justify="left")
        self.comparison_status.pack(fill="x",padx=8,pady=5)
        ttk.Button(wrapper,text="Refresh literature comparison",command=self._refresh_comparison).pack(anchor="e",padx=8,pady=5)
        # Preserve the former Compare criteria matrix beneath the new literature comparison.
        self._build_existing_compare_panel()

    def _build_existing_compare_panel(self):
        # Only the parent directive changes; calculations and table are inherited.
        super().build_compare_tab()

    def _model_snapshot(self):
        mixed=calc_mixed_hsp(self.case_data["solvents"],[self.solvent_a,self.solvent_b],[self.ratio_a,100-self.ratio_a])
        domains=evaluate_domains(self.case_data["polymer"],mixed)
        red=", ".join(f"{name}: {data['red']:.3f}" for name,data in domains.items())
        return {"system":(self.solvent_a,self.solvent_b),"ratio_a":self.ratio_a,
            "red":red,"hsp":mixed,"domains":domains,"model_polymer":self.case_data["polymer"]["name"]}

    def _comparison_rows(self):
        if self.no_lit_var.get():return []
        m=self._model_snapshot(); rows=[]
        student=self.polymer_var.get()
        selected=("Nafion" if self.current_case_id=="nafion" else "Cellulose acetate" if self.current_case_id=="cellulose_acetate" else "Polystyrene (PS)")
        model_supported=student==selected
        for rec in self.literature_records:
            sa,sb=_canonical_solvent(rec["solvent_a"]),_canonical_solvent(rec["solvent_b"]) if rec["solvent_b"] else None
            model_a,model_b=m["system"]
            same_poly=rec["polymer"].strip().lower()==student.strip().lower() and model_supported
            pure_model=(model_a if self.ratio_a==100 else model_b if self.ratio_a==0 else None)
            same_set=(same_poly and {sa,sb}=={model_a,model_b}) if sb else (same_poly and sa==pure_model)
            if not model_supported:note="N/A — Model calibrated for a different polymer"
            elif rec["polymer"].strip().lower()!=student.strip().lower():note="N/A — paper investigates a different polymer"
            elif not same_set:note="Different solvent system; no direct comparison"
            else:
                fractions=_volume_ratio(rec["ratio"])
                if fractions and sb:
                    fraction_model=m["ratio_a"] if sa==model_a else 100-m["ratio_a"]
                    note=f"Same system; volume ratio difference {abs(fractions[0]-fraction_model):.1f} percentage points (solvent A basis)"
                elif not sb:note="Pure-solvent record; confirm model is set to this pure solvent"
                else:note="Same system; ratio basis missing or different; N/A"
            rows.append({"paper":rec["title"],"polymer":rec["polymer"],
                "paper_system":rec["solvent_a"]+(" / "+rec["solvent_b"] if rec["solvent_b"] else ""),
                "paper_ratio":rec["ratio"] or "N/A","model_system":f"{model_a} / {model_b}",
                "model_ratio":f"{m['ratio_a']:.0f}/{100-m['ratio_a']:.0f} v/v",
                "red":m["red"] if model_supported else "N/A", "comparison":note})
        return rows

    def _refresh_comparison(self):
        self.lit_compare_tree.delete(*self.lit_compare_tree.get_children())
        for r in self._comparison_rows():
            self.lit_compare_tree.insert("","end",values=tuple(r[k] for k in
                ("paper","polymer","paper_system","paper_ratio","model_system","model_ratio","red","comparison")))
        self.comparison_status.config(text=("No suitable literature data found; model results are separate." if self.no_lit_var.get() else
            f"Compared {len(self.literature_records)} student-entered paper(s). RED is model output, not experimental dissolution." if self.literature_records else
            "No paper records to compare. Add sources in INVESTIGATE or select the no-data option."))

    def _add_navigation(self):
        pages=(self.tab_define,self.tab_investigate,self.tab_model,self.tab_compare,self.tab_decide,self.tab_reflect)
        for i,page in enumerate(pages):
            bar=tk.Frame(page,bg=DARK["figure"])
            bar.pack(side="bottom",fill="x",padx=10,pady=4)
            if i:ttk.Button(bar,text="← Back",command=lambda:self._step(-1)).pack(side="left")
            if i<5:ttk.Button(bar,text="Next →",style="Primary.TButton",command=lambda:self._step(1)).pack(side="right")

    def _save_responses(self):
        self.response_state={"problem_definition":self.get_student_field("problem_definition"),
            "polymer":self.polymer_var.get(),"polarity":self.polarity_var.get(),
            "class":self.class_var.get(),"suggestion":self.suggest_var.get(),
            "literature_records":copy.deepcopy(self.literature_records),"no_literature":self.no_lit_var.get(),
            "model":self._model_snapshot(),"comparison":self._comparison_rows(),
            "recommendation":self.recommendation_combo.get()}

    def _step(self,direction):
        self._save_responses()
        pages=(self.tab_define,self.tab_investigate,self.tab_model,self.tab_compare,self.tab_decide,self.tab_reflect)
        current=pages.index(self.notebook.nametowidget(self.notebook.select()))
        next_page=pages[max(0,min(5,current+direction))]
        self.notebook.select(next_page)
        if next_page==self.tab_compare:self._refresh_comparison()

    def _tab_changed(self,event=None):
        if self.notebook.select()==str(self.tab_compare):self._refresh_comparison()

    def refresh_all(self):
        case=self.case_data
        self.define_context.delete("1.0","end")
        self.define_context.insert("end",f"Engineering problem:\n{case['problem_statement']}\n\nApplication:\n{case['application']}\n\nPolymer specification:\n{case['polymer']['material_specification']}")
        names=list(case["solvents"])
        self.solvent_a_combo["values"]=names;self.solvent_b_combo["values"]=names
        self.solvent_a_combo.set(self.solvent_a);self.solvent_b_combo.set(self.solvent_b)
        self.ratio_scale.set(self.ratio_a)
        self.recommendation_combo["values"]=names+[f"Custom blend: {self.solvent_a}/{self.solvent_b}","Shortlist / insufficient evidence"]
        # The legacy weighted comparison table was removed; refresh the active
        # student literature-versus-Model table instead.
        # Do not call legacy populate_compare()/populate_conflict(): their widgets
        # were removed by the redesigned Compare and DECIDE pages.
        self.populate_investigate()
        self.update_model()
        self.check_recommendation_gate()
        self._update_suggestions()
        self._refresh_comparison()

    def collect_student_fields(self):
        values=super().collect_student_fields()
        self._save_responses()
        values.update({"define_polymer":self.polymer_var.get(),"define_polarity":self.polarity_var.get(),
            "define_class":self.class_var.get(),"define_suggested":self.suggest_var.get(),
            "literature_records":copy.deepcopy(self.literature_records),"no_literature":self.no_lit_var.get(),
            "literature_comparison":self._comparison_rows()})
        return values


# UI revision: the legacy weighted conditional-index panel is deliberately omitted.
# Compare now contains only student literature-to-current-Model comparison.
class PolymerSolvPages(PolymerSolvPages):
    def build_compare_tab(self):
        self.directive(
            self.tab_compare,
            "COMPARE — Compare student-entered papers with current Model inputs and outputs. "
            "Non-comparable values are shown as N/A."
        )
        wrapper = tk.Frame(self.tab_compare, bg=DARK["axes"], bd=1, relief="solid")
        wrapper.pack(fill="both", expand=True, padx=10, pady=6)
        tk.Label(
            wrapper, text="Literature evidence versus current Model", bg=DARK["axes"],
            fg=DARK["blue"], font=("Segoe UI", 10, "bold")
        ).pack(anchor="w", padx=10, pady=(10, 2))
        tk.Label(
            wrapper,
            text="This table is based on student-entered paper records. Model RED is thermodynamic screening evidence, not a reported experimental outcome.",
            bg=DARK["axes"], fg=DARK["subtext"], wraplength=1280, justify="left",
            font=("Segoe UI", 8)
        ).pack(anchor="w", padx=10, pady=(0, 8))
        columns = ("Paper", "Paper polymer", "Paper solvents", "Paper ratio", "Model solvents", "Model ratio", "Model RED", "Comparison")
        self.lit_compare_tree = ttk.Treeview(wrapper, columns=columns, show="headings", height=15)
        widths = (180, 130, 175, 120, 175, 110, 190, 365)
        for column, width in zip(columns, widths):
            self.lit_compare_tree.heading(column, text=column)
            self.lit_compare_tree.column(column, width=width, anchor="center")
        self.lit_compare_tree.pack(fill="both", expand=True, padx=10, pady=(0, 8), side="left")
        scrollbar = ttk.Scrollbar(wrapper, orient="vertical", command=self.lit_compare_tree.yview)
        scrollbar.pack(side="right", fill="y", pady=(0, 8))
        self.lit_compare_tree.configure(yscrollcommand=scrollbar.set)
        self.comparison_status = tk.Label(
            self.tab_compare, text="", bg=DARK["figure"], fg=DARK["amber"],
            wraplength=1280, justify="left", font=("Segoe UI", 8)
        )
        self.comparison_status.pack(fill="x", padx=12, pady=(0, 5))
        ttk.Button(
            self.tab_compare, text="Refresh literature comparison", style="Primary.TButton",
            command=self._refresh_comparison
        ).pack(anchor="e", padx=12, pady=(0, 9))


# DECIDE revision: scientific question-and-answer workflow based on student literature records and current Model output.
class PolymerSolvPages(PolymerSolvPages):
    DECIDE_TEXT_FIELDS = (
        ("decide_model_prediction", "1. What does the Model predict for the selected solvent or blend, based on RED and any domain-specific results?", 4),
        ("decide_literature_finding", "2. What did the literature report under its stated polymer grade, composition, concentration, temperature, and processing conditions?", 4),
        ("decide_comparability", "3. Are the literature system and Model system directly comparable? Identify differences in polymer, solvent identity, blend ratio, concentration, temperature, or method.", 4),
        ("decide_agreement_analysis", "4. Where do Model and literature evidence agree, and where do they conflict?", 4),
        ("decide_conflict_mechanisms", "5. If they conflict, what scientific mechanisms or limitations could explain the difference? Consider polymer grade, kinetics, morphology, water content, mixing history, evaporation, or HSP-model limitations.", 4),
        ("decide_recommendation_reason", "6. Which candidate or shortlist do you recommend, and what evidence supports that choice?", 4),
        ("decide_missing_evidence", "7. What evidence is still missing, and what experiment would verify the decision? State a measurable endpoint.", 4),
    )

    def build_decide_tab(self):
        self.directive(
            self.tab_decide,
            "DECIDE — Answer every scientific question using the student-entered literature records and the current Model. "
            "A Model result is screening evidence, not proof of experimental suitability."
        )
        outer = tk.Frame(self.tab_decide, bg=DARK["figure"])
        outer.pack(fill="both", expand=True, padx=10, pady=4)
        canvas = tk.Canvas(outer, bg=DARK["figure"], highlightthickness=0)
        scroll = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        self.decide_form = tk.Frame(canvas, bg=DARK["figure"])
        self.decide_form.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.decide_form, anchor="nw")
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")
        self.decide_form.bind("<Enter>", lambda e: canvas.bind_all("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1*(event.delta/120)), "units")))
        self.decide_form.bind("<Leave>", lambda e: canvas.unbind_all("<MouseWheel>"))

        evidence = tk.Frame(self.decide_form, bg=DARK["axes"], bd=1, relief="solid")
        evidence.pack(fill="x", padx=2, pady=(0, 7))
        tk.Label(evidence, text="Student-entered literature records and current Model", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=(9, 3))
        self.decide_evidence_text = tk.Text(evidence, height=10, bg="#0F172A", fg="white", state="disabled", wrap="word", relief="flat", font=("Segoe UI", 8))
        self.decide_evidence_text.pack(fill="x", padx=12, pady=(0, 9))

        questions = tk.Frame(self.decide_form, bg=DARK["axes"], bd=1, relief="solid")
        questions.pack(fill="both", expand=True, padx=2, pady=(0, 7))
        tk.Label(questions, text="Scientific decision questions — all answers are required", bg=DARK["axes"], fg=DARK["amber"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=(9, 2))
        tk.Label(questions, text="Use the evidence above. Do not treat a low RED value as experimental proof.", bg=DARK["axes"], fg=DARK["subtext"], font=("Segoe UI", 8)).pack(anchor="w", padx=12, pady=(0, 7))
        for key, label, height in self.DECIDE_TEXT_FIELDS:
            tk.Label(questions, text=label, bg=DARK["axes"], fg=DARK["subtext"], wraplength=1260, justify="left", font=("Segoe UI", 8, "bold")).pack(anchor="w", padx=12, pady=(7, 2))
            widget = tk.Text(questions, height=height, bg="#0F172A", fg="white", insertbackground="white", wrap="word", relief="flat", font=("Segoe UI", 8))
            widget.pack(fill="x", padx=12)
            self.student_widgets[key] = widget

        tk.Label(questions, text="4a. Overall relationship between Model and literature evidence (required dropdown)", bg=DARK["axes"], fg=DARK["subtext"], font=("Segoe UI", 8, "bold")).pack(anchor="w", padx=12, pady=(8, 2))
        self.decide_evidence_relationship = ttk.Combobox(questions, state="readonly", values=("", "Agree", "Partly agree", "Conflict", "Not directly comparable"))
        self.decide_evidence_relationship.pack(fill="x", padx=12, pady=(0, 4))
        self.student_widgets["decide_evidence_relationship"] = self.decide_evidence_relationship

        tk.Label(questions, text="Recommendation / shortlist (required dropdown)", bg=DARK["axes"], fg=DARK["subtext"], font=("Segoe UI", 8, "bold")).pack(anchor="w", padx=12, pady=(8, 2))
        self.recommendation_combo = ttk.Combobox(questions, state="readonly")
        self.recommendation_combo.pack(fill="x", padx=12, pady=(0, 4))
        self.recommendation_combo.bind("<<ComboboxSelected>>", self.check_recommendation_gate)
        self.recommendation_status = ttk.Label(questions, text="Select an option to inspect its documented gate.", style="Info.TLabel", wraplength=1240)
        self.recommendation_status.pack(fill="x", padx=12, pady=(0, 10))

    def _decide_evidence_summary(self):
        snapshot = self._model_snapshot()
        lines = [
            "CURRENT MODEL OUTPUT",
            f"Model system: {snapshot['system'][0]} / {snapshot['system'][1]} at {snapshot['ratio_a']:.0f}/{100-snapshot['ratio_a']:.0f} v/v",
            f"Model RED: {snapshot['red']}",
            "Domain results: " + "; ".join(f"{name}: RED {row['red']:.3f} ({row['screening']})" for name, row in snapshot['domains'].items()),
            "Model RED is a thermodynamic screening calculation, not a reported experimental outcome.",
            "",
            "STUDENT-ENTERED LITERATURE RECORDS",
        ]
        if self.no_lit_var.get():
            lines.append("Student selected: no suitable literature data found.")
        elif not self.literature_records:
            lines.append("No student-entered literature records yet.")
        else:
            for i, rec in enumerate(self.literature_records, 1):
                lines.extend((
                    f"Paper {i}: {rec.get('title') or 'Untitled record'}",
                    f"  DOI/URL: {rec.get('doi') or 'Not entered'} | Authors: {rec.get('authors') or 'Not entered'} | Year: {rec.get('year') or 'Not entered'}",
                    f"  Polymer: {rec.get('polymer') or 'Not entered'} | Solvents: {rec.get('solvent_a') or 'Not entered'}" + (f" / {rec.get('solvent_b')}" if rec.get('solvent_b') else ""),
                    f"  Ratio: {rec.get('ratio') or 'Not entered'} | Concentration: {rec.get('concentration') or 'Not entered'} | Temperature: {rec.get('temperature') or 'Not entered'} | Time: {rec.get('time') or 'Not entered'}",
                    f"  Reported observation: {rec.get('observation') or 'Not entered'}",
                    f"  Student notes/caveats: {rec.get('notes') or 'Not entered'}",
                ))
        rows = self._comparison_rows()
        lines.extend(("", "CURRENT COMPARE-PAGE OUTCOME"))
        if rows:
            for row in rows:
                lines.append(f"{row['paper']}: {row['comparison']} | Model RED: {row['red']}")
        else:
            lines.append("No direct literature-to-Model comparison row is available.")
        return "\n".join(lines)

    def refresh_decide_evidence(self):
        if not hasattr(self, "decide_evidence_text"):
            return
        text = self._decide_evidence_summary()
        self.decide_evidence_text.configure(state="normal")
        self.decide_evidence_text.delete("1.0", "end")
        self.decide_evidence_text.insert("1.0", text)
        self.decide_evidence_text.configure(state="disabled")

    def _required_decide_answers_missing(self):
        missing = []
        for key, label, _ in self.DECIDE_TEXT_FIELDS:
            if not self.get_student_field(key):
                missing.append(label.split(". ", 1)[0])
        if not self.get_student_field("decide_evidence_relationship"):
            missing.append("4a. evidence relationship dropdown")
        if not self.recommendation_combo.get().strip():
            missing.append("recommendation / shortlist")
        return missing

    def _step(self, direction):
        pages = (self.tab_define, self.tab_investigate, self.tab_model, self.tab_compare, self.tab_decide, self.tab_reflect)
        current = pages.index(self.notebook.nametowidget(self.notebook.select()))
        next_page = pages[max(0, min(5, current + direction))]
        if current == 4 and next_page == self.tab_reflect:
            missing = self._required_decide_answers_missing()
            if missing:
                messagebox.showwarning("Complete DECIDE", "Answer all DECIDE questions before continuing to REFLECT.\n\nMissing: " + ", ".join(missing))
                return
        self._save_responses()
        self.notebook.select(next_page)
        if next_page == self.tab_compare:
            self._refresh_comparison()
        elif next_page == self.tab_decide:
            self.refresh_decide_evidence()

    def _tab_changed(self, event=None):
        selected = self.notebook.select()
        if selected == str(self.tab_compare):
            self._refresh_comparison()
        elif selected == str(self.tab_decide):
            self.refresh_decide_evidence()

    def _save_responses(self):
        super()._save_responses()
        self.response_state["decide_scientific_answers"] = self._decide_answers()

    def _decide_answers(self):
        answers = {key: self.get_student_field(key) for key, _, _ in self.DECIDE_TEXT_FIELDS}
        answers["decide_evidence_relationship"] = self.get_student_field("decide_evidence_relationship")
        answers["recommendation"] = self.recommendation_combo.get()
        return answers

    def collect_student_fields(self):
        values = super().collect_student_fields()
        values.update(self._decide_answers())
        values["decide_scientific_questions"] = {label: self.get_student_field(key) for key, label, _ in self.DECIDE_TEXT_FIELDS}
        values["decide_scientific_questions"]["4a. Overall relationship between Model and literature evidence"] = self.get_student_field("decide_evidence_relationship")
        return values


# Formal PDF and 300-dpi figure helpers. Imports are local to keep the standard GUI startup lightweight.
def _safe_report_value(value):
    value = str(value).strip() if value is not None else ""
    return value if value else "Not entered"


def _pdf_page_number(canvas, document):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColorRGB(0.25, 0.30, 0.37)
    canvas.drawRightString(A4[0] - 42, 24, f"Page {document.page}")
    canvas.drawString(42, 24, "PolymerSolv Studio — Material Selection Report")
    canvas.restoreState()


def create_material_selection_pdf(filename, report, figures):
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib import colors
    from reportlab.lib.units import cm
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Image, Table, TableStyle
    from xml.sax.saxutils import escape
    from tempfile import TemporaryDirectory
    from pathlib import Path
    import os
    styles = getSampleStyleSheet()
    title = ParagraphStyle("MSRTitle", parent=styles["Title"], textColor=colors.HexColor("#0D4B6B"), fontSize=19, leading=23)
    section = ParagraphStyle("MSRSection", parent=styles["Heading2"], textColor=colors.HexColor("#0D4B6B"), spaceBefore=12, spaceAfter=6)
    body = ParagraphStyle("MSRBody", parent=styles["BodyText"], fontSize=8.5, leading=12, spaceAfter=4)
    small = ParagraphStyle("MSRSmall", parent=body, fontSize=7.5, leading=10)
    story = [Paragraph("PolymerSolv Studio", title), Paragraph("Material Selection Report", title), Spacer(1, 10)]
    metadata = [["Student", escape(_safe_report_value(report["student"]))], ["Date", escape(report["date"])], ["Status", "DRAFT — missing fields are shown as ‘Not entered’."], ["Case", escape(_safe_report_value(report["case"]))]]
    table = Table(metadata, colWidths=(3.0*cm, 13.0*cm))
    table.setStyle(TableStyle([("BACKGROUND", (0,0), (0,-1), colors.HexColor("#E6F3F8")), ("GRID", (0,0), (-1,-1), 0.25, colors.HexColor("#94A3B8")), ("VALIGN", (0,0), (-1,-1), "TOP"), ("FONTNAME", (0,0), (-1,-1), "Helvetica"), ("FONTSIZE", (0,0), (-1,-1), 8), ("LEFTPADDING", (0,0), (-1,-1), 6), ("RIGHTPADDING", (0,0), (-1,-1), 6), ("TOPPADDING", (0,0), (-1,-1), 5), ("BOTTOMPADDING", (0,0), (-1,-1), 5)]))
    story.extend([table, Spacer(1, 11)])
    for heading, content in report["sections"]:
        story.append(Paragraph(escape(heading), section))
        for paragraph in content.split("\n"):
            if paragraph.strip():
                story.append(Paragraph(escape(paragraph), body))
            else:
                story.append(Spacer(1, 4))
    story.append(Paragraph("Figures", section))
    story.append(Paragraph("The following figures use the application's dark Model styling. They are generated from the current Model state and are thermodynamic-screening evidence, not experimental proof of suitability.", body))
    with TemporaryDirectory() as tempdir:
        for index, (caption, figure) in enumerate(figures, 1):
            path = Path(tempdir) / f"figure_{index}.png"
            figure.savefig(path, format="png", dpi=300, facecolor=figure.get_facecolor(), bbox_inches="tight")
            story.append(Spacer(1, 8))
            story.append(Paragraph(f"Figure {index}. {escape(caption)}", small))
            image = Image(str(path), width=16.5*cm, height=9.2*cm, kind="proportional")
            story.append(image)
            if index < len(figures): story.append(PageBreak())
        doc = SimpleDocTemplate(filename, pagesize=A4, leftMargin=1.45*cm, rightMargin=1.45*cm, topMargin=1.45*cm, bottomMargin=1.45*cm, title="PolymerSolv Studio Material Selection Report")
        doc.build(story, onFirstPage=_pdf_page_number, onLaterPages=_pdf_page_number)


# REFLECT revision: automatically generated all-stage report, formal PDF, and publication-resolution figure export.
class PolymerSolvPages(PolymerSolvPages):
    FIGURE_SPECS = (
        ("HSP projections", "HSP projections showing adopted polymer-domain spheres and the current solvent mixture."),
        ("Domain-specific RED screening evidence", "Domain-specific RED values for the candidate solvents and current mixture."),
        ("Binary-mixture sensitivity", "RED sensitivity across the current binary solvent-mixture composition range."),
    )

    def build_reflect_tab(self):
        self.directive(self.tab_reflect, "REFLECT — Review the automatically generated evidence summary, then export a formal draft report or 300 dpi Model figures.")
        top = tk.Frame(self.tab_reflect, bg=DARK["axes"], bd=1, relief="solid")
        top.pack(fill="x", padx=10, pady=(0, 6))
        tk.Label(top, text="Student name", bg=DARK["axes"], fg=DARK["subtext"], font=("Segoe UI", 8, "bold")).grid(row=0, column=0, sticky="w", padx=12, pady=(9, 2))
        name = tk.Entry(top, bg="#0F172A", fg="white", insertbackground="white", relief="flat", font=("Segoe UI", 9))
        name.grid(row=1, column=0, sticky="ew", padx=12, pady=(0, 9))
        self.student_widgets["student_name"] = name
        top.columnconfigure(0, weight=1)
        actions = tk.Frame(top, bg=DARK["axes"]); actions.grid(row=0, column=1, rowspan=2, sticky="e", padx=12)
        ttk.Button(actions, text="Refresh report", style="Primary.TButton", command=self.refresh_reflect_report).pack(side="left", padx=3)
        ttk.Button(actions, text="Download PDF report", style="Primary.TButton", command=self.download_reflect_pdf).pack(side="left", padx=3)
        ttk.Button(actions, text="Export all figures", style="Primary.TButton", command=self.export_all_model_figures).pack(side="left", padx=3)
        report_box = tk.Frame(self.tab_reflect, bg=DARK["axes"], bd=1, relief="solid")
        report_box.pack(fill="both", expand=True, padx=10, pady=(0, 6))
        tk.Label(report_box, text="Automatically generated Material Selection Report preview", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=(9, 3))
        tk.Label(report_box, text="Draft report: ‘Not entered’ identifies missing student information. Literature records are student-entered and unverified by the software.", bg=DARK["axes"], fg=DARK["amber"], font=("Segoe UI", 8), wraplength=1280, justify="left").pack(anchor="w", padx=12, pady=(0, 5))
        report_scroll = ttk.Scrollbar(report_box, orient="vertical")
        self.reflect_report_text = tk.Text(report_box, bg="#0F172A", fg="white", state="disabled", wrap="word", relief="flat", font=("Consolas", 8), yscrollcommand=report_scroll.set)
        report_scroll.config(command=self.reflect_report_text.yview)
        report_scroll.pack(side="right", fill="y", padx=(0, 10), pady=(0, 8))
        self.reflect_report_text.pack(fill="both", expand=True, padx=(12, 0), pady=(0, 8))
        figures = tk.Frame(self.tab_reflect, bg=DARK["axes"], bd=1, relief="solid")
        figures.pack(fill="x", padx=10, pady=(0, 4))
        tk.Label(figures, text="300 dpi publication figure export", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=(8, 3))
        tk.Label(figures, text="Each export keeps the dark Model styling. Choose PNG, JPG/JPEG, TIFF, PDF, or SVG in the save dialog.", bg=DARK["axes"], fg=DARK["subtext"], font=("Segoe UI", 8)).pack(anchor="w", padx=12, pady=(0, 7))
        buttons = tk.Frame(figures, bg=DARK["axes"]); buttons.pack(fill="x", padx=12, pady=(0, 9))
        for index, (name, _) in enumerate(self.FIGURE_SPECS):
            ttk.Button(buttons, text=f"Export {name}", style="Primary.TButton", command=lambda i=index: self.export_one_model_figure(i)).pack(side="left", padx=(0, 6))

    def _report_data(self):
        snap = self._model_snapshot()
        remaining = 100 - snap["ratio_a"]
        problem = self.get_student_field("problem_definition")
        selected_polymer = self.polymer_var.get()
        if selected_polymer == "Other / Custom": selected_polymer = self.custom_var.get().strip() or "Other / Custom (not specified)"
        lit = self.literature_records
        literature_lines = ["Evidence status: Student-entered literature records are unverified by the software."]
        if self.no_lit_var.get(): literature_lines.append("Student selected: no suitable literature data found.")
        elif not lit: literature_lines.append("No literature records entered; search status is not established.")
        else:
            for i, rec in enumerate(lit, 1):
                literature_lines.append(f"Paper {i}: {_safe_report_value(rec.get('title'))}; {_safe_report_value(rec.get('authors'))} ({_safe_report_value(rec.get('year'))}); DOI/URL: {_safe_report_value(rec.get('doi'))}.")
                literature_lines.append(f"  Polymer: {_safe_report_value(rec.get('polymer'))}; system: {_safe_report_value(rec.get('solvent_a'))}" + (f" / {_safe_report_value(rec.get('solvent_b'))}" if rec.get('solvent_b') else "") + f"; ratio: {_safe_report_value(rec.get('ratio'))}; concentration: {_safe_report_value(rec.get('concentration'))}; temperature: {_safe_report_value(rec.get('temperature'))}; time: {_safe_report_value(rec.get('time'))}.")
                literature_lines.append(f"  Reported observation: {_safe_report_value(rec.get('observation'))}; caveats: {_safe_report_value(rec.get('notes'))}.")
        domains = "\n".join(f"{name}: Ra {row['ra']:.3f}; R0 {row['R0']:.3f}; RED {row['red']:.3f}; screening: {row['screening']}." for name, row in snap['domains'].items())
        comparisons = self._comparison_rows()
        compare_lines = ["Model RED is a thermodynamic-screening indicator and is not an experimental observation."]
        if comparisons:
            compare_lines += [f"{row['paper']}: Model system {row['model_system']} at {row['model_ratio']}; model RED {row['red']}; comparison: {row['comparison']}." for row in comparisons]
        else: compare_lines.append("No direct literature-to-Model comparison row is available.")
        decision = []
        for key, label, _ in self.DECIDE_TEXT_FIELDS:
            decision.append(f"{label}\n{_safe_report_value(self.get_student_field(key))}")
        decision.append("Overall relationship between Model and literature evidence\n" + _safe_report_value(self.get_student_field("decide_evidence_relationship")))
        decision.append("Selected recommendation / shortlist\n" + _safe_report_value(self.recommendation_combo.get()))
        gate = self.recommendation_status.cget("text") if hasattr(self, "recommendation_status") else "Not entered"
        decision.append("Documented gate message\n" + _safe_report_value(gate))
        sections = [
            ("1. Define", f"Engineering problem: {_safe_report_value(problem)}\nDecision criteria / constraints: {_safe_report_value(self.get_student_field('criteria_definition'))}\nSelected polymer: {_safe_report_value(selected_polymer)}\nPolarity / class: {_safe_report_value(self.polarity_var.get())} / {_safe_report_value(self.class_var.get())}\nCatalog suggestion: {_safe_report_value(self.suggest_var.get())}. Catalog information is not evidence of dissolution."),
            ("2. Investigate", "\n".join(literature_lines)),
            ("3. Model", f"Teaching model case: {snap['model_polymer']}\nCurrent system: {snap['system'][0]} / {snap['system'][1]} at {snap['ratio_a']:.0f}/{remaining:.0f} v/v.\nMixed HSP: δD {snap['hsp'][0]:.3f}, δP {snap['hsp'][1]:.3f}, δH {snap['hsp'][2]:.3f} MPa^0.5.\n{domains}"),
            ("4. Compare", "\n".join(compare_lines)),
            ("5. Decide", "\n\n".join(decision)),
            ("6. Reflect", "The final report is a draft evidence record. It does not certify material suitability, verify literature sources, or replace experimental validation."),
        ]
        return {"student": self.get_student_field("student_name"), "date": date.today().isoformat(), "case": self.case_data["name"], "sections": sections}

    def refresh_reflect_report(self):
        if not hasattr(self, "reflect_report_text"): return
        report = self._report_data()
        preview = ["POLYMERSOLV STUDIO — MATERIAL SELECTION REPORT", f"Student: {_safe_report_value(report['student'])}", f"Date: {report['date']}", "Status: DRAFT — missing fields are shown as ‘Not entered’.", ""]
        for heading, content in report["sections"]: preview += [heading.upper(), content, ""]
        self.reflect_report_text.configure(state="normal")
        self.reflect_report_text.delete("1.0", "end")
        self.reflect_report_text.insert("1.0", "\n".join(preview))
        self.reflect_report_text.configure(state="disabled")

    def _build_model_figures(self):
        snap = self._model_snapshot(); a, b = snap['system']; ratio = snap['ratio_a']; label = f"{ratio:.0f}/{100-ratio:.0f} {a}/{b}"
        return [
            (self.FIGURE_SPECS[0][1], plot_hsp_projections(self.case_data, snap['hsp'], label)),
            (self.FIGURE_SPECS[1][1], plot_red_evidence(self.case_data, snap['hsp'], label)),
            (self.FIGURE_SPECS[2][1], plot_mixture_sensitivity(self.case_data, a, b, ratio)),
        ]

    def _image_filetypes(self):
        return [("PNG image (300 dpi)", "*.png"), ("JPEG image (300 dpi)", "*.jpg"), ("TIFF image (300 dpi)", "*.tiff"), ("PDF vector graphic", "*.pdf"), ("SVG vector graphic", "*.svg")]

    def _save_figure(self, figure, filename):
        suffix = Path(filename).suffix.lower()
        fmt = {".jpg":"jpeg", ".jpeg":"jpeg", ".tif":"tiff", ".tiff":"tiff", ".png":"png", ".pdf":"pdf", ".svg":"svg"}.get(suffix)
        if not fmt: raise ValueError("Choose PNG, JPG/JPEG, TIFF, PDF, or SVG.")
        figure.savefig(filename, format=fmt, dpi=300, facecolor=figure.get_facecolor(), bbox_inches="tight")

    def export_one_model_figure(self, index):
        from tkinter import filedialog
        name, _ = self.FIGURE_SPECS[index]
        extension = ".png"
        filename = filedialog.asksaveasfilename(title=f"Export {name}", initialfile=name.lower().replace(" ", "_") + extension, defaultextension=extension, filetypes=self._image_filetypes())
        if not filename: return
        figure = self._build_model_figures()[index][1]
        try:
            self._save_figure(figure, filename)
        except Exception as exc:
            messagebox.showerror("Figure export failed", str(exc))
        else:
            messagebox.showinfo("Figure exported", f"Saved 300 dpi figure to:\n{filename}")
        finally:
            plt.close(figure)

    def export_all_model_figures(self):
        from tkinter import filedialog
        from pathlib import Path
        folder = filedialog.askdirectory(title="Choose folder for all 300 dpi Model figures")
        if not folder: return
        choice = filedialog.asksaveasfilename(title="Choose figure format (the filename is used only to choose format)", initialfile="figures.png", defaultextension=".png", filetypes=self._image_filetypes())
        if not choice: return
        suffix = Path(choice).suffix.lower()
        if suffix not in (".png", ".jpg", ".jpeg", ".tif", ".tiff", ".pdf", ".svg"):
            messagebox.showerror("Figure export failed", "Choose a supported figure format."); return
        figures = self._build_model_figures(); paths = []
        try:
            for i, (caption, figure) in enumerate(figures, 1):
                filename = str(Path(folder) / f"PolymerSolv_Model_Figure_{i}{suffix}")
                self._save_figure(figure, filename); paths.append(filename)
        except Exception as exc:
            messagebox.showerror("Figure export failed", str(exc))
        else:
            messagebox.showinfo("Figures exported", "Saved three 300 dpi figures:\n" + "\n".join(paths))
        finally:
            for _, figure in figures: plt.close(figure)

    def download_reflect_pdf(self):
        from tkinter import filedialog
        filename = filedialog.asksaveasfilename(title="Save Material Selection Report", initialfile="PolymerSolv_Material_Selection_Report.pdf", defaultextension=".pdf", filetypes=[("PDF document", "*.pdf")])
        if not filename: return
        figures = self._build_model_figures()
        try:
            create_material_selection_pdf(filename, self._report_data(), figures)
        except Exception as exc:
            messagebox.showerror("PDF export failed", f"Could not create report:\n{exc}")
        else:
            messagebox.showinfo("PDF report saved", f"Draft Material Selection Report saved to:\n{filename}")
        finally:
            for _, figure in figures: plt.close(figure)

    def _tab_changed(self, event=None):
        super()._tab_changed(event)
        if self.notebook.select() == str(self.tab_reflect): self.refresh_reflect_report()

    def collect_student_fields(self):
        values = super().collect_student_fields()
        values["reflect_report"] = self._report_data()
        return values


# DEFINE revision: polarity-controlled class and same-polarity single/binary solvent suggestions.
class PolymerSolvPages(PolymerSolvPages):
    def build_define_tab(self):
        self.directive(self.tab_define, "DEFINE — State your own problem. Suggested solvents and solvent systems are polarity-based catalog matches, not literature evidence or proof of dissolution.")
        container = tk.Frame(self.tab_define, bg=DARK["figure"])
        container.pack(fill="both", expand=True, padx=10, pady=5)
        left = tk.Frame(container, bg=DARK["axes"], bd=1, relief="solid")
        left.pack(side="left", fill="both", expand=True, padx=(0, 5))
        right = tk.Frame(container, bg=DARK["axes"], bd=1, relief="solid")
        right.pack(side="right", fill="both", expand=True, padx=(5, 0))
        tk.Label(left, text="Teaching case for the existing Model", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=8)
        self.case_var = tk.StringVar(value="nafion")
        for cid, case in CASES_DATABASE.items():
            tk.Radiobutton(left, text=case["name"], variable=self.case_var, value=cid, command=self.change_case, bg=DARK["axes"], fg="white", selectcolor="#0F172A", activebackground=DARK["axes"], activeforeground="white").pack(anchor="w", padx=10)
        self.define_context = tk.Text(left, bg="#0F172A", fg="white", wrap="word", relief="flat", font=("Segoe UI", 9))
        self.define_context.pack(fill="both", expand=True, padx=12, pady=10)
        tk.Label(right, text="Your problem definition", bg=DARK["axes"], fg=DARK["blue"], font=("Segoe UI", 10, "bold")).pack(anchor="w", padx=12, pady=10)
        self.student_widgets["problem_definition"] = tk.Text(right, height=7, bg="#0F172A", fg="white", insertbackground="white", wrap="word", relief="flat", font=("Segoe UI", 10))
        self.student_widgets["problem_definition"].pack(fill="x", padx=12, pady=(0, 12))
        tk.Label(right, text="Polymer type (USER-EDITABLE catalog above)", bg=DARK["axes"], fg=DARK["subtext"]).pack(anchor="w", padx=12)
        self.polymer_var = tk.StringVar(value="Nafion")
        polymer_box = ttk.Combobox(right, textvariable=self.polymer_var, values=POLYMER_CHOICES, state="readonly")
        polymer_box.pack(fill="x", padx=12, pady=5); polymer_box.bind("<<ComboboxSelected>>", self._polymer_changed)
        self.custom_var = tk.StringVar(); self.custom_line = tk.Frame(right, bg=DARK["axes"])
        tk.Label(self.custom_line, text="Custom polymer name", bg=DARK["axes"], fg="white").pack(anchor="w")
        tk.Entry(self.custom_line, textvariable=self.custom_var, bg="#0F172A", fg="white", insertbackground="white").pack(fill="x")
        tk.Label(right, text="Solvent polarity", bg=DARK["axes"], fg=DARK["subtext"]).pack(anchor="w", padx=12)
        self.polarity_var = tk.StringVar(value="Any / Not specified")
        self.polarity_combo = ttk.Combobox(right, textvariable=self.polarity_var, values=("Polar", "Nonpolar", "Any / Not specified"), state="readonly")
        self.polarity_combo.pack(fill="x", padx=12, pady=5); self.polarity_combo.bind("<<ComboboxSelected>>", self._polarity_changed)
        tk.Label(right, text="Solvent class", bg=DARK["axes"], fg=DARK["subtext"]).pack(anchor="w", padx=12)
        self.class_var = tk.StringVar(value="Any / Not specified")
        self.class_combo = ttk.Combobox(right, textvariable=self.class_var, values=("Protic", "Aprotic", "Any / Not specified"), state="readonly")
        self.class_combo.pack(fill="x", padx=12, pady=5); self.class_combo.bind("<<ComboboxSelected>>", self._update_suggestions)
        tk.Label(right, text="Suggested single solvents and same-polarity binary solvent systems (catalog matches; NOT validated dissolving solvents)", bg=DARK["axes"], fg=DARK["subtext"], wraplength=500, justify="left").pack(anchor="w", padx=12, pady=(4,0))
        self.suggest_var = tk.StringVar()
        self.suggest_combo = ttk.Combobox(right, textvariable=self.suggest_var, state="readonly")
        self.suggest_combo.pack(fill="x", padx=12, pady=6)
        self.define_note = tk.Label(right, text="", bg=DARK["axes"], fg=DARK["amber"], wraplength=500, justify="left")
        self.define_note.pack(fill="x", padx=12, pady=6)
        self._polarity_changed()

    def _polarity_changed(self, event=None):
        if self.polarity_var.get() == "Nonpolar":
            self.class_var.set("Not applicable for nonpolar")
            self.class_combo.configure(values=("Not applicable for nonpolar",), state="disabled")
        else:
            self.class_combo.configure(values=("Protic", "Aprotic", "Any / Not specified"), state="readonly")
            if self.class_var.get() == "Not applicable for nonpolar":
                self.class_var.set("Any / Not specified")
        self._update_suggestions()

    def _polymer_changed(self, event=None):
        if self.polymer_var.get() == "Other / Custom":
            self.custom_line.pack(fill="x", padx=12, pady=6)
        else:
            self.custom_line.pack_forget()
        self._update_suggestions()

    def _update_suggestions(self, event=None):
        polarity = self.polarity_var.get()
        solvent_class = self.class_var.get()
        eligible = []
        for name, (solvent_polarity, kind) in SOLVENT_POLARITY_CATALOG.items():
            if polarity != "Any / Not specified" and solvent_polarity != polarity:
                continue
            if polarity != "Nonpolar" and solvent_class not in ("Any / Not specified", "Not applicable for nonpolar") and kind != solvent_class:
                continue
            eligible.append(name)
        singles = eligible[:]
        systems = [f"{eligible[i]} / {eligible[j]}" for i in range(len(eligible)) for j in range(i + 1, len(eligible))]
        values = singles + systems
        self.suggest_combo["values"] = values
        if self.suggest_var.get() not in values:
            self.suggest_var.set(values[0] if values else "")
        message = f"{len(singles)} single-solvent and {len(systems)} same-polarity binary-system catalog match(es)."
        if polarity == "Nonpolar":
            message += " Solvent class is not applicable to nonpolar selections; no protic/aprotic label is shown."
        message += f" The Model still uses {self.case_data['polymer']['name']}; changing polymer selection does not make RED a prediction for another polymer."
        self.define_note.config(text=message)


if __name__ == "__main__":
    main()
