"""
All 47 OligoToxDB toxicity endpoints, organized by assay system and mechanism.

This is the single source of truth for endpoint naming across the entire codebase.
Every endpoint has: name, assay system, unit, mechanism, and IC50 thresholds.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Endpoint:
    name: str                  # canonical key used in all tables
    display: str               # human-readable label for UI
    assay_system: str          # PHH | KidneyOrganoid | PBMC | Platelet | MPS
    unit: str                  # measurement unit
    mechanism: str             # hepatotoxicity | nephrotoxicity | immunotoxicity | complement | coagulopathy | thrombocytopenia
    higher_is_toxic: bool      # True if higher value = more toxic (e.g., LDH); False for viability
    ic50_inactive_um: float = 50.0   # IC50 threshold for "inactive" classification
    ic50_high_um: float = 1.0        # IC50 threshold for "high" classification


# ─── All 47 endpoints ────────────────────────────────────────────────────────

ENDPOINTS: list[Endpoint] = [

    # ── Hepatotoxicity — Primary Human Hepatocytes ────────────────────────
    Endpoint("Cell_viability_ATPLite",   "Cell Viability (ATP-lite)",       "PHH",          "% vehicle", "hepatotoxicity",   False),
    Endpoint("LDH_release",              "LDH Release",                     "PHH",          "% vehicle", "hepatotoxicity",   True),
    Endpoint("ALT_secretion",            "ALT Secretion",                   "PHH",          "U/L",       "hepatotoxicity",   True),
    Endpoint("AST_secretion",            "AST Secretion",                   "PHH",          "U/L",       "hepatotoxicity",   True),
    Endpoint("Caspase_3_7",              "Caspase 3/7 Activity",            "PHH",          "RLU",       "hepatotoxicity",   True),
    Endpoint("HMGB1_release",            "HMGB1 Release (Necrosis)",        "PHH",          "ng/mL",     "hepatotoxicity",   True),
    Endpoint("ROS_CellROX",              "Reactive Oxygen Species (ROS)",   "PHH",          "AU",        "hepatotoxicity",   True),
    Endpoint("Lipid_OilRedO",            "Lipid Accumulation (Oil Red O)",  "PHH",          "AU",        "hepatotoxicity",   True),
    Endpoint("Mitochondrial_JC1",        "Mitochondrial Potential (JC-1)",  "PHH",          "ratio",     "hepatotoxicity",   False),
    Endpoint("Bile_transport_CDFDA",     "Bile Transport (CDFDA)",          "PHH",          "% control", "hepatotoxicity",   False),
    Endpoint("Gene_module_hepatotox",    "Hepatotoxicity Gene Module",      "PHH",          "score",     "hepatotoxicity",   True),
    Endpoint("RNAseq_hepatotoxicity",    "RNA-seq (Hepatotoxicity PC)",     "PHH",          "log2FC",    "hepatotoxicity",   True),

    # ── Nephrotoxicity — Kidney Proximal Tubule Organoids ─────────────────
    Endpoint("Organoid_viability_3D",    "Organoid Viability (CellTiter-Glo 3D)", "KidneyOrganoid", "% vehicle", "nephrotoxicity", False),
    Endpoint("KIM1_secretion",           "KIM-1 Secretion",                 "KidneyOrganoid", "ng/mL",   "nephrotoxicity",   True),
    Endpoint("NGAL_secretion",           "NGAL Secretion",                  "KidneyOrganoid", "ng/mL",   "nephrotoxicity",   True),
    Endpoint("Organoid_caspase_3_7",     "Caspase 3/7 (Organoid)",         "KidneyOrganoid", "RLU",     "nephrotoxicity",   True),
    Endpoint("Lysosomal_LysoTracker",    "Lysosomal Accumulation",          "KidneyOrganoid", "AU",      "nephrotoxicity",   True),
    Endpoint("TightJunction_ZO1",        "Tight Junction Integrity (ZO-1)", "KidneyOrganoid", "AU",      "nephrotoxicity",   False),
    Endpoint("Mitochondrial_Mito",       "Mitochondrial Stress (Mitotracker)", "KidneyOrganoid", "AU",  "nephrotoxicity",   True),
    Endpoint("Gene_module_nephrotox",    "Nephrotoxicity Gene Module",      "KidneyOrganoid", "score",   "nephrotoxicity",   True),

    # ── Immunotoxicity — PBMC Co-culture ─────────────────────────────────
    Endpoint("IFNa",                     "IFN-α",                           "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("IFNg",                     "IFN-γ",                           "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("IL1b",                     "IL-1β",                           "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("IL6",                      "IL-6",                            "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("IL8",                      "IL-8 (CXCL8)",                    "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("IL10",                     "IL-10",                           "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("IL12p70",                  "IL-12p70",                        "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("TNFa",                     "TNF-α",                           "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("MCP1",                     "MCP-1 (CCL2)",                    "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("IP10",                     "IP-10 (CXCL10)",                  "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("GM_CSF",                   "GM-CSF",                          "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("IL4",                      "IL-4",                            "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("IL17A",                    "IL-17A",                          "PBMC",         "pg/mL",     "immunotoxicity",   True),
    Endpoint("TLR9_activation",          "TLR9 Activation (NF-κB)",        "PBMC",         "RLU",       "immunotoxicity",   True),
    Endpoint("NK_CD69",                  "NK Cell Activation (CD69+)",      "PBMC",         "% cells",   "immunotoxicity",   True),

    # ── Complement — PBMC ────────────────────────────────────────────────
    Endpoint("C3a",                      "Complement C3a",                  "PBMC",         "ng/mL",     "complement",       True),
    Endpoint("C5a",                      "Complement C5a",                  "PBMC",         "ng/mL",     "complement",       True),
    Endpoint("C1q_binding",              "C1q Binding",                     "PBMC",         "AU",        "complement",       True),

    # ── Coagulopathy — Platelet-Rich Plasma ──────────────────────────────
    Endpoint("Platelet_aggregation",     "Platelet Aggregation (LTA)",      "Platelet",     "% max",     "thrombocytopenia", True),
    Endpoint("P_selectin_CD62P",         "P-Selectin (CD62P+)",             "Platelet",     "% cells",   "thrombocytopenia", True),
    Endpoint("CD63_degran",              "CD63 (Platelet Degranulation)",   "Platelet",     "% cells",   "thrombocytopenia", True),
    Endpoint("aPTT",                     "aPTT (Clotting Time)",            "Platelet",     "seconds",   "coagulopathy",     True),
    Endpoint("PT",                       "PT (Prothrombin Time)",           "Platelet",     "seconds",   "coagulopathy",     True),
    Endpoint("FactorXa_inhib",           "Factor Xa Inhibition",            "Platelet",     "% inhibit", "coagulopathy",     True),
    Endpoint("Thrombin_inhib",           "Thrombin Inhibition",             "Platelet",     "% inhibit", "coagulopathy",     True),

    # ── Multi-organ MPS ───────────────────────────────────────────────────
    Endpoint("MPS_ALT",                  "Chip Effluent ALT (MPS)",         "MPS",          "U/L",       "hepatotoxicity",   True),
    Endpoint("MPS_KIM1",                 "Chip Effluent KIM-1 (MPS)",       "MPS",          "ng/mL",     "nephrotoxicity",   True),
]

# Quick lookup dictionaries
ENDPOINT_BY_NAME: dict[str, Endpoint] = {ep.name: ep for ep in ENDPOINTS}
ENDPOINTS_BY_SYSTEM: dict[str, list[Endpoint]] = {}
for ep in ENDPOINTS:
    ENDPOINTS_BY_SYSTEM.setdefault(ep.assay_system, []).append(ep)
ENDPOINTS_BY_MECHANISM: dict[str, list[Endpoint]] = {}
for ep in ENDPOINTS:
    ENDPOINTS_BY_MECHANISM.setdefault(ep.mechanism, []).append(ep)

# Flat name lists for ML use
ALL_ENDPOINT_NAMES: list[str] = [ep.name for ep in ENDPOINTS]
PRIMARY_ENDPOINT_NAMES: list[str] = [
    ep.name for ep in ENDPOINTS
    if ep.assay_system in ("PHH", "PBMC", "Platelet")
    and ep.mechanism != "gene_module"
    and "RNAseq" not in ep.name
]

ASSAY_SYSTEMS = ["PHH", "KidneyOrganoid", "PBMC", "Platelet", "MPS"]
MECHANISMS = ["hepatotoxicity", "nephrotoxicity", "immunotoxicity", "complement", "coagulopathy", "thrombocytopenia"]
