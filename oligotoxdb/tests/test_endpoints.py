"""Tests for the endpoint registry."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from oligotoxdb.endpoints import (
    ENDPOINTS, ENDPOINT_BY_NAME, ENDPOINTS_BY_SYSTEM,
    ALL_ENDPOINT_NAMES, PRIMARY_ENDPOINT_NAMES, ASSAY_SYSTEMS,
)


def test_endpoint_count():
    assert len(ENDPOINTS) == 47


def test_all_endpoints_have_required_fields():
    for ep in ENDPOINTS:
        assert ep.name, f"Empty name in {ep}"
        assert ep.display, f"Empty display in {ep}"
        assert ep.assay_system in ASSAY_SYSTEMS, f"Invalid system: {ep.assay_system}"
        assert ep.unit, f"Empty unit in {ep}"
        assert ep.mechanism, f"Empty mechanism in {ep}"


def test_endpoint_names_unique():
    names = [ep.name for ep in ENDPOINTS]
    assert len(names) == len(set(names)), "Duplicate endpoint names found"


def test_lookup_by_name():
    assert "Cell_viability_ATPLite" in ENDPOINT_BY_NAME
    assert "aPTT" in ENDPOINT_BY_NAME
    assert "IFNa" in ENDPOINT_BY_NAME
    assert ENDPOINT_BY_NAME["aPTT"].assay_system == "Platelet"


def test_endpoints_by_system_coverage():
    all_in_systems = []
    for system_eps in ENDPOINTS_BY_SYSTEM.values():
        all_in_systems.extend(ep.name for ep in system_eps)
    # every endpoint should appear in exactly one system
    assert set(all_in_systems) == set(ALL_ENDPOINT_NAMES)


def test_primary_endpoints_subset():
    primary_set = set(PRIMARY_ENDPOINT_NAMES)
    all_set = set(ALL_ENDPOINT_NAMES)
    assert primary_set.issubset(all_set), "Primary endpoints must be a subset of all endpoints"
    assert len(PRIMARY_ENDPOINT_NAMES) > 0


def test_mechanisms_covered():
    mechanisms = {ep.mechanism for ep in ENDPOINTS}
    required = {"hepatotoxicity", "nephrotoxicity", "immunotoxicity", "complement", "coagulopathy"}
    assert required.issubset(mechanisms)


def test_higher_is_toxic_consistent():
    # Viability endpoints should have higher_is_toxic = False
    viability = [ep for ep in ENDPOINTS if "viability" in ep.name.lower()]
    for ep in viability:
        assert not ep.higher_is_toxic, f"{ep.name} should have higher_is_toxic=False"

    # LDH should be toxic when high
    ldh = ENDPOINT_BY_NAME["LDH_release"]
    assert ldh.higher_is_toxic
