"""Core data structures and ML pipeline for OligoTox."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Dict


class BackboneClass(str, Enum):
    PS = "PS"
    PO = "PO"
    PMO = "PMO"
    LNA = "LNA"
    PNA = "PNA"
    F2_ANA = "F2_ANA"
    MORPHOLINO = "MORPHOLINO"


class SugarModification(str, Enum):
    MOE = "MOE"
    OMe = "OMe"
    LNA_SUGAR = "LNA_SUGAR"
    NONE = "NONE"


class ToxicityEndpoint(str, Enum):
    HEPATOTOXICITY = "HEPATOTOXICITY"
    NEPHROTOXICITY = "NEPHROTOXICITY"
    IMMUNOTOXICITY = "IMMUNOTOXICITY"
    COMPLEMENT = "COMPLEMENT"
    COAGULOPATHY = "COAGULOPATHY"


def compute_gc_content(sequence: str) -> float:
    """Compute GC content of a DNA/RNA sequence."""
    if not sequence:
        raise ValueError("Sequence must not be empty.")
    seq = sequence.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


def compute_cpg_ratio(sequence: str) -> float:
    """Compute CpG dinucleotide ratio."""
    seq = sequence.upper()
    denom = max(len(seq) - 1, 1)
    return seq.count("CG") / denom


@dataclass
class Oligonucleotide:
    oligo_id: str
    sequence: str
    backbone: BackboneClass
    modifications: List[str] = field(default_factory=list)
    gc_content: float = field(init=False)

    def __post_init__(self) -> None:
        self.gc_content = compute_gc_content(self.sequence)


@dataclass
class ToxicityRecord:
    oligo_id: str
    endpoint: ToxicityEndpoint
    value: float  # [0, 1]
    cell_system: str
    confidence: float  # [0, 1]

    def __post_init__(self) -> None:
        if not (0.0 <= self.value <= 1.0):
            raise ValueError(f"value must be in [0,1], got {self.value}")
        if not (0.0 <= self.confidence <= 1.0):
            raise ValueError(f"confidence must be in [0,1], got {self.confidence}")


@dataclass
class FeatureVector:
    oligo_id: str
    features: Dict[str, float]


def extract_features(oligo: Oligonucleotide) -> FeatureVector:
    """Extract a fixed feature vector from an Oligonucleotide."""
    feats: Dict[str, float] = {
        "gc_content": oligo.gc_content,
        "cpg_ratio": compute_cpg_ratio(oligo.sequence),
        "sequence_length": float(len(oligo.sequence)),
        "is_ps": 1.0 if oligo.backbone == BackboneClass.PS else 0.0,
        "is_lna": 1.0 if oligo.backbone == BackboneClass.LNA else 0.0,
        "n_modifications": float(len(oligo.modifications)),
    }
    return FeatureVector(oligo_id=oligo.oligo_id, features=feats)


class OligotoxPipeline:
    """Lightweight ML pipeline for oligonucleotide toxicity prediction."""

    def __init__(self, model_type: str = "xgb") -> None:
        self.model_type = model_type
        self._model = None
        self._feature_order: List[str] = []

    def feature_matrix(self, feature_vectors: List[FeatureVector]) -> List[List[float]]:
        """Convert FeatureVectors to a 2-D list (rows = samples)."""
        if not self._feature_order and feature_vectors:
            self._feature_order = sorted(feature_vectors[0].features.keys())
        return [
            [fv.features[k] for k in self._feature_order]
            for fv in feature_vectors
        ]

    def fit(self, X: List[FeatureVector], y: List[float]) -> None:
        """Fit a LogisticRegression model (sklearn imported lazily)."""
        from sklearn.linear_model import LogisticRegression  # noqa: PLC0415

        X_mat = self.feature_matrix(X)
        y_bin = [int(v > 0.5) for v in y]
        self._model = LogisticRegression(max_iter=1000)
        self._model.fit(X_mat, y_bin)

    def predict_proba(self, X: List[FeatureVector]) -> List[float]:
        """Return predicted probabilities for the positive class."""
        if self._model is None:
            raise RuntimeError("Pipeline has not been fitted yet.")
        X_mat = self.feature_matrix(X)
        probs = self._model.predict_proba(X_mat)
        return [float(p[1]) for p in probs]
