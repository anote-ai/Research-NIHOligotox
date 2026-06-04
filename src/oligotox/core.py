"""Core data structures and ML pipeline for OligoTox."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional


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


# Motif weights for known toxic sequence patterns.
# Keys are uppercase motif strings; values are risk scores in [0, 1].
_MOTIF_WEIGHTS: Dict[str, float] = {
    "CCGG": 0.30,   # CpG-dense, complement activation risk
    "GGGG": 0.25,   # G-quadruplex forming, coagulopathy
    "TTTT": 0.10,   # poly-T, protein-binding artefact
    "GCGC": 0.20,   # high-GC palindrome, immune stimulation
    "CGCG": 0.20,   # CpG-rich palindrome
    "AAAA": 0.05,   # mild non-specific binding
}


def compute_motif_score(sequence: str) -> float:
    """Compute a motif-based toxicity risk score for a sequence.

    Scans for known problematic k-mer patterns and sums their weighted
    contributions, capped at 1.0.

    Args:
        sequence: DNA/RNA sequence (case-insensitive).

    Returns:
        A float in [0, 1] representing cumulative motif risk.
    """
    seq = sequence.upper()
    total = 0.0
    for motif, weight in _MOTIF_WEIGHTS.items():
        count = sum(1 for i in range(len(seq) - len(motif) + 1) if seq[i:i + len(motif)] == motif)
        total += count * weight
    return min(total, 1.0)


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
        "motif_score": compute_motif_score(oligo.sequence),
    }
    return FeatureVector(oligo_id=oligo.oligo_id, features=feats)


class ToxicityPredictor:
    """Rule-based heuristic toxicity predictor.

    Applies domain-knowledge rules derived from oligonucleotide pharmacology:

    - PS backbone increases hepatotoxicity risk via protein binding.
    - High GC content (>0.60) increases complement activation risk.
    - CpG-rich sequences trigger immunotoxicity (TLR9 pathway).
    - G-quadruplex-prone sequences (GGGG motif) raise coagulopathy risk.
    - Motif score contributes across all endpoints.
    """

    # Hepatotoxicity thresholds
    _HEPATOTOX_PS_BOOST: float = 0.25
    _HEPATOTOX_BASE: float = 0.10

    # Complement activation threshold for GC content
    _COMPLEMENT_GC_THRESHOLD: float = 0.60
    _COMPLEMENT_GC_BOOST: float = 0.30

    # CpG immunotoxicity threshold
    _IMMUNO_CPG_THRESHOLD: float = 0.05
    _IMMUNO_CPG_BOOST: float = 0.30

    def predict(
        self,
        oligo: Oligonucleotide,
        endpoint: ToxicityEndpoint,
    ) -> float:
        """Return a predicted toxicity score in [0, 1] for the given endpoint.

        Args:
            oligo: The oligonucleotide to assess.
            endpoint: The toxicity endpoint of interest.

        Returns:
            Predicted toxicity probability in [0, 1].
        """
        motif = compute_motif_score(oligo.sequence)
        gc = oligo.gc_content
        cpg = compute_cpg_ratio(oligo.sequence)

        if endpoint == ToxicityEndpoint.HEPATOTOXICITY:
            score = self._HEPATOTOX_BASE + motif * 0.20
            if oligo.backbone == BackboneClass.PS:
                score += self._HEPATOTOX_PS_BOOST
            if gc > 0.55:
                score += 0.10

        elif endpoint == ToxicityEndpoint.COMPLEMENT:
            score = 0.05 + motif * 0.15
            if gc >= self._COMPLEMENT_GC_THRESHOLD:
                score += self._COMPLEMENT_GC_BOOST
            if oligo.backbone in (BackboneClass.PS, BackboneClass.LNA):
                score += 0.10

        elif endpoint == ToxicityEndpoint.IMMUNOTOXICITY:
            score = 0.05 + motif * 0.25
            if cpg >= self._IMMUNO_CPG_THRESHOLD:
                score += self._IMMUNO_CPG_BOOST
            if oligo.backbone == BackboneClass.PS:
                score += 0.10

        elif endpoint == ToxicityEndpoint.NEPHROTOXICITY:
            score = 0.05 + motif * 0.15
            if len(oligo.sequence) > 25:
                score += 0.10
            if oligo.backbone == BackboneClass.PS:
                score += 0.10

        elif endpoint == ToxicityEndpoint.COAGULOPATHY:
            score = 0.05 + motif * 0.30
            seq = oligo.sequence.upper()
            g4_count = sum(1 for i in range(len(seq) - 3) if seq[i:i + 4] == "GGGG")
            score += min(g4_count * 0.15, 0.45)

        else:
            score = 0.05 + motif * 0.10

        return float(min(max(score, 0.0), 1.0))

    def predict_all_endpoints(
        self,
        oligo: Oligonucleotide,
    ) -> Dict[str, float]:
        """Return predicted toxicity for all ToxicityEndpoints."""
        return {
            ep.value: self.predict(oligo, ep)
            for ep in ToxicityEndpoint
        }

    def risk_category(self, score: float) -> str:
        """Map a toxicity score to a qualitative risk band."""
        if score < 0.20:
            return "low"
        if score < 0.50:
            return "moderate"
        if score < 0.75:
            return "high"
        return "very_high"


class OligotoxPipeline:
    """Lightweight ML pipeline for oligonucleotide toxicity prediction."""

    def __init__(self, model_type: str = "xgb") -> None:
        self.model_type = model_type
        self._model: Optional[object] = None
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
        self._model.fit(X_mat, y_bin)  # type: ignore[union-attr]

    def predict_proba(self, X: List[FeatureVector]) -> List[float]:
        """Return predicted probabilities for the positive class."""
        if self._model is None:
            raise RuntimeError("Pipeline has not been fitted yet.")
        X_mat = self.feature_matrix(X)
        probs = self._model.predict_proba(X_mat)  # type: ignore[union-attr]
        return [float(p[1]) for p in probs]
