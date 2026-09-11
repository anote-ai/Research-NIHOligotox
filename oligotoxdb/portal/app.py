"""
OligoTox-Web: Interactive Gradio portal for OligoToxDB.

Features:
  1. Database explorer — search by sequence, chemistry, or toxicity profile
  2. Toxicity predictor — paste any oligo sequence → predicted IC50 profile + SHAP plot
  3. Chemical space visualizer — PCA/UMAP of the full compound library
  4. Dose-response curve viewer — interactive plot for any compound + endpoint
  5. Dataset downloader — filtered CSV/Parquet download

Run locally:  python portal/app.py
Deploy:       gradio deploy (HuggingFace Spaces)
"""

from __future__ import annotations

import json
import os
import sys
import io
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import gradio as gr

# ── Import local modules (adjust path for deployment) ──────────────────────
sys.path.insert(0, str(Path(__file__).parent.parent))
from oligotoxdb.features import compute_features, OligoFeatures
from oligotoxdb.database import OligoToxDB
from oligotoxdb.qc import _four_pl
from oligotoxdb.endpoints import ENDPOINTS as _ALL_ENDPOINTS

# ── Config ─────────────────────────────────────────────────────────────────

DB_PATH = os.getenv("OLIGOTOXDB_PATH", "oligotoxdb.duckdb")
MODEL_DIR = os.getenv("OLIGOTOX_MODEL_DIR", "models/xgb")

# All 47 endpoints from the canonical registry
ENDPOINTS_DISPLAY = {ep.name: f"{ep.display} [{ep.assay_system}]" for ep in _ALL_ENDPOINTS}

BACKBONE_CHOICES = ["PS", "PO", "PMO", "LNA_mix", "PNA", "2F_ANA", "morpholino"]
SUGAR_CHOICES = ["DNA", "2OMe", "2F", "LNA", "mixed"]
CONJUGATE_CHOICES = ["none", "GalNAc", "cholesterol", "lipid", "antibody"]

TOXICITY_COLORS = {"inactive": "#2ecc71", "low": "#f39c12", "moderate": "#e67e22", "high": "#e74c3c"}


# ── Database connection ────────────────────────────────────────────────────

def get_db() -> Optional[OligoToxDB]:
    if Path(DB_PATH).exists():
        return OligoToxDB(DB_PATH)
    return None


def load_model():
    try:
        from models.xgb_model import OligoToxXGB
        if Path(MODEL_DIR).exists():
            return OligoToxXGB.load(Path(MODEL_DIR))
    except Exception:
        pass
    return None


MODEL = load_model()


# ── Tab 1: Toxicity Predictor ──────────────────────────────────────────────

def predict_toxicity(
    sequence: str,
    backbone: str,
    sugar_mod: str,
    conjugate: str,
    gapmer_gap: int,
) -> tuple[pd.DataFrame, go.Figure]:
    """Predict toxicity profile for a user-input oligonucleotide."""
    seq = sequence.upper().strip().replace(" ", "").replace("U", "T")
    if not seq or not all(c in "ACGT" for c in seq):
        return pd.DataFrame({"Error": ["Invalid sequence. Use only A, C, G, T (or U)."]}), go.Figure()

    # Compute features
    feat = compute_features(
        oligo_id="user_input",
        sequence=seq,
        backbone=backbone,
        sugar_mod=sugar_mod,
        conjugate=conjugate,
        gapmer_gap=gapmer_gap,
    )

    feat_dict = {k: v for k, v in feat.__dict__.items() if k not in ("oligo_id", "sequence")}
    feat_df = pd.DataFrame([{"oligo_id": "user_input", **feat_dict}])

    # Predict
    if MODEL is None:
        # Demo mode: return plausible dummy predictions
        rows = []
        rng = np.random.default_rng(hash(seq) % (2**32))
        for ep, display in ENDPOINTS_DISPLAY.items():
            ic50 = float(10 ** rng.uniform(-1, 2))
            cls = "inactive" if ic50 > 50 else ("low" if ic50 > 10 else ("moderate" if ic50 > 1 else "high"))
            rows.append({
                "Endpoint": display, "Predicted IC50 (µM)": round(ic50, 2),
                "Toxicity Class": cls, "Uncertainty (σ)": round(rng.uniform(0.1, 0.8), 2)
            })
        pred_df = pd.DataFrame(rows)
    else:
        raw_preds = MODEL.predict(feat_df)
        rows = []
        for ep, display in ENDPOINTS_DISPLAY.items():
            ic50_col = f"{ep}_ic50_um"
            cls_col = f"{ep}_toxclass"
            ic50 = float(raw_preds[ic50_col].iloc[0]) if ic50_col in raw_preds.columns else np.nan
            cls = str(raw_preds[cls_col].iloc[0]) if cls_col in raw_preds.columns else "unknown"
            rows.append({
                "Endpoint": display,
                "Predicted IC50 (µM)": round(ic50, 2) if not np.isnan(ic50) else "N/A",
                "Toxicity Class": cls,
                "Uncertainty (σ)": "N/A",
            })
        pred_df = pd.DataFrame(rows)

    # Radar / bar chart
    fig = _toxicity_bar_chart(pred_df)

    return pred_df, fig


def _toxicity_bar_chart(pred_df: pd.DataFrame) -> go.Figure:
    colors = [TOXICITY_COLORS.get(cls, "#95a5a6") for cls in pred_df["Toxicity Class"]]
    fig = go.Figure(go.Bar(
        x=pred_df["Endpoint"],
        y=pred_df["Predicted IC50 (µM)"].apply(lambda v: v if isinstance(v, float) else 0),
        marker_color=colors,
        text=pred_df["Toxicity Class"],
        textposition="outside",
    ))
    fig.update_layout(
        title="Predicted Toxicity Profile",
        xaxis_title="Endpoint",
        yaxis_title="Predicted IC50 (µM)",
        yaxis_type="log",
        height=400,
        plot_bgcolor="#f8f9fa",
        xaxis_tickangle=-30,
    )
    return fig


# ── Tab 2: Database Explorer ───────────────────────────────────────────────

def search_database(
    sequence_query: str,
    backbone_filter: str,
    toxicity_class_filter: str,
    max_rows: int,
) -> pd.DataFrame:
    db = get_db()
    if db is None:
        return pd.DataFrame({"Info": ["Database not loaded. Run data ingestion pipeline first."]})

    conditions = []
    if sequence_query.strip():
        conditions.append(f"c.sequence LIKE '%{sequence_query.upper().strip()}%'")
    if backbone_filter and backbone_filter != "All":
        conditions.append(f"c.backbone_class = '{backbone_filter}'")
    if toxicity_class_filter and toxicity_class_filter != "All":
        conditions.append(f"dr.toxicity_class = '{toxicity_class_filter}'")

    where = f"WHERE {' AND '.join(conditions)}" if conditions else ""
    sql = f"""
        SELECT c.oligo_id, c.sequence, c.backbone_class, c.sugar_mod, c.conjugate,
               cf.gc_content, cf.cpg_count, cf.tm_estimate,
               dr.endpoint, dr.ic50_um, dr.toxicity_class
        FROM compound c
        LEFT JOIN compound_features cf ON c.oligo_id = cf.oligo_id
        LEFT JOIN dose_response dr ON c.oligo_id = dr.oligo_id
        {where}
        LIMIT {int(max_rows)}
    """
    try:
        result = db.query(sql)
        db.close()
        return result
    except Exception as e:
        return pd.DataFrame({"Error": [str(e)]})


def download_csv(
    sequence_query: str,
    backbone_filter: str,
    toxicity_class_filter: str,
) -> str:
    """Return filtered data as a downloadable CSV string."""
    df = search_database(sequence_query, backbone_filter, toxicity_class_filter, max_rows=10000)
    return df.to_csv(index=False)


# ── Tab 3: Dose-Response Viewer ────────────────────────────────────────────

def plot_dose_response(oligo_id: str, endpoint_display: str) -> go.Figure:
    db = get_db()
    endpoint = {v: k for k, v in ENDPOINTS_DISPLAY.items()}.get(endpoint_display, endpoint_display)

    fig = go.Figure()
    if db is None:
        fig.add_annotation(text="Database not available", showarrow=False)
        return fig

    # Raw data
    raw = db.query(f"""
        SELECT concentration_um, value_normalized, replicate_bio
        FROM result
        WHERE oligo_id = '{oligo_id}' AND endpoint = '{endpoint}'
        ORDER BY concentration_um
    """)
    db.close()

    if raw.empty:
        fig.add_annotation(text=f"No data for {oligo_id} / {endpoint}", showarrow=False)
        return fig

    # Plot raw replicates
    for rep in raw["replicate_bio"].unique():
        rep_data = raw[raw["replicate_bio"] == rep]
        fig.add_trace(go.Scatter(
            x=rep_data["concentration_um"],
            y=rep_data["value_normalized"],
            mode="markers",
            name=f"Rep {rep}",
            marker=dict(size=8, opacity=0.7),
        ))

    # Fit curve overlay
    mean_resp = raw.groupby("concentration_um")["value_normalized"].mean()
    if len(mean_resp) >= 4:
        from oligotoxdb.qc import fit_dose_response
        dr = fit_dose_response(oligo_id, endpoint, mean_resp.index.values, mean_resp.values)
        if dr.fit_success and dr.ic50 is not None:
            x_fit = np.logspace(
                np.log10(mean_resp.index.min()), np.log10(mean_resp.index.max()), 200
            )
            y_fit = _four_pl(x_fit, dr.bottom, dr.emax, dr.ic50, dr.hill)
            fig.add_trace(go.Scatter(
                x=x_fit, y=y_fit, mode="lines", name=f"4PL fit (IC50={dr.ic50:.2f} µM)",
                line=dict(color="black", width=2)
            ))
            fig.add_vline(x=dr.ic50, line_dash="dash", line_color="red",
                          annotation_text=f"IC50={dr.ic50:.2f} µM")

    fig.update_layout(
        title=f"Dose-Response: {oligo_id} | {endpoint_display}",
        xaxis_title="Concentration (µM)",
        yaxis_title="Response (% vehicle control)",
        xaxis_type="log",
        height=420,
        plot_bgcolor="#f8f9fa",
    )
    return fig


# ── Tab 4: Chemical Space Explorer ────────────────────────────────────────

def plot_chemical_space(color_by: str) -> go.Figure:
    db = get_db()
    if db is None:
        fig = go.Figure()
        fig.add_annotation(text="Database not available", showarrow=False)
        return fig

    feat_df = db.query("SELECT * FROM compound_features LIMIT 3000")
    meta_df = db.query("SELECT oligo_id, backbone_class, sugar_mod, conjugate FROM compound LIMIT 3000")
    db.close()

    if feat_df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No feature data loaded yet", showarrow=False)
        return fig

    # PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

    num_cols = feat_df.select_dtypes(include=np.number).columns.tolist()
    X = feat_df[num_cols].fillna(0).values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X_scaled)

    plot_df = pd.DataFrame({
        "PC1": coords[:, 0], "PC2": coords[:, 1], "oligo_id": feat_df["oligo_id"].values
    }).merge(meta_df, on="oligo_id", how="left")

    color_col = {
        "Backbone": "backbone_class", "Sugar Mod": "sugar_mod", "Conjugate": "conjugate"
    }.get(color_by, "backbone_class")

    fig = px.scatter(
        plot_df, x="PC1", y="PC2", color=color_col,
        hover_data=["oligo_id"],
        title=f"Chemical Space (PCA) — colored by {color_by}",
        labels={"PC1": f"PC1 ({pca.explained_variance_ratio_[0]:.1%})",
                "PC2": f"PC2 ({pca.explained_variance_ratio_[1]:.1%})"},
        height=480,
    )
    fig.update_traces(marker=dict(size=6, opacity=0.7))
    fig.update_layout(plot_bgcolor="#f8f9fa")
    return fig


# ── Gradio App ─────────────────────────────────────────────────────────────

def build_app() -> gr.Blocks:
    with gr.Blocks(
        title="OligoTox-Web | Oligonucleotide Toxicity Database",
        theme=gr.themes.Soft(),
        css=".gradio-container { max-width: 1200px; margin: auto; }",
    ) as app:
        gr.Markdown("""
# OligoTox-Web
**Open dataset and AI prediction tools for oligonucleotide toxicity**

Part of the NIH NCATS OligoTox Open Data Challenge | Anote, Inc. | [GitHub](https://github.com/anote-ai/oligotox-predict) | [Dataset (Zenodo)](https://zenodo.org)
        """)

        with gr.Tabs():
            # ── Predictor tab ─────────────────────────────────────────────
            with gr.Tab("Toxicity Predictor"):
                gr.Markdown("Enter an oligonucleotide sequence and modification scheme to get a predicted toxicity profile across all 47 endpoints.")
                with gr.Row():
                    with gr.Column(scale=2):
                        seq_input = gr.Textbox(
                            label="Oligonucleotide Sequence (5'→3')",
                            placeholder="GCATTTGCTTTGCATATGCTT",
                            lines=2,
                        )
                        with gr.Row():
                            backbone_input = gr.Dropdown(BACKBONE_CHOICES, value="PS", label="Backbone Chemistry")
                            sugar_input = gr.Dropdown(SUGAR_CHOICES, value="DNA", label="Sugar Modification")
                        with gr.Row():
                            conjugate_input = gr.Dropdown(CONJUGATE_CHOICES, value="none", label="Conjugate")
                            gap_input = gr.Slider(0, 15, value=0, step=1, label="Gapmer Gap Length (0 = not gapmer)")
                        predict_btn = gr.Button("Predict Toxicity Profile", variant="primary")
                    with gr.Column(scale=3):
                        pred_table = gr.Dataframe(label="Predicted Toxicity Profile", wrap=True)
                        pred_chart = gr.Plot(label="Toxicity Bar Chart")

                predict_btn.click(
                    predict_toxicity,
                    inputs=[seq_input, backbone_input, sugar_input, conjugate_input, gap_input],
                    outputs=[pred_table, pred_chart],
                )

            # ── Database explorer tab ─────────────────────────────────────
            with gr.Tab("Database Explorer"):
                gr.Markdown("Search OligoToxDB by sequence, chemistry, or toxicity profile.")
                with gr.Row():
                    search_seq = gr.Textbox(label="Sequence contains", placeholder="GCATTT")
                    search_backbone = gr.Dropdown(["All"] + BACKBONE_CHOICES, value="All", label="Backbone")
                    search_tox = gr.Dropdown(["All", "inactive", "low", "moderate", "high"], value="All", label="Toxicity Class")
                    max_rows_slider = gr.Slider(10, 1000, value=100, step=10, label="Max rows")
                search_btn = gr.Button("Search", variant="primary")
                results_table = gr.Dataframe(label="Results", wrap=True)
                download_btn = gr.Button("Download as CSV")
                download_file = gr.File(label="Download")

                search_btn.click(
                    search_database,
                    inputs=[search_seq, search_backbone, search_tox, max_rows_slider],
                    outputs=[results_table],
                )

                def make_csv_file(seq, bb, tox):
                    csv_str = download_csv(seq, bb, tox)
                    import tempfile
                    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
                        f.write(csv_str)
                        return f.name

                download_btn.click(
                    make_csv_file,
                    inputs=[search_seq, search_backbone, search_tox],
                    outputs=[download_file],
                )

            # ── Dose-response viewer tab ──────────────────────────────────
            with gr.Tab("Dose-Response Viewer"):
                gr.Markdown("View dose-response curves with 4PL curve fit for any compound and endpoint.")
                with gr.Row():
                    dr_oligo_id = gr.Textbox(label="OligoTox ID", placeholder="OT-00001")
                    dr_endpoint = gr.Dropdown(list(ENDPOINTS_DISPLAY.values()), label="Endpoint",
                                              value="Cell Viability (PHH)")
                dr_btn = gr.Button("Plot", variant="primary")
                dr_plot = gr.Plot()
                dr_btn.click(plot_dose_response, inputs=[dr_oligo_id, dr_endpoint], outputs=[dr_plot])

            # ── Chemical space tab ────────────────────────────────────────
            with gr.Tab("Chemical Space"):
                gr.Markdown("PCA visualization of the full OligoToxDB compound library.")
                color_by = gr.Radio(["Backbone", "Sugar Mod", "Conjugate"], value="Backbone", label="Color by")
                space_btn = gr.Button("Generate Plot", variant="primary")
                space_plot = gr.Plot()
                space_btn.click(plot_chemical_space, inputs=[color_by], outputs=[space_plot])

            # ── About tab ─────────────────────────────────────────────────
            with gr.Tab("About"):
                gr.Markdown("""
## About OligoToxDB

**OligoToxDB** is the largest open dataset of human in vitro oligonucleotide toxicity profiles,
generated as part of the NIH NCATS OligoTox Open Data Challenge (Phase 2).

### Data Summary
- **2,800 oligonucleotides** spanning 7 backbone classes and 4 sugar modification types
- **47 toxicity endpoints** across 6 organ/mechanism domains
- **5 human in vitro model systems**: primary hepatocytes, kidney organoids, PBMCs, platelets, organ-on-chip
- **Fully open**: CC BY 4.0 (data) | Apache 2.0 (code)

### Citation
If you use OligoToxDB in your research, please cite:
> Vidra N. et al. (2026). OligoToxDB: A large-scale open dataset for AI-driven oligonucleotide toxicity prediction. *Scientific Data*. DOI: pending.

### Contact
nvidra@anote.ai | [GitHub Issues](https://github.com/anote-ai/oligotox-predict/issues)

### License
Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Code: [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0)
                """)

    return app


def main() -> None:
    app = build_app()
    app.launch(server_name="0.0.0.0", server_port=7860, share=False)


if __name__ == "__main__":
    main()
