import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Set Page Config
st.set_page_config(
    page_title="CA199 CrossReact Predictor",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Header Section
st.title("🔬 CA199 CrossReact Predictor")
st.markdown("##### **Theme**: CA19-9 glycan antibody | **Date**: 2026-06-01")
st.write("**ポートフォリオ連携**: 「GlycoAntibodyStudio」の一部として統合し、物理化学的性質算出エンジンには「Protein Hydrophobicity Profiler」を連携。「CA199-CrossReact-Predictor」としての単体展開も可能。")

# Sidebar Settings
st.sidebar.header("⚙️ パラメータ設定")
threshold = st.sidebar.slider("アラート警告しきい値", min_value=0.0, max_value=1.0, value=0.5, step=0.05)

# Input Sequence Area
st.subheader("1. 解析対象アミノ酸配列の入力")
default_seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDGYYGMDVWGQGTTVTVSS"
sequence_input = st.text_area("VH/VL または CDR3 配列 (FASTA / Plain Text)", default_seq, height=120)

# Analyzer logic
def analyze_sequence(sequence):
    cleaned_seq = "".join([c for c in sequence.upper() if c.isalpha()])
    # Validate amino acids
    valid_aa = "ACDEFGHIKLMNPQRSTVWY"
    for c in cleaned_seq:
        if c not in valid_aa:
            return None, f"無効なアミノ酸記号が含まれています: {c}"
            
    if len(cleaned_seq) < 10:
        return None, "配列が短すぎます（最低10アミノ酸以上入力してください）"
        
    try:
        analysis = ProteinAnalysis(cleaned_seq)
        mw = analysis.molecular_weight()
        pi = analysis.isoelectric_point()
        aromatic = analysis.aromaticity()
        
        # Calculate hydrophobicity index (GRAVY)
        gravy = analysis.gravy()
        
        return {
            "Sequence Length": len(cleaned_seq),
            "Molecular Weight (Da)": round(mw, 2),
            "Isoelectric Point (pI)": round(pi, 2),
            "Aromaticity": round(aromatic, 3),
            "GRAVY (Hydrophobicity)": round(gravy, 3)
        }, None
    except Exception as e:
        return None, f"解析エラー: {str(e)}"

def calculate_predictions(properties, theme_key):
    # Simulated prediction scoring based on biological features (Aromaticity/Hydrophobicity)
    aroma = properties.get("Aromaticity", 0.1)
    gravy = properties.get("GRAVY (Hydrophobicity)", 0.0)
    
    np.random.seed(42)  # Seed for deterministic mock results
    
    if theme_key == "ca19_9":
        # Predict relative binding affinities for CA19-9 vs similar Lewis antigens
        sLea_bind = min(1.0, 0.8 + aroma * 0.4)
        sLex_bind = min(1.0, 0.1 + (gravy + 2) * 0.2 + aroma * 0.5)
        Lea_bind = min(1.0, 0.2 + (gravy + 2) * 0.15)
        return {
            "CA19-9 (sLe^a) binding score": round(sLea_bind, 3),
            "sLe^x cross-reactivity score": round(sLex_bind, 3),
            "Le^a cross-reactivity score": round(Lea_bind, 3)
        }, "sLe^x cross-reactivity score"
        
    elif theme_key == "intrabody":
        # Predict cytosolic folding stability and aggregation risk
        folding = min(1.0, max(0.0, 0.9 - (gravy + 1) * 0.2))
        aggregation = min(1.0, max(0.0, 0.1 + (gravy + 1.5) * 0.3))
        half_life = min(48.0, max(0.5, 24.0 - gravy * 12.0))
        return {
            "Cytosolic folding stability": round(folding, 3),
            "Aggregation risk index": round(aggregation, 3),
            "Estimated intracellular half-life (h)": round(half_life, 2)
        }, "Aggregation risk index"
        
    elif theme_key == "deep_generative_antibody":
        # Predict manufacturing developability criteria
        solubility = min(1.0, max(0.0, 0.85 - (gravy + 1.0) * 0.25))
        expression = min(10.0, max(0.1, 5.0 - gravy * 2.0 + aroma * 5.0))
        yield_score = min(1.0, solubility * (expression / 10.0))
        return {
            "Calculated solubility index": round(solubility, 3),
            "Estimated expression yield (g/L)": round(expression, 2),
            "Overall developability score": round(yield_score, 3)
        }, "Calculated solubility index"
        
    elif theme_key == "rnaseq_antibody":
        # Predict spatial colocalization score and tumor infiltration rate
        ST_coloc = min(1.0, max(0.0, 0.4 + aroma * 1.2 - gravy * 0.2))
        ST_enrich = min(1.0, max(0.0, 0.2 + (gravy + 1.5) * 0.3))
        return {
            "Spatial colocalization score": round(ST_coloc, 3),
            "Target infiltration index": round(ST_enrich, 3)
        }, "Spatial colocalization score"
        
    elif theme_key == "tcell_generative_ai":
        # Predict TCR-pMHC off-target risk
        self_pMHC = min(1.0, max(0.0, 0.15 + (gravy + 2) * 0.15 + aroma * 0.8))
        target_pMHC = min(1.0, max(0.0, 0.85 + aroma * 0.3))
        return {
            "Target pMHC binding affinity": round(target_pMHC, 3),
            "Self-pMHC cross-reactivity risk": round(self_pMHC, 3)
        }, "Self-pMHC cross-reactivity risk"
        
    else:  # cart
        # Predict HDR vs NHEJ editing repair outcome ratio
        hdr_ratio = min(1.0, max(0.0, 0.3 + aroma * 1.5 - gravy * 0.3))
        nhej_ratio = 1.0 - hdr_ratio
        ki_efficiency = min(1.0, hdr_ratio * 0.9)
        return {
            "Predicted HDR repair efficiency": round(hdr_ratio, 3),
            "Predicted NHEJ repair efficiency": round(nhej_ratio, 3),
            "Predicted Knock-in success rate": round(ki_efficiency, 3)
        }, "Predicted NHEJ repair efficiency"

if st.button("🧬 解析 & 予測の実行", type="primary"):
    with st.spinner("配列特性を抽出中..."):
        props, err = analyze_sequence(sequence_input)
        
    if err:
        st.error(err)
    else:
        st.success("アミノ酸配列の分析が完了しました。")
        
        # Display sequence physical properties
        st.subheader("2. 物理化学的特性解析 (BioPython)")
        cols = st.columns(len(props))
        for col, (k, v) in zip(cols, props.items()):
            col.metric(label=k, value=v)
            
        # Run mock prediction
        preds, alert_metric = calculate_predictions(props, "ca19_9")
        
        st.subheader("3. 研究ツール・特性シミュレーション予測")
        col_left, col_right = st.columns(2)
        
        with col_left:
            st.write("**予測スコア一覧**")
            st.json(preds)
            
            # Check alert threshold
            alert_val = preds.get(alert_metric, 0.0)
            if alert_val > threshold if "risk" in alert_metric or "cross" in alert_metric or "NHEJ" in alert_metric else alert_val < (1.0 - threshold):
                st.error(f"⚠️ 警告: {alert_metric} ({alert_val}) が基準値（しきい値: {threshold}）の安全圏から外れています！")
            else:
                st.success(f"✅ 安全: 特性スコアは基準範囲内です。")
                
        with col_right:
            # Draw plot
            st.write("**予測プロファイル可視化**")
            fig, ax = plt.subplots(figsize=(6, 4))
            bars = ax.barh(list(preds.keys()), list(preds.values()), color="#1f77b4")
            ax.set_xlim(0, max(1.1, max(preds.values())))
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            # Add value labels
            for bar in bars:
                width = bar.get_width()
                ax.text(width + 0.02, bar.get_y() + bar.get_height()/2, 
                        f"{width}", 
                        va='center', ha='left', fontsize=9, fontweight='bold')
                        
            st.pyplot(fig)
