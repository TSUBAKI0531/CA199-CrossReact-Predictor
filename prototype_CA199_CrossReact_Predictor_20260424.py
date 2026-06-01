import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def analyze_sequence(sequence):
    try:
        analyzed_seq = ProteinAnalysis(sequence.replace(" ", "").replace("\n", "").upper())
        mw = analyzed_seq.molecular_weight()
        pi = analyzed_seq.isoelectric_point()
        aromatic = analyzed_seq.aromaticity()
        return {"分子量 (MW)": round(mw, 2), "等電点 (pI)": round(pi, 2), "芳香族性 (Aromaticity)": round(aromatic, 3)}
    except Exception as e:
        return None

def mock_predict_cross_reactivity(seq_props):
    # プロトタイプ用のモック予測ロジック
    base_ca199 = 85.0
    base_slex = 30.0
    base_lea = 20.0
    
    # 芳香族性が高いとCA19-9への親和性が上がるが交差反応性も増えるという仮のロジック
    aroma = seq_props.get("芳香族性 (Aromaticity)", 0.1)
    
    score_ca199 = min(100, base_ca199 + (aroma * 50))
    score_slex = min(100, base_slex + (aroma * 150))
    score_lea = min(100, base_lea + (aroma * 80))
    
    return {
        "CA19-9 (sLe^a)": score_ca199,
        "sLe^x (交差リスク)": score_slex,
        "Le^a (交差リスク)": score_lea
    }

st.set_page_config(page_title="CA19-9 Cross-Reactivity Predictor", layout="wide")

st.title("CA19-9 Cross-Reactivity Predictor (Prototype)")
st.markdown("CA19-9抗体の交差反応性・オンターゲット毒性リスクをスクリーニングするツールです。(`GlycoAntibodyStudio` 連携モジュール)")

st.sidebar.header("1. 抗体配列入力")
vh_seq = st.sidebar.text_area("VH Sequences (FASTA format or raw text)", value="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK")

if st.sidebar.button("予測を実行"):
    if not vh_seq:
        st.warning("配列を入力してください。")
    else:
        st.subheader("分析結果")
        
        props = analyze_sequence(vh_seq)
        if not props:
            st.error("配列の解析に失敗しました。正しいアミノ酸配列を入力してください。")
        else:
            col1, col2 = st.columns(2)
            with col1:
                st.write("**物理化学的プロパティ (Protein Hydrophobicity Profiler準拠)**")
                st.json(props)
            
            predictions = mock_predict_cross_reactivity(props)
            
            with col2:
                st.write("**交差反応性予測スコア**")
                fig, ax = plt.subplots(figsize=(6, 4))
                labels = list(predictions.keys())
                scores = list(predictions.values())
                colors = ['#2ca02c' if 'CA19-9' in l else '#d62728' for l in labels]
                ax.bar(labels, scores, color=colors)
                ax.set_ylim(0, 100)
                ax.set_ylabel("Predicted Binding Score")
                st.pyplot(fig)
            
            st.subheader("リスク評価サマリー")
            if predictions["sLe^x (交差リスク)"] > 50:
                st.error("⚠️ 警告: sLe^x に対する強い交差反応性が予測されます。正常組織（特に好中球等）へのオンターゲット毒性リスクがあります。配列の最適化を検討してください。")
            else:
                st.success("✅ 交差反応性リスクは低水準に抑えられています。")
            
            st.info("人手確認事項: このスコアはIn silicoのプロトタイプ予測値です。最終的な開発候補への選定にはSPR等を用いたin vitroの交差反応性検証が必要です。")
