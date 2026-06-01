# SPEC NOTE [矛盾①]: 仕様書「保存形式」に CSV/JSON ダウンロード機能が要求されているが未実装。
# SPEC NOTE [矛盾②]: 仕様書「判定結果 Pass/Alert/Fail の3段階」に対し、本実装は Pass/Alert の2段階のみ。
import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

st.set_page_config(
    page_title="CA199 CrossReact Predictor",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("🔬 CA199 CrossReact Predictor")
st.markdown("##### **Theme**: CA19-9 glycan antibody | **Date**: 2026-06-01")
st.write("**ポートフォリオ連携**: 「GlycoAntibodyStudio」の一部として統合し、物理化学的性質算出エンジンには「Protein Hydrophobicity Profiler」を連携。「CA199-CrossReact-Predictor」としての単体展開も可能。")

st.sidebar.header("⚙️ パラメータ設定")
threshold = st.sidebar.slider("アラート警告しきい値", min_value=0.0, max_value=1.0, value=0.5, step=0.05)

st.subheader("1. 解析対象アミノ酸配列の入力")
default_seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDGYYGMDVWGQGTTVTVSS"
sequence_input = st.text_area("VH/VL または CDR3 配列 (FASTA / Plain Text)", default_seq, height=120)


def analyze_sequence(sequence: str) -> tuple[dict[str, float | int] | None, str | None]:
    """アミノ酸配列を解析して物理化学的特性を返す。

    Args:
        sequence: 解析対象のアミノ酸配列文字列（FASTA またはプレーンテキスト）。

    Returns:
        (properties, error) のタプル。
        解析成功時は properties に各特性値の辞書、error は None。
        失敗時は properties が None、error にエラーメッセージ文字列を返す。

    Raises:
        なし（例外は内部でキャッチしてエラーメッセージとして返す）。
    """
    logger.info("配列解析を開始します（入力文字数: %d）", len(sequence))
    cleaned_seq = "".join([c for c in sequence.upper() if c.isalpha()])
    valid_aa = "ACDEFGHIKLMNPQRSTVWY"
    for c in cleaned_seq:
        if c not in valid_aa:
            msg = f"無効なアミノ酸記号が含まれています: {c}"
            logger.error(msg)
            return None, msg

    if len(cleaned_seq) < 10:
        msg = "配列が短すぎます（最低10アミノ酸以上入力してください）"
        logger.error(msg)
        return None, msg

    try:
        analysis = ProteinAnalysis(cleaned_seq)
        mw = analysis.molecular_weight()
        pi = analysis.isoelectric_point()
        aromatic = analysis.aromaticity()
        gravy = analysis.gravy()
        logger.info(
            "配列解析完了: length=%d, MW=%.2f, pI=%.2f, Aromaticity=%.3f, GRAVY=%.3f",
            len(cleaned_seq), mw, pi, aromatic, gravy,
        )
        return {
            "Sequence Length": len(cleaned_seq),
            "Molecular Weight (Da)": round(mw, 2),
            "Isoelectric Point (pI)": round(pi, 2),
            "Aromaticity": round(aromatic, 3),
            "GRAVY (Hydrophobicity)": round(gravy, 3),
        }, None
    except ValueError as e:
        msg = f"配列値エラー: {e}"
        logger.exception("ProteinAnalysis で ValueError が発生しました")
        return None, msg
    except Exception as e:
        msg = f"解析エラー: {e}"
        logger.exception("ProteinAnalysis で予期しないエラーが発生しました")
        return None, msg


def calculate_predictions(
    properties: dict[str, float | int],
    theme_key: str,
) -> tuple[dict[str, float], str]:
    """物理化学的特性からテーマ別の予測スコアを算出する。

    Args:
        properties: analyze_sequence() が返す特性値辞書。
            "Aromaticity" と "GRAVY (Hydrophobicity)" キーを使用する。
        theme_key: 予測モデルを選択するテーマキー。
            {"ca19_9", "intrabody", "deep_generative_antibody",
             "rnaseq_antibody", "tcell_generative_ai", "cart"} のいずれか。

    Returns:
        (predictions, alert_metric) のタプル。
        predictions: 各スコア名をキー、スコア値を値とする辞書。
        alert_metric: アラート判定に使用するスコア名。
    """
    aroma = properties.get("Aromaticity", 0.1)
    gravy = properties.get("GRAVY (Hydrophobicity)", 0.0)
    logger.info(
        "予測計算を開始します (theme_key=%s, aromaticity=%.3f, gravy=%.3f)",
        theme_key, aroma, gravy,
    )

    np.random.seed(42)

    if theme_key == "ca19_9":
        sLea_bind = min(1.0, 0.8 + aroma * 0.4)
        sLex_bind = min(1.0, 0.1 + (gravy + 2) * 0.2 + aroma * 0.5)
        Lea_bind = min(1.0, 0.2 + (gravy + 2) * 0.15)
        preds: dict[str, float] = {
            "CA19-9 (sLe^a) binding score": round(sLea_bind, 3),
            "sLe^x cross-reactivity score": round(sLex_bind, 3),
            "Le^a cross-reactivity score": round(Lea_bind, 3),
        }
        alert_metric = "sLe^x cross-reactivity score"
    elif theme_key == "intrabody":
        folding = min(1.0, max(0.0, 0.9 - (gravy + 1) * 0.2))
        aggregation = min(1.0, max(0.0, 0.1 + (gravy + 1.5) * 0.3))
        half_life = min(48.0, max(0.5, 24.0 - gravy * 12.0))
        preds = {
            "Cytosolic folding stability": round(folding, 3),
            "Aggregation risk index": round(aggregation, 3),
            "Estimated intracellular half-life (h)": round(half_life, 2),
        }
        alert_metric = "Aggregation risk index"
    elif theme_key == "deep_generative_antibody":
        solubility = min(1.0, max(0.0, 0.85 - (gravy + 1.0) * 0.25))
        expression = min(10.0, max(0.1, 5.0 - gravy * 2.0 + aroma * 5.0))
        yield_score = min(1.0, solubility * (expression / 10.0))
        preds = {
            "Calculated solubility index": round(solubility, 3),
            "Estimated expression yield (g/L)": round(expression, 2),
            "Overall developability score": round(yield_score, 3),
        }
        alert_metric = "Calculated solubility index"
    elif theme_key == "rnaseq_antibody":
        ST_coloc = min(1.0, max(0.0, 0.4 + aroma * 1.2 - gravy * 0.2))
        ST_enrich = min(1.0, max(0.0, 0.2 + (gravy + 1.5) * 0.3))
        preds = {
            "Spatial colocalization score": round(ST_coloc, 3),
            "Target infiltration index": round(ST_enrich, 3),
        }
        alert_metric = "Spatial colocalization score"
    elif theme_key == "tcell_generative_ai":
        self_pMHC = min(1.0, max(0.0, 0.15 + (gravy + 2) * 0.15 + aroma * 0.8))
        target_pMHC = min(1.0, max(0.0, 0.85 + aroma * 0.3))
        preds = {
            "Target pMHC binding affinity": round(target_pMHC, 3),
            "Self-pMHC cross-reactivity risk": round(self_pMHC, 3),
        }
        alert_metric = "Self-pMHC cross-reactivity risk"
    else:  # cart
        hdr_ratio = min(1.0, max(0.0, 0.3 + aroma * 1.5 - gravy * 0.3))
        nhej_ratio = 1.0 - hdr_ratio
        ki_efficiency = min(1.0, hdr_ratio * 0.9)
        preds = {
            "Predicted HDR repair efficiency": round(hdr_ratio, 3),
            "Predicted NHEJ repair efficiency": round(nhej_ratio, 3),
            "Predicted Knock-in success rate": round(ki_efficiency, 3),
        }
        alert_metric = "Predicted NHEJ repair efficiency"

    logger.info("予測計算完了: %s", preds)
    return preds, alert_metric


if st.button("🧬 解析 & 予測の実行", type="primary"):
    with st.spinner("配列特性を抽出中..."):
        props, err = analyze_sequence(sequence_input)

    if err:
        st.error(err)
    else:
        st.success("アミノ酸配列の分析が完了しました。")

        st.subheader("2. 物理化学的特性解析 (BioPython)")
        cols = st.columns(len(props))
        for col, (k, v) in zip(cols, props.items()):
            col.metric(label=k, value=v)

        preds, alert_metric = calculate_predictions(props, "ca19_9")

        st.subheader("3. 研究ツール・特性シミュレーション予測")
        col_left, col_right = st.columns(2)

        with col_left:
            st.write("**予測スコア一覧**")
            st.json(preds)

            if not preds:
                logger.error("予測結果が空です")
                st.error("予測結果が空です。入力を確認してください。")
            else:
                alert_val = preds.get(alert_metric, 0.0)
                is_alert = (
                    alert_val > threshold
                    if "risk" in alert_metric or "cross" in alert_metric or "NHEJ" in alert_metric
                    else alert_val < (1.0 - threshold)
                )
                if is_alert:
                    logger.warning(
                        "アラート閾値超過: %s = %.3f (threshold=%.2f)",
                        alert_metric, alert_val, threshold,
                    )
                    st.error(f"⚠️ 警告: {alert_metric} ({alert_val}) が基準値（しきい値: {threshold}）の安全圏から外れています！")
                else:
                    st.success("✅ 安全: 特性スコアは基準範囲内です。")

        with col_right:
            st.write("**予測プロファイル可視化**")
            fig, ax = plt.subplots(figsize=(6, 4))
            bars = ax.barh(list(preds.keys()), list(preds.values()), color="#1f77b4")
            ax.set_xlim(0, max(1.1, max(preds.values())))
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            for bar in bars:
                width = bar.get_width()
                ax.text(
                    width + 0.02, bar.get_y() + bar.get_height() / 2,
                    f"{width}",
                    va='center', ha='left', fontsize=9, fontweight='bold',
                )

            st.pyplot(fig)
            plt.close(fig)
