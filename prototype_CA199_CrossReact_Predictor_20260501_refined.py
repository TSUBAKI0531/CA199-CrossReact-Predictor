import logging
from typing import Dict, Optional

import matplotlib.figure
import matplotlib.pyplot as plt
import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def analyze_sequence(sequence: str) -> Optional[Dict[str, float]]:
    """VH配列からタンパク質物性を算出する。

    Args:
        sequence: アミノ酸一文字表記のVH配列。空白・改行を含んでもよい。

    Returns:
        分子量(MW)、等電点(pI)、芳香族性(Aromaticity)を格納したdict。
        入力が空または無効なアミノ酸を含む場合はNone。
    """
    cleaned = sequence.replace(" ", "").replace("\n", "").upper()
    if not cleaned:
        logger.warning("空の配列が入力されました")
        return None
    try:
        analyzed_seq = ProteinAnalysis(cleaned)
        mw = analyzed_seq.molecular_weight()
        pi = analyzed_seq.isoelectric_point()
        aromatic = analyzed_seq.aromaticity()
        logger.info("配列解析完了: MW=%.2f, pI=%.2f, Aromaticity=%.3f", mw, pi, aromatic)
        return {
            "MW": round(mw, 2),
            "pI": round(pi, 2),
            "Aromaticity": round(aromatic, 3),
        }
    except ValueError as e:
        logger.error("無効なアミノ酸配列です: %s", e)
        return None
    except Exception as e:
        logger.error("配列解析中に予期しないエラーが発生しました: %s", e)
        return None


def mock_predict(seq_props: Dict[str, float]) -> Dict[str, float]:
    """配列物性から各ルイス抗原エピトープへの交差反応性スコアをモック予測する。

    芳香族性(Aromaticity)を主要なパラメータとしてスコアを算出する。
    実装はモックであり、実際の結合親和性モデルへの置き換えを前提とする。

    Args:
        seq_props: analyze_sequence() が返すdict。"Aromaticity" キーを使用する。

    Returns:
        エピトープ名をキー、予測スコア(0〜100)を値とするdict。
        キーは "CA19-9 (sLe^a)", "sLe^x", "Le^a"。
    """
    aroma = seq_props.get("Aromaticity", 0.1)
    scores: Dict[str, float] = {
        "CA19-9 (sLe^a)": min(100.0, 85 + aroma * 50),
        "sLe^x": min(100.0, 30 + aroma * 150),
        "Le^a": min(100.0, 20 + aroma * 80),
    }
    logger.info("交差反応性スコア算出完了: %s", scores)
    return scores


def render_chart(preds: Dict[str, float]) -> matplotlib.figure.Figure:
    """交差反応性スコアを棒グラフとして描画する。

    CA19-9 (sLe^a) を緑、その他のルイス抗原を赤で表示する。

    Args:
        preds: mock_predict() が返すスコアdict。

    Returns:
        描画済みのmatplotlib Figureオブジェクト。
    """
    colors = ["green" if k == "CA19-9 (sLe^a)" else "red" for k in preds]
    fig, ax = plt.subplots()
    ax.bar(preds.keys(), preds.values(), color=colors)
    ax.set_ylabel("交差反応性スコア (%)")
    ax.set_ylim(0, 110)
    ax.set_title("Cross-Reactivity Profile")
    logger.info("グラフ描画完了")
    return fig


def main() -> None:
    """CA19-9 Cross-Reactivity PredictorのStreamlitアプリエントリポイント。"""
    st.title("CA19-9 Cross-Reactivity Predictor")
    st.write("連携: GlycoAntibodyStudio, Protein Hydrophobicity Profiler")

    seq = st.text_area("VH Sequence", "EVQLVESGGGLVQPGGSLRLSCAASGFTFS")

    if st.button("予測"):
        logger.info("予測ボタンが押されました")

        props = analyze_sequence(seq)
        if props is None:
            st.error(
                "配列の解析に失敗しました。有効なアミノ酸一文字配列を入力してください。"
            )
            return

        st.subheader("配列物性")
        st.json(props)

        preds = mock_predict(props)

        st.subheader("交差反応性スコア")
        fig = render_chart(preds)
        st.pyplot(fig)

        if preds["sLe^x"] > 50:
            st.error(
                "⚠ 交差反応性リスク高: sLe^x スコアが閾値(50)を超えています。"
                " 配列の再設計を検討してください。"
            )
        else:
            st.success("交差反応性リスクは許容範囲内です。")


main()
