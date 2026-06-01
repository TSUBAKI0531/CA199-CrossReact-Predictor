import logging
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
)
logger = logging.getLogger(__name__)


def analyze_sequence(sequence: str) -> dict[str, float] | None:
    """アミノ酸配列の物理化学的プロパティを解析する。

    BioPythonの ProteinAnalysis を用いて、分子量・等電点・芳香族性を計算する。
    入力配列は空白・改行を除去し、大文字に正規化してから解析する。

    Args:
        sequence: 解析対象のアミノ酸配列。FASTA形式のヘッダーを含まない
            生の配列文字列、またはスペース・改行を含む文字列を受け付ける。

    Returns:
        解析結果を格納した辞書。キーは "分子量 (MW)"、"等電点 (pI)"、
        "芳香族性 (Aromaticity)" で、値はそれぞれ丸め済みの float。
        解析に失敗した場合は None を返す。

    Raises:
        解析中に発生した例外はキャッチしてログに記録し、呼び出し元には
        None を返すことで通知する（例外は再送出しない）。
    """
    logger.info("配列解析を開始します (長さ: %d文字)", len(sequence))
    try:
        cleaned = sequence.replace(" ", "").replace("\n", "").upper()
        analyzed_seq = ProteinAnalysis(cleaned)
        mw = analyzed_seq.molecular_weight()
        pi = analyzed_seq.isoelectric_point()
        aromatic = analyzed_seq.aromaticity()
        result = {
            "分子量 (MW)": round(mw, 2),
            "等電点 (pI)": round(pi, 2),
            "芳香族性 (Aromaticity)": round(aromatic, 3),
        }
        logger.info("配列解析成功: MW=%.2f, pI=%.2f, Aromaticity=%.3f", mw, pi, aromatic)
        return result
    except Exception as e:
        logger.error("配列解析に失敗しました: %s", e, exc_info=True)
        return None


def mock_predict_cross_reactivity(seq_props: dict[str, float]) -> dict[str, float]:
    """物理化学的プロパティから交差反応性スコアをモック予測する。

    プロトタイプ用の簡易予測ロジック。芳香族性（Aromaticity）を主要因として
    CA19-9・sLe^x・Le^a に対する結合スコアを算出する。
    芳香族性が高いほど CA19-9 への親和性が上昇するが、同時に交差反応性も
    増加するという仮定に基づいている。

    Args:
        seq_props: ``analyze_sequence`` が返す物理化学的プロパティの辞書。
            "芳香族性 (Aromaticity)" キーの値を使用する。キーが存在しない場合は
            デフォルト値 0.1 を使用する。

    Returns:
        各ターゲットの予測結合スコア（0〜100）を格納した辞書。
        キーは "CA19-9 (sLe^a)"、"sLe^x (交差リスク)"、"Le^a (交差リスク)"。

    Raises:
        ValueError: スコア計算中に数値変換が失敗した場合。
    """
    logger.info("交差反応性モック予測を開始します")
    try:
        base_ca199 = 85.0
        base_slex = 30.0
        base_lea = 20.0

        aroma = seq_props.get("芳香族性 (Aromaticity)", 0.1)

        score_ca199 = min(100, base_ca199 + (aroma * 50))
        score_slex = min(100, base_slex + (aroma * 150))
        score_lea = min(100, base_lea + (aroma * 80))

        predictions = {
            "CA19-9 (sLe^a)": score_ca199,
            "sLe^x (交差リスク)": score_slex,
            "Le^a (交差リスク)": score_lea,
        }
        logger.info(
            "予測スコア: CA19-9=%.1f, sLe^x=%.1f, Le^a=%.1f",
            score_ca199,
            score_slex,
            score_lea,
        )
        return predictions
    except ValueError as e:
        logger.error("スコア計算中に数値エラーが発生しました: %s", e, exc_info=True)
        raise


st.set_page_config(page_title="CA19-9 Cross-Reactivity Predictor", layout="wide")

st.title("CA19-9 Cross-Reactivity Predictor (Prototype)")
st.markdown(
    "CA19-9抗体の交差反応性・オンターゲット毒性リスクをスクリーニングするツールです。"
    "(`GlycoAntibodyStudio` 連携モジュール)"
)

st.sidebar.header("1. 抗体配列入力")
vh_seq = st.sidebar.text_area(
    "VH Sequences (FASTA format or raw text)",
    value="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK",
)

if st.sidebar.button("予測を実行"):
    logger.info("予測実行ボタンが押されました")
    if not vh_seq.strip():
        st.warning("配列を入力してください。")
        logger.warning("空の配列が入力されました")
    else:
        st.subheader("分析結果")

        props = analyze_sequence(vh_seq)
        if props is None:
            st.error(
                "配列の解析に失敗しました。正しい1文字表記のアミノ酸配列を入力してください。"
                "（例: EVQLVES... ）\n\n"
                "詳細はターミナルのログを確認してください。"
            )
        else:
            col1, col2 = st.columns(2)
            with col1:
                st.write("**物理化学的プロパティ (Protein Hydrophobicity Profiler準拠)**")
                st.json(props)

            try:
                predictions = mock_predict_cross_reactivity(props)
            except ValueError:
                st.error("交差反応性スコアの計算中にエラーが発生しました。ログを確認してください。")
                st.stop()

            with col2:
                st.write("**交差反応性予測スコア**")
                fig, ax = plt.subplots(figsize=(6, 4))
                labels = list(predictions.keys())
                scores = list(predictions.values())
                colors = ["#2ca02c" if "CA19-9" in l else "#d62728" for l in labels]
                ax.bar(labels, scores, color=colors)
                ax.set_ylim(0, 100)
                ax.set_ylabel("Predicted Binding Score")
                st.pyplot(fig)

            st.subheader("リスク評価サマリー")
            if predictions["sLe^x (交差リスク)"] > 50:
                st.error(
                    "⚠️ 警告: sLe^x に対する強い交差反応性が予測されます。"
                    "正常組織（特に好中球等）へのオンターゲット毒性リスクがあります。"
                    "配列の最適化を検討してください。"
                )
                logger.warning(
                    "高い交差反応性リスクを検出: sLe^x スコア=%.1f",
                    predictions["sLe^x (交差リスク)"],
                )
            else:
                st.success("✅ 交差反応性リスクは低水準に抑えられています。")
                logger.info(
                    "交差反応性リスク低: sLe^x スコア=%.1f",
                    predictions["sLe^x (交差リスク)"],
                )

            st.info(
                "人手確認事項: このスコアはIn silicoのプロトタイプ予測値です。"
                "最終的な開発候補への選定にはSPR等を用いたin vitroの交差反応性検証が必要です。"
            )
