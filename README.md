# CA199 CrossReact Predictor

**テーマ**: CA19-9 glycan antibody（CA19-9グリカン抗体）

## 概要

CA19-9グリカン抗体の研究において、正常組織に発現する類似ルイス抗原（sLe^x, Le^a等）との交差反応性をアミノ酸配列から予測するツールです。BioPythonによる物理化学的特性解析とシミュレーションモデルにより、設計段階でのオンターゲット毒性リスクをスクリーニングします。結果はPass/Alert 2段階で判定し、棒グラフとともに可視化します。

## 入力

- アミノ酸配列（VH/VL または CDR3）
- 形式: FASTAまたはプレーンテキスト

## 出力

| スコア | 説明 |
|--------|------|
| CA19-9 (sLe^a) binding score | CA19-9抗原に対する結合スコア |
| sLe^x cross-reactivity score | sLe^x抗原への交差反応リスクスコア |
| Le^a cross-reactivity score | Le^a抗原への交差反応リスクスコア |

判定: **Pass / Alert** 2段階（アラート閾値はサイドバーで調整可能）

## 使用技術

- Python
- Streamlit
- BioPython（ProteinAnalysis: MW / pI / Aromaticity / GRAVY 算出）

## ローカル起動方法

```bash
pip install -r requirements.txt
streamlit run app.py
```
