# 仕様書: CA19-9 Cross-Reactivity Predictor

## 1. 概要
CA19-9抗体開発において、正常ルイス抗原（sLe^x, Le^aなど）との交差反応性を予測・スクリーニングするStreamlitアプリ。

## 2. 目的
- 抗体の配列（FASTA）等から結合性スコアを比較し、特異性の高いリード抗体を選別する。
- `GlycoAntibodyStudio` エコシステムの一部として機能させる。

## 3. 主要機能
- FASTA配列入力（VH/VL）
- BioPythonを用いた配列特性抽出と交差反応性モックスコアリング
- PandasとStreamlitを用いた結合プロファイルの可視化とアラート

## 4. 技術スタック
- Python, Streamlit, pandas, BioPython, matplotlib
