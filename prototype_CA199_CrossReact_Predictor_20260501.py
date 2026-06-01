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
        return {"MW": round(mw, 2), "pI": round(pi, 2), "Aromaticity": round(aromatic, 3)}
    except:
        return None

def mock_predict(seq_props):
    aroma = seq_props.get("Aromaticity", 0.1)
    return {"CA19-9 (sLe^a)": min(100, 85 + aroma*50), "sLe^x": min(100, 30 + aroma*150), "Le^a": min(100, 20 + aroma*80)}

st.title("CA19-9 Cross-Reactivity Predictor")
st.write("連携: GlycoAntibodyStudio, Protein Hydrophobicity Profiler")
seq = st.text_area("VH Sequence", "EVQLVESGGGLVQPGGSLRLSCAASGFTFS")
if st.button("予測"):
    props = analyze_sequence(seq)
    if props:
        st.json(props)
        preds = mock_predict(props)
        fig, ax = plt.subplots()
        ax.bar(preds.keys(), preds.values(), color=['green', 'red', 'red'])
        st.pyplot(fig)
        if preds["sLe^x"] > 50: st.error("交差反応性リスク高")
