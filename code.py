import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors

st.set_page_config(page_title="Molecular Info", page_icon="🧪")

st.title("Molecule Information")

smiles = st.text_input("Enter SMILES code:", "")

if smiles:
    try:
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            st.error("Invalid SMILES code. Please try again.")
        else:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            st.write(f"Molar Mass: {mw:.2f} g/mol")
            st.write(f"LogP: {logp:.2f}")
            st.write(f"Hydrogen Bond Donors: {hbd}")
            st.write(f"Hydrogen Bond Acceptors: {hba}")
    except Exception as e:
        st.error(f"Error: {e}")
