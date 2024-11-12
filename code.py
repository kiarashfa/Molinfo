import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors

st.set_page_config(page_title="Molecular Info", page_icon="🧪")

st.title("Molecule Information")
st.write("Enter a SMILES code to get molecular information")

# Example SMILES
st.write("Example: Aspirin = CC(=O)OC1=CC=CC=C1C(=O)O")

smiles = st.text_input("Enter SMILES code:", "")

if smiles:
    try:
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            st.error("Invalid SMILES code. Please try again.")
        else:
            # Calculate basic properties
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            
            # Display properties
            st.write("### Molecular Properties")
            st.write(f"Molar Mass: {mw:.2f} g/mol")
            st.write(f"LogP: {logp:.2f}")

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
