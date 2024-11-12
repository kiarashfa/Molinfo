import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

st.set_page_config(page_title="Molecular Info", page_icon="🧪", layout="wide")

st.title("Molecule Information")
st.write("Enter a SMILES code to get molecular information")

# Create two columns for the layout
col1, col2 = st.columns(2)

with col1:
    # Example SMILES with expandable section
    with st.expander("See example SMILES codes"):
        st.write("**Aspirin:** CC(=O)OC1=CC=CC=C1C(=O)O")
        st.write("**Caffeine:** CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        st.write("**Paracetamol:** CC(=O)NC1=CC=C(O)C=C1")

    # Input section
    smiles = st.text_input("Enter SMILES code:", "")

if smiles:
    try:
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            st.error("Invalid SMILES code. Please try again.")
        else:
            # Calculate properties
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            
            # Display molecule image
            with col1:
                try:
                    img = Draw.MolToImage(mol)
                    st.image(img, caption='Molecule Structure', use_column_width=True)
                except Exception as e:
                    st.warning("Could not generate molecule image")
            
            # Display properties
            with col2:
                st.write("### Molecular Properties")
                
                # Create a more organized display of properties
                properties = {
                    "Molar Mass": f"{mw:.2f} g/mol",
                    "LogP": f"{logp:.2f}",
                    "H-Bond Donors": str(hbd),
                    "H-Bond Acceptors": str(hba),
                    "Topological Polar Surface Area": f"{tpsa:.2f} Å²"
                }
                
                for prop, value in properties.items():
                    st.write(f"**{prop}:** {value}")

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")

# Add some information at the bottom
st.markdown("""---
This app uses RDKit to analyze molecular properties from SMILES codes.
- Molar Mass: The mass of one mole of the substance
- LogP: Measure of lipophilicity (higher values = more lipophilic)
- H-Bond Donors/Acceptors: Number of hydrogen bond donors and acceptors
- TPSA: Topological Polar Surface Area, useful for drug-likeness prediction
""")
