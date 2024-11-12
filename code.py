import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import IPythonConsole

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
            
            # Display properties first (moving this up since we might have issues with visualization)
            with col2:
                st.write("### Molecular Properties")
                
                properties = {
                    "Molar Mass": f"{mw:.2f} g/mol",
                    "LogP": f"{logp:.2f}",
                    "H-Bond Donors": str(hbd),
                    "H-Bond Acceptors": str(hba),
                    "Topological Polar Surface Area": f"{tpsa:.2f} Å²"
                }
                
                for prop, value in properties.items():
                    st.write(f"**{prop}:** {value}")
                
                # Additional calculated properties
                ro5_violations = sum([
                    mw > 500,
                    logp > 5,
                    hbd > 5,
                    hba > 10
                ])
                
                st.write("### Drug-likeness")
                st.write(f"Lipinski Rule of 5 violations: {ro5_violations}")
                if ro5_violations <= 1:
                    st.success("Molecule passes Lipinski's Rule of 5! ✅")
                else:
                    st.warning("Molecule violates Lipinski's Rule of 5")

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")

# Add information at the bottom
st.markdown("""---
### About the Properties
- **Molar Mass:** The mass of one mole of the substance (g/mol)
- **LogP:** Measure of lipophilicity (higher values = more lipophilic)
- **H-Bond Donors/Acceptors:** Number of hydrogen bond donors and acceptors
- **TPSA:** Topological Polar Surface Area, useful for drug-likeness prediction

### Lipinski's Rule of 5
A molecule is more likely to be orally active if it meets these criteria:
- Molecular weight ≤ 500
- LogP ≤ 5
- H-bond donors ≤ 5
- H-bond acceptors ≤ 10
""")
