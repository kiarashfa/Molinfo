import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor

st.set_page_config(page_title="Molecular Info", page_icon="🧪", layout="wide")

def mol_to_svg(mol):
    """Convert molecule to SVG string"""
    rdDepictor.Compute2DCoords(mol)
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(400, 400)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg

st.title("Molecule Information")
st.write("Enter a SMILES code to get molecular information")

# Example SMILES with descriptions
with st.expander("See example SMILES codes"):
    st.write("🌟 **Common Molecules:**")
    st.code("Aspirin: CC(=O)OC1=CC=CC=C1C(=O)O")
    st.code("Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    st.code("Paracetamol: CC(=O)NC1=CC=C(O)C=C1")

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
            
            # Create two columns
            col1, col2 = st.columns(2)
            
            # Display molecule SVG
            with col1:
                st.write("### Molecule Structure")
                svg = mol_to_svg(mol)
                st.write(svg, unsafe_allow_html=True)
            
            # Display properties
            with col2:
                st.write("### Molecular Properties")
                st.write(f"🏋️ Molar Mass: {mw:.2f} g/mol")
                st.write(f"💧 LogP: {logp:.2f}")
                st.write(f"🔵 H-Bond Donors: {hbd}")
                st.write(f"🔴 H-Bond Acceptors: {hba}")
                st.write(f"📊 TPSA: {tpsa:.2f} Å²")
                
                # Calculate Lipinski's Rule of 5 violations
                ro5_violations = sum([
                    mw > 500,
                    logp > 5,
                    hbd > 5,
                    hba > 10
                ])
                
                st.write("\n### Drug-likeness")
                st.write(f"Lipinski Rule of 5 violations: {ro5_violations}")
                if ro5_violations <= 1:
                    st.success("Molecule passes Lipinski's Rule of 5! ✅")
                else:
                    st.warning("Molecule violates Lipinski's Rule of 5")

    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        st.write("Please check your SMILES code and try again.")
