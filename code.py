import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

st.set_page_config(page_title="Molecular Info", page_icon="🧪", layout="wide")

st.title("Molecule Information")
st.write("Enter a SMILES code to get detailed molecular information")

# Example SMILES with descriptions
with st.expander("See example SMILES codes"):
    st.write("""
    - Aspirin: `CC(=O)OC1=CC=CC=C1C(=O)O`
    - Caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
    - Paracetamol: `CC(=O)NC1=CC=C(O)C=C1`
    - Ibuprofen: `CC(C)CC1=CC=C(C=C1)[C@H](C)C(=O)O`
    """)

smiles = st.text_input("Enter SMILES code:", "")

if smiles:
    try:
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            st.error("Invalid SMILES code. Please try again.")
        else:
            # Calculate properties
            mw = Descriptors.ExactMolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            rings = Descriptors.RingCount(mol)
            aromatic_rings = sum(1 for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
            
            # Create columns for layout
            col1, col2 = st.columns(2)
            
            # Basic Properties
            with col1:
                st.write("### Basic Properties")
                st.write(f"🏋️ Molecular Weight: {mw:.2f} g/mol")
                st.write(f"💧 LogP: {logp:.2f}")
                st.write(f"🔵 H-Bond Donors: {hbd}")
                st.write(f"🔴 H-Bond Acceptors: {hba}")
                st.write(f"📊 TPSA: {tpsa:.2f} Å²")
                st.write(f"🔄 Rotatable Bonds: {rotatable_bonds}")
                st.write(f"⭕ Total Rings: {rings}")
                st.write(f"💫 Aromatic Rings: {aromatic_rings}")

            # Drug-likeness Analysis
            with col2:
                st.write("### Drug-likeness Analysis")
                
                # Lipinski's Rule of 5
                ro5_violations = sum([
                    mw > 500,
                    logp > 5,
                    hbd > 5,
                    hba > 10
                ])
                
                st.write("#### Lipinski's Rule of 5")
                violations = []
                if mw > 500: violations.append("Molecular weight > 500")
                if logp > 5: violations.append("LogP > 5")
                if hbd > 5: violations.append("H-bond donors > 5")
                if hba > 10: violations.append("H-bond acceptors > 10")
                
                if ro5_violations == 0:
                    st.success("✅ Passes all Lipinski's Rules!")
                elif ro5_violations == 1:
                    st.warning(f"⚠️ Has {ro5_violations} violation: {violations[0]}")
                else:
                    st.error(f"❌ Has {ro5_violations} violations:")
                    for v in violations:
                        st.write(f"- {v}")
                
                # Veber's Rules
                st.write("#### Veber's Rules")
                veber_violations = []
                if rotatable_bonds > 10:
                    veber_violations.append("Rotatable bonds > 10")
                if tpsa > 140:
                    veber_violations.append("TPSA > 140")
                
                if not veber_violations:
                    st.success("✅ Passes Veber's Rules!")
                else:
                    st.error("❌ Violations:")
                    for v in veber_violations:
                        st.write(f"- {v}")
                
                # Lead-likeness
                st.write("#### Lead-likeness")
                if 200 <= mw <= 350 and -1 <= logp <= 3:
                    st.success("✅ Compound is in lead-like space!")
                else:
                    st.warning("⚠️ Compound is outside lead-like space")
            
            # Chemical Formula
            formula = AllChem.rdMolDescriptors.CalcMolFormula(mol)
            st.write(f"### Chemical Formula: {formula}")

    except Exception as e:
        st.error(f"An error occurred while analyzing the molecule: {str(e)}")
        st.write("Please check your SMILES code and try again.")
