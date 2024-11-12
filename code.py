import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

# Configure the page
st.set_page_config(
    page_title="Molecular Info",
    page_icon="🧪",
    layout="wide"
)

def main():
    st.title("Molecule Information")
    st.write("Enter a SMILES code to get molecular information")

    # Example SMILES
    with st.expander("See example SMILES codes"):
        st.code("Aspirin: CC(=O)OC1=CC=CC=C1C(=O)O")
        st.code("Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    smiles = st.text_input("Enter SMILES code:", "")

    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            
            if mol is None:
                st.error("Invalid SMILES code. Please check your input and try again.")
            else:
                # Calculate properties
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                hbd = Descriptors.NumHDonors(mol)
                hba = Descriptors.NumHAcceptors(mol)
                tpsa = Descriptors.TPSA(mol)
                
                # Create two columns
                col1, col2 = st.columns(2)
                
                # Display molecule image
                with col1:
                    try:
                        img = Draw.MolToImage(mol)
                        st.image(img, caption='Molecule Structure')
                    except Exception as e:
                        st.warning("Could not generate molecule image")
                
                # Display properties
                with col2:
                    st.write("### Molecular Properties")
                    st.write(f"🏋️ Molar Mass: {mw:.2f} g/mol")
                    st.write(f"💧 LogP: {logp:.2f}")
                    st.write(f"🔵 H-Bond Donors: {hbd}")
                    st.write(f"🔴 H-Bond Acceptors: {hba}")
                    st.write(f"📊 TPSA: {tpsa:.2f} Å²")

        except Exception as e:
            st.error(f"An error occurred: {str(e)}")
            st.write("Please check your SMILES code and try again.")

if __name__ == "__main__":
    main()
