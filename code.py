import streamlit as st

st.set_page_config(page_title="Molecular Info", page_icon="🧪")

st.title("Molecule Information")
st.write("Basic test to verify Streamlit is working")

test_input = st.text_input("Type something to test:", "")
if test_input:
    st.write(f"You typed: {test_input}")
