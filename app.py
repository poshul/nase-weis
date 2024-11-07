import streamlit as st
from pathlib import Path

if __name__ == '__main__':
    pages = [
        st.Page(Path("content", "quickstart.py"), title="Quickstart", icon="👋"),
        st.Page(Path("content", "naseweis.py"), title="NASEWEIS", icon="👃"),
    ]

    pg = st.navigation(pages)
    pg.run()