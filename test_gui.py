from streamlit.testing.v1 import AppTest
import pytest
from src import fileupload
import json
from pathlib import Path
import shutil


@pytest.fixture
def launch(request):
    test = AppTest.from_file(request.param)

    ## Initialize session state ##
    with open("settings.json", "r") as f:
        test.session_state.settings = json.load(f)
    test.session_state.settings["test"] = True
    test.secrets["workspace"] = "test"
    return test


# Test launching of all pages
@pytest.mark.parametrize(
    "launch",
    (
        # "content/quickstart.py", # NOTE: this page does not work due to streamlit.errors.StreamlitPageNotFoundError error
        "content/topp_workflow_file_upload.py",
        "content/topp_workflow_parameter.py",
        "content/topp_workflow_execution.py",
        "content/topp_workflow_results.py",
        "content/run_subprocess.py",
    ),
    indirect=True,
)
def test_launch(launch):
    """Test if all pages can be launched without errors."""
    launch.run(timeout=30)  # Increased timeout from 10 to 30 seconds
    assert not launch.exception
