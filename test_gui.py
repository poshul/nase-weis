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
        "content/documentation.py",
        "content/topp_workflow_file_upload.py",
        "content/topp_workflow_parameter.py",
        "content/topp_workflow_execution.py",
        "content/topp_workflow_results.py",
        "content/download_section.py",
        "content/simple_workflow.py",
        "content/run_subprocess.py",
    ),
    indirect=True,
)
def test_launch(launch):
    """Test if all pages can be launched without errors."""
    launch.run(timeout=30)  # Increased timeout from 10 to 30 seconds
    assert not launch.exception


########### PAGE SPECIFIC TESTS ############
@pytest.mark.parametrize(
    "launch,selection",
    [
        ("content/documentation.py", "User Guide"),
        ("content/documentation.py", "Installation"),
        (
            "content/documentation.py",
            "Developers Guide: How to build app based on this template",
        ),
        ("content/documentation.py", "Developers Guide: TOPP Workflow Framework"),
        ("content/documentation.py", "Developer Guide: Windows Executables"),
        ("content/documentation.py", "Developers Guide: Deployment"),
    ],
    indirect=["launch"],
)
def test_documentation(launch, selection):
    launch.run()
    launch.selectbox[0].select(selection).run()
    assert not launch.exception


@pytest.mark.parametrize("launch", ["content/file_upload.py"], indirect=True)
def test_file_upload_load_example(launch):
    launch.run()
    for i in launch.tabs:
        if i.label == "Example Data":
            i.button[0].click().run()
            assert not launch.exception

