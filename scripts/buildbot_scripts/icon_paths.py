from pathlib import Path

buildbot_script_path = Path(__file__).parent.absolute()
buildbot_list_path = (buildbot_script_path / "experiment_lists").absolute()
base_path = buildbot_script_path.parents[1].absolute()
run_path = (base_path / "run").absolute()
