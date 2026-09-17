import argparse
import sys
from pathlib import Path
from streamlit.web import cli as stcli

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default=None, help="Path to .yaml config file")
    args = parser.parse_args()

    main_script = Path(__file__).parent/"RiboMARTA_main.py"
    command = ["streamlit", "run", str(main_script)]
    if args.config is not None:
        command += ["--", "--config", args.config]
    sys.argv = command

    stcli.main()

if __name__ == "__main__":
    main()