import re
import argparse
from pathlib import Path
import sys

TYPE_HINT_REPLACEMENTS = {
    "calc.constants.": "constants.",
    "calc.interp.": "interp.",
    "calc.layer.": "layer.",
    "calc.params.": "params.",
    "calc.parcel.": "parcel.",
    "calc.thermo.": "thermo.",
    "calc.winds.": "winds.",
}

IMPORT_REPLACEMENTS = {
    "import calc.constants": "import constants",
    "import calc.interp": "import .interp",
    "import calc.layer": "import .layer",
    "import calc.params": "import .params",
    "import calc.parcel": "import .parcel",
    "import calc.thermo": "import .thermo",
    "import calc.winds": "import .winds",
}

def fix_stubs(file_path: Path):
    """Reads a stub file and performs a two-stage fix on its content."""
    try:
        original_content = file_path.read_text(encoding='utf-8')
        modified_content = original_content

        for old_hint, new_hint in TYPE_HINT_REPLACEMENTS.items():
            modified_content = modified_content.replace(old_hint, new_hint)

        for bad_import, good_import in IMPORT_REPLACEMENTS.items():
            pattern = re.compile(f"^{re.escape(bad_import)}$", re.MULTILINE)
            modified_content = pattern.sub(good_import, modified_content)

        if original_content != modified_content:
            file_path.write_text(modified_content, encoding='utf-8')
            print(f"   ... Patched stubs in {file_path.name}")

    except Exception as e:
        print(f"Error processing {file_path}: {e}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="A script to fix circular imports in auto-generated nanobind stub files."
    )
    parser.add_argument(
        "stubs_directory",
        help="The root directory where the .pyi stubs were generated."
    )
    args = parser.parse_args()

    stubs_dir = Path(args.stubs_directory)
    if not stubs_dir.is_dir():
        print(f"Error: The specified directory does not exist: {stubs_dir}", file=sys.stderr)
        sys.exit(1)
        
    print(f"Scanning for stub files in: {stubs_dir.resolve()}")

    found_files = list(stubs_dir.glob("*.pyi"))
    if not found_files:
        print("Warning: No submodule '__init__.pyi' files were found to patch.")
        
    for pyi_file in found_files:
        fix_stubs(pyi_file)

    print("\nStub import patching complete.")


if __name__ == "__main__":
    main()
