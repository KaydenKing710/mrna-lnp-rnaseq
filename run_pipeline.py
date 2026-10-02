from pathlib import Path
import shutil
import subprocess
import sys
import time

# PROJECT CONFIGURATION
PROJECT_ROOT = Path(__file__).resolve().parent


PYTHON_STEPS = [
    (
        "01 - Preprocessing",
        "analysis/01_preprocessing/prepare_counts.py",
    ),
    (
        "02 - Quality Control",
        "analysis/02_quality_control/qc_analysis.py",
    ),
    (
        "03 - PyDESeq2 Differential Expression",
        "analysis/03_pydeseq2/run_pydeseq2.py",
    ),
    (
        "04 - PyDESeq2 Volcano Plot",
        "analysis/03_pydeseq2/volcano_plot.py",
    ),
    (
        "05 - Gene Annotation",
        "analysis/04_biological_analysis/annotate_genes.py",
    ),
    (
        "06 - Selected-Gene Heatmap",
        "analysis/04_biological_analysis/heatmap_selected_genes.py",
    ),
    (
        "07 - GO Enrichment",
        "analysis/04_biological_analysis/go_enrichment.py",
    ),
]


COMPARISON_STEPS = [
    (
        "09 - DEG Set Comparison",
        "analysis/06_method_comparison/compare_deg_sets.py",
    ),
    (
        "10 - Fold-Change Comparison",
        "analysis/06_method_comparison/logfc_scatter.py",
    ),
    (
        "11 - Discordant DEG Analysis",
        "analysis/06_method_comparison/analyze_discordant_genes.py",
    ),
]


EDGER_SCRIPT = (
    "analysis/05_edger/run_edger.R"
)

# HELPER FUNCTIONS
def print_header(title):

    print("\n")
    print("=" * 70)
    print(title)
    print("=" * 70)


def run_python_step(title, relative_path):

    print_header(title)

    script_path = (
        PROJECT_ROOT
        / relative_path
    )

    if not script_path.exists():

        raise FileNotFoundError(
            f"Script not found: "
            f"{script_path}"
        )

    subprocess.run(
        [
            sys.executable,
            str(script_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )


def find_rscript():

    # First try PATH
    rscript = shutil.which(
        "Rscript"
    )

    if rscript:
        return rscript


    # Common Windows locations
    possible_paths = [
        Path(
            r"C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
        ),
        Path(
            r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
        ),
        Path(
            r"C:\Program Files\R\R-4.5.3\bin\Rscript.exe"
        ),
    ]

    for path in possible_paths:

        if path.exists():
            return str(path)


    return None


def run_edger():

    print_header(
        "08 - edgeR Differential Expression"
    )

    rscript = find_rscript()

    if rscript is None:

        raise RuntimeError(
            "Rscript was not found. "
            "Install R or add Rscript "
            "to the system PATH."
        )


    script_path = (
        PROJECT_ROOT
        / EDGER_SCRIPT
    )

    if not script_path.exists():

        raise FileNotFoundError(
            f"edgeR script not found: "
            f"{script_path}"
        )


    subprocess.run(
        [
            rscript,
            "--vanilla",
            str(script_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )

# MAIN PIPELINE
def main():

    start_time = time.time()

    print_header(
        "mRNA-LNP RNA-seq Analysis Pipeline"
    )

    print(
        f"Project directory:\n"
        f"{PROJECT_ROOT}"
    )

    print(
        f"\nPython executable:\n"
        f"{sys.executable}"
    )

    # Python analysis before edgeR
    for title, script in PYTHON_STEPS:

        run_python_step(
            title,
            script,
        )

    # edgeR


    run_edger()


    # Method comparison


    for title, script in COMPARISON_STEPS:

        run_python_step(
            title,
            script,
        )

    # Final summary
    elapsed = (
        time.time()
        - start_time
    )

    print_header(
        "PIPELINE COMPLETED SUCCESSFULLY"
    )

    print(
        f"Total runtime: "
        f"{elapsed / 60:.2f} minutes"
    )

    print(
        "\nMain output directories:"
    )

    print(
        "results/figures/"
    )

    print(
        "results/tables/"
    )

# ENTRY POINT
if __name__ == "__main__":

    try:

        main()

    except subprocess.CalledProcessError as error:

        print_header(
            "PIPELINE FAILED"
        )

        print(
            "A pipeline step returned "
            "an error."
        )

        print(
            f"Exit code: "
            f"{error.returncode}"
        )

        sys.exit(
            error.returncode
        )

    except Exception as error:

        print_header(
            "PIPELINE FAILED"
        )

        print(
            f"{type(error).__name__}: "
            f"{error}"
        )

        sys.exit(1)