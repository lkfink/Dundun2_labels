"""
Spectral quality-control analysis for the dùndún stimuli.
For AP&P Revision
Lauren Fink + chatGPT 5.6 Sol

Outputs
-------
1. spectral_metrics.csv
2. individual_spectra/<stimulus>.png
3. stimulus_spectra_all.pdf
4. stimulus_spectra_heatmap.png/pdf
5. stimulus_spectral_rolloff.png/pdf

Edit AUDIO_DIR below, then run:

    python plot_stimulus_spectra.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import librosa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import welch


# =============================================================================
# Configuration
# =============================================================================

AUDIO_DIR = Path("/Users/laurenfink/Documents/myProjects/Dundun/fromMPIEA_laptop/dundun/Madita_stimuli/Final_recordings_wav")

OUTPUT_DIR = AUDIO_DIR / "spectral_qc"
INDIVIDUAL_DIR = OUTPUT_DIR / "individual_spectra"

AUDIO_EXTENSIONS = {
    ".wav",
    ".mp3",
    ".flac",
    ".ogg",
    ".aif",
    ".aiff",
    ".m4a",
}

# Upper frequency shown in the plots.
MAX_PLOT_HZ = 20_000

# Frequencies specifically mentioned by the reviewer.
REFERENCE_FREQUENCIES_HZ = (7_000, 15_000)

# Display range for spectra normalized to each recording's maximum.
MIN_RELATIVE_DB = -100
MAX_RELATIVE_DB = 3

# Welch parameters.
TARGET_SEGMENT_LENGTH = 16_384
TARGET_NFFT = 32_768

# Heatmap frequency spacing.
HEATMAP_STEP_HZ = 25

# "filename", "rolloff_95", or "rolloff_99"
HEATMAP_SORT = "filename"
HEATMAP_SORT = "power_above_7khz"


# =============================================================================
# Helper functions
# =============================================================================

def find_audio_files(directory: Path) -> list[Path]:
    """Find supported audio files recursively."""

    if not directory.exists():
        raise FileNotFoundError(
            f"Audio directory does not exist:\n{directory}"
        )

    files = sorted(
        path
        for path in directory.rglob("*")
        if path.is_file()
        and path.suffix.lower() in AUDIO_EXTENSIONS
        and OUTPUT_DIR not in path.parents
    )

    if not files:
        raise FileNotFoundError(
            f"No supported audio files were found in:\n{directory}"
        )

    return files


def load_audio(path: Path) -> tuple[np.ndarray, int]:
    """
    Load at the file's native sample rate and mix channels to mono.
    """

    try:
        signal, sample_rate = librosa.load(
            path,
            sr=None,
            mono=True,
            dtype=np.float64,
        )
    except Exception as exc:
        raise RuntimeError(
            f"Could not decode {path.name}.\n"
            "For MP3 files, make sure FFmpeg is installed."
        ) from exc

    signal = np.asarray(signal, dtype=np.float64)
    signal = signal[np.isfinite(signal)]

    if signal.size < 256:
        raise ValueError(
            f"{path.name} contains too few usable samples."
        )

    # Remove any constant DC offset.
    signal = signal - np.mean(signal)

    if np.max(np.abs(signal)) == 0:
        raise ValueError(
            f"{path.name} appears to contain silence only."
        )

    return signal, int(sample_rate)


def calculate_spectrum(
    signal: np.ndarray,
    sample_rate: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Estimate the one-sided power spectral density using Welch's method.

    Returns
    -------
    frequencies:
        Frequency values in Hz.
    psd:
        Linear power spectral density.
    relative_db:
        PSD normalized so that the strongest frequency bin is 0 dB.
    """

    segment_length = min(
        TARGET_SEGMENT_LENGTH,
        signal.size,
    )

    # Use a power-of-two FFT at least as long as the Welch segment.
    minimum_nfft = 2 ** math.ceil(
        math.log2(segment_length)
    )

    nfft = max(
        TARGET_NFFT,
        minimum_nfft,
    )

    frequencies, psd = welch(
        signal,
        fs=sample_rate,
        window="hann",
        nperseg=segment_length,
        noverlap=segment_length // 2,
        nfft=nfft,
        detrend="constant",
        return_onesided=True,
        scaling="density",
    )

    psd = np.maximum(
        psd,
        np.finfo(float).tiny,
    )

    relative_db = 10 * np.log10(
        psd / np.max(psd)
    )

    return frequencies, psd, relative_db


def spectral_rolloff(
    frequencies: np.ndarray,
    psd: np.ndarray,
    proportion: float,
) -> float:
    """
    Frequency below which the requested proportion of total power occurs.
    """

    if not 0 < proportion < 1:
        raise ValueError(
            "Rolloff proportion must be between zero and one."
        )

    cumulative_power = np.cumsum(psd)
    total_power = cumulative_power[-1]

    if total_power <= 0:
        return float("nan")

    target = proportion * total_power
    index = int(
        np.searchsorted(
            cumulative_power,
            target,
            side="left",
        )
    )

    index = min(index, frequencies.size - 1)

    return float(frequencies[index])


def power_proportion_above(
    frequencies: np.ndarray,
    psd: np.ndarray,
    threshold_hz: float,
) -> float:
    """Proportion of total spectral power above a frequency."""

    total_power = float(np.sum(psd))

    if total_power <= 0:
        return float("nan")

    high_frequency_power = float(
        np.sum(psd[frequencies >= threshold_hz])
    )

    return high_frequency_power / total_power


def save_individual_spectrum(
    file_name: str,
    frequencies: np.ndarray,
    relative_db: np.ndarray,
    metrics: dict[str, float | str],
    output_path: Path,
) -> plt.Figure:
    """Create and save an individual spectrum figure."""

    max_frequency = min(
        MAX_PLOT_HZ,
        float(frequencies[-1]),
    )

    figure, axis = plt.subplots(
        figsize=(8, 4.5)
    )

    axis.plot(
        frequencies / 1000,
        relative_db,
        linewidth=0.8,
    )

    # for reference_hz in REFERENCE_FREQUENCIES_HZ:
    #     if reference_hz <= max_frequency:
    #         axis.axvline(
    #             reference_hz / 1000,
    #             linestyle="--",
    #             linewidth=0.8,
    #             alpha=0.65,
    #         )

    axis.axvline(
        float(metrics["rolloff_95_hz"]) / 1000,
        linestyle=":",
        linewidth=1,
        alpha=0.9,
        label="95% roll-off",
    )

    axis.axvline(
        float(metrics["rolloff_99_hz"]) / 1000,
        linestyle="-.",
        linewidth=1,
        alpha=0.9,
        label="99% roll-off",
    )

    axis.set_xlim(
        0,
        max_frequency / 1000,
    )

    axis.set_ylim(
        MIN_RELATIVE_DB,
        MAX_RELATIVE_DB,
    )

    axis.set_xlabel("Frequency (kHz)")
    axis.set_ylabel("Power spectral density (dB relative to peak)")
    axis.set_title(file_name)

    summary_text = (
        f"95% roll-off: "
        f"{float(metrics['rolloff_95_hz']) / 1000:.2f} kHz\n"
        f"99% roll-off: "
        f"{float(metrics['rolloff_99_hz']) / 1000:.2f} kHz\n"
        f"Power >7 kHz: "
        f"{100 * float(metrics['power_above_7khz']):.3f}%\n"
        f"Power >15 kHz: "
        f"{100 * float(metrics['power_above_15khz']):.3f}%"
    )

    axis.text(
        0.985,
        0.96,
        summary_text,
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.85,
        },
    )

    axis.legend(
        loc="lower left",
        frameon=False,
        fontsize=8,
    )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    figure.tight_layout()

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    return figure


def make_heatmap(
    spectra: list[dict],
    metrics: pd.DataFrame,
) -> None:
    """
    Create a normalized spectral heatmap across all stimuli.
    """

    maximum_common_hz = min(
        MAX_PLOT_HZ,
        max(item["frequencies"][-1] for item in spectra),
    )

    common_frequencies = np.arange(
        0,
        maximum_common_hz + HEATMAP_STEP_HZ,
        HEATMAP_STEP_HZ,
    )

    rows = []

    for item in spectra:
        interpolated = np.interp(
            common_frequencies,
            item["frequencies"],
            item["relative_db"],
            left=np.nan,
            right=np.nan,
        )

        rows.append(interpolated)

    matrix = np.vstack(rows)
    names = np.array(
        [item["file_name"] for item in spectra]
    )

    # if HEATMAP_SORT == "rolloff_95":
    #     order = np.argsort(
    #         metrics["rolloff_95_hz"].to_numpy()
    #     )
    # elif HEATMAP_SORT == "rolloff_99":
    #     order = np.argsort(
    #         metrics["rolloff_99_hz"].to_numpy()
    #     )
    # else:
    #     order = np.argsort(
    #         np.char.lower(names.astype(str))
    #     )


    if HEATMAP_SORT == "power_above_7khz":
        order = np.argsort(
            metrics["percent_power_above_7khz"].to_numpy()
        )
    elif HEATMAP_SORT == "rolloff_95":
        order = np.argsort(
            metrics["rolloff_95_hz"].to_numpy()
        )
    elif HEATMAP_SORT == "rolloff_99":
        order = np.argsort(
            metrics["rolloff_99_hz"].to_numpy()
        )
    else:
        order = np.argsort(
            np.char.lower(names.astype(str))
        )

    matrix = matrix[order]
    names = names[order]

    figure_height = max(
        7,
        0.28 * len(names) + 2,
    )

    figure, axis = plt.subplots(
        figsize=(11, figure_height)
    )

    image = axis.imshow(
        matrix,
        aspect="auto",
        origin="upper",
        interpolation="nearest",
        extent=[
            common_frequencies[0] / 1000,
            common_frequencies[-1] / 1000,
            len(names) - 0.5,
            -0.5,
        ],
        vmin=MIN_RELATIVE_DB,
        vmax=0,
        cmap="viridis",
    )

    axis.set_yticks(
        np.arange(len(names))
    )
    axis.set_yticklabels(
        names,
        fontsize=8,
    )

    for reference_hz in REFERENCE_FREQUENCIES_HZ:
        if reference_hz <= maximum_common_hz:
            axis.axvline(
                reference_hz / 1000,
                linestyle="--",
                linewidth=1,
            )

    axis.set_xlabel("Frequency (kHz)")
    axis.set_ylabel("Stimulus")
    axis.set_title(
        "Long-term average spectra across stimuli"
    )

    colorbar = figure.colorbar(
        image,
        ax=axis,
        pad=0.015,
    )
    colorbar.set_label(
        "Power spectral density (dB relative to each stimulus peak)"
    )

    figure.tight_layout()

    for suffix in ("png", "pdf"):
        figure.savefig(
            OUTPUT_DIR / f"stimulus_spectra_heatmap.{suffix}",
            dpi=300 if suffix == "png" else None,
            bbox_inches="tight",
            facecolor="white",
        )

    plt.close(figure)


def make_rolloff_plot(metrics: pd.DataFrame) -> None:
    """Plot 95% and 99% roll-off frequencies for every stimulus."""

    plot_data = metrics.sort_values(
        "rolloff_95_hz",
        ascending=True,
    ).reset_index(drop=True)

    y_positions = np.arange(
        len(plot_data)
    )

    figure_height = max(
        6,
        0.27 * len(plot_data) + 2,
    )

    figure, axis = plt.subplots(
        figsize=(9, figure_height)
    )

    axis.scatter(
        plot_data["rolloff_95_hz"] / 1000,
        y_positions,
        marker="o",
        label="95% roll-off",
    )

    axis.scatter(
        plot_data["rolloff_99_hz"] / 1000,
        y_positions,
        marker="x",
        label="99% roll-off",
    )

    for y_position, row in plot_data.iterrows():
        axis.plot(
            [
                row["rolloff_95_hz"] / 1000,
                row["rolloff_99_hz"] / 1000,
            ],
            [y_position, y_position],
            linewidth=0.6,
            alpha=0.5,
        )

    for reference_hz in REFERENCE_FREQUENCIES_HZ:
        axis.axvline(
            reference_hz / 1000,
            linestyle="--",
            linewidth=1,
            alpha=0.75,
        )

    axis.set_yticks(y_positions)
    axis.set_yticklabels(
        plot_data["file_name"],
        fontsize=8,
    )

    axis.set_xlabel("Spectral roll-off frequency (kHz)")
    axis.set_ylabel("Stimulus")
    axis.set_title("Spectral bandwidth of experimental stimuli")

    axis.legend(
        loc="lower right",
        frameon=False,
    )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    figure.tight_layout()

    for suffix in ("png", "pdf"):
        figure.savefig(
            OUTPUT_DIR / f"stimulus_spectral_rolloff.{suffix}",
            dpi=300 if suffix == "png" else None,
            bbox_inches="tight",
            facecolor="white",
        )

    plt.close(figure)


# =============================================================================
# Main analysis
# =============================================================================

def main() -> None:
    """Run the complete spectral QC analysis."""

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    INDIVIDUAL_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    audio_files = find_audio_files(
        AUDIO_DIR
    )

    print(
        f"Found {len(audio_files)} audio files."
    )

    all_spectra: list[dict] = []
    metrics_rows: list[dict] = []
    failed_files: list[tuple[str, str]] = []

    pdf_path = OUTPUT_DIR / "stimulus_spectra_all.pdf"

    with PdfPages(pdf_path) as pdf:
        for file_number, path in enumerate(
            audio_files,
            start=1,
        ):
            print(
                f"[{file_number}/{len(audio_files)}] "
                f"Processing {path.name}"
            )

            try:
                signal, sample_rate = load_audio(path)

                frequencies, psd, relative_db = (
                    calculate_spectrum(
                        signal,
                        sample_rate,
                    )
                )

                rolloff_95 = spectral_rolloff(
                    frequencies,
                    psd,
                    0.95,
                )

                rolloff_99 = spectral_rolloff(
                    frequencies,
                    psd,
                    0.99,
                )

                power_above_7khz = power_proportion_above(
                    frequencies,
                    psd,
                    7_000,
                )

                power_above_15khz = power_proportion_above(
                    frequencies,
                    psd,
                    15_000,
                )

                row = {
                    "file_name": path.name,
                    "file_path": str(path),
                    "sample_rate_hz": sample_rate,
                    "duration_seconds": signal.size / sample_rate,
                    "nyquist_hz": sample_rate / 2,
                    "rolloff_95_hz": rolloff_95,
                    "rolloff_99_hz": rolloff_99,
                    "power_above_7khz": power_above_7khz,
                    "percent_power_above_7khz": (
                        100 * power_above_7khz
                    ),
                    "power_above_15khz": power_above_15khz,
                    "percent_power_above_15khz": (
                        100 * power_above_15khz
                    ),
                }

                metrics_rows.append(row)

                all_spectra.append(
                    {
                        "file_name": path.name,
                        "frequencies": frequencies,
                        "psd": psd,
                        "relative_db": relative_db,
                    }
                )

                safe_stem = path.stem.replace(
                    "/",
                    "_",
                )

                figure = save_individual_spectrum(
                    file_name=path.name,
                    frequencies=frequencies,
                    relative_db=relative_db,
                    metrics=row,
                    output_path=(
                        INDIVIDUAL_DIR
                        / f"{safe_stem}_spectrum.png"
                    ),
                )

                pdf.savefig(
                    figure,
                    bbox_inches="tight",
                    facecolor="white",
                )

                plt.close(figure)

            except Exception as exc:
                failed_files.append(
                    (path.name, str(exc))
                )
                print(
                    f"WARNING: skipped {path.name}: {exc}",
                    file=sys.stderr,
                )

    if not metrics_rows:
        raise RuntimeError(
            "No audio files were analyzed successfully."
        )

    metrics = pd.DataFrame(
        metrics_rows
    ).sort_values(
        "file_name",
        key=lambda column: column.str.lower(),
    )

    metrics.to_csv(
        OUTPUT_DIR / "spectral_metrics.csv",
        index=False,
    )

    if failed_files:
        pd.DataFrame(
            failed_files,
            columns=["file_name", "error"],
        ).to_csv(
            OUTPUT_DIR / "failed_audio_files.csv",
            index=False,
        )

    # Match the spectrum list order to the metrics table.
    spectrum_lookup = {
        item["file_name"]: item
        for item in all_spectra
    }

    ordered_spectra = [
        spectrum_lookup[file_name]
        for file_name in metrics["file_name"]
    ]

    make_heatmap(
        ordered_spectra,
        metrics.reset_index(drop=True),
    )

    make_rolloff_plot(
        metrics
    )

    print()
    print("Analysis complete.")
    print(f"Outputs saved to:\n{OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()