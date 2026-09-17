import streamlit as st
from module0 import global_aggregate
from app_helpers import AnalysisStatus, TranscriptIDMode, GenomicCoordinatesMode

from bokeh.plotting import figure
from bokeh.models import (
    ColumnDataSource, HoverTool, PanTool, Range1d, Span, Label, Legend
)
import pandas as pd
from bokeh.models import LinearAxis
from bokeh.layouts import column
from bokeh.models import BoxAnnotation
from streamlit_bokeh import streamlit_bokeh

# Frame colors used everywhere — single source of truth
FRAME_COLORS = {1: "#48cf4b", 2: "#68ace3", 3: "#e66343"}
CODON_COLORS = {"ATG": "white", "stop codon": "black", "CTG": "#e3e3e3"}

# Introducing smoothing function:
def plot_smoother(df: pd.DataFrame, smoothing:bool = True):
    """
    Calculates the frame over the global aggregate, then smoothes by it.
    Also adds column "X" used by plotting function "subplot_single_transcript_plot"

    :param df: pandas.DataFrame
        DataFrame containing a 'Global Aggregate' column
    :param smoothing: logical. Default is true. If false, data is not smoothed.
    :return df: pandas.DataFrame
        DataFrame with the smoothed Global Aggregate
    :return k: int. The k-mer value (k) applied to smooth peaks, if smoothing == True.
    """
    # Copies the position to column "X". Used to calculate Frame and for visualisation later
    df["X"] = df.index

    df["X"] = df["X"].astype(int)
    # Assign frame to each position
    df["Frame"] = (df["X"] % 3) + 1

    if smoothing == True:
        # Code to smooth cov by frame. Adjust based on transcript length
        transcript_len = len(df.index)
        if transcript_len < 3500:
            k = 3
        elif transcript_len < 6000:
            k = 5
        elif transcript_len < 10000:
            k = 9
        else:
            k = 13

        k2 = k * 3
        for i in range(0, len(df), k2):
            # Select the slice of 9 rows using iloc
            chunk = df.iloc[i:i + k2]
            # out of those, groupby Frame, calculate the mean
            chunk_averages = chunk.groupby('Frame')['Global Aggregate'].transform("mean")
            # store the new values as a new column
            df.loc[chunk.index, "Global Aggregate smoothed"] = chunk_averages

        return df, k

    elif smoothing == False:
        # No smoothing, indicated as k = 1
        k=1

        return df, k


def _map_type(codon):
    if codon == "ATG":
        return "ATG"
    elif codon in ["CTG", "GTG", "TTG"] : #, "GTG", "TTG"]:
        return "CTG"
    elif codon in ["TGA", "TAA", "TAG"]:
        return "stop codon"
    else:
        return None

# Helper - finds codon in a given frame
def _find_codons(seq, frame, codons):
    rows = []
    for i in range(frame - 1, len(seq) - 2, 3):
        codon = str(seq[i:i + 3])
        # if codon in codons:  # This way everything is mapped, regardless of given list
        rows.append({
            "codon": codon,
            "position": i,
            "Frame": int(frame),
            "type": _map_type(codon)
        })
    locations = pd.DataFrame(rows)

    return locations


def compute_codons_by_frame(sequence):
    codons = ["ATG", "CTG", "GTG", "TTG", "TGA", "TAA", "TAG"]
    frame_1 = _find_codons(sequence, 1, codons)
    frame_2 = _find_codons(sequence, 2, codons)
    frame_3 = _find_codons(sequence, 3, codons)

    temp = pd.concat([frame_1, frame_2, frame_3], ignore_index=True)
    all_reading_frames = temp.sort_values(by="position", ascending = True).reset_index(drop = True)

    return all_reading_frames


def _add_vline_and_label(p, x, max_y, text):
    """Add a dashed vertical line and a label to mark start/stop codons."""
    vline = Span(
        location=x,
        dimension="height",
        line_color="black",
        line_width=2,
        line_dash="dashed",
    )
    label = Label(
        x=x,
        y=max_y,
        text=text,
        text_font_size="11pt",
        text_color="black",
        text_align="center",
        y_offset=-20,
        y_units="data",
        background_fill_color="white",
        background_fill_alpha=1.0,
        padding=4,
    )
    p.add_layout(vline)
    p.add_layout(label)


def subplot_single_transcript_plot2(
    df: pd.DataFrame,
    x_range: Range1d,
    y_axis: str,
    utr5_length: int,
    cds_length: int,
    txid: str = "",
    k: int = 1,
    mappability_df: pd.DataFrame = None,
) -> figure:
    """
    Frame-coloured coverage plot. Uses vbar glyphs (one per frame) so hover
    snaps cleanly to nucleotide positions.
    """
    df = df.copy()
    df["Frame"] = df["Frame"].astype(int)

    y_upper = df[y_axis].max()
    y_upper = max(y_upper, 10) # If a transcript has 0 counts, then we set the max y on the plot to 10 (for visualisation)
    y_lower = min(-1 * df[y_axis].quantile(0.10), -0.05 * y_upper)
    y_lower = min(y_lower, -0.05)  # Ading safety line to not encour into "-0.0" bug. We keep some margin, but avoid -0.0

    title = (
        f"Global Aggregate for {txid}, smoothed by k = {k}"
        if k > 1 else f"Global Aggregate for {txid}"
    )

    p = figure(
        x_range=x_range,
        y_range=Range1d(start=y_lower, end=y_upper * 1.2),
        sizing_mode="stretch_width",
        height=700,
        background_fill_color="white",
        border_fill_color="white",
        x_axis_label="Position on the transcript (nt)",
        y_axis_label="Aggregate Counts",
        title=title,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        active_scroll="wheel_zoom",
    )
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.outline_line_color = None

    # We plot tas area the various frames, but use an invisible bar for the hover functionality
    legend_items = []
    hover_renderers = []
    for frame in (1, 2, 3):
        sub = df[df["Frame"] == frame].sort_values("X")
        if sub.empty:
            continue
        source = ColumnDataSource(sub)

        if k == 1:
            # Histogram mode: one bar per nucleotide (easier to observe triplet periodicity)
            renderer = p.vbar(
                x="X", top=y_axis, bottom=0, width=1,
                source=source,
                fill_color=FRAME_COLORS[frame],
                fill_alpha=0.6,
                line_color=None,
            )
            # The renderer itself is hoverable, no need for a separate invisible bar
            hover_renderers.append(renderer)
        else:

            # Smoothed mode: area plot
            renderer = p.varea(
                x="X", y1=0, y2=y_axis,
                source=source,
                fill_color=FRAME_COLORS[frame],
                fill_alpha=0.5,
            )

            # Invisible hover bars
            hover_bar = p.vbar(
                x="X", top=y_axis, bottom=0, width=1,
                source=source,
                fill_color=FRAME_COLORS[frame],
                fill_alpha=0,
                hover_fill_color=FRAME_COLORS[frame],
                hover_fill_alpha=0.4,
                line_color=None,
            )
            hover_renderers.append(hover_bar)

        legend_items.append((f"Frame {frame}", [renderer]))

    # Optional mappability overlay on secondary y-axis
    if mappability_df is not None:

        p.extra_y_ranges = {"mappability": Range1d(start=0, end=1)}

        mapp_line = p.line(
            x="X",
            y="Mappability",
            source=ColumnDataSource(mappability_df),
            line_color="gray",
            line_dash="dotted",
            line_width=1.5,
            y_range_name="mappability",
        )

        p.add_layout(
            LinearAxis(
                y_range_name="mappability",
                axis_label="Mappability",
                axis_label_text_color="gray",
                major_label_text_color="gray",
                axis_line_color="gray",
            ),
            "right",
        )

        legend_items.append(("Mappability", [mapp_line]))

    # Start/Stop labels
    start_pos = utr5_length
    stop_pos = utr5_length + cds_length

    _add_vline_and_label(p, start_pos, y_upper * 1.1, "Start")
    _add_vline_and_label(p, stop_pos, y_upper * 1.1, "Stop")

    # Legend OUTSIDE the plot area, above the title
    legend = Legend(items=legend_items, orientation="horizontal", click_policy="hide")
    p.add_layout(legend, "above")

    # Hover — vline mode now works because we have vbar glyphs with x positions
    hover = HoverTool(
        tooltips=[
            ("Position", "@X"),
            ("Frame", "@Frame"),
            ("Value", f"@{{{y_axis}}}{{0.00}}"),
        ],
        mode="vline",
        renderers = hover_renderers,
    )
    p.add_tools(hover)

    # Typography - Size
    p.title.text_font_size = "14pt"
    p.xaxis.axis_label_text_font_size = "13pt"
    p.xaxis.major_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "13pt"
    p.yaxis.major_label_text_font_size = "14pt"
    # Typography - Colour
    p.title.text_color = "black"
    p.xaxis.axis_label_text_color = "black"
    p.yaxis.axis_label_text_color = "black"
    p.xaxis.major_label_text_color = "black"
    p.yaxis.major_label_text_color = "black"
    p.legend.label_text_color = "black"
    p.legend.background_fill_color = "white"

    # Limit pan and zoom
    p.y_range.bounds = (y_lower, y_upper * 1.5)
    p.x_range.bounds = (x_range.start, x_range.end)
    p.extra_y_ranges = {"mappability": Range1d(start=0, end=1, bounds=(0, 1))}

    return p


def ORF_plotter2(
    df: pd.DataFrame,
    frame: int,
    x_range: Range1d,
) -> figure:
    """
    Single-frame ORF track. Background coloured by frame, with ATG (white),
    stop codons (black), and CTG (light gray) highlighted at their positions.

    Args:
        df: DataFrame with columns 'position', 'Frame', and 'type'
        frame: Reading frame to display (1, 2, or 3)
        x_range: Shared Range1d for x-axis (links to main plot)

    Returns:
        Bokeh figure object for use as a stacked track.
    """
    plot_df = df.copy()

    # Mask: rows not belonging to this frame lose their type (will get bg colour)
    mask = plot_df["Frame"] != frame
    plot_df.loc[mask, "type"] = None
    plot_df["type"] = plot_df["type"].astype(str)

    bg_colour = FRAME_COLORS[frame]

    codon_color_map = {
        "ATG": "white",
        "stop codon": "black",
        "CTG": "#e3e3e3",
        "None": bg_colour,
        "nan": bg_colour,
    }
    plot_df["color"] = plot_df["type"].map(lambda x: codon_color_map.get(x, bg_colour))

    source = ColumnDataSource(plot_df)

    p = figure(
        x_range=x_range,
        y_range=Range1d(start=0, end=1),
        height=75,
        sizing_mode="stretch_width",
        title=f"Frame {frame}",
        background_fill_color="white",
        border_fill_color="white",
        tools="save",
        toolbar_location="right",
    )

    bars = p.vbar(
        x="position",
        top=1,
        bottom=0,
        width=1,
        color="color",
        source=source,
    )

    hover = HoverTool(
        tooltips=[
            ("Position", "@position"),
            ("Codon Type", "@type"),
        ],
        renderers=[bars],
    )
    p.add_tools(hover)

    # Hide axes and gridlines — decorative track
    p.yaxis.visible = False
    p.xaxis.visible = False
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.outline_line_color = None

    # Title styling
    p.title.text_color = "black"
    p.title.text_font_size = "10pt"

    # Bound pan to data range
    p.x_range.bounds = (x_range.start, x_range.end)
    p.y_range.bounds = (0, 1)

    return p

def plot_ga_with_frames2(
    current_run: AnalysisStatus,
    all_studies_df: pd.DataFrame,
    y_column: str,
    smoothing: bool = True,
    normalize: bool = False,
) -> column:
    """
    Frame-aware global aggregate plot with ORF tracks beneath, all sharing
    a linked x-axis. Mappability (if available) is overlaid on the main plot.
    """
    # Prepare data
    ga_df = global_aggregate(all_studies_df, normalize = normalize)
    processed_df, k = plot_smoother(ga_df, smoothing)

    # Build shared numeric x-range
    x_min = int(processed_df["X"].min())
    x_max = int(processed_df["X"].max())
    margin = max(1, int((x_max - x_min) * 0.01))
    shared_x_range = Range1d(start=x_min - margin, end=x_max + margin)

    # Main transcript plot (with optional mappability overlay)
    transcript_plot = subplot_single_transcript_plot2(
        df=processed_df,
        x_range=shared_x_range,
        y_axis=y_column,
        utr5_length=current_run.utr5_length,
        cds_length=current_run.cds_length,
        txid=current_run.txid,
        k=k,
        mappability_df=current_run.mappability_df,
    )

    # ORF tracks for each frame
    orf_info = compute_codons_by_frame(current_run.transcript_sequence)
    orf_plots = [ORF_plotter2(orf_info, frame, shared_x_range) for frame in (1, 2, 3)]

    # Stack everything
    combined_plot = column(
        transcript_plot,
        *orf_plots,
        sizing_mode="stretch_width",
        height=925,
    )

    return combined_plot

# Makes transcript plot for genco mode
def colorless_transcript_plot(
    df: pd.DataFrame,
    y_axis: str,
    k: int = 1,
    mappability_df: pd.DataFrame = None,
) -> figure:
    """
    Simple coverage plot. Uses vbar glyphs (one per frame) so hover
    snaps cleanly to nucleotide positions.
    """
    df = df.copy()

    y_upper = max(1, df[y_axis].max())
    y_lower = min(-1 * df[y_axis].quantile(0.10), -0.05 * y_upper)

    title = (
        f"Global Aggregate, smoothed by k = {k}"
        if k > 1 else f"Global Aggregate"
    )

    x_range = Range1d(start=df["Position"].min(), end=df["Position"].max())

    p = figure(
        x_range=x_range,
        y_range=Range1d(start=y_lower, end=y_upper * 1.2),
        sizing_mode="stretch_width",
        height=700,
        background_fill_color="white",
        border_fill_color="white",
        x_axis_label="Position on the region (nt)",
        y_axis_label="Aggregate Counts",
        title=title,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        active_scroll="wheel_zoom",
    )
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.outline_line_color = None

    # We plot tas area the various frames, but use an invisible bar for the hover functionality
    legend_items = []
    hover_renderers = []

    source = ColumnDataSource(df)

    area = p.varea(
            x="Position", y1=0, y2=y_axis,
            source=source,
            fill_color="gray",
            fill_alpha=0.5,
        )

        # Invisible hover bars
    hover_bar = p.vbar(
            x="Position", top=y_axis, bottom=0, width=1,
            source=source,
            fill_color="gray",
            fill_alpha=0,
            hover_fill_color="gray",
            hover_fill_alpha=0.4,
            line_color=None,
        )
    hover_renderers.append(hover_bar)

    legend_items.append((f"Coverage", [area]))

    # Optional mappability overlay on secondary y-axis
    if mappability_df is not None:

        p.extra_y_ranges = {"mappability": Range1d(start=0, end=1)}

        mapp_line = p.line(
            x="Position",
            y="Mappability",
            source=ColumnDataSource(mappability_df),
            line_color="black",
            line_dash="dotted",
            line_width=1.5,
            y_range_name="mappability",
        )

        p.add_layout(
            LinearAxis(
                y_range_name="mappability",
                axis_label="Mappability",
                axis_label_text_color="gray",
                major_label_text_color="gray",
                axis_line_color="gray",
            ),
            "right",
        )

        legend_items.append(("Mappability", [mapp_line]))


    # Legend OUTSIDE the plot area, above the title
    legend = Legend(items=legend_items, orientation="horizontal", click_policy="hide")
    p.add_layout(legend, "above")

    # Hover — vline mode now works because we have vbar glyphs with x positions
    hover = HoverTool(
        tooltips=[
            ("Position", "@Position"),
            ("Value", f"@{{{y_axis}}}{{0.00}}"),
        ],
        mode="vline",
        renderers = hover_renderers,
    )
    p.add_tools(hover)

    # Typography - Size
    p.title.text_font_size = "14pt"
    p.xaxis.axis_label_text_font_size = "13pt"
    p.xaxis.major_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "13pt"
    p.yaxis.major_label_text_font_size = "14pt"
    # Typography - Colour
    p.title.text_color = "black"
    p.xaxis.axis_label_text_color = "black"
    p.yaxis.axis_label_text_color = "black"
    p.xaxis.major_label_text_color = "black"
    p.yaxis.major_label_text_color = "black"
    p.legend.label_text_color = "black"
    p.legend.background_fill_color = "white"

    # Limit pan and zoom
    p.y_range.bounds = (y_lower, y_upper * 1.5)
    p.x_range.bounds = (x_range.start, x_range.end)
    p.extra_y_ranges = {"mappability": Range1d(start=0, end=1, bounds=(0, 1))}

    return p



#### NEW HELPERS FOR COLOURLESS STP (Used in page 5 - tab5) ###

def plot_stp_with_frames(
    current_run: AnalysisStatus,
    all_studies_df: pd.DataFrame,
    run_option: list = None,
    normalize: bool = False,
) -> column:
    """
    Per-study global aggregate plot with ORF tracks beneath, sharing a
    linked x-axis. Raw counts, no smoothing.
    """

    df = all_studies_df.copy()
    if "Position" in df.columns:
        df = df.drop(columns = "Position")
    df["X"] = df.index

    x_min = int(df["X"].min())
    x_max = int(df["X"].max())
    margin = max(1, int((x_max - x_min) * 0.01))
    shared_x_range = Range1d(start=x_min - margin, end=x_max + margin)

    coordinates = {"x_slice": current_run.x_slice, "y_slice": current_run.y_slice, "n_slice": current_run.n_slice}

    transcript_plot = subplot_single_transcript_plot(
        df=df,
        x_range=shared_x_range,
        utr5_length=current_run.utr5_length,
        cds_length=current_run.cds_length,
        txid=current_run.txid,
        run_option=run_option,
        mappability_df=current_run.mappability_df,
        coordinates = coordinates
    )

    orf_info = compute_codons_by_frame(current_run.transcript_sequence)
    orf_plots = [ORF_plotter2(orf_info, frame, shared_x_range) for frame in (1, 2, 3)]

    combined_plot = column(
        transcript_plot, *orf_plots,
        sizing_mode="stretch_width", height=925,
    )
    return combined_plot

from bokeh.palettes import Category10, Category20

def _add_region_box(p, start, end, label=None, color="#d4d4d4", alpha=0.35):
    """Shade a background region spanning the full plot height (like the
    sandpaper regions in the old Plotly plot)."""
    box = BoxAnnotation(
        left=start,
        right=end,
        fill_color=color,
        fill_alpha=alpha,
        line_width=0,
        level="underlay",  # sits behind the glyphs
    )
    p.add_layout(box)
    if label is not None:
        lbl = Label(
            x=(start + end) / 2,
            y=0,
            y_units="screen",
            y_offset=-4,
            text=label,
            text_font_size="10pt",
            text_color="black",
            text_align="center",
            background_fill_color=color,
            background_fill_alpha=0.3,
            padding=3,
        )
        p.add_layout(lbl)

def subplot_single_transcript_plot(
    df: pd.DataFrame,
    x_range: Range1d,
    utr5_length: int,
    cds_length: int,
    txid: str = "",
    run_option: list = None,
    mappability_df: pd.DataFrame = None,
    coordinates: dict = None,
) -> figure:
    """
    Per-study coverage plot. Each SRR/study column is drawn as its own
    coloured bar trace (raw counts, no smoothing, no frame colouring).
    Expects wide format: an 'X' column plus one count column per study.
    """
    df = df.copy()

    # Identify study columns (everything except X / Mappability)
    reserved = {"X", "Mappability"}
    study_cols = [c for c in df.columns if c not in reserved]
    if run_option is not None:
        study_cols = [c for c in study_cols if c in run_option]

    # Colour palette sized to the number of studies
    n = len(study_cols)
    if n <= 10:
        palette = Category10[max(3, n)]
    else:
        palette = (Category20[20] * (n // 20 + 1))[:n]

    # Y bounds across all plotted studies
    if study_cols:
        y_upper = max(df[study_cols].max().max(), 10)
    else:
        y_upper = 10
    y_lower = min(-0.05 * y_upper, -0.05)  # avoid the -0.0 bounds bug

    p = figure(
        x_range=x_range,
        y_range=Range1d(start=y_lower, end=y_upper * 1.2),
        sizing_mode="stretch_width",
        height=700,
        background_fill_color="white",
        border_fill_color="white",
        x_axis_label="Position on the transcript (nt)",
        y_axis_label="Raw Counts",
        title=f"Raw Counts for {txid}",
        tools="pan,wheel_zoom,box_zoom,reset,save",
        active_scroll="wheel_zoom",
    )
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.outline_line_color = None

    # One bar trace per study
    legend_items = []
    hover_renderers = []
    for i, col in enumerate(study_cols):
        sub = df[["X", col]].sort_values("X")
        source = ColumnDataSource(sub)
        renderer = p.vbar(
            x="X", top=col, bottom=0, width=1,
            source=source,
            fill_color=palette[i % len(palette)],
            fill_alpha=0.6,
            hover_fill_alpha=0.9,
            line_color=None,
        )
        # per-renderer hover so the tooltip reads the right column
        p.add_tools(HoverTool(
            tooltips=[("Position", "@X"), ("Study", col), ("Counts", f"@{{{col}}}{{0}}")],
            mode="vline",
            renderers=[renderer],
        ))
        legend_items.append((col, [renderer]))
        hover_renderers.append(renderer)

    # Optional mappability overlay on secondary y-axis
    if mappability_df is not None:
        p.extra_y_ranges = {"mappability": Range1d(start=0, end=1, bounds=(0, 1))}
        mapp_line = p.line(
            x="X", y="Mappability",
            source=ColumnDataSource(mappability_df),
            line_color="gray", line_dash="dotted", line_width=1.5,
            y_range_name="mappability",
        )
        p.add_layout(
            LinearAxis(
                y_range_name="mappability",
                axis_label="Mappability",
                axis_label_text_color="gray",
                major_label_text_color="gray",
                axis_line_color="gray",
            ),
            "right",
        )
        legend_items.append(("Mappability", [mapp_line]))

    # Start/Stop labels
    start_pos = utr5_length
    stop_pos = utr5_length + cds_length
    _add_vline_and_label(p, start_pos, y_upper * 1.1, "Start")
    _add_vline_and_label(p, stop_pos, y_upper * 1.1, "Stop")

    if coordinates is not None:
        region_labels = {
            "x_slice": "Baseline",
            "y_slice": "Target",
            "n_slice": "Noise Region",
        }
        for key, (start, end) in coordinates.items():
            _add_region_box(p, start, end, label=region_labels.get(key))


    legend = Legend(items=legend_items, orientation="horizontal", click_policy="hide")
    p.add_layout(legend, "above")

    # Typography
    p.title.text_font_size = "14pt"
    p.xaxis.axis_label_text_font_size = "13pt"
    p.xaxis.major_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "13pt"
    p.yaxis.major_label_text_font_size = "14pt"
    p.title.text_color = "black"
    p.xaxis.axis_label_text_color = "black"
    p.yaxis.axis_label_text_color = "black"
    p.xaxis.major_label_text_color = "black"
    p.yaxis.major_label_text_color = "black"
    p.legend.label_text_color = "black"
    p.legend.background_fill_color = "white"

    # Limit pan and zoom
    p.y_range.bounds = (y_lower, y_upper * 1.5)
    p.x_range.bounds = (x_range.start, x_range.end)

    return p


def genco_single_transcript_plot(
    df: pd.DataFrame,
    x_range: Range1d,
    run_option: list = None,
    mappability_df: pd.DataFrame = None,
    genomic_merged_ranges: dict = None,    # needs info on baseline/target/noise
) -> figure:
    """
    Per-study coverage plot. Each SRR/study column is drawn as its own
    coloured bar trace (raw counts, no smoothing, no frame colouring).
    Expects wide format: an 'X' column plus one count column per study.
    """
    df = df.copy()

    # Identify study columns (everything except X / Mappability)
    reserved = {"X", "Mappability"}
    study_cols = [c for c in df.columns if c not in reserved]
    if run_option is not None:
        study_cols = [c for c in study_cols if c in run_option]

    # Colour palette sized to the number of studies
    n = len(study_cols)
    if n <= 10:
        palette = Category10[max(3, n)]
    else:
        palette = (Category20[20] * (n // 20 + 1))[:n]

    # Y bounds across all plotted studies
    if study_cols:
        y_upper = max(df[study_cols].max().max(), 10)
    else:
        y_upper = 10
    y_lower = min(-0.05 * y_upper, -0.05)  # avoid the -0.0 bounds bug

    p = figure(
        x_range=x_range,
        y_range=Range1d(start=y_lower, end=y_upper * 1.2),
        sizing_mode="stretch_width",
        height=700,
        background_fill_color="white",
        border_fill_color="white",
        x_axis_label="Position on the transcript (nt)",
        y_axis_label="Raw Counts",
        title=f"Raw Counts",
        tools="pan,wheel_zoom,box_zoom,reset,save",
        active_scroll="wheel_zoom",
    )
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.outline_line_color = None

    # One bar trace per study
    legend_items = []
    hover_renderers = []
    for i, col in enumerate(study_cols):
        sub = df[["X", col]].sort_values("X")
        source = ColumnDataSource(sub)
        renderer = p.vbar(
            x="X", top=col, bottom=0, width=1,
            source=source,
            fill_color=palette[i % len(palette)],
            fill_alpha=0.6,
            hover_fill_alpha=0.9,
            line_color=None,
        )
        # per-renderer hover so the tooltip reads the right column
        p.add_tools(HoverTool(
            tooltips=[("Position", "@X"), ("Study", col), ("Counts", f"@{{{col}}}{{0}}")],
            mode="vline",
            renderers=[renderer],
        ))
        legend_items.append((col, [renderer]))
        hover_renderers.append(renderer)

    # Optional mappability overlay on secondary y-axis
    if mappability_df is not None:
        p.extra_y_ranges = {"mappability": Range1d(start=0, end=1, bounds=(0, 1))}
        mapp_line = p.line(
            x="X", y="Mappability",
            source=ColumnDataSource(mappability_df),
            line_color="gray", line_dash="dotted", line_width=1.5,
            y_range_name="mappability",
        )
        p.add_layout(
            LinearAxis(
                y_range_name="mappability",
                axis_label="Mappability",
                axis_label_text_color="gray",
                major_label_text_color="gray",
                axis_line_color="gray",
            ),
            "right",
        )
        legend_items.append(("Mappability", [mapp_line]))

    # Code to add boundaries regions (baseline/target/noise)
    delim_1 = genomic_merged_ranges["y_slice"][0]
    delim_2 = genomic_merged_ranges["y_slice"][1]
    vline = Span(
        location=delim_1,
        dimension="height",
        line_color="black",
        line_width=2,
        line_dash="dashed",
    )
    vline2 = Span(
        location=delim_2,
        dimension="height",
        line_color="black",
        line_width=2,
        line_dash="dashed",
    )
    p.add_layout(vline)
    p.add_layout(vline2)
    # Code to add textbox over regions
    _add_region_box(p, genomic_merged_ranges["x_slice"][0],
                    genomic_merged_ranges["x_slice"][1], label="Baseline", color="#d4d4d4", alpha=0)

    _add_region_box(p, genomic_merged_ranges["y_slice"][0],
                    genomic_merged_ranges["y_slice"][1], label="Target", color="#d4d4d4", alpha=0.35)
    _add_region_box(p, genomic_merged_ranges["n_slice"][0],
                    genomic_merged_ranges["n_slice"][1], label="Noise", color="#d4d4d4", alpha=0)

    legend = Legend(items=legend_items, orientation="horizontal", click_policy="hide")
    p.add_layout(legend, "above")

    # Typography
    p.title.text_font_size = "14pt"
    p.xaxis.axis_label_text_font_size = "13pt"
    p.xaxis.major_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "13pt"
    p.yaxis.major_label_text_font_size = "14pt"
    p.title.text_color = "black"
    p.xaxis.axis_label_text_color = "black"
    p.yaxis.axis_label_text_color = "black"
    p.xaxis.major_label_text_color = "black"
    p.yaxis.major_label_text_color = "black"
    p.legend.label_text_color = "black"
    p.legend.background_fill_color = "white"

    # Limit pan and zoom
    p.y_range.bounds = (y_lower, y_upper * 1.5)
    p.x_range.bounds = (x_range.start, x_range.end)

    return p