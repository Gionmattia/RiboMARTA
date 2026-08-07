# mappability_helpers.py
from transcriptID_mode_helpers import *
from bokeh.transform import linear_cmap
from bokeh.models import ColumnDataSource, HoverTool, PanTool
from bokeh.plotting import figure

# Extracts mappability_for_txid as df
def txid_map_extraction(txid, mappability_path, gr):
    """
    Function loads coverage information into a dataframe (all_studies_df) which is stored into current_run.
    Also extract features lengths (5', CDs, 3') and stores them in current_run object.

    :param current_run: current run - object where txid mode variables and results are stored. See app_helpers.py
    :param mappability_path: absolute path to the mappability bigwig file.
    :param gr: annotation object, as loaded by pyranges.
    :return: all_ranges, strandedness
    """
    # 1 Retrieve transcript transcript_ranges
    transcript_ranges = extract_tx_pyranges(gr, txid)
    # 2  Retrieve transcript features (CDS, UTRs) transcript_ranges and strandedness
    all_ranges, strandedness = regions_df_formatter(transcript_ranges)
    # FIX TO work with UCSC format
    all_ranges["Chromosome"] = "chr" + all_ranges["Chromosome"].astype(str)

    # 4 Extract coverage from given bigwigs
    path_as_list = [] # Workaround, original function works on lists.
    path_as_list.append(mappability_path)
    mappability_df = extract_from_multiple_bigwigs(path_as_list, all_ranges, strandedness, n_jobs=1)

    return mappability_df


# Makes plot (bokeh) under main plot
def bokeh_mappability_track(mappability_df, x_range=None):
    # Masking as in your code
    plot_df = mappability_df.copy()
    plot_df['Position'] = plot_df['Position'].astype(str)

    if x_range is None:
        x_range = plot_df['Position'].tolist()

    source = ColumnDataSource(plot_df)
    palette = [f"#{i:02x}{i:02x}{i:02x}" for i in range(255, -1, -1)]# 0 is lighter, 1 is darker

    mapper = linear_cmap(
        field_name='Mappability',  # your 0-1 column
        palette=palette,
        low=0,
        high=1
    )

    p = figure(
        x_range=x_range,
        y_range=(0, 0.1),
        height= 35,
        width=700,
        title=None,
        toolbar_location= None,
        tools="wheel_zoom, reset",
        min_border_top=0,
        min_border_bottom=0,
    )

    bars = p.vbar(
        x='Position',
        top=0.1,
        width=0.9,
        color=mapper,
        source=source,
    )

    hover = HoverTool(
        tooltips=[
            ("Position", "@Position"),
            ("Mappability", "@Mappability{0.3f}"),
        ],
        renderers=[bars]
    )
    pan = PanTool(dimensions="width")
    p.add_tools(hover, pan)
    p.toolbar.active_drag = pan

    p.yaxis.visible = False
    p.xaxis.visible = False
    p.ygrid.grid_line_color = None

    return p