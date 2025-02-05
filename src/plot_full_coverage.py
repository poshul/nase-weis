# Authors: Hendrik Weisser, Samuel Wein
# License: BSD 3-clause

import csv
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import importlib
from io import StringIO
import os


#' Read an mzTab file containing NucleicAcidSearchEngine results
#'
#' @param path Path to the mzTab file
#'
#' @return List of data frames for the different mzTab parts
def read_mzTab(path):
    def process_buffer(buffer, col_names, result):
        data = pd.read_csv(StringIO('\n'.join(buffer)), sep="\t", header=None, names=col_names, na_values="null")
        result[data.iloc[0, 0]] = data.iloc[:, 1:]
        return result

    with open(path, 'r') as file:
        content = file.readlines()

    prefixes = [line[:3] for line in content]
    result = {}
    buffer = []
    col_names = ["MTH", "name", "value"]

    for i, line in enumerate(content):
        if not line.strip():
            continue

        if prefixes[i][2] == "H" and (i == 0 or content[i - 1].strip() == ""):
            if buffer:
                result = process_buffer(buffer, col_names, result)
            col_names = content[i].strip().split("\t")
            buffer = []
        else:
            buffer.append(line.strip())

    process_buffer(buffer, col_names, result)
    return result

#' Split a (modified) oligonucleotide sequence given as a string into a list of
#' (modified) nucleotides
def split_sequence(seq):
    split = list(seq)
    bracket_open = [i for i, char in enumerate(split) if char == "["]
    bracket_close = [i for i, char in enumerate(split) if char == "]"]

    assert len(bracket_open) == len(bracket_close), "Mismatch between open and close brackets"

    mods = []
    if bracket_open:
        for i in range(len(bracket_open)):
            mods.append("".join(split[bracket_open[i] + 1:bracket_close[i]]))
        for i in range(len(bracket_open)):
            del split[bracket_open[i] + 1:bracket_close[i]]

    remove = [] #FIXME
    if split[0] == "p":
        split[1] = "p" + split[1]
        split = split[1:]

    if split[-1] == "p":
        split[-2] = split[-2] + "p"
        split = split[:-1]

    return split

#' Resolve ambiguities in sequence assignments due to different localisations
#' of modifications
# TODO: This is untested use with caution
def resolve_ambiguities(mztab, groups={"m5C?": ["C", "m5C"]}):
    """
    Resolve ambiguities in sequence assignments due to different localizations of modifications.
    """
    mztab["OSM"]["sequence"] = mztab["OSM"]["sequence"].astype(str)
    mztab["OLI"]["sequence"] = mztab["OLI"]["sequence"].astype(str)
    parts = {key: df for key, df in mztab["OSM"].groupby("spectra_ref")}
    nrows = {key: len(df) for key, df in parts.items()}
    ambig = {key: df for key, df in parts.items() if nrows[key] > 1}
    unexpected = set()
    
    resolved = []
    for spectra_ref, part in ambig.items():
        split_seqs = [list(seq) for seq in part["sequence"]]
        lengths = [len(seq) for seq in split_seqs]
        assert all(length == lengths[0] for length in lengths), "Length mismatch in sequences."
        
        mat = np.array(split_seqs)
        consensus = []
        
        for column in mat.T:
            unique_nucs = sorted(set(column))
            if len(unique_nucs) > 1:
                matches = [group for group, members in groups.items() if set(unique_nucs).issubset(members)]
                if matches:
                    consensus.append(matches[0])
                else:
                    signature = ", ".join(unique_nucs)
                    if signature not in unexpected:
                        print(f"Warning: unexpected ambiguity: {signature}")
                        unexpected.add(signature)
                    consensus.append("|".join(unique_nucs))
            else:
                consensus.append(unique_nucs[0])
        
        con_seq = "".join(
            f"[{nuc}]" if len(nuc) > 2 or (len(nuc) == 2 and not (nuc.startswith("p") or nuc.endswith("p"))) else nuc
            for nuc in consensus
        )
        
        if con_seq not in mztab["OLI"]["sequence"].values:
            rows = mztab["OLI"][mztab["OLI"]["sequence"].isin(part["sequence"])].copy()
            rows["sequence"] = con_seq
            mztab["OLI"] = pd.concat([mztab["OLI"], rows.drop_duplicates()], ignore_index=True)
        
        resolved.append(part.iloc[[0]].assign(sequence=con_seq))
    
    unambig = [df for key, df in parts.items() if nrows[key] == 1]
    mztab["OSM"] = pd.concat(unambig + resolved, ignore_index=True)
    mztab["OLI"] = mztab["OLI"][mztab["OLI"]["sequence"].isin(mztab["OSM"]["sequence"])]
    mztab["OSM"]["sequence"] = mztab["OSM"]["sequence"].astype("category")
    mztab["OLI"]["sequence"] = mztab["OLI"]["sequence"].astype("category")
    
    return mztab

#' Read in RNA modification definitions (MODOMICS format)
#'
#' @param directory Directory containing modification data file(s)
#' @param include_custom Include custom mod. definitions?
#'
#' @return Data frame containing modification definitions
def read_rna_modifications(directory=None, include_custom=True):
    """
    Read in RNA modification definitions (MODOMICS format).

    Parameters:
        directory (str): Directory containing modification data file(s).
        include_custom (bool): Whether to include custom modification definitions.

    Returns:
        pd.DataFrame: Data frame containing modification definitions.
    """
    print(directory)
    if directory is None:
        pyopenms_loc = importlib.util.find_spec("pyopenms").submodule_search_locations[0]
        if pyopenms_loc is not None:
            directory = os.path.join(pyopenms_loc, "share/OpenMS/CHEMISTRY/")
        else:
            raise ModuleNotFoundError("pyopenms not found")

    # Define column types
    col_types = {i: str for i in range(9)}
    col_types[3] = "category"
    col_types[6] = "category"
    col_types[7] = float
    col_types[8] = float

    # Read Modomics.tsv
    mod_file = os.path.join(directory, "Modomics.tsv")
    mod_info = pd.read_csv(mod_file, sep="\t", quoting=csv.QUOTE_NONE, encoding="utf-8", skiprows=4, dtype=col_types, na_values="None")

    if include_custom:
        # Add extra column type for alternatives
        col_types[len(col_types)] = str  
        custom_file = os.path.join(directory, "Custom_RNA_modifications.tsv")
        custom_mods = pd.read_csv(custom_file, sep="\t", quoting=csv.QUOTE_NONE, encoding="utf-8", skiprows=4, dtype=col_types)

        # Add missing "alternatives" column if not present
        if "alternatives" not in mod_info.columns:
            mod_info["alternatives"] = ""

        # Append custom modifications
        mod_info = pd.concat([mod_info, custom_mods], ignore_index=True)

    # Use just "Q" instead of "QtRNA" for queuosine
    mod_info["short_name"] = mod_info["short_name"].str.replace("QtRNA$", "Q", regex=True)

    return mod_info


## Coverage plotting functions:

#' Generate a discrete color scale for a range of values
#'
#' @param nmax Maximum value
#' @param nmin Minimum value
#' @param n_colors Number of colors to use
#'
#' @return Vector of break points named with color codes
def get_color_scale(nmax, nmin=1, n_colors=8):
    assert nmin <= nmax
    nmax = nmax - nmin + 1

    if nmax <= n_colors:
        scale = np.arange(1, nmax + 1)
    else:
        exp = n_colors - 1
        if nmax < np.ceil(1.5 ** exp):
            scale = np.linspace(1, nmax, n_colors)
        else:
            base = nmax ** (1 / exp)
            scale = np.ceil(np.round(base ** np.arange(0, exp + 1), 1))

    scale = scale + nmin - 1
    colors = plt.cm.magma(np.linspace(0, 1, len(scale)))
    scale_colors = {str(i): mcolors.to_hex(colors[i]) for i in range(len(scale))}
    return scale, scale_colors


#' Generate a table of "stacked bars" representing oligonucleotides
#'
#' @param oligo_data Oligonucleotide data (from mzTab)
#' @param osm_data Spectrum-match data (from mzTab)
#' @param nuc_length Length of the full RNA sequence
#' @param color_scale Color scale (from [get.color.scale()])
#' @param mod_info Table with RNA modification data (from MODOMICS)
#' @param mods User-defined mapping of modifications to symbols
#'
#' @return Matrix of HTML table cells
def make_coverage_table(oligo_data, osm_data, nuc_length, color_scale, colors, mod_info, mods={}):
    """
    Generate a table of "stacked bars" representing oligonucleotides.
    
    :param oligo_data: Oligonucleotide data (from mzTab)
    :param osm_data: Spectrum-match data (from mzTab)
    :param nuc_length: Length of the full RNA sequence
    :param color_scale: Color scale (from get_color_scale())
    :param mod_info: Table with RNA modification data (from MODOMICS)
    :param mods: User-defined mapping of modifications to symbols
    :return: Matrix of HTML table cells
    """
    html_esc = {"&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;"}
    
    oligo_data = oligo_data.sort_values(by=["start", "end", "sequence"], ascending=[True, False, True])
    oligo_seqs = oligo_data["sequence"].unique()
    split_seqs = {seq: split_sequence(seq) for seq in oligo_seqs}  # Simulating split_sequence function
    
    table = np.full((1, nuc_length), "", dtype=object)
    
    for _, row in oligo_data.iterrows():
        seq = row["sequence"]
        count = (osm_data["sequence"] == seq).sum()
        bucket = next((i for i, val in enumerate(color_scale) if val >= count), None)
        color = colors[str(int(bucket))] if bucket is not None else "#FFFFFF"
        split_seq = split_seqs[seq]
        start, end = row["start"] - 1, row["end"] - 1
        oligo_length = end - start + 1 #FIXME
        
        if len(split_seq) != oligo_length:
            print(row)
            raise ValueError("Length mismatch")
        
        row_ind = next((i for i in range(len(table)) if all(cell == "" for cell in table[i, start:end+1])), None)
        
        if row_ind is None:
            table = np.vstack([table, np.full((1, nuc_length), "", dtype=object)])
            row_ind = len(table) - 1

        if start == end:
            table[row_ind, start] = '<td class="left right'
        else:
            table[row_ind, start] = '<td class="left"'
            table[row_ind, end] = '<td class="right"'
        
        if oligo_length > 2:
            table[row_ind, start+1:end] = '<td class="inner"'
        
        for j, nucleotide in enumerate(split_seq):
            if re.match(r"^p?[ACGU]p?$", split_seq[j]):
                mod = ""
            else:
                if nucleotide in mods:
                    char = mods[nucleotide]
                elif "|" in nucleotide:
                    char = "&nbsp;"
                else:
                    nucleotide = nucleotide.rstrip("p")
                    pos = mod_info[mod_info["short_name"] == nucleotide].index
                    char = "@" if pos.empty else mod_info.loc[pos[0], "html_abbrev"]
                mod = f'<div class="mod" title="{nucleotide}">{char}</div>' if char else ""
            table[row_ind, start + j] += (f' style="background-color:{color}" title="spectral count: {count}">'
                                         f'{mod}</td>\n')
    
    return table

#' Return the HTML page header for the coverage plot
def get_html_header():
    """Return the HTML page header for the coverage plot."""
    return """<!doctype html>
<html>
<head>
<style>
body {
  font-family: sans-serif;
}
table, th, td {
  border: none;
  border-spacing: 0px 4px;
}
.scale {
  height: 10px;
  width: 40px;
}
.nums {
  writing-mode: sideways-lr;
  font-size: 50%;
}
.left {
  border-left: 2px solid white;
}
.right {
  border-right: 2px solid white;
}
.left, .inner, .right {
  height: 10px;
}
.mod {
  max-height: 10px;
  text-align: center;
  background-color: green;
  font-weight: bold;
  color: white;
  font-size: 10px;
  line-height: 10px;
}
</style>
</head>
<body>
"""

#' Return HTML code for the color scale (from [get_color_scale()])
def get_scale_html(color_scale, colors):
    """Return HTML code for the color scale."""
    html = "<table style=\"text-align:center\">\n<tr>\n<td>Spectral counts:</td>\n"
    lower = [color_scale[0]] + [x + 1 for x in color_scale[:-1]]
    for i, color in enumerate(color_scale):
        lo = int(lower[i])
        hi = int(color)
        text = f"{lo}-{hi}" if lo != hi else str(hi)
        html += f"<td>{text}</td>\n"
    html += "</tr>\n<tr>\n<td/>\n"
    for color in range(len(color_scale)):
        html += f"<td class=\"scale\" style=\"background-color:{colors[str(color)]}\"></td>\n"
    html += "</tr>\n</table>\n"
    return html

#' Helper function to generate HTML code for the coverage plot of a single RNA
#'
#' @param accession Database accession of the RNA
#' @param mztab1 First mzTab file containing RNA data
#' @param mztab2 Optional second mzTab file containing RNA data
#' @param labels Labels to show for two mzTab files
#' @param break_at Add line breaks after this many bases in the sequence
#' @param color_scale Color scale (from [get_color_scale()])
#' @param mod_info Table with RNA modification data (from MODOMICS)
#' @param description Description of the RNA
#'
#' @return HTML code
def make_coverage_html_single(accession, mztab1, mztab2=None, labels=None,
                               break_at=float('inf'), color_scale=None, colors=None, mod_info=None, description=None):
    """Helper function to generate HTML code for the coverage plot of a single RNA."""
    if accession in mztab1['NUC']['accession'].values:
        seq = mztab1['NUC'].loc[mztab1['NUC']['accession'] == accession, "opt_sequence"].values[0]
        split_seq = split_sequence(seq)
    else:
        seq = mztab2['NUC'].loc[mztab2['NUC']['accession'] == accession, "opt_sequence"].values[0]
        split_seq = split_sequence(seq)
    
    nums = [''] * len(split_seq)
    fives = list(range(5, len(split_seq) + 1, 5))
    for f in fives:
        nums[f - 1] = str(f)
    
    # Generate coverage table
    oligo_data1 = mztab1['OLI'][mztab1['OLI']['accession'] == accession]
    coverage_table1 = make_coverage_table(oligo_data1, mztab1['OSM'], len(split_seq), color_scale, colors, mod_info)
    covered = sum([any(col != "" for col in coverage_table1[:, i]) for i in range(len(split_seq))])
    coverage1 = covered / len(split_seq)

    if mztab2 is not None:
        oligo_data2 = mztab2['OLI'][mztab2['OLI']['accession'] == accession]
        coverage_table2 = make_coverage_table(oligo_data2, mztab2['OSM'], len(split_seq), color_scale, colors, mod_info)
        covered = sum([any(col != "" for col in coverage_table2[:, i]) for i in range(len(split_seq))])
        coverage2 = covered / len(split_seq)
        row_labels = [f"{label}:" for label in labels]
    else:
        coverage_table2 = coverage_table1
        row_labels = ["", ""]

    # write header
    html = f"<h1>{accession}</h1>\n"
    if description is not None and not '' and not np.isnan(description):
        print(description)
        html += f"<p>{description}</p>\n"
    html += f"<p>Sequence coverage: {round(coverage1 * 100, 2)}%"
    if mztab2 is not None:
        html += f" ({labels[0]}) / {round(coverage2 * 100, 2)}% ({labels[1]})"
    html += "</p>\n<table>\n"
    
    # Handle line breaks
    break_at = min(break_at, len(split_seq))
    parts = (len(split_seq) + break_at - 1) // break_at  # Ceiling division
    row_label_style = "max-height:10px; line-height:10px; white-space:nowrap"
    
    for part in range(parts):
        current_cols = list(range(part * break_at, min((part + 1) * break_at, len(split_seq))))
        if mztab2 is not None:
            # Add oligonucleotides (1)
            html += f"<tr>\n<td style=\"{row_label_style}\" rowspan=\"{len(coverage_table1)}\"><b>{row_labels[0]}</b></td>\n"
            for i in range(len(coverage_table1), 0, -1):
                for j in current_cols:
                    s = coverage_table1[i - 1, j]
                    html += "<td/>\n" if s == "" else s
                html += "</tr>\n"

            # Add position numbers (1)
            html += f"<tr class=\"nums\">\n<td/>\n{''.join([f'<td>{nums[i]}</td>' for i in current_cols])}\n</tr>\n"
            html += "<tr>\n<td>Sequence:</td>\n"
        
        # Add RNA sequence
        html += f"<tr><td/>\n{''.join([f'<th>{split_seq[i]}</th>' for i in current_cols])}\n</tr>\n"
        
        # Add position numbers (2)
        html += f"<tr class=\"nums\">\n<td/>\n{''.join([f'<td>{nums[i]}</td>' for i in current_cols])}\n</tr>\n"
        
        # Add oligonucleotides (2)
        html += f"<tr>\n<td style=\"{row_label_style}\" rowspan=\"{len(coverage_table2)}\"><b>{row_labels[1]}</b></td>\n"
        for i in range(len(coverage_table2)):
            for j in current_cols:
                s = coverage_table2[i, j]
                html += "<td/>\n" if s == "" else s
            html += "</tr>\n"
        
        html += f"<tr>\n<td colspan=\"{len(current_cols) + 1}\"><hr/></td>\n</tr>\n"
    
    html += "</table>\n"
    return html

#' Generate a full HTML coverage plot for RNA data in one or two mzTab files
#'
#' @param mztab1 First mzTab file containing RNA data
#' @param mztab2 Optional second mzTab file containing RNA data
#' @param labels Labels to show for two mzTab files
#' @param break_at Add line breaks after this many bases in the sequence
#' @param mod_info Table with RNA modification data (from MODOMICS)
#'
#' @return HTML code
def make_coverage_html(mztab1, mztab2=None, labels=None, break_at=100, mod_info=read_rna_modifications()):
    """Generate a full HTML coverage plot for RNA data in one or two mzTab files."""
    html = get_html_header()
    count_range = [min(mztab1['OSM']['sequence'].value_counts()), max(mztab1['OSM']['sequence'].value_counts())]
    accessions = sorted(mztab1['NUC']['accession'].unique())
    
    if mztab2 is not None:
        if labels is None:
            raise ValueError("Labels must be provided when comparing two mzTab files")
        assert len(labels) == 2
        count_range = [min(count_range[0], min(mztab2['OSM']['sequence'].value_counts())),
                       max(count_range[1], max(mztab2['OSM']['sequence'].value_counts()))]
        accessions = sorted(set(accessions).union(mztab2['NUC']['accession'].unique()))
    
    color_scale, colors = get_color_scale(count_range[1], count_range[0])
    html += get_scale_html(color_scale,colors)
    
    for accession in accessions:
        desc = mztab1['NUC'][mztab1['NUC']['accession'] == accession]["description"].values[0]
        acc_html = make_coverage_html_single(accession, mztab1, mztab2, labels, break_at, color_scale, colors, mod_info, desc)
        html += acc_html
    
    html += "</body>\n</html>\n"
    return html