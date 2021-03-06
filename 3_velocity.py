import numpy as np
import pandas as pd
import scvelo as scv
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import string

scv.settings.verbosity = 1  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

# Load spliced and unspliced abundances data
adata = scv.read_loom("velocity_files/Puck_190926_03_C0IH1.loom")
scv.pl.proportions(adata, save='info.png', show=False)

# Correct the names of bead barcodes of adata
adata.obs_names = pd.Index([barcode_id[21:35] for barcode_id in adata.obs_names])

# Add location data
bead_locations = pd.read_csv("velocity_files/Puck_190926_03_bead_locations.csv", index_col='barcodes')
list_coordinates = [[bead_locations['xcoord'][barcode], bead_locations['ycoord'][barcode]] for barcode in adata.obs_names]
adata.obsm['X_spatial'] = np.array(list_coordinates)

# Make variable names unique
adata.var_names_make_unique()

# Prepare data about subsets
markers = pd.read_csv("velocity_files/all_markers_050.csv")
beads = pd.read_csv("velocity_files/results_rctd_table.csv", index_col = 0)
value_counts = markers['types'].value_counts()
to_keep = value_counts[value_counts >= 5].index
types = set(markers[markers['types'].isin(to_keep)]['types'])


for type_combo in types:
    # Start processing data
    print(f'----------------------------- Processing {type_combo} -----------------------------')
    type_combo_sep = type_combo.split('/')
    type_set = [type_combo_sep[0], type_combo, type_combo_sep[1]]

    # Filter data to markers and beads for current types
    type_markers = markers[markers['types'] == type_combo]
    type_beads = beads[beads['types'].isin(type_set)]
    
    # Check for matching beads and variable names
    if any(marker not in adata.var_names for marker in type_markers['ID'].values):
        print(f"{sum(marker not in adata.var_names for marker in type_markers['ID'])} out of {len(type_markers['ID'])} {type_combo} markers can't be found in the velocyto data")
        # print([marker for marker in markers['ID'] if marker not in adata.var_names])
    if any(barcode not in adata.obs_names for barcode in type_beads['barcode'].values):
        print(f"{sum(barcode not in adata.obs_names for barcode in type_beads['barcode'])} out of {len(type_beads['barcode'])} {type_combo} beads can't be found in the velocyto data")
        # print([barcode for barcode in beads['barcode'] if barcode not in adata.obs_names])

    # Ensure type beads and markers are in velocyto data
    type_beads = type_beads[type_beads['barcode'].isin(adata.obs_names)]
    type_markers = type_markers[type_markers['ID'].isin(adata.var_names)]

    # Subset adata to current types
    type_data = None
    type_data = adata[type_beads['barcode'].values, type_markers['ID'].values].copy()
    
    # Add bead type metadata
    type_data.obs = pd.concat([type_data.obs, type_beads[['barcode', 'spot_class', 'type1', 'type2', 'types', 'prop1', 'prop2']].set_index('barcode')], axis = 1)
    color_mapping = dict(zip(type_set, list(range(1, len(type_set)+1))))
    type_data.obs['color_map'] = [color_mapping[type_val] for type_val in type_data.obs['types']]
    
    # Get beads and markers to modify (bead doublets and markers split by dominance)
    doublet_beads = type_beads[type_beads['spot_class'].isin(['doublet_certain', 'doublet_uncertain', 'reject'])]
    doublet_indices = [i for i, barcode in enumerate(type_data.obs_names) if barcode in doublet_beads['barcode'].values]
    type1_markers = type_markers[type_markers['pct.1'] > type_markers['pct.2']]
    type2_markers = type_markers[type_markers['pct.1'] < type_markers['pct.2']]
    
    # Normalise counts
    spl = type_data.layers['spliced'].toarray()
    unspl = type_data.layers['unspliced'].toarray()

    for i, mrk in enumerate(type_data.var_names):
        if mrk in type1_markers['ID'].values:
            spl[doublet_indices, i] = spl[doublet_indices, i]/doublet_beads['prop1']
            unspl[doublet_indices, i] = unspl[doublet_indices, i]/doublet_beads['prop1']
        elif mrk in type2_markers['ID'].values:
            spl[doublet_indices, i] = spl[doublet_indices, i]/doublet_beads['prop2']
            unspl[doublet_indices, i] = unspl[doublet_indices, i]/doublet_beads['prop2']
        else:
            print(f"Could not find {mrk} in the marker data. Something went wrong, check whether all data sources align.")

    type_data.layers['spliced'] = csr_matrix(spl)
    type_data.layers['unspliced'] = csr_matrix(unspl)

    # Create heterotypic beads-only object
    type_data_het = type_data[type_data.obs['types'] == type_combo, :].copy()

    type_data_a = type_data[type_data.obs['types'].isin(type_set[:2]), type_data.var_names.isin(type1_markers['ID'].values)]
    type_data_b = type_data[type_data.obs['types'].isin(type_set[1:]), type_data.var_names.isin(type2_markers['ID'].values)]

    # Set the stage for plots
    colour_params = [{'color': '#35b779'}, {'color': 'types', 'palette': ['#fde725', '#35b779']}, {'color': 'types', 'palette': ['#35b779', '#31688e']}]
    fig, axs = plt.subplots(2, 3, frameon = False)


    # Try calculating velocities for the full set (heterotypic and both homotypic) and heterotypic only
    index_to_text = {0: 'heterotypic', 1: 'heterotypic and type 1', 2: 'heterotypic and type 2'}
    for extent_marker, type_data_obj in enumerate([type_data_het, type_data_a, type_data_b]):
        try:
            print(f'------ Processing {index_to_text[extent_marker]} velocities for {type_combo}')
            # Preprocess data
            print('--- Removing duplicates')
            scv.pp.remove_duplicate_cells(type_data_obj)
            print('--- Calculating moments')
            scv.pp.moments(type_data_obj, n_pcs=30, n_neighbors=30, use_highly_variable = False)

            # Show markers and beads that survived scv preprocessing
            # if any(marker not in adata.var_names for marker in type_markers['ID'].values):
            #     print(f"{sum(marker not in adata.var_names for marker in type_markers['ID'].values)} out of {len(type_markers['ID'].values)} {type_combo} markers were filtered out by scvelo preprocessing")
            if any(barcode not in adata.obs_names for barcode in type_beads['barcode'].values):
                print(f"{sum(barcode not in adata.obs_names for barcode in type_beads['barcode'].values)} out of {len(type_beads['barcode'].values)} {type_combo} beads were filtered out by scvelo preprocessing")

            # Estimate RNA velocity
            print('--- Calculating velocity')
            scv.tl.velocity(type_data_obj, mode = 'stochastic', filter_genes = False, min_r2 = -1, use_highly_variable = False)
            print('--- Calculating velocity graph')
            scv.tl.velocity_graph(type_data_obj)

            # Project the velocities
            print('--- Calculating umap')
            scv.tl.umap(type_data_obj)

            # Add subplots to grid (in axs)
            print('--- Plot velocity on umap')
            axs[1, extent_marker] = scv.pl.velocity_embedding_stream(type_data_obj, basis = 'umap', size = 50, alpha = 0.8, title = '', show = False, **colour_params[extent_marker], ax = axs[1, extent_marker], legend_loc = 'right margin')

            # Project the velocities on spatial locations (with adata scatter base to show overall plot)
            print('--- Plot velocity on spatial')
            axs[0, extent_marker] = scv.pl.scatter(adata, basis = 'spatial', xlim = (0, 6000), ylim = (0, 6000), size = 100, alpha = 0.05, show = False, ax = axs[0, extent_marker])
            axs[0, extent_marker] = scv.pl.velocity_embedding_stream(type_data_obj, basis = 'spatial', xlim = (0,6000), ylim = (0,6000), size = 100, title = '', show = False, **colour_params[extent_marker], ax = axs[0, extent_marker], legend_loc = 'right margin')

        except ValueError as e:
            # Remove failed attempts from plot
            axs[0, extent_marker].remove()
            axs[1, extent_marker].remove()
            # Print error catch message
            print(f"ERROR: {index_to_text[extent_marker].capitalize()} {type_combo} failed.")
            print(e)
            if "all the input array dimensions" in str(e):
                print(
                    """This error is likely due to a vstack call in scv.tl.velocity -> compute_stochastic -> leastsq_generalised. 
                    This happens when scvelo keeps only one gene for velocity because the velocities are too similar,
                    and this single gene becomes a mismatched shape matrix due to unexpected conversions between (n,) vectors and (n,1) single-row matrices.
                    scvelo does not check for single genes at that stage, which should be raised with them."""
                    )

    print('--- Finish and save plot')
    # Add subplot letters
    for n, ax in enumerate(axs.flat):
        ax = ax.text(-0.1, 1.1, string.ascii_uppercase[n], transform=ax.transAxes, size=14, weight='bold')
    # Save composite plot
    fig.savefig(f'figures/{type_combo_sep[0]}_{type_combo_sep[1]}_ab.png', bbox_inches = 'tight')
    
    print(f'------------------------------ Finished {type_combo} ------------------------------')

