import numpy as np
import pandas as pd
import scvelo as scv
from anndata import AnnData
from scipy.sparse import csr_matrix

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

# Load spliced and unspliced abundances data
adata = scv.read_loom("velocity_files/Puck_190926_03_C0IH1.loom")
scv.pl.proportions(adata, save='info.png', show=False)

# Preprocess data
# scv.pp.filter_and_normalize(adata, min_shared_counts=10, n_top_genes=None) # Disabled n_top_genes because we manually subset markers
# scv.pp.remove_duplicate_cells(adata)
# scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Correct the names of bead barcodes of adata
adata.obs_names = pd.Index([barcode_id[21:35] for barcode_id in adata.obs_names])

# Add location data
bead_locations = pd.read_csv("velocity_files/Puck_190926_03_bead_locations.csv", index_col='barcodes')
list_coordinates = [[bead_locations['xcoord'][barcode], bead_locations['ycoord'][barcode]] for barcode in adata.obs_names]
adata.obsm['X_spatial'] = np.array(list_coordinates)

# Make variable names unique
adata.var_names_make_unique()

# Prepare data about subsets
markers = pd.read_csv("velocity_files/all_markers.csv")
beads = pd.read_csv("velocity_files/results_rctd_table.csv", index_col = 0)
value_counts = markers['types'].value_counts()
to_keep = value_counts[value_counts >= 10].index
types = set(markers[markers['types'].isin(to_keep)]['types'])

spatial_params = {'basis': 'spatial', 'xlim': (0,6000), 'ylim': (0,6000), 'size': 100}
base_plot = scv.pl.scatter(adata, **spatial_params, alpha = 0.2)

for type_combo in types:
    try:
        print(f'-----------------------------Processing {type_combo}-----------------------------')
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
    
        # type_marker_ids = [marker for marker in type_marker_ids if marker in adata.var_names]
        # type_bead_ids = [barcode for barcode in type_bead_ids if barcode in adata.obs_names]
    
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
    
        # Preprocess data
        print('--- Normalising')
        scv.pp.filter_and_normalize(type_data) # Removed shared counts because it was filtering out all genes
        scv.pp.filter_and_normalize(type_data_het)
        print('--- Removing duplicates')
        scv.pp.remove_duplicate_cells(type_data)
        scv.pp.remove_duplicate_cells(type_data_het)
        print('--- Calculating moments')
        scv.pp.moments(type_data, n_pcs=30, n_neighbors=30)
        scv.pp.moments(type_data_het, n_pcs=30, n_neighbors=30)
    
        # Show markers and beads that survived scv preprocessing
        if any(marker not in adata.var_names for marker in type_markers['ID'].values):
            print(f"{sum(marker not in adata.var_names for marker in type_markers['ID'].values)} out of {len(type_markers['ID'].values)} {type_combo} markers were filtered out by scvelo preprocessing")
        if any(barcode not in adata.obs_names for barcode in type_beads['barcode'].values):
            print(f"{sum(barcode not in adata.obs_names for barcode in type_beads['barcode'].values)} out of {len(type_beads['barcode'].values)} {type_combo} beads were filtered out by scvelo preprocessing")
    
        # Estimate RNA velocity
        print('--- Calculating velocity')
        scv.tl.velocity(type_data, mode = 'stochastic')
        scv.tl.velocity(type_data_het, mode = 'stochastic')
        print('--- Calculating velocity graph')
        scv.tl.velocity_graph(type_data)
        scv.tl.velocity_graph(type_data_het)
    
        # Project the velocities
        print('--- Calculating umap')
        scv.tl.umap(type_data)
        scv.tl.umap(type_data_het)
        print('--- Plot velocity on umap')
        scv.pl.velocity_embedding_stream(type_data, basis='umap', color='color_map', title = None, save=f'umap_{type_combo_sep[0]}_{type_combo_sep[1]}.png', show = False)
        scv.pl.velocity_embedding_stream(type_data_het, basis='umap', save=f'umap_{type_combo_sep[0]}_{type_combo_sep[1]}_het.png', show = False)
    
        # Project the velocities on spatial locations
        print('--- Plot velocity on spatial')
        velocity_plot = scv.pl.velocity_embedding_stream(type_data, color='types', save=f'spatial_{type_combo_sep[0]}_{type_combo_sep[1]}.png', show = False, **spatial_params)
        scv.pl.velocity_embedding_stream(type_data_het, save=f'spatial_{type_combo_sep[0]}_{type_combo_sep[1]}_het.png', show = False, **spatial_params)
        # base_plot + velocity_plot
    
        # # Store results
        # # Test whether running things on type_data is propagated back to adata
        # # otherwise, merge like this (or initialise with the first and merge everything afterwards)
        # scv.utils.merge(comb_data, type_data)
    
        # Save data to file
        # Read in again with scv.read('filename.loom', X_name = '')
        print('--- Write file')
        type_data.write_loom(f"velocity_files/velocity__{type_combo_sep[0]}_{type_combo_sep[1]}.loom", write_obsm_varm=True)
        print(f'-----------------------------Finished {type_combo}-----------------------------')
    except ValueError as e:
        print(f"{type_combo} failed, presumably due to the vstack bug.")
        print(e)
