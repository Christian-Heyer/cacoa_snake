def get_input_files(wildcards, f_type):
    return config["samples"][wildcards.sample][f_type]


get_mtx_file = partial(get_input_files, f_type="mtx_file")
get_barcodes_file = partial(get_input_files, f_type="barcodes_file")
get_features_file = partial(get_input_files, f_type="features_file")
