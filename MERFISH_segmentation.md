# MERFISH Data Analysis
Last updated by Chan (12/7/2023)
Anyone working on cell segmentation of this data should feel free to modify this document.
## Changelog:

## General Workflow:
1. Generate MERSCOPE data
2. Use vizgen-postprocessing tool to resegment data
3. Apply modification by Baysor segmentation
4. Reload segments into .vzg file and into MERSCOPE Visualizer

## Useful Links:
 - https://github.com/kharchenkolab/Baysor#news-060--2023-04-13
 - https://www.nature.com/articles/s41587-021-01044-w
 - https://github.com/Vizgen/vizgen-postprocessing
 - https://vizgen.github.io/vizgen-postprocessing/index.html

## Notes
Sample dataset can be found in /lab/solexa_reddien/MERFISH_cellpose_Baysor_Sample/ which contains all necessary files and scripts for analysis.

## Getting Started
### Downloading Vizgen Post-processing tool
Notes: Currently (12/7/2023) VPT requires >=python 3.9 and <python3.11
1. Set up environment: (requires conda)
```sh
conda activate py310
pip install vpt[all]
```
2. Modify parameter file if necessary (sample in /lab/solexa_reddien/Chan/Baysor_Test/cellpose_default_1_ZLevel.csv)
```sh

  "experiment_properties": {
    "all_z_indexes": [0, 1, 2, 3, 4, 5, 6],
    "z_positions_um": [1.5, 3, 4.5, 6, 7.5, 9, 10.5]
  },
  "segmentation_tasks": [
    {
      "task_id": 0,
      "segmentation_family": "Cellpose",
      "entity_types_detected": [
        "cell"
      ],
      "z_layers": [
        3
      ],
      "segmentation_properties": {
        "model": "cyto2",
        "model_dimensions": "2D",
        "custom_weights": null,
        "version": "latest"
      },
      "task_input_data": [
        {
          "image_channel": "PolyT",
          "image_preprocessing": [
            {
              "name": "normalize",
              "parameters": {
                "type": "CLAHE",
                "clip_limit": 0.01,
                "filter_size": [
                  100,
                  100
                ]
              }
            }
          ]
        },
        {
          "image_channel": "DAPI",
          "image_preprocessing": [
            {
              "name": "normalize",
              "parameters": {
                "type": "CLAHE",
                "clip_limit": 0.01,
                "filter_size": [
                  100,
                  100
                ]
              }
            }
          ]
        }
      ],
      "segmentation_parameters": {
        "nuclear_channel": "DAPI",
        "entity_fill_channel": "PolyT",
        "diameter": 70,
        "flow_threshold": 0.95,
        "mask_threshold": -5.5,
        "minimum_mask_size": 500
      },
      "polygon_parameters": {
        "simplification_tol": 2,
        "smoothing_radius": 10,
        "minimum_final_area": 500
      }
    },
    {
      "task_id": 1,
      "segmentation_family": "Cellpose",
      "entity_types_detected": [
        "cell"
      ],
      "z_layers": [
        3
      ],
      "segmentation_properties": {
        "model": "nuclei",
        "model_dimensions": "2D",
        "custom_weights": null,
        "version": "latest"
      },
      "task_input_data": [
        {
          "image_channel": "DAPI",
          "image_preprocessing": [
            {
              "name": "normalize",
              "parameters": {
                "type": "CLAHE",
                "clip_limit": 0.01,
                "filter_size": [
                  100,
                  100
                ]
              }
            }
          ]
        }
      ],
      "segmentation_parameters": {
        "nuclear_channel": "DAPI",
        "entity_fill_channel": null,
        "diameter": 55,
        "flow_threshold": 0.8,
        "mask_threshold": -3,
        "minimum_mask_size": 500
      },
      "polygon_parameters": {
        "simplification_tol": 2,
        "smoothing_radius": 10,
        "minimum_final_area": 500
      }
    }
  ],
  "segmentation_task_fusion": {
    "entity_fusion_strategy": "harmonize",
    "fused_polygon_postprocessing_parameters": {
      "min_distance_between_entities": 1,
      "min_final_area": 500
    }
  },
  "output_files": [
    {
      "entity_types_output": [
        "cell"
      ],
      "files": {
        "run_on_tile_dir": "result_tiles",
        "mosaic_geometry_file": "cellpose_mosaic_space.parquet",
        "micron_geometry_file": "cellpose_micron_space.parquet",
        "cell_metadata_file": "cellpose_cell_metadata.csv"
      }
    }
  ]
}
```
3. Start segmentation - I ran the following script:
```sh
sbatch ./SCRIPT.sh
```
where SCRIPT.sh contains the details about running the code on slurm:
```sh
#!/bin/bash
# Configuration values for SLURM job submission.
# One hash before SBATCH is not a comment, but two hashes are.
#SBATCH --job-name=cellposeMERFISH       # friendly name for job.
#SBATCH --nodes=1                     # ensure cpus are on one node
#SBATCH --ntasks=1                    # run a single task
#SBATCH --cpus-per-task=24             # number of cpus/threads requested.
#SBATCH --mem=120gb                     # memory requested.
#SBATCH --partition=nvidia-t4-20      # partition (queue) to use
#SBATCH --output error-%j.out  # name of output file.  %j inserts jobid
#SBATCH --gres=gpu:1                  # This is needed to actually access a gpu
#SBATCH --mail-type=ALL              # enable this to have slurm email you when a job starts / finishes
#SBATCH --mail-user=YOUR_EMAIL@wi.mit.edu
vpt run-segmentation --segmentation-algorithm /PATH_TO_PARAMETER_FILE.json --input-images /PATH_TO_IMAGE_FOLDER/ --input-micron-to-mosaic /PATH_TO/micron_to_mosaic_pixel_transform.csv --output-path /PATH_TO_OUTPUT/
```
This process will take a while - also modify cpu, allocated memory, and partition if necessary. 

For my samples, this seems to take around ~5 hours. 
For region1: 135 tiles (from 742277x35385 pixel images), used ~4gb memory but gave it 24 cores on GPU (T4) - took around 5:50hrs
For region3: 96 tiles (from 53905x33533 pixel images), used ~4gb memory but gave it 12 cores on GPU (A4000) - took around 4:22hrs. 

4. Partition transcripts into cells
```sh
vpt partition-transcripts --input-boundaries /PATH/TO/CELLPOSE/RESULT/cellpose_micron_space.parquet --input-transcripts /PATH/TO/detected_transcripts.csv --output-entity-by-gene /PATH/TO/OUTPUT/cell_by_gene.csv --output-transcripts /PATH/TO/OUTPUT/detected_transcripts.csv
```
This process generates a cell matrix and modifies the transcript list to associate each with a cell_id
This process took around ~3mins with 8gb memory and 4 CPU

5. OPTIONAL: Sum signal intensities for each cell for all co-detections - can be useful when using sequential rounds potentially.
```sh
vpt derive-entity-metadata --input-boundaries /PATH/TO/CELLPOSE/RESULT/cellpose_micron_space.parquet --input-entity-by-gene /PATH/TO/OUTPUT/cell_by_gene.csv --output-metadata /PATH/TO/OUTPUT/cell_metadata.csv

vpt sum-signals --input-images="/PATH/TO/IMAGES/mosaic(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif" --input-boundaries /PATH/TO/CELLPOSE/RESULTS/cellpose_micron_space.parquet --input-micron-to-mosaic /PATH/TO/FILE/micron_to_mosaic_pixel_transform.csv --output-csv /PATH/TO/OUTPUT/sum_signals.csv
```
The metadata will contain information about location, size, and shape of cells in the data if you want to perform e.g. neighbor analysis downstream or size filter cells.
Took around ~5 mins with 4 cores and ~1gb memory for region1

The sum-signals will have information about each stains brightness for all cells including DAPI and polyT.
This took around ~2:21 hrs with 4gb and 4 cores

6. Update  .vzg file for visualization on MERSCOPE Visualizer

I would recommend generating a file even when performing Baysor segmentation below to generate two files to compare segmentations. This will be useful in fine-tuning parameters down the line. 

Note: This might require >50GB of free space since it will need to unpack all of the data, so make sure there is adequate storage whereever you do this.

```sh
vpt update-vzg --input-vzg /PATH/TO/ORIGINAL/experiment.vzg --input-boundaries /PATH/TO/CELLPOSE/cellpose_micron_space.parquet --input-entity-by-gene /PATH/TO/OUTPUT/cell_by_gene.csv --input-metadata /PATH/TO/OUTPUT/cell_metadata.csv --output-vzg /PATH/TO/OUTPUT/experiment_cellpose.vzg
```

### Setting up Baysor
1. Install Baysor
```sh
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git")); Pkg.build()'
```
2. Add to PATH
```sh
echo $PATH
export PATH="/home/USERNAME/.julia/bin:$PATH"
```
3. Generate config.toml. Example:
```sh
[data]
x = "global_x" # Name of the x column in the input data. Default: "x"
y = "global_y" # Name of the y column in the input data. Default: "y"
z = "global_z" # Name of the y column in the input data. Default: "z"
gene = "transcript_id" # Name of gene column in the input data. Default: "gene"
force_2d = true # Ignores z-column in the data if it is provided
min_molecules_per_gene = 1 # Minimal number of molecules per gene. Default: 1
exclude_genes = "" # Comma-separated list of genes or regular expressions to ignore during segmentation. Example: 'Blank*,MALAT1'
min_molecules_per_cell = 1 # Minimal number of molecules for a cell to be considered as real. It's an important parameter, as it's used to infer several other parameters. Default: 3
min_molecules_per_segment = 0 # Minimal number of molecules in a segmented region, required for this region to be considered as a possible cell. Default: min-molecules-per-cell / 4
confidence_nn_id = 0 # Number of nearest neighbors to use for confidence estimation. Default: min-molecules-per-cell / 2 + 1
[segmentation]
scale = -1.0 # Negative values mean it must be estimated from `min_molecules_per_cell`
scale_std = "25%" # Standard deviation of scale across cells. Can be either number, which means absolute value of the std, or string ended with "%" to set it relative to scale. Default: "25%"
estimate_scale_from_centers = true # Use scale estimate from DAPI if provided. Default: true
n_clusters = 9 # Number of clusters to use for cell type segmentation. Default: 4
prior_segmentation_confidence = 0.2 # Confidence of the prior segmentation. Default: 0.2
iters = 500 # Number of iterations for the cell segmentation algorithm. Default: 500
n_cells_init = 0 # Initial number of cells
nuclei_genes = "" # Comma-separated list of nuclei-specific genes. If provided, `cyto-genes` has to be set, as well.
cyto_genes = "" # Comma-separated list of cytoplasm-specific genes. If provided, `nuclei-genes` has to be set, as well.
# The parameters below are not supposed to be changed normally
new_component_weight = 0.2 # Prior weight of assignment a molecule to new component. Default: 0.2
new_component_fraction = 0.3 # Fraction of distributions, sampled at each stage. Default: 0.3
[plotting]
gene_composition_neigborhood = 0 # Number of neighbors (i.e. 'k' in k-NN), which is used for gene composition visualization. Larger numbers leads to more global patterns. Default: estimate from min-molecules-per-cell
min_pixels_per_cell = 15 # Number of pixels per cell of minimal size, used to estimate size of the final plot. For most protocols values around 7-30 give enough visualization quality. Default: 15
```
4. Run preview experiment. This will give a sense of transcript spread and noise within the sample. 
```sh
baysor preview -c /PATH_TO/config.toml /PATH_TO_TRANSCRIPTS/detected_transcripts.csv
```
This will generate a html file with some useful statistics.

5. Run Baysor segmentation:

To run a de novo segmentation experiment, run the following:
```sh
baysor run -s 50 -x global_x -y global_y -g gene -m 5 -c /PATH_TO/config.toml --save-polygons=geojson -p --count-matrix-format tsv /PATH/TO/detected_transcripts.csv
```
I set -s to 50 since it would correspond roughly to ~5um in size (since 1 pixel ~ .108um). 
To use a previous cellpose segmentation result:
```sh
baysor run [-s 50 -x global_x -y global_y --gene transcript_id] --prior-segmentation-confidence 0.7 --save-polygons=geojson -p -c config.toml /PATH_TO/detected_transcripts.csv :cell_id
```
For prior segmentation confidence, default is 0.2, but for sparse data, >0.7 can work better. Higher results in less modification of prior segmentation. 

In cell_id column, cells should have integer assignments, with 0 as unassigned transcripts. 

6. Convert segmentation into VPT compatible parquet file
```sh
vpt convert-geometry --input-boundaries /PATH/TO/BAYSOR/segmentation_polygons.json --output-boundaries /PATH/TO/OUTPUT/
```

7. Update .vzg file (Can also get metadata? Untested)

```sh
vpt update-vzg --input-vzg /PATH/TO/ORIGINAL/experiment.vzg --input-boundaries /PATH/TO/BAYSOR/RESULT/baysor_cellpose_micron_space.parquet --input-entity-by-gene /PATH/TO/OUTPUT/cell_by_gene.csv --output-vzg /PATH/TO/OUTPUT/experiment_cellpose.vzg
```
## Helpful stuff

```sh
baysor run --help
```
The following list should be helpful in setting up a de novo segmentation or utilizing the cellpose output from VPT. 
```sh
 run <args> [options] [flags]
Args
  <coordinates>                                             CSV file with coordinates of molecules and gene type
  [prior_segmentation]                                      Image or a MAT file with segmentation mask (either boolean
                                                            or component indexing) or CSV column with integer
                                                            segmentation labels. If it's the column name, it should be
                                                            preceded ':' symbol (e.g. :cell)
Options
  -c, --config <config.toml>                                TOML file with a config
  -x, --x-column <x>                                        Name of x column. Overrides the config value.
  -y, --y-column <y>                                        Name of y column. Overrides the config value.
  -z, --z-column <z>                                        Name of z column. Overrides the config value.
  -g, --gene-column <gene>                                  Name of gene column. Overrides the config value.
  -m, --min-molecules-per-cell <m>                          Minimal number of molecules for a cell to be considered as
                                                            real. It's an important parameter, as it's used to infer
                                                            several other parameters. Overrides the config value.
  -s, --scale <s>                                           Scale parameter, which suggest approximate cell radius for
                                                            the algorithm. Must be in the same units as x and
                                                            y molecule coordinates. Overrides the config
                                                            value. Sets estimate-scale-from-centers to
                                                            false.
  --scale-std <ss>                                          Standard deviation of scale across cells. Can be either
                                                            number, which means absolute value of the std, or string
                                                            ended with '%' to set it relative to scale (default: "25%")
  --n-clusters <nc>                                         Number of molecule clusters, i.e. major cell types. Depends
                                                            on protocol resolution, but should not be too high. In most
                                                            cases something between 3 and 15 should work well. (default:
                                                            4)
  --prior-segmentation-confidence <p>                       Confidence of the prior_segmentation results.
                                                            Value in [0; 1]. If you want the final segmentation not
                                                            contradicting to prior_segmentation, set it to 1.
                                                            Otherwise, if you assume errors in prior_segmentation,
                                                            values in [0.2-0.7] allow flexibility for the algorithm.
                                                            (default: 0.2)
  -o, --output <path>                                       Name of the output file or path to the output directory
                                                            (default: "segmentation.csv")
  --save-polygons <format>                                  Save estimated cell boundary polygons to a file with a
                                                            specified format. Only 'GeoJSON' format is
                                                            currently supported.
  --count-matrix-format <format>                            Storage format of the segmentec cell count matrix. Either
                                                            'loom' or 'tsv' (default: 'loom')
Flags
  -p, --plot                                                Save pdf with plot of the segmentation
  --no-ncv-estimation                                       Turns off neighborhood composition vectors estimation
  -h, --help                                                Print this help message.
```


