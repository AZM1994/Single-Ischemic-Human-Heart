# SigNet

SigNet is a package to study genetic mutational processes.
Check out our [theoretical background page](documentation/theoretical_background.md) for further information on this topic.
As of now, it contains 3 solutions:

- **[SigNet Refitter](documentation/signet_refitter.md)**: Tool for signature decomposition.
- **[SigNet Generator](documentation/signet_generator.md)**: Tool for realistic mutational data generation.
- **[SigNet Detector](documentation/signet_detector.md)**: Tool for mutational vector out-of-distribution detection.


This is the official code implementation of the paper: [Mutational signature decomposition with deep neural networks reveals origins of clock-like processes and hypoxia dependencies](https://www.biorxiv.org/content/10.1101/2023.12.06.570467v1).
By [Claudia Serrano](https://www.linkedin.com/in/claudia-serrano-colome-440aa9182/), [Oleguer Canal](https://github.com/OleguerCanal), et all.

## Readme contents

You can use SigNet in 3 different ways depending on your workflow:

1. **[Python Package](#python-package)**
   1. Python Package Installation
   2. Python Package Usage

2. **[Command Line Interface](#command-line-interface)** (CLI)

3. **[Source Code](#source-code)**
   1. Downloading Source Code
   2. Code-Basics


## Python Package
Recommended if you want to integrate SigNet as part of your python workflow, or intending to re-train models on custom data with limited ANN architectural changes.
You can install the python package running:

```BASH
pip install signaturesnet
```

Once installed, you can run Signet Refitter like so:

```python
import pandas as pd
from signaturesnet.modules.signet_module import SigNet

# Read your mutational data
mutations = pd.read_csv("your_input", header=0, index_col=0)

# Load & Run signet
signet = SigNet(opportunities_name_or_path="your_normalization_file")
results = signet(mutation_dataset=mutations)

# Extract results
w, u, l, c, _ = results.get_output()

# Store results
results.save(path='Output', name="this_experiment_filename")

# Plot figures
results.plot_results(save=True)
```

For a more usage examples: Check out the [examples folder](examples/):
   - [refitter_example.py](examples/refitter_example.py) for a usage example.
   - [generator_example.py](examples/generator_example.py) for a usage example.
   - [detector_example.py](examples/detector_example.py) for a usage example.

**NOTE**: _It is recommended that you work on a [custom python virtualenvironment](https://virtualenv.pypa.io/en/latest/) to avoid package version mismatches._


## Command Line Interface

Recommended if only interested in running SigNet modules independently and **not** willing to retrain models or change the source code.<br>
**NOTE**: _This option is only tested on Debian-based Linux distributions_. Steps:

1. Download the [signaturesnet exectuable](TODOlink_to_executable)
2. Change directory to wherever you downloaded it: `cd <wherever/you/downloaded/the/executable/>` 
3. Make it executable by your user: `sudo chmod u+x signaturesnet`

__Refitter:__

The following example shows how to use [SigNet Refitter](documentation/signet_refitter.md).


```BASH
cd <wherever/you/downloaded/the/executable/>
./signaturesnet refitter  [--input_format {counts, bed, vcf}]
                   [--input_data INPUTFILE]
                   [--reference_genome REFGENOME]
                   [--normalization {None, exome, genome, PATH_TO_ABUNDANCES}] 
                   [--only_nnls ONLYNNLS]
                   [--cutoff CUTOFF]
                   [--output OUTPUT]
                   [--plot_figs False]
```

- `--input_format`: Name of the format of the input. The default is 'counts'. Please refer to [Mutations Input](documentation/input_output_formats.md##Mutations-Input) for further details.

- `--input_data`: Path to the file containing the mutational counts. Please refer to [Mutations Input](documentation/input_output_formats.md##Mutations-Input) for further details.

- `--reference_genome`: Name or path to the reference genome. Needed when input_format is bed or vcf.

- `--normalization`: As the INPUTFILE contain counts, we need to normalize them according to the abundances of each trinucleotide on the genome region we are counting the mutations.
  - Choose `None` (default): If you don't want any normalization.
  - Choose `exome`:  If the data that is being input comes from Whole Exome Sequencing. This will normalize the counts according to the trinucleotide abundances in the exome.
  - Choose `genome`: If the data comes from Whole Genome Sequencing.
  - Set a `PATH_TO_ABUNDANCES` to use a custom normalization file. Please refer to [Normalization Input](documentation/input_output_formats.md##Mutations-Input) for further details on the input format.

- `--only_nnls`: Whether to use NNLS mode only (the finetuner is not run). Default: `False`.

- `--cutoff`: Cutoff to be applied to the final weights. Default: 0.01.

- `--output` Path to the folder where all the output files (weights guesses and figures) will be stored. By default, this folder will be called "Output" and will be created in the current directory. Please refer to [SigNet Refitter Output](documentation/input_output_formats.md##Signet-Refitter-Output) for further details on the output format.

- `--plot_figs` Whether to generate output plots or not. Possible options are `True` or `False`.


__Detector:__

```BASH
cd <wherever/you/downloaded/the/executable/>
./signaturesnet detector  [--input_data INPUTFILE]
                   [--normalization {None, exome, genome, PATH_TO_ABUNDANCES}] 
                   [--output OUTPUT]
```

(Same arguments as before)

__Generator:__

```BASH
cd <wherever/you/downloaded/the/executable/>
./signaturesnet generator  [--n_datapoints INT]
                    [--output OUTPUT]
```

- `--n_datapoints`: Number of signature weight combinations to generate.


## Source Code

Is the option which gives more flexibility.
Recommended if you want to play around with the code, re-train custom models or [do contributions](documentation/).

### Downloading Source Code

Clone the repo and install it as an editable pip package like so:

```BASH
git clone git@github.com:weghornlab/SigNet.git
cd SigNet
pip install -e .
```

Refer [here](documentation/code_structure.md) for the project code organization.
