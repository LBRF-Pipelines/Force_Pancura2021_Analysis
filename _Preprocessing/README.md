# CFS Batch-Conversion Script

This subfolder contains Python code for batch-converting the raw Signal `.cfs` files into formats easily imported into R.

Each CFS file is converted into two separate files: a `.csv` file containing the metadata for each Signal frame (in the `Converted/frames` subfolder) and an Apache Feather (`.arrow`) file containing the actual signal data for each frame (in the `Converted/signals` subfolder).

## Dependencies

To run this script, you will need Python 3.7 or newer installed on your computer. You will also need to install the `pycfs-signal` and `pyarrow` Python packages, which you can do by running the following commands in a terminal window:

```bash
pip3 install git+https://github.com/a-hurst/pycfs-signal
pip3 install pyarrow
```


## Usage

To use the script, first place the raw CFS files to convert in the `Raw` subdirectory. Then, open a terminal window in this folder and run `python3 convert_cfs.py`. After a few seconds, the converted data will be available in a new folder named `Converted`.

