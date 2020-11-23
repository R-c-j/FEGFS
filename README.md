## Article title

## Table of contents
* [Install](#Install)
	* [Install using `pip`](#Install using `pip`)
* [Usage](#Usage)
	* [example](#example)
	* [real data](#real data)
## Install
### Install using `pip`

Firstly, we suggest to build a new virtual environment through `Anaconda`:
```
conda create -n FFM python=3.7
```
Create and activate the virtual environment environment `FFM`:
```
conda activate FFM
```

## Usage
All functions of FFM can be found in the script folder `FFMC.py` Is running:
```
python FFMC.py -d <data_name> -c <num_clusters> -n <retention_ratio> -i <GO_Term_path> -e　<expression_matrix_path> -o <outputpath> -l <label_path>
```
### example
Run `FFMC` as an example in script:
```
python FFMC.py -d test -c 5 -n 0.4 -i ./example/GO_Term.xlsx -e ./example/test_count_matrix.csv -o .example/Term_matrix -l ./example/test_label.csv
```
### real data
Here we take `Pollen` as an example:
```
python FFMC.py -d test -c 5 -n 0.4 -i ./example/GO_Term.xlsx -e ./example/test_count_matrix.csv -o .example/Term_matrix -l ./example/test_label.csv
```

