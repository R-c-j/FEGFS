## Article title

## Table of contents
* [Install](#Install)
	* [Install using `pip`](#Install)
	* [package](#Package)
* [Usage](#Usage)
	* [Example](#Example)
	* [Real data](#Example) 
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
### Package

| package | version |
| ---- | ---- |
| pandas | 0.25.1 |
| numpy | 1.16.5 |
| scikit-learn | 0.21.3 |
| matplotlib | 3.3.1 |

## Usage
All functions of FFM can be found in the script folder `FFMC.py` Is running:
```
python FFMC.py -d <data_name> -c <num_clusters> -n <retention_ratio> -i <GO_Term_path> -e　<expression_matrix_path> -o <outputpath> -l <label_path>
```
### Example
Run `FFMC` as an example in script:
```
python FFMC.py -d test -c 5 -n 0.4 -i ./example/GO_Term.xlsx -e ./example/test_count_matrix.csv -o ./example/Term_matrix -l ./example/test_label.csv
```
### Real data
Here we take `Pollen` as an example:
```
python FFMC.py -d pollen -c 11 -n 0.4 -i ./pollen/pollen_GO_Term.xlsx -e ./pollen/pollen_count_matrix.csv -o ./pollen/pollen_Term_matrix -l ./pollen/pollen_label.csv
```

