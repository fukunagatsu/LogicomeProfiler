# Logicome Profiler
Logicome Profiler exhaustively detects statistically significant triplet logic relationships from a binary matrix dataset.

## Version
Version 1.1.0 (2020/01/30)

## Usage
./LogicomeProfiler <input_file> <output_file> <input_file> <output_file> <controlling procedure> <familywise error rate or false discovery rate>

If you would like to use the FWER and the FDR for the multiple testing criteria, please set the controlling procedure to 0 and 1, respectively. For example, when you would like to use the FWER criteria and set the significance level to 0.05, the command line is as follows:

    ./LogicomeProfiler test_input.txt test_output.txt 0 0.05

## Input file format
The number of samples and items are described in the first line. In the second line, the name of each sample is listed. After the third line, the presence of absence of each item in each sample is described.

## Example of input file format
    3 5
    item_name sample1 sample2 sample3
    item1 1 0 0
    item2 1 0 1
    item3 0 1 0
    item4 0 1 1
    item5 1 1 0

## Reference
Tsukasa Fukunaga and Wataru Iwasaki. "Logicome Profiler: Exhaustive detection of statistically significant logic relationships from comparative omics data." PLoS One, 15, e0232106 (2020)
