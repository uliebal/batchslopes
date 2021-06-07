batchslopes

Author: Ulf Liebal
Institution: iAMB, RWTH Aachen
Contact: ulf.liebal@rwth-aachen.de

batchslopes helps in identifying growth rates in growth curves from the 
growth profiler. Multiple experiments are bundled in a single csv file 
with a common time column. The tool bins the data and selects bin regions
with high correlation coefficient. These regions are further binned until
the correlation coefficient decreases and the top correlation coefficient
and the corresponding slope of the logarithmic OD are reported.
