# Region SPECific Tool for continuous directed evolution (RESPECTevo)

Codes and processed data associated with RESPECTevo

You can run the code by either using git-clone or downloading the files directly as a ZIP archive from GitHub, then extracting them to your desired directory.
Once downloaded, navigate to the appropriate folder and run the scripts as needed. Installation requires less than one minute on a typical computer.

The processed data used for plotting and statistical analysis are stored in folders located near the corresponding code files. These output files are used to generate figures included in the article. If the code does not function as expected, please ensure that data paths are correctly specified. By default, relative paths are used, defined by the basicpath argument.

Demo files for testing NGS_process.py are included in folders adjacent to the script. In the readcount_to_xlsx function, the demo files are processed outputs from Varscan readcount applied to lacZ sequencing data. The resulting .xlsx file is saved in the plotting folder.

For the barcode_splitter_final function, the demo input is an artificially generated dataset. The function identifies and uses the last 4-letter barcode to split the input file. After execution, the function generates separate output files containing sequences with valid barcodes, correct lengths and primer binding site.

Running the Plotting and Statistics functions on the full deposited dataset takes approximately 10–20 minutes. The NGS_process function completes in under one minute on the demo file, but requires about 2 hours when processing the Real dataset.

Version of the code has been tested on
Python v.3.9.18
-pandas v 1.5.1
-numpy v 1.24.3
-matplotlib v 3.5.3
-openpyxl v 3. 0. 10
-seaborn v 0. 13. 2
-scipy v 1. 10. 1
-biopython v 1. 78
