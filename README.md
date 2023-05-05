# TMP-ALIGN
T&MP-Align computes some alignments and similarity scores between two combinations of organism+pathway.
Furthermore, it outputs information regarding individual pathways that can be useful for it's analysis.

There tool can be launched in some different ways depending of user's criteria:

#1 - Analysis of individual pathways. You have to run:
                python main organism pathway
For example:    python main ggo 00020

Output is saved in 'Results/Path_Info' folder by default. You'll be able to study T-invariants, graphs, etc

#2 - Comparation and alignment between two pathways. You have to run:
                python main organism1 organism2 pathway
For example:    python main fri fre 00020

Alignment info is saved in 'Results/Align' folder by default. You can see T-invariant alignment, paths and
reaction alignment, etc.

#3 - Batch alignment. You can compute a large set of alignments without having to enter data every time
using a CSV file. CSV location has to be entered in settings.py file. In addition, you can specify
CSV's separator, column from where the tool has to extract the org. code, etc. It's all explained in settings.py

To execute a batch alignment, run the tool without specifying any argument, like so:
                python main
                
->  T&MP-Align automatically retrieves all the files and information that it needs every time so you don't have
    to worry about downloading and putting files into the adequate directory. All you need is an active internet
    connection with accessible port 80.
    
---------------------------------------------------------------------------------------------------------------

SETTINGS.PY FILE

This file includes the tool's configuration and useful parameters. You can set:
-Path's limit to compute
-Aligment weights
-Directories to save the results
-Paths computation parameters
-CSV file for batch alignment
...

Note that this file is in pure Python language, check it's syntax after editing it.

