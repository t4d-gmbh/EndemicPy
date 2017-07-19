# UIU Dynamics on Ant Colonies

**Please make sure that you use endemicPy version [v0.1.0]() if you want to run this example**
Click [here]() to download this version.


This is a minimal working example to simulate UIU dynamics on connection data of an ant colony.
The connection data is present in form of a [raw text file](./data_example.txt) and corresponds only to a small sample
of the data of a single colony used in our paper [Short-term activity cycles impede information transmission in ant colonies](https://doi.org/10.1371/journal.pcbi.1005527).

[uiu_single_trans.py](./uiu_single_trans.py) is the script to load the connection data and run UIU-porpagation on those connections.
It can be run from the command line.
For further details on how to run it open a console, navigate to the folder containing the
script and type:

    python uiu_single_trans.py --help

An exemplary command to run UIU-propagation on the example data could look like this:

    python uiu_single_trans.py -i data_example.txt -b 10. -m 0.01 -t 2. -s random -o simulation_output.txt -r 2

which will perform 2 runs, `-r 2` i.e. independent simulations, on the temporal network computed from the contact data in `data_example.txt`.
For each run a seed, i.e. initial carrier of the information will be chosen at random, `-s random`.
The output will be written in the file `simulation_output.txt` and will look similar to this:

```csv
#time, count
-- rep 0 --
1390257057.18, 0
1390257076.45, 3
1390257093.43, 4
1390257130.91, 3
1390257133.22, 4
...
-- rep 1 --
...
```

The content should be fairly easy to understand.
The line `-- rep x --` indicates that a new simulation starts.
Each column then reports the time and the count of informed individuals at that moment.

If more detailed output is needed, e.g. each transmission event should be reported, the option `-d` can be added along with specifying a file to write the detailed output in:

    python uiu_single_trans.py -i data_example.txt -b 10. -m 0.01 -t 2. -s random -o simulation_output.txt -r 2 -d detailed_output.txt

The content of `detailed_output.txt` will then look similar to this:

```csv
1390257057.18, 24, seed
1390257120.07, 24
1390257057.18, 25, seed
1390257057.21, 215, 25
1390257076.32, 592, 25
1390257093.41, 44, 215
1390257130.91, 44
1390257133.2, 44, 25
1390257140.61, 592
1390257189.71, 592, 25
1390257207.26, 25
1390257207.29, 215
1390257211.7, 25, 44
1390257252.1, 25
1390257264.79, 44
1390257282.97, 592
```

Each line has either 2 or 3 columns.
If only 2 columns are present, then this describes a recover event (so from I to U). 
As an example the line `1390257120.07, 24` writes:
At time `1390257120.07` individual `24` lost the information.
If 3 columns are given, then the line describes a transmission event.
So the line `1390257057.21, 215, 25` writes: At time `1390257057.21` individual `215` got a transmission from individual `25`.
An exception are the lines containing the `seed`.
Such a line basically indicates that a new simulation is reported now and gives the id of the initial carrier of the information.

If we consider the exemplary content of the `detailed_output.txt` from above, we see that in a first simulation individual `24` was the initial carrier but did not transmit to any other ant.
In the 2nd simulation individual `25` was the initial carrier, managed to transmit to `215`, `592`, `44` and again to `592` before it lost the information, and so on...

<!---    
Note that you can also explicitly specify the column names holding the start and stop times and the node ids/names such that you can load 
connection data with slightly different format. For the here presented example the column names were `Starttime`, `Stoptime`, `Tag1`, `Tag2` and the delimiter between
columns is `,`. These are also the default values for the [uiu_single_trans.py](./uiu_single_trans.py) script. If your column names deviate or you use another delimiter, like TAB you need
to provide the names of your columns and/or the delimiter as arguments to the scripts. 

As an example suppose you have a data file looking somehow like so:

```csv
Node1 Node2 T_start T_stop
a b 2 3
a c 2 4
...
```

The command to carry out 2 UIU simulations on this data could look like:

    python uiu_single_trans.py 
-->
