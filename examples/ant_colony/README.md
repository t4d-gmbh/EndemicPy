# UIU dynamics on an ant colony

This is a minimal working example to simulate UIU dynamics on connection data of an ant colony.
The connection data is present in form of a [raw text file](./data_example.txt) and corresponds only to a small sample
of the data of a single colony used in our paper [Short-term activity cycles impede information transmission in ant colonies](https://doi.org/10.1371/journal.pcbi.1005527).

[uiu_single_trans.py](./uiu_single_trans.py) is the script to load the connection data and run UIU-porpagation on those connections.
It can be run from the command line. For further details on how to run it open a console, navigate to the folder containing the
script and type:

    python uiu_single_trans.py --help

An exemplary command could look like so:

    python uiu_single_trans.py ...
    
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
