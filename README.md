# EndemicPy
## Transmission dynamics in various structures

EndemicPy, or simply _endemic_ is a python package under development aiming to simulate a vast range of transmission dynamics on various host structure models.


## How to use _endemic_
_endemic_ can either be installed or one can simply copy the [endemic](endemic/) folder into the same folder as the python script that is using the package.
If you decide to simply copy the folder, you can ignore the installation step below.

### Versions
_endemic_ has different releases that were (and are) developed for various projects. 

*You need to use the correct release for the type of simulations you plan to do.*
Using the wrong version might result in errors or incoherent output.

**If you are using EndemicPy for your own work, please make sure that you correctly refer to this package and to the research article related to the specific release you are using.**

To learn which release is the right one for you, here is a list of the published projects along with the used release:

- **[v0.1.0](https://github.com/j-i-l/EndemicPy/releases/tag/v0.1.0)** for [Short-term activity cycles impede information transmission in ant colonies](https://doi.org/10.1371/journal.pcbi.1005527)

- **[v0.3.0](https://github.com/j-i-l/EndemicPy/releases/tag/v0.3.0)** for [Host population structure impedes reversion to drug sensitivity after discontinuation of treatment]( https://doi.org/10.1371/journal.pcbi.1005704)

### Dependencies
The only non-standard python package _endemic_ uses in [numpy](http://www.numpy.org/), so please make sure that you have a recent numpy version installed.

### Installation
To install the package you might want to setup a [virtualenv](https://virtualenv.pypa.io/en/stable/) which is not a requirement but a recommendation. 
To install the package, simply open a console, `cd` into the EndemicPy folder and type:

    python setup.py install

If everything works fine, you can now simply import _endemic_ as a package in your python scripts.

### Examples
Check out the [examples](examples/) for more information on how EndemicPy is used in different projects.
