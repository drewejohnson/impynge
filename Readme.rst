=======
impynge
=======
Simple neutron displacement damage calculator

Usage
-----
::

  usage: impynge.py [-h] [-e E_MIN] [-E E_MAX] [-S] [-N] {iron,graphite}      
  
  positional arguments:
    {iron,graphite}     Name of the target

  optional arguments:
    -h, --help          show this help message and exit
    -e E_MIN, --e-min E_MIN
                        Minimum neutron energy [eV]
    -E E_MAX, --e-max E_MAX
                        Maximum neutron energy [eV]
    -S, --save          Save figures to this directory
    -N, --no-show       Do not show figures.
  
