# OTIS: Order Through Informed Swapping

OTIS generates FCC and BCC atomic lattices with short-range order (SRO). 

OTIS is run on the command line as 

    python otis.py input_file
    
Python3 is required as is the NumPy module.

The input file is structured as

    lattice_type
    N1 N2 N3
    lattice_parameter
    n
    elem1 elem2 elem3 ... 
    c1 c2 c3 ...
    input_type
    a11 a12 a13 ...
    a21 a22 a23 ...
    a31 a32 a33 ...
    ...
    output_file
    
The first line specifies the lattice type. The currently available options are `BCC` or `FCC`. The simulation cell is created to align with the primitive directions for that lattice type. For `BCC`, the primitive directions are e<sub>1</sub> = \[111̅], e<sub>2</sub> = \[1̅11], and e<sub>3</sub> = \[11̅1]. For `FCC`, the primitive directions are e<sub>1</sub> = \[110], e<sub>2</sub> = \[011], and e<sub>3</sub> = \[101].

The second line gives the number of primitive cells in primitive direction. There are `N1` cells in the e<sub>1</sub> direction, `N2` cells in the e<sub>2</sub> direction, and `N3` cells in the e<sub>3</sub> direction

`lattice_parameter` is the lattice parameter of the *primitive* cell i.e. the distance between nearest neighbors. This is only used when outputting the lattice positions. The units are not important and can be whatever you need for your use case. 

`n` is the number of component elements.

`elem1`, `elem2`, and so on are the element names. There should be `n` entries on this line.

`c1`, `c2`, and so on are the concentrations of each element. There should be `n` entries on this line.

`input_type` can be either `wc` or `pij`. `wc` indicates that the subsequent lines will contain Warren-Cowley paramters, while `pij` indicates the subsequent lines will contain a probability matrix. 

The next `n` lines contain an `n` x `n` matrix of either Warren-Cowley parameters or bond probabilities. The Warren-Cowley parameters are defined as 

α<sub>ij</sub> = (p<sub>ij</sub> - c<sub>j</sub>) / (δ<sub>ij</sub> - c<sub>j</sub>)

where p<sub>ij</sub> is the probability of the nearest neighbor of an i-type atom being j-type. The ordering of rows and columns should correspond with the order of element names given previously. The code will check that the probabilities add to one for each element.

The final line contains the filename for the lattice to be output after swapping is completed. The output file has an `xyz` format in which each line contains the element name followed by the x, y, and z position. 

