This is a makefile C code for creating the cone shape speckle pattern. 

In order to run the program a basic `Scatterers.txt` file is needed.  This
`Scatterers.txt` file is produced by another C code,
`random_scatterers_creater.c`, which will give us a position matrix of a
fixed number of scatterers positions.

There are several output files for this program.  

- A .h5 file which contains the cone speckle pattern. Using the command
	'h5topng' this file can be changed into .png file to view.
	
-	A `visting counting.txt` file which contains the information of the
	radiation counting for every scatterer. Using this .txt file, an
	localization drawing can be done by the matlab codes.  

-	A `free_path.txt` file which contains the accumulated path for each
	iteration in the program. Using the matlab code, an exponenatial decade
	fuction with the information of mean-free-path as the expectation. 
