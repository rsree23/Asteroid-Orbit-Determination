# For all files
Please see the attached paper, which will help with understanding the purpose of the programs and the
goal of the whole project. This was all done at The Summer Science Program in Astrophysics, 2018. www.ssp.org.

# For Orbit_Determination.py specifically

My code will output all permutations (1,2,5), (1,3,5), (1,4,5) with the roots for each.
On my input file and in the code I refer to the observations as 1, 2, 3, 4, 5.  
The actual correspondance to the 7 times my team went observing are as follows:  
Obs Night: 1, 2, 3, 4, 5, 6, 7  
Obs Number: 1, 2, 3, Cloudy Night (none), 4, Cloudy Night (none), 5

Also, in my input file, if a time was 06:00:00, I wrote it as 6 0 0, since having the two zeroes or a zero in front
of a number confuses python and causes errors. Also, due to the way I read from a file in my code, I took away any colons (I need
spaces to split values, not colons).

With regard to the print output, my code will output orbital elements not only for all observation night permutations,
but also for all positive, real Scalar Equation of Lagrange roots, even if their respective orbital elements far off (the
range is negative). You will see that this is the case for two roots in my 1, 2, 5 permutation, but the third has a positive
range and works.

I also have included the option to differentially correct orbital elements.

You will see a user input yes or no (y/n) for each observation night to differentially correct.

If the range was negative, there is no option to differentially correct those orbital elements (can already disregard
knowing that they had a negative rho scalar).
