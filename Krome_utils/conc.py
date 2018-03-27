filenames = ['outstromgrenSi0.dat','outstromgrenSi1.dat','outstromgrenSi2.dat','outstromgrenSi3.dat']
with open('outstromgrenSi.dat', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
