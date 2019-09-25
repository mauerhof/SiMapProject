# Usage: make -j N name='NAME'
# N is the number of processors
# optional: NAME is a generic name, all models named NAME*.in will be run
# C. Morisset

# change the following line so that the string to the right of the equals sign 
# gives the path to the cloudy executable that you want to use for this run
CLOUDY = /home/evol/mauerhof/c17.00/source/cloudy.exe

# the following wildcard search returns the set of input scripts you want to run
SRC = $(wildcard ${name}*.in)

OBJ = $(SRC:.in=.out)
all: $(OBJ)

%.out: %.in
	$(CLOUDY) < $< > $@
# Notice the previous line has TAB in first column
