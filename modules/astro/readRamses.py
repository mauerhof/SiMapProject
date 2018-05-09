def look_for(directory, parameter):
	
	f = open(directory, 'r')
	parameters = f.read().split()
	
	for i in range(len(parameters)):
		if parameters[i] == parameter:                  # looks for the word inside the look_for parenthesis
			number = eval(parameters[i+2])           # equalizes the variable parameter to the value after the equal sign in parameters.txt
			break
		else:
			continue

	return number



def read_units(directory, timestep):
	
	tsString = '_' + format(timestep, "05")
	directory = directory+'output'+tsString+'/info'+tsString+'.txt'
	
	f = open(directory, 'r')
	parameters = f.read().split()
	
	for i in range(len(parameters)):
		if parameters[i] == 'unit_l':                 
			unit_length = eval(parameters[i+2])
			unit_d = eval(parameters[i+5]) 
			unit_time = eval(parameters[i+8])   
			break

	return [unit_length,unit_d, unit_time]
	
	
def read_groups(directory, timestep):
	
	tsString = '_' + format(timestep, "05")
	directory = directory+'output'+tsString+'/info_rt'+tsString+'.txt'
	
	f = open(directory, 'r')
	parameters = f.read().split()
	
	groups = []
	
	for i in range(len(parameters)):
		if parameters[i] == 'groupL0':
			break
	i += 3
	while(parameters[i] != 'groupL1'):
		groups.append(eval(parameters[i]))
		i += 1
	
	return groups        
			  


