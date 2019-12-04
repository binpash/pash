#! /home/users/cordier/.linuxbrew/bin/python3

class StateManager (object):
	"""
	Pipeline State Management Class (not yet implemented)
	Kind of like a flexible git 'light'. 
	"""
	def __init__ (self):
		return self

	#
	# Module Specific State Management
	#

	def create (self, args, status):
		"""
		Create State - Record Individual State to $module.state File
			@params
				args 	- Str - Required - List of state arguments
				status 	- Int - Required - C Status Code
			@return
				status  - Int -          - C Status Code
		"""
		return status

	def read (self, args, status):
		"""
		Read State - Check Whether State is Recorded in $module.state File
			@params
				args 	- Str - Required - List of state arguments
				status 	- Int - Required - C Status Code
			@return
				status  - Int -          - C Status Code
		"""
		return status

	def destroy (self, args, status):
		"""
		Destroy State - Delte Recorded State in $module.state File
			@params
				args 	- Str - Required - List of state arguments
				status 	- Int - Required - C Status Code
			@return
				status  - Int -          - C Status Code
		"""
		return status

	def reset (self, args, status):
		"""
		Reset State - Clear Module State(s) from $module.state File
			@params
				args 	- Str - Required - List of state arguments
				status 	- Int - Required - C Status Code
			@return
				status  - Int -          - C Status Code
		"""
		return status

	#
	# Global Pipeline State Management
	#

	def commit (self, args, status):
		"""
		Commit State - Commit Module State to Global pipeline.state File
			@params
				args 	- Str - Required - List of state arguments
				status 	- Int - Required - C Status Code
			@return
				status  - Int -          - C Status Code
		"""
		return status

	def committed (self, args, status):
		"""
		Committed State - Check Whether Module State is Committed to Global pipeline.state File
			@params
				args 	- Str - Required - List of state arguments
				status 	- Int - Required - C Status Code
			@return
				status  - Int -          - C Status Code
		"""
		return status

	def rollback (self, args, status):
		"""
		Rollback State - Rollback Module State in Global pipeline.state File
			@params
				args 	- Str - Required - List of state arguments
				status 	- Int - Required - C Status Code
			@return
				status  - Int -          - C Status Code
		"""
		return status

if __name__ == "__main__":

	# Imports
	import argparse

	# Parse Arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--method", type = str, help = "State Method ('put', 'has', 'commit')")
	parser.add_argument("-a", "--args", type = str, help = "State Arguments")
	parser.add_argumen("-s", "--status", type = int, help = "C Status Code")
	# Arg References
	argsDict = vars(parser.parse_args())
	method = argsDict["method"]
	args = argsDict["args"]
	status = argsDict["status"]

	# Init State Manager & Call Method with Args & Status Code
	state = StateManager()
	return state[method](args, status)

else:
	
	pass
