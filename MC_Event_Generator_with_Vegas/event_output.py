#==========================================================#
# Process: e+e- -> Z/gamma -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.16
# Version: 1.0
#==========================================================#

class Event_Output(object):

	def output(self, i, p1, p2, p3, p4):
		with open("event.txt","a") as events: 
		    events.write("===================="+"event "+str(i)+"====================")
		    events.write("\n")
		    events.write("pem:	"+str(p1))
		    events.write("\n")
		    events.write("pep:	"+str(p2))
		    events.write("\n")
		    events.write("pmm:	"+str(p3))
		    events.write("\n")
		    events.write("pmp:	"+str(p4))
		    events.write("\n")
		    events.write("\n")
