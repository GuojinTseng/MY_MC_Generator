#==========================================================#
# Process: pp -> Z -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.23
# Version: 1.0
#==========================================================#

#import guojin's module
import run_card, param_card

class Event_Output(object):

	def output_headers(self):
		with open("event.lhe", "a") as events:
			events.write("#************************************************************")
			events.write("\n")
			events.write("#*                 Guojin's Event Generator                 *")
			events.write("\n")
			events.write("#*                                                          *")
			events.write("\n")
			events.write("#*                *                       *                 *")
			events.write("\n")
			events.write("#*                  *                   *                   *")
			events.write("\n")
			events.write("#*                    * * G-J Tseng * *                     *")
			events.write("\n")
			events.write("#*                  *                   *                   *")
			events.write("\n")
			events.write("#*                *                       *                 *")
			events.write("\n")
			events.write("#*                                                          *")
			events.write("\n")
			events.write("#*                                                          *")
			events.write("\n")
			events.write("#*         VERSION 1.0                 2018-07-23           *")
			events.write("\n")
			events.write("#*                                                          *")
			events.write("\n")
			events.write("#*                       Find me at                         *")
			events.write("\n")
			events.write("#*                  guojintseng@pku.edu.cn                  *")
			events.write("\n")
			events.write("#*                                                          *")
			events.write("\n")
			events.write("#************************************************************")
			events.write("\n")
			events.write("#*                                                          *")
			events.write("\n")
			events.write("#*             Copyright (C) 2018 Guojin Tseng              *")
			events.write("\n")
			events.write("#*                                                          *")
			events.write("\n")
			events.write("#************************************************************")
			events.write("\n")
			events.write("\n")
			events.write("#**************************run_card**************************")
			events.write("\n")
			events.write("COM Energy in GeV: ")
			events.write(str(run_card.ECM))
			events.write("\n")
			events.write("Random Seed: ")
			events.write(str(run_card.seed))
			events.write("\n")
			events.write("Number of Sprinkled Points: ")
			events.write(str(run_card.N))
			events.write("\n")
			events.write("Number of Generated events: ")
			events.write(str(run_card.Nevents))
			events.write("\n")
			events.write("\n")
			events.write("#************************param_card**************************")
			events.write("\n")
			events.write("Z Boson Mass: ")
			events.write(str(param_card.MZ))
			events.write("\n")
			events.write("Z Boson Width: ")
			events.write(str(param_card.GAMMAZ))
			events.write("\n")
			events.write("ALPHA_EM @MZ: ")
			events.write(str(param_card.alpha))
			events.write("\n")
			events.write("G_Fermi Constant: ")
			events.write(str(param_card.G_Fermi))
			events.write("\n")
			events.write("sin^2(weinberg angle)")
			events.write(str(param_card.sw2))
			events.write("\n")
			events.write("\n")

	def output(self, i, p1, p2, p3, p4):
		with open("event.lhe","a") as events:
			events.write("<event>")
			events.write("\n")
			events.write("pq1:	"+str(p1))
			events.write("\n")
			events.write("pq1:	"+str(p2))
			events.write("\n")
			events.write("pmm:	"+str(p3))
			events.write("\n")
			events.write("pmp:	"+str(p4))
			events.write("\n")
			events.write("<\event>")
			events.write("\n")
			events.write("\n")