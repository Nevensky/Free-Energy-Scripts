#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
#plt.style.use("ggplot")
matplotlib.rcParams.update({'font.size':7 })

import click

#sendgrid
import sendgrid
import os
from sendgrid.helpers.mail import *
import base64

umbrella_dir=os.getcwd()

@click.command()
@click.option('-mdp',default=umbrella_dir+'/MDP',help='MDP file dirs.')

def inputs(mdp):
	'''Analysis script for Umbrella sampling. \n

Developed by Neven Golenic | neven.golenic@gmail.com'''
	umbrella_title="DHK/MeCN 0.8nm window sep in range (0.6 - 1.8)nm"

	histogram = np.genfromtxt("analysis/histogram.xvg",delimiter="\t",skip_header=17)[:,:-1]
	pmf = np.genfromtxt("analysis/profile.xvg",delimiter="\t",skip_header=17)
	#print(histogram)


	fig, (ax1,ax2)= plt.subplots(2, 1)



	ax1.plot(pmf[:,0],pmf[:,1])
	ax1.set_title("Potential of Mean Force")
	ax1.set_ylabel("E [kJ/mol]")
	ax1.set_xlabel("Reaction coordinate [nm]")

	group_hist =[]
	len_hist = np.shape(histogram)[1]
	for k in range(len_hist):
		if k!=0:
			ax2.plot(histogram[:,0],histogram[:,k])

	ax2.set_title("Umbrella histograms")
	ax2.set_ylabel("count")
	ax2.set_xlabel("Reaction coordinate [nm]")

	plt.tight_layout()

	plt.savefig("analysis/wham.pdf")


	def encodemail():
		#encode pull_distances.pdf as pull_distances_attachment in base64
		with open("analysis/pull_distances.pdf", "rb") as f:
			pull_distances_attachment = base64.b64encode(f.read())
		#adapt for sending
		pull_distances_attachment = str(pull_distances_attachment)[2:-1]

		# encode wham.pdf as wham_attachment in base64
		with open("analysis/wham.pdf", "rb") as f:
			wham_attachment = base64.b64encode(f.read())
		#adapt for sending
		wham_attachment = str(wham_attachment)[2:-1]

		# encode wham.pdf as wham_attachment in base64
		with open("system.top", "rb") as f:
			topology_attachment = base64.b64encode(f.read())
		#adapt for sending
		topology_attachment = str(topology_attachment)[2:-1]
		return wham_attachment,pull_distances_attachment,topology_attachment

	def sendmail(wham_attachment,pull_distances_attachment,topology_attachment):
		sg = sendgrid.SendGridAPIClient(apikey='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')
		from_email = Email("umbrella@mail.nevensky.net")
		to_email = Email("neven.golenic@gmail.com")
		subject = "Umbrella sampling report - "+umbrella_title
		content = Content("text/plain", "Umbrella sampling report. Attachment contains plots of PMF and histograms.")
		mail = Mail(from_email, subject, to_email, content)
		personalization = Personalization()
		personalization.add_to(Email("ghorvat@chem.pmf.hr","Gordan Horvat"))
		personalization.add_cc(Email("neven.golenic@chem.pmf.hr","Neven Golenic"))
		mail.add_personalization(personalization)
		attachment = Attachment()
		attachment.set_content(wham_attachment)
		attachment.set_type("application/pdf")
		attachment.set_filename("wham.pdf")
		attachment.set_disposition("attachment")
		attachment.set_content_id("WHAM report")
		mail.add_attachment(attachment)
		attachment2 = Attachment()
		attachment2.set_content(pull_distances_attachment)
		attachment2.set_type("application/pdf")
		attachment2.set_filename("pull_distances.pdf")
		attachment2.set_disposition("attachment")
		attachment2.set_content_id("Pull CONF distances")
		mail.add_attachment(attachment2)
		attachment3 = Attachment()
		attachment3.set_content(topology_attachment)
		attachment3.set_type("application/top")
		attachment3.set_filename("system.top")
		attachment3.set_disposition("attachment")
		attachment3.set_content_id("Topology")
		mail.add_attachment(attachment3)
		response = sg.client.mail.send.post(request_body=mail.get())
		print(response.status_code)
		print(response.body)
		print(response.headers)


	wham_attachment,pull_distances_attachment,topology_attachment=encodemail()
	sendmail(wham_attachment,pull_distances_attachment,topology_attachment)

if __name__ == '__main__':
	inputs()
