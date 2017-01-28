#!/usr/bin/env python3
# Updates BAR scripts from master branch on github repo.
import requests
import json
import sys
root = sys.path[0]
current_version ='bar_v1.0.0'

r = requests.get('https://api.github.com/repos/Nevensky/Free-Energy-Skripte/releases')
print("Checking for new version of BAR script.")
if(r.ok):
	repoItem = json.loads(r.text or r.content)
	for item in repoItem:
		if 'bar' in item['tag_name']:
			latest_version = item['tag_name']
			break
	if latest_version != current_version:
		print("Updating to: "+latest_version)
		os.system('wget https://github.com/Nevensky/Free-Energy-Skripte/releases/tag/'+latest_version)
	else:
		print("Nothing to update.")