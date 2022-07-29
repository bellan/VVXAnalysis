from os.path import expanduser
from os import environ, path

user = environ.get('USER', None)
home = expanduser("~")
personalFolder = home+"/VVXStuff" #Change your directory 

if(user is not None):
    site_location = "/eos/home-%s/%s/www/Analysis" %(user[0], user)
    if(path.exists(site_location)):
        personalFolder = site_location
    # eg. Analysis = "VVXAnalyzer", region = "CR3P1F"
