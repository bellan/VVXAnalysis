from os.path import expanduser
from os import environ

user = environ.get('User', None)
home = expanduser("~")
personalFolder = home+"/VVXStuff" #Change your directory 

if(user is not None):
    personalFolder = "/eos/home-%s/%s/www/Analysis/last" %(user[0], user, "EXT")
    # eg. Analysis = "VVXAnalyzer", region = "CR3P1F"
