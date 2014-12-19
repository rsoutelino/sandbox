import os
import urllib2

top_level_url = "http://metocean.fugrogeos.com/marinha/"
target_url = "http://metocean.fugrogeos.com/marinha/Members/Data_month.csv"

# create a password manager
password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
password_mgr.add_password(None, top_level_url, "rsoutelino", "@#Upwelling")
handler = urllib2.HTTPBasicAuthHandler(password_mgr)
opener = urllib2.build_opener(handler)

# use the opener to fetch a URL
f = opener.open(target_url)
urllib2.install_opener(opener)

with open(os.path.basename(target_url), "wb") as local_file:
    local_file.write(f.read())
